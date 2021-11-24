"""Module for representation and analysis of MS measurements"""

from ..measurements import Measurement
from ..spectra import Spectrum
from ..plotters.ms_plotter import MSPlotter, STANDARD_COLORS
from ..exceptions import SeriesNotFoundError, QuantificationError
from ..constants import (
    AVOGADROS_CONSTANT,
    BOLTZMAN_CONSTANT,
    STANDARD_TEMPERATURE,
    STANDARD_PRESSURE,
    DYNAMIC_VISCOSITIES,
    MOLECULAR_DIAMETERS,
    MOLAR_MASSES,
)
from ..data_series import TimeSeries, ValueSeries
from ..db import Saveable
import re
import numpy as np


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    extra_column_attrs = {
        "ms_meaurements": {
            "mass_aliases",
            "signal_bgs",
        },
    }

    def __init__(
        self,
        name,
        mass_aliases=None,
        signal_bgs=None,
        tspan_bg=None,
        calibration=None,
        **kwargs,
    ):
        """Initializes a MS Measurement

        Args:
            name (str): The name of the measurement
            calibration (dict): calibration constants whereby the key
                corresponds to the respective signal name.
            mass_aliases (dict): {mass: mass_name} for any masses that
                do not have the standard 'M<x>' format used by ixdat.
            signal_bgs (dict): {mass: S_bg} where S_bg is the background signal
                in [A] for the mass (typically set with a timespan by `set_bg()`)
            calibration (ECMSCalibration): A calibration for the MS signals
            tspan_bg (timespan): background time used to set masses
        """
        super().__init__(name, **kwargs)
        self.calibration = calibration  # TODO: Not final implementation
        self.mass_aliases = mass_aliases or {}
        self.signal_bgs = signal_bgs or {}
        self.tspan_bg = tspan_bg

    def __getitem__(self, item):
        """Try standard lookup, then check if item is a flux or alias for a mass"""
        try:
            return super().__getitem__(item)
        except SeriesNotFoundError:
            if item in self.mass_aliases:
                return self[self.mass_aliases[item]]
            if item.startswith("n_"):  # it's a flux!
                mol = item.split("_")[-1]
                return self.get_flux_series(mol)
            else:
                raise

    def set_bg(self, tspan_bg=None, mass_list=None):
        """Set background values for mass_list to the average signal during tspan_bg."""
        mass_list = mass_list or self.mass_list
        tspan_bg = tspan_bg or self.tspan_bg
        for mass in mass_list:
            t, v = self.grab(mass, tspan_bg)
            self.signal_bgs[mass] = np.mean(v)

    def reset_bg(self, mass_list=None):
        """Reset background values for the masses in mass_list"""
        mass_list = mass_list or self.mass_list
        for mass in mass_list:
            if mass in self.signal_bgs:
                del self.signal_bgs[mass]

    def grab_signal(
        self,
        signal_name,
        tspan=None,
        t_bg=None,
        removebackground=False,
        include_endpoints=False,
    ):
        """Returns t, S where S is raw signal in [A] for a given signal name (ie mass)

        Args:
            signal_name (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            t_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            removebackground (bool): Whether to remove a pre-set background if available
                Defaults to False. (Note in grab_flux it defaults to True.)
            include_endpoints (bool): Whether to ensure tspan[0] and tspan[-1] are in t
        """
        time, value = self.grab(
            signal_name, tspan=tspan, include_endpoints=include_endpoints
        )

        if t_bg is None:
            if removebackground and signal_name in self.signal_bgs:
                return time, value - self.signal_bgs[signal_name]
            return time, value

        else:
            _, bg = self.grab(signal_name, tspan=t_bg)
            return time, value - np.average(bg)

    def grab_cal_signal(self, signal_name, tspan=None, t_bg=None):
        """Returns a calibrated signal for a given signal name. Only works if
        calibration dict is not None.

        Args:
            signal_name (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            t_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
        """
        # TODO: Not final implementation.
        # FIXME: Depreciated! Use grab_flux instead!
        if self.calibration is None:
            print("No calibration dict found.")
            return

        time, value = self.grab_signal(signal_name, tspan=tspan, t_bg=t_bg)

        return time, value * self.calibration[signal_name]

    def grab_flux(
        self,
        mol,
        tspan=None,
        tspan_bg=None,
        removebackground=True,
        include_endpoints=False,
    ):
        """Return the flux of mol (calibrated signal) in [mol/s]

        Args:
            mol (str or MSCalResult): Name of the molecule or a calibration thereof
            tspan (list): Timespan for which the signal is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            removebackground (bool): Whether to remove a pre-set background if available
                Defaults to True.
        """
        if isinstance(mol, str):
            if not self.calibration or mol not in self.calibration:
                raise QuantificationError(
                    f"Can't quantify {mol} in {self}: "
                    f"Not in calibration={self.calibration}"
                )
            mass, F = self.calibration.get_mass_and_F(mol)
        elif isinstance(mol, MSCalResult):
            mass = mol.mass
            F = mol.F
        else:
            raise TypeError("mol must be str or MSCalResult")
        x, y = self.grab_signal(
            mass,
            tspan=tspan,
            t_bg=tspan_bg,
            removebackground=removebackground,
            include_endpoints=include_endpoints,
        )
        n_dot = y / F
        return x, n_dot

    def grab_flux_for_t(
        self,
        mol,
        t,
        tspan_bg=None,
        removebackground=False,
        include_endpoints=False,
    ):
        """Return the flux of mol (calibrated signal) in [mol/s] for a given time vec

        Args:
            mol (str): Name of the molecule.
            t (np.array): The time vector along which to give the flux
            tspan_bg (tspan): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            removebackground (bool): Whether to remove a pre-set background if available
        """
        t_0, y_0 = self.grab_flux(
            mol,
            tspan_bg=tspan_bg,
            removebackground=removebackground,
            include_endpoints=include_endpoints,
        )
        y = np.interp(t, t_0, y_0)
        return y

    def get_flux_series(self, mol, tspan=None):
        """Return a ValueSeries with the calibrated flux of mol during tspan"""
        t, n_dot = self.grab_flux(mol, tspan=tspan)
        tseries = TimeSeries(
            name="n_dot_" + mol + "-t", unit_name="s", data=t, tstamp=self.tstamp
        )
        vseries = ValueSeries(
            name="n_dot_" + mol, unit_name="mol/s", data=n_dot, tseries=tseries
        )
        return vseries

    def integrate_signal(self, mass, tspan, tspan_bg, ax=None):
        """Integrate a ms signal with background subtraction and evt. plotting

        TODO: Should this, like grab_signal does now, have the option of using a
            background saved in the object rather than calculating a new one?

        Args:
            mass (str): The mass for which to integrate the signal
            tspan (tspan): The timespan over which to integrate
            tspan_bg (tspan): Timespan at which the signal is at its background value
            ax (Axis): axis to plot on. Defaults to None
        """
        t, S = self.grab_signal(mass, tspan=tspan, include_endpoints=True)
        if tspan_bg:
            t_bg, S_bg_0 = self.grab_signal(
                mass, tspan=tspan_bg, include_endpoints=True
            )
            S_bg = np.mean(S_bg_0) * np.ones(t.shape)
        else:
            S_bg = np.zeros(t.shape)
        if ax:
            if ax == "new":
                fig, ax = self.plotter.new_ax()
            ax.fill_between(t, S_bg, S, color=STANDARD_COLORS[mass], alpha=0.2)
        return np.trapz(S - S_bg, t)

    @property
    def mass_list(self):
        """List of the masses for which ValueSeries are contained in the measurement"""
        return [self.as_mass(col) for col in self.series_names if self.is_mass(col)]

    def is_mass(self, item):
        if re.search("^M[0-9]+$", item):
            return True
        if item in self.mass_aliases.values():
            return True
        return False

    def as_mass(self, item):
        if re.search("^M[0-9]+$", item):
            return item
        else:
            try:
                return next(k for k, v in self.mass_aliases.items() if v == item)
            except StopIteration:
                raise TypeError(f"{self} does not recognize '{item}' as a mass.")

    @property
    def plotter(self):
        """The default plotter for MSMeasurement is MSPlotter"""
        if not self._plotter:
            self._plotter = MSPlotter(measurement=self)
        return self._plotter


class MSCalResult(Saveable):
    """A class for a mass spec calibration result.

    TODO: How can we generalize calibration? I think that something inheriting directly
        from saveable belongs in a top-level module and not in a technique module
    """

    column_attrs = {"name", "mol", "mass", "cal_type", "F"}

    def __init__(
        self,
        name=None,
        mol=None,
        mass=None,
        cal_type=None,
        F=None,
    ):
        super().__init__()
        self.name = name or f"{mol} at {mass}"
        self.mol = mol
        self.mass = mass
        self.cal_type = cal_type
        self.F = F

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name={self.name}, mol={self.mol}, "
            f"mass={self.mass}, F={self.F})"
        )

    @property
    def color(self):
        return STANDARD_COLORS[self.mass]


class MSInlet:
    """A class for describing the inlet to the mass spec

    Every MSInlet describes the rate and composition of the gas entering a mass
    spectrometer. The default is a Spectro Inlets EC-MS chip.
    """

    def __init__(
        self,
        *,
        l_cap=1e-3,
        w_cap=6e-6,
        h_cap=6e-6,
        gas="He",
        T=STANDARD_TEMPERATURE,
        p=STANDARD_PRESSURE,
        verbose=True,
    ):
        """Create a Chip object given its properties

        Args:
            l_cap (float): capillary length [m]. Defaults to design parameter.
            w_cap (float): capillary width [m]. Defaults to design parameter.
            h_cap (float): capillary height [m]. Defaults to design parameter.
            p (float): system pressure in [Pa] (if to change from that in medium)
            T (float): system temperature in [K] (if to change from that in medium)
            gas (str): the gas at the start of the inlet.
            verbose (bool): whether to print stuff to the terminal
        """
        self.verbose = verbose
        self.l_cap = l_cap
        self.w_cap = w_cap
        self.h_cap = h_cap
        self.p = p
        self.T = T
        self.gas = gas  # TODO: Gas mixture class. This must be a pure gas now.

    def calc_n_dot_0(
        self, gas=None, w_cap=None, h_cap=None, l_cap=None, T=None, p=None
    ):
        """Calculate the total molecular flux through the capillary in [s^-1]

        Uses Equation 4.10 of Daniel's Thesis.

        Args:
            w_cap (float): capillary width [m], defaults to self.w_cap
            h_cap (float): capillary height [m], defaults to self.h_cap
            l_cap (float): capillary length [m], defaults to self.l_cap
            gas (dict or str): the gas in the chip, defaults to self.gas
            T (float): Temperature [K], if to be updated
            p (float): pressure [Pa], if to be updated
        Returns:
            float: the total molecular flux in [s^-1] through the capillary
        """

        if w_cap is None:
            w_cap = self.w_cap  # capillary width in [m]
        if h_cap is None:
            h_cap = self.h_cap  # capillary height in [m]
        if l_cap is None:
            l_cap = self.l_cap  # effective capillary length in [m]
        if T is None:
            T = self.T
        if p is None:
            p = self.p

        pi = np.pi
        eta = DYNAMIC_VISCOSITIES[gas]  # dynamic viscosity in [Pa*s]
        s = MOLECULAR_DIAMETERS[gas]  # molecule diameter in [m]
        m = MOLAR_MASSES[gas] * 1e-3 / AVOGADROS_CONSTANT  # molecule mass in [kg]

        d = ((w_cap * h_cap) / pi) ** 0.5 * 2
        # d = 4.4e-6  #used in Henriksen2009
        a = d / 2
        p_1 = p
        lambda_ = d  # defining the transitional pressure
        # ...from setting mean free path equal to capillary d
        p_t = BOLTZMAN_CONSTANT * T / (2 ** 0.5 * pi * s ** 2 * lambda_)
        p_2 = 0
        p_m = (p_1 + p_t) / 2  # average pressure in the transitional flow region
        v_m = (8 * BOLTZMAN_CONSTANT * T / (pi * m)) ** 0.5
        # a reciprocal velocity used for short-hand:
        nu = (m / (BOLTZMAN_CONSTANT * T)) ** 0.5

        # ... and now, we're ready for the capillary equation.
        #   (need to turn of black and flake8 for tolerable format)
        # fmt: off
        #   Equation 4.10 of Daniel Trimarco's PhD Thesis:
        N_dot = (                                                               # noqa
            1 / (BOLTZMAN_CONSTANT * T) * 1 / l_cap * (                         # noqa
                (p_t - p_2) * a**3 * 2 * pi / 3 * v_m + (p_1 - p_t) * (         # noqa
                    a**4 * pi / (8 * eta) * p_m  + a**3 * 2 * pi / 3 * v_m * (  # noqa
                        (1 + 2 * a * nu * p_m / eta) / (                        # noqa
                        1 + 2.48 * a * nu * p_m / eta                           # noqa
                        )                                                       # noqa
                    )                                                           # noqa
                )                                                               # noqa
            )                                                                   # noqa
        )                                                                       # noqa
        # fmt: on
        n_dot = N_dot / AVOGADROS_CONSTANT
        return n_dot

    def gas_flux_calibration(
        self,
        measurement,
        mol,
        mass,
        tspan=None,
        tspan_bg=None,
        ax=None,
    ):
        """
        Args:
            measurement (MSMeasurement): The measurement with the calibration data
            mol (str): The name of the molecule to calibrate
            mass (str): The mass to calibrate at
            tspan (iter): The timespan to average the signal over. Defaults to all
            tspan_bg (iter): Optional timespan at which the signal is at its background.
            ax (matplotlib axis): the axis on which to indicate what signal is used
                with a thicker line. Defaults to none

        Returns MSCalResult: a calibration result containing the sensitivity factor for
            mol at mass
        """
        t, S = measurement.grab_signal(mass, tspan=tspan, t_bg=tspan_bg)
        if ax:
            ax.plot(t, S, color=STANDARD_COLORS[mass], linewidth=5)

        n_dot = self.calc_n_dot_0(gas=mol)
        F = np.mean(S) / n_dot
        return MSCalResult(
            name=f"{mol}_{mass}",
            mol=mol,
            mass=mass,
            cal_type="gas_flux_calibration",
            F=F,
        )


class MSSpectrum(Spectrum):
    """Nothing to add to normal spectrum yet.
    TODO: Methods for co-plotting ref spectra from a database
    """

    pass
