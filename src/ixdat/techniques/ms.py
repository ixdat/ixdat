"""Module for representation and analysis of MS measurements"""

import re
import numpy as np
import json  # FIXME: This is for MSCalibration.export, but shouldn't have to be here.

from ..measurements import Measurement, Calibration
from ..spectra import Spectrum
from ..plotters.ms_plotter import MSPlotter, STANDARD_COLORS
from ..exceptions import QuantificationError
from ..constants import (
    AVOGADROS_CONSTANT,
    BOLTZMAN_CONSTANT,
    STANDARD_TEMPERATURE,
    STANDARD_PRESSURE,
    DYNAMIC_VISCOSITIES,
    MOLECULAR_DIAMETERS,
    MOLAR_MASSES,
)
from ..data_series import ValueSeries
from ..db import Saveable


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    extra_column_attrs = {"ms_measurement": ("tspan_bg",)}
    default_plotter = MSPlotter

    def __init__(self, name, **kwargs):
        tspan_bg = kwargs.pop("tspan_bg", None)
        super().__init__(name, **kwargs)
        self.tspan_bg = tspan_bg

    @property
    def ms_calibration(self):
        ms_cal_list = []
        tspan_bg = None
        signal_bgs = {}
        for cal in self.calibration_list:
            ms_cal_list = ms_cal_list + getattr(cal, "ms_cal_list", [])
            for mass, bg in getattr(cal, "signal_bgs", {}).items():
                if mass not in signal_bgs:
                    signal_bgs[mass] = bg
            tspan_bg = tspan_bg or getattr(cal, "tspan_bg", None)
        return MSCalibration(ms_cal_results=ms_cal_list, signal_bgs=signal_bgs)

    @property
    def signal_bgs(self):
        return self.ms_calibration.signal_bgs

    def set_bg(self, tspan_bg=None, mass_list=None):
        """Set background values for mass_list to the average signal during tspan_bg."""
        mass_list = mass_list or self.mass_list
        tspan_bg = tspan_bg or self.tspan_bg
        signal_bgs = {}
        for mass in mass_list:
            t, v = self.grab(mass, tspan_bg)
            signal_bgs[mass] = np.mean(v)
        self.add_calibration(MSCalibration(signal_bgs=signal_bgs))

    def reset_bg(self, mass_list=None):
        """Reset background values for the masses in mass_list"""
        mass_list = mass_list or self.mass_list
        for mass in mass_list:
            if mass in self.signal_bgs:
                del self.signal_bgs[mass]

    def grab(
        self,
        item,
        tspan=None,
        tspan_bg=None,
        include_endpoints=False,
        remove_background=False,
    ):
        """Returns t, S where S is raw signal in [A] for a given signal name (ie mass)

        Args:
            item (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if
                available. This is special to MSMeasurement.
                Defaults to False, but in grab_flux it defaults to True.
            include_endpoints (bool): Whether to ensure tspan[0] and tspan[-1] are in t
        """
        time, value = super().grab(
            item, tspan=tspan, include_endpoints=include_endpoints
        )

        if tspan_bg:
            _, bg = self.grab(item, tspan=tspan_bg)
            return time, value - np.average(bg)
        elif remove_background:
            if item in self.signal_bgs:
                return time, value - self.signal_bgs[item]
            elif self.tspan_bg:
                _, bg = self.grab(item, tspan=self.tspan_bg)
                return time, value - np.average(bg)
        return time, value

    def grab_for_t(self, item, t, tspan_bg=None, remove_background=False):
        """Return a numpy array with the value of item interpolated to time t

        Args:
            item (str): The name of the value to grab
            t (np array): The time vector to grab the value for
            tspan_bg (iterable): Optional. A timespan defining when `item` is at its
                baseline level. The average value of `item` in this interval will be
                subtracted from what is returned.
            remove_background (bool): Whether to remove a pre-set background if
                available. This is special to MSMeasurement.
                Defaults to False, but in grab_flux it defaults to True.
        """
        t_0, v_0 = self.grab(
            item, tspan_bg=tspan_bg, remove_background=remove_background
        )
        v = np.interp(t, t_0, v_0)
        return v

    def grab_signal(self, *args, **kwargs):
        """Alias for grab()"""
        return self.grab(*args, **kwargs)

    def grab_flux(
        self,
        mol,
        tspan=None,
        tspan_bg=None,
        remove_background=True,
        include_endpoints=False,
    ):
        """Return the flux of mol (calibrated signal) in [mol/s]

        Note:
        `grab_flux(mol, ...)` is identical to `grab(f"n_dot_{mol}", ...)` with
        removebackround=True by default. An MSCalibration does the maths.

        Args:
            mol (str or MSCalResult): Name of the molecule or a ms_calibration thereof
            tspan (list): Timespan for which the signal is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if available
                Defaults to True.
        """
        return self.grab(
            # grab() invokes __getitem__, which invokes the `Calibration`. Specifically,
            # `MSCalibration.calibrate_series()` interprets item names starting with
            # "n_" as molecule fluxes, and checks itself for a sensitivity factor.
            f"n_dot_{mol}",
            tspan=tspan,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
            include_endpoints=include_endpoints,
        )

    def grab_flux_for_t(
        self,
        mol,
        t,
        tspan_bg=None,
        remove_background=False,
        include_endpoints=False,
    ):
        """Return the flux of mol (calibrated signal) in [mol/s] for a given time vec

        Args:
            mol (str): Name of the molecule.
            t (np.array): The time vector along which to give the flux
            tspan_bg (tspan): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if available
        """
        t_0, y_0 = self.grab_flux(
            mol,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
            include_endpoints=include_endpoints,
        )
        y = np.interp(t, t_0, y_0)
        return y

    def get_flux_series(self, mol):
        """Return a ValueSeries with the calibrated flux of mol"""
        return self[f"n_dot_{mol}"]

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
            t_bg, S_bg_0 = self.grab_signal(mass, tspan=tspan_bg, include_endpoints=True)
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
        if item in self.reverse_aliases and self.is_mass(self.reverse_aliases[item][0]):
            return True
        return False

    def as_mass(self, item):
        if re.search("^M[0-9]+$", item):
            return item
        new_item = self.reverse_aliases[item][0]
        if self.is_mass(new_item):
            return self.as_mass(new_item)
        raise TypeError(f"{self} does not recognize '{item}' as a mass.")


class MSCalResult(Saveable):
    """A class for a mass spec ms_calibration result.

    FIXME: I think that something inheriting directly from Saveable does not belong in
        a technique module.
    """

    table_name = "ms_cal_results"
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


class MSCalibration(Calibration):
    """Class for mass spec calibrations. TODO: replace with powerful external package"""

    extra_linkers = {"ms_calibration_results": ("ms_cal_results", "ms_cal_result_ids")}
    # FIXME: signal_bgs are not saved at present. Should they be a separate table
    #   of Saveable objects like ms_cal_results or should they be a single json value?
    child_attrs = [
        "ms_cal_results",
    ]

    def __init__(
        self,
        name=None,
        date=None,
        tstamp=None,  # FIXME: No need to have both a date and a tstamp?
        setup=None,
        ms_cal_results=None,
        signal_bgs=None,
        technique="MS",
        measurement=None,
    ):
        """
        Args:
            name (str): Name of the ms_calibration
            date (str): Date of the ms_calibration
            setup (str): Name of the setup where the ms_calibration is made
            ms_cal_results (list of MSCalResult): The mass spec calibrations
            measurement (MSMeasurement): The measurement
        """
        super().__init__(
            name=name or f"EC-MS ms_calibration for {setup} on {date}",
            technique=technique,
            tstamp=tstamp,
            measurement=measurement,
        )
        self.date = date
        self.setup = setup
        self.ms_cal_results = ms_cal_results or []
        self.signal_bgs = signal_bgs or {}

    @property
    def ms_cal_result_ids(self):
        return [cal.id for cal in self.ms_cal_results]

    @property
    def mol_list(self):
        return list({cal.mol for cal in self.ms_cal_results})

    @property
    def mass_list(self):
        return list({cal.mass for cal in self.ms_cal_results})

    @property
    def name_list(self):
        return list({cal.name for cal in self.ms_cal_results})

    def __contains__(self, mol):
        return mol in self.mol_list or mol in self.name_list

    def __iter__(self):
        yield from self.ms_cal_results

    def calibrate_series(self, key, measurement=None):
        """Return a calibrated series for `key` if possible.

        If key starts with "n_", it is interpreted as a molecule flux. This method then
        searches the calibration for a sensitivity factor for that molecule uses it to
        divide the relevant mass signal from the measurement. Example acceptable keys:
        "n_H2", "n_dot_H2".
        If the key does not start with "n_", or the calibration can't find a relevant
        sensitivity factor and mass signal, this method returns None.
        """
        measurement = measurement or self.measurement
        if key.startswith("n_"):  # it's a flux!
            mol = key.split("_")[-1]
            try:
                mass, F = self.get_mass_and_F(mol)
            except QuantificationError:
                # Calibrations just return None when they can't get what's requested.
                return
            signal_series = measurement[mass]
            y = signal_series.data
            if mass in measurement.signal_bgs:
                # FIXME: How to make this optional to user of MSMeasuremt.grab()?
                y = y - measurement.signal_bgs[mass]
            n_dot = y / F
            return ValueSeries(
                name=f"n_dot_{mol}",
                unit_name="mol/s",
                data=n_dot,
                tseries=signal_series.tseries,
            )

    def get_mass_and_F(self, mol):
        """Return the mass and sensitivity factor to use for simple quant. of mol"""
        cal_list_for_mol = [cal for cal in self if cal.mol == mol or cal.name == mol]
        Fs = [cal.F for cal in cal_list_for_mol]
        if not Fs:
            raise QuantificationError(f"{self} has no sensitivity factor for {mol}")
        index = np.argmax(np.array(Fs))

        the_good_cal = cal_list_for_mol[index]
        return the_good_cal.mass, the_good_cal.F

    def get_F(self, mol, mass):
        """Return the sensitivity factor for mol at mass"""
        cal_list_for_mol_at_mass = [
            cal
            for cal in self
            if (cal.mol == mol or cal.name == mol) and cal.mass == mass
        ]
        F_list = [cal.F for cal in cal_list_for_mol_at_mass]
        if not F_list:
            raise QuantificationError(
                f"{self} has no sensitivity factor for {mol} at {mass}"
            )
        return np.mean(np.array(F_list))

    def scaled_to(self, ms_cal_result):
        """Return a new ms_calibration w scaled sensitivity factors to match one given"""
        F_0 = self.get_F(ms_cal_result.mol, ms_cal_result.mass)
        scale_factor = ms_cal_result.F / F_0
        calibration_as_dict = self.as_dict()
        new_cal_list = []
        for cal in self.ms_cal_results:
            cal = MSCalResult(
                name=cal.name,
                mass=cal.mass,
                mol=cal.mol,
                F=cal.F * scale_factor,
                cal_type=cal.cal_type + " scaled",
            )
            new_cal_list.append(cal)
        calibration_as_dict["ms_cal_results"] = new_cal_list
        del calibration_as_dict["ms_cal_result_ids"]
        # ^ FIXME: ms_cal_result_ids via MemoryBackend
        calibration_as_dict["name"] = calibration_as_dict["name"] + " scaled"
        return self.__class__.from_dict(calibration_as_dict)

    @classmethod
    def read(cls, path_to_file):
        """Read an MSCalibration from a json-formatted text file"""
        with open(path_to_file) as f:
            obj_as_dict = json.load(f)
        # put the MSCalResults (exported as dicts) into objects:
        obj_as_dict["ms_cal_results"] = [
            MSCalResult.from_dict(ms_cal_as_dict)
            for ms_cal_as_dict in obj_as_dict["ms_cal_results"]
        ]
        return cls.from_dict(obj_as_dict)

    def export(self, path_to_file=None):
        """Export an ECMSCalibration as a json-formatted text file"""
        path_to_file = path_to_file or (self.name + ".ix")
        self_as_dict = self.as_dict()
        # replace the ms_cal_result ids with the dictionaries of the results themselves:
        del self_as_dict["ms_cal_result_ids"]
        self_as_dict["ms_cal_results"] = [cal.as_dict() for cal in self.ms_cal_results]
        with open(path_to_file, "w") as f:
            json.dump(self_as_dict, f, indent=4)


class MSInlet:
    """A class for describing the inlet to the mass spec

    Every MSInlet describes the rate and composition of the gas entering a mass
    spectrometer. The default is a Spectro Inlets EC-MS chip.
    TODO: Replace with powerful external package.
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
        """Create an MSInlet object given its properties.

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
        self.l_cap_eff = {}
        self.w_cap = w_cap
        self.h_cap = h_cap
        self.p = p
        self.T = T
        self.gas = gas  # TODO: Gas mixture class. This must be a pure gas now.

    def calc_l_cap_eff(
        self, n_dot_measured, gas=None, w_cap=None, h_cap=None, T=None, p=None
    ):
        """Calculate gas specific effective length of the capillary in [m]
        and add {gas:value} to l_cap_eff (dict)

        Args:
            w_cap (float): Capillary width [m], defaults to self.w_cap
            h_cap (float): Capillary height [m], defaults to self.h_cap
            n_dot_measured (float): Measured flux of gas [mol/s]
            gas (dict or str): The gas in the chip, defaults to self.gas
            T (float): Temperature [K], if to be updated
            p (float): Pressure [Pa], if to be updated
        Returns:
            float: Gas specific effective length in [m]
        """

        n_dot_predicted = self.calc_n_dot_0(gas=gas, w_cap=w_cap, h_cap=h_cap, T=T, p=p)

        l_cap_gas_specific_eff = self.l_cap * n_dot_predicted / n_dot_measured
        self.l_cap_eff[
            gas
        ] = l_cap_gas_specific_eff  # add effective l_cap for specific gas

        return l_cap_gas_specific_eff

    def update_l_cap(self, gases=[]):
        """Update self.l_cap from average of values in dict l_cap_eff

        Args:
            gases (list): List of gases to average l_cap, default all
        Returns:
            float: Averaged effective capilllary length in [m]
        """
        if self.l_cap_eff and not gases:
            self.l_cap = np.mean(list(self.l_cap_eff.values()))
        elif self.l_cap_eff and gases:
            _l_cap = 0
            for gas in gases:
                _l_cap += self.l_cap_eff[gas]
            self.l_cap = _l_cap / len(gases)

        return self.l_cap

    def calc_n_dot_0(self, gas=None, w_cap=None, h_cap=None, l_cap=None, T=None, p=None):
        """Calculate the total molecular flux through the capillary in [s^-1]

        Uses Equation 4.10 of Trimarco, 2017. "Real-time detection of sub-monolayer
        desorption phenomena during electrochemical reactions: Instrument development
        and applications." PhD Thesis, Technical University of Denmark.

        Args:
            w_cap (float): Capillary width [m], defaults to self.w_cap
            h_cap (float): Capillary height [m], defaults to self.h_cap
            l_cap (float): Capillary length [m], defaults to self.l_cap
            gas (dict or str): The gas in the chip, defaults to self.gas
            T (float): Temperature [K], if to be updated
            p (float): Pressure [Pa], if to be updated
        Returns:
            float: The total molecular flux in [s^-1] through the capillary
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
        p_t = BOLTZMAN_CONSTANT * T / (2**0.5 * pi * s**2 * lambda_)
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
        carrier_mol=None,
        mol_conc_ppm=None,
    ):
        """
        Args:
            measurement (MSMeasurement): The measurement with the ms_calibration data
            mol (str): The name of the molecule to calibrate
            mass (str): The mass to calibrate at
            tspan (iter): The timespan to average the signal over. Defaults to all
            tspan_bg (iter): Optional timespan at which the signal is at its background.
            ax (matplotlib axis): The axis on which to indicate what signal is used
                with a thicker line. Defaults to none
            carrier_mol (str): The name of the molecule of the carrier gas if
                a dilute analyte is used. Calibration assumes total flux of the
                capillary is the same as the flux of pure carrier gas. Defaults
                to None.
            mol_conc_ppm (float): Concentration of the dilute analyte in the carrier gas
                in ppm. Defaults to None.

        Returns MSCalResult: a ms_calibration result containing the sensitivity factor
            for mol at mass
        """
        t, S = measurement.grab_signal(mass, tspan=tspan, tspan_bg=tspan_bg)
        if ax:
            ax.plot(t, S, color=STANDARD_COLORS[mass], linewidth=5)
        if carrier_mol:
            if mol_conc_ppm:
                cal_type = "carrier_gas_flux_calibration"
            else:
                raise QuantificationError(
                    "Cannot use carrier gas calibration without analyte"
                    " concentration. mol_conc_ppm is missing."
                )
        elif mol_conc_ppm:
            raise QuantificationError(
                "Cannot use carrier gas calibration without carrier"
                " gas definition. carrier_mol is missing."
            )
        else:
            cal_type = "gas_flux_calibration"
            mol_conc_ppm = 10**6
            carrier_mol = mol

        n_dot = self.calc_n_dot_0(gas=carrier_mol) * mol_conc_ppm / 10**6
        F = np.mean(S) / n_dot
        return MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type=cal_type,
            F=F,
        )


class MSSpectrum(Spectrum):
    """Nothing to add to normal spectrum yet.
    TODO: Methods for co-plotting ref spectra from a database
    """

    pass
