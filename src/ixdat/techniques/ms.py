"""Module for representation and analysis of MS measurements"""

import re
import numpy as np

from ..measurement_base import Measurement
from ..spectra import Spectrum, SpectrumSeries, SpectroMeasurement
from ..plotters import MSPlotter, MSSpectroPlotter
from ..plotters.ms_plotter import STANDARD_COLORS
from ..exporters import MSExporter, MSSpectroExporter
from ..tools import deprecate
from ..plugins import plugins
from ..calculators.ms_calculators import (
    MSCalResult,
    MSBackgroundSet,
    MSCalibration,
    MSConstantBackground,
)

# and, for back-compatibility until 0.3.1:
from ..calculators.ms_calculators import MSInlet  # noqa: F401


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    default_plotter = MSPlotter
    default_exporter = MSExporter
    background_calculator_types = [MSBackgroundSet]

    def __init__(self, name, **kwargs):
        tspan_bg = kwargs.pop("tspan_bg", None)
        super().__init__(name, **kwargs)
        self.tspan_bg = tspan_bg
        self._siq_quantifier = None  # Used with external quantification package

    @property
    def ms_calibration(self):
        return self.calculators["MS calibration"]

    @property
    def siq_calculator(self):
        return self.calculators["siq calculator"]

    @property
    def signal_bgs(self):
        return self.calculators["MS background"]

    def set_bg(self, tspan=None, tspan_bg=None, mass_list=None):
        """Set background values for mass_list to the average signal during tspan_bg."""
        tspan = tspan or tspan_bg
        background = MSBackgroundSet.from_measurement_point(
            measurement=self, tspan=tspan, mass_list=mass_list
        )
        self.add_calculator(background)

    def reset_bg(self, mass_list=None):
        """Reset background values for all masses or the masses in mass_list"""
        if mass_list is None:
            new_calculator_list = [
                cal
                for cal in self.calculator_list
                if type(cal) not in self.background_calculator_types
            ]
            self._calculator_list = new_calculator_list
            self.consolidate_calculators()
        else:
            self.add_calculator(
                MSBackgroundSet(
                    bg_list=[MSConstantBackground(mass=mass, bg=0) for mass in mass_list]
                )
            )

    def grab_signal(self, *args, **kwargs):
        """Alias for grab()"""
        return self.grab(*args, **kwargs)

    @deprecate(
        "0.1", "Use `remove_background` instead.", "0.3.1", kwarg_name="removebackground"
    )
    def grab_flux(
        self,
        mol,
        tspan=None,
        tspan_bg=None,
        remove_background=True,
        removebackground=None,
        include_endpoints=False,
    ):
        """Return the flux of mol (calibrated signal) in [mol/s]

        .. note::
          * With native ixdat quantification (use_siq=False),
            `grab_flux(mol, ...)` is identical to `grab(f"n_dot_{mol}", ...)` with
            remove_background=True by default. An MSCalibration does the maths.
          * With an external quantification package (use_siq=True), the maths are done
            here with the help of self.quantifier

        Args:
            mol (str or MSCalResult): Name of the molecule or a ms_calibration thereof
            tspan (list): Timespan for which the signal is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if available
                Defaults to True.
            removebackground (bool): DEPRECATED. Use `remove_background`.
            include_endpoints (bool): Whether to interpolate for tspan[0] and tspan[-1]
        """
        if removebackground is not None:
            remove_background = removebackground

        if isinstance(mol, MSCalResult):
            t, signal = self.grab(
                mol.mass,
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                include_endpoints=include_endpoints,
            )
            return t, signal / mol.F
        return self.grab(
            # grab() invokes __getitem__, which invokes the `Calibration`. Specifically,
            # `MSCalibration.calculate_series()` interprets item names starting with
            # "n_" as molecule fluxes, and checks itself for a sensitivity factor.
            f"n_dot_{mol}",
            tspan=tspan,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
            include_endpoints=include_endpoints,
        )

    @deprecate(
        "0.1", "Use `remove_background` instead.", "0.3.1", kwarg_name="removebackground"
    )
    def grab_flux_for_t(
        self,
        mol,
        t,
        tspan_bg=None,
        remove_background=False,
        removebackground=None,
    ):
        """Return the flux of mol (calibrated signal) in [mol/s] for a given time vec

        Args:
            mol (str): Name of the molecule.
            t (np.array): The time vector along which to give the flux
            tspan_bg (tspan): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if available
            removebackground (bool): DEPRECATED. Use `remove_background`.
        """
        if removebackground is not None:
            remove_background = removebackground
        t_0, y_0 = self.grab_flux(
            mol,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
        )
        y = np.interp(t, t_0, y_0)
        return y

    def get_flux_series(self, mol):
        """Return a ValueSeries with the calibrated flux of mol"""
        return self[f"n_dot_{mol}"]

    def integrate_signal(self, mass, tspan, tspan_bg, ax=None):
        """Integrate a ms signal with background subtraction and evt. plotting

        TODO:
            Should this, like grab_signal does now, have the option of using a
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

    def integrate_flux(self, mol, tspan, tspan_bg, ax=None):

        """Integrate a calibrated ms signal with background subtraction and evt.
        plotting (copy of integrate_signal method)

        TODO:
          Should this, like grab_signal does now, have the option of using a
          background saved in the object rather than calculating a new one?

        TODO:
            Ensure fill_between considers the non-standard unit in the figure

        Args:
            mol (str): The molecule name for which to integrate the signal
            tspan (tspan): The timespan over which to integrate
            tspan_bg (tspan): Timespan at which the signal is at its background value
            ax (Axis): axis to plot on. Defaults to None
        """
        t, S = self.grab_flux(mol, tspan=tspan, include_endpoints=True)
        if tspan_bg:
            t_bg, S_bg_0 = self.grab_flux(mol, tspan=tspan_bg, include_endpoints=True)
            S_bg = np.mean(S_bg_0) * np.ones(t.shape)
        else:
            S_bg = np.zeros(t.shape)
        if ax:
            if ax == "new":
                fig, ax = self.plotter.new_ax()
            ax.fill_between(t, S_bg, S, color=STANDARD_COLORS[mol], alpha=0.2)
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
        raise TypeError(f"{self!r} does not recognize '{item}' as a mass.")

    # --- METHODS WHICH HAVE BEEN MOVED TO `Calculator` CLASSES ---- #

    @deprecate(
        "0.2.13",
        "Use `MSCalibration.gas_flux_calibration` instead.",
        "0.3.1",
    )
    def gas_flux_calibration(self, *args, **kwargs):
        return MSCalibration.gas_flux_calibration(measurement=self, *args, **kwargs)

    @deprecate(
        "0.2.13",
        "Use `MSCalibration.gas_flux_calibration_curve` instead.",
        "0.3.1",
    )
    def gas_flux_calibration_curve(self, *args, **kwargs):
        return MSCalibration.gas_flux_calibration_curve(
            measurement=self, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.gas_flux_calibration` instead.",
        "0.3.1",
    )
    def siq_gas_flux_calibration(self, mol, mass, tspan, chip=None):
        return plugins.siq.Calculator.gas_flux_calibration(
            measurement=self, mol=mol, mass=mass, tspan=tspan, chip=None
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.gas_flux_calibration_curve` instead.",
        "0.3.1",
    )
    def siq_gas_flux_calibration_curve(
        self,
        *args,
        **kwargs,
    ):
        return plugins.siq.Calculator.gas_flux_calibration_curve(
            measurement=self, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.multicomp_gas_flux_calibration` instead.",
        "0.3.1",
    )
    def siq_multicomp_gas_flux_calibration(self, *args, **kwargs):
        return plugins.siq.Calculator.siq_multicomp_gas_flux_calibration(
            measurement=self, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `siq_calculator.grab_fluxes` instead.",
        "0.3.1",
    )
    def grab_siq_fluxes(
        self,
        tspan=None,
        tspan_bg=None,
        remove_background=False,
        include_endpoints=False,
    ):
        return self.siq_calculator.grab_fluxes(
            tspan=None, tspan_bg=None, remove_background=False, include_endpoints=False
        )


class MSSpectrum(Spectrum):
    """Nothing to add to normal Spectrum yet.
    TODO: Methods for co-plotting ref spectra from a database
    """

    pass


class MSSpectrumSeries(SpectrumSeries):
    """Nothing to add to normal SpectrumSeries yet."""

    pass


class MSSpectroMeasurement(MSMeasurement, SpectroMeasurement):
    extra_column_attrs = SpectroMeasurement.extra_column_attrs
    default_plotter = MSSpectroPlotter
    default_exporter = MSSpectroExporter

    # FIXME: https://github.com/ixdat/ixdat/pull/166#discussion_r1486023530
