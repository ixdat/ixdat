"""Module for representation and analysis of MS measurements"""

import re
import numpy as np
import json  # FIXME: This is for MSCalibration.export, but shouldn't have to be here.
import warnings

from ..measurements import Measurement, Calculator
from ..spectra import Spectrum, SpectrumSeries, SpectroMeasurement
from ..plotters import MSPlotter, MSSpectroPlotter
from ..plotters.ms_plotter import STANDARD_COLORS
from ..exporters import MSExporter, MSSpectroExporter
from ..exceptions import QuantificationError, SeriesNotFoundError
from ..constants import (
    AVOGADRO_CONSTANT,
    BOLTZMANN_CONSTANT,
    STANDARD_TEMPERATURE,
    STANDARD_PRESSURE,
    DYNAMIC_VISCOSITIES,
    MOLECULAR_DIAMETERS,
    MOLAR_MASSES,
)
from ..data_series import ValueSeries
from ..db import Saveable
from ..tools import deprecate
from ..config import plugins


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    # FIXME: tspan_bg should be column of a Calculator
    #   (see https://github.com/ixdat/ixdat/issues/164)
    extra_column_attrs = {"ms_measurement": ("tspan_bg",)}
    default_plotter = MSPlotter
    default_exporter = MSExporter

    def __init__(self, name, **kwargs):
        tspan_bg = kwargs.pop("tspan_bg", None)
        super().__init__(name, **kwargs)
        self.tspan_bg = tspan_bg
        self._siq_quantifier = None  # Used with external quantification package

    @property
    def ms_calibration(self):
        return self.calculators["MS calibration"]

    @property
    def background_calculator_types(self):
        # FIXME: This should be a static class attribute, but the problem is that
        #    it is defined below.
        return [MSBackgroundSet]

    @property
    def siq_calculator(self):
        if not plugins.use_si_quant:
            raise QuantificationError("siq is not activated.")
        for calculator in self.calculator_list:
            if isinstance(calculator, plugins.siq.Calculator):
                return calculator
        raise QuantificationError("measurement has no siq Calculator.")

    @property
    def signal_bgs(self):
        return self.calculators["MS background"]

    def set_bg(self, tspan=None, tspan_bg=None, mass_list=None):
        """Set background values for mass_list to the average signal during tspan_bg."""
        tspan = tspan or tspan_bg
        background = MSBackgroundSet.from_measurement_point(
            measurement=self, tspan=tspan_bg, mass_list=mass_list
        )
        self.add_calculator(background)

    def reset_bg(self, mass_list=None):
        """Reset background values for the masses in mass_list"""
        raise NotImplementedError

    def grab_signal(self, *args, **kwargs):
        """Alias for grab()"""
        return self.grab(*args, **kwargs)

    @deprecate(
        "0.1", "Use `remove_background` instead.", "0.3", kwarg_name="removebackground"
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

        Note:
        - With native ixdat quantification (use_siq=False),
          `grab_flux(mol, ...)` is identical to `grab(f"n_dot_{mol}", ...)` with
          remove_background=True by default. An MSCalibration does the maths.
        - With an external quantification package (use_siq=True), the maths are done
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
        "0.1", "Use `remove_background` instead.", "0.3", kwarg_name="removebackground"
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
        raise TypeError(f"{self!r} does not recognize '{item}' as a mass.")

    # --- METHODS WHICH HAVE BEEN MOVED TO `Calculator` CLASSES ---- #

    @deprecate(
        "0.2.13",
        "Use `MSCalibration.gas_flux_calibration` instead.",
        "0.3",
    )
    def gas_flux_calibration(self, *args, **kwargs):
        return MSCalibration.gas_flux_calibration(measurement=self, *args, **kwargs)

    @deprecate(
        "0.2.13",
        "Use `MSCalibration.gas_flux_calibration_curve` instead.",
        "0.3",
    )
    def gas_flux_calibration_curve(self, *args, **kwargs):
        return MSCalibration.gas_flux_calibration_curve(
            measurement=self, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.gas_flux_calibration` instead.",
        "0.3",
    )
    def siq_gas_flux_calibration(self, mol, mass, tspan, chip=None):
        return plugins.siq.Calculator.gas_flux_calibration(
            measurement=self, mol=mol, mass=mass, tspan=tspan, chip=None
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.gas_flux_calibration_curve` instead.",
        "0.3",
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
        "0.3",
    )
    def siq_multicomp_gas_flux_calibration(self, *args, **kwargs):
        return plugins.siq.Calculator.siq_multicomp_gas_flux_calibration(
            measurement=self, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `siq_calculator.grab_fluxes` instead.",
        "0.3",
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


class MSConstantBackground(Saveable):
    """A small Saveable wrapper around a constant MS background

    FIXME: There's nothing about this data structure unique to MS. Should there be a
    general ConstantBackground class for all Measurement types?
    """

    table_name = "ms_constant_backgrounds"
    column_attrs = {"mass", "bg"}

    def __init__(
        self,
        mass,
        bg,
        name=None,
    ):
        self.mass = mass
        self.bg = bg
        self.name = name or f"{mass} bg"

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name={self.name}, mol={self.mol}, bg={self.bg})"
        )

    def get_bg_for_t(self, t):
        return self.bg


class MSBackgroundSet(Calculator):
    calculator_name = "MS background"

    extra_linkers = {
        "ms_background_constants": ("ms_constant_backgrounds", "ms_constant_bg_ids")
        # Additional types of backgrounds would go here
    }
    child_attrs = [
        "constant_bg_list",
        # FIXME: This repeats info in extra_linkers. Should be possible to combine. #75.
    ]

    def __init__(
        self,
        name=None,
        date=None,
        bg_list=None,
        technique="MS",
        tstamp=None,
        measurement=None,
    ):
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.constant_bg_list = []
        for bg in bg_list:
            if type(bg) is MSConstantBackground:
                self.constant_bg_list.append(bg)
            # elif type(bg) is OtherBackgroundClassWhenItIsReady:
            #   self.other_bacground_list.append(bg)
            else:
                warnings.warn(
                    f"Initiating an MSBackgoundSet with a background list including {bg}"
                    ", which is not a recognized MS background type. Skipping that one."
                )

    @classmethod
    def from_measurement_point(cls, measurement, tspan, mass_list=None):
        mass_list = mass_list or measurement.mass_list
        bg_list = []
        for mass in mass_list:
            _, S = measurement.grab(mass, remove_background=False, tspan=tspan)
            bg_list.append(MSConstantBackground(mass=mass, bg=np.mean(S)))
        return cls(bg_list=bg_list, measurement=measurement)

    @property
    def ms_constant_bg_ids(self):
        return [cal.id for cal in self.constant_bg_list]

    @property
    def bg_list(self):
        return self.constant_bg_list  # + self.other_bg_lists_when_they_are_here

    @property
    def mass_list(self):
        return list({bg.mass for bg in self.bg_list})

    @property
    def available_series_names(self):
        return set([f"{mass}" for mass in self.mass_list])

    def calculate_series(self, mass, measurement=None):
        """Return a calibrated series for `key` if possible.

        If key starts with "n_", it is interpreted as a molecule flux. This method then
        searches the calibration for a sensitivity factor for that molecule uses it to
        divide the relevant mass signal from the measurement. Example acceptable keys:
        "n_H2", "n_dot_H2".
        If the key does not start with "n_", or the calibration can't find a relevant
        sensitivity factor and mass signal, this method returns None.
        """
        measurement = measurement or self.measurement

        raw_series = measurement[mass + "-raw"]
        bg_object = next(obj for obj in self.bg_list if obj.mass == mass)
        bg = bg_object.get_bg_for_t(raw_series.t)

        return ValueSeries(
            name=mass + "_bg_corrected",
            unit_name="mol/s",
            data=raw_series.data - bg,
            tseries=raw_series.tseries,
        )


class MSCalResult(Saveable):
    """A small saveable wrapper around an MS calibration."""

    table_name = "ms_cal_results"
    column_attrs = {"name", "mol", "mass", "cal_type", "F"}

    def __init__(
        self,
        mol,
        mass,
        F,
        name=None,
        cal_type=None,
    ):
        super().__init__()
        self.name = name or f"{mol}@{mass}"
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

    @classmethod
    def from_siq(cls, siq_cal_point):
        return cls(
            name=siq_cal_point.mol + "@" + siq_cal_point.mass,
            mol=siq_cal_point.mol,
            mass=siq_cal_point.mass,
            cal_type=siq_cal_point.F_type,
            F=siq_cal_point.F,
        )

    def to_siq(self):
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSCalPoint.to_siq` only works when using "
                "`spectro_inlets_quantification` "
                "(`ixdat.options.activate_siq()`). "
                "For native ixdat MS quantification, use `MSMeasurement.calibrate`"
            )
        return plugins.siq.CalPoint(
            mol=self.mol,
            mass=self.mass,
            F=self.F,
            F_type=self.cal_type,
        )


class MSCalibration(Calculator):
    """Class for mass spec calibrations.

    This is a simple version, just implementing sensitivity facotrs.
    Consider using the more powerful spectro_inlets_quantification package instead.
    """

    calculator_name = "MS calibration"

    extra_linkers = {"ms_calibration_results": ("ms_cal_results", "ms_cal_result_ids")}
    child_attrs = [
        "ms_cal_results",
    ]

    def __init__(
        self,
        mol=None,
        mass=None,
        F=None,
        name=None,
        date=None,
        tstamp=None,  # FIXME: No need to have both a date and a tstamp?
        setup=None,
        ms_cal_results=None,
        technique="MS",
        measurement=None,
    ):
        """
        Args:
            name (str): Name of the ms_calibration
            date (str): Date of the ms_calibration
            setup (str): Name of the setup where the ms_calibration is made
            ms_cal_results (list of MSCalResult): The mass spec calibrations
            measurement (MSMeasurement): The measurement used to make the calibration
        """
        if plugins.use_siq:
            warnings.warn(
                "spectro_inlets_quantification is active but you are making a native "
                "ixdat MSCalibration. It's probably best not to mix the two!"
            )
        super().__init__(
            name=name,
            technique=technique,
            tstamp=tstamp,
            measurement=measurement,
        )
        if mol and mass and (F is not None):
            if ms_cal_results:
                raise QuantificationError(
                    "Attempt to initiate MSCalibration with both"
                    "mol, mass, and F, *and* a ms_cal_results list. Choose one."
                )
            ms_cal_results = [MSCalResult(mol, mass, F)]
        self.date = date
        self.setup = setup
        self.ms_cal_results = ms_cal_results or []

    @classmethod
    @deprecate(
        "0.2.6",
        "Use `inlet` instead. Or consider using `siq_gas_flux_calibration` "
        "with the `spectro_inlets_quantification` package.",
        "0.3",
        kwarg_name="chip",
    )
    def gas_flux_calibration(
        cls,
        measurement,
        mol,
        mass,
        inlet=None,
        chip=None,
        tspan=None,
        tspan_bg=None,
        ax=None,
        carrier_mol=None,
        mol_conc_ppm=None,
    ):
        """
        Fit mol's sensitivity at mass based on one period with steady gas composition.

        Args:
            mol (str): The name of the molecule to calibrate
            mass (str): The mass to calibrate at
            inlet (MSInlet): An object with a `calc_n_dot_0` method for total flux to
                the vacuum chamber containing the mass spectrometer.
            chip (MSInlet): DEPRECATED. Old name for `inlet`.
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

        Returns MSCalResult: a MS calibration result containing the sensitivity factor
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
        n_dot = inlet.calc_n_dot_0(gas=carrier_mol) * mol_conc_ppm / 10**6
        F = np.mean(S) / n_dot
        cal_result = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type=cal_type,
            F=F,
        )
        return cls(ms_cal_results=[cal_result], measurement=measurement)

    @classmethod
    def gas_flux_calibration_curve(
        cls,
        measurement,
        mol,
        mass,
        inlet=None,
        chip=None,
        tspan_list=None,
        carrier_mol=None,
        mol_conc_ppm=None,
        p_inlet=None,
        tspan_bg=None,
        ax="new",
        axis_measurement=None,
        remove_bg_on_axis_measurement=True,
        return_ax=False,
    ):
        """Fit mol's sensitivity at mass from 2+ periods of steady gas composition.

        Args:
            inlet (MSInlet): An object with a `calc_n_dot_0` method for total flux to
                the vacuum chamber containing the mass spectrometer.
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            inlet (MSInlet): An object with a `calc_n_dot_0` method for total flux to
                the vacuum chamber containing the mass spectrometer.
            chip (MSInlet): DEPRECATED. Old name for `inlet`.
            tspan_list (list of tspan): The timespans of steady concentration
                or pressure
            carrier_mol (str): The name of the molecule of the carrier gas if
                a dilute analyte is used. Calibration assumes total flux of the
                capillary is the same as the flux of pure carrier gas. Defaults
                to None.
            mol_conc_ppm (float or list): Concentration of the dilute analyte in
                the carrier gas in ppm. Defaults to None. Accepts float (for pressure
                calibration) or list for concentration calibration. If list needs
                to be same length as tspan_list or selector_list.
            p_inlet (float, list): Pressure at the inlet (Pa). Overwrites the pressure
                inherent to self (i.e. the MSInlet object). Accepts float (for conc.
                calibration) or list for pressure calibration. If list, then
                needs to be same length as tspan_list or selector_list.
            tspan_bg (tspan): The time to use as a background
            ax (Axis): The axis on which to plot the ms_calibration curve result.
                Defaults to a new axis.
            axis_measurement (Axis): The MS plot axes to highlight the
                ms_calibration on. Defaults to None.
            remove_bg_on_axis_measurement (bool):
                Whether the plot on axis_measurement is showing raw data or bg
                subtracted data. Defaults to True, i.e. plotting data with the
                same bg subtraction as used for the calibration.
            return_ax (bool): Whether to return the axis on which the calibration
                curve is plotted together with the MSCalResult. Defaults to False.

        Return MSCalResult(, Axis): The result of the MS calibration (and calibration
            curve axis if requested) based on flux calculation during selected time
            periods.
        TODO: automatically recognize the pressure from measurement (if available)
        """
        # prepare three lists to loop over to determine molecule flux in the
        # different periods of steady gas composition
        if not isinstance(mol_conc_ppm, list):
            mol_conc_ppm_list = [mol_conc_ppm for x in tspan_list]
        else:
            mol_conc_ppm_list = mol_conc_ppm
        if isinstance(p_inlet, list):
            p_list = p_inlet
        else:
            p_list = [p_inlet for _ in tspan_list]
        if not len(mol_conc_ppm_list) == len(p_list) == len(tspan_list):
            raise ValueError(
                "Length of input lists for concentrations"
                " and tspan or pressures and tspan is not equal"
            )
        S_list = []
        n_dot_list = []
        if carrier_mol:
            if None not in mol_conc_ppm_list:
                cal_type = "carrier_gas_flux_calibration_curve"
            else:
                raise QuantificationError(
                    "Cannot use carrier gas calibration without analyte"
                    " concentration. 'mol_conc_ppm' is missing. For a pure gas,"
                    "use 'mol' instead of 'carrier_mol' and don't give a 'mol_conc_ppm'"
                )
        elif None not in mol_conc_ppm_list:
            raise QuantificationError(
                "Cannot use carrier gas calibration without carrier"
                " gas definition. 'carrier_mol' is missing. For a pure gas,"
                "use 'mol' instead of 'carrier_mol' and don't give a 'mol_conc_ppm'"
            )
        else:
            cal_type = "gas_flux_calibration_curve"
            # redefine mol_conc_ppm_list to compensate for unit correction done
            # in the calculation of n_dot below
            mol_conc_ppm = 10**6
            mol_conc_ppm_list = [mol_conc_ppm for x in tspan_list]
            # specify that the gas given as mol is now the carrier_mol
            carrier_mol = mol
        for tspan, mol_conc_ppm, pressure in zip(tspan_list, mol_conc_ppm_list, p_list):
            t, S = measurement.grab_signal(mass, tspan=tspan, tspan_bg=tspan_bg)
            if axis_measurement:
                if remove_bg_on_axis_measurement:
                    t_plot, S_plot = t, S
                else:
                    t_plot, S_plot = measurement.grab_signal(mass, tspan=tspan)
                axis_measurement.plot(
                    t_plot, S_plot, color=STANDARD_COLORS[mass], linewidth=5
                )
            n_dot = (
                inlet.calc_n_dot_0(gas=carrier_mol, p=pressure) * mol_conc_ppm / 10**6
            )
            S_list.append(np.mean(S))
            n_dot_list.append(n_dot)
        n_dot_vec = np.array(n_dot_list)
        S_vec = np.array(S_list)
        pfit = np.polyfit(n_dot_vec, S_vec, deg=1)
        F = pfit[0]
        if ax:
            color = STANDARD_COLORS[mass]
            if ax == "new":
                ax = measurement.plotter.new_ax(
                    xlabel="molecule flux / [nmol/s]", ylabel="signal / [nA]"
                )
            ax.plot(n_dot_vec * 1e9, S_vec * 1e9, "o", color=color)
            n_dot_fit = np.array([0, max(n_dot_vec)])
            S_fit = n_dot_fit * pfit[0] + pfit[1]
            ax.plot(n_dot_fit * 1e9, S_fit * 1e9, "--", color=color)
        cal_result = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type=cal_type,
            F=F,
        )
        return cls(ms_cal_results=[cal_result], measurement=measurement)

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

    def __getattr__(self, attr):
        """A MSCalibration with one cal result can be dropped in for that result:"""
        ms_cal_results = super().__getattr__("ms_cal_results")
        if len(ms_cal_results) == 1:
            ms_cal_result = self.ms_cal_results[0]
            if hasattr(ms_cal_result, attr):
                return getattr(ms_cal_result, attr)
        return super().__getattr__(attr)

    @property
    def available_series_names(self):
        return set([f"n_dot_{mol}" for mol in self.mol_list])

    def calculate_series(self, key, measurement=None):
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
                # Calculators just return None when they can't get what's requested.
                return
            try:
                signal_series = measurement[mass]
            except SeriesNotFoundError:
                # Calculators just return None when they can't get what's requested.
                # It could be that another calibration contains a sensitivity factor
                # for the desired mol at an available mass.
                return
            y = signal_series.data
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
            raise QuantificationError(f"{self!r} has no sensitivity factor for {mol}")
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
                f"{self!r} has no sensitivity factor for {mol} at {mass}"
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

    @classmethod
    def from_siq(cls, siq_calculator):

        # A complication is that it can be either a Calibration or a SensitivityList.
        # Either way, the sensitivity factors are in `sf_list`:
        ms_cal_results = [MSCalResult.from_siq(cal) for cal in siq_calculator.sf_list]
        # if it's a Calibration, we want the metadata:
        try:
            calibration = cls(
                name=siq_calculator.name,
                date=siq_calculator.date,
                setup=siq_calculator.setup,
                ms_cal_results=ms_cal_results,
            )
        # if not, we just want the data:
        except AttributeError:
            calibration = cls(ms_cal_results=ms_cal_results)
        return calibration

    def to_siq(self):
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSCalibration.to_siq` only works when using "
                "`spectro_inlets_quantification` "
                "(`ixdat.plugins.activate_siq()`). "
            )
        cal_list = [cal.to_siq() for cal in self.ms_cal_results]
        return plugins.siq.Calculator(
            name=self.name,
            date=self.date,
            setup=self.setup,
            cal_list=cal_list,
        )


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

    def update_l_cap(self, gases=None):
        """Update self.l_cap from average of values in dict l_cap_eff

        Args:
            gases (list): List of gases to average l_cap, default all
        Returns:
            float: Averaged effective capilllary length in [m]
        """
        gases = gases or []
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

        # TODO: make it so that DYNAMIC_VISCOSITIES[gas] can just be a float if someone
        #   enters it without having access to the temperature-dependent values.
        if T < DYNAMIC_VISCOSITIES[gas][0, 0] or T > DYNAMIC_VISCOSITIES[gas][-1, 0]:
            warnings.warn(
                "Insufficient data in constants.py to appropriately estimate "
                f"the dynamic viscosity for {gas} at temperature: {T}K",
                stacklevel=2,
            )
        _eta_v = DYNAMIC_VISCOSITIES[gas][:, 1]  # list of known eta(T) for 'gas'
        _eta_T = DYNAMIC_VISCOSITIES[gas][:, 0]  # list of paired Ts for eta(T)

        eta = np.interp(T, _eta_T, _eta_v)  # dynamic viscosity of gas at T in [Pa*s]

        s = MOLECULAR_DIAMETERS[gas]  # molecule diameter in [m]
        m = MOLAR_MASSES[gas] * 1e-3 / AVOGADRO_CONSTANT  # molecule mass in [kg]

        d = ((w_cap * h_cap) / pi) ** 0.5 * 2
        # d = 4.4e-6  #used in Henriksen2009
        a = d / 2
        p_1 = p
        lambda_ = d  # defining the transitional pressure
        # ...from setting mean free path equal to capillary d
        p_t = BOLTZMANN_CONSTANT * T / (2**0.5 * pi * s**2 * lambda_)
        p_2 = 0
        p_m = (p_1 + p_t) / 2  # average pressure in the transitional flow region
        v_m = (8 * BOLTZMANN_CONSTANT * T / (pi * m)) ** 0.5
        # a reciprocal velocity used for short-hand:
        nu = (m / (BOLTZMANN_CONSTANT * T)) ** 0.5

        # ... and now, we're ready for the capillary equation.
        #   (need to turn of black and flake8 for tolerable format)
        # fmt: off
        #   Equation 4.10 of Daniel Trimarco's PhD Thesis:
        N_dot = (                                                               # noqa
            1 / (BOLTZMANN_CONSTANT * T) * 1 / l_cap * (                        # noqa
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
        n_dot = N_dot / AVOGADRO_CONSTANT
        return n_dot

    # --- METHODS WHICH HAVE BEEN MOVED TO `Calculator` CLASSES ---- #

    @deprecate(
        last_supported_release="0.2.5",
        update_message=("`gas_flux_calibration` is now a method of `MSCalibration`"),
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def gas_flux_calibration(self, measurement, mol, mass, *args, **kwargs):
        return MSCalibration.gas_flux_calibration(
            measurement=measurement, mol=mol, mass=mass, inlet=self, *args, **kwargs
        )

    @deprecate(
        last_supported_release="0.2.5",
        update_message=(
            "`gas_flux_calibration_curve` is now a method of `MSCalibration`"
        ),
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def gas_flux_calibration_curve(
        self,
        measurement,
        mol,
        mass,
        *args,
        **kwargs,
    ):
        return MSCalibration.gas_flux_calibration_curve(
            measurement=measurement, mol=mol, mass=mass, inlet=self, *args, **kwargs
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
    extra_column_attrs = {
        **MSMeasurement.extra_column_attrs,
        **SpectroMeasurement.extra_column_attrs,
    }
    default_plotter = MSSpectroPlotter
    default_exporter = MSSpectroExporter

    # FIXME: https://github.com/ixdat/ixdat/pull/166#discussion_r1486023530
