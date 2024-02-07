"""Module for representation and analysis of MS measurements"""

import re
import numpy as np
import json  # FIXME: This is for MSCalibration.export, but shouldn't have to be here.
import warnings

from ..measurements import Measurement, Calibration
from ..spectra import Spectrum, SpectrumSeries, SpectroMeasurement
from ..plotters.ms_plotter import MSPlotter, MSSpectroPlotter, STANDARD_COLORS
from ..exporters.ms_exporter import MSExporter, MSSpectroExporter
from ..exceptions import QuantificationError
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


def _with_siq_quantifier(method):
    """Decorate a Measurement-copying method to copy the quantifier to the new one.

    This means that the quantifier doesn't need to be re-applied when getting a
    subset of the measurement.

    Method is a method of a Measurement which returns a Measurement of the same type,
    such as cut() and multicut() and the methods that use these ones, such as select,
    and indexing for CyclicVoltammogram objects.
    """

    def method_with_siq_quantifier(*args, **kwargs):
        new_measurement = method(*args, **kwargs)
        old_measurement = args[0]
        if old_measurement.siq_quantifier:
            new_measurement.set_siq_quantifier(old_measurement.siq_quantifier)
        return new_measurement

    return method_with_siq_quantifier


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    # FIXME: tspan_bg should be column of a DataTreater, not MSMeasurement:
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
        new_signal_bgs = {}
        for mass in mass_list:
            if mass in self.signal_bgs:
                new_signal_bgs[mass] = 0
        self.add_calibration(MSCalibration(signal_bgs=new_signal_bgs))

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
            item (str): Name of the signal. If `item` has the form f"n_dot_{mol}", then
                grab_flux(mol) is returned.
            tspan (list): Timespan for which the signal is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if
                available. This is special to MSMeasurement.
                Defaults to False, but in grab_flux it defaults to True.
            include_endpoints (bool): Whether to ensure tspan[0] and tspan[-1] are in t
        """
        if plugins.use_siq and item.startswith("n_dot_"):
            return self.grab_flux(
                item.removeprefix("n_dot_"),
                tspan=tspan,
                tspan_bg=tspan_bg,
                include_endpoints=include_endpoints,
            )
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

        if plugins.use_siq:
            # We have to calculate the fluxes of all the mols and masses in the
            # quantifier's sensitivity matrix. But this method only returns one.
            # TODO: The results should therefore be cached. But how to know when they
            #   need to be recalculated?
            t, n_dots = self.grab_siq_fluxes(
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                include_endpoints=include_endpoints,
            )
            return t, n_dots[mol]

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
            # `MSCalibration.calibrate_series()` interprets item names starting with
            # "n_" as molecule fluxes, and checks itself for a sensitivity factor.
            f"n_dot_{mol}",
            tspan=tspan,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
            include_endpoints=include_endpoints,
        )

    def grab_siq_fluxes(
        self, tspan=None, tspan_bg=None, remove_background=False, include_endpoints=False
    ):
        """Return a time vector and a dictionary with all the quantified fluxes

        Args:
            tspan (list): Timespan for which the signal is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            remove_background (bool): Whether to remove a pre-set background if available
                Defaults to True..
            include_endpoints (bool): Whether to interpolate for tspan[0] and tspan[-1]
        """
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSMeasurement.grab_siq_fluxes` only works when using "
                "`spectro_inlets_quantification` "
                "(`ixdat.plugins.activate_siq()`). "
            )
        sm = self._siq_quantifier.sm
        signals = {}
        t = None
        for mass in sm.mass_list:
            if t is None:
                t, S = self.grab(
                    mass,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=include_endpoints,
                )
            else:
                S = self.grab_for_t(
                    mass,
                    t=t,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                )
            signals[mass] = S
        n_dots = sm.calc_n_dot(signals=signals)
        return t, n_dots

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

    @deprecate(
        "0.2.6",
        "Use `inlet` instead. Or consider using `siq_gas_flux_calibration` "
        "with the `spectro_inlets_quantification` package.",
        "0.3",
        kwarg_name="chip",
    )
    def gas_flux_calibration(
        self,
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
        if plugins.use_siq:
            warnings.warn(
                "spectro_inlets_quantification is active but you are using the native "
                "ixdat version of `MSMeasurement.gas_flux_calibration`"
            )
        t, S = self.grab_signal(mass, tspan=tspan, tspan_bg=tspan_bg)
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
        return MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type=cal_type,
            F=F,
        )

    def gas_flux_calibration_curve(
        self,
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
        if plugins.use_siq:
            warnings.warn(
                "spectro_inlets_quantification is active but you are using the native "
                "ixdat version of `MSMeasurement.siq_gas_flux_calibration_curve`"
            )
        return self._gas_flux_calibration_curve(
            mol=mol,
            mass=mass,
            inlet=inlet,
            chip=chip,
            tspan_list=tspan_list,
            carrier_mol=carrier_mol,
            mol_conc_ppm=mol_conc_ppm,
            p_inlet=p_inlet,
            tspan_bg=tspan_bg,
            ax="new",
            axis_measurement=axis_measurement,
            remove_bg_on_axis_measurement=remove_bg_on_axis_measurement,
            return_ax=return_ax,
        )

    @deprecate(
        "0.2.6",
        "Use `inlet` instead. Or consider using `siq_gas_flux_calibration_curve` "
        "with the `spectro_inlets_quantification` package.",
        "0.3",
        kwarg_name="chip",
    )
    def _gas_flux_calibration_curve(
        self,
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
        """Helper function. See gas_flux_calibraiton_curve for argument descriptions."""

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
            t, S = self.grab_signal(mass, tspan=tspan, tspan_bg=tspan_bg)
            if axis_measurement:
                if remove_bg_on_axis_measurement:
                    t_plot, S_plot = t, S
                else:
                    t_plot, S_plot = self.grab_signal(mass, tspan=tspan)
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
                ax = self.plotter.new_ax(
                    xlabel="molecule flux / [nmol/s]", ylabel="signal / [nA]"
                )
            ax.plot(n_dot_vec * 1e9, S_vec * 1e9, "o", color=color)
            n_dot_fit = np.array([0, max(n_dot_vec)])
            S_fit = n_dot_fit * pfit[0] + pfit[1]
            ax.plot(n_dot_fit * 1e9, S_fit * 1e9, "--", color=color)
        cal = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type=cal_type,
            F=F,
        )
        if return_ax:
            return cal, ax
        return cal

    def siq_gas_flux_calibration(self, mol, mass, tspan, chip=None):
        """Simple pure-gas flux calibration, using `spectro_inlets_quantification`

        Args:
            mol (str): Name of molecule to be calibrated (e.g. "He")
            mass (str): Mass at which to calibrate it (e.g. "M4")
            tspan (timespan): A timespan during which the pure gas is in the chip
            chip (Chip, optional): An object defining the capillary inlet, if different
                than the standard chip assumed by the external package.

        Returns CalPoint: An object from `spectro_inlets_quantification`,
           representing the calibration result
        """
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSMeasurement.siq_gas_flux_calibration` only works when using "
                "`spectro_inlets_quantification` "
                "(`ixdat.options.activate_siq()`). For native ixdat MS quantification, "
                "use `gas_flux_calibration` instead."
            )
        Chip = plugins.siq.Chip
        CalPoint = plugins.siq.CalPoint

        chip = chip or Chip()
        n_dot = chip.calc_n_dot_0(gas=mol)
        S = self.grab_signal(mass, tspan=tspan)[1].mean()
        F = S / n_dot
        return CalPoint(mol=mol, mass=mass, F=F, F_type="capillary", date=self.yyMdd)

    def siq_gas_flux_calibration_curve(
        self,
        mol,
        mass,
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
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            tspan_list (list of tspan): The timespans of steady concentration
                or pressure
            carrier_mol (str): The name of the molecule of the carrier gas if
                a dilute analyte is used. Calibration assumes total flux of the
                capillary is the same as the flux of pure carrier gas. Defaults
                to None.
            mol_conc_ppm (float, list): Concentration of the dilute analyte in
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
            axis_measurement (Axes): The MS plot axes to highlight the
                ms_calibration on. Defaults to None. These axes are not returned.
            remove_bg_on_axis_measurement (bool):
                Whether the plot on axis_measurement is showing raw data or bg
                subtracted data. Defaults to True, i.e. plotting data with the
                same bg subtraction as used for the calibration.
            return_ax (bool): Whether to return the axis on which the calibration
                curve is plotted together with the MSCalResult. Defaults to False.

        Returns CalPoint: An object from `spectro_inlets_quantification`,
           representing the calibration result

        TODO: automatically recognize the pressure from measurement (if available)
        """
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSMeasurement.siq_gas_flux_calibration` only works when using "
                "`spectro_inlets_quantification`"
                "(`ixdat.options.activate_siq()`). "
                "For native ixdat MS quantification, use `gas_flux_calibration`"
                "instead."
            )
        Chip = plugins.siq.Chip

        chip = chip or Chip()

        cal, ax = self._gas_flux_calibration_curve(
            inlet=chip,
            mol=mol,
            mass=mass,
            tspan_list=tspan_list,
            carrier_mol=carrier_mol,
            mol_conc_ppm=mol_conc_ppm,
            p_inlet=p_inlet,
            tspan_bg=tspan_bg,
            ax=ax,
            axis_measurement=axis_measurement,
            remove_bg_on_axis_measurement=remove_bg_on_axis_measurement,
            return_ax=True,
        )

        cal = cal.to_siq()
        if return_ax:
            return cal, ax
        return cal

    def siq_multicomp_gas_flux_calibration(
        self, mol_list, mass_list, gas, tspan, gas_bg=None, tspan_bg=None, chip=None
    ):
        """Calibration of multiple components of a calibration gas simultaneously

        Uses a matrix equation and the reference spectra in the molecule data files.

        The results are only as accurate as the reference spectrum used. For this reason,
        this method is a last resort and it is recommended *not* to use a multicomponent
        calibration gas. Instead, get a separate calibration gas for each molecule to
        be calibrated.

        Here is an explanation of the math used in this method:

        The fundamental matrix equation is:
          S_vec = F_mat @ n_dot_vec
        Elementwise, this is:
         S_M = sum_i ( F^i_M * n_dot^i )
        Rewrite to show that sensitivity factors follow each molecule's spectrum:
         S_M = sum_i (F_weight_i * spectrum^i_M * n_dot^i)
        And regroup the parts that only depend on the molecule (^i):
         S_M = sum_i (spectrum^i_M * (F_weight^i * n_dot^i))
         S_M = sum_i (spectrum^i_M * sensitivity_flux^i)
        Change back into a matrix equation, and solve it:
         S_vec = spectrum_mat @ sensitivity_flux_vec
         sensitivity_flux_vec = spectrum_mat^-1 @ S_vec   # eq. 1
        Ungroup the part we grouped before (the "sensitivity_flux"):
         F_weight^i = sensitivity_flux^i / n_dot^i        # eq. 2
        And, in the end, each sensitivity factor is:
         F_M^i = F_weight^i * spectrum^i_M                # eq. 3

        Equations 1, 2, and 3 are implemented in the code of this method.

        Args:
            mol_list (list of str): List of the names of the molecules to calibrate
            mass_list (list of str): List of the masses to calibrate
            gas (Gas, dict, or str): Composition of the calibration gas, e.g.
               {"Ar": 0.95, "H2": 0.05} for 5% H2 in Ar
            tspan (Timespan): Timespan during which the calibration gas is in the chip
            gas_bg (Gas, dict, or str): Composition of the background gas
            tspan_bg (Timespan): Timespan during which the background gas is in the chip
            chip (Chip, optional): object describing the MS capillary, if different than
               the standard chip in the MS quantification package

        Returns Calibration: An object from `spectro_inlets_quantification`,
           representing all the calibration results from the calibration.
        """
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSMeasurement.siq_multicomp_gas_flux_calibration` "
                "only works when using `spectro_inlets_quantification` "
                "(`ixdat.plugins.activate_siq()`). "
            )
        Chip = plugins.siq.Chip
        CalPoint = plugins.siq.CalPoint
        Calibration = plugins.siq.Calibration

        chip = chip or Chip()
        chip.gas = gas
        flux = chip.calc_n_dot()

        chip_bg = chip or Chip()
        chip_bg.gas = gas_bg
        flux_bg = chip_bg.calc_n_dot()

        delta_flux_list = []
        for mol in mol_list:
            delta_flux = flux.get(mol, 0) - flux_bg.get(mol, 0)
            delta_flux_list.append(delta_flux)
        delta_flux_vec = np.array(delta_flux_list)

        delta_signal_list = []
        for mass in mass_list:
            S = self.grab_signal(mass, tspan=tspan)[1].mean()
            S_bg = self.grab_signal(mass, tspan=tspan_bg)[1].mean()
            delta_S = S - S_bg
            delta_signal_list.append(delta_S)
        delta_signal_vec = np.array(delta_signal_list)

        spectrum_vec_list = []
        for mol in mol_list:
            spectrum = chip.gas.mdict[mol].norm_spectrum
            spectrum_vec = np.array([spectrum.get(mass, 0) for mass in mass_list])
            spectrum_vec_list.append(spectrum_vec)
        spectrum_mat = np.stack(spectrum_vec_list).transpose()

        inverse_spectrum_mat = np.linalg.inv(spectrum_mat)
        sensitivity_flux_vec = inverse_spectrum_mat @ delta_signal_vec  # eq. 1
        F_weight_vec = sensitivity_flux_vec / delta_flux_vec  # eq. 2

        cal_list = []
        for i, mol in enumerate(mol_list):
            for M, mass in enumerate(mass_list):
                F = F_weight_vec[i] * spectrum_mat[M, i]  # eq. 3
                if F:
                    cal = CalPoint(
                        mol=mol, mass=mass, F=F, F_type="capillary", date=self.yyMdd
                    )
                    cal_list.append(cal)

        return Calibration(cal_list=cal_list)

    def set_siq_quantifier(
        self,
        quantifier=None,
        calibration=None,
        mol_list=None,
        mass_list=None,
        carrier="He",
    ):
        """Set the `spectro_inlets_quantification` quantifier.

        The Quantifier is an object with the method `calc_n_dot`, which takes a
        dictionary of signals or signal vectors in [A] and return a dictionary of
        molecular fluxes in [mol/s].
        The quantifier typically does this by solving the linear equations of
        S_M = sum_i ( F_M^i * n_dot^i )
        Where n_dot^i is the flux to the vacuum chamber of molecule i in [mol/s], S_M
        is the signal at mass M in [A], and F_M^i is the *sensitivity factor* of molecule
        i at mass M.
        The quantifier thus needs access to a set of sensitivity factors.

        The quantifier can be built in this method (avoiding explicit import of the
        `spectro_inlets_quantification` package) by providing the sensitivity factors
        in the form of a `Calibration` (which can be obtained from e.g.
        MSMeasurement.multicomp_gas_flux_cal) and the specification of which ones to
        use by `mol_list` and `mass_list`.
        The quantifier will always use all the masses in `mass_list` to solve for the
        flux of all the mols in `mol_list`.

        The argument `carrier` is required by some quantifiers but only used if
        partial pressures before the MS inlet are required (`quantifier.calc_pp`)

        Quantification is only as accurate as your sensitivity factors!

        Args:
            quantifier (Quantifier): The quantifier, if prepared before method call.
               No additional arguments needed. Otherwise, the following three are needed:
            calibration (Calibration): The calibration to build the quantifier with
            mol_list (list of str): The list of molecules to use in flux calculations.
               These should all be represented in the Calibration. If not provided,
               we'll use all the mols in the Calibration.
            mass_list (list of str): The list of masses to use in flux calculations.
               These should all be represented in the Calibration. If not provided,
               we'll use all the masses in the Calibration.
            carrier (optional, str): The carrier gas in the experiment. Defaults to "He".
        """
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSMeasurement.set_siq_quantifier` only works when using "
                "`spectro_inlets_quantification` "
                "(`ixdat.options.activate_siq()`). "
                "For native ixdat MS quantification, use `MSMeasurement.calibrate`"
            )
        Quantifier = plugins.siq.Quantifier

        if quantifier:
            self._siq_quantifier = quantifier
        else:
            mol_list = mol_list or calibration.mol_list
            mass_list = mass_list or calibration.mass_list
            self._siq_quantifier = Quantifier(
                calibration=calibration,
                mol_list=mol_list,
                mass_list=mass_list,
                carrier=carrier,
            )

    @property
    def siq_quantifier(self):
        return self._siq_quantifier

    cut = _with_siq_quantifier(Measurement.cut)
    multicut = _with_siq_quantifier(Measurement.multicut)


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
    def from_siq(cls, siq_calibration):

        # A complication is that it can be either a Calibration or a SensitivityList.
        # Either way, the sensitivity factors are in `sf_list`:
        ms_cal_results = [MSCalResult.from_siq(cal) for cal in siq_calibration.sf_list]
        # if it's a Calibration, we want the metadata:
        try:
            calibration = cls(
                name=siq_calibration.name,
                date=siq_calibration.date,
                setup=siq_calibration.setup,
                ms_cal_results=ms_cal_results,
            )
        # if not, we just want the data:
        except AttributeError:
            calibration = cls(ms_cal_results=ms_cal_results)
        return calibration

    def to_siq(self):
        if not plugins.use_siq:
            raise QuantificationError(
                "`MSCalPoint.to_siq` only works when using "
                "`spectro_inlets_quantification` "
                "(`ixdat.options.activate_siq()`). "
                "For native ixdat MS quantification, use `MSMeasurement.calibrate`"
            )
        cal_list = [cal.to_siq() for cal in self.ms_cal_results]
        return plugins.siq.Calibration(
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

    @deprecate(
        last_supported_release="0.2.5",
        update_message=("`gas_flux_calibration` is now a method of `MSMeasurement`"),
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
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
        return measurement.gas_flux_calibration(
            mol=mol,
            inlet=self,
            mass=mass,
            tspan=tspan,
            tspan_bg=tspan_bg,
            ax=ax,
            carrier_mol=carrier_mol,
            mol_conc_ppm=mol_conc_ppm,
        )

    @deprecate(
        last_supported_release="0.2.5",
        update_message=(
            "`gas_flux_calibration_curve` is now a method of `MSMeasurement`"
        ),
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def gas_flux_calibration_curve(
        self,
        measurement,
        mol,
        mass,
        tspan_list=None,
        carrier_mol=None,
        mol_conc_ppm=None,
        p_inlet=None,
        tspan_bg=None,
        ax="new",
        axis_measurement=None,
        return_ax=False,
    ):
        return measurement.gas_flux_calibration_curve(
            mol=mol,
            inlet=self,
            mass=mass,
            tspan_list=tspan_list,
            tspan_bg=tspan_bg,
            ax=ax,
            carrier_mol=carrier_mol,
            mol_conc_ppm=mol_conc_ppm,
            p_inlet=p_inlet,
            axis_measurement=axis_measurement,
            return_ax=return_ax,
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

    # FIXME: automate this in inheritance of hyphenated techniques:
    extra_column_attrs = {
        "spectro_measurements": {"spectrum_id"},
        "ms_measurement": {"tspan_bg"},
    }
    default_plotter = MSSpectroPlotter
    default_exporter = MSSpectroExporter

    # FIXME: this shouldn't be necessary. See #164.
    cut = _with_siq_quantifier(SpectroMeasurement.cut)
    multicut = _with_siq_quantifier(SpectroMeasurement.multicut)
