"""Module for representation and analysis of EC measurements"""

import numpy as np

from ..measurements import Measurement, append_series, time_shifted
from ..data_series import ValueSeries, ConstantValue
from ..exceptions import SeriesNotFoundError
from ..exporters.ec_exporter import ECExporter


class ECMeasurement(Measurement):
    """Class implementing electrochemistry measurements

    TODO: Implement a unit library for current and potential, A_el and RE_vs_RHE
    TODO:   so that e.g. current can be seamlessly normalized to mass OR area.

    The main job of this class is making sure that the ValueSeries most essential for
    visualizing and normal electrochemistry measurements (i.e. excluding impedance spec.,
    RRDE, etc, which would need new classes) are always available in the correct form as
    the measurement is added with others, reduced to a selection, calibrated and
    normalized, etc. These most important ValueSeries are:

    - `potential`: The working-electrode potential typically in [V].
      If `ec_meas` is an `ECMeasurement`, then `ec_meas["potential"]` always returns a
      `ValueSeries` characterized by:

        - calibrated and/or corrected, if the measurement has been calibrated with the
          reference electrode potential (`RE_vs_RHE`, see `calibrate`) and/or corrected
          for ohmic drop (`R_Ohm`, see `correct_ohmic_drop`).
        - A name that makes clear any calibration and/or correction
        - Data which spans the entire timespan of the measurement - i.e. whenever EC data
          is being recorded, `potential` is there, even the name of the raw
          `ValueSeries` (what the acquisition software calls it) changes. Indeed
          `ec_meas["potential"].tseries` is the measurement's definitive time variable.

    - `current`: The working-electrode current typically in [mA] or [mA/cm^2].
      `ec_meas["current"]` always returns a `ValueSeries` characterized by:

        - normalized if the measurement has been normalized with the electrode area
          (`A_el`, see `normalize`)
        - A name that makes clear whether it is normalized
        - Data which spans the entire timespan of the measurement

    - `selector`: A counter series distinguishing sections of the measurement program.
      This is essential for analysis of complex measurements as it allows for
      corresponding parts of experiments to be isolated and treated identically.
      `selector` in `ECMeasurement` is defined to incriment each time one or more of
      the following changes:

        - `loop_number`: A parameter saved by some potentiostats (e.g. BioLogic) which
          allow complex looped electrochemistry programs.
        - `file_number`: The id of the component measurement from which each section of
          the data (the origin of each `ValueSeries` concatenated to `potential`)
        - `cycle_number`: An incrementer within a file saved by a potentiostat.

    The names of these ValueSeries, which can also be used to index the measurement, are
    conveniently available as properties:

    - `ec_meas.t_str` is the name of the definitive time, which corresponds to potential.
    - `ec_meas.E_str` is the name of the raw potential
    - `ec_meas.V_str` is the name to the calibrated and/or corrected potential
    - `ec_meas.I_str` is the name of the raw current
    - `ec_meas.J_str` is the name of the normalized current.
    - `ec_meas.sel_str` is the name of the default selector, i.e. "selector"

    Numpy arrays from important `DataSeries` are also directly accessible via attributes:

    - `ec_meas.t` for `ec_meas["potential"].t`
    - `ec_meas.v` for `ec_meas["potential"].data`
    - `ec_meas.j` for `ec_meas["current"].data`

    `ECMeasurement` comes with an `ECPlotter` which either plots `potential` and
    `current` against time (`ec_meas.plot_measurement()`) or plots `current` against
    `potential (`ec_meas.plot_vs_potential()`).

    It turns out that keeping track of current, potential, and selector when combining
    datasets is enough of a job to fill a class. Thus, the more exciting
    electrochemistry-related functionality should be implemented in inheriting classes
    such as `CyclicVoltammagram`.
    """

    extra_column_attrs = {
        "ec_meaurements": {
            "ec_technique",
            "RE_vs_RHE",
            "R_Ohm",
            "raw_potential_names",
            "A_el",
            "raw_current_names",
        }
    }

    def __init__(
        self,
        name,
        *,
        technique=None,
        metadata=None,
        s_ids=None,
        series_list=None,
        m_ids=None,
        component_measurements=None,
        reader=None,
        plotter=None,
        exporter=None,
        sample=None,
        lablog=None,
        tstamp=None,
        ec_technique=None,
        t_str="time / [s]",
        E_str="raw potential / [V]",
        V_str="$U_{RHE}$ / [V]",
        RE_vs_RHE=None,
        R_Ohm=None,
        raw_potential_names=("Ewe/V", "<Ewe>/V"),  # TODO: reader must define this
        I_str="raw current / [mA]",
        J_str="J / [mA cm$^{-2}$]",
        A_el=None,
        raw_current_names=("I/mA", "<I>/mA"),  # TODO: reader must define this
        cycle_names=("cycle number",),
    ):
        """initialize an electrochemistry measurement

        Args:
            name (str): The name of the measurement
                TODO: Decide if metadata needs the json string option.
                TODO:   See: https://github.com/ixdat/ixdat/pull/1#discussion_r546436991
            metadata (dict): Free-form measurement metadata
            technique (str): The measurement technique
            s_ids (list of int): The id's of the measurement's DataSeries, if
                to be loaded (instead of given directly in series_list)
            series_list (list of DataSeries): The measurement's DataSeries
            m_ids (list of int): The id's of the component measurements, if to be
                loaded. None unless this is a combined measurement (typically
                corresponding to more than one file).
            component_measurements (list of Measurements): The measurements of which
                this measurement is a combination
            reader (Reader): The file reader (None unless read from a file)
            plotter (Plotter): The visualization tool for the measurement
            exporter (Exporter): The exporting tool for the measurement
            sample (Sample): The (already loaded) sample being measured
            lablog (LabLog): The log entry with e.g. notes taken during the measurement
            tstamp (float): The nominal starting time of the measurement, used for
                data selection, visualization, and exporting.
            t_str (str): Name of the main time variable (corresponding to potential)
            E_str (str): Name of raw potential (so called because potential is saved
                as "Ewe/V" in biologic .mpt files)
            V_str (str): Name of calibrated potential
            RE_vs_RHE (float): Reference electrode potential in [V] on the RHE scale.
                If RE_vs_RHE is not None, the measurement is considered *calibrated*,
                and will use the calibrated potential `self[self.V_str]` by default
                TODO: Unit
            R_Ohm (float): Ohmic drop in [Ohm]. If R_Ohm is not None, the ohmic drop
                is corrected for when returning potential.
                TODO: Unit
            raw_potential_names (tuple of str): The names of the VSeries which
                represent raw working electrode current. This is typically how the data
                acquisition software saves potential
            I_str (str): Name of raw current
            J_str (str): Name of normalized current
            A_el (float): Area of electrode in [cm^2].
                If A_el is not None, the measurement is considered *normalized*,
                and will use the calibrated current `self[self.J_str]` by default
                TODO: Unit
            raw_current_names (tuple of str): The names of the VSeries which represent
                raw working electrode current. This is typically how the data
                acquisition software saves current.
        """

        calibration = self.calibration if hasattr(self, "calibration") else None
        super().__init__(
            name,
            technique=technique,
            metadata=metadata,
            s_ids=s_ids,
            series_list=series_list,
            m_ids=m_ids,
            component_measurements=component_measurements,
            reader=reader,
            plotter=plotter,
            exporter=exporter,
            sample=sample,
            lablog=lablog,
            tstamp=tstamp,
        )  # FIXME: The super init sets self.calibration to None but I can't see why.
        self.calibration = calibration or ECCalibration(RE_vs_RHE, A_el)
        if RE_vs_RHE is not None:  # if given as an arg RE_vs_RHE trumps calibration
            self.RE_vs_RHE = RE_vs_RHE
        if A_el is not None:
            self.A_el = A_el
        self.ec_technique = ec_technique
        self.t_str = t_str
        self.E_str = E_str
        self.V_str = V_str
        self.R_Ohm = R_Ohm
        self.raw_potential_names = raw_potential_names
        self.I_str = I_str
        self.J_str = J_str
        self.raw_current_names = raw_current_names
        self.cycle_names = cycle_names

        self.sel_str = "selector"
        self.cycle_str = "cycle_number"

        self.plot_vs_potential = self.plotter.plot_vs_potential

        self._raw_potential = None
        self._raw_current = None
        self._selector = None
        self._file_number = None
        if self.potential:
            if all(
                [
                    (current_name not in self.series_names)
                    for current_name in self.raw_current_names
                ]
            ):
                self.series_list.append(
                    ConstantValue(
                        name=self.raw_current_names[0],
                        unit_name="mA",
                        value=0,
                    )
                )
                self._populate_constants()  # So that OCP currents are included as 0.
                # TODO: I don't like this. ConstantValue was introduced to facilitate
                #   ixdat's laziness, but I can't find anywhere else to put the call to
                #   _populate_constants() that can find the right tseries. This is a
                #   violation of laziness as bad as what it was meant to solve.
            if all(
                [
                    (cycle_name not in self.series_names)
                    for cycle_name in self.cycle_names
                ]
            ):
                self.series_list.append(
                    ConstantValue(
                        name=self.cycle_names[0],
                        unit_name=None,
                        value=0,
                    )
                )
                self._populate_constants()  # So that everything has a cycle number

    @property
    def A_el(self):
        return self.calibration.A_el

    @A_el.setter
    def A_el(self, A_el):
        self.calibration.A_el = A_el

    @property
    def RE_vs_RHE(self):
        return self.calibration.RE_vs_RHE

    @RE_vs_RHE.setter
    def RE_vs_RHE(self, RE_vs_RHE):
        self.calibration.RE_vs_RHE = RE_vs_RHE

    def _populate_constants(self):
        """Replace any ConstantValues with ValueSeries on potential's tseries

        TODO: This function flagrantly violates laziness. Not only does it fill up all
            the ConstantValue's with long vectors before they're needed, it also forces
            raw_potential to be built before it is needed.
            A lazier solution is needed.
        """
        for (i, s) in enumerate(self.series_list):
            if isinstance(s, ConstantValue):
                tseries = self.potential.tseries
                series = s.get_vseries(tseries=tseries)
                self.series_list[i] = series

    def __getitem__(self, item):
        """Return the (concatenated) (time-shifted) `DataSeries` with name `item`

        If `item` matches one of the strings for managed series described in the class
        docstring, item retrieval will still first look in `series_list` (see
        `ixdat.Measurement.__getitem__()` for this inherited behavior), but will then
        return the corresponding managed attribute of this `ECMeasurement`. (See
        class docstring.) Only if `item` matches neither these strings nor the names of
        the `DataSeries` in `series_list` is a `SeriesNotFoundError` raised.

        TODO: I would like to decorate this with a with_time_shifted() decorator to
            enforce here (rather than all over as now) that the time is referenced
            to the measurement tstamp. But not obvious to me how the decorator would have
            access to self.tstamp.
        """
        try:
            return super().__getitem__(item)
        except SeriesNotFoundError:
            if item == self.t_str:  # master time (potential's tseries)
                return self.potential.tseries
            if item == self.E_str:  # raw potential
                return self.raw_potential
            elif item == self.V_str:  # (calibrated) (corrected) potential
                return self.potential
            elif item == self.I_str:  # raw current
                return self.raw_current
            elif item == self.J_str:  # (normalized) current
                return self.current
            elif item == self.sel_str:  # selector
                return self.selector
            elif item == "potential":
                return self.potential
            elif item == "current":
                return self.current
            elif item == "raw_potential":
                return self.raw_potential
            elif item == "raw_current":
                return self.raw_current
            elif item == self.potential.name:
                return self.potential
            elif item == self.current.name:
                return self.current
            raise SeriesNotFoundError(f"{self} doesn't have item '{item}'")

    @property
    def raw_potential(self):
        """Return a time-shifted ValueSeries for the raw potential, built first time."""
        if not self._raw_potential:
            try:
                self._find_or_build_raw_potential()
            except SeriesNotFoundError as e:
                print(f"Warning!!! {self} encountered: {e}")
                return
        return time_shifted(self._raw_potential, tstamp=self.tstamp)
        # FIXME. Hidden attributes not scaleable cache'ing

    def _find_or_build_raw_potential(self):
        """Build the raw potential and store it data_series and as self._raw_potential()
        # TODO, it should instead be stored in a `cached_series_list`.

        This works by finding all the series that have names matching the raw potential
        names list `self.raw_potential_names` (which should be provided by the Reader).
        If there is only one, it just shifts it to t=0 at self.tstamp.
        FIXME:
            If there are multiple it appends them with t=0 at self.tstamp. In this
            case it also appends the `TimeSeries` to `series_list` since *bad things
            might happen?* if the `TimeSeries` of a `ValueSeries` in `series_list` is
            not itself in `series_list`. But this results in redundant TimeSeries.
        """
        potential_series_list = [
            s for s in self.series_list if s.name in self.raw_potential_names
        ]
        if len(potential_series_list) == 1:
            self._raw_potential = time_shifted(
                potential_series_list[0], tstamp=self.tstamp
            )
        elif len(potential_series_list) > 1:
            raw_potential = append_series(potential_series_list, tstamp=self.tstamp)
            if self._raw_current and self._raw_current.tseries == raw_potential.tseries:
                # Then we can re-use the tseries from raw_current rather than
                # saving a new one :D
                potential_tseries = self._raw_current.tseries
            else:
                potential_tseries = raw_potential.tseries
                self.series_list.append(potential_tseries)
            self._raw_potential = ValueSeries(
                name=self.E_str,
                data=raw_potential.data,
                unit_name=raw_potential.unit_name,
                tseries=potential_tseries,
            )
            self[
                self.E_str
            ] = self._raw_potential  # TODO: Better cache'ing. This saves.
        else:
            raise SeriesNotFoundError(
                f"{self} does not have a series corresponding to raw potential."
                f" Looked for series with names in {self.raw_potential_names}"
            )

    @property
    def raw_current(self):
        """Return a time-shifted ValueSeries for the raw current, built first time."""
        if not self._raw_current:
            self._find_or_build_raw_current()
        return time_shifted(self._raw_current, tstamp=self.tstamp)
        # FIXME. Hidden attributes not scaleable cache'ing

    def _find_or_build_raw_current(self):
        """Build the raw current and store it data_series and as self._raw_current()

        This works the way as `_find_or_build_raw_potential`. See the docstring there.
        FIXME: it also has the same problems.
        """
        current_series_list = [
            s for s in self.series_list if s.name in self.raw_current_names
        ]
        if len(current_series_list) == 1:
            self._raw_current = time_shifted(current_series_list[0], tstamp=self.tstamp)
        elif len(current_series_list) > 1:
            raw_current = append_series(current_series_list, tstamp=self.tstamp)
            if self._raw_potential and (
                self._raw_potential.tseries == raw_current.tseries
            ):
                # Then we can re-use the tseries from raw_potential rather than
                # saving a new one :D
                current_tseries = self._raw_potential.tseries
            else:
                current_tseries = raw_current.tseries
                self.series_list.append(current_tseries)
            self._raw_current = ValueSeries(
                name=self.I_str,
                data=raw_current.data,
                unit_name=raw_current.unit_name,
                tseries=current_tseries,
            )
            self[
                self.I_str
            ] = self._raw_current  # TODO: better cache'ing. This is saved
        else:
            raise SeriesNotFoundError(
                f"{self} does not have a series corresponding to raw current."
                f" Looked for series with names in {self.raw_current_names}"
            )

    def calibrate(self, RE_vs_RHE=None, A_el=None, R_Ohm=None):
        """Calibrate the EC measurement (all args optional)

        Args:
            RE_vs_RHE (float): reference electode potential on RHE scale in [V]
            A_el (float): electrode area in [cm^2]
            R_Ohm (float): ohmic drop resistance in [Ohm]
        """
        if RE_vs_RHE is not None:  # it can be 0!
            self.calibrate_RE(RE_vs_RHE=RE_vs_RHE)
        if A_el:
            self.normalize_current(A_el=A_el)
        if R_Ohm:
            self.correct_ohmic_drop(R_Ohm=R_Ohm)

    def calibrate_RE(self, RE_vs_RHE):
        """Calibrate the reference electrode by providing `RE_vs_RHE` in [V].

        Return string: The name of the calibrated potential
        """
        self.RE_vs_RHE = RE_vs_RHE
        return self.V_str

    def normalize_current(self, A_el):
        """Normalize current to electrod surface area by providing `A_el` in [cm^2].

        Return string: The name of the normalized current
        """
        self.A_el = A_el
        return self.J_str

    def correct_ohmic_drop(self, R_Ohm):
        """Correct for ohmic drop by providing `R_Ohm` in [Ohm].

        Return string: The name of the corrected potential
        """
        self.R_Ohm = R_Ohm
        return self.V_str

    @property
    def potential(self):
        """The ValueSeries with the ECMeasurement's potential.

        This is result of the following:
        - Starts with `self.raw_potential`
        - if the measurement is "calibrated" i.e. `RE_vs_RHE` is not None: add
        `RE_vs_RHE` to the potential data and change its name from `E_str` to `V_str`
        - if the measurement is "corrected" i.e. `R_Ohm` is not None: subtract
        `R_Ohm` times the raw current from the potential and add " (corrected)" to
        its name.
        """

        if self.V_str in self.series_names:
            return self[self.V_str]
        raw_potential = self.raw_potential
        if self.RE_vs_RHE is None and self.R_Ohm is None:
            return raw_potential
        fixed_V_str = raw_potential.name
        fixed_potential_data = raw_potential.data
        fixed_unit_name = raw_potential.unit_name
        if self.RE_vs_RHE is not None:
            fixed_V_str = self.V_str
            fixed_potential_data = fixed_potential_data + self.RE_vs_RHE
            fixed_unit_name = "V <RHE>"
        if self.R_Ohm:
            fixed_V_str += " (corrected)"
            fixed_potential_data = (
                fixed_potential_data
                - self.R_Ohm * self.grab_for_t("raw_current", raw_potential.t) * 1e-3
            )  # TODO: Units. The 1e-3 here is to bring raw_current.data from [mA] to [A]
        return ValueSeries(
            name=fixed_V_str,
            data=fixed_potential_data,
            unit_name=fixed_unit_name,
            tseries=raw_potential.tseries,
        )  # TODO: Better cache'ing. This is not cached at all.

    @property
    def current(self):
        """The ValueSeries with the ECMeasurement's current.

        This is result of the following:
        - Starts with `self.raw_current`
        - if the measurement is "normalized" i.e. `A_el` is not None: divide the current
        data by `A_el`, change its name from `I_str` to `J_str`, and add `/cm^2` to
        its unit.
        """
        if self.J_str in self.series_names:
            return self[self.J_str]
        raw_current = self.raw_current
        if self.A_el is None:
            return raw_current
        else:
            return ValueSeries(
                name=self.J_str,
                data=raw_current.data / self.A_el,
                unit_name=raw_current.unit_name + "/cm^2",
                tseries=raw_current.tseries,
            )

    def grab_potential(self, tspan=None, cal=True):
        """Return t and potential (if cal else raw_potential) [V] vectors cut by tspan"""
        if cal:
            return self.grab("potential", tspan=tspan)
        else:
            return self.grab("raw_potential", tspan=tspan)

    def grab_current(self, tspan=None, norm=True):
        """Return t [s] and current (if cal else raw_current) [V] vectors cut by tspan"""
        if norm:
            return self.grab("current", tspan=tspan)
        else:
            return self.grab("raw_current", tspan=tspan)

    @property
    def t(self):
        """The definitive time np array of the measurement, corresponding to potential"""
        return self.potential.t.copy()

    @property
    def v(self):
        """The potential [V] numpy array of the measurement"""
        return self.potential.data.copy()

    @property
    def j(self):
        """The current ([mA] or [mA/cm^2]) numpy array of the measurement"""
        return self.current.data.copy()

    @property
    def plotter(self):
        """The default plotter for ECMeasurement is ECPlotter"""
        if not self._plotter:
            from ..plotters.ec_plotter import ECPlotter

            self._plotter = ECPlotter(measurement=self)

        return self._plotter

    @property
    def exporter(self):
        """The default plotter for ECMeasurement is ECExporter"""
        if not self._exporter:
            self._exporter = ECExporter(measurement=self)
        return self._exporter

    @property
    def selector(self):
        """The ValuSeries which is used by default to select parts of the measurement.

        See the class docstring for details.
        """
        if self.sel_str not in self.series_names:
            self._build_selector()
        return time_shifted(self[self.sel_str], tstamp=self.tstamp)

    def _build_selector(self, sel_str=None):
        """Build `selector` from `cycle number`, `loop_number`, and `file_number`

        See the class docstring for details.
        """
        sel_str = sel_str or self.sel_str
        changes = np.tile(False, self.t.shape)
        col_list = ["cycle number", "loop_number", "file_number"]
        for col in col_list:
            if col in self.series_names:
                values = self[col].data
                if len(values) == 0:
                    print("WARNING: " + col + " is empty")
                    continue
                elif not len(values) == len(changes):
                    print("WARNING: " + col + " has an unexpected length")
                    continue
                n_down = np.append(
                    values[0], values[:-1]
                )  # comparing with n_up instead puts selector a point ahead
                changes = np.logical_or(changes, n_down < values)
        selector = np.cumsum(changes)
        selector_series = ValueSeries(
            name=sel_str,
            unit_name="",
            data=selector,
            tseries=self.potential.tseries,
        )
        self[self.sel_str] = selector_series  # TODO: Better cache'ing. This gets saved.

    @property
    def cycle_number(self):
        """The cycle number ValueSeries, requires building from component measurements"""
        cycle_series_list = [s for s in self.series_list if s.name in self.cycle_names]
        if len(cycle_series_list) == 1:
            return time_shifted(cycle_series_list[0], tstamp=self.tstamp)
        elif len(cycle_series_list) > 1:
            cycle = append_series(cycle_series_list, tstamp=self.tstamp)
            return ValueSeries(
                name=self.cycle_str,
                data=cycle.data,
                unit_name=cycle.unit_name,
                tseries=cycle.tseries,
            )
        else:
            return
            raise SeriesNotFoundError(
                f"{self} does not have a series corresponding to cycle number."
            )
        # TODO: better cache'ing. This one is not cache'd at all

    @property
    def file_number(self):
        """The file number ValueSeries, requires building from component measurements."""
        if "file_number" in self.series_names:
            self._build_file_number()
        return time_shifted(self["file_number"], tstamp=self.tstamp)

    def _build_file_number(self):
        """Build the """
        file_number_series_list = []
        for m in self.component_measurements:
            vseries = m.potential
            file_number_series = ValueSeries(
                name="file_number",
                unit_name="",
                data=np.tile(m.id, vseries.t.shape),
                tseries=vseries.tseries,
            )
            file_number_series_list.append(file_number_series)
        file_number = append_series(file_number_series_list, tstamp=self.tstamp)
        self[
            "file_number"
        ] = file_number  # TODO: better cache'ing. This one gets saved.

    def as_cv(self):
        """Convert self to a CyclicVoltammagram"""
        from .cv import CyclicVoltammagram

        self_as_dict = self.as_dict()
        self_as_dict["series_list"] = self.series_list
        self_as_dict["technique"] = "CV"
        del self_as_dict["s_ids"]
        # Note, this works perfectly! All needed information is in self_as_dict :)
        return CyclicVoltammagram.from_dict(self_as_dict)


class ECCalibration:
    """A small container for RHE_vs_RE and A_el"""

    def __init__(self, RE_vs_RHE=None, A_el=None):
        self.RE_vs_RHE = RE_vs_RHE
        self.A_el = A_el
