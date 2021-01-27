"""Module for representation and analysis of EC measurements"""

import numpy as np

from ..measurements import Measurement, append_series
from ..data_series import ValueSeries, ConstantValue
from ..exceptions import SeriesNotFoundError


class ECMeasurement(Measurement):
    """Class implementing raw electrochemistry measurements

    TODO:
        Implement a unit library for current and potential, A_el and RE_vs_RHE
        so that e.g. current can be seamlessly normalized to mass OR area.
    """

    extra_column_attrs = {
        "ec_meaurements": {
            "ec_technique",
            "RE_vs_RHE",
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
        metadata_json_string=None,
        s_ids=None,
        series_list=None,
        m_ids=None,
        component_measurements=None,
        reader=None,
        plotter=None,
        exporter=None,
        sample=None,
        sample_name=None,
        lablog=None,
        tstamp=None,
        ec_technique=None,
        E_str="raw potential / [V]",
        V_str="$U_{RHE}$ / [V]",
        RE_vs_RHE=None,
        raw_potential_names=("Ewe/V", "<Ewe>/V"),
        I_str="raw current / [mA]",
        J_str="J / [mA cm$^{-2}$]",
        A_el=None,
        raw_current_names=("I/mA", "<I>/mA"),
        cycle_names=("cycle number",),
    ):
        """initialize an electrochemistry measurement

        Args:
            name (str): The name of the measurement
                TODO: Decide if metadata needs the json string option.
                    See: https://github.com/ixdat/ixdat/pull/1#discussion_r546436991
            metadata (dict): Free-form measurement metadata
            metadata_json_string (str): Free-form measurement metadata as json string
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
            sample_name (str): The name of the sample being measured (will be loaded)
            lablog (LabLog): The log entry with e.g. notes taken during the measurement
            tstamp (float): The nominal starting time of the measurement, used for
                data selection, visualization, and exporting.
            E_str (str): Name of raw potential (so called because potential is saved
                as "Ewe/V" in biologic .mpt files)
            V_str (str): Name of calibrated potential
            RE_vs_RHE (float): Reference electrode potential in [V] on the RHE scale.
                If RE_vs_RHE is not None, the measurement is considered *calibrated*,
                and will use the calibrated potential `self[self.V_str]` by default
                TODO: Unit
            raw_potential_names (tuple of str): The names of the VSeries which
                represent raw working electrode current. This is typically how the data
                acquisition software saves potential
            I_str (str): Name of raw current
            J_str (str): Name of normalized current
            A_el (float): Area of electrode in [cm^2].
                If A_el is not None, the measurement is considered *normalized*,
                and will use the calibrated potential `self[self.V_str]` by default
                TODO: Unit
            raw_current_names (tuple of str): The names of the VSeries which represent
                raw working electrode current. This is typically how the data
                acquisition software saves current.
        """
        super().__init__(
            name,
            technique=technique,
            metadata=metadata,
            metadata_json_string=metadata_json_string,
            s_ids=s_ids,
            series_list=series_list,
            m_ids=m_ids,
            component_measurements=component_measurements,
            reader=reader,
            plotter=plotter,
            exporter=exporter,
            sample=sample,
            sample_name=sample_name,
            lablog=lablog,
            tstamp=tstamp,
        )
        self.ec_technique = ec_technique
        self.E_str = E_str
        self.V_str = V_str
        self.RE_vs_RHE = RE_vs_RHE
        self.raw_potential_names = raw_potential_names
        self.I_str = I_str
        self.J_str = J_str
        self.A_el = A_el
        self.raw_current_names = raw_current_names
        self.cycle_names = cycle_names

        self.sel_str = "selector"
        self.cycle_str = "cycle_number"

        self._raw_potential = None
        self._raw_current = None
        self._selector = None
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
            # TODO: I don't like this. The ConstantValue was introduced to facilitate
            #   ixdat's laziness, but I can't find anywhere else to put the call to
            #   _populate_constants() that can find the right tseries. That's because
            #   once measurements are added together, it's not completely easy to
            #   tell which DataSeries come form the same component measurements.
        if all(
            [(cycle_name not in self.series_names) for cycle_name in self.cycle_names]
        ):
            self.series_list.append(
                ConstantValue(
                    name=self.cycle_names[0],
                    unit_name=None,
                    value=0,
                )
            )
            self._populate_constants()  # So that everything has a cycle number

    def _populate_constants(self):
        for (i, s) in enumerate(self.series_list):
            if isinstance(s, ConstantValue):
                tseries = self.potential.tseries
                current_series = s.get_vseries(tseries=tseries)
                self.series_list[i] = current_series

    def __getitem__(self, item):
        if item == self.E_str:
            return self.raw_potential
        elif item == self.V_str:
            return self.potential
        elif item == self.I_str:
            return self.raw_current
        elif item == self.J_str:
            return self.current
        elif item == self.sel_str:
            return self.selector
        return super().__getitem__(item)

    @property
    def raw_potential(self):
        if not self._raw_potential:
            self._find_or_build_raw_potential()
        return self._raw_potential

    def _find_or_build_raw_potential(self):
        potential_series_list = [
            s for s in self.series_list if s.name in self.raw_potential_names
        ]
        if len(potential_series_list) == 1:
            self._raw_potential = potential_series_list[0]
        elif len(potential_series_list) > 1:
            raw_potential = append_series(potential_series_list)
            if self._raw_current and self._raw_current.tseries == raw_potential.tseries:
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
            self[self.E_str] = self._raw_potential
        else:
            raise SeriesNotFoundError(
                f"{self} does not have a series corresponding to raw potential."
                f" Looked for series with names in {self.raw_potential_names}"
            )

    @property
    def raw_current(self):
        if not self._raw_current:
            self._find_or_build_raw_current()
        return self._raw_current

    def _find_or_build_raw_current(self):
        current_series_list = [
            s for s in self.series_list if s.name in self.raw_current_names
        ]
        if len(current_series_list) == 1:
            self._raw_current = current_series_list[0]
        elif len(current_series_list) > 1:
            raw_current = append_series(current_series_list)
            if self._raw_potential and (
                self._raw_potential.tseries == raw_current.tseries
            ):
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
            self[self.I_str] = self._raw_current
        else:
            raise SeriesNotFoundError(
                f"{self} does not have a series corresponding to raw current."
                f" Looked for series with names in {self.raw_current_names}"
            )

    def calibrate(self, RE_vs_RHE, A_el=None):
        self.RE_vs_RHE = RE_vs_RHE if RE_vs_RHE is not None else self.RE_vs_RHE
        self.A_el = A_el if A_el is not None else self.A_el

    def normalize(self, A_el, RE_vs_RHE=None):
        self.A_el = A_el if A_el is not None else self.A_el
        self.RE_vs_RHE = RE_vs_RHE if RE_vs_RHE is not None else self.RE_vs_RHE

    @property
    def potential(self):
        raw_potential = self.raw_potential
        if self.RE_vs_RHE is None:
            return raw_potential
        else:
            return ValueSeries(
                name=self.V_str,
                data=raw_potential.data + self.RE_vs_RHE,
                unit_name=raw_potential.unit_name,
                tseries=raw_potential.tseries,
            )

    @property
    def current(self):
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

    def get_potential(self, tspan=None):
        t = self.potential.t.copy()
        v = self.potential.data.copy()
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[mask]
            v = v[mask]
        return t, v

    def get_current(self, tspan=None):
        t = self.current.t.copy()
        j = self.current.data.copy()
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[mask]
            j = j[mask]
        return t, j

    @property
    def t(self):
        return self.potential.t.copy()

    @property
    def v(self):
        return self.potential.data.copy()

    @property
    def j(self):
        return self.current.data.copy()

    @property
    def plotter(self):
        if not self._plotter:
            from ..plotters.ec_plotter import ECPlotter

            self._plotter = ECPlotter(measurement=self)
        return self._plotter

    @property
    def selector(self):
        if self.sel_str not in self.series_names:
            self.build_selector()
        return self._selector

    def build_selector(self, sel_str=None):
        sel_str = sel_str or self.sel_str
        changes = np.tile(False, self.t.shape)
        self["file_number"] = self.file_number
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
        self._selector = selector_series

    @property
    def cycle_number(self):
        cycle_series_list = [s for s in self.series_list if s.name in self.cycle_names]
        if len(cycle_series_list) == 1:
            return cycle_series_list[0]
        elif len(cycle_series_list) > 1:
            cycle = append_series(cycle_series_list)
            return ValueSeries(
                name=self.cycle_str,
                data=cycle.data,
                unit_name=cycle.unit_name,
                tseries=cycle.tseries,
            )
        else:
            raise SeriesNotFoundError(
                f"{self} does not have a series corresponding to raw current."
                f" Looked for series with names in {self.raw_current_names}"
            )

    @property
    def file_number(self):
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
        file_number = append_series(file_number_series_list)
        return file_number
