"""Module for representation and analysis of EC measurements"""

import numpy as np

from ..measurements import Measurement, append_series
from ..data_series import ValueSeries
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

        self._raw_potential = None
        self._raw_current = None

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
            self.series_list.append(self._raw_potential)
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
            self.series_list.append(self._raw_current)

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
