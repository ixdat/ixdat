# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 21:13:50 2024

@author: Søren
"""
import warnings
import numpy as np
from ..measurement_base import Calculator
from ..data_series import ValueSeries, ConstantValue, append_series
from ..exceptions import SeriesNotFoundError, QuantificationError


class Indexer(Calculator):

    calculator_type = "indexer"
    available_series_names = {"file_number", "selector"}

    def __init__(
        self,
        name=None,
        technique=None,
        tstamp=None,
        measurement=None,
        selector_name=None,
        columns=None,
        extra_columns=None,
    ):
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.measurement = measurement
        self.selector_name = selector_name
        self.columns = columns
        self.extra_columns = extra_columns

    def calculate_series(self, key, measurement):
        if key == "file_number":
            return self._build_file_number_series(measurement=measurement)
        elif key == "selector":
            return self._build_selector_series(measurement=measurement)
        raise QuantificationError(f"{self} cannot calculate {key}")

    def _build_file_number_series(self, measurement):
        """Build a `file_number` series based on component measurements times."""
        series_to_append = []
        for i, m in enumerate(measurement.component_measurements or [measurement]):
            if (
                measurement.control_technique_name
                and not m.technique == measurement.control_technique_name
            ):
                continue
            if not measurement.control_series_name:
                tseries = m.time_series[0]
            else:
                try:
                    tseries = m[measurement.control_series_name].tseries
                except SeriesNotFoundError:
                    continue
            series_to_append.append(
                ConstantValue(name="file_number", unit_name="", data=i, tseries=tseries)
            )
        return append_series(
            series_to_append, name="file_number", tstamp=measurement.tstamp
        )

    def _build_selector_series(
        self, measurement, selector_name=None, columns=None, extra_columns=None
    ):
        """Build a `selector` series which demarcates the data.

        The `selector` is a series which can be used to conveniently and powerfully
        grab sections of the data. It is built up from less powerful demarcation series
        in the raw data (like `cycle_number`, `step_number`, `loop_number`, etc) and
        `file_number` by counting the cumulative changes in those series.
        See slide 3 of:
        https://www.dropbox.com/s/sjxzr52fw8yml5k/21E18_DWS3_cont.pptx?dl=0

        Args:
            selector_name (str): The name to use for the selector series
            columns (list): The list of demarcation series. The demarcation series have
                to have equal-length tseries, which should be the one pointed to by the
                meausrement's `control_series_name`.
            extra_columns (list): Extra demarcation series to include if needed.
        """
        # the name of the selector series:
        selector_name = selector_name or self.selector_name or measurement.selector_name
        # a vector that will be True at the points where a series changes:
        changes = np.tile(False, measurement.t.shape)
        # the names of the series which help demarcate the data
        columns = columns or self.columns or measurement.selection_series_names
        extra_columns = extra_columns or self.extra_columns
        if extra_columns:
            columns += extra_columns
        for col in columns:
            try:
                vseries = measurement[col]
            except SeriesNotFoundError:
                continue
            values = vseries.data
            if len(values) == 0:
                warnings.warn("WARNING: " + col + " is empty")
                continue
            elif not len(values) == len(changes):
                warnings.warn("WARNING: " + col + " has an unexpected length")
                continue
            # a vector which is shifted one.
            last_value = np.append(values[0], values[:-1])
            # comparing value and last_value shows where in the vector changes occur:
            changes = np.logical_or(changes, last_value != values)
        # taking the cumsum makes a vector that increases 1 each time one of the
        #   original demarcation vector changes
        selector_data = np.cumsum(changes)
        selector_series = ValueSeries(
            name=selector_name,
            unit_name="",
            data=selector_data,
            tseries=measurement[measurement.control_series_name].tseries,
        )
        return selector_series
