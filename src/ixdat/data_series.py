"""This module defines the DataSeries class, the elementary data structure of ixdat

An ixdat DataSeries is a wrapper around a numpy array containing the metadata needed
to combine it with other DataSeries. Typically this means a reference to the time
variable corresponding to the rows of the array. The time variable itself is a special
case, TimeSeries, which must know its absolute (unix) timestamp.
"""

import numpy as np
from .db import Saveable
from .units import Unit
from .exceptions import TimeError, AxisError


class DataSeries(Saveable):
    """The base class for all numerical data representation in ixdat.

    These class's objects are saved and loaded as rows in the data_series table
    """

    table_name = "data_series"
    column_attrs = {
        "name",
        "unit_name",
        "data",
    }

    def __init__(self, name, unit_name, data):
        """initialize a data series with its name, unit, and data (id handled by parent)

        Args:
            name (str): The name of the data series
            unit_name (str): The name of the unit in which the data is stored
            data (np.array): The numerical data
        """
        super().__init__()
        self.name = name
        self.unit = Unit(unit_name)
        self._data = data

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return the right type of DataSeries based on the info in its serialization"""
        if "tstamp" in obj_as_dict:
            return TimeSeries(**obj_as_dict)
        elif "t_ids" in obj_as_dict:
            return ValueSeries(**obj_as_dict)
        elif "a_ids" in obj_as_dict:
            return Field(**obj_as_dict)
        elif "value" in obj_as_dict:
            return ConstantValue(**obj_as_dict)
        return cls(**obj_as_dict)

    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.id}, name='{self.name}')"

    @property
    def data(self):
        """The data as a np.array, loaded the first time it is needed."""
        if self._data is None:
            self._data = self.load_data()  # inherited from Saveable.
        return self._data

    @property
    def unit_name(self):
        """The name of the data series' unit"""
        return self.unit.name


class TimeSeries(DataSeries):
    """Class to store time data. These are characterized by having a tstamp"""

    extra_column_attrs = {"tstamps": {"tstamp"}}

    def __init__(self, name, unit_name, data, tstamp):
        """Initiate a TimeSeries with name, unit_name, data, and a tstamp (float)

        Args (in addition to those of parent):
            tstamp (float): The unix timestamp of the time at which t=0 in the data
        """
        super().__init__(name, unit_name, data)
        self.tstamp = tstamp

    @property
    def tseries(self):
        """Trivially, a TimeSeries is its own TimeSeries"""
        return self


class ValueSeries(DataSeries):
    """Class to store scalar values that are measured over time.

    Characterized by a reference to the corresponding time series. This reference is
    represented in relational databases as a row in an auxiliary linker table
    """

    extra_linkers = {"value_time": ("data_series", "t_ids")}

    def __init__(self, name, unit_name, data, t_id=None, t_ids=None, tseries=None):
        """Initiate a ValueSeries with a TimeSeries or a reference thereto

        Args (in addition to those of parent):
            t_id (int): The id of the corresponding TimeSeries, if not given directly
            t_ids (list of int): [t_id], only so that a backend can pass t_id as a list
            tseries (TimeSeries): The corresponding TimeSeries, if available
        """
        super().__init__(name, unit_name, data)
        self._tseries = tseries
        # TODO: This could probably be handled more nicely with PlaceHolderObjects
        #   see: Measurement and
        #   https://github.com/ixdat/ixdat/pull/1#discussion_r551518461
        if t_ids and not t_id:
            t_id = t_ids[0]
        self._t_id = t_id
        if tseries and t_id:
            if not t_id == tseries.id:
                raise TimeError(f"{self} initiated with non-matching t_id and tseries")
        if tseries is None and t_id is None:
            raise TimeError(f"{self} initiated without t_id or tseries.")

    @property
    def t_id(self):
        """int: the id of the TimeSeries"""
        if self._tseries:
            return self._tseries.id
        return self._t_id

    @property
    def t_ids(self):
        """list: the id of the TimeSeries, in a list for consistent linker table def."""
        return [self.t_id]

    @property
    def tseries(self):
        """The TimeSeries describing when the data in the ValueSeries was recorded"""
        if not self._tseries:
            self._tseries = TimeSeries.get(i=self.t_id)
            self._t_id = None  # to avoid any confusion of two t_id's
        return self._tseries

    @property
    def v(self):
        """The value as a 1-d np array"""
        return self.data

    @property
    def t(self):
        """The measurement times as a 1-d np array"""
        return self.tseries.data

    @property
    def tstamp(self):
        """The timestamp, from the TimeSeries of the ValueSeries"""
        return self.tseries.tstamp


class Field(DataSeries):
    """Class for storing multi-dimensional data spanning 'axes'

    Characterized by a list of references to these axes, which are themselves also
    DataSeries. This is represented in the extra linkers.
    """

    extra_linkers = {"field_axes": ("data_series", "a_ids")}

    def __init__(self, name, unit_name, data, a_ids=None, axes_series=None):
        """Initiate the Field and check that the supplied axes make sense.

        Args (in addition to those of parent):
            a_ids (list of int): The ids of the corresponding axes DataSeries, if not
                the series are not given directly as `axes_series`
            axes_series (list of DataSeries): The DataSeries describing the axes which
                the field's data spans, if available
        """
        super().__init__(name, unit_name, data)
        N = len(a_ids) if a_ids is not None else len(axes_series)
        self.N_dimensions = N
        self._a_ids = a_ids if a_ids is not None else ([None] * N)
        # TODO: This could probably be handled more nicely with PlaceHolderObjects
        #   see: Measurement and
        #   https://github.com/ixdat/ixdat/pull/1#discussion_r551518461
        self._axes_series = axes_series if axes_series is not None else ([None] * N)
        self._check_axes()  # raises an AxisError if something's wrong

    def get_axis_id(self, axis_number):
        """Return the id of the `axis_number`'th axis of the data"""
        if self._axes_series[axis_number]:
            return self._axes_series[axis_number].id
        return self._a_ids[axis_number]

    def get_axis_series(self, axis_number):
        """Return the DataSeries of the `axis_number`'th axis of the data"""
        if not self._axes_series[axis_number]:
            self._axes_series[axis_number] = DataSeries.get(i=self._a_ids[axis_number])
            # And so as not have two id's for the axis_number'th axis:
            self._a_ids[axis_number] = None
        return self._axes_series[axis_number]

    @property
    def a_ids(self):
        """List of the id's of the axes spanned by the field"""
        return [self.get_axis_id(n) for n in range(self.N_dimensions)]

    @property
    def axes_series(self):
        """List of the DataSeries defining the axes spanned by the field"""
        return [self.get_axis_series(n) for n in range(self.N_dimensions)]

    def _check_axes(self):
        """Check that there are no contradictions in the Field's axes_series and id's"""
        N = self.N_dimensions
        if len(self._a_ids) != N:
            raise AxisError(
                f"{self} is {N}-D but initiated with {len(self._a_ids)} axis id's"
            )
        if len(self._axes_series) != N:
            raise AxisError(
                f"{self} is {N}-D but initiated with {len(self._axes_series)} axes"
            )
        for n, (a_id, axis_series) in enumerate(zip(self._a_ids, self._axes_series)):
            if a_id is not None and axis_series is not None and a_id != axis_series.id:
                raise AxisError(
                    f"{self} initiated with contradicting id's for {n}'th axis"
                )
            elif a_id is None and axis_series is None:
                raise AxisError(
                    f"{self} has no axis id for series or id for its {n}'th axis"
                )

    @property
    def data(self):
        """When loading data, Field checks that its dimensions match its # of axes"""
        if self._data is None:
            self._data = self.load_data()
            if len(self._data.shape) != self.N_dimensions:
                raise AxisError(
                    f"{self} has {self.N_dimensions} axes but its data is "
                    f"{len(self._data.shape)}-dimensional."
                )
        return self._data


class ConstantValue(DataSeries):
    """This is a stand-in for a VSeries for when we know the value is constant"""

    extra_column_attrs = {"constants": {"value"}}

    def __init__(self, name, unit_name, data=None, value=None):
        super().__init__(name=name, unit_name=unit_name, data=np.array([]))
        if not np.array(value).size == 1:
            raise AxisError(
                f"Can't initiate {self} with data={self.value}. Data must have size 1."
            )
        self.value = value

    def get_vseries(self, tseries):
        data = self.value * np.ones(tseries.data.shape)
        return ValueSeries(
            name=self.name, unit_name=self.unit_name, data=data, tseries=tseries
        )
