"""This module defines the DataSeries class, the elementary data structure of ixdat

An ixdat DataSeries is a wapper around a numpy array containing the metadata needed
to combine it with other DataSeries. Typically this means a reference to the time
variable corresponding to the rows of the array. The time variable itself is a special
case, TimeSeries, which must know its absolute (unix) time.
"""

from .db import Saveable
from .units import Unit
from .exceptions import TimeError, AxisError


class DataSeries(Saveable):
    """The DataSeries class"""

    table_name = "data_series"
    column_attrs = {"id": "i", "name": "name", "unit": "unit_name", "data": "data"}

    def __init__(self, name, unit, data):
        """initialize a data series"""
        super().__init__()
        self.name = name
        self.unit_name = unit
        self._data = data

    @classmethod
    def from_dict(cls, obj_as_dict):
        if "tstamp" in obj_as_dict:
            return TimeSeries(**obj_as_dict)
        elif "t_id" in obj_as_dict:
            return ValueSeries(**obj_as_dict)
        elif "a_ids" in obj_as_dict:
            return Field(**obj_as_dict)

    def __repr__(self):
        return f"{self.__class__}(id={self.id}, name='{self.name}')"

    @property
    def data(self):
        if self._data is None:
            self._data = self.load_data()
        return self._data

    @property
    def unit(self):
        return Unit(self.unit_name)


class TimeSeries(DataSeries):

    extra_column_attrs = {"tstamps": {"tstamp": "tstamp"}}

    def __init__(self, name, unit, data, tstamp):
        super().__init__(name, unit, data)
        self.tstamp = tstamp


class ValueSeries(DataSeries):

    extra_linkers = {"value_time": ("data_series", {"t_ids": "t_ids"})}

    def __init__(self, name, unit, data, t_id=None, tseries=None):
        super().__init__(name, unit, data)
        self._tseries = tseries
        self._t_id = t_id
        if tseries and t_id:
            if not t_id == tseries.id:
                raise TimeError(f"{self} initiated with non-matching t_id and tseries")
        if tseries is None and t_id is None:
            raise TimeError(f"{self} initiated without t_id or tseries.")

    @property
    def t_id(self):
        if self._tseries:
            return self._tseries.id
        return self._t_id

    @property
    def a_ids(self):
        return [self.t_id]

    @property
    def tseries(self):
        if not self._tseries:
            self._tseries = TimeSeries.open(i=self.t_id)
        return self._tseries

    @property
    def v(self):
        return self.data

    @property
    def t(self):
        return self.tseries.data

    @property
    def tstamp(self):
        return self.tseries.tstamp


class Field(DataSeries):

    extra_linkers = {"field_axes": ("data_series", {"a_ids": "a_ids"})}

    def __init__(self, name, unit, data, a_ids=None, axes_series=None):
        super().__init__(name, unit, data)
        N = len(a_ids) if a_ids is not None else len(axes_series)
        self.N_dimensions = N
        self._a_ids = a_ids if a_ids is not None else ([None] * N)
        self._axes_series = axes_series if axes_series is not None else ([None] * N)
        self._check_axes()  # raises an AxisError if something's wrong

    def get_axis_id(self, axis_number):
        if self._axes_series[axis_number]:
            return self._axes_series[axis_number].id
        return self._a_ids[axis_number]

    def get_axis_series(self, axis_number):
        if not self._axes_series[axis_number]:
            self._axes_series[axis_number] = DataSeries.open(i=self._a_ids[axis_number])
        return self._axes_series[axis_number]

    @property
    def a_ids(self):
        return [self.get_axis_id(n) for n in range(self.N_dimensions)]

    @property
    def axes_series(self):
        return [self.get_axis_series(n) for n in range(self.N_dimensions)]

    def _check_axes(self):
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
        if self._data is None:
            self._data = self.load_data()
            if len(self._data.shape) != self.N_dimensions:
                raise AxisError(
                    f"{self} has {self.N_dimensions} axes but its data is "
                    f"{len(self._data.shape)}-dimensional."
                )
        return self._data
