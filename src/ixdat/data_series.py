"""This module defines the DataSeries class, the elementary data structure of ixdat

An ixdat DataSeries is a wrapper around a numpy array containing the metadata needed
to combine it with other DataSeries. Typically this means a reference to the time
variable corresponding to the rows of the array. The time variable itself is a special
case, TimeSeries, which must know its absolute (unix) timestamp.
"""

import numpy as np
from .db import Saveable
from .tools import tstamp_to_string
from .units import Unit
from .exceptions import AxisError, BuildError


class DataSeries(Saveable):
    """The base class for all numerical data representation in ixdat.

    These class's objects are saved and loaded as rows in the data_series table
    """

    table_name = "data_series"
    column_attrs = {"name", "unit_name", "data", "series_type"}
    series_type = "series"

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
        series_type = obj_as_dict.pop("series_type")
        series_class = SERIES_CLASSES[series_type]
        return series_class(**obj_as_dict)

    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.id}, name='{self.name}')"

    @property
    def data(self):
        """The data as a np.array, loaded the first time it is needed."""
        if self._data is None:
            self._data = self.load_data()  # inherited from Savable.
        return self._data

    @property
    def unit_name(self):
        """The name of the data series' unit"""
        return self.unit.name

    @property
    def shape(self):
        return self.data.shape

    @property
    def size(self):
        return self.data.size


class TimeSeries(DataSeries):
    """Class to store time data. These are characterized by having a tstamp"""

    extra_column_attrs = {"tstamps": {"tstamp"}}
    series_type = "tseries"

    def __init__(self, name, unit_name, data, tstamp):
        """Initiate a TimeSeries with name, unit_name, data, and a tstamp (float)

        Args (in addition to those of parent, :class:`.DataSeries`):
            tstamp (float): The unix timestamp of the time at which t=0 in the data
        """
        super().__init__(name, unit_name, data)
        self.tstamp = tstamp

    def __str__(self):
        """Return TimeSeries string representation"""
        # On the form: TimeSeries: 'NAME'. Min, max: 12, 4000 [s] @ 22E18 14:34:55
        return (
            f"{self.__class__.__name__}: '{self.name}'. "
            f"Min, max: {min(self.data):.0f}, {max(self.data):.0f} [{self.unit.name}] "
            f"@ {tstamp_to_string(self.tstamp)}"
        )

    @property
    def t(self):
        return self.data

    @property
    def tseries(self):
        """Trivially, a TimeSeries is its own TimeSeries"""
        return self


class Field(DataSeries):
    """Class for storing multi-dimensional data spanning 'axes'

    Characterized by a list of references to these axes, which are themselves also
    DataSeries. This is represented in the extra linkers.
    """

    extra_linkers = {"field_axes": ("data_series", "a_ids")}
    child_attrs = ["axes_series"]
    series_type = "field"

    def __init__(self, name, unit_name, data, a_ids=None, axes_series=None):
        """Initiate the Field and check that the supplied axes make sense.

        Args (in addition to those of parent, :class:`.DataSeries`):
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
                f"{self!r} is {N}-D but initiated with {len(self._a_ids)} axis id's"
            )
        if len(self._axes_series) != N:
            raise AxisError(
                f"{self!r} is {N}-D but initiated with {len(self._axes_series)} axes"
            )
        for n, (a_id, axis_series) in enumerate(zip(self._a_ids, self._axes_series)):
            if a_id is not None and axis_series is not None and a_id != axis_series.id:
                raise AxisError(
                    f"{self!r} initiated with contradicting id's for {n}'th axis"
                )
            elif a_id is None and axis_series is None:
                raise AxisError(
                    f"{self!r} has no axis id for series or id for its {n}'th axis"
                )

    @property
    def data(self):
        """When loading data, Field checks that its dimensions match its # of axes"""
        if self._data is None:
            self._data = self.load_data()
            if len(self._data.shape) != self.N_dimensions:
                raise AxisError(
                    f"{self!r} has {self.N_dimensions} axes but its data is "
                    f"{len(self._data.shape)}-dimensional."
                )
        return self._data.copy()  # TODO: make data series data immutable with numpy flag
        # see: https://github.com/ixdat/ixdat/pull/101/files#r1126172936

    @property
    def tstamp(self):
        """The unix time corresponding to t=0 for the time-resolved axis of the Field

        The timestamp of a Field is the timestamp of its TimeSeries or ValueSeries
        """
        for s in self.axes_series:
            if isinstance(s, (ValueSeries, TimeSeries)):
                return s.tstamp


class ValueSeries(Field):
    """Class to store scalar values that are measured over time.

    Characterized by a reference to the corresponding time series. This reference is
    represented in relational databases as a row in an auxiliary linker table
    """

    series_type = "vseries"

    def __init__(
        self,
        name,
        unit_name,
        data,
        t_id=None,
        tseries=None,
        a_ids=None,
        axes_series=None,
    ):
        """Initiate a ValueSeries with a TimeSeries or a reference thereto

        Args (in addition to those of :class:`.Field`):
            t_id (int): The id of the corresponding TimeSeries, if not given directly
                (can also be supplied as `a_ids[0]`)
            tseries (TimeSeries): The corresponding TimeSeries, if available
                (can also be supplied as `axes_series[0]`)
        """
        a_ids = a_ids or [t_id]
        axes_series = axes_series or [tseries]
        super().__init__(name, unit_name, data, a_ids, axes_series)
        # TODO: This could probably be handled more nicely with PlaceHolderObjects
        #   see: Measurement and
        #   https://github.com/ixdat/ixdat/pull/1#discussion_r551518461

    def __str__(self):
        """Return string representation"""
        # Return a string representation on the form:
        # ValueSeries: 'NAME'. Min, max: -1.23, 4.56 [V]
        return (
            f"{self.__class__.__name__}: '{self.name}'. "
            f"Min, max: {min(self.data):.1e}, {max(self.data):.1e} [{self.unit.name}]"
        )

    @property
    def tseries(self):
        return self.axes_series[0]

    @property
    def t_id(self):
        """int: the id of the TimeSeries"""
        if self._axes_seriess:
            return self.tseries.id
        return self.a_ids[0]

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

    def __hash__(self):
        return super().__hash__()


class ConstantValue(ValueSeries):
    """This is a stand-in for a VSeries for when we know the value is constant"""

    series_type = "constantvalue"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._expanded_data = None

    @property
    def data(self):
        if self._expanded_data is None:
            if self._data is None:
                self._data = self.load_data()  # inherited from Savable.
            self._expanded_data = np.ones(self.t.shape) * self._data
        return self._expanded_data


SERIES_CLASSES = {
    cls.series_type: cls
    for cls in [DataSeries, TimeSeries, Field, ValueSeries, ConstantValue]
}


def append_series(series_list, sorted=True, name=None, tstamp=None):
    """Return series appending series_list relative to series_list[0].tseries.tstamp

    Args:
        series_list (list of Series): The series to append (must all be of same type)
        sorted (bool): Whether to sort the data so that time only goes forward
        name (str): Name to give the appended series. Defaults to series_list[0].name
        tstamp (unix tstamp): The t=0 of the returned series or its TimeSeries.
    """
    s0 = series_list[0]
    if isinstance(s0, TimeSeries):
        return append_tseries(series_list, sorted=sorted, name=name, tstamp=tstamp)
    elif isinstance(s0, ValueSeries):
        return append_vseries_by_time(
            series_list, sorted=sorted, name=name, tstamp=tstamp
        )
    raise BuildError(
        f"An algorithm of append_series for series like {s0!r} is not yet implemented"
    )


def append_vseries_by_time(series_list, sorted=True, name=None, tstamp=None):
    """Return new ValueSeries with the data in series_list appended

    Args:
        series_list (list of ValueSeries): The value series to append
        sorted (bool): Whether to sort the data so that time only goes forward
        name (str): Name to give the appended series. Defaults to series_list[0].name
        tstamp (unix tstamp): The t=0 of the returned ValueSeries' TimeSeries.
    """
    name = name or series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    data = np.array([])
    tseries_list = [s.tseries for s in series_list]
    if not all(isinstance(ts, TimeSeries) for ts in tseries_list):
        raise BuildError(
            f"can't append {series_list} w incompatible tseries list = {tseries_list}"
        )
    tseries, sort_indeces = append_tseries(
        tseries_list, sorted=sorted, return_sort_indeces=True, tstamp=tstamp
    )

    for s in series_list:
        data = np.append(data, s.data)
    if sorted:
        data = data[sort_indeces]

    return cls(name=name, unit_name=unit.name, data=data, tseries=tseries)


def append_tseries(
    series_list, sorted=True, return_sort_indeces=False, name=None, tstamp=None
):
    """Return new TimeSeries with the data appended.

    Args:
        series_list (list of TimeSeries): The time series to append
        sorted (bool): Whether to sort the data so that time only goes forward
        return_sort_indeces (bool): Whether to return the indeces that sort the data
        name (str): Name to give the appended series. Defaults to series_list[0].name
        tstamp (unix tstamp): The t=0 of the returned TimeSeries.
    """
    name = name or series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    tstamp = tstamp or series_list[0].tstamp
    data = np.array([])

    for s in series_list:
        if not (s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        data = np.append(data, s.data + s.tstamp - tstamp)

    if sorted:
        sort_indices = np.argsort(data)
        data = data[sort_indices]
    else:
        sort_indices = None

    tseries = cls(name=name, unit_name=unit.name, data=data, tstamp=tstamp)
    if return_sort_indeces:
        return tseries, sort_indices
    return tseries


def time_shifted(series, tstamp=None):
    """Return a series with the time shifted to be relative to tstamp"""
    if tstamp is None:
        return series
    if tstamp == series.tstamp:
        return series
    cls = series.__class__
    if isinstance(series, TimeSeries):
        new_data = series.data + series.tstamp - tstamp  # shift the time.
        return cls(
            name=series.name,
            unit_name=series.unit.name,
            data=new_data,
            tstamp=tstamp,
        )
    elif isinstance(series, ValueSeries):
        series = cls(
            name=series.name,
            unit_name=series.unit.name,
            data=series.data,
            tseries=time_shifted(series.tseries, tstamp=tstamp),
        )
    return series


def get_tspans_from_mask(t, mask):
    """Return a list of tspans for time intervals remaining when mask is applied to t

    FIXME: This is pure numpy manipulation and probably belongs somewhere else.
    """
    mask_prev = np.append(False, mask[:-1])  # the mask shifted right by one
    mask_next = np.append(mask[1:], False)  # the mask shifted left by one
    # An array that is True where intervals meeting the criteria start:
    #   (This includes at [0] if mask[0] is True.)
    interval_starts_here = np.logical_and(np.logical_not(mask_prev), mask)
    # An array that is True where intervals meeting the criteria finish:
    #   (This includes at [-1] if mask[-1] is True.)
    interval_ends_here = np.logical_and(mask, np.logical_not(mask_next))

    t_starts = list(t[interval_starts_here])  # the start times implied by the mask
    t_ends = list(t[interval_ends_here])  # the finish times implied by the mask
    tspans = zip(t_starts, t_ends)  # and, the timespans where the criteria is met!
    return tspans
