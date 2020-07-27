"""This module defines the DataSeries class, the elementary data structure of ixdat

An ixdat DataSeries is a wapper around a numpy array containing the metadata needed
to combine it with other DataSeries. Typically this means a reference to the time
variable corresponding to the rows of the array. The time variable itself is a special
case, TimeSeries, which must know its absolute (unix) time.
"""


class DataSeries:
    """The DataSeries class"""

    def __init__(self):
        """initialize a data series"""
        pass
