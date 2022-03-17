"""Module defining the custom exceptions used in ixdat"""


class DataBaseError(Exception):
    """ixdat errors having with the DataBase"""


class TimeError(Exception):
    """ixdat errors having to do with time"""


class AxisError(Exception):
    """ixdat errors having to do with axes of multi-dimensional data"""


class SeriesNotFoundError(Exception):
    """ixdat errors having to do with a desired data series not available"""


class BuildError(Exception):
    """ixdat errors having to do with series not being able to be appended"""


class ReadError(Exception):
    """ixdat errors having to do with file reading"""


class TechniqueError(Exception):
    """ixdat errors having to do with techniques and their limitations"""


class QuantificationError(Exception):
    """ixdat errors having to do with techniques and their limitations"""


class DeprecationError(Exception):
    """Error raised when a deprecated component has reached hard deprecation"""
