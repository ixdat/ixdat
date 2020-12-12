"""Module defining the custom exceptions used in ixdat"""


class DataBaseError(Exception):
    """ixdat errors having with the DataBase"""


class TimeError(Exception):
    """ixdat errors having to do with time"""


class AxisError(Exception):
    """ixdat errors having to do with axes of multi-dimensional data"""
