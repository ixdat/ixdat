"""Configuration variables including stnadard data directory

constants:
    STANDARD_DATA_DIRECTORY (Path): the directory in which to save by default with
        the default database backend, which is to save files in a directory.
        ixdat will make the directory if it does not exist.
    STANDARD_SUFFIX (str): The file extension for JSON format metadata files
    STANDARD_DATA_SUFFIX (str): The file extension for numpy.save format data files
"""
from pathlib import Path

STANDARD_METADATA_SUFFIX = "ix"
STANDARD_DATA_SUFFIX = "ixdata"
STANDARD_DATA_DIRECTORY = Path.home() / "ixdat"

if not STANDARD_DATA_DIRECTORY.exists():
    STANDARD_DATA_DIRECTORY.mkdir()
