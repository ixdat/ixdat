"""Import readers and build the READER_CLASSES dictionary for direct import

Constants:
    READER_CLASSES (dict): Dictionary of {reader_name: ReaderClass} where
        reader_name is the name of the backend (like "directory") and ReaderClass
        is the reader class for parsing files.
"""
from ..techniques import TECHNIQUE_CLASSES

# ixdat
from .ixdat_csv import IxdatCSVReader

# potentiostats
from .biologic import BiologicMPTReader
from .autolab import AutolabTXTReader

# ec-ms
from .zilien import ZilienTSVReader
from .ec_ms_pkl import EC_MS_CONVERTER

READER_CLASSES = {
    "ixdat": IxdatCSVReader,
    "biologic": BiologicMPTReader,
    "autolab": AutolabTXTReader,
    "zilien": ZilienTSVReader,
    "EC_MS": EC_MS_CONVERTER,
}
