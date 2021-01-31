"""Import readers and build the READER_CLASSES dictionary for direct import

Constants:
    READER_CLASSES (dict): Dictionary of {reader_name: ReaderClass} where
        reader_name is the name of the backend (like "directory") and ReaderClass
        is the reader class for parsing files.
"""
from ..techniques import TECHNIQUE_CLASSES
from .ec_ms_pkl import EC_MS_CONVERTER
from .zilien import ZilienTSVReader
from .biologic import BiologicMPTReader

READER_CLASSES = {
    "EC_MS": EC_MS_CONVERTER,
    "zilien": ZilienTSVReader,
    "biologic": BiologicMPTReader,
}
