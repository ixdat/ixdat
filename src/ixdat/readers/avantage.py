from pathlib import Path
from ..spectra import Spectrum

class AvantageAVGReader:
    """A class for importing a .avg file exported by DataSpace_BatchDump.exe"""
    def __init__(self, path_to_file=None):
        self.path_to_file = path_to_file

    def read(self, path_to_file, cls=Spectrum):
        path_to_file = Path(path_to_file or self.path_to_file)

        # This is just for testing. When finished it'll return a Spectrum object.
        return path_to_file