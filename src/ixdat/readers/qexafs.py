"""Readers for 'qexafs' files exported by Diamond's B18-Core"""

from pathlib import Path
from .. import Spectrum


class QexafsDATReader:
    def read(self, path_to_file, cls=Spectrum):
        with open(path_to_file, "r"):
            pass
        return path_to_file
