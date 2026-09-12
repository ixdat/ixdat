"""Sandbox script to aid in development and demo of the Avantage reader"""

from tools_for_demos import DEMO_DATA_DIR
from ixdat import Spectrum

path_to_file = DEMO_DATA_DIR / "avantage/XPS Survey.avg"

spectrum = Spectrum.read(path_to_file, reader="avantage")

spectrum.plot()
