"""Sandbox script to aid in development and demo of the XRDML reader"""

from tools_for_demos import DEMO_DATA_DIR

from ixdat import Spectrum

path_to_file = DEMO_DATA_DIR / "xrdml/GI-XRD Path 2_1 omega 0p5 step 10s.xrdml"


spectrum = Spectrum.read(path_to_file, reader="xrdml")

spectrum.plot()
