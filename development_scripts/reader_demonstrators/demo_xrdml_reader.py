"""Sandbox script to aid in development and demo of the XRDML reader"""

from pathlib import Path

from ixdat import Spectrum

path_to_file = (
    Path.home()
    / "Dropbox/ixdat_resources/test_data/xrdml"
    / "GI-XRD Path 2_1 omega 0p5 step 10s.xrdml"
)

spectrum = Spectrum.read(path_to_file, reader="xrdml")

spectrum.plot()
