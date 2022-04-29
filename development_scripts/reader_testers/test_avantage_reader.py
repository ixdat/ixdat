from pathlib import Path

from ixdat import Spectrum

path_to_file = (
    Path.home() / "Dropbox/ixdat_resources/test_data/avantage" / "XPS Survey.avg"
)

spectrum = Spectrum.read(path_to_file, reader="avantage")

spectrum.plot()
