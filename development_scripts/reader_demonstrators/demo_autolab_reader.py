"""For use in development of the autolab reader. Requires access to sample data."""

from pathlib import Path

from ixdat import Measurement

path_to_file = Path.home() / (
    "Dropbox/ixdat_resources/test_data/autolab/autolab_test_file.txt"
)

meas = Measurement.read(path_to_file, reader="autolab")

meas.plot()
