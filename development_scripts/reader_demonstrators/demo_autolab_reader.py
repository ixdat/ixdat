"""For use in development of the autolab reader. Requires access to sample data."""

from tools_for_demos import DEMO_DATA_DIR
from ixdat import Measurement

path_to_file = DEMO_DATA_DIR / "autolab/autolab_test_file.txt"

meas = Measurement.read(path_to_file, reader="autolab")

meas.plot()
