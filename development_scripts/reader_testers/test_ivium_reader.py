"""For use in development of the ivium reader. Requires access to sample data."""

from pathlib import Path
import pandas as pd

from ixdat import Measurement
from ixdat.techniques import CyclicVoltammogram

path_to_file = Path.home() / (
    "Dropbox/ixdat_resources/test_data/ivium/ivium_test_dataset"
)
path_to_single_file = path_to_file.parent / (path_to_file.name + "_1")
df = pd.read_csv(path_to_single_file, sep=r"\s+", header=1)

meas = Measurement.read(path_to_file, reader="ivium")

meas.save()

meas.plot_measurement()

meas_cv = CyclicVoltammogram.read(path_to_file, reader="ivium")

meas_cv.save()

meas_cv.plot_measurement()
meas_cv.redefine_cycle(start_potential=0.4, redox=False)
for i in range(4):
    meas_cv[i].plot()
