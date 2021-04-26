from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt

from ixdat.techniques import CyclicVoltammagram

path_to_file = Path.home() / (
    "Dropbox/ixdat_resources/test_data/ivium/ivium_test_dataset"
)
path_to_single_file = path_to_file.parent / (path_to_file.name + "_1")
df = pd.read_csv(path_to_single_file, sep=r"\s+", header=1)

meas = CyclicVoltammagram.read(path_to_file, reader="ivium")

meas.save()

meas.plot_measurement()
meas.redefine_cycle(start_potential=0.4, redox=False)
for i in range(4):
    meas[i].plot()
