from pathlib import Path
from ixdat import Measurement


data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/nordic_tdms/24B07_0_Pt"

c1 = Measurement.read_set(data_dir / "CV_101448_ 1.tdms", reader="nordic", suffix=".tdms")

meas = Measurement.read_set(data_dir, reader="nordic", suffix=".tdms")

meas.plot()
