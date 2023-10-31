from pathlib import Path
from ixdat import Measurement

# Data here: https://www.dropbox.com/scl/fo/z1599uzzi0avmtco67bsd/h?rlkey=ozhkgtq5cfp467j8kavpnoe3c&dl=0

data_dir = Path.home() / "Dropbox/ixdat_resources/debugging/23J31_Issue_93"
# this dataset has an empty file, 02_Rude4A_03_CVA_C01.mpt


meas = Measurement.read_set(data_dir / "02", suffix=".mpt")

meas.plot()

part = meas.cut(tspan=[20, 30])

part.plot()
