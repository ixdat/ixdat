from tools_for_demos import DEMO_DATA_DIR
from ixdat import Measurement

data_dir = DEMO_DATA_DIR / "nordic_tdms/24B07_0_Pt"

c1 = Measurement.read(data_dir / "CV_101448_ 1.tdms", reader="nordic")

c1.plot()

meas = Measurement.read_set(data_dir, reader="nordic", suffix=".tdms")

meas.plot()

cv = meas.as_cv()
cv.redefine_cycle(start_potential=0.4, redox=True)

if False:
    # FIXME AttributeError: module 'matplotlib.cm' has no attribute 'get_cmap'
    cv[15:30].plot_cycles()
