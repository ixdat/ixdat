from tools_for_demos import DEMO_DATA_DIR
from ixdat import Measurement

data_dir = DEMO_DATA_DIR / "cinfdata/Krabbe"

tpms = Measurement.read(
    data_dir / "baratron_temp_measurement.txt.csv",
    reader="ixdat",
    technique="reactor",
    aliases={"pressure": ["Reactor pressure"], "temperature": ["TC temperature"]},
)

axes = tpms.plot()

axes[0].get_figure().savefig("tpms_plot.png")
