from pathlib import Path
from ixdat import Measurement

data_dir = Path("~/Dropbox/ixdat_resources/test_data/cinfdata/Krabbe").expanduser()

tpms = Measurement.read(
    data_dir / "baratron_temp_measurement.txt.csv",
    reader="ixdat",
    technique="reactor",
    aliases={"pressure": ["Reactor pressure"], "temperature": ["TC temperature"]},
)

axes = tpms.plot()

axes[0].get_figure().savefig("tpms_plot.png")
