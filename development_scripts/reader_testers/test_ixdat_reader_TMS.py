from pathlib import Path
from ixdat import Measurement

data_dir = Path("~/Dropbox/ixdat_resources/test_data/cinfdata/Krabbe").expanduser()

tms = Measurement.read(
    data_dir / "baratron_temp_measurement.txt.csv",
    reader="ixdat",
    technique="reactor",
)

tms.plot()
