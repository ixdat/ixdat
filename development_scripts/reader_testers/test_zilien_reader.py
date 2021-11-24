from pathlib import Path

from ixdat import Measurement
from ixdat.techniques import MSMeasurement

data_dir = Path(r"C:\Users\scott\Dropbox\ixdat_resources\test_data\zilien_with_ec")

path_to_file = data_dir / "2021-02-01 17_44_12.tsv"

# This imports it with the EC data
ecms = Measurement.read(path_to_file, reader="zilien")
ecms.plot_measurement()

# This imports it without the EC data:
ms = MSMeasurement.read(path_to_file, reader="zilien")
ms.plot_measurement()  # nice. one panel, no MS :)

# This adds in the EC data from Biologic:

ec = Measurement.read_set(
    data_dir / "2021-02-01 17_44_12", reader="biologic", suffix=".mpt"
)
ecms_2 = ec + ms
ecms_2.plot_measurement()
