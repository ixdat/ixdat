"""For use in development of the zilien reader(s). Requires access to sample data."""

from tools_for_demos import DEMO_DATA_DIR
from ixdat import Measurement
from ixdat.techniques import MSMeasurement, ECMeasurement

data_dir = DEMO_DATA_DIR / "zilien_with_ec"

path_to_file = data_dir / "2021-02-01 17_44_12.tsv"

# This imports it with the EC data
ecms = Measurement.read(path_to_file, reader="zilien")
ecms.calibrate_RE(0)
ecms.plot_measurement()

# This imports it as just an MS measurement.
ms = MSMeasurement.read(path_to_file, reader="zilien")
ms.plot_measurement()  # nice. one panel, no MS :)

# This adds in the EC data from Biologic:

ec = Measurement.read_set(
    data_dir / "2021-02-01 17_44_12", reader="biologic", suffix=".mpt"
)
ecms_2 = ec + ms
ecms_2.plot_measurement()

# This imports it as just an EC measurement
ec_2 = ECMeasurement.read(path_to_file, reader="zilien")
ec_2.plot()


# This plots just the EC data from the first EC-MS measurement:

ecms.ec_plotter.plot_measurement()
ecms.ec_plotter.plot_vs_potential()
ecms.ms_plotter.plot_measurement()

# This plots it as a cyclic voltammagram
if False:
    # FIXME: Measurement.__init__() got an unexpected keyword argument 'spectrum_id'
    ecms_cv = ecms.as_cv()
    ecms_cv.ec_plotter.plot_vs_potential()
    ecms_cv.plot()
