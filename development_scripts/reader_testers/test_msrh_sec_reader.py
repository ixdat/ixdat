from pathlib import Path
from ixdat.readers.msrh_sec import MsrhSECReader

data_dir = Path(r"C:\Users\scott\Dropbox\ixdat_resources\test_data\sec")

path_to_sec = data_dir / "test-7SEC.csv"
path_to_wl = data_dir / "WL.csv"
path_to_jv = data_dir / "test-7_JV.csv"

sec_meas = MsrhSECReader().read(
    path_to_sec, path_to_wl, path_to_jv, scan_rate=1, tstamp=1
)

sec_meas.plot()
sec_meas.whitelight_spectrum.plot()
print(sec_meas["spectra"].data.shape)
