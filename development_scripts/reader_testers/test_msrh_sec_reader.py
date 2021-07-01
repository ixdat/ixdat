from pathlib import Path
from ixdat import Measurement

data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/sec"

path_to_sec = data_dir / "test-7SEC.csv"
path_to_wl = data_dir / "WL.csv"
path_to_jv = data_dir / "test-7_JV.csv"

sec_meas = Measurement.read(
    path_to_sec,
    path_to_wl_file=path_to_wl,
    path_to_jv_file=path_to_jv,
    scan_rate=1,
    tstamp=1,
    reader="msrh_sec",
)

axes = sec_meas.plot_measurement(V_ref=0.66, cmap_name="jet", make_colorbar=True)
axes[0].get_figure().savefig("sec_example_with_colorbar.png")
ax = sec_meas.plot_waterfall(V_ref=0.66, cmap_name="jet", make_colorbar=True)
ax.get_figure().savefig("sec_waterfall_example.png")
