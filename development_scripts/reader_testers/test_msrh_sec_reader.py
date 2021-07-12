"""Demonstrate simple importing and plotting SEC data"""

from pathlib import Path
from ixdat import Measurement

from matplotlib import pyplot as plt

plt.close("all")


data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/sec"

sec_meas = Measurement.read(
    data_dir / "test-7SEC.csv",
    path_to_ref_spec_file=data_dir / "WL.csv",
    path_to_V_J_file=data_dir / "test-7_JV.csv",
    scan_rate=1,
    tstamp=1,
    reader="msrh_sec",
)

axes = sec_meas.plot_measurement(
    V_ref=0.4,
    cmap_name="jet",
    make_colorbar=True,
)

ax = sec_meas.plot_waterfall(
    V_ref=0.4,
    cmap_name="jet",
    make_colorbar=True,
)

axes2 = sec_meas.plot_vs_potential(V_ref=0.66, cmap_name="jet", make_colorbar=False)
axes2 = sec_meas.plot_vs_potential(
    V_ref=0.4, vspan=[0.5, 2], cmap_name="jet", make_colorbar=False
)

ax = sec_meas.get_dOD_spectrum(V_ref=0.4, V=0.8).plot(color="b", label="species 1")
sec_meas.get_dOD_spectrum(V_ref=0.8, V=1.2).plot(color="g", ax=ax, label="species 2")
sec_meas.get_dOD_spectrum(V_ref=1.2, V=1.49).plot(color="r", ax=ax, label="species 3")
ax.legend()
