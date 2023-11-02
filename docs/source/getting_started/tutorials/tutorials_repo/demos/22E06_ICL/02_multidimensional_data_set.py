from pathlib import Path
from ixdat import Measurement

data_dir = Path("../../data/02")

sec_meas = Measurement.read(
    data_dir / "test-7SEC.csv",
    path_to_ref_spec_file=data_dir / "WL.csv",
    path_to_U_J_file=data_dir / "test-7_JV.csv",
    scan_rate=1,
    tstamp=1,
    reader="msrh_sec",
)

sec_meas.plot(cmap_name="jet")
sec_meas.calibrate_RE(RE_vs_RHE=0.2)
sec_meas.set_reference_spectrum(V_ref=0.6)

axes = sec_meas.plot_vs_potential(cmap_name="jet")
axes[1].set_yscale("log")

spectrum_1 = sec_meas.get_dOD_spectrum(V=1.3, V_ref=0.6)
spectrum_2 = sec_meas.get_dOD_spectrum(V=1.5, V_ref=1.3)

ax = spectrum_1.plot(color="k", label="before onset")
spectrum_2.plot(color="r", label="change around onset", ax=ax)
ax.set_ylim(bottom=0)

sec_meas.plot_waterfall()
