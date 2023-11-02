from pathlib import Path
from ixdat import Measurement
from matplotlib import pyplot as plt


data_dir = Path("../data/02")

sec_meas = Measurement.read(
    data_dir / "test-7SEC.csv",
    path_to_ref_spec_file=data_dir / "WL.csv",
    path_to_U_J_file=data_dir / "test-7_JV.csv",
    scan_rate=1,
    tstamp=1,
    reader="msrh_sec",
)

sec_meas.calibrate_RE(RE_vs_RHE=0.2)
sec_meas.set_reference_spectrum(V_ref=0.6)

ax = sec_meas.plot_waterfall()

r =  0.85
fig = ax.get_figure()
fig.set_figwidth(6 * r)
fig.set_figheight(4.5 * r)
fig.tight_layout()
fig.savefig("demo_02.png", dpi=600)

axes2 = sec_meas.plot_measurement()
axes2[0].set_xlabel("time / [s]")
fig = axes2[0].get_figure()
fig.set_figwidth(6 * r)
fig.set_figheight(4.5 * r)
fig.tight_layout()
fig.savefig("demo_02_measurement.png", dpi=600)

r = 0.75
ax3 = sec_meas.get_dOD_spectrum(t=300).plot()
fig = ax3.get_figure()
fig.set_figwidth(6 * r)
fig.set_figheight(4.5 * r)
fig.tight_layout()
fig.savefig("demo_02_spectrum.png", dpi=600)

r = 0.85
sec_meas.plot_wavelengths_vs_potential(wavelengths=["w800", "w580", "w480"])
fig = plt.gcf()  # FIXME: plot_wavelengths_vs_potential doesn't return an axis!
fig.set_figwidth(6 * r)
fig.set_figheight(4.5 * r)
fig.tight_layout()
fig.savefig("demo_02_wavelengths.png", dpi=600)