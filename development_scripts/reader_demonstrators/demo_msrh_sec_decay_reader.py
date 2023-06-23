"""For use in development of the MSRH SEC reader. Requires access to sample data.
MSRH = molecular science research hub, at Imperial College London.
"""

from pathlib import Path
from ixdat import Measurement

data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/sec"

sec_meas = Measurement.read(
    # data_dir / "decay/PDtest-1.35-1OSP-SP.csv",
    data_dir / "decay/PDtest-1.33-1OSP-SP.csv",
    path_to_ref_spec_file=data_dir / "WL.csv",
    # path_to_t_V_file=data_dir / "decay/PDtest-1.35-1OSP-E-t.csv",
    # path_to_t_J_file=data_dir / "decay/PDtest-1.35-1OSP-J-t.csv",
    path_to_t_U_file=data_dir / "decay/PDtest-1.33-1OSP-E-t.csv",
    path_to_t_J_file=data_dir / "decay/PDtest-1.33-1OSP-J-t.csv",
    tstamp=1,
    reader="msrh_sec_decay",
)
# Suggestion: command-line switching for development scripts.
#  https://github.com/ixdat/ixdat/pull/30/files#r810014299

sec_meas.calibrate_RE(RE_vs_RHE=0.26)

sec_meas.set_reference_spectrum(t_ref=5)

axes = sec_meas.plot_measurement(
    # V_ref=0.66,  # can't do a V_ref for this as can't interpolate on potential..
    # So OD will be calculated using the reference spectrum in WL.csv
    # cmap_name="jet",
    cmap_name="inferno",
    make_colorbar=False,
)
# axes[0].get_figure().savefig("decay_vs_t.png")

axes = sec_meas.plot_wavelengths(wavelengths=["w500", "w600", "w700", "w800"])

ax_w = sec_meas.plot_waterfall()

# exit()
# ax_w.get_figure().savefig("decay_waterfall.png")

ref_spec = sec_meas.reference_spectrum
resting_spec = sec_meas.get_spectrum(t=5)  # 5 seconds in, i.e. before the pulse
working_spec = sec_meas.get_spectrum(t=20)  # during the pulse.
decaying_spec = sec_meas.get_spectrum(t=40)  # after the pulse.

ax = resting_spec.plot(color="k", label="resting")
working_spec.plot(color="r", label="working", ax=ax)
decaying_spec.plot(color="b", label="decaying", ax=ax)
ref_spec.plot(color="0.5", linestyle="--", label="reference", ax=ax)
ax.legend()
# ax.get_figure().savefig("select raw spectra.png")

resting_OD_spec = sec_meas.get_dOD_spectrum(t=5)  # 5 seconds in, i.e. before the pulse
working_OD_spec = sec_meas.get_dOD_spectrum(t=20)  # during the pulse.
decaying_OD_spec = sec_meas.get_dOD_spectrum(t=40)  # after the pulse.

ax_OD = resting_OD_spec.plot(color="k", label="resting")
working_OD_spec.plot(color="r", label="working", ax=ax_OD)
decaying_OD_spec.plot(color="b", label="decaying", ax=ax_OD)
ax_OD.legend()
# ax_OD.get_figure().savefig("select OD spectra.png")
