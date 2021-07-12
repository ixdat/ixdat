"""Demonstrate simple importing and plotting SEC decay data"""

from pathlib import Path
from ixdat import Measurement

data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/sec"

sec_meas = Measurement.read(
    data_dir / "decay/PDtest-1.35-1OSP-SP.csv",
    path_to_ref_spec_file=data_dir / "WL.csv",
    path_to_t_V_file=data_dir / "decay/PDtest-1.35-1OSP-E-t.csv",
    path_to_t_J_file=data_dir / "decay/PDtest-1.35-1OSP-J-t.csv",
    tstamp=1,
    reader="msrh_sec_decay",
)

if True:  # Replace reference with spectrum at t=5
    from ixdat.data_series import Field

    ref_spec = sec_meas.get_spectrum(t=5)
    reference = Field(
        name="reference",
        unit_name="counts",
        data=ref_spec.field.data,
        axes_series=ref_spec.field.axes_series,
    )
    sec_meas["reference"] = reference

axes = sec_meas.plot_measurement(
    # V_ref=0.66,  # can't do a V_ref for this as can't interpolate on potential..
    # So OD will be calculated using the reference spectrum in WL.csv
    cmap_name="jet",
    make_colorbar=False,
)
# axes[0].get_figure().savefig("decay_vs_t.png")

ax_w = sec_meas.plot_waterfall()

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
