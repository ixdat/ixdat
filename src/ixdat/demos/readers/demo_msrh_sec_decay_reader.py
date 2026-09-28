"""Demo of the MSRH SEC decay reader. Requires the demo data.

MSRH = molecular science research hub, at Imperial College London.
"""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_sec_decay(data_dir=None):
    """Return a calibrated SpectroECMeasurement of a potential pulse and decay."""
    root = get_demo_data_dir(data_dir)
    sec_meas = Measurement.read(
        # root / "sec/decay/PDtest-1.35-1OSP-SP.csv",
        root / "sec/decay/PDtest-1.33-1OSP-SP.csv",
        path_to_ref_spec_file=root / "sec/WL.csv",
        # path_to_t_V_file=root / "sec/decay/PDtest-1.35-1OSP-E-t.csv",
        # path_to_t_J_file=root / "sec/decay/PDtest-1.35-1OSP-J-t.csv",
        path_to_t_U_file=root / "sec/decay/PDtest-1.33-1OSP-E-t.csv",
        path_to_t_J_file=root / "sec/decay/PDtest-1.33-1OSP-J-t.csv",
        tstamp="now",
        reader="msrh_sec_decay",
    )
    sec_meas.calibrate_RE(RE_vs_RHE=0.26)
    sec_meas.set_reference_spectrum(t_ref=5)
    return sec_meas


EXAMPLES = {"sec_decay": load_sec_decay}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if not show:
        return results
    sec_meas = results["sec_decay"]

    sec_meas.plot_measurement(
        # V_ref=0.66,  # can't do a V_ref for this as can't interpolate on potential..
        # So OD will be calculated using the reference spectrum in WL.csv
        cmap_name="inferno",
        make_colorbar=False,
    )

    sec_meas.plot_wavelengths(wavelengths=["w500", "w600", "w700", "w800"])
    sec_meas.plot_waterfall()

    ax = sec_meas.get_spectrum(t=5).plot(color="k", label="resting")  # before pulse
    sec_meas.get_spectrum(t=20).plot(color="r", label="working", ax=ax)  # in pulse
    sec_meas.get_spectrum(t=40).plot(color="b", label="decaying", ax=ax)  # after
    sec_meas.reference_spectrum.plot(
        color="0.5", linestyle="--", label="reference", ax=ax
    )
    ax.legend()

    ax_OD = sec_meas.get_dOD_spectrum(t=5).plot(color="k", label="resting")
    sec_meas.get_dOD_spectrum(t=20).plot(color="r", label="working", ax=ax_OD)
    sec_meas.get_dOD_spectrum(t=40).plot(color="b", label="decaying", ax=ax_OD)
    ax_OD.legend()
    plt.show()
    return results


if __name__ == "__main__":
    results = main()
