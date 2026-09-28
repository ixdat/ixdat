"""Demo of the MSRH SEC reader. Requires the demo data.

MSRH = molecular science research hub, at Imperial College London.
"""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_sec(data_dir=None):
    """Return a calibrated SpectroECMeasurement read from MSRH SEC files."""
    root = get_demo_data_dir(data_dir)
    sec_meas = Measurement.read(
        root / "sec/test-7SEC.csv",
        path_to_ref_spec_file=root / "sec/WL.csv",
        path_to_U_J_file=root / "sec/test-7_JV.csv",
        scan_rate=1,
        tstamp="now",
        reader="msrh_sec",
    )
    sec_meas.calibrate_RE(RE_vs_RHE=0.26)  # provide RE potential in [V] vs RHE
    sec_meas.normalize_current(A_el=1)  # provide electrode area in [cm^2]
    sec_meas.set_reference_spectrum(V_ref=0.66)
    return sec_meas


EXAMPLES = {"sec": load_sec}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    sec_meas = results["sec"]

    # Test export and reload
    export_name = "exported_sec.csv"
    sec_meas.export(export_name)
    sec_reloaded = Measurement.read(export_name, reader="ixdat")
    sec_reloaded.set_reference_spectrum(V_ref=0.66)
    results["reloaded"] = sec_reloaded

    if show:
        ax = sec_meas.get_dOD_spectrum(V=1.0, V_ref=0.66).plot(
            color="b", label="species 1"
        )
        sec_meas.get_dOD_spectrum(V=1.4, V_ref=1.0).plot(
            color="g", label="species 2", ax=ax
        )
        sec_meas.get_dOD_spectrum(V=1.7, V_ref=1.4).plot(
            color="r", label="species 3", ax=ax
        )
        ax.legend()

        sec_reloaded.plot_vs_potential(cmap_name="jet")
        sec_reloaded.continuous = False
        sec_reloaded.plot_vs_potential(cmap_name="jet")

        sec_meas.plot_measurement(V_ref=0.4, cmap_name="jet", make_colorbar=True)
        ax = sec_meas.plot_waterfall(V_ref=0.4, cmap_name="jet", make_colorbar=True)
        ax.get_figure().savefig("sec_waterfall.png")

        sec_meas.plot_vs_potential(V_ref=0.66, cmap_name="jet", make_colorbar=False)
        sec_meas.plot_vs_potential(
            V_ref=0.66,
            vspan=[1.0, 1.5],
            wlspan=[500, 700],
            cmap_name="jet",
            make_colorbar=False,
        )

        ax = sec_meas.get_dOD_spectrum(V_ref=0.66, V=1.0).plot(
            color="b", label="species 1"
        )
        sec_meas.get_dOD_spectrum(V_ref=1.0, V=1.45).plot(
            color="g", ax=ax, label="species 2"
        )
        sec_meas.get_dOD_spectrum(V_ref=1.45, V=1.75).plot(
            color="r", ax=ax, label="species 3"
        )
        ax.legend()

        axes = sec_meas.plot_wavelengths_vs_potential(
            wavelengths=["w460", "w600", "w850"]
        )
        axes[0].set_ylabel("intense!")
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
