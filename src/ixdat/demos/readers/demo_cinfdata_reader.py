"""Demo of the cinfdata reader, reproducing Trimarco et al. 2018, Fig. 3.

Requires the demo data.
"""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_ms(data_dir=None):
    """Return an MSMeasurement read from a cinfdata .txt file."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / "cinfdata/Trimarco2018_fig3/QMS_1.txt", reader="cinfdata"
    )


def load_ec(data_dir=None):
    """Return a calibrated ECMeasurement read from a set of biologic files."""
    root = get_demo_data_dir(data_dir)
    ec_meas = Measurement.read_set(
        root / "cinfdata/Trimarco2018_fig3/09_fig4", reader="biologic"
    )
    ec_meas.calibrate(RE_vs_RHE=0.65, A_el=0.196)
    return ec_meas


def load_ecms(data_dir=None):
    """Return an ECMSMeasurement combining the cinfdata MS and biologic EC data."""
    return load_ec(data_dir=data_dir) + load_ms(data_dir=data_dir)


EXAMPLES = {"ms": load_ms, "ec": load_ec, "ecms": load_ecms}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    ecms_meas = results["ecms"]
    ecms_meas.export("trimarco2018_fig3_data.csv", tspan=[0, 600])

    if show:
        results["ms"].plot_measurement()
        results["ec"].plot_measurement()

        axes = ecms_meas.plot_measurement(
            mass_lists=[["M4", "M28"], ["M44", "M2"]],
            tspan_bg=[None, [30, 40]],
            legend=False,
            unit="pA",
        )
        axes[2].set_ylim([-7, 70])
        axes[0].set_ylim([-1.8e3, 18e3])
        axes[0].get_figure().tight_layout()

    ecms_meas.set_bg(tspan_bg=[0, 10])
    cv = ecms_meas.as_cv()
    cv.redefine_cycle(start_potential=0.39, redox=False)
    results["cv"] = cv

    if show:
        axes_cv = cv[2].plot(mass_list=["M2", "M44"], logplot=False)
        axes_cv = cv[1].plot(
            mass_list=["M2", "M44"], linestyle="--", axes=axes_cv, logplot=False
        )
        axes_cv[0].get_figure().savefig("Trimarco2018_ixdat.png")
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
