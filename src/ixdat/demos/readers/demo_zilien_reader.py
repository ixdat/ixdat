"""Demo of the zilien reader(s). Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir
from ixdat.techniques import MSMeasurement, ECMeasurement

FILE_NAME = "zilien_with_ec/2021-02-01 17_44_12.tsv"


def load_ecms(data_dir=None):
    """Return an ECMSMeasurement with both the EC and MS data of a zilien file."""
    root = get_demo_data_dir(data_dir)
    ecms = Measurement.read(root / FILE_NAME, reader="zilien")
    ecms.calibrate_RE(0)
    return ecms


def load_ms(data_dir=None):
    """Return an MSMeasurement with only the MS data of a zilien file."""
    root = get_demo_data_dir(data_dir)
    return MSMeasurement.read(root / FILE_NAME, reader="zilien")


def load_ec(data_dir=None):
    """Return an ECMeasurement with only the EC data of a zilien file."""
    root = get_demo_data_dir(data_dir)
    return ECMeasurement.read(root / FILE_NAME, reader="zilien")


def load_ecms_biologic(data_dir=None):
    """Return an ECMSMeasurement combining zilien MS data and biologic EC data."""
    root = get_demo_data_dir(data_dir)
    ec = Measurement.read_set(
        root / "zilien_with_ec/2021-02-01 17_44_12", reader="biologic", suffix=".mpt"
    )
    return ec + load_ms(data_dir=data_dir)


EXAMPLES = {
    "ecms": load_ecms,
    "ms": load_ms,
    "ec": load_ec,
    "ecms_biologic": load_ecms_biologic,
}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        ecms = results["ecms"]
        ecms.plot_measurement()
        results["ms"].plot_measurement()  # one panel, no EC
        results["ecms_biologic"].plot_measurement()
        results["ec"].plot()

        # Plot only the EC or only the MS data of the EC-MS measurement
        ecms.ec_plotter.plot_measurement()
        ecms.ec_plotter.plot_vs_potential()
        ecms.ms_plotter.plot_measurement()

        ecms_cv = ecms.as_cv()
        ecms_cv.ec_plotter.plot_vs_potential()
        ecms_cv.plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
