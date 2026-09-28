"""Demo of the Pfeiffer reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_air_mid(data_dir=None):
    """Return an MSMeasurement of air read from a Pfeiffer MID .dat file."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root
        / "pfeiffer"
        / "MID_air, Position 1, RGA PrismaPro 200 44526001, 003-02-2021 17'41'12 - Bin.dat",  # noqa: E501
        reader="pfeiffer",
    )


EXAMPLES = {"air_mid": load_air_mid}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    meas = results["air_mid"]
    if show:
        meas.set_bg(tspan_bg=[180, 200])
        meas.plot_measurement(logplot=False)
        meas.reset_bg()
        meas.plot_measurement(logplot=False)
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
