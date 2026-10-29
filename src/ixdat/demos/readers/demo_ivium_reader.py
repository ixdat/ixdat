"""Demo of the ivium reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir
from ixdat.techniques import CyclicVoltammogram


def load_ec(data_dir=None):
    """Return an ECMeasurement read from a set of ivium files."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(root / "ivium/ivium_test_dataset", reader="ivium")


def load_cv(data_dir=None):
    """Return a CyclicVoltammogram read from a set of ivium files."""
    root = get_demo_data_dir(data_dir)
    cv = CyclicVoltammogram.read(root / "ivium/ivium_test_dataset", reader="ivium")
    cv.redefine_cycle(start_potential=0.4, redox=False)
    return cv


EXAMPLES = {"ec": load_ec, "cv": load_cv}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["ec"].plot_measurement()
        results["cv"].plot_measurement()
        for i in range(4):
            results["cv"][i].plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
