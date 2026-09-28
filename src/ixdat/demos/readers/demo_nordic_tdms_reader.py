"""Demo of the nordic .tdms reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_single_file(data_dir=None):
    """Return an ECMeasurement read from a single nordic .tdms file."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / "nordic_tdms/24B07_0_Pt/CV_101448_ 1.tdms", reader="nordic"
    )


def load_cv(data_dir=None):
    """Return a CyclicVoltammogram combining all the .tdms files in a folder."""
    root = get_demo_data_dir(data_dir)
    meas = Measurement.read_set(
        root / "nordic_tdms/24B07_0_Pt", reader="nordic", suffix=".tdms"
    )
    cv = meas.as_cv()
    cv.redefine_cycle(start_potential=0.4, redox=True)
    return cv


EXAMPLES = {"single_file": load_single_file, "cv": load_cv}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["single_file"].plot()
        results["cv"].plot()
        results["cv"][15:30].plot_cycles()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
