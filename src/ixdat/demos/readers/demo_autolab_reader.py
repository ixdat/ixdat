"""Demo of the autolab reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_single_file(data_dir=None):
    """Return an ECMeasurement read from a single autolab test file."""
    path_to_file = get_demo_data_dir(data_dir) / "autolab/autolab_test_file.txt"
    # The file has no timestamp. Without `tstamp`, the reader prompts for one.
    return Measurement.read(
        path_to_file, reader="autolab", tstamp=path_to_file.stat().st_mtime
    )


EXAMPLES = {"single_file": load_single_file}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["single_file"].plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
