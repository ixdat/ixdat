"""Demo of the Avantage XPS reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Spectrum
from ixdat.demos import get_demo_data_dir


def load_xps_survey(data_dir=None):
    """Return an XPS survey Spectrum read from an Avantage .avg file."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(root / "avantage/XPS Survey.avg", reader="avantage")


EXAMPLES = {"xps_survey": load_xps_survey}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["xps_survey"].plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
