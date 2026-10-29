"""Demo of the ixdat .csv reader on a reactor pressure and temperature measurement.

Requires the demo data.
"""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def load_tpms(data_dir=None):
    """Return a reactor Measurement of pressure and temperature from an ixdat .csv."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / "cinfdata/Krabbe/baratron_temp_measurement.txt.csv",
        reader="ixdat",
        technique="reactor",
        aliases={"pressure": ["Reactor pressure"], "temperature": ["TC temperature"]},
    )


EXAMPLES = {"tpms": load_tpms}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        axes = results["tpms"].plot()
        axes[0].get_figure().savefig("tpms_plot.png")
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
