"""Demo of the ixdat .csv reader and exporter.

Reads a biologic file from the repository's ``test_data/``, exports it as an
ixdat .csv and reads it back. The tutorial loader reads a .csv written by ixdat
v0.1 from GitHub, so it requires an internet connection.
"""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_test_data_dir

TUTORIALS_URL = "https://raw.githubusercontent.com/ixdat/tutorials/"


def load_cv(data_dir=None):
    """Return a calibrated CyclicVoltammogram read from a biologic .mpt file."""
    root = get_test_data_dir(data_dir)
    meas = Measurement.read(root / "biologic/Pt_poly_cv_CUT.mpt", reader="biologic")
    meas.calibrate(RE_vs_RHE=0.01, R_Ohm=100, A_el=0.196)
    return meas.as_cv()


def load_tutorial_v0p1(data_dir=None):
    """Return an ECMeasurement read from an ixdat v0.1 .csv on the tutorials page.

    `data_dir` is not used. It is there to match the other loaders.
    """
    return Measurement.read_url(
        TUTORIALS_URL + "ixdat_v0p1/loading_appending_and_saving/co_strip.csv",
        reader="ixdat",
        aliases={
            "t": ["time/s"],
            "raw_current": ["raw current / [mA]"],
            "raw_potential": ["raw potential / [V]"],
        },
    )


EXAMPLES = {
    "cv": load_cv,
    "tutorial_v0p1": load_tutorial_v0p1,
}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    results["cv"].export("test.csv")
    results["reloaded"] = Measurement.read("test.csv", reader="ixdat")
    if show:
        results["reloaded"].plot()
        results["tutorial_v0p1"].plot_measurement()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
