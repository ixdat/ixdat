"""Demo of the XRDML reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Spectrum
from ixdat.demos import get_demo_data_dir


def load_gi_xrd(data_dir=None):
    """Return a grazing-incidence XRD Spectrum read from an .xrdml file."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(
        root / "xrdml/GI-XRD Path 2_1 omega 0p5 step 10s.xrdml", reader="xrdml"
    )


EXAMPLES = {"gi_xrd": load_gi_xrd}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["gi_xrd"].plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
