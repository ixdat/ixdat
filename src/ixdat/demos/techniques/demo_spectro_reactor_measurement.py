"""Demo of combining a reactor measurement with a series of mass spectra.

Requires the demo data.
"""

from matplotlib import pyplot as plt

from ixdat import Measurement, Spectrum
from ixdat.demos import get_demo_data_dir

FOLDER = "cinfdata/Krabbe"


def load_spectra(data_dir=None):
    """Return a SpectrumSeries of mass spectra read from an ixdat .txt export."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(root / FOLDER / "spectrumseries_mass_spec.txt", reader="ixdat")


def load_reactor(data_dir=None):
    """Return a reactor Measurement of pressure, temperature and mass signals."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / FOLDER / "masstime.csv",
        reader="ixdat",
        technique="reactor",
        aliases={"pressure": ["Reactor pressure"], "temperature": ["TC temperature"]},
    )


def load_spectro_reactor(data_dir=None):
    """Return the reactor Measurement combined with the mass spectra."""
    return load_reactor(data_dir=data_dir) + load_spectra(data_dir=data_dir)


EXAMPLES = {
    "spectra": load_spectra,
    "reactor": load_reactor,
    "spectro_reactor": load_spectro_reactor,
}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["spectra"].heat_plot()
        axes = results["spectro_reactor"].plot()
        axes[0].get_figure().savefig("spectro_tpms_plot.png")
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
