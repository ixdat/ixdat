"""Demo of the qexafs reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Spectrum, Measurement
from ixdat.demos import get_demo_data_dir

FOLDER = "qexafs/constant potential"


def load_xas(data_dir=None):
    """Return an XAS Spectrum read from a single qexafs .dat file."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(
        root / FOLDER / "540117_IrO2_crys_0.60V_1.dat", reader="qexafs", technique="XAS"
    )


def load_multi_spec(data_dir=None):
    """Return a MultiSpectrum with all the columns of a single qexafs .dat file."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(root / FOLDER / "540117_IrO2_crys_0.60V_1.dat", reader="qexafs")


def load_xas_series(data_dir=None):
    """Return an XAS SpectrumSeries read from a set of qexafs .dat files."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read_set(
        part=root / FOLDER / "IrO2_crys", suffix=".dat", reader="qexafs", technique="XAS"
    )


def load_ec(data_dir=None):
    """Return an ECMeasurement read from a biologic .mpt file."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(root / FOLDER / "IrO2_CA_0p60_C02.mpt", reader="biologic")


def load_ec_xas(data_dir=None):
    """Return an EC-XAS measurement combining the EC data and the XAS series."""
    ec_xas = load_ec(data_dir=data_dir) + load_xas_series(data_dir=data_dir)
    ec_xas.tstamp += ec_xas.t[0]
    return ec_xas


EXAMPLES = {
    "xas": load_xas,
    "multi_spec": load_multi_spec,
    "xas_series": load_xas_series,
    "ec": load_ec,
    "ec_xas": load_ec_xas,
}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["xas"].plot()

        ax = results["multi_spec"]["QexafsFFI0"].plot()
        ax2 = ax.twinx()
        results["multi_spec"]["I0"].plot(ax=ax2, color="r")
        ax2.set_ylim([-1e6, 2e6])

        ax = results["xas_series"].heat_plot()
        ax.get_figure().savefig("xas_heat_plot.png")

        results["ec"].plot()
        results["ec_xas"].plot(xspan=[11200, 11300])
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
