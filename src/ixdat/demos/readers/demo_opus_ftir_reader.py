"""Demo of the OPUS FTIR reader combined with biologic EC data.

Requires the demo data.
"""

from matplotlib import pyplot as plt

from ixdat import Spectrum, Measurement
from ixdat.demos import get_demo_data_dir

FOLDER = "opus_ftir/dpt_from_Matthew/231205 DME 3% EtOH"


def load_ftir(data_dir=None):
    """Return an FTIR SpectrumSeries read from a set of OPUS .dpt files."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(
        root / FOLDER,
        time_first="05/12/2023 15:20:33.696 (GMT+0)",  # %d/%m/%Y %H:%M:%S.%f
        time_last="05/12/2023 17:59:57.725 (GMT+0)",
        reader="opus_ftir",
    )


def load_ec(data_dir=None):
    """Return an ohmic-drop corrected ECMeasurement read from biologic files."""
    root = get_demo_data_dir(data_dir)
    ec = Measurement.read_set(root / FOLDER, suffix=".mpt", reader="biologic")
    ec.calibrate(R_Ohm=30)
    return ec


def load_ec_ftir(data_dir=None):
    """Return an EC-FTIR measurement combining the EC data and the FTIR spectra."""
    return load_ec(data_dir=data_dir) + load_ftir(data_dir=data_dir)


EXAMPLES = {"ftir": load_ftir, "ec": load_ec, "ec_ftir": load_ec_ftir}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if not show:
        return results
    ftir = results["ftir"]

    ftir.heat_plot()
    ftir.plot_waterfall()
    ftir.plot(  # stacked spectra plot
        dt=1000,
        xspan=[1000, 1500],
        xspan_bg=[1000, 1020],
        color="k",
        y_values="n",
        average=False,
    )

    results["ec"].plot()
    results["ec_ftir"].plot_measurement()
    results["ec_ftir"].plot_stacked_spectra(
        dn=40,
        xspan=[1050, 1400],
        xspan_bg=[1050, 1060],
        scale_factor=1.5,
        # average=False,  # no averaging
        average=5,  # average 5 spectra each side
        # average=True,  # most averaging
    )
    plt.show()
    return results


if __name__ == "__main__":
    results = main()
