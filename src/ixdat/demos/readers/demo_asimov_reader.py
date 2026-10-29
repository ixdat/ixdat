"""Demo of the Asimov reader with an NMR spectrum and a Biologic CV.

The first run opens the normal Asimov login flow.
"""

import matplotlib.pyplot as plt

from ixdat import Measurement, Spectrum

NMR_SPECTRUM_ID = "300f4dcf-1585-51a6-81bc-f867910efe1a"
BIOLOGIC_CV_ID = "fbbf4edb-f288-5e85-bed0-e3bc89edc2e5"


def load_nmr_spectrum(data_dir=None):
    """Return the permanent NMR spectrum example from Asimov.

    `data_dir` is not used. It is there to match the other loaders.
    """
    return Spectrum.read(NMR_SPECTRUM_ID, reader="asimov")


def load_biologic_cv(data_dir=None):
    """Return the permanent Biologic cyclic voltammetry example from Asimov.

    `data_dir` is not used. It is there to match the other loaders.
    """
    return Measurement.read(BIOLOGIC_CV_ID, reader="asimov")


EXAMPLES = {"nmr_spectrum": load_nmr_spectrum, "biologic_cv": load_biologic_cv}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        results["nmr_spectrum"].plot()
        results["biologic_cv"].plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
