"""Read and plot two examples from asimov

The first run opens the normal asimov login flow. The examples are an NMR
spectrum and a Biologic cyclic voltammetry measurement.
"""

import matplotlib.pyplot as plt

from ixdat import Measurement, Spectrum

NMR_SPECTRUM_ID = "300f4dcf-1585-51a6-81bc-f867910efe1a"
BIOLOGIC_CV_ID = "fbbf4edb-f288-5e85-bed0-e3bc89edc2e5"


def load_demo_objects():
    """Return the permanent NMR spectrum and Biologic CV Asimov examples."""
    nmr_spectrum = Spectrum.read(NMR_SPECTRUM_ID, reader="asimov")
    biologic_cv = Measurement.read(BIOLOGIC_CV_ID, reader="asimov")
    return nmr_spectrum, biologic_cv


def main():
    nmr_spectrum, biologic_cv = load_demo_objects()
    nmr_spectrum.plot()
    biologic_cv.plot()

    plt.show()


if __name__ == "__main__":
    main()
