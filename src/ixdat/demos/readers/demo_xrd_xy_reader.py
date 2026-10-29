"""Demo of the XRDXYReader on two powder diffraction datasets read from URLs.

Requires an internet connection.

  1. ZSM-5 zeolite -- synchrotron .xye (3-column: 2-theta, intensity, error), no header
     Source: https://github.com/stefsmeets/lines (MIT license)

  2. Rutile TiO2 — Debye-calculated I(Q) .dat (comma-separated, header "Q,I(Q)")
     ICSD entry 001504, tetragonal, from DebyeCalculator
     Source: https://github.com/FrederikLizakJohansen/DebyeCalculator (MIT license)
"""

import matplotlib.pyplot as plt

from ixdat import Spectrum

ZSM5_URL = "https://raw.githubusercontent.com/stefsmeets/lines/main/testing/zsm-5.xye"
RUTILE_URL = (
    "https://raw.githubusercontent.com/FrederikLizakJohansen/DebyeCalculator"
    "/main/debyecalculator/unittests_files"
    "/icsd_001504_cc_r6_lc_2.85_6_tetragonal_Iq.dat"
)


def load_zsm5(data_dir=None):
    """Return an XRD Spectrum of ZSM-5 zeolite vs 2-theta, with an error column.

    The file has no header, so the reader defaults to "two theta / degree".
    `data_dir` is not used. It is there to match the other loaders.
    """
    return Spectrum.read_url(
        ZSM5_URL, reader="xrdxy", name="ZSM-5 zeolite (synchrotron)"
    )


def load_rutile(data_dir=None):
    """Return a calculated XRD Spectrum of rutile TiO2 vs Q.

    The header line "Q,I(Q)" tells the reader that the x-axis is Q.
    `data_dir` is not used. It is there to match the other loaders.
    """
    return Spectrum.read_url(RUTILE_URL, reader="xrdxy", name="Rutile TiO2 I(Q)")


EXAMPLES = {"zsm5": load_zsm5, "rutile": load_rutile}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    for spec in results.values():
        print(spec)
    if show:
        for spec in results.values():
            spec.plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
