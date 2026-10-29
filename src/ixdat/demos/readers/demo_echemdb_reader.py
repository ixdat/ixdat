"""Demo of the EChemDB reader, comparing a measured CV with an EChemDB reference.

Reads the measured CV from the repository's ``test_data/`` and the reference
from echemdb.org, so it requires an internet connection.

@author: Søren
@contributor: Frederik
"""

import matplotlib.pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_test_data_dir

REF_ID = "briega-martos_2021_cation_48_f1cs_black"


def load_measured_cv(data_dir=None):
    """Return a CyclicVoltammogram read from a biologic .mpt file."""
    root = get_test_data_dir(data_dir)
    return Measurement.read(root / "biologic/Pt_poly_cv.mpt", reader="biologic").as_cv()


def load_reference_cv(data_dir=None):
    """Return a reference CyclicVoltammogram from EChemDB, calibrated to RHE.

    `data_dir` is not used. It is there to match the other loaders.
    """
    ref_cv = Measurement.read(REF_ID, reader="echemdb").as_cv()
    # Nernst shift SHE to RHE
    ref_cv.calibrate(RE_vs_RHE=0.060 * ref_cv.metadata["electrolyte"]["ph"]["value"])
    return ref_cv


EXAMPLES = {"measured_cv": load_measured_cv, "reference_cv": load_reference_cv}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if not show:
        return results

    fig, ax = plt.subplots(figsize=(6, 4))
    results["measured_cv"][2].plot(ax=ax, color="C0", label="Measured cycle #3")
    results["reference_cv"].plot(ax=ax, color="C2", label="EchemDB ref vs RHE")

    ax.set_xlabel("Potential vs. RHE (V)")
    ax.set_ylabel("Current density (A/m²)")
    ax.legend(loc="best")
    ax.grid(True)
    plt.tight_layout()
    plt.show()
    return results


if __name__ == "__main__":
    results = main()
