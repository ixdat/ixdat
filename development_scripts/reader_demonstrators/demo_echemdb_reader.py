# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 10:48:47 2025

@author: Søren
@contributor: Frederik
"""

from pathlib import Path
from ixdat import Measurement
import matplotlib.pyplot as plt
import numpy as np

# load CV from a Biologic .mpt file and convert to EC technique
mpt_path = Path(__file__).parent / "../../test_data/biologic/Pt_poly_cv_CUT.mpt"
my_cv = Measurement.read(mpt_path, reader="biologic").as_cv()

# pull out the 3rd cycle and plot
my_cycle = my_cv[2]
fig, ax = plt.subplots(figsize=(6, 4))
my_cycle.plot(ax=ax, color="tab:red", label="Measured cycle #3")

# load reference from EChemDB
ref_id = "briega-martos_2021_cation_48_f1Cs_black"
try:
    ref_cv = Measurement.read(ref_id, reader="echemdb")
    print(ref_cv)
    ref_cycle = ref_cv.as_cv()[1]
    ref_cycle.plot(ax=ax, color="C1", linestyle=":", label="EchemDB ref (latest)")

    # calibrate and re-plot
    ref_cycle.calibrate_RE(RE_vs_RHE=0.7)  # shift potential by +0.7 V
    ref_cycle.plot(ax=ax, color="k", linewidth=1, alpha=0.6, label="EchemDB ref vs RHE")

except Exception as e:
    print(f"[Warning] Could not load or plot latest EchemDB ref: {e}")
    raise e

# fetch an older reference (v0.4.1) and overlay, should be the same
try:
    ref_old = Measurement.read(ref_id, reader="echemdb", version="0.4.1").as_cv()
    ref_cycle_old = ref_old[1]
    ref_cycle_old.plot(ax=ax, color="C2", linestyle="-.", label="EchemDB ref (v0.4.1)")
    
    # calibrate and re-plot
    ref_cycle_old.calibrate_RE(RE_vs_RHE=0.7)  # shift potential by +0.7 V
    ref_cycle_old.plot(ax=ax, color="k", linewidth=1, alpha=0.6, label="EchemDB ref (v0.4.1) vs RHE")
except Exception as e:
    print(f"[Warning] Could not load v0.4.1 ref: {e}")

# Finalize plot
ax.set_xlabel("Potential vs. SHE (V)")
ax.set_ylabel("Current density (A/m²)")
ax.set_title("Comparison of Measured vs. Reference CVs")
ax.legend(loc="best")
ax.grid(True)

plt.tight_layout()
plt.show()
