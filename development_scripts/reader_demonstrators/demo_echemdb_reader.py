# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 10:48:47 2025

@author: Søren
@contributor: Frederik
"""

from pathlib import Path
from ixdat import Measurement
import matplotlib.pyplot as plt


# Load and parse biologic .mpt file
mpt_path = Path(__file__).parent / "../../test_data/biologic/Pt_poly_cv_CUT.mpt",
my_cv = Measurement.read(mpt_path, reader="biologic").as_cv()

# Select a single CV cycle
my_cv_cycle = my_cv[2]  # a CyclicVoltammogram object
ax = my_cv_cycle.plot(color="r")

# Try to load a reference CV from echemdb
try:
    ref_path = "briega_martos_2018_understanding_j3045_f1c_black", 
    ref_cv = Measurement.read(ref_path, reader="echemdb")
    ref_cv.plot(ax=ax, color="k", linestyle="--")
except Exception as e:
    print(f"[Warning] Failed to load reference CV from echemdb: {e}")
finally:
    # QT rendering workaround [Frederik]
    plt.show()
