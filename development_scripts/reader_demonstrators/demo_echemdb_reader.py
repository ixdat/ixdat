# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 10:48:47 2025

@author: Søren
"""

from pathlib import Path
from ixdat import Measurement



my_cv = Measurement.read(
    Path.home() / "git/ixdat/test_data/biologic/Pt_poly_cv_CUT.mpt",
    reader="biologic"
).as_cv()

my_cv_cycle = my_cv[2]  # a CyclicVoltammogram object
ax = my_cv_cycle.plot(color="r")



ref_cv = Measurement.read(
    "briega_martos_2018_understanding_j3045_f1c_black", reader="echemdb"
)  # also a CyclicVoltammogram object
ref_cv.plot(ax=ax, color="k", linestyle="--")





