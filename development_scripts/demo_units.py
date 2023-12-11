# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:15:11 2023

@author: Søren
"""
from pathlib import Path
from ixdat import Measurement
from matplotlib import pyplot as plt


PATH_TO_DATAFILE = Path(__file__).parent / "../test_data/biologic/Pt_poly_cv_CUT.mpt"

ec_measurement = Measurement.read(PATH_TO_DATAFILE, reader="biologic")

t_hr, I_uA = ec_measurement.grab("raw_current", unit_name="uA", t_unit_name="hour")

fig, ax = plt.subplots()

ax.plot(t_hr, I_uA, "k")
