# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
# from pathlib import Path
from matplotlib import pyplot as plt

from ixdat import Measurement

plt.close("all")

my_id = 4
loaded_meas = Measurement.get(my_id)

loaded_meas.plot()
loaded_meas.plot_vs_potential()

cv = loaded_meas.as_cv()

cv.plot()
