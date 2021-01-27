# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
# from pathlib import Path
from matplotlib import pyplot as plt

from ixdat import Measurement

plt.close("all")

my_id = 5
loaded_meas = Measurement.get(my_id)

cv = loaded_meas.as_cv()
cv.plot_measurement(J_str="cycle")
cv_selection = cv[10:16]

cv_selection.plot_measurement(J_str="cycle")
cv_selection.redefine_cycle(start_potential=0.4, redox=1)
cv_selection.plot_measurement(J_str="cycle")

ax = cv_selection[1].plot(label="cycle 1")
cv_selection[2].plot(ax=ax, linestyle="--", label="cycle 2")
ax.legend()
