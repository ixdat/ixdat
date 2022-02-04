# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
from pathlib import Path
from matplotlib import pyplot as plt

from ixdat.techniques import CyclicVoltammogram

plt.close("all")

stripping_cycle = CyclicVoltammogram.get(1)
base_cycle = CyclicVoltammogram.get(2)

diff = stripping_cycle.diff_with(base_cycle)

diff.plot_diff()
diff.plot_measurement()
diff.plot()
