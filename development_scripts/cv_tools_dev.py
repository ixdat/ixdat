# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
from pathlib import Path
from matplotlib import pyplot as plt

from ixdat.techniques import CyclicVoltammagram

plt.close("all")

stripping_cycle = CyclicVoltammagram.get(1)
base_cycle = CyclicVoltammagram.get(2)

diff = stripping_cycle.diff_with(base_cycle)

diff.plot_diff()
diff.plot_measurement()
diff.plot()
