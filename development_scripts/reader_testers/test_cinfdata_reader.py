# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
from pathlib import Path
from matplotlib import pyplot as plt

from ixdat import Measurement

path_to_file = (
    Path.home()
    / "Dropbox/ixdat_resources/test_data/cinfdata/Trimarco2018_fig3/QMS_1.txt"
)

meas = Measurement.read(path_to_file, reader="cinfdata")

meas.plot_measurement()
