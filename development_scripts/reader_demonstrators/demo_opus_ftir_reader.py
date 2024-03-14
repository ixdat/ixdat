# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 16:20:40 2023

@author: Søren
"""
from pathlib import Path
from ixdat import Spectrum, Measurement
from matplotlib import pyplot as plt

plt.close("all")


data_dir = (
    Path.home() / "Dropbox/WORKSPACES/ICL people/Matthew/opus_ftir/DPT files for Soren"
)


ftir = Spectrum.read(
    data_dir / "231205 DME 3% EtOH",
    time_first="05/12/2023 15:20:33.696 (GMT+0)",  #  %d/%m/%Y %H:%M:%S.%f  23 chars
    time_last="05/12/2023 17:59:57.725 (GMT+0)",
    reader="opus_ftir",
)

ftir.heat_plot()  # heat plot
ftir.plot_waterfall()  # waterfall plot
ftir.plot(
    dn=20,
    xspan=[1000, 1500],
    xspan_bg=[1000, 1020],
    scale_factor=2,
    color="k",
    y_values="n",
    average=False,
)  # stacked spectra plot


ec = Measurement.read_set(
    data_dir / "231205 DME 3% EtOH",
    suffix=".mpt",
    reader="biologic",
)
ec.calibrate(R_Ohm=30)
ec.plot()

ecftir = ec + ftir

ecftir.plot_measurement()
ecftir.plot_stacked_spectra(
    dn=20,
    xspan=[800, 1800],
    xspan_bg=[800, 820],
)
