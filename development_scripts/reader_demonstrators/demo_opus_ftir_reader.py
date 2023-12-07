# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 16:20:40 2023

@author: Søren
"""
from pathlib import Path
from ixdat import Spectrum, Measurement


data_dir = Path.home() / "Dropbox/WORKSPACES/ICL people/Matthew/opus_ftir/DPT files for Soren"


ftir = Spectrum.read(
   data_dir / "231205 DME 3% EtOH",
   time_first="05/12/2023 15:20:33.696 (GMT+0)",  #  %d/%m/%Y %H:%M:%S.%f  23 chars
   time_last="05/12/2023 17:59:57.725 (GMT+0)",
   reader="opus_ftir"
)

# ftir.plot()  # waterfall plot


ec = Measurement.read_set(
   data_dir / "231205 DME 3% EtOH",
   suffix=".mpt",
   reader="biologic",
)
ec.calibrate(R_Ohm=30)

ecftir = ec + ftir
ecftir.plot(
   V_step=0.1,  # plot a spectrum for every 0.1 V
)  # stacked with potential on y-axis
