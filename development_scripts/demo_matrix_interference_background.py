# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:26:19 2023

@author: Søren
"""
from pathlib import Path
from ixdat import Measurement


data_dir = Path.home() / (
    "Dropbox/WORKSPACES/ICL people/Daisy/data/21D_LIB/2021-04-28 14_56_39 21D28_NMC_Li_mesh"
)
export_file = data_dir / "exported.csv"
if False:  # load raw data
    ms = Measurement.read(
        data_dir / "2021-04-28 14_56_39 21D28_NMC_Li_mesh.tsv",
        technique="MS",
    )
    ms.export(export_file, tspan=ms.tspan, time_step=120)
else:  # load exported data
    ms = Measurement.read(export_file, reader="ixdat")


ms.plot()


ms.plot(mass_list=["M15", "M44", "M26", "M2"], logplot=False)

background = ms.remove_matrix_interference(
    mass_ref="M15", mass_list=["M44", "M26"], tspan_norm=[16500, 17000]
)

diff1 = background.remove_bg_from_series("M44-bg-subtracted").data - ms["M44"].data

print("diff1:")
# print(diff1)

# sanity check:
ratio = diff1 / ms.grab_for_t("M15", ms["M44"].t)
print("ratio:")
# print(ratio)  # should be constant


diff2 = (
    ms.grab("M44", remove_background=True)[1]
    - ms.grab("M44", remove_background=False)[1]
)
print("diff2 (should equal diff1):")
# print(diff2)

ms.plot(mass_list=["M15", "M44", "M26", "M2"], logplot=False, remove_background=True)
