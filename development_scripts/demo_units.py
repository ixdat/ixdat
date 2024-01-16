# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:15:11 2023

@author: Søren
"""
from pathlib import Path
from ixdat import Measurement
from matplotlib import pyplot as plt
from ixdat.units import ureg


ureg.setup_matplotlib(True)


PATH_TO_DATAFILE = Path(__file__).parent / "../test_data/biologic/Pt_poly_cv_CUT.mpt"

ec = Measurement.read(PATH_TO_DATAFILE, reader="biologic")

t_hr, I_uA = ec.grab(
    "raw_current",
    unit_name="uA",
    t_unit_name="hour",
    return_quantity=False,
)

fig, ax = plt.subplots()
ax.yaxis.set_units(I_uA.u)
ax.xaxis.set_units(t_hr.u)

ax.plot(t_hr, I_uA, "k")

ax.plot(t_hr, I_uA.to("mA") * 0.7, "b")


print(ec.potential)
