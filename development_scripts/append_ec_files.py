# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
from pathlib import Path
from matplotlib import pyplot as plt

from ixdat import Measurement

plt.close("all")

data_dir = Path.home() / ("Dropbox/ixdat_resources/tutorials_data/extended_platinum_ec/")

if True:  # test reading from directory
    full_measurement = Measurement.read_set(data_dir, suffix="mpt")
    full_measurement.plot(J_name="cycle number")
    # test selection:
    cycle = full_measurement.select([1, 2])
    cycle.plot()
    full_measurement.cut([1000, 2000]).plot()
    full_measurement.select_values(loop_number=0, **{"cycle number": 1}).plot()
    # # The following line raises a TypeError:
    # full_measurement.select_values(loop_number=1, **{"loop_number": 0}).plot()

if True:

    ocp_file = data_dir / "01_demo_02_OCV_C01.mpt"
    cv_file = data_dir / "01_demo_03_CVA_C01.mpt"
    cp_file = data_dir / "01_demo_04_CP_C01.mpt"

    ocp_meas = Measurement.read(ocp_file, name="Pt_demo_ocp")
    print("read ocp file!")
    cv_meas = Measurement.read(cv_file, reader="biologic", name="Pt_demo_cv")
    print("read cv file!")
    cp_meas = Measurement.read(cp_file, reader="biologic", name="Pt_demo_cp")
    print("read cp file!")

    ocp_id = ocp_meas.save()
    cv_id = cv_meas.save()
    cp_id = cp_meas.save()
else:
    ocp_meas = Measurement.get(1)
    cv_meas = Measurement.get(2)
    cp_meas = Measurement.get(3)

combined_meas = ocp_meas + cv_meas + cp_meas

t_start = 13700
tspan = [0, 5000]

combined_meas.calibrate(RE_vs_RHE=0.715, A_el=0.196)

combined_meas.tstamp += t_start

combined_meas.plot(tspan=tspan)
combined_meas.plot(tspan=tspan, J_name="raw_CE_potential")

cut_meas = combined_meas.cut(tspan=tspan)
cut_meas.plot(J_name="selector")
cut_meas.rebuild_selector(columns=["cycle number"])
cut_meas.plot(J_name="selector")

select_meas = cut_meas.select_values(selector=range(4, 8))
select_meas.correct_ohmic_drop(R_Ohm=200)
select_meas.plot_vs_potential()

if False:  # check a save-load round trip
    if True:
        select_meas.name = "selected_measurement"
        select_meas.save()  # this changes its ID!
        my_id = select_meas.id
    else:
        my_id = 4

    # check the loaded measurement

    loaded_meas = Measurement.get(my_id)

    loaded_meas.plot_vs_potential()
    loaded_meas.plot(J_str="selector")

# check the cyclic voltammatry selection stuff

cv = combined_meas.as_cv()

cv.calibrate(RE_vs_RHE=0.72)

cv_selection = cv[10:16]

cv_selection.plot_measurement(J_name="cycle")
cv_selection.redefine_cycle(start_potential=0.45, redox=1)
cv_selection.plot_measurement(J_name="cycle")

ax = cv_selection[1].plot(label="cycle 1")
cv_selection[2].plot(ax=ax, linestyle="--", label="cycle 2")
ax.legend()

cv_selection.export(path_to_file=Path.home() / "test.csv")
