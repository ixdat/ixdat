# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 22:02:57 2021

@author: scott
"""
from pathlib import Path
from matplotlib import pyplot as plt

from ixdat import Measurement

plt.close("all")

if False:
    data_dir = Path.home() / (
        "Dropbox/ixdat_resources/20B12_Data_Analysis_Workshop/example data set/"
    )

    ocp_file = data_dir / "01_Trimi1_cont_02_OCV_C01.mpt"
    cv_file = data_dir / "01_Trimi1_cont_03_CVA_C01.mpt"
    cp_file = data_dir / "01_Trimi1_cont_04_CP_C01.mpt"

    ocp_meas = Measurement.read(ocp_file, reader="biologic", name="Pt_demo_ocp")
    print("read ocp file!")
    cv_meas = Measurement.read(cv_file, reader="biologic", name="Pt_demo_cv")
    print("read cv file!")
    cp_meas = Measurement.read(cp_file, reader="biologic", name="Pt_demo_cp")
    print("read cp file!")

    ocp_id = ocp_meas.save()
    cv_id = cv_meas.save()
    cp_id = cp_meas.save()
else:
    ocp_meas = Measurement.get(7)
    cv_meas = Measurement.get(8)
    cp_meas = Measurement.get(9)

combined_meas = ocp_meas + cv_meas + cp_meas

t_start = 13700
tspan = [0, 3000]

combined_meas.calibrate(RE_vs_RHE=0.715, A_el=0.196)

combined_meas.tstamp += t_start

combined_meas.plot(tspan=tspan)

cut_meas = combined_meas.cut(tspan=tspan)
cut_meas.plot(J_str="selector")

select_meas = cut_meas.select_value(selector=7)
select_meas.correct_ohmic_drop(R_Ohm=200)
select_meas.plot()

select_meas.name = "selected_measurement"
select_meas.save()  # this changes its ID!
my_id = select_meas.id

del select_meas

loaded_meas = Measurement.get(my_id)

loaded_meas.plot()
