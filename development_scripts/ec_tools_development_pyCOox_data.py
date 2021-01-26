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
    data_dir = Path(
        "~/Dropbox/ixdat_resources/20B12_Data_Analysis_Workshop/example data set/"
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
    ocp_meas = Measurement.get(3)
    cv_meas = Measurement.get(4)
    cp_meas = Measurement.get(5)

combined_meas = ocp_meas + cv_meas + cp_meas

v_list = ["Ewe/V", "<I>/mA", "I/mA"]
t_start = 13700
tspan = [0, 3000]

if False:
    for meas in [ocp_meas, cv_meas, cp_meas, combined_meas]:
        meas.tstamp += t_start

        meas.plot(v_list=v_list, tspan=tspan)
else:
    combined_meas.tstamp += t_start

combined_meas.plot(tspan=tspan)
