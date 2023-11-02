# -*- coding: utf-8 -*-
"""
This module is a worked example which answers the challenge at the end of the
ipython notebook "difference_between_two_cvs" here:
https://github.com/ixdat/tutorials/blob/main/simple_ec_analysis/difference_between_two_cvs.ipynb

I suggest to go through that notebook and give the exercise a try yourself
before taking a look at the suggested solutions here.
"""
import numpy as np
from matplotlib import pyplot as plt
from ixdat.techniques import CyclicVoltammogram as CV

plt.close("all")

meas = CV.read_url(
    "https://raw.githubusercontent.com/ixdat/tutorials/main/electrochemistry/data/oxide_reduction.csv",
    reader="ixdat",
)

meas.plot_measurement()
meas.tstamp = meas.tstamp + meas.t[0]


if True:  # Method 1
    tspan_red = [362, 377]
    t_red, I_red = meas.grab("raw_current", tspan=tspan_red)
    v_red = meas.grab_for_t("potential", t_red)

    tspan_base = [470, 485]
    t_base, I_base = meas.grab("raw_current", tspan=tspan_base)
    v_base = meas.grab_for_t("potential", t_base)

    fig, ax = plt.subplots()
    ax.plot(v_red, I_red, "g")
    ax.plot(v_base, I_base, "k")

    Q_reduction_cycle = np.trapz(I_red, t_red) * 1e-3
    Q_base_cycle = np.trapz(I_base, t_base) * 1e-3

    Q_red = Q_reduction_cycle - Q_base_cycle
    print(f"----------\nMethod 1: reduced with {Q_red*1e6} uC passed.\n---------")


#  this is used for both methods 2 and 3
meas = meas.cut([340, 700])
meas.redefine_cycle(start_potential=1.05, redox=False)
meas.plot_measurement(J_str="cycle")

if True:  # Method 2
    Q_reduction_cycle = (
        meas[1].integrate(
            "raw_current",
            vspan=[1.0, 0.5],
            ax="new",
        )
        * 1e-3
    )
    Q_base_cycle = (
        meas[2].integrate(
            "raw_current",
            vspan=[1.0, 0.5],
            ax="new",
        )
        * 1e-3
    )
    Q_red = Q_reduction_cycle - Q_base_cycle

    print(f"----------\nMethod 2: reduced with {Q_red*1e6} uC passed.\n---------")

if True:  # Method 3
    diff = meas[1].diff_with(meas[2])
    diff.plot()
    Q_red = diff.integrate("raw_current", tspan=[360, 380], ax="new") * 1e-3

    print(f"----------\nMethod 3:  reduced with {Q_red*1e6} uC passed.\n---------")
