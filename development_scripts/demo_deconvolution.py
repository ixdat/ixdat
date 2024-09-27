# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:43:48 2024

Script demonstrating the mass transport deconvolution class & methods

email a.winiwarter@imperial.ac.uk for the data

@author: awiniwar
"""

import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path
from ixdat.calculators.ecms_calculators import MSCalResult, ECMSImpulseResponse
from ixdat import Measurement
from ixdat import plugins

plugins.activate_siq()

# ------------------ Impulse response fitting of measured data + comparison with model
# insert your folder here:
data_folder = Path(
    r"C:\Users\Søren\Dropbox\ixdat_resources\test_data\deconvolution"  # noqa
)
# data_folder = Path.home() / "Dropbox/ixdat_resources/test_data/deconvolution"
data_reg_cell_pulses_o2 = Measurement.read(
    data_folder / "O2_pulses.csv",
    reader="ixdat",
)

data_reg_cell_pulses_o2.plot_measurement()


# ------------ O2 calibration
# First, we mannually make a siq calibration point:
F_o2 = MSCalResult(
    name="O2", mol="O2", mass="M32", cal_type="ECMS calibration", F=0.11755
).to_siq()
# Then put it in a siq Calculator:
calibration = plugins.siq.Calculator(cal_list=[F_o2])
# And get the calculator ready to use on measurements:
calibration.set_quantifier(mol_list=["O2"], mass_list=["M32"], carrier="He")

# Add the siq Calculator it to the measurement:
data_reg_cell_pulses_o2.add_calculator(calibration)
# This makes the flux of O2 available as `"n_dot_O2"`.

# prepare a figure to co-plot the measured and modelled data:
fig2, ax2 = plt.subplots(nrows=1, ncols=1)
ax2.set_xlabel("time / [s]")
ax2.set_ylabel("norm. intensity")

# pulse length: duration of the pulse, both for measured and modeled pulse. this is
# considered when normalzing the area of the pulse to 1
pulse_length = 100
# pulse inf: defines the end of the pulse. Needs to be equal or larger than pulse_length.
# This is considered for the baseline correction and the automatic label generation for
# the plot. Can be beneficial to be longer than pulse_length if there are signal spikes
# which impact normalizsation.
pulse_inf = 100
dt = 0.1  # (default is 0.1, 0.6 is as measured)

for t_impulse in [18172, 18472, 18772, 19072]:  # list of pulse start times. accuracy
    # to the second seems to be necessary
    # select the data for the above pulse start times
    pulse = data_reg_cell_pulses_o2.cut(tspan=[t_impulse - 10, t_impulse + pulse_inf])
    # shift timestamp so all pulses start at 0
    pulse.tstamp += 0.9355739 + t_impulse  # the number is from the beginning of the
    # full dataset which has been cut from the sample data to decrease size of file
    # determine the integrated flux during the pulse for automatic labelling
    pulse_int = (
        pulse.integrate_flux(mol="O2", tspan=[0, pulse_inf], tspan_bg=[-5, 0]) * 1e9
    )
    # extract the impulse response from data
    imp_resp = ECMSImpulseResponse.from_measurement(
        mol="O2",
        measurement=pulse,
        tspan=[-5, pulse_length],
        tspan_bg=[[-5, 0], [pulse_inf - 5, pulse_inf]],
    )
    # add the extracted impulse response to the figure
    ax2.plot(
        imp_resp.t_kernel,
        imp_resp.kernel,
        marker="o",
        ms=3,
        markerfacecolor="w",
        linestyle="",
        label="n$_{O2}$" + " = {:.2f} nmol".format(pulse_int),
    )

# model impulse response, manually vary the working distance (wd) to find the best fit
# for the data. in the example case, with HOR in the same cell, a working distance
# of 230um was measured (data not included)
wd = 220e-6
imp_resp_model = ECMSImpulseResponse.from_parameters(
    mol="O2",
    working_distance=wd,
    A_el=0.197,
    D=2.1e-9,  # optional if using siq
    H_v_cc=33,  # optional if using siq
    n_dot=None,  # optional if using siq
    carrier_gas=None,  # defaults to He
    gas_volume=1e-10,  # optional if using siq
    duration=pulse_length,  # needs to
    # be the same as experimental data
    # if normalized to peak area
    # (default!)
    dt=dt,
)
# x_model=np.arange(0, pulse_length, dt)
ax2.plot(
    imp_resp_model.t_kernel,
    imp_resp_model.kernel,
    linestyle="--",
    label="model= {:.0f} $\mu$m".format(wd * 1e6),  # noqa
)
ax2.legend()

# ----------- Deconvolution of ECMS data
# import data
ca_dark_day1 = Measurement.read(
    data_folder / "to_deconvolute_data.csv",
    reader="ixdat",
)
ca_dark_day1.calibrate(A_el=0.35**2 * np.pi)
ca_dark_day1.plot(mass_list=["M32"], logplot=False)

# Define model parameters and deconvolute the oxygen signal
wd = 180e-6  # different working distance as different day and cell
imp_resp_model = ECMSImpulseResponse.from_parameters(
    mol="O2",
    working_distance=wd,
    A_el=0.35**2 * np.pi,
    D=2.1e-9,  # optional if using siq
    H_v_cc=33,  # optional if using siq
    carrier_gas=None,  # default is He
    gas_volume=1e-10,  # optional if using siq
)


# now let's grab the deconvoluted current for a tspan of interest
# first, need to calibrate the measurement, we'll use siq here
ca_dark_day1.add_calculator(calibration)
# now let's deconvolute
t_deconvoluted, v_deconvoluted = imp_resp_model.calc_deconvoluted_signal(
    mol="O2",
    measurement=ca_dark_day1,
    tspan=[-10, 210],
    tspan_bg=[-10, -1],
)
# and check what we got
fig3, ax3 = plt.subplots(nrows=1, ncols=1)
ax3.plot(t_deconvoluted, v_deconvoluted)
print(type(t_deconvoluted))
print(type(v_deconvoluted))


# Often we have more tspans that are interesting, so let's loop through all the data
# we want to deconvolute
tspan_list_1_dark = [[-10, 210], [208, 580]]
# each tspan needs to include 0 if no t_zero given
t_zero_list_1_dark = [6.6, 227]
# now deconvolute
imp_resp_model.deconvolute_for_tspans(
    tspan_list=tspan_list_1_dark,
    t_zero_list=t_zero_list_1_dark,
    measurement=ca_dark_day1,
    mol="O2",
    name="CA_hematite_dark_day1_deconvoluted_180um",  # this will
    # automatically save the plots under this name
    export_data=False,
)

# Another way to use the deconvolution is to attach it to the measurement as
# as a calculator. This makes the deconvoluted flux of the deconvolution's molecule
# ´mol´ available to grab or look up with the key `f"n_dot_{mol}-deconvoluted"`
ca_dark_day1.add_calculator(imp_resp_model)

# look up the deconovluted signal:
t, n_dot_O2_decon = ca_dark_day1.grab("n_dot_O2-deconvoluted", tspan_bg=[0, 10])

# plot it together with raw signal:
axes = ca_dark_day1.plot(mol_list=["O2"], tspan_bg=[0, 10], logplot=False)
axes[0].plot(t, n_dot_O2_decon, color="0.5", label="O2 deconvoluted")
axes[0].legend()
