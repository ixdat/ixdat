# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 18:16:13 2024

@author: Søren

This script demonstrates that the challenges of chaining calculators have been solved.
See https://github.com/ixdat/ixdat/issues/183
"""
from pathlib import Path
from matplotlib import pyplot as plt
from ixdat import Measurement
from ixdat.exceptions import SeriesNotFoundError
from ixdat.techniques.ms import MSCalibration


plt.close("all")

data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/ixdat_exports"

ecms = Measurement.read(data_dir / "trimarco2018_fig3_data.csv", reader="ixdat")

try:
    t, n_dot_H2 = ecms.grab("n_dot_H2", tspan=[100, 200])
except SeriesNotFoundError as e:
    print("got error message, as intended:")
    print(e)


# Chaining MS quantification and background subraction

ecms.add_calculator(MSCalibration("H2", "M2", F=0.21))

t, n_dot_H2 = ecms.grab("n_dot_H2", tspan=[100, 200])

ecms.set_bg(mass_list=["M2"], tspan=[20, 30])

n_dot_H2_bg_removed = ecms.grab_for_t("n_dot_H2", t=t)


fig, ax = plt.subplots()

ax.plot(t, n_dot_H2, "b--")
ax.plot(t, n_dot_H2_bg_removed, "b")
ax.set_xlabel("time / [s]")
ax.set_ylabel("flux / [nmol/s]")


n_dot_H2_raw = ecms.grab_for_t("n_dot_H2", t=t, remove_background=False)
ax.plot(t, n_dot_H2_raw, "k--")

M2_bg_removed = ecms.grab_for_t("M2", t=t, remove_background=True)
M2_raw = ecms.grab_for_t("M2", t=t, remove_background=False)

F_bg_implied = (M2_raw - M2_bg_removed) / (n_dot_H2_raw - n_dot_H2_bg_removed)
assert F_bg_implied[0] == 0.21


t, U_raw = ecms.grab("potential")
print(f"max(U_raw) = {max(U_raw)}")
# Finds no calculators,
# but finds that "potential" can be an alias for "raw_potential"

cal1 = ecms.calibrate(RE_vs_RHE=0.72)

# ECMeasurement.calibrate() raises a warning and implements the last RE_vs_RHE:
cal1a = ecms.calibrate(RE_vs_RHE=0.715)

# This raises a Warning because it makes a third calculator responding to "potential":
U = ecms.grab_for_t("potential", t=t)
print(f"max(U) = {max(U)}")
# ECMeasurement.calibrate() adds R_Ohm to the last
cal2 = ecms.calibrate(R_Ohm=100)
print(f"max(U) = {max(U)}")
U_corr = ecms.grab_for_t("potential", t=t)
print(f"max(U_corr) = {max(U_corr)}")
# The "-raw" suffix ensures that no calculators are applied:
U_raw_again = ecms.grab_for_t("potential-raw", t=t)
print(f"max(U_raw_again) = {max(U_raw_again)}")
U_again = ecms.grab_for_t("potential", calculator_list=[cal1a], t=t)

assert max(U_again) == max(U)
assert max(U_corr) != max(U)
assert max(U_raw_again) == max(U_raw)
assert max(U_raw) != max(U)
