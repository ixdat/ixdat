"""For use in development of the cinfdata reader. Requires access to sample data."""

from pathlib import Path
from ixdat import Measurement, plugins
from ixdat.techniques.ms import MSCalResult, MSCalibration

data_dir = (
    Path.home() / "Dropbox/ixdat_resources/tutorials_data/"
    "23J02_ec_ms_quantification/Zenodo_8400063"
)
ms = Measurement.read(
    data_dir / "2022-07-19 10_02_32 HER_OER_calibration.tsv", reader="zilien"
)
# ms_meas.plot_measurement()
ec = Measurement.read_set(data_dir / "HER_OER", suffix=".mpt")
ec.calibrate(RE_vs_RHE=0, A_el=0.196)

ecms = ec + ms

ecms.tstamp += 9000
axes = ecms.plot(tspan=[0, 3000])

# this has background subtraction and goes through zero:
cal_result = ecms.ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan_list=[[600, 700], [1150, 1250], [1800, 1900], [2350, 2450]],
    tspan_bg=[400, 450],
    axes_measurement=axes,
    force_through_zero=True,
)
# this is forced through zero but without background subtraction gives the wrong answer:
cal_result_2 = ecms.ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan_list=[[600, 700], [1150, 1250], [1800, 1900], [2350, 2450]],
    force_through_zero=True,
)
# this doesn't go through zero but still gives the right answer:
cal_result_3 = ecms.ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan_list=[[600, 700], [1150, 1250], [1800, 1900], [2350, 2450]],
)

plugins.activate_siq()

# The following issues a warning and returns an ixdat MSCalResult :)
cal_1 = ecms.ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan_list=[[600, 700], [1150, 1250], [1800, 1900], [2350, 2450]],
)
# The following returns a siq CalPoint :)
siq_cal_2 = ecms.siq_ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan_list=[[600, 700], [1150, 1250], [1800, 1900], [2350, 2450]],
    tspan_bg=[400, 450],
    force_through_zero=True,
)

# And here we demonstrate all the interconversions between
# siq and native ixdat calibration objects:
# (this is not a natural workflow, just some code to show that the methods work.)

cal_2 = MSCalResult.from_siq(siq_cal_2)
print(cal_2)

siq_cal_1 = cal_1.to_siq()
print(siq_cal_1)

# You can't directly add MSCalResults, but that operation *is* implemented in siq :)
siq_calibration = cal_1.to_siq() + cal_2.to_siq()
print(siq_calibration)

# The following dowsn't work because it's a SensitivityList instead of a Calibration :(
# siq_calibration.plot_as_spectrum()

calibration = MSCalibration.from_siq(siq_calibration)
print(calibration)

siq_calibration_again = calibration.to_siq()

siq_calibration_again.plot_as_spectrum()  # works! :)
