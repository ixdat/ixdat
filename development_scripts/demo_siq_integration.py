from pathlib import Path
import ixdat
from ixdat import Measurement
from ixdat.calculators.ms_calculators import MSInlet, MSCalibration

PATH_TO_DATA_FOLDER = (
    Path(__file__).parent.parent / "submodules/ixdat-large-test-files/zilien_version_1"
)

ms = Measurement.read(
    PATH_TO_DATA_FOLDER / "2022-04-06 16_17_23 full set.tsv", technique="MS"
)  # an MSMeasurement
# ms.plot()

# ---- Spectro Inlets calibration inaccessible without activating siq ----- #
print(ixdat.config.plugins.use_siq)  # False


# FIXME: now, accessing a siq class through ixdat.config.plugins.siq gives `None` if
#   siq has not been activated. An error message would be more appropriate.
try:
    siqCalculator0 = ixdat.config.plugins.siq.Calculator
except Exception as e:
    print(e)  # should explain that siq has not been activated.
else:
    print("siqCalculator0 = " + str(siqCalculator0))
    print("No error when getting a siq Calculator class despite not activating siq!")


# ---- Native calibration ----- #
native_cal = MSCalibration.gas_flux_calibration(
    measurement=ms,
    inlet=MSInlet(),
    mol="He",
    mass="M4",
    tspan=[100, 200],
)
print(native_cal)  # An ixdat MSCalResult object

# ---- Spectro Inlets calibration ----- #
ixdat.config.plugins.activate_siq()
siqCalculator = ixdat.config.plugins.siq.Calculator

quant_cal = siqCalculator.gas_flux_calibration(
    measurement=ms, mol="He", mass="M4", tspan=[100, 200]
)
# siqCalculator objects require their method `set_quantifier` to be called before
# they provide
quant_cal.set_quantifier(carrier="He", mol_list=["He"], mass_list=["M4"])
print(quant_cal)  # A CalPoint object of the external package

# The following works but the calculator that results from adding quant_cal to itself
# does not provide any series because it has not yet had a quantifier set
print(quant_cal + quant_cal)
