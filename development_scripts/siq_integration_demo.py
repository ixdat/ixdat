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


try:
    siqCalculator = ixdat.config.plugins.siq.Calculator  # gives an error
except Exception as e:
    print(e)  # explains that native ixdat requires an MSInlet to be specifically defined

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
quant_cal.set_quantifier(carrier="He", mol_list=["He"], mass_list=["M4"])
print(quant_cal)  # A CalPoint object of the external package

print(quant_cal + quant_cal)
