from pathlib import Path
import ixdat
from ixdat import Measurement
from ixdat.techniques.ms import MSInlet

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
    native_cal = ms.siq_gas_flux_calibration(mol="He", mass="M4", tspan=[100, 200])
except Exception as e:
    print(e)  # explains that native ixdat requires an MSInlet to be specifically defined

# ---- Native calibration ----- #
native_cal = ms.gas_flux_calibration(
    inlet=MSInlet(), mol="He", mass="M4", tspan=[100, 200]
)
print(native_cal)  # An ixdat MSCalResult object

# ---- Spectro Inlets calibration ----- #
ixdat.config.plugins.activate_siq()

quant_cal = ms.siq_gas_flux_calibration(mol="He", mass="M4", tspan=[100, 200])
print(quant_cal)  # A CalPoint object of the external package

calibration = ixdat.plugins.siq.Calibration(cal_list=[quant_cal])

print(calibration + quant_cal)
