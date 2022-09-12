from pathlib import Path
import ixdat
from ixdat import Measurement
from ixdat.techniques.ms import MSInlet


PATH_TO_DATA_FOLDER = Path(__file__).parent.parent / "test_data/Zilien version 1"
ms = Measurement.read(
    PATH_TO_DATA_FOLDER / "2022-04-06 16_17_23 full set.tsv", technique="MS"
)  # an MSMeasurement
# ms.plot()

# ---- Spectro Inlets calibration inaccessible without use_si_quant ----- #
print(ixdat.config.plugins.use_si_quant)  # False
try:
    native_cal = ms.gas_flux_calibration(mol="He", mass="M4", tspan=[100, 200])
except Exception as e:
    print(e)  # explains that native ixdat requires an MSInlet to be specifically defined

# ---- Native calibration ----- #
native_cal = MSInlet().gas_flux_calibration(
    mol="He", mass="M4", measurement=ms, tspan=[100, 200]
)
print(native_cal)  # An ixdat MSCalResult object

# ---- Spectro Inlets calibration ----- #
ixdat.config.plugins.use_si_quant = True

quant_cal = ms.gas_flux_calibration(mol="He", mass="M4", tspan=[100, 200])
print(quant_cal)  # A CalPoint object of the external package
