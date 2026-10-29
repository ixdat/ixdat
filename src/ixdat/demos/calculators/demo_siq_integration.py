"""Demo of ixdat's native MS calibration next to the Spectro Inlets (siq) one.

Requires the demo data and the ``spectro_inlets_quantification`` (siq) package.
"""

import ixdat
from ixdat import Measurement
from ixdat.calculators.ms_calculators import MSInlet, MSCalibration
from ixdat.demos import get_demo_data_dir


def load_ms(data_dir=None):
    """Return an MSMeasurement read from a zilien (version 1) .tsv file."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / "zilien_version_1/2022-04-06 16_17_23 full set.tsv", technique="MS"
    )


EXAMPLES = {"ms": load_ms}


def main(show=True):
    """Print the native and the siq calibrations. `show` has no plots to switch."""
    results = {name: load() for name, load in EXAMPLES.items()}
    ms = results["ms"]

    # ---- Spectro Inlets calibration is inaccessible without activating siq
    print(ixdat.plugins.use_siq)
    # FIXME: accessing a siq class through ixdat.plugins.siq gives `None` if siq has
    #   not been activated. An error message would be more appropriate.
    try:
        siq_calculator = ixdat.plugins.siq.Calculator
    except Exception as e:
        print(e)  # should explain that siq has not been activated.
    else:
        print(f"siq Calculator without activating siq: {siq_calculator}")

    # ---- Native calibration
    native_cal = MSCalibration.gas_flux_calibration(
        measurement=ms, inlet=MSInlet(), mol="He", mass="M4", tspan=[100, 200]
    )
    print(native_cal)  # an ixdat MSCalResult
    results["native_cal"] = native_cal

    # ---- Spectro Inlets calibration
    ixdat.plugins.activate_siq()
    quant_cal = ixdat.plugins.siq.Calculator.gas_flux_calibration(
        measurement=ms, mol="He", mass="M4", tspan=[100, 200]
    )
    # siq Calculators only provide series after `set_quantifier` has been called
    quant_cal.set_quantifier(carrier="He", mol_list=["He"], mass_list=["M4"])
    print(quant_cal)  # a CalPoint of the siq package
    results["siq_cal"] = quant_cal

    # Adding quant_cal to itself works, but the resulting calculator provides no
    # series, since it has no quantifier set.
    print(quant_cal + quant_cal)
    return results


if __name__ == "__main__":
    results = main()
