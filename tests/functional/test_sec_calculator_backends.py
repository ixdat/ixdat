"""Saving and loading spectroelectrochemistry calculators on the database backends"""

import sys
from pathlib import Path

import numpy as np

from ixdat import Measurement
from ixdat.calculators.sec_calculators import (
    DeltaODCalculator,
    Denoiser,
    SpectrumDeconvolutor,
)

sys.path.insert(0, str(Path(__file__).parents[1] / "unit"))
from test_sec_calculators import two_transition_sec  # noqa: E402


def test_round_trip_of_sec_calculators(fresh_backend):
    sec, _ = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1, xspan=[450, 850]
    )
    calculators = [
        DeltaODCalculator(normalise=False, xspan=[400, 600], Vstep=0.01),
        Denoiser(denoise_method="Savitzky-Golay", sg_window=7),
        decon,
    ]
    for calculator in calculators:
        sec.add_calculator(calculator)

    loaded = Measurement.get(sec.save())

    assert loaded.spectra_type == "Intensity"
    for calculator in calculators:
        loaded_calculator = loaded.calculators[calculator.calculator_type]
        assert type(loaded_calculator) is type(calculator)
        assert repr(loaded_calculator) == repr(calculator)
    loaded_decon = loaded.calculators["spectrum_deconvolutor"]
    np.testing.assert_allclose(loaded_decon.converged_data, decon.converged_data)
    np.testing.assert_allclose(loaded_decon.wavelength, decon.wavelength)
    assert len(loaded_decon.potentialspans) == len(decon.potentialspans)
