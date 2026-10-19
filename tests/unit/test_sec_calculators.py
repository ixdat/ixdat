"""Tests for the spectroelectrochemistry calculators."""

import numpy as np
import pytest

from ixdat.calculators import CALCULATOR_CLASSES
from ixdat.calculators.sec_calculators import DeltaODCalculator
from ixdat.data_series import DataSeries, Field, TimeSeries, ValueSeries
from ixdat.techniques import ECMeasurement
from ixdat.techniques.spectroelectrochemistry import (
    ECOpticalMeasurement,
    OpticalSpectrumSeries,
)

TSTAMP = 1756386155.0
N_SPECTRA = 11
WL = np.array([400.0, 500.0, 600.0])


def make_sec(spectra_type="Intensity", data=None, potential=None):
    """Return an ECOpticalMeasurement with spectra at t = 0, 1, ..., 10 s.

    By default the potential sweeps linearly from 0 to 1 V, and the intensity at
    each wavelength decays as 10**(-k * U), so that dOD = k * U with k = 1, 2, 3.
    """
    t = np.arange(float(N_SPECTRA))
    potential = np.linspace(0, 1, N_SPECTRA) if potential is None else potential
    if data is None:
        data = 10 ** (-np.outer(potential, [1.0, 2.0, 3.0]))
    ec_tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=TSTAMP)
    ec = ECMeasurement(
        name="ec",
        technique="EC",
        tstamp=TSTAMP,
        series_list=[
            ec_tseries,
            ValueSeries("raw_potential", "V", potential, tseries=ec_tseries),
            ValueSeries("raw_current", "mA", np.zeros_like(t), tseries=ec_tseries),
        ],
    )
    opt_tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=TSTAMP)
    field = Field(
        name="intensity",
        unit_name="counts",
        data=data,
        axes_series=[opt_tseries, DataSeries("wavelength", "nm", WL)],
    )
    optical = OpticalSpectrumSeries(
        name="optical",
        technique="Optical",
        tstamp=TSTAMP,
        field=field,
        continuous=True,
        spectra_type=spectra_type,
    )
    return ec + optical


def expected_dOD(potential=None):
    potential = np.linspace(0, 1, N_SPECTRA) if potential is None else potential
    return np.outer(potential, [1.0, 2.0, 3.0])


def test_dOD_calculator_is_built_in_and_registered():
    sec = make_sec()
    assert isinstance(sec, ECOpticalMeasurement)
    assert isinstance(sec.dOD_calculator, DeltaODCalculator)
    assert CALCULATOR_CLASSES["dOD_calculator"] is DeltaODCalculator


def test_dOD_lookup_uses_the_reference_spectrum():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    np.testing.assert_allclose(sec["dOD"].data, expected_dOD())


def test_dOD_lookup_is_cached_until_the_reference_changes():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    first = sec["dOD"]
    assert sec["dOD"] is first
    sec.set_reference_spectrum(t_ref=10)
    assert sec["dOD"] is not first
    np.testing.assert_allclose(sec["dOD"].data, expected_dOD() - expected_dOD()[-1])


def test_intensity_without_reference_raises():
    sec = make_sec()
    with pytest.raises(ValueError, match="reference"):
        sec["dOD"]


def test_calc_dOD_with_explicit_reference_matches_lookup():
    sec = make_sec()
    np.testing.assert_allclose(sec.calc_dOD(index_ref=0).data, expected_dOD())


def test_absorption_with_reference_subtracts():
    absorbance = expected_dOD() + 0.5
    sec = make_sec(spectra_type="Absorption", data=absorbance)
    np.testing.assert_allclose(sec.calc_dOD(index_ref=0).data, expected_dOD())


def test_absorption_without_reference_warns_and_passes_through():
    absorbance = expected_dOD() + 0.5
    sec = make_sec(spectra_type="Absorption", data=absorbance)
    with pytest.warns(UserWarning):
        np.testing.assert_allclose(sec["dOD"].data, absorbance)


def test_transmission_in_percent_without_reference():
    transmission = 100 * 10 ** (-expected_dOD())
    sec = make_sec(spectra_type="Transmission", data=transmission)
    with pytest.warns(UserWarning):
        np.testing.assert_allclose(sec["dOD"].data, expected_dOD())


@pytest.mark.parametrize(
    "step, n_diff",
    [({"idx_step": 2}, 5), ({"tstep": 2.0}, 5), ({"Vstep": 0.2}, 5)],
)
def test_diff_lookup_with_each_step_type(step, n_diff):
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    for attr, value in step.items():
        setattr(sec.dOD_calculator, attr, value)
    sec.dOD_calculator.normalise = False
    diff = sec["dOD>diff"]
    assert diff.data.shape == (n_diff, 3)
    # dOD rises by k * 0.2 V per step of 2 spectra
    np.testing.assert_allclose(diff.data, np.tile([0.2, 0.4, 0.6], (n_diff, 1)))
    np.testing.assert_allclose(diff.axes_series[0].t, [2, 4, 6, 8, 10])


def test_diff_is_normalised_to_max_absolute_value():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    diff = sec.dOD_calculator.calc_diff(sec, idx_step=2, normalise=True)
    np.testing.assert_allclose(diff.data, np.tile([1 / 3, 2 / 3, 1.0], (5, 1)))


def test_diff_normalisation_keeps_sign_of_negative_changes():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=10)
    diff = sec.dOD_calculator.calc_diff(sec, idx_step=2)
    assert np.all(diff.data > 0)
    sec_down = make_sec(potential=np.linspace(1, 0, N_SPECTRA))
    sec_down.set_reference_spectrum(t_ref=0)
    diff_down = sec_down.dOD_calculator.calc_diff(sec_down, idx_step=2)
    np.testing.assert_allclose(diff_down.data, np.tile([-1 / 3, -2 / 3, -1], (5, 1)))


def test_Vstep_follows_a_cathodic_sweep():
    sec = make_sec(potential=np.linspace(1, 0, N_SPECTRA))
    sec.set_reference_spectrum(t_ref=0)
    diff = sec.dOD_calculator.calc_diff(sec, Vstep=0.2, normalise=False)
    np.testing.assert_allclose(diff.data, np.tile([-0.2, -0.4, -0.6], (5, 1)))


def test_Vstep_rejects_non_monotonic_potential():
    potential = np.concatenate([np.linspace(0, 1, 6), np.linspace(0.8, 0, 5)])
    sec = make_sec(potential=potential)
    sec.set_reference_spectrum(t_ref=0)
    with pytest.raises(ValueError, match="monotonic"):
        sec.dOD_calculator.calc_diff(sec, Vstep=0.2)


def test_a_new_step_replaces_the_old_one():
    calc = DeltaODCalculator(idx_step=2)
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    calc.calc_diff(sec, tstep=2.0)
    assert calc.idx_step is None and calc.tstep == 2.0


def test_diff_requires_a_step():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    with pytest.raises(ValueError, match="exactly one"):
        sec["dOD>diff"]


def test_diff_works_on_a_selected_cut():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    cut = sec.cut(tspan=[3.5, 10.5])
    cut.dOD_calculator.idx_step = 1
    cut.dOD_calculator.normalise = False
    np.testing.assert_allclose(cut["dOD>diff"].data, np.tile([0.1, 0.2, 0.3], (6, 1)))


def test_cut_without_reference_works():
    sec = make_sec(spectra_type="Absorption", data=expected_dOD())
    cut = sec.cut(tspan=[3.5, 10.5])
    assert cut.reference_spectrum is None


def test_round_trip_through_as_dict():
    calc = DeltaODCalculator(normalise=False, xspan=[400, 600], Vstep=0.01)
    copy = DeltaODCalculator.from_dict(calc.as_dict())
    assert isinstance(copy, DeltaODCalculator)
    assert repr(copy) == repr(calc)
