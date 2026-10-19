"""Tests for the spectroelectrochemistry calculators."""

import numpy as np
import pytest

from ixdat.calculators import CALCULATOR_CLASSES
from ixdat.calculators.sec_calculators import (
    DeltaODCalculator,
    Denoiser,
    SpectrumDeconvolutor,
    calc_cycle_difference,
)
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


def noisy_sec(seed=0):
    """Return a sec with many wavelengths: smooth spectra plus noise."""
    rng = np.random.default_rng(seed)
    wl = np.linspace(400, 900, 101)
    smooth = np.outer(np.linspace(0.1, 1, N_SPECTRA), np.exp(-((wl - 650) / 80) ** 2))
    return smooth, smooth + 0.01 * rng.standard_normal(smooth.shape), wl


def field_of(data, wl, name="dOD"):
    tseries = TimeSeries("time", "s", np.arange(float(len(data))), tstamp=TSTAMP)
    return Field(
        name=name,
        unit_name="",
        data=data,
        axes_series=[tseries, DataSeries("wavelength", "nm", wl)],
    )


def test_savitzky_golay_reduces_noise():
    smooth, noisy, wl = noisy_sec()
    denoised = Denoiser(denoise_method="Savitzky-Golay", sg_window=11).denoise_spectra(
        field_of(noisy, wl)
    )
    assert denoised.name == "dOD>smooth"
    assert np.std(denoised.data - smooth) < 0.5 * np.std(noisy - smooth)


def test_PCA_reduces_noise():
    pytest.importorskip("sklearn")
    smooth, noisy, wl = noisy_sec()
    denoised = Denoiser(PCA_explained_variance=0.9).denoise_spectra(
        field_of(noisy, wl)
    )
    assert np.std(denoised.data - smooth) < 0.5 * np.std(noisy - smooth)


def test_unknown_denoise_method_raises():
    _, noisy, wl = noisy_sec()
    with pytest.raises(ValueError, match="denoise_method"):
        Denoiser(denoise_method="magic").denoise_spectra(field_of(noisy, wl))


def test_normalise_scales_to_unit_range():
    _, noisy, wl = noisy_sec()
    denoised = Denoiser(
        denoise_method="Savitzky-Golay", sg_window=11, normalise=True
    ).denoise_spectra(field_of(noisy, wl))
    np.testing.assert_allclose(np.min(denoised.data, axis=1), 0, atol=1e-12)
    np.testing.assert_allclose(np.max(denoised.data, axis=1), 1)


def test_smooth_suffix_chains_with_diff():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    sec.dOD_calculator.idx_step = 2
    sec.add_calculator(
        Denoiser(denoise_method="Savitzky-Golay", sg_window=3, sg_poly_order=1)
    )
    result = sec["dOD>diff>smooth"]
    assert result.data.shape == sec["dOD>diff"].data.shape
    # A linear filter leaves spectra linear in wavelength unchanged
    np.testing.assert_allclose(result.data, sec["dOD>diff"].data)


def test_measurement_denoise_spectra_uses_defaults_for_None():
    _, noisy, wl = noisy_sec()
    sec = make_sec()
    denoised = sec.denoise_spectra(
        spectra_field=field_of(noisy, wl),
        denoise_method="Savitzky-Golay",
        sg_window=None,
        sg_poly_order=None,
    )
    assert denoised.data.shape == noisy.shape


def test_denoiser_round_trip_through_as_dict():
    calc = Denoiser(denoise_method="Savitzky-Golay", sg_window=7, xspan=[400, 800])
    copy = Denoiser.from_dict(calc.as_dict())
    assert isinstance(copy, Denoiser)
    assert repr(copy) == repr(calc)


def two_transition_sec(n=81, noise=1e-4, seed=0, n_dead=0):
    """Return a sec with two redox transitions, at 0.3 V and 0.7 V, and their shapes.

    Each transition has its own Gaussian spectrum. The intensity is referenced to the
    first spectrum, where both transitions are (practically) not started.
    """
    rng = np.random.default_rng(seed)
    wl = np.linspace(400, 900, 101)
    shapes = np.array(
        [np.exp(-((wl - 500) / 40) ** 2), np.exp(-((wl - 750) / 40) ** 2)]
    )
    potential = np.linspace(0, 1, n)
    completion = np.array(
        [1 / (1 + np.exp(-(potential - U0) / 0.03)) for U0 in (0.3, 0.7)]
    )
    dOD = completion.T @ shapes
    intensity = 10 ** (-dOD) * (1 + noise * rng.standard_normal(dOD.shape))
    # The first n_dead wavelengths read negative counts in one spectrum, as dead or
    # over-subtracted pixels do
    intensity[n // 2, :n_dead] = -1
    t = np.arange(float(n))
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
    optical = OpticalSpectrumSeries(
        name="optical",
        technique="Optical",
        tstamp=TSTAMP,
        field=Field(
            name="intensity",
            unit_name="counts",
            data=intensity,
            axes_series=[opt_tseries, DataSeries("wavelength", "nm", wl)],
        ),
        continuous=True,
        spectra_type="Intensity",
    )
    sec = ec + optical
    sec.set_reference_spectrum(t_ref=0)
    return sec, shapes


def test_deconvolutor_finds_both_converged_spectra():
    sec, shapes = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1
    )
    assert len(decon.converged_data) == 2
    for converged, shape in zip(decon.converged_data, shapes):
        assert np.corrcoef(converged, shape)[0, 1] > 0.999
    lo, hi = decon.potentialspans
    assert abs(np.mean(lo) - 0.3) < 0.05 and abs(np.mean(hi) - 0.7) < 0.05
    assert len(decon.adjacent_distances) == 79
    assert len(decon.spectrum_list) == 2


def test_deconvolutor_lookup_reconstructs_the_dOD():
    sec, _ = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1
    )
    sec.add_calculator(decon)
    coefficient_0 = sec["coefficient_0"]
    assert isinstance(coefficient_0, ValueSeries)
    assert coefficient_0.data.shape == (81,)
    # The residuals are in mOD, on a signal of up to 1000 mOD
    assert np.max(np.abs(sec["residuals"].data)) < 5
    np.testing.assert_allclose(sec["reconstructed"].data, sec["dOD"].data, atol=5e-3)
    # The first transition completes before the second starts
    assert sec["coefficient_0"].data[40] > 10 * sec["coefficient_1"].data[40]


def test_deconvolutor_fits_a_different_measurement():
    sec, _ = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1
    )
    other, _ = two_transition_sec(n=41)
    coefficients = decon.deconvolute(other)
    assert len(coefficients) == 2 and coefficients[0].data.shape == (41,)


def test_deconvolutor_fits_on_its_own_wavelengths():
    sec, _ = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1, xspan=[450, 850]
    )
    assert decon.wavelength.min() >= 450 and decon.wavelength.max() <= 850
    decon.deconvolute(sec)
    assert decon.reconstructed_field.data.shape == (81, len(decon.wavelength))


def test_deconvolutor_bounds_keep_coefficients_increasing():
    sec, _ = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1, bounds_bool=True
    )
    for coefficient in decon.deconvolute(sec):
        assert np.all(np.diff(coefficient.data) >= -1e-9)


def test_no_convergence_raises():
    sec, _ = two_transition_sec()
    with pytest.raises(ValueError, match="No convergent regions"):
        SpectrumDeconvolutor.from_spectromeasurement(sec, idx_step=1, conv_limit=0)


def test_deconvolutor_round_trip_through_as_dict():
    sec, _ = two_transition_sec()
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sec, idx_step=1, min_region_width=1
    )
    copy = SpectrumDeconvolutor.from_dict(decon.as_dict())
    np.testing.assert_allclose(copy.converged_data, decon.converged_data)
    np.testing.assert_allclose(copy.wavelength, decon.wavelength)
    assert copy.potentialspans == decon.potentialspans


@pytest.fixture
def deconvoluted_sec():
    import matplotlib

    matplotlib.use("Agg")
    sec, _ = two_transition_sec()
    sec.dOD_calculator.idx_step = 1
    decon = SpectrumDeconvolutor.from_spectromeasurement(sec, min_region_width=1)
    sec.add_calculator(decon)
    return sec


def test_plot_dOD_difference_spectra(deconvoluted_sec):
    ax = deconvoluted_sec.plotter.plot_dOD_difference_spectra()
    assert len(ax.get_lines()) == 80


def test_plot_converged_spectra(deconvoluted_sec):
    ax = deconvoluted_sec.plotter.plot_converged_spectra()
    assert len(ax.get_lines()) == 2


def test_plot_fit_and_residuals(deconvoluted_sec):
    axes = deconvoluted_sec.plotter.plot_fit_and_residuals()
    assert len(axes[0].get_lines()) == 2
    assert len(axes[1].get_lines()) == 81


def test_plot_fit_reconstruction(deconvoluted_sec):
    ax = deconvoluted_sec.plotter.plot_fit_reconstruction()
    assert len(ax.get_lines()) == 81


def test_plotting_without_deconvolutor_explains_what_to_do():
    sec, _ = two_transition_sec()
    with pytest.raises(ValueError, match="add_calculator"):
        sec.plotter.plot_converged_spectra()


def two_cycle_sec(drift=0.0):
    """Return a sec with two triangle-wave cycles of 21 spectra each.

    The dOD follows the potential, plus `drift` OD in the second cycle.
    """
    sweep = np.concatenate([np.linspace(0, 1, 11), np.linspace(0.9, 0, 10)])
    potential = np.concatenate([sweep, sweep])
    cycle = np.repeat([1.0, 2.0], len(sweep))
    dOD = np.outer(potential + drift * (cycle - 1), [1.0, 2.0, 3.0])
    t = np.arange(float(len(potential)))
    ec_tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=TSTAMP)
    ec = ECMeasurement(
        name="ec",
        technique="EC",
        tstamp=TSTAMP,
        series_list=[
            ec_tseries,
            ValueSeries("raw_potential", "V", potential, tseries=ec_tseries),
            ValueSeries("raw_current", "mA", np.zeros_like(t), tseries=ec_tseries),
            ValueSeries("cycle", "", cycle, tseries=ec_tseries),
        ],
    )
    opt_tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=TSTAMP)
    optical = OpticalSpectrumSeries(
        name="optical",
        technique="Optical",
        tstamp=TSTAMP,
        field=Field(
            name="intensity",
            unit_name="counts",
            data=10 ** (-dOD),
            axes_series=[opt_tseries, DataSeries("wavelength", "nm", WL)],
        ),
        continuous=True,
        spectra_type="Intensity",
    )
    sec = ec + optical
    sec.set_reference_spectrum(t_ref=0)
    return sec


def test_identical_cycles_give_zero_difference():
    diff = calc_cycle_difference(two_cycle_sec())
    assert diff.unit_name == "mOD"
    assert diff.data.shape[1] == 3
    np.testing.assert_allclose(diff.data, 0, atol=1e-9)


def test_cycle_difference_shows_drift_in_mOD():
    diff = calc_cycle_difference(two_cycle_sec(drift=0.01))
    expected = np.tile([10.0, 20.0, 30.0], (len(diff.data), 1))
    np.testing.assert_allclose(diff.data, expected)


def test_plot_dOD_cycle_diff():
    import matplotlib

    matplotlib.use("Agg")
    sec = two_cycle_sec(drift=0.01)
    ax = sec.plotter.plot_dOD_cycle_diff(
        upper_noise_bound=1, lower_noise_bound=-1
    )
    assert len(ax.get_lines()) == len(calc_cycle_difference(sec).data)


@pytest.mark.parametrize("width", [0, 150])
def test_track_wavelength_matches_dOD(width):
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    w500 = sec.track_wavelength(500, width=width)
    np.testing.assert_allclose(w500.data, expected_dOD()[:, 1])
    assert "w500" in sec.tracked_wavelengths
    np.testing.assert_allclose(sec["w500"].data, w500.data)


def test_track_wavelength_with_explicit_reference_on_absorption():
    absorbance = expected_dOD() + 0.5
    sec = make_sec(spectra_type="Absorption", data=absorbance)
    w600 = sec.track_wavelength(600, width=0, index_ref=0)
    np.testing.assert_allclose(w600.data, expected_dOD()[:, 2])


def test_denoiser_leaves_nan_wavelengths_alone():
    smooth, noisy, wl = noisy_sec()
    noisy[3, :5] = np.nan
    denoised = Denoiser(denoise_method="Savitzky-Golay", sg_window=11).denoise_spectra(
        field_of(noisy, wl)
    )
    assert np.all(np.ma.getmaskarray(denoised.data)[:, :5])
    assert not np.any(np.ma.getmaskarray(denoised.data)[:, 5:])


def test_step_smaller_than_spacing_warns():
    sec = make_sec()
    sec.set_reference_spectrum(t_ref=0)
    with pytest.warns(UserWarning, match="smaller than the spacing"):
        sec.dOD_calculator.calc_diff(sec, Vstep=0.05)


def test_deconvolutor_leaves_out_nan_wavelengths():
    sec, _ = two_transition_sec(n_dead=10)
    sec.dOD_calculator.idx_step = 1
    decon = SpectrumDeconvolutor.from_spectromeasurement(sec, min_region_width=1)
    assert len(decon.wavelength) == 91
    sec.add_calculator(decon)
    assert not np.any(np.ma.getmaskarray(sec["residuals"].data))
    assert np.max(np.abs(sec["residuals"].data)) < 5
