"""Tests for the cycle/direction-handling fixes made while reviewing PR #3.

These cover:
- CyclicVoltammogram.redefine_cycle(turning_point=True), including the
  same_sign NameError that hit the default redox=None case with more than
  one turning point.
- ECOpticalMeasurement._check_direction, shared by get_dOD_cycle and
  get_dOD_difference_spectra.
- get_dOD_cycle's two crash fixes: a cycle with no turning point now raises
  the intended ValueError instead of NameError, and anodic_mask/cathodic_mask
  no longer go unbound when the sample right before a turning point has an
  exactly-zero gradient.
- The end-to-end path PR #3 depends on (as_cv -> redefine_cycle -> + optical)
  on real fixture data: get_dOD_cycle, get_dOD_difference_spectra,
  denoise_spectra, get_convergent_spectra, plot_convergent_spectra.
"""

from pathlib import Path

import numpy as np
import pytest

from ixdat import Measurement, Spectrum
from ixdat.data_series import DataSeries, Field, TimeSeries, ValueSeries
from ixdat.techniques import ECMeasurement
from ixdat.techniques.spectroelectrochemistry import ECOpticalMeasurement
from ixdat.techniques.spectroelectrochemistry import OpticalSpectrumSeries


DATA_DIR = Path(__file__).resolve().parents[2] / "test_data" / "oceanview_sec"
OPTICAL_FILE = DATA_DIR / "mini_oceanview__0__15-02-35-123.txt"
EC_FILE = DATA_DIR / "demo_oceanview_ecoptical_sec.mpt"


# --- helpers -----------------------------------------------------------


def make_triangular_cv(v, tstamp=1000.0):
    """Return a CyclicVoltammogram whose potential is the given array."""
    t = np.arange(len(v), dtype=float)
    tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=tstamp)
    potential = ValueSeries(
        name="raw_potential", unit_name="V", data=np.asarray(v), tseries=tseries
    )
    current = ValueSeries(
        name="raw_current", unit_name="mA", data=np.ones(len(v)), tseries=tseries
    )
    ec = ECMeasurement(
        name="synthetic EC",
        technique="EC",
        series_list=[tseries, potential, current],
        tstamp=tstamp,
    )
    return ec.as_cv()


def make_ec_optical_with_cycle(v, cycle, tstamp=1000.0, n_wavelengths=2):
    """Return an ECOpticalMeasurement with a hand-set 'cycle' series.

    Every point belongs to a single, arbitrary spectrum so the shapes line
    up; this is for testing get_dOD_cycle's masking logic directly, not for
    representing a realistic acquisition.
    """
    t = np.arange(len(v), dtype=float)
    tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=tstamp)
    potential = ValueSeries(
        name="raw_potential", unit_name="V", data=np.asarray(v), tseries=tseries
    )
    current = ValueSeries(
        name="raw_current", unit_name="mA", data=np.ones(len(v)), tseries=tseries
    )
    ec = ECMeasurement(
        name="synthetic EC",
        technique="EC",
        series_list=[tseries, potential, current],
        tstamp=tstamp,
    )

    wl = DataSeries(
        name="wavelength",
        unit_name="nm",
        data=np.arange(400.0, 400.0 + n_wavelengths * 100, 100),
    )
    spectra = np.random.default_rng(0).random((len(v), n_wavelengths)) + 1
    spec_field = Field(
        name="intensity", unit_name="a.u.", data=spectra, axes_series=[tseries, wl]
    )
    optical = OpticalSpectrumSeries(
        name="synthetic optical",
        reader=None,
        technique="Optical",
        tstamp=tstamp,
        field=spec_field,
        continuous=True,
        spectra_type=None,
    )

    ec_optical = ec + optical
    ec_optical.set_reference_spectrum(t_ref=0)
    cycle_series = ValueSeries(
        name="cycle",
        unit_name="",
        data=np.asarray(cycle, dtype=float),
        tseries=ec_optical.potential.tseries,
    )
    ec_optical.replace_series("cycle", cycle_series)
    return ec_optical


def read_ec_optical_fixture():
    """Read the real EC + optical fixtures and define cycles by turning point."""
    optical = Spectrum.read(OPTICAL_FILE, reader="oceanview")
    ec = Measurement.read(EC_FILE, reader="biologic")
    cv = ec.as_cv()
    cv.redefine_cycle(turning_point=True, redox=None, N_points=5, N_sep=10)
    ec_optical = cv + optical
    ec_optical.set_reference_spectrum(t_ref=55)
    return ec_optical


# --- CyclicVoltammogram.redefine_cycle(turning_point=True) -------------


def test_redefine_cycle_turning_point_counts_reversals():
    """Four sweep legs (3 reversals) should register as three extra cycles."""
    v = np.concatenate(
        [
            np.linspace(0, 1, 10),
            np.linspace(1, 0, 10),
            np.linspace(0, 1, 10),
            np.linspace(1, 0, 10),
        ]
    )
    cv = make_triangular_cv(v)
    cv.redefine_cycle(turning_point=True, redox=None, N_points=3, N_sep=5)

    assert sorted(set(cv["cycle"].data.tolist())) == [0.0, 1.0, 2.0, 3.0]


def test_redefine_cycle_turning_point_redox_none_default_does_not_crash():
    """Regression test: redox=None (the default) with multiple turning points
    used to raise NameError on an unassigned same_sign."""
    v = np.concatenate(
        [
            np.linspace(0, 1, 10),
            np.linspace(1, 0, 10),
            np.linspace(0, 1, 10),
        ]
    )
    cv = make_triangular_cv(v)
    # redox left at its default (None) -- this used to crash.
    cv.redefine_cycle(turning_point=True, N_points=3, N_sep=5)

    assert len(set(cv["cycle"].data.tolist())) > 1


def test_redefine_cycle_turning_point_redox_true_counts_only_upward_turns():
    v = np.concatenate(
        [
            np.linspace(0, 1, 10),
            np.linspace(1, 0, 10),
            np.linspace(0, 1, 10),
            np.linspace(1, 0, 10),
        ]
    )
    cv = make_triangular_cv(v)
    cv.redefine_cycle(turning_point=True, redox=True, N_points=3, N_sep=5)

    # only the negative-to-positive reversal at index 20 counts
    assert sorted(set(cv["cycle"].data.tolist())) == [0.0, 1.0]


# --- ECOpticalMeasurement._check_direction ------------------------------


@pytest.mark.parametrize("direction", [0, 1, None])
def test_check_direction_accepts_valid_values(direction):
    ECOpticalMeasurement._check_direction(direction)  # must not raise


@pytest.mark.parametrize("direction", [2, -1, "anodic"])
def test_check_direction_rejects_invalid_values(direction):
    with pytest.raises(ValueError, match="direction must be"):
        ECOpticalMeasurement._check_direction(direction)


# --- get_dOD_cycle crash fixes -------------------------------------------


def test_get_dOD_cycle_no_turning_point_raises_value_error():
    """A cycle whose potential never reverses has no turning point.

    This used to raise NameError (valid_indices was never assigned);
    it must raise the documented ValueError instead.
    """
    v = np.arange(10.0)  # monotonic: zero turning points
    ec_optical = make_ec_optical_with_cycle(v, cycle=np.zeros(10))

    with pytest.raises(ValueError, match="No valid turning point"):
        ec_optical.get_dOD_cycle(cycle_number=0, direction=0, N_points=2)


def test_get_dOD_cycle_handles_zero_gradient_before_turning_point():
    """A flat sample right before a turning point used to leave
    anodic_mask/cathodic_mask unbound (NameError)."""
    v = np.array([0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 2.0, 1.0, 0.0, -1.0])
    # dUdt via central differences is exactly 0 one sample before the turn.
    assert np.gradient(v)[4] == 0.0
    ec_optical = make_ec_optical_with_cycle(v, cycle=np.zeros(10))

    anodic = ec_optical.get_dOD_cycle(cycle_number=0, direction=0, N_points=2)
    cathodic = ec_optical.get_dOD_cycle(cycle_number=0, direction=1, N_points=2)

    assert anodic.data.shape[0] > 0
    assert cathodic.data.shape[0] > 0


def test_get_dOD_cycle_rejects_invalid_direction():
    v = np.array([0.0, 1.0, 2.0, 1.0, 0.0])
    ec_optical = make_ec_optical_with_cycle(v, cycle=np.zeros(5))

    with pytest.raises(ValueError, match="direction must be"):
        ec_optical.get_dOD_cycle(cycle_number=0, direction=5)


@pytest.mark.parametrize(
    "call",
    [
        lambda m: m.get_dOD_difference_spectra(cycle_number=0, direction=0),
        lambda m: m.get_dOD_cycle_noise(cycle_1=0, cycle_2=0, direction=0),
    ],
)
def test_direction_split_does_not_need_an_explicit_n_points(call):
    """Regression test: get_dOD_difference_spectra and get_dOD_cycle_noise used
    to default N_points to None and forward it straight into get_dOD_cycle,
    which raised TypeError ('<' not supported between int and NoneType) the
    moment direction was 0 or 1 without an explicit N_points. The array needs
    to be long enough either side of the turning point for get_dOD_cycle's own
    default (N_points=10) to find it.
    """
    v = np.concatenate([np.linspace(0, 1, 15), np.linspace(1, 0, 15)])
    ec_optical = make_ec_optical_with_cycle(v, cycle=np.zeros(30))

    call(ec_optical)  # must not raise TypeError


# --- end-to-end path on real fixture data ---------------------------------


def test_sec_cycle_pipeline_end_to_end_on_real_fixture():
    """Exercises the as_cv -> redefine_cycle -> + optical path PR #3 needs,
    and the functions fixed while reviewing it: get_dOD_cycle,
    get_dOD_difference_spectra, denoise_spectra, get_convergent_spectra,
    and plot_convergent_spectra (which used to always crash on
    smooth_spectra)."""
    ec_optical = read_ec_optical_fixture()

    dod_cycle = ec_optical.get_dOD_cycle(cycle_number=0)
    assert dod_cycle.data.shape[0] > 0

    diff = ec_optical.get_dOD_difference_spectra(cycle_number=0)
    assert diff.data.shape[0] > 0

    denoised = ec_optical.denoise_spectra(
        spectra_field=diff,
        denoise_method="Savitzky-Golay",
        sg_window=5,
        sg_poly_order=2,
    )
    assert denoised.data.shape == diff.data.shape

    converge = ec_optical.get_convergent_spectra(
        cycle_number=0, conv_limit=0.9, min_region_width=1
    )
    axes = ec_optical.plot_convergent_spectra(converge_output=converge)
    assert axes is not None


def test_plot_convergent_spectra_no_longer_takes_smooth_spectra():
    """smooth_spectra was removed; passing it must fail loudly, not be
    silently ignored, and the default call must not raise TypeError."""
    ec_optical = read_ec_optical_fixture()
    converge = ec_optical.get_convergent_spectra(
        cycle_number=0, conv_limit=0.9, min_region_width=1
    )

    with pytest.raises(TypeError):
        ec_optical.plot_convergent_spectra(converge_output=converge, smooth_spectra=5)


def test_get_convergent_spectra_handles_no_convergent_region():
    """Regression test: with normalise=True (the default), zero convergent
    regions used to raise AxisError (np.min(empty_array, axis=1)) instead of
    returning an empty result."""
    ec_optical = read_ec_optical_fixture()

    converge = ec_optical.get_convergent_spectra(
        cycle_number=0, conv_limit=-1  # impossible to satisfy -> zero regions
    )

    assert len(converge[1]) == 0
