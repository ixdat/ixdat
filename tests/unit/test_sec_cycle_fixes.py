"""Tests for cycle definition by turning point and the EC-optical path on real data.

These cover:
- CyclicVoltammogram.redefine_cycle(turning_point=True): counting reversals, only
  upward reversals with redox=True, numbering after start_potential, and the error
  for a sweep without reversals.
- The EC-optical path on the committed fixture files: reading, combining, dOD,
  difference spectra, denoising, and plotting.
"""

from pathlib import Path

import numpy as np
import pytest

from ixdat import Measurement, Spectrum
from ixdat.calculators.sec_calculators import Denoiser
from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.techniques import ECMeasurement


DATA_DIR = Path(__file__).resolve().parents[2] / "test_data" / "oceanview_sec"
OPTICAL_FILE = DATA_DIR / "mini_oceanview__0__15-02-35-123.txt"
EC_FILE = DATA_DIR / "demo_oceanview_ecoptical_sec.mpt"

# Four sweep legs, so three reversals: up at 0-1, down, up, down
FOUR_LEGS = np.concatenate(
    [
        np.linspace(0, 1, 10),
        np.linspace(1, 0, 10),
        np.linspace(0, 1, 10),
        np.linspace(1, 0, 10),
    ]
)


# --- helpers -----------------------------------------------------------


def make_triangular_cv(v, tstamp=1000.0):
    """Return a CyclicVoltammogram with the given potential, one point per second."""
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


def read_ec_optical_fixture():
    """Read the EC and optical fixtures, one anodic sweep, and combine them."""
    optical = Spectrum.read(OPTICAL_FILE, reader="oceanview")
    ec = Measurement.read(EC_FILE, reader="biologic")
    ec_optical = ec.as_cv() + optical
    ec_optical.set_reference_spectrum(t_ref=55)
    return ec_optical


# --- CyclicVoltammogram.redefine_cycle(turning_point=True) -------------
# The synthetic sweeps have only 10 points per leg, so res_points is set to 2.


def test_redefine_cycle_turning_point_counts_reversals():
    """Four sweep legs (3 reversals) give cycles 0 to 3."""
    cv = make_triangular_cv(FOUR_LEGS)
    cv.redefine_cycle(turning_point=True, N_points=3, N_sep=5, res_points=2)

    cycle = cv["cycle"].data
    assert sorted(set(cycle.tolist())) == [0.0, 1.0, 2.0, 3.0]
    # each new cycle starts at a reversal, 10 points apart
    np.testing.assert_allclose(np.diff(np.flatnonzero(np.diff(cycle))), 10, atol=2)


def test_redefine_cycle_turning_point_redox_true_counts_only_upward_turns():
    cv = make_triangular_cv(FOUR_LEGS)
    cv.redefine_cycle(
        turning_point=True, redox=True, N_points=3, N_sep=5, res_points=2
    )

    # only the reversal from down to up, around index 20, starts a new cycle
    cycle = cv["cycle"].data
    assert sorted(set(cycle.tolist())) == [0.0, 1.0]
    assert abs(np.flatnonzero(cycle)[0] - 20) <= 2


def test_redefine_cycle_turning_point_after_start_potential():
    """Cycle 1 starts where the potential first crosses start_potential."""
    cv = make_triangular_cv(FOUR_LEGS)
    cv.redefine_cycle(
        turning_point=True, start_potential=0.5, N_points=3, N_sep=5, res_points=2
    )

    cycle = cv["cycle"].data
    assert cycle[0] == 0
    first_of_cycle_1 = np.flatnonzero(cycle == 1)[0]
    assert cv.U[first_of_cycle_1 - 1] < 0.5 <= cv.U[first_of_cycle_1]
    assert sorted(set(cycle.tolist())) == [0.0, 1.0, 2.0, 3.0, 4.0]


def test_redefine_cycle_turning_point_without_reversal_raises():
    cv = make_triangular_cv(np.linspace(0, 1, 40))
    with pytest.raises(ValueError, match="No valid turning points"):
        cv.redefine_cycle(turning_point=True, res_points=2)


# --- the EC-optical path on the committed fixtures ------------------------


def test_sec_pipeline_on_real_fixture():
    """The fixture is one anodic sweep from 0.38 V to 2 V followed by a hold."""
    ec_optical = read_ec_optical_fixture()
    assert ec_optical.technique == "EC-Optical"

    sweep = ec_optical.select(tspan=[55, 95])
    assert sweep["dOD"].data.shape == sweep.spectra.data.shape

    sweep.dOD_calculator.idx_step = 1
    diff = sweep["dOD>diff"]
    assert diff.data.shape[0] == sweep.spectra.data.shape[0] - 1

    sweep.add_calculator(
        Denoiser(denoise_method="Savitzky-Golay", sg_window=5, sg_poly_order=2)
    )
    assert sweep["dOD>diff>smooth"].data.shape == diff.data.shape


def test_plot_wavelengths_vs_cv_accepts_caller_supplied_axes():
    """Passing axes= explicitly used to raise UnboundLocalError on ax_current."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ec_optical = read_ec_optical_fixture()
    fig, ax = plt.subplots()
    axes = ec_optical.plot_wavelengths_vs_cv(
        wavelengths=["w650"], axes=(ax, ax.twinx())
    )
    assert axes[0] is ax
    plt.close(fig)
