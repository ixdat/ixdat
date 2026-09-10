"""Demonstrate OceanView EC-optical data in ixdat.

The script reads a small OceanView optical spectrum-series fixture and a BioLogic EC
fixture, combines them into an ECOpticalMeasurement, and shows the basic EC-optical
plots. It then reads a real multi-cycle CV and pairs it with synthesized optical data
(no real spectrometer recording exists for that CV) to show the cycle-based dOD
analysis -- turning-point cycles, anodic/cathodic splitting, denoising, and
convergent-spectra detection -- which needs more than one full sweep to demonstrate.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ixdat import Measurement, Spectrum
from ixdat.data_series import DataSeries, Field, TimeSeries
from ixdat.techniques.spectroelectrochemistry import OpticalSpectrumSeries


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DATA_DIR = REPO_ROOT / "test_data" / "oceanview_sec"
BIOLOGIC_DIR = REPO_ROOT / "test_data" / "biologic"

OPTICAL_FILE = DATA_DIR / "mini_oceanview__0__15-02-35-123.txt"
EC_FILE = DATA_DIR / "demo_oceanview_ecoptical_sec.mpt"
MULTI_CYCLE_CV_FILE = BIOLOGIC_DIR / "Pt_poly_cv.mpt"


def title_current_figure(title):
    """Title the most recently created figure, so its window is identifiable."""
    plt.gcf().suptitle(title)


def read_ec_optical_measurement():
    """Read the EC and optical files and combine them into ECOpticalMeasurement."""
    optical = Spectrum.read(OPTICAL_FILE, reader="oceanview")
    ec = Measurement.read(EC_FILE, reader="biologic")
    ec_optical = ec + optical
    ec_optical.set_reference_spectrum(t_ref=55)
    return ec_optical


def demo_average_every():
    """Show how average_every reduces the reader's output without changing it."""
    unaveraged = Spectrum.read(OPTICAL_FILE, reader="oceanview")
    averaged = Spectrum.read(OPTICAL_FILE, reader="oceanview", average_every=4)
    print(
        f"average_every=1 (default): {unaveraged.field.data.shape[0]} spectra; "
        f"average_every=4: {averaged.field.data.shape[0]} spectra"
    )


def synthesize_optical_series(cv, tspan, tau=3.0, seed=0):
    """Fabricate an OpticalSpectrumSeries whose absorption peak relaxes toward a
    potential-dependent equilibrium, for the cycle-based demo below.

    No real spectrometer recording exists for MULTI_CYCLE_CV_FILE, so this makes
    up a plausible signal: a Gaussian peak whose height lags the CV's real
    potential with relaxation time `tau`, the way many optical signals actually
    behave after a potential step. That relaxation is what get_convergent_spectra
    is designed to detect, so a signal that tracked potential instantaneously
    would never show a convergent region.
    """
    t_full, v_full = cv.grab(cv.U_name)
    t = np.arange(tspan[0], tspan[1], 1.0)
    v = np.interp(t, t_full, v_full)
    wavelengths = np.arange(400.0, 901.0, 50.0)

    equilibrium_height = 0.3 * (v - v_full.min()) / (v_full.max() - v_full.min())
    peak_height = np.zeros_like(equilibrium_height)
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        peak_height[i] = peak_height[i - 1] + (
            equilibrium_height[i] - peak_height[i - 1]
        ) * (1 - np.exp(-dt / tau))

    rng = np.random.default_rng(seed)
    peak_center, peak_width = 550, 40
    spectra = [
        2.0
        - h * np.exp(-0.5 * ((wavelengths - peak_center) / peak_width) ** 2)
        + rng.normal(scale=0.0005, size=len(wavelengths))
        for h in peak_height
    ]

    tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=cv.tstamp)
    wl_series = DataSeries(name="wavelength", unit_name="nm", data=wavelengths)
    field = Field(
        name="intensity",
        unit_name="a.u.",
        data=np.array(spectra),
        axes_series=[tseries, wl_series],
    )
    return OpticalSpectrumSeries(
        name="synthesized optical (demo only, see synthesize_optical_series)",
        reader=None,
        technique="Optical",
        tstamp=cv.tstamp,
        field=field,
        continuous=True,
        spectra_type=None,
    )


def demo_cycle_based_analysis():
    """Demonstrate cycle-based SEC analysis on a real multi-cycle CV.

    Cycle 1 of this CV is a full anodic-then-cathodic sweep (one turning point in
    the middle), so it's enough to show the direction split as well as comparing
    two whole cycles to each other.
    """
    print("reading a real 10-cycle CV")
    ec = Measurement.read(MULTI_CYCLE_CV_FILE, reader="biologic")
    cv = ec.as_cv()
    cv.redefine_cycle(turning_point=True, redox=True, N_points=5, N_sep=10)

    tspan = (12.4, 79.8)  # cycles 1 and 2
    print(f"synthesizing optical data over t={tspan} (cycles 1 and 2)")
    optical = synthesize_optical_series(cv, tspan=tspan)
    ec_optical = cv + optical
    ec_optical.set_reference_spectrum(t_ref=tspan[0])

    print("plotting the anodic and cathodic waterfalls for cycle 1")
    ec_optical.plot_waterfall_cycle(cycle_number=1, direction=0)
    title_current_figure("Cycle 1, anodic")
    ec_optical.plot_waterfall_cycle(cycle_number=1, direction=1)
    title_current_figure("Cycle 1, cathodic")

    print("plotting the difference between cycle 1 and cycle 2 (anodic)")
    ec_optical.plot_dOD_cycle_diff(cycle_1=1, cycle_2=2, direction=0)
    title_current_figure("Cycle 1 vs cycle 2, anodic difference")

    print("tracking a wavelength vs potential across cycles 1-2, current overlaid")
    ec_optical.plot_wavelengths_vs_cv(wavelengths=["w550"], tspan=tspan)
    title_current_figure("Tracked wavelength vs potential (cycles 1-2)")

    print("getting and plotting the dOD difference spectra for cycle 1 (anodic)")
    diff_spectra = ec_optical.get_dOD_difference_spectra(cycle_number=1, direction=0)
    ec_optical.plot_dOD_difference_spectra(diff_spectra=diff_spectra)
    title_current_figure("dOD difference spectra, cycle 1 anodic")

    print("denoising the difference spectra")
    denoised = ec_optical.denoise_spectra(
        spectra_field=diff_spectra,
        denoise_method="Savitzky-Golay",
        sg_window=5,
        sg_poly_order=2,
    )
    ec_optical.plotter.plot_waterfall_vs(
        measurement=ec_optical, field=denoised, vs=ec_optical.U_name
    )
    title_current_figure("Denoised dOD difference spectra, cycle 1 anodic")

    print("getting and plotting the convergent spectra for cycle 1 (anodic)")
    # conv_limit is looser than the 0.01 default because this demo's optical data
    # is only sampled once a second -- coarser than a real spectrometer.
    converge = ec_optical.get_convergent_spectra(
        cycle_number=1, direction=0, conv_limit=0.05, min_region_width=1
    )
    ec_optical.plot_convergent_spectra(converge_output=converge)
    title_current_figure("Convergent spectra, cycle 1 anodic")


def main():
    demo_average_every()

    print("reading EC and optical data")
    ec_optical = read_ec_optical_measurement()
    print(f"combined measurement technique: {ec_optical.technique}")
    print(
        "optical spectra: "
        f"{len(ec_optical.spectrum_series)} spectra x {len(ec_optical.wavelength.data)} wavelengths"
    )

    print("plotting EC-optical heat map")
    ec_optical.plot_measurement(tspan=[55, 95], wlspan=[400, 900], t_ref=55)
    title_current_figure("EC-optical heat map")

    print("plotting waterfall")
    ec_optical.plot_waterfall(t_ref=55)
    title_current_figure("Waterfall")

    print("plotting wavelength tracking")
    ec_optical.plot_wavelengths(wavelengths=["w650", "w800"], tspan=[55, 95])
    title_current_figure("Wavelength tracking vs time")

    print("plotting wavelength tracking vs potential")
    ec_optical.plot_wavelengths_vs_potential(
        wavelengths=["w650", "w800"], tspan=[55, 95]
    )
    title_current_figure("Wavelength tracking vs potential")

    print("plotting wavelength tracking vs potential, with current overlaid")
    ec_optical.plot_wavelengths_vs_cv(wavelengths=["w650", "w800"], tspan=[55, 95])
    title_current_figure("Wavelength tracking vs potential, with current")

    demo_cycle_based_analysis()

    plt.show()


if __name__ == "__main__":
    main()
