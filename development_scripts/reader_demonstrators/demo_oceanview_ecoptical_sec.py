"""Demonstrate OceanView EC-optical data in ixdat.

With no further setup, the script reads a small OceanView optical spectrum-series
fixture and a BioLogic EC fixture, combines them into an ECOpticalMeasurement, and
shows the basic EC-optical plots. This part is self-contained and needs no external
data.

Set IXDAT_DEMO_REAL_DATA_DIR to the unzipped IrOx_SEC_benchmark.zip (a real
multi-cycle EC-optical recording on IrOx, too large to commit to the repo)
to instead run the full analysis on one anodic sweep, as subplots of one figure:
one cycle per sweep, the dOD waterfall, smoothed difference spectra, converged
spectra and the fit, and wavelength tracking.
"""

import os
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from ixdat import Measurement, Spectrum
from ixdat.calculators.sec_calculators import Denoiser, SpectrumDeconvolutor


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DATA_DIR = REPO_ROOT / "test_data" / "oceanview_sec"

OPTICAL_FILE = DATA_DIR / "mini_oceanview__0__15-02-35-123.txt"
EC_FILE = DATA_DIR / "demo_oceanview_ecoptical_sec.mpt"

REAL_DATA_DIR = os.environ.get("IXDAT_DEMO_REAL_DATA_DIR")


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


def demo_small_fixture():
    """The basic EC-optical demo: real, small, committed fixture data."""
    demo_average_every()

    print("reading EC and optical data")
    ec_optical = read_ec_optical_measurement()
    print(f"combined measurement technique: {ec_optical.technique}")
    print(
        "optical spectra: "
        f"{len(ec_optical.spectrum_series)} spectra x "
        f"{len(ec_optical.wavelength.data)} wavelengths"
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


def read_real_measurement(data_dir):
    """Read a real multi-cycle CV on IrOx and its optical spectra.

    Expects, as unzipped from IrOx_SEC_benchmark.zip:
    data_dir/Data files/Data/6_CVs_repeat_new_baseline_new_alignment_03_CVA_C01.mpt
    and data_dir/Data files/Data/
    Ir_Ox_CV_benchmark_realigned_test_signal_Relative__0__14-52-00-835.txt

    The spectra were saved in OceanView's "Relative" mode, i.e. as % transmission
    relative to a background spectrum, so they are read as "Transmission".
    average_every=2 halves the number of spectra while reading.
    """
    data_dir = Path(data_dir) / "Data files" / "Data"
    optical = Spectrum.read(
        data_dir
        / "Ir_Ox_CV_benchmark_realigned_test_signal_Relative__0__14-52-00-835.txt",
        reader="oceanview",
        spectra_type="Transmission",
        average_every=2,
    )
    ec = Measurement.read(
        data_dir / "6_CVs_repeat_new_baseline_new_alignment_03_CVA_C01.mpt",
        reader="biologic",
    )
    return ec.as_cv(), optical


def demo_real_data_analysis():
    """Run the SEC analysis on one anodic sweep of the real data, as one figure."""
    cv, optical = read_real_measurement(REAL_DATA_DIR)

    fig = plt.figure(figsize=(30, 20))
    gs = gridspec.GridSpec(4, 4, figure=fig)

    def cell(row, col):
        return fig.add_subplot(gs[row, col])

    ax = cell(0, 0)
    cv.plot_measurement(J_name="cycle", axes=[ax, ax.twinx()])
    ax.set_title("Echem, before turning-point cycles")

    # A new cycle each time the sweep turns, so each cycle is one sweep
    cv.redefine_cycle(turning_point=True, res_points=10)
    ax = cell(1, 0)
    cv.plot_measurement(J_name="cycle", axes=[ax, ax.twinx()])
    ax.set_title("Echem, one cycle per sweep")

    sec = cv + optical

    ax_heat, ax_pot = cell(2, 0), cell(3, 0)
    sec.plot_measurement(
        wlspan=[400, 900],
        min_threshold=0,
        max_threshold=0.25,
        axes=[ax_heat, ax_pot, None, ax_pot.twinx()],
    )
    ax_heat.set_title("Full EC-optical dataset")

    print("selecting one anodic sweep")
    sweep = sec.select(cycle=5)
    # Reference the sweep to its start, at 0.66 V
    sweep.set_reference_spectrum(spectrum=sweep.spectrum_series[0])

    ax = cell(0, 1)
    sweep.plot_waterfall(ax=ax, xlim=(400, 900), ylim=(-0.02, 0.25))
    ax.set_title("Anodic sweep, dOD waterfall")

    print("difference spectra every 20 mV, smoothed")
    sweep.dOD_calculator.Vstep = 0.02
    sweep.dOD_calculator.xspan = [400, 900]
    sweep.add_calculator(Denoiser(denoise_method="Savitzky-Golay", sg_window=21))
    diff_smooth = sweep["dOD>diff>smooth"]
    ax = cell(1, 1)
    sweep.plot_dOD_difference_spectra(field=diff_smooth, xlim=(400, 900), ax=ax)
    ax.set_title("Anodic sweep, smoothed dOD difference spectra")

    print("finding converged spectra and fitting the sweep")
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sweep,
        field_name="dOD>diff>smooth",
        xspan=[400, 900],
        conv_limit=0.1,
        min_region_separation=6,
    )
    print(f"converged spectra at {decon.potentialspans}")
    sweep.add_calculator(decon)

    ax = cell(2, 1)
    ax.plot(decon.adjacent_distances)
    ax.axhline(0.1, color="k", linestyle="--")
    ax.set_xlabel("difference spectrum")
    ax.set_ylabel("distance to next")
    ax.set_title("Anodic sweep, adjacent-spectra distances")

    ax = cell(3, 1)
    sweep.plot_converged_spectra(ax=ax)
    ax.set_title("Anodic sweep, converged spectra")

    ax = cell(0, 2)
    sweep.plot_fit_reconstruction(ax=ax, xlim=(400, 900))
    ax.set_title("Anodic sweep, fit reconstruction")

    ax_top, ax_bottom = cell(1, 2), cell(2, 2)
    sweep.plot_fit_and_residuals(axes=[ax_top, ax_bottom])
    ax_top.set_title("Anodic sweep, fit coefficients and residuals")

    ax = cell(3, 2)
    sweep.plot_wavelengths_vs_cv(
        wavelengths=["w450", "w600", "w800"], axes=(ax, ax.twinx())
    )
    ax.set_title("Anodic sweep, wavelengths vs potential")

    print("tracking wavelengths across the full dataset")
    ax_top, ax_bottom = cell(0, 3), cell(1, 3)
    sec.plot_wavelengths(
        wavelengths=["w450", "w600", "w800"],
        axes=[ax_top, ax_bottom, None, ax_bottom.twinx()],
    )
    ax_top.set_title("Tracked wavelengths vs time")

    fig.suptitle("SEC analysis of one anodic sweep on IrOx", fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.97))


def main():
    if REAL_DATA_DIR:
        demo_real_data_analysis()
    else:
        demo_small_fixture()
        print(
            "\nSet IXDAT_DEMO_REAL_DATA_DIR to the unzipped IrOx_SEC_benchmark.zip to"
            " also run the full analysis on real multi-cycle data (one cycle per"
            " sweep, difference spectra, converged spectra and the fit)."
        )

    plt.show()


if __name__ == "__main__":
    main()
