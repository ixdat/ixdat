"""Demonstrate OceanView EC-optical data in ixdat.

With no further setup, the script reads a small OceanView optical spectrum-series
fixture and a BioLogic EC fixture, combines them into an ECOpticalMeasurement, and
shows the basic EC-optical plots. This part is self-contained and needs no external
data.

Set IXDAT_DEMO_REAL_DATA_DIR to the unzipped IrOx_SEC_benchmark.zip (a real
multi-cycle EC-optical recording on IrOx, too large to commit to the repo)
to instead run the full analysis on one anodic sweep, as one figure read row by
row: the measurement, one sweep, finding the transitions, and checking the fit.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

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
    ec_optical.plotter.plot_wavelengths(wavelengths=["w650", "w800"], tspan=[55, 95])
    title_current_figure("Wavelength tracking vs time")

    print("plotting wavelength tracking vs potential")
    ec_optical.plotter.plot_wavelengths_vs_potential(
        wavelengths=["w650", "w800"], tspan=[55, 95]
    )
    title_current_figure("Wavelength tracking vs potential")

    print("plotting wavelength tracking vs potential, with current overlaid")
    ec_optical.plotter.plot_wavelengths_vs_potential(
        wavelengths=["w650", "w800"], tspan=[55, 95], coplot=True
    )
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
    """Run the SEC analysis on the real data, as one figure read row by row:

    1. The measurement: dOD, potential, and tracked wavelengths vs time.
    2. One anodic sweep: its dOD vs potential, and its difference spectra.
    3. Finding the transitions: where the difference spectra stop changing shape,
       the converged spectra there, and how much of each the sweep contains.
    4. Checking the fit: the reconstructed spectra, the residuals, and measured vs
       fitted dOD at three wavelengths.
    """
    cv, optical = read_real_measurement(REAL_DATA_DIR)
    wavelengths = ["w450", "w600", "w800"]

    # A new cycle each time the sweep turns, so each cycle is one sweep
    cv.redefine_cycle(turning_point=True)
    sec = cv + optical

    # One subfigure per row, so each row gets its own heading. Constrained layout
    # makes room for the labels, twin axes, and colorbars.
    fig = plt.figure(figsize=(21, 24), layout="constrained")
    rows = fig.subfigures(4, 1, hspace=0.04)
    axes = np.array([row.subplots(1, 3) for row in rows])
    for row, title in zip(
        rows,
        [
            "1. The measurement",
            "2. One anodic sweep",
            "3. Finding the transitions",
            "4. Checking the fit",
        ],
    ):
        row.suptitle(title, x=0.01, ha="left", fontsize=17, fontweight="bold")

    # The same colors for each wavelength as plot_wavelengths_vs_potential uses
    cmap = plt.get_cmap("jet")
    norm = plt.Normalize(vmin=min(sec.wl), vmax=max(sec.wl))
    colors = [cmap(norm(float(wl_str[1:]))) for wl_str in wavelengths]

    print("1. the measurement")
    sec.plot_measurement(
        wlspan=[400, 900],
        zspan=[0, 0.25],
        J_name="cycle",
        axes=[axes[0, 0], axes[0, 1], None, axes[0, 1].twinx()],
    )
    axes[0, 0].set_title("1a. dOD vs time (no reference)")
    axes[0, 1].set_title("1b. Potential and cycle number vs time")
    for wl_str, color in zip(wavelengths, colors):
        sec.track_wavelength(float(wl_str[1:]))
        t, y = sec.grab(wl_str)
        axes[0, 2].plot(t, y, color=color, label=wl_str)
    axes[0, 2].set_xlabel("time / [s]")
    axes[0, 2].set_ylabel("dOD")
    axes[0, 2].legend()
    axes[0, 2].set_title("1c. Tracked wavelengths vs time")

    print("2. one anodic sweep, referenced to its start at 0.66 V")
    sweep = sec.select(cycle=6)
    sweep.set_reference_spectrum(spectrum=sweep.spectrum_series[0])
    sweep.plot_waterfall(ax=axes[1, 0])
    axes[1, 0].set_xlim(400, 900)
    axes[1, 0].set_ylim(-0.02, 0.25)
    axes[1, 0].set_title("2a. dOD spectra, colored by potential")
    sweep.plotter.plot_wavelengths_vs_potential(
        wavelengths=wavelengths, axes=(axes[1, 1], axes[1, 1].twinx())
    )
    axes[1, 1].set_title("2b. Tracked wavelengths and current vs potential")

    # Difference spectra every 20 mV, smoothed over 21 wavelengths
    sweep.dOD_calculator.Vstep = 0.02
    sweep.dOD_calculator.xspan = [400, 900]
    sweep.add_calculator(Denoiser(denoise_method="Savitzky-Golay", sg_window=21))
    diff_smooth = sweep["dOD>diff>smooth"]
    sweep.plotter.plot_dOD_difference_spectra(field=diff_smooth, ax=axes[1, 2])
    axes[1, 2].set_xlim(400, 900)
    axes[1, 2].set_ylim(-1.1, 1.1)
    axes[1, 2].set_title("2c. Smoothed difference spectra (normalised)")

    print("3. finding the transitions")
    conv_limit = 0.1
    decon = SpectrumDeconvolutor.from_spectromeasurement(
        sweep,
        field_name="dOD>diff>smooth",
        xspan=[400, 900],
        conv_limit=conv_limit,
        min_region_separation=6,
    )
    print(f"converged spectra at {decon.potentialspans}")
    sweep.add_calculator(decon)

    tseries = diff_smooth.axes_series[0]
    U_diff = sweep.grab_for_t(
        "potential", t=tseries.t + tseries.tstamp - sweep.tstamp
    )
    axes[2, 0].plot(U_diff[:-1], decon.adjacent_distances)
    axes[2, 0].axhline(conv_limit, color="k", linestyle="--", label="conv_limit")
    for span in decon.potentialspans:
        axes[2, 0].axvspan(*span, color="gainsboro")
    axes[2, 0].set_xlabel(sweep.U_name)
    axes[2, 0].set_ylabel("distance to next difference spectrum")
    axes[2, 0].legend()
    axes[2, 0].set_title("3a. Convergent regions (grey)")
    sweep.plotter.plot_converged_spectra(ax=axes[2, 1])
    axes[2, 1].set_title("3b. Converged spectra (normalised)")
    sweep.plotter.plot_fit_and_residuals(axes=[axes[2, 2], axes[3, 1]])
    axes[3, 1].set_xlim(400, 900)
    axes[2, 2].set_title("3c. Fit coefficient of each converged spectrum")

    print("4. checking the fit")
    sweep.plotter.plot_fit_reconstruction(ax=axes[3, 0])
    axes[3, 0].set_xlim(400, 900)
    axes[3, 0].set_ylim(-0.02, 0.25)
    axes[3, 0].set_title("4a. Reconstructed dOD spectra (compare 2a)")
    axes[3, 1].set_title("4b. Residuals / [mOD]")
    measured = sweep["dOD"]
    fitted = sweep["reconstructed"]
    U = sweep.grab_for_t(
        "potential",
        t=measured.axes_series[0].t + measured.axes_series[0].tstamp - sweep.tstamp,
    )
    for wl_str, color in zip(wavelengths, colors):
        wl = float(wl_str[1:])
        i_meas = np.argmin(np.abs(measured.axes_series[1].data - wl))
        i_fit = np.argmin(np.abs(fitted.axes_series[1].data - wl))
        axes[3, 2].plot(U, measured.data[:, i_meas], color=color, label=wl_str)
        axes[3, 2].plot(U, fitted.data[:, i_fit], color="k", linestyle="--")
    axes[3, 2].plot([], [], color="k", linestyle="--", label="fit")
    axes[3, 2].set_xlabel(sweep.U_name)
    axes[3, 2].set_ylabel("dOD")
    axes[3, 2].legend()
    axes[3, 2].set_title("4c. Measured and fitted dOD vs potential")

    fig.suptitle("SEC analysis of IrOx: one anodic sweep", fontsize=20)


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
