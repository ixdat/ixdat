"""Demonstrate OceanView EC-optical data in ixdat.

With no further setup, the script reads a small OceanView optical spectrum-series
fixture and a BioLogic EC fixture, combines them into an ECOpticalMeasurement, and
shows the basic EC-optical plots. This part is self-contained and needs no external
data.

Set IXDAT_DEMO_REAL_DATA_DIR to a directory with echem/ and uv-vis/ subfolders (a
real multi-cycle Ni(OH)2/KOH EC-optical recording, too large to commit to the repo)
to instead run the full cycle-based analysis this feature was built for, all as
subplots of one figure for cycle 4: turning-point cycles, anodic/cathodic
waterfalls, wavelength tracking, and the convergent-spectra/fitting pipeline, using
the same ROIs and parameters as the analysis this was developed against.
"""

import os
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from ixdat import Measurement, Spectrum


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
    """Read the real multi-cycle CV and its optical data, and calibrate the CV.

    Expects data_dir/echem/01_CVs_04_CVA_C01.mpt and
    data_dir/uv-vis/01_CVs_25msinteg_av4_QEP065181__0__09-54-15-000.txt.

    boxcar_width and average_every reduce the spectra to about 1s (1 mV, for this
    scan's 1 mV/s rate) resolution -- the file itself is at 100 ms resolution
    (25 ms integration, averaged 4x by OceanView), too fine to be worth keeping.
    """
    data_dir = Path(data_dir)
    optical = Spectrum.read(
        data_dir / "uv-vis" / "01_CVs_25msinteg_av4_QEP065181__0__09-54-15-000.txt",
        reader="oceanview",
        spectra_type="Intensity",
        boxcar_width=50,
        average_every=10,
    )
    ec = Measurement.read(
        data_dir / "echem" / "01_CVs_04_CVA_C01.mpt", reader="biologic"
    )
    return ec.as_cv(), optical


def demo_real_data_analysis():
    """Reproduce the cycle-based SEC analysis this feature was built for, as one
    figure of subplots -- one representative cycle (4) rather than a cycle-to-cycle
    comparison, so this stays a single cycle's worth of panels.
    """
    cv, optical = read_real_measurement(REAL_DATA_DIR)

    fig = plt.figure(figsize=(30, 20))
    gs = gridspec.GridSpec(4, 5, figure=fig)

    def cell(row, col):
        return fig.add_subplot(gs[row, col])

    cv.calibrate(RE_vs_RHE=0.8663)
    ax = cell(0, 0)
    cv.plot_measurement(J_name="cycle", axes=[ax, ax.twinx()])
    ax.set_title("Echem, before turning-point cycles")

    cv.redefine_cycle(turning_point=True, redox=True, N_points=5, N_sep=900)
    ax = cell(1, 0)
    cv.plot_measurement(J_name="cycle", axes=[ax, ax.twinx()])
    ax.set_title("Echem, turning-point cycles")

    cv.calibrate(R_Ohm=16.2)
    ax = cell(2, 0)
    cv.plot_measurement(axes=[ax, ax.twinx()])
    ax.set_title("Echem, iR-compensated")

    ax = cell(3, 0)
    cv[4].plot(ax=ax)
    ax.set_title("Cycle 4 alone")

    ec_optical = cv + optical
    ec_optical.set_reference_spectrum(t_ref=ec_optical.t[0])

    wl = ec_optical.wavelength.data
    wl_mask = (wl >= 300) & (wl <= 900)

    print("finding the start time of cycle 4, anodic and cathodic")
    anodic_dod = ec_optical.get_dOD_cycle(cycle_number=4, direction=0, N_points=200)
    anodic_start = anodic_dod.axes_series[0].data.min()
    cathodic_dod = ec_optical.get_dOD_cycle(cycle_number=4, direction=1, N_points=200)
    cathodic_start = cathodic_dod.axes_series[0].data.min()

    print("plotting the anodic and cathodic waterfalls for cycle 4")
    for row, direction, t_ref, label in (
        (0, 0, anodic_start, "anodic"),
        (1, 1, cathodic_start, "cathodic"),
    ):
        cycle_data = ec_optical.get_dOD_cycle(
            t_ref=t_ref, cycle_number=4, direction=direction, N_points=200
        ).data[:, wl_mask]
        ax = cell(row, 1)
        ec_optical.plot_waterfall_cycle(
            t_ref=t_ref,
            cycle_number=4,
            xlim=(300, 900),
            ylim=(cycle_data.min() - 0.01, cycle_data.max() + 0.01),
            direction=direction,
            N_points=200,
            ax=ax,
        )
        ax.set_title(f"Cycle 4, {label} waterfall")

    print("cycle 4: dOD difference spectra, convergent spectra, and fit")
    diff_4 = ec_optical.get_dOD_difference_spectra(
        cycle_number=4,
        N_points=100,
        wlmin=300,
        wlmax=900,
        step=5,
        direction=0,
        normalise=True,
        t_ref=anodic_start,
    )
    ax = cell(2, 1)
    ec_optical.plot_dOD_difference_spectra(
        diff_spectra=diff_4, xlim=(350, 900), ylim=(0, 1), ax=ax
    )
    ax.set_title("Cycle 4, dOD difference spectra")

    converge_4 = ec_optical.get_convergent_spectra(
        diff_spectra_field=diff_4,
        wlmin=300,
        wlmax=900,
        conv_limit=0.2,
        converged_spectrum_form="average",
        min_region_width=1,
        min_region_separation=10,
        smooth_distances=True,
        window_length=20,
    )
    ax = cell(3, 1)
    ax.plot(converge_4[2])
    ax.set_title("Cycle 4, adjacent-spectra distances")

    ax = cell(0, 2)
    ec_optical.plot_convergent_spectra(converge_output=converge_4, ax=ax)
    ax.set_title("Cycle 4, convergent spectra")

    cycle_4 = ec_optical.get_dOD_cycle(
        t_ref=anodic_start, cycle_number=4, direction=0, N_points=100
    )
    ax = cell(1, 2)
    ec_optical.plot_fit_reconstruction(
        cycle_data=cycle_4,
        converged_data=converge_4,
        noise_bound=10,
        bounds_bool=True,
        wlmin=300,
        ax=ax,
    )
    ax.set_title("Cycle 4, fit reconstruction")

    ax = cell(2, 2)
    ec_optical.plot_wavelengths_vs_cv(
        wavelengths=["w404", "w500", "w345"],
        tspan=(anodic_start, cathodic_start),
        axes=(ax, ax.twinx()),
    )
    ax.set_title("Cycle 4, wavelengths vs potential")

    ax_heat, ax_pot = cell(0, 3), cell(1, 3)
    ec_optical.plot(
        wlspan=[300, 900],
        min_threshold=-1.5,
        max_threshold=1.5,
        axes=[ax_heat, ax_pot, None, ax_pot.twinx()],
    )
    ax_heat.set_title("Full EC-optical dataset")

    ax_top, ax_bottom = cell(2, 3), cell(3, 3)
    ec_optical.plot_fit_and_residuals(
        cycle_data=cycle_4,
        converged_data=converge_4,
        noise_bound=10,
        bounds_bool=True,
        wlmin=300,
        axes=[ax_top, ax_bottom],
    )
    ax_top.set_title("Cycle 4, fit and residuals")

    print("tracking wavelengths across the full dataset")
    ax_top, ax_bottom = cell(0, 4), cell(1, 4)
    ec_optical.plot_wavelengths(
        wavelengths=["w404", "w500", "w345"],
        axes=[ax_top, ax_bottom, None, ax_bottom.twinx()],
    )
    ax_top.set_title("Tracked wavelengths vs time")

    fig.suptitle("Cycle 4 SEC analysis", fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.97))


def main():
    if REAL_DATA_DIR:
        demo_real_data_analysis()
    else:
        demo_small_fixture()
        print(
            "\nSet IXDAT_DEMO_REAL_DATA_DIR to a directory with echem/ and uv-vis/ "
            "subfolders to also run the full cycle-based analysis on real "
            "multi-cycle data (turning-point cycles, waterfalls, cycle-to-cycle "
            "comparison, and the convergent-spectra/fitting pipeline)."
        )

    plt.show()


if __name__ == "__main__":
    main()
