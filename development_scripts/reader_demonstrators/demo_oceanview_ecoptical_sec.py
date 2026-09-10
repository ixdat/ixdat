"""Demonstrate OceanView EC-optical data in ixdat.

With no further setup, the script reads a small OceanView optical spectrum-series
fixture and a BioLogic EC fixture, combines them into an ECOpticalMeasurement, and
shows the basic EC-optical plots. This part is self-contained and needs no external
data.

Set IXDAT_DEMO_REAL_DATA_DIR to a directory with echem/ and uv-vis/ subfolders (a
real multi-cycle Ni(OH)2/KOH EC-optical recording, too large to commit to the repo)
to instead run the full cycle-based analysis this feature was built for: turning-point
cycles, anodic/cathodic waterfalls, cycle-to-cycle comparison, wavelength tracking,
and the convergent-spectra/fitting pipeline, using the same ROIs and parameters as
the analysis this was developed against.
"""

import os
from pathlib import Path

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
    cv = ec.as_cv()
    cv.calibrate(RE_vs_RHE=0.8663)

    print("plotting the full echem before cycles are redefined")
    cv.plot_measurement(J_name="cycle")
    title_current_figure("Full echem, before turning-point cycles")

    cv.redefine_cycle(turning_point=True, redox=True, N_points=5, N_sep=900)

    print("plotting the full echem with turning-point cycles")
    cv.plot_measurement(J_name="cycle")
    title_current_figure("Full echem, turning-point cycles")

    cv.calibrate(R_Ohm=16.2)
    print("plotting the full, iR-compensated echem")
    cv.plot_measurement()
    title_current_figure("Full echem, iR-compensated")

    print("plotting cycle 4 alone")
    cv[4].plot()
    title_current_figure("Cycle 4")

    return cv, optical


def demo_real_data_analysis():
    """Reproduce the cycle-based SEC analysis this feature was built for.

    Cycles 4 and 5 are the least-changing pair in this dataset, so they're the
    ones examined in detail: waterfalls, a cycle-to-cycle comparison, and the
    convergent-spectra/fitting pipeline (with cycle 5 re-run over a narrower
    wavelength range, since below 425 nm is where the two cycles diverge).
    """
    cv, optical = read_real_measurement(REAL_DATA_DIR)

    ec_optical = cv + optical
    ec_optical.set_reference_spectrum(t_ref=ec_optical.t[0])

    print("plotting the full EC-optical dataset (350-900 nm is where signal lives)")
    ec_optical.plot(wlspan=[300, 900], min_threshold=-1.5, max_threshold=1.5)
    title_current_figure("Full EC-optical dataset")

    wl = ec_optical.wavelength.data
    wl_mask = (wl >= 300) & (wl <= 900)

    print("finding the start time of cycles 4 and 5, anodic and cathodic")
    anodic_start = {}
    cathodic_start = {}
    for cycle in (4, 5):
        dod = ec_optical.get_dOD_cycle(cycle_number=cycle, direction=0, N_points=200)
        anodic_start[cycle] = dod.axes_series[0].data.min()
        dod = ec_optical.get_dOD_cycle(cycle_number=cycle, direction=1, N_points=200)
        cathodic_start[cycle] = dod.axes_series[0].data.min()

    print("plotting the anodic and cathodic waterfalls for cycles 4 and 5")
    for cycle in (4, 5):
        for direction, start_times, label in (
            (0, anodic_start, "anodic"),
            (1, cathodic_start, "cathodic"),
        ):
            cycle_data = ec_optical.get_dOD_cycle(
                t_ref=start_times[cycle],
                cycle_number=cycle,
                direction=direction,
                N_points=200,
            ).data[:, wl_mask]
            ec_optical.plot_waterfall_cycle(
                t_ref=start_times[cycle],
                cycle_number=cycle,
                xlim=(300, 900),
                ylim=(cycle_data.min() - 0.01, cycle_data.max() + 0.01),
                direction=direction,
                N_points=200,
            )
            title_current_figure(f"Cycle {cycle}, {label}")

    print("plotting the difference between cycle 4 and cycle 5")
    ec_optical.plot_dOD_cycle_diff(
        direction=0,
        N_points=200,
        cycle_1=4,
        cycle_2=5,
        t_ref1=anodic_start[4],
        t_ref2=anodic_start[5],
        upper_noise_bound=1,
        lower_noise_bound=-1,
        xlim=(300, 900),
    )
    title_current_figure("Cycle 4 vs cycle 5, anodic difference")
    ec_optical.plot_dOD_cycle_diff(
        direction=1,
        N_points=200,
        cycle_1=4,
        cycle_2=5,
        t_ref1=cathodic_start[4],
        t_ref2=cathodic_start[5],
        upper_noise_bound=1,
        lower_noise_bound=-1,
        xlim=(300, 900),
    )
    title_current_figure("Cycle 4 vs cycle 5, cathodic difference")

    print("tracking wavelengths across the full dataset and across cycles 4-5")
    ec_optical.plot_wavelengths(wavelengths=["w404", "w500", "w345"])
    title_current_figure("Tracked wavelengths vs time")
    ec_optical.plot_wavelengths_vs_cv(
        wavelengths=["w404", "w500", "w345"],
        tspan=(anodic_start[4], anodic_start[5]),
    )
    title_current_figure("Tracked wavelengths vs potential (cycles 4-5)")

    print("cycle 4: dOD difference spectra, convergent spectra, and fit")
    diff_4 = ec_optical.get_dOD_difference_spectra(
        cycle_number=4,
        N_points=100,
        wlmin=300,
        wlmax=900,
        step=5,
        direction=0,
        normalise=True,
        t_ref=anodic_start[4],
    )
    ec_optical.plot_dOD_difference_spectra(
        diff_spectra=diff_4, xlim=(350, 900), ylim=(0, 1)
    )
    title_current_figure("Cycle 4, dOD difference spectra")

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
    plt.figure()
    plt.plot(converge_4[2])
    title_current_figure("Cycle 4, adjacent-spectra distances")

    ec_optical.plot_convergent_spectra(converge_output=converge_4)
    title_current_figure("Cycle 4, convergent spectra")

    cycle_4 = ec_optical.get_dOD_cycle(
        t_ref=anodic_start[4], cycle_number=4, direction=0, N_points=100
    )
    ec_optical.plot_fit_and_residuals(
        cycle_data=cycle_4,
        converged_data=converge_4,
        noise_bound=10,
        bounds_bool=True,
        wlmin=300,
    )
    title_current_figure("Cycle 4, fit and residuals")
    ec_optical.plot_fit_reconstruction(
        cycle_data=cycle_4,
        converged_data=converge_4,
        noise_bound=10,
        bounds_bool=True,
        wlmin=300,
    )
    title_current_figure("Cycle 4, fit reconstruction")

    print("cycle 5: same pipeline, narrower wavelength range (425-900 nm)")
    diff_5 = ec_optical.get_dOD_difference_spectra(
        cycle_number=5,
        N_points=200,
        wlmin=425,
        wlmax=900,
        step=10,
        direction=0,
        normalise=True,
        t_ref=anodic_start[4],
    )
    ec_optical.plot_dOD_difference_spectra(
        diff_spectra=diff_5, xlim=(425, 900), ylim=(0, 1)
    )
    title_current_figure("Cycle 5, dOD difference spectra")

    converge_5 = ec_optical.get_convergent_spectra(
        diff_spectra_field=diff_5,
        wlmin=425,
        wlmax=900,
        conv_limit=0.05,
        converged_spectrum_form="average",
        min_region_width=1,
        min_region_separation=5,
        smooth_distances=True,
        window_length=10,
    )
    plt.figure()
    plt.plot(converge_5[2])
    title_current_figure("Cycle 5, adjacent-spectra distances")

    ec_optical.plot_convergent_spectra(converge_output=converge_5)
    title_current_figure("Cycle 5, convergent spectra")

    cycle_5 = ec_optical.get_dOD_cycle(
        t_ref=anodic_start[5], cycle_number=5, direction=0, N_points=200
    )
    ec_optical.plot_fit_and_residuals(
        cycle_data=cycle_5,
        converged_data=converge_5,
        noise_bound=0,
        bounds_bool=True,
        wlmin=425,
    )
    title_current_figure("Cycle 5, fit and residuals")
    ec_optical.plot_fit_reconstruction(
        cycle_data=cycle_5,
        converged_data=converge_5,
        noise_bound=0,
        bounds_bool=True,
        wlmin=425,
    )
    title_current_figure("Cycle 5, fit reconstruction")


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
