"""Demonstrate OceanView EC-optical data in ixdat.

The script reads a small OceanView optical spectrum-series fixture and a BioLogic EC
fixture, combines them into an ECOpticalMeasurement, and shows the basic EC-optical
plots, the cycle-based dOD analysis, and the denoising/convergent-spectra pipeline.
"""

from pathlib import Path

import matplotlib.pyplot as plt

from ixdat import Measurement, Spectrum


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DATA_DIR = REPO_ROOT / "test_data" / "oceanview_sec"

OPTICAL_FILE = DATA_DIR / "mini_oceanview__0__15-02-35-123.txt"
EC_FILE = DATA_DIR / "demo_oceanview_ecoptical_sec.mpt"


def read_ec_optical_measurement():
    """Read the EC and optical files and combine them into ECOpticalMeasurement.

    The EC data is read as a CyclicVoltammogram first, so redefine_cycle can define
    cycles from turning points in the potential sweep rather than a fixed potential.
    `cv + optical` then combines it with the optical data directly, into the same
    ECOpticalMeasurement that `ec + optical` would give for a plain ECMeasurement.
    """
    optical = Spectrum.read(OPTICAL_FILE, reader="oceanview")
    ec = Measurement.read(EC_FILE, reader="biologic")
    cv = ec.as_cv()
    cv.redefine_cycle(turning_point=True, redox=None, N_points=5, N_sep=10)
    ec_optical = cv + optical
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


def main():
    demo_average_every()

    print("reading EC and optical data")
    ec_optical = read_ec_optical_measurement()
    print(f"combined measurement technique: {ec_optical.technique}")
    print(
        "optical spectra: "
        f"{len(ec_optical.spectrum_series)} spectra x {len(ec_optical.wavelength.data)} wavelengths"
    )
    # This fixture's optical window (55-95 s) only spans one CV cycle, so the
    # direction-based (anodic/cathodic) split needs more data than it has -- a
    # turning point has to fall inside the optical window. Real datasets covering
    # several cycles can pass direction=0 (anodic) or direction=1 (cathodic).

    print("plotting EC-optical heat map")
    ec_optical.plot_measurement(tspan=[55, 95], wlspan=[400, 900], t_ref=55)

    print("plotting waterfall")
    ec_optical.plot_waterfall(t_ref=55)

    print("plotting wavelength tracking")
    ec_optical.plot_wavelengths(wavelengths=["w650", "w800"], tspan=[55, 95])

    print("plotting wavelength tracking vs potential")
    ec_optical.plot_wavelengths_vs_potential(
        wavelengths=["w650", "w800"], tspan=[55, 95]
    )

    print("plotting wavelength tracking vs potential, with current overlaid")
    ec_optical.plot_wavelengths_vs_cv(wavelengths=["w650", "w800"], tspan=[55, 95])

    print("getting and plotting the dOD difference spectra for cycle 0")
    diff_spectra = ec_optical.get_dOD_difference_spectra(cycle_number=0)
    ec_optical.plot_dOD_difference_spectra(diff_spectra=diff_spectra)

    print("denoising the difference spectra")
    denoised = ec_optical.denoise_spectra(
        spectra_field=diff_spectra,
        denoise_method="Savitzky-Golay",
        sg_window=5,
        sg_poly_order=2,
    )
    ec_optical.plotter.plot_waterfall_vs(
        measurement=ec_optical,
        field=denoised,
        vs=ec_optical.U_name,
    )

    print("getting and plotting the convergent spectra for cycle 0")
    converge = ec_optical.get_convergent_spectra(
        cycle_number=0, conv_limit=0.9, min_region_width=1
    )
    ec_optical.plot_convergent_spectra(converge_output=converge)

    plt.show()


if __name__ == "__main__":
    main()
