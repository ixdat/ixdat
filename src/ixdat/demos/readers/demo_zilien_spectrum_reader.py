"""Demo of the zilien reader with mass spectra. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Spectrum, Measurement
from ixdat.demos import get_demo_data_dir

FOLDER = "zilien_with_spectra"
MEAS_NAME = "2023-05-16 11_34_16 mix_cal_gas_glass_slide.tsv"
SPEC_NAME = (
    "mix_cal_gas_glass_slide mass scans/"
    "mass scan started at measurement time 0000066.tsv"
)


def load_spectrum(data_dir=None):
    """Return a mass Spectrum read from a single zilien mass scan file."""
    root = get_demo_data_dir(data_dir)
    return Spectrum.read(root / FOLDER / SPEC_NAME, reader="zilien")


def load_ms_with_spectra(data_dir=None):
    """Return an MS measurement which includes the mass scans as spectra."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / FOLDER / MEAS_NAME,
        reader="zilien",
        technique="MS-MS_spectra",  # include MS spectra
        # technique="MS",  # do not include MS spectra by default
        # include_mass_scans=True,  # include spectra! (overrides technique)
        # include_mass_scans=False,  # don't include spectra! (overrides technique)
    )


def load_ms_no_spectra(data_dir=None):
    """Return an MS measurement which leaves out the mass scans."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / FOLDER / MEAS_NAME,
        reader="zilien",
        technique="MS",
        include_mass_scans=False,
    )


EXAMPLES = {
    "spectrum": load_spectrum,
    "ms_with_spectra": load_ms_with_spectra,
    "ms_no_spectra": load_ms_no_spectra,
}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    spec = results["spectrum"]
    meas = results["ms_with_spectra"]
    meas_no_spec = results["ms_no_spectra"]

    # Test Spectrum exporting and re-reading
    spec.export("./my_spectrum.csv")
    results["reloaded_spectrum"] = Spectrum.read("./my_spectrum.csv", reader="ixdat")

    if show:
        spec.plot(color="k")
        ax = results["reloaded_spectrum"].plot(color="g")
        ax.set_yscale("log")

        meas.plot(mass_list=["M40", "M18"])
        meas[1].plot()  # plots a spectrum

    if False:  # test SpectroMSMeasurement exporting and re-reading. Works! :)
        meas.export("./my_spectro_ms_measurement.csv")
        loaded_meas = Measurement.read("./my_spectro_ms_measurement.csv", reader="ixdat")
        ax = loaded_meas.plot(mass_list=["M2", "M18", "M28", "M32", "M40"])

    meas.spectrum_series.continuous = True
    if show:
        meas.plot(mass_list=["M40", "M18"])
        meas_no_spec.plot(mass_list=["M40", "M18"])

    # Cutting and adding MS measurements with and without spectra
    meas.spectrum_series.continuous = False
    meas_p1 = meas.cut(tspan=[0, 3000])
    meas_p2 = meas.cut(tspan=[3000, 4000])
    meas_p3 = meas.cut(tspan=[4500, 5000])  # no spectra here!
    meas_p1_no_spec = meas_no_spec.cut(tspan=[0, 3000])
    meas_p2_no_spec = meas_no_spec.cut(tspan=[3000, 4000])

    results["p1"] = meas_p1
    results["joined"] = meas_p1 + meas_p2
    results["joined_p2_no_spec"] = meas_p1 + meas_p2_no_spec
    results["joined_p1_no_spec"] = meas_p1_no_spec + meas_p2
    print(f"spectra in 4500-5000 s: {len(meas_p3.spectrum_series)}")  # 0
    results["joined_all"] = meas_p1 + meas_p2 + meas_p3  # order doesn't matter!
    n_joined = len(results["joined_all"].spectrum_series)
    print(f"spectra in the joined measurement: {n_joined}")  # 4

    if show:
        meas_p1.plot()
        results["joined"].plot()
        results["joined_p2_no_spec"].plot()
        results["joined_p1_no_spec"].plot()
        meas_p3.plot()  # bottom panel is empty
        results["joined_all"].plot()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
