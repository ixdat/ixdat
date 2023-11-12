from pathlib import Path
from ixdat import Measurement, Spectrum


data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/cinfdata/Krabbe"


spectra = Spectrum.read(data_dir / "spectrumseries_mass_spec.txt", reader="ixdat")
spectra.heat_plot()


if True:
    meas = Measurement.read(
        data_dir / "masstime.csv",
        reader="ixdat",
        technique="reactor",
        aliases={"pressure": ["Reactor pressure"], "temperature": ["TC temperature"]},
    )
    # meas.plot()

    spec_meas = meas + spectra
    axes = spec_meas.plot()

    axes[0].get_figure().savefig("spectro_tpms_plot.png")
