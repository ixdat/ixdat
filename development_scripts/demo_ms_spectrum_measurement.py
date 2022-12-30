from pathlib import Path
from ixdat import Measurement, Spectrum


data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/cinfdata/Krabbe"


spectra = Spectrum.read(
    data_dir / "spectrumseries_mass_spec.txt", reader="ixdat_spectrum"
)
spectra.heat_plot()


if True:
    meas = Measurement.read(data_dir / "masstime.csv", reader="ixdat")
    # meas.plot()

    spec_meas = meas + spectra
    spec_meas.plot()
