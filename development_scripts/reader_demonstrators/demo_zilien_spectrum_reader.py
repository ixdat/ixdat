"""For use in development of zilien spectrum reader. Requires access to sample data."""

from pathlib import Path
from ixdat import Spectrum, Measurement


data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/zilien_with_spectra"

path_to_meas = data_dir / "2023-05-16 11_34_16 mix_cal_gas_glass_slide.tsv"
path_to_spec = (
    data_dir
    / "mix_cal_gas_glass_slide mass scans/mass scan started at measurement time 0000066.tsv"
)


spec = Spectrum.read(
    path_to_spec,
    reader="zilien",
)

spec.plot(color="k")


if False:  # test Spectrum saving and loading.
    # Works!
    s_id = spec.save()
    loaded = Spectrum.get(s_id)
    ax = loaded.plot(color="g")
    ax.set_yscale("log")

if False:  # test Spectrum exporting and re-reading
    # FIXME: Doesn't work yet. MSSpectrum has no default exporter.
    s_id = spec.export("./my_spectrum.csv")
    loaded = Spectrum.read("./my_spectrum.csv", reader="ixdat")
    ax = loaded.plot(color="g")
    ax.set_yscale("log")


meas = Measurement.read(path_to_meas, reader="zilien", technique="MS")
meas.plot(
    mass_list=["M40", "M18"],
)

meas[1].plot()  # plots a spectrum

if False:  # test SpectroMSMeasurement saving and loading.
    # FIXME: broken :(  Object-relational mapping should ensure this always works.
    m_id = meas.save()
    loaded = Measurement.get(m_id)
    ax = loaded.plot(mass_list=["M2", "M18", "M32", "M40"])

if False:  # test SpectroMSMeasurement exporting and re-reading
    # FIXME: Doesn't work yet. SpectroMSMeasurement.export doesn't save spectra.
    meas.export("./my_spectro_ms_measurement.csv")
    loaded = Measurement.read("./my_spectro_ms_measurement.csv", reader="ixdat")
    ax = loaded.plot(mass_list=["M2", "M18", "M28", "M32", "M40"])

meas.spectrum_series.continuous = True
meas.plot(
    mass_list=["M40", "M18"],
)

meas_no_spec = Measurement.read(
    path_to_meas, reader="zilien", technique="MS", include_mass_scans=False
)
meas_no_spec.plot(
    mass_list=["M40", "M18"],
)
