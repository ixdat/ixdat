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


if True:  # test Spectrum saving and loading.
    # Works!
    s_id = spec.save()
    loaded = Spectrum.get(s_id)
    ax = loaded.plot(color="g")
    ax.set_yscale("log")

if True:  # test Spectrum exporting and re-reading
    # FIXME: Doesn't work yet. MSSpectrum has no default exporter.
    spec.export("./my_spectrum.csv")
    loaded_spec = Spectrum.read("./my_spectrum.csv", reader="ixdat")
    ax = loaded_spec.plot(color="g")
    ax.set_yscale("log")


meas = Measurement.read(path_to_meas, reader="zilien", technique="MS")
meas.plot(
    mass_list=["M40", "M18"],
)

meas[1].plot()  # plots a spectrum

if False:  # test SpectroMSMeasurement saving and loading.
    # Woohoo, this works now! (Must have been cae8cf1)... but very slow...
    m_id = meas.save()
    loaded = Measurement.get(m_id)
    ax = loaded.plot(mass_list=["M2", "M18", "M32", "M40"])

if True:  # test SpectroMSMeasurement exporting and re-reading. Works! :)
    meas.export("./my_spectro_ms_measurement.csv")
    loaded_meas = Measurement.read("./my_spectro_ms_measurement.csv", reader="ixdat")
    ax = loaded_meas.plot(mass_list=["M2", "M18", "M28", "M32", "M40"])

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


meas.spectrum_series.continuous = False
meas_p1 = meas.cut(tspan=[0, 3000])
meas_p1.plot()
meas_p2 = meas.cut(tspan=[3000, 4000])

meas_joined = meas_p1 + meas_p2  # tests adding of two MSSpectroMeasurement objects
meas_joined.plot()

meas_p1_no_spec = meas_no_spec.cut(tspan=[0, 3000])
meas_p2_no_spec = meas_no_spec.cut(tspan=[3000, 4000])

meas_joined_p2_no_spec = meas_p1 + meas_p2_no_spec
meas_joined_p2_no_spec.plot()

meas_joined_p1_no_spec = meas_p1_no_spec + meas_p2
meas_joined_p1_no_spec.plot()
