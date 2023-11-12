"""For use in development of zilien spectrum reader. Requires access to sample data."""

from pathlib import Path
from ixdat import Spectrum, Measurement

data_dir = (
    Path.home() / "/home/scott/Dropbox/ixdat_resources/test_data/zilien_with_spectra"
)

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


meas = Measurement.read(path_to_meas, reader="zilien", technique="MS")

meas.plot()


if True:  # test Spectrum saving and loading
    s_id = spec.save()
    loaded = Spectrum.get(s_id)
    ax = loaded.plot(color="g")
    ax.set_yscale("log")
