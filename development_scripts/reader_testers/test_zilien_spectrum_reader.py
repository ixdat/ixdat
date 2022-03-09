"""For use in development of zilien spectrum reader. Requires access to sample data."""

from pathlib import Path
from ixdat import Spectrum

path_to_file = (
    Path(r"C:\Users\scott\Dropbox\ixdat_resources\test_data\zilien_spectra")
    / "mass scan started at measurement time 0001700.tsv"
)

spec = Spectrum.read(
    path_to_file,
    reader="zilien_spec",
)

spec.plot(color="k")

s_id = spec.save()

loaded = Spectrum.get(s_id)
loaded.plot(color="g")
