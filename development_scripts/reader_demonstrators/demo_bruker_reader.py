"""Demo for the Bruker TopSpin NMR reader.

Reads the vendored 1D 1H NMR experiment under
``test_data/bruker/MTBLS1_ADG19007u_162_10`` (MetaboLights MTBLS1 human
urine, 700 MHz, ``noesypr1d``, 128 scans), prints the metadata that the
reader has lifted out of ``acqus`` / ``procs``, and plots the chemical
shift spectrum.

Requires the optional ``nmrglue`` dependency: ``pip install nmrglue``.
"""

from pathlib import Path

from ixdat import Spectrum
from ixdat.readers.bruker import ACQUS_KEYS, PROCS_KEYS


DATA_DIR = (
    Path(__file__).parent.parent.parent
    / "test_data"
    / "bruker"
    / "MTBLS1_ADG19007u_162_10"
)


# 1. Read the experiment folder.
spec = Spectrum.read(DATA_DIR, reader="bruker")

# 2. Print the basics.
print(f"class       : {type(spec).__name__}")
print(f"technique   : {spec.technique}")
print(f"name        : {spec.name}")
print(f"tstamp      : {spec.tstamp}")
print(
    f"x ({spec.xseries.unit_name}) : "
    f"{spec.x.min():.2f} .. {spec.x.max():.2f}  ({spec.x.size} points)"
)
print(f"y           : {spec.y.min():.3g} .. {spec.y.max():.3g}")
print()

# 3. Print the lifted acquisition + processing parameters.
md = spec.metadata
acq_keys = ACQUS_KEYS
proc_keys = tuple(f"proc_{k}" for k in PROCS_KEYS)

print("acquisition parameters (from acqus):")
for k in acq_keys:
    if k in md:
        print(f"  {k:10s} = {md[k]}")
print()
print("processing parameters (from procs):")
for k in proc_keys:
    if k in md:
        print(f"  {k:14s} = {md[k]}")
print()
print(
    f"full acqus dict has {len(md['acqus'])} keys; "
    f"full procs dict has {len(md['procs'])} keys."
)

# 4. Plot the spectrum.
ax = spec.plot(color="k", linewidth=0.6)
ax.set_title(
    f"{spec.name}\n"
    f"{md.get('NUC1', '?')} NMR, {md.get('PULPROG', '?')}, "
    f"BF1 = {md.get('BF1', '?')} MHz, "
    f"NS = {md.get('NS', '?')}, solvent = {md.get('SOLVENT', '?')}"
)
ax.get_figure().tight_layout()
ax.get_figure().show()
