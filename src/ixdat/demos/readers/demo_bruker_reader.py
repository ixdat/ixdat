"""Demo of the Bruker TopSpin NMR reader.

Reads the vendored 1D 1H NMR experiment under
``test_data/bruker/MTBLS1_ADG19007u_162_10`` (MetaboLights MTBLS1 human
urine, 700 MHz, ``noesypr1d``, 128 scans), prints the metadata that the
reader has lifted out of ``acqus`` / ``procs``, and plots both the processed
chemical-shift spectrum and the raw FID.

Requires the optional ``nmrglue`` dependency: ``pip install nmrglue``.
"""

import matplotlib.pyplot as plt

from ixdat import Spectrum
from ixdat.demos import get_test_data_dir
from ixdat.readers.bruker import ACQUS_KEYS, PROCS_KEYS

FOLDER = "bruker/MTBLS1_ADG19007u_162_10"


def load_spectrum(data_dir=None):
    """Return the processed 1H NMR Spectrum of a Bruker experiment folder."""
    root = get_test_data_dir(data_dir)
    return Spectrum.read(root / FOLDER, reader="bruker")


def load_fid(data_dir=None):
    """Return the raw FID (time-domain signal) of a Bruker experiment folder."""
    root = get_test_data_dir(data_dir)
    return Spectrum.read(root / FOLDER, reader="bruker", processed=False)


EXAMPLES = {"spectrum": load_spectrum, "fid": load_fid}


def print_summary(spec):
    """Print the basics and the acquisition and processing parameters of `spec`."""
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

    md = spec.metadata
    print("acquisition parameters (from acqus):")
    for k in ACQUS_KEYS:
        if k in md:
            print(f"  {k:10s} = {md[k]}")
    print()
    print("processing parameters (from procs):")
    for k in (f"proc_{k}" for k in PROCS_KEYS):
        if k in md:
            print(f"  {k:14s} = {md[k]}")
    print()
    print(
        f"full acqus dict has {len(md['acqus'])} keys; "
        f"full procs dict has {len(md['procs'])} keys."
    )


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    spec = results["spectrum"]
    fid = results["fid"]
    print_summary(spec)

    if show:
        md = spec.metadata
        ax = spec.plot(color="k", linewidth=0.6)
        ax.set_title(
            f"{spec.name}\n"
            f"{md.get('NUC1', '?')} NMR, {md.get('PULPROG', '?')}, "
            f"BF1 = {md.get('BF1', '?')} MHz, "
            f"NS = {md.get('NS', '?')}, solvent = {md.get('SOLVENT', '?')}"
        )
        ax_fid = fid.plot(color="k", linewidth=0.6)
        ax_fid.set_title(f"{fid.name} — raw FID")
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
