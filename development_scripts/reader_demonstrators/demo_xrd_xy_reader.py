"""Sandbox script to aid in development and demo of the XRDXYReader

Downloads two real powder diffraction datasets:

  1. ZSM-5 zeolite — synchrotron .xye (3-column: 2θ, intensity, error), no header
     Source: https://github.com/stefsmeets/lines (MIT license)

  2. Rutile TiO2 — Debye-calculated I(Q) .dat (comma-separated, header "Q,I(Q)")
     ICSD entry 001504, tetragonal, from DebyeCalculator
     Source: https://github.com/FrederikLizakJohansen/DebyeCalculator (MIT license)

Run from repo root:
    python development_scripts/reader_demonstrators/demo_xrd_xy_reader.py
"""

import urllib.request
from pathlib import Path

from ixdat import Spectrum

import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).parent / "xrd_xy_data"
DATA_DIR.mkdir(exist_ok=True)

ZSM5_URL = (
    "https://raw.githubusercontent.com/stefsmeets/lines/main/testing/zsm-5.xye"
)
RUTILE_URL = (
    "https://raw.githubusercontent.com/FrederikLizakJohansen/DebyeCalculator"
    "/main/debyecalculator/unittests_files"
    "/icsd_001504_cc_r6_lc_2.85_6_tetragonal_Iq.dat"
)

zsm5_path = DATA_DIR / "zsm-5.xye"
rutile_path = DATA_DIR / "icsd_001504_rutile_Iq.dat"

for url, path in [(ZSM5_URL, zsm5_path), (RUTILE_URL, rutile_path)]:
    if not path.exists():
        print(f"Downloading {path.name} ...")
        urllib.request.urlretrieve(url, path)

# --- 1. ZSM-5 synchrotron .xye -------------------------------------------
# No header: reader defaults to "two theta / degree" + stores error column
zsm5 = Spectrum.read(zsm5_path, reader="xy", name="ZSM-5 zeolite (synchrotron)")
print(zsm5)
zsm5.plot()
plt.show()

# --- 2. Rutile TiO2 I(Q) from DebyeCalculator ----------------------------
# Header line "Q,I(Q)", comma-separated: reader detects Q-space automatically
rutile = Spectrum.read(rutile_path, reader="xy", name="Rutile TiO2 I(Q)")
print(rutile)
rutile.plot()
plt.show()
