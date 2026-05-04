"""Sandbox script to aid in development and demo of the XRDXYReader

Reads two real powder diffraction datasets directly from URLs:

  1. ZSM-5 zeolite -- synchrotron .xye (3-column: 2-theta, intensity, error), no header
     Source: https://github.com/stefsmeets/lines (MIT license)

  2. Rutile TiO2 — Debye-calculated I(Q) .dat (comma-separated, header "Q,I(Q)")
     ICSD entry 001504, tetragonal, from DebyeCalculator
     Source: https://github.com/FrederikLizakJohansen/DebyeCalculator (MIT license)

Run from repo root:
    python development_scripts/reader_demonstrators/demo_xrd_xy_reader.py
"""

from ixdat import Spectrum

import matplotlib.pyplot as plt

ZSM5_URL = "https://raw.githubusercontent.com/stefsmeets/lines/main/testing/zsm-5.xye"
RUTILE_URL = (
    "https://raw.githubusercontent.com/FrederikLizakJohansen/DebyeCalculator"
    "/main/debyecalculator/unittests_files"
    "/icsd_001504_cc_r6_lc_2.85_6_tetragonal_Iq.dat"
)

# --- 1. ZSM-5 synchrotron .xye -------------------------------------------
# No header: reader defaults to "two theta / degree" + stores error column
zsm5 = Spectrum.read_url(ZSM5_URL, reader="xrdxy", name="ZSM-5 zeolite (synchrotron)")
print(zsm5)
zsm5.plot()
plt.show()

# --- 2. Rutile TiO2 I(Q) from DebyeCalculator ----------------------------
# Header line "Q,I(Q)", comma-separated: reader detects Q-space automatically
rutile = Spectrum.read_url(RUTILE_URL, reader="xrdxy", name="Rutile TiO2 I(Q)")
print(rutile)
rutile.plot()
plt.show()
