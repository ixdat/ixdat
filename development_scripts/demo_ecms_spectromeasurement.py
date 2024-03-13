# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 22:02:26 2024

@author: Søren
"""

from pathlib import Path
from ixdat import Measurement


data_dir = Path.home() / (
    r"Dropbox\ixdat_resources\debugging\24B05_ECMS_with_mass_scans\23K15"
)
ms_no_spec = Measurement.read(
    data_dir / "MS/2023-11-15 12_20_46 MS.tsv",
    include_mass_scans=False,
)
ms = Measurement.read(
    data_dir / "MS/2023-11-15 12_20_46 MS.tsv",
    # include_mass_scans=True,  # not needed, they are included by default :)
)

ec = Measurement.read_set(data_dir / "EC/", suffix=".mpr")
ec.calibrate(RE_vs_RHE=0.5)

ecms_no_spec = ec + ms_no_spec
ecms_no_spec.plot()

ecms = ec + ms

ecms.plot()

ecms[1].plot()  # a mass scan!
