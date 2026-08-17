#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:40:30 2026

@author: soren
"""

from pathlib import Path
from matplotlib import pyplot as plt
from ixdat.db import change_database
from tools_for_demos import view_tables

plt.close("all")

sqlite_file = Path(".") / "all_demo_data.sqlite"
if sqlite_file.exists():
    sqlite_file.unlink()

change_database("sqlite", db_path=sqlite_file)


# demo_asimov_reader has the objects in a main() method so can't be import from

if False:
    # FIXME: ValueError: Object arrays cannot be saved when allow_pickle=False
    from demo_autolab_reader import meas as meas_autolab

    meas_autolab.save()


if True:
    from demo_avantage_reader import spectrum as spec_avantage

    spec_avantage.save()


if True:
    from demo_biologic_mpr_reader import combined_measurement_list as mm_mpr

    meas_mpr = mm_mpr[0]  # One is enough.
    meas_mpr.save()


if True:
    from demo_bruker_reader import fid as spec_bruker

    spec_bruker.save()


if True:
    from demo_cinfdata_reader import ecms_meas as meas_cinfdata

    meas_cinfdata.save()


if True:
    from demo_cinfdata_reader_TPMS import tpms as meas_tpms

    meas_tpms.save()


if False:  # requires internet
    from demo_echemdb_reader import ref_cycle as meas_echemdb

    meas_echemdb.save()


if True:
    from demo_ivium_reader import meas_cv as meas_ivium

    meas_ivium.save()


if True:
    from demo_ixdat_csv_reader import meas_loaded as meas_ixdat_csv

    meas_ixdat_csv.save()


if True:
    from demo_msrh_sec_reader import sec_meas as meas_sec

    meas_sec.save()


if True:
    from demo_msrh_sec_decay_reader import sec_meas as meas_sec_decay

    meas_sec_decay.save()


if True:
    from demo_nordic_tdms_reader import cv as meas_nordic

    meas_nordic.save()


if True:
    from demo_opus_ftir_reader import ftir as spec_opus
    from demo_opus_ftir_reader import ecftir as meas_opus_biologic

    spec_opus.save()
    meas_opus_biologic.save()  # doesn't duplicate the spectra data :)


if False:
    # FIXME: ValueError: Object arrays cannot be saved when allow_pickle=False
    from demo_pfeiffer_reader import meas as meas_pfeiffer

    meas_pfeiffer.save()


if True:
    from demo_qexafs_reader import xas_series as spec_qexafs
    from demo_qexafs_reader import ec_xas as meas_qexafs_biologic

    spec_qexafs.save()
    meas_qexafs_biologic.save()  # doesn't duplicate the spectra data :)


if False:  # requires internet
    from demo_xrd_xy_reader import zsm5 as spec_xrd_xye
    from demo_xrd_xy_reader import rutile as spec_xrd_xy

    spec_xrd_xye.save()
    spec_xrd_xy.save()


if True:
    from demo_xrdml_reader import spectrum as spec_xrdml

    spec_xrdml.save()


if True:
    from demo_zilien_reader import ecms as meas_zilien

    meas_zilien.save()


if True:
    from demo_zilien_spectrum_reader import meas_p1 as meas_zilien_spec_p1
    from demo_zilien_spectrum_reader import meas_joined as meas_zilien_spec

    meas_zilien_spec.save()
    meas_zilien_spec_p1.save()  # doesn't re-save the already-saved spectra.

view_tables(sqlite_file)
