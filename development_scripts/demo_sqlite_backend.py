# -*- coding: utf-8 -*-
"""A tour of ixdat's SQLite database backend.

This demo saves real measurements to a relational database in a single local
file and shows off everything you gain by doing so:

1.  saving is one line, and the whole schema is created for you
2.  everything lives in one portable file with referential integrity
3.  loading is lazy - metadata first, numerical data only on demand
4.  objects come back by name or id, as the class they were saved as,
    with their calculators (e.g. calibrations) intact
5.  measurements built from other measurements share rows - no data duplication
6.  every ixdat data family works: measurements, spectra, calculators, ...
7.  the database speaks SQL: query it with pandas or any SQLite tool
8.  rows can be updated in place
9.  the full relational schema of ixdat's data model is generated from the
    classes themselves (what PR #75 drew by hand: https://dbdiagram.io/d/601b...)

It only uses data files shipped in the ixdat repository, so it runs anywhere:

    python demo_sqlite_backend.py
"""

import shutil
import sqlite3
import tempfile
import time
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from ixdat import Measurement, Spectrum
from ixdat.db import DB, change_database

plt.close("all")

TEST_DATA_DIR = Path(__file__).absolute().parent.parent / "test_data"

DEMO_DIR = Path(tempfile.gettempdir()) / "ixdat_sqlite_demo"
shutil.rmtree(DEMO_DIR, ignore_errors=True)  # start fresh on every run
DB_FILE = DEMO_DIR / "demo_project.sqlite"


def part(title):
    print("\n" + "=" * 75 + f"\n{title}\n" + "=" * 75)


# --------------------------------------------------------------------------- #
part("1. One line to point ixdat at a database file")
# --------------------------------------------------------------------------- #

change_database("sqlite", db_path=DB_FILE)
print(f"Active backend: {DB.backend}")
print("(No new dependencies - this is python's built-in sqlite3.)")

# --------------------------------------------------------------------------- #
part("2. Saving a measurement (and its calibration) is one line each")
# --------------------------------------------------------------------------- #

ec = Measurement.read(TEST_DATA_DIR / "biologic/Pt_poly_cv_CUT.mpt", reader="biologic")
ec.calibrate(RE_vs_RHE=0.72, A_el=0.196)  # a Calculator, saved along with it

tic = time.perf_counter()
ec_id = ec.save()
toc = time.perf_counter()
print(f"Saved {ec!r} with all its data series in {(toc - tic) * 1000:.1f} ms")

# --------------------------------------------------------------------------- #
part("3. What that created: a real relational database in a single file")
# --------------------------------------------------------------------------- #

print(f"Everything is in the one file {DB_FILE}")
print(f"  ({DB_FILE.stat().st_size / 1024:.0f} kB - open it with any SQLite tool,")
print("   e.g. DB Browser for SQLite, DBeaver, datasette, or the sqlite3 CLI)\n")

# The tables were created on demand, derived from the ixdat classes involved:
con = sqlite3.connect(DB_FILE)  # a plain SQL connection, to peek behind the scenes
tables = [
    row[0]
    for row in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    )
]
for table in tables:
    count = con.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
    print(f"  table {table:.<28} {count:>3} rows")

print("\nWith honest-to-goodness foreign keys, e.g.:\n")
print(
    con.execute(
        "SELECT sql FROM sqlite_master WHERE name='measurement_series'"
    ).fetchone()[0]
)
violations = con.execute("PRAGMA foreign_key_check").fetchall()
print(f"\nForeign key violations in the database: {violations or 'none'}")

# --------------------------------------------------------------------------- #
part("4. Loading is lazy: metadata now, numbers only when you need them")
# --------------------------------------------------------------------------- #

tic = time.perf_counter()
loaded_ec = Measurement.get(ec_id)
toc = time.perf_counter()
print(f"Measurement.get({ec_id}) took {(toc - tic) * 1000:.2f} ms")
print("  ... because it read *no* numerical data:")
series = loaded_ec.series_list[0]
print(f"  {series.name!r}._data is {series._data}")

tic = time.perf_counter()
n_points = len(series.data)  # NOW the array is fetched from the database
toc = time.perf_counter()
print(
    f"  first access to .data fetched {n_points} points "
    f"in {(toc - tic) * 1000:.2f} ms"
)

# --------------------------------------------------------------------------- #
part("5. Objects come back by name, as what they were, calibration included")
# --------------------------------------------------------------------------- #

by_name = Measurement.load("Pt_poly_cv_CUT.mpt")
print(f"Measurement.load('Pt_poly_cv_CUT.mpt') -> {by_name!r}")
print(f"  class: {type(by_name).__name__} (not just a generic Measurement)")
print(f"  calculators: {by_name.calculator_list}")
print(f"  equal to what was saved:  by_name == ec  is  {by_name == ec}")

t, v_raw = by_name.grab("raw_potential")
t, v = by_name.grab("potential")  # calibrated, thanks to the reloaded ECCalibration
print(
    f"  mean(potential - raw_potential) = {(v - v_raw).mean():.2f} V"
    "  <- the RE_vs_RHE=0.72 calibration survived the round trip"
)

# --------------------------------------------------------------------------- #
part("6. Composed measurements share rows - data is never duplicated")
# --------------------------------------------------------------------------- #

composed = ec.select(cycle=1) + ec.select(cycle=3)
composed_id = composed.save()
print(f"Saved composed measurement (id={composed_id}) built from two cycles.\n")

n_references, n_series = con.execute(
    "SELECT COUNT(*), COUNT(DISTINCT data_series_id) FROM measurement_series"
).fetchone()
n_measurements = con.execute("SELECT COUNT(*) FROM measurement").fetchone()[0]
print(f"  {n_measurements} measurements now reference series "
      f"{n_references} times ...")
print(f"  ... but only {n_series} distinct data series (arrays) are stored:")
print("  the composed measurement and its components *share* rows.")

reloaded_composed = Measurement.get(composed_id)
print(f"\nRound trip of the composed measurement: "
      f"{reloaded_composed == composed}")

# --------------------------------------------------------------------------- #
part("7. Not just measurements: every ixdat data family gets its tables")
# --------------------------------------------------------------------------- #

xrd = Spectrum.read(TEST_DATA_DIR / "xrd/twotheta_2th.xy", reader="xrdxy")
xrd_id = xrd.save()
loaded_xrd = type(xrd).get(xrd_id)
print(f"Saved and reloaded {loaded_xrd!r}")
print(f"  class: {type(loaded_xrd).__name__}, equal: {loaded_xrd == xrd}")
print("  (spectra live in their own tables: multispectrum, data_series, "
      "field_axes, ...)")

# --------------------------------------------------------------------------- #
part("8. The database speaks SQL: instant analytics over everything saved")
# --------------------------------------------------------------------------- #

print("All measurements, joined with their EC details and series counts:\n")
# NOTE: this join only finds pure ECMeasurement rows, like the ones saved above.
# An ECMSMeasurement writes its ec_technique to "ecms_measurements" instead of
# "ec_measurements" - see the "Known limitation" section of
# docs/source/diving_deeper/backend.rst - so a query meant to cover every EC-ish
# measurement should filter on measurement.technique or UNION both tables.
print(
    pd.read_sql_query(
        """
        SELECT m.id, m.name, m.technique, e.ec_technique,
               COUNT(ms.data_series_id) AS n_series
        FROM measurement m
        LEFT JOIN ec_measurements e ON e.id = m.id
        LEFT JOIN measurement_series ms ON ms.measurement_id = m.id
        GROUP BY m.id ORDER BY m.id
        """,
        con,
    ).to_string(index=False)
)

print("\nAll calculators, joined with their calibration constants:\n")
print(
    pd.read_sql_query(
        """
        SELECT c.id, c.name, c.calculator_type, cal.RE_vs_RHE, cal.A_el, cal.R_Ohm
        FROM calculator c
        LEFT JOIN ec_calibration cal ON cal.id = c.id
        """,
        con,
    ).to_string(index=False)
)
print("\n(Anything - pandas, R, LibreOffice, grafana - can query this file.)")

# --------------------------------------------------------------------------- #
part("9. Rows can be updated in place")
# --------------------------------------------------------------------------- #

loaded_ec.name = "Pt_poly_cv_CUT.mpt (analyzed 2026-07-02)"
DB.backend.save(loaded_ec, force=True)
print(f"Renamed and updated row id={loaded_ec.id}. The database now says:")
print("  " + str(con.execute(
    "SELECT name FROM measurement WHERE id = ?", (loaded_ec.id,)
).fetchone()[0]))

# --------------------------------------------------------------------------- #
part("10. The whole schema of ixdat's data model, generated from the classes")
# --------------------------------------------------------------------------- #

report = DB.backend.schema_report()
schema_file = DEMO_DIR / "ixdat_schema.sql"
schema_file.write_text(report)
n_tables = report.count("CREATE TABLE")
print(f"schema_report() derived {n_tables} tables from the imported ixdat classes")
print(f"  (full DDL written to {schema_file})")
print("\nA taste - inheritance becomes extension tables joined on id:\n")
for statement in report.split(";"):
    if '"ecms_measurements"' in statement or '"tstamps"' in statement:
        print(statement.strip() + ";\n")

# --------------------------------------------------------------------------- #
part("Finally: plots, straight from the database")
# --------------------------------------------------------------------------- #

axes = Measurement.load("Pt_poly_cv_CUT.mpt (analyzed 2026-07-02)").plot_measurement()
axes[0].set_title("EC measurement reloaded from SQLite")

ax = type(xrd).get(xrd_id).plot()
ax.set_title("XRD spectrum reloaded from SQLite")

print(f"\nDone. The database file is yours to explore: {DB_FILE}")

if not plt.isinteractive():
    plt.show()

con.close()
DB.backend.close()
