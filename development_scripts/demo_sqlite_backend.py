"""A tour of ixdat's SQLite database backend.

This demo saves measurements and spectra to a relational database in a single
local file and walks through the main persistence features:
"""

import shutil
import sqlite3
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ixdat import Measurement, Spectrum
from ixdat.data_series import DataSeries, Field, TimeSeries
from ixdat.db import DB, change_database
from ixdat.spectra import SpectrumSeries
from ixdat.techniques.spectroelectrochemistry import ECOpticalMeasurement

plt.close("all")

TEST_DATA_DIR = Path(__file__).absolute().parent.parent / "test_data"

DEMO_DIR = Path(tempfile.gettempdir()) / "ixdat_sqlite_demo"
# Start with an empty folder so every run produces the same ids and row counts.
shutil.rmtree(DEMO_DIR, ignore_errors=True)
DB_FILE = DEMO_DIR / "demo_project.sqlite"


def part(title):
    print("\n" + "=" * 75 + f"\n{title}\n" + "=" * 75 + "\n")


part("1. One line to point ixdat at a database file")

change_database("sqlite", db_path=DB_FILE)
print(f"Active backend: {DB.backend}")

part("2. Saving a measurement (and its calibration) is one line each")

ec = Measurement.read(TEST_DATA_DIR / "biologic/Pt_poly_cv_CUT.mpt", reader="biologic")
ec.calibrate(RE_vs_RHE=0.72, A_el=0.196)  # a Calculator, saved along with it

tic = time.perf_counter()
ec_id = ec.save()
toc = time.perf_counter()
print(
    f"Saved {ec!r} in {DB_FILE} with all its data series in {(toc - tic) * 1000:.1f} ms"
)
print(f"  ({DB_FILE.stat().st_size / 1024:.0f} kB)")

# A database table is much like one sheet in a spreadsheet. Each row holds one
# saved item. ixdat reads the class definitions and makes the needed sheets when
# an item is first saved. The table plan stays beside the Python code.
#
# sqlite3 is Python's built-in way to look directly inside the database file.
con = sqlite3.connect(DB_FILE)
tables = [
    row[0]
    for row in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    )
]
for table in tables:
    count = con.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
    print(f"  table {table:.<28} {count:>3} rows")

# A foreign key is an id that points to a row in another table. For example, a
# measurement stores the ids of its data series. SQLite checks these pointers,
# which prevents a measurement from pointing to a row that does not exist.
print("\nWith honest-to-goodness foreign keys, e.g.:\n")
print(
    con.execute(
        "SELECT sql FROM sqlite_master WHERE name='measurement_series'"
    ).fetchone()[0]
)
violations = con.execute("PRAGMA foreign_key_check").fetchall()
print(f"\nForeign key violations in the database: {violations or 'none'}")

# The schema version says which table layout the file uses. ixdat checks it
# before reading the file. An index is a small lookup guide that helps SQLite
# find common values, such as a name or technique, without scanning every row.
schema_version = con.execute(
    'SELECT value FROM "_ixdat_metadata" WHERE key = ?', ("schema_version",)
).fetchone()[0]
lookup_indexes = con.execute(
    """
    SELECT COUNT(*) FROM sqlite_master
    WHERE type = 'index' AND name LIKE 'ixdat_%'
    """
).fetchone()[0]
print(
    f"Storage schema version: {schema_version}; lookup indexes created: {lookup_indexes}"
)

part("3. Loading is lazy: metadata now, numbers only when you need them")

# Names and other small details are cheap to load. Numerical arrays can be much
# larger, so ixdat leaves them in the database until code asks for `.data`.
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

part("4. Objects come back by name, as what they were, calibration included")

# The saved technique tells ixdat which kind of Measurement to rebuild. Links
# stored in the database also reconnect the measurement to its calibration.
by_name = Measurement.load("Pt_poly_cv_CUT.mpt")
print(f"Measurement.load('Pt_poly_cv_CUT.mpt') -> {by_name!r}")
print(f"  concrete class: {type(by_name).__name__}")
print(f"  calculators: {by_name.calculator_list}")
print(f"  equal to what was saved:  by_name == ec  is  {by_name == ec}")

t, v_raw = by_name.grab("raw_potential")
t, v = by_name.grab("potential")  # calibrated, thanks to the reloaded ECCalibration
print(
    f"  mean(potential - raw_potential) = {(v - v_raw).mean():.2f} V"
    "  <- the RE_vs_RHE=0.72 calibration survived the round trip"
)

part("5. Composed measurements reuse their component data rows")

composed = ec.select(cycle=1) + ec.select(cycle=3)
composed_id = composed.save()
print(f"Saved composed measurement (id={composed_id}) built from two cycles.\n")

# The measurement_series table is a list of pointers from measurements to data
# series. Several measurements can point to the same series row. Large arrays
# are then stored once even when they are used in several places.
n_references, n_series = con.execute(
    "SELECT COUNT(*), COUNT(DISTINCT data_series_id) FROM measurement_series"
).fetchone()
n_measurements = con.execute("SELECT COUNT(*) FROM measurement").fetchone()[0]
print(f"  {n_measurements} measurements now reference series {n_references} times ...")
print(f"  ... but only {n_series} distinct data series (arrays) are stored:")
print("  the composed measurement and its components *share* rows.")

reloaded_composed = Measurement.get(composed_id)
print(f"\nRound trip of the composed measurement: {reloaded_composed == composed}")

part("6. Spectra and other ixdat data families use the same backend")

# A MultiSpectrum groups several spectral fields that share one x-axis. The
# database keeps pointers to those fields in the same order as the Python list.
xrd = Spectrum.read(TEST_DATA_DIR / "xrd/twotheta_2th.xy", reader="xrdxy")
xrd_id = xrd.save()
loaded_xrd = type(xrd).get(xrd_id)
print(
    f"XRD: {type(loaded_xrd).__name__}, {len(loaded_xrd.fields)} field(s), "
    f"equal after reload: {loaded_xrd == xrd}"
)

# A SpectrumSeries is a grid of numbers. One axis is time and the other is
# wavelength. The Field holds the grid and pointers to those two axes.
wavelength = DataSeries(
    name="wavelength / nm",
    unit_name="nm",
    data=np.linspace(400, 700, 4),
)
optical_time = TimeSeries(
    name="optical time / s",
    unit_name="s",
    data=np.array([0.0, 1.0, 2.0]),
    tstamp=1.6e9,
)
optical_series = SpectrumSeries(
    name="time-resolved optical spectra",
    technique="EC-Optical",
    tstamp=optical_time.tstamp,
    field=Field(
        name="intensity",
        unit_name="counts",
        data=np.array(
            [
                [1.0, 2.0, 3.0, 4.0],
                [2.0, 3.0, 4.0, 5.0],
                [3.0, 4.0, 5.0, 6.0],
            ]
        ),
        axes_series=[optical_time, wavelength],
    ),
)
reference = Spectrum(
    name="optical reference",
    technique="optical",
    tstamp=optical_time.tstamp,
    field=Field(
        name="reference intensity",
        unit_name="counts",
        data=np.ones(4),
        axes_series=[wavelength],
    ),
)

# This measurement owns the changing SpectrumSeries above and points to one
# separate Spectrum used as its reference. The database saves the whole group
# and remembers how its pieces connect.
optical = ECOpticalMeasurement(
    name="synthetic EC-Optical measurement",
    technique="EC-Optical",
    ec_technique="synthetic potential sweep",
    tstamp=optical_time.tstamp,
    series_list=[optical_time],
    spectrum_series=optical_series,
    reference_spectrum=reference,
)
optical_id = optical.save()
loaded_optical = Measurement.get(optical_id)
lazy_optical_field = loaded_optical.spectrum_series.field
print(
    f"EC-Optical: {type(loaded_optical).__name__}, "
    f"linked {type(loaded_optical.spectrum_series).__name__} and "
    f"{type(loaded_optical.reference_spectrum).__name__}"
)
print(f"  spectrum array before access: {lazy_optical_field._data}")
print(f"  spectrum array after access:  shape={loaded_optical.spectra.data.shape}")

# These tables hold the connections described above. The list relationships store a
# position for each item. The single reference stores one spectrum id.
print(
    "  relational paths: multispectrum_fields for XRD fields, field_axes for "
    "spectral axes, and ec_optical_measurements for the reference"
)

part("7. A complete object graph is one transaction")

# A transaction makes a save all-or-nothing. This example has metadata that
# cannot be saved as JSON. SQLite removes the child series it had started to
# save, so no loose or half-saved rows remain.
rows_before = con.execute('SELECT COUNT(*) FROM "data_series"').fetchone()[0]
invalid = Measurement(
    name="this graph will roll back",
    technique="simple",
    metadata={"not JSON": object()},
    series_list=[DataSeries(name="unsaved child", unit_name="V", data=np.array([1.0]))],
)
try:
    invalid.save()
except TypeError:
    pass
rows_after = con.execute('SELECT COUNT(*) FROM "data_series"').fetchone()[0]
print(
    "An invalid parent rolled back itself and its child: "
    f"data-series rows stayed at {rows_before} ({rows_after == rows_before})."
)

part("8. The database speaks SQL: instant analytics over everything saved")

# JOIN lines up rows from different tables using their ids. LEFT JOIN keeps a
# measurement in the result even when it has no matching EC detail or series.
# GROUP BY gathers all series links for one measurement so COUNT can total them.
print("All measurements, joined with their EC details and series counts:\n")
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

part("9. Rows can be updated in place")

# `force=True` writes the changed values back to the existing row with this id.
loaded_ec.name = "Pt_poly_cv_CUT.mpt (analyzed 2026-07-02)"
DB.backend.save(loaded_ec, force=True)
print(f"Renamed and updated row id={loaded_ec.id}. The database now says:")
updated_name = con.execute(
    "SELECT name FROM measurement WHERE id = ?", (loaded_ec.id,)
).fetchone()[0]
print(f"  {updated_name}")

part("10. The whole schema of ixdat's data model, generated from the classes")

# A schema is the database's blueprint: its tables, columns, and links. ixdat
# builds this blueprint from the save information already written on its Python
# classes. A plugin class can describe its own saved values in the same way.
report = DB.backend.schema_report()
schema_file = DEMO_DIR / "ixdat_schema.sql"
schema_file.write_text(report)
n_tables = report.count("CREATE TABLE")
print(f"schema_report() derived {n_tables} tables from the imported ixdat classes")
print(f"  (full DDL written to {schema_file})")
print("\nA taste - inheritance becomes extension tables joined on id:\n")

# In Python, a TimeSeries is a DataSeries with an extra timestamp. SQL stores
# the shared DataSeries values in one table and the timestamp in a small extra
# table. Both rows use the same id, so SQLite knows they describe one object.
# EC-MS measurements follow the same pattern for their extra values.
for statement in report.split(";"):
    if '"ecms_measurements"' in statement or '"tstamps"' in statement:
        print(statement.strip() + ";\n")

part("11. Close and reopen the complete project")

# Closing ends the current connection to the file. The saved rows stay on disk.
# Loading from the same file rebuilds each requested object and its links.
con.close()
DB.backend.close()
change_database("sqlite", db_path=DB_FILE)

reopened_ec = Measurement.load("Pt_poly_cv_CUT.mpt (analyzed 2026-07-02)")
reopened_xrd = type(xrd).get(xrd_id)
reopened_optical = Measurement.get(optical_id)
print(
    f"Reopened {type(reopened_ec).__name__}, {type(reopened_xrd).__name__}, "
    f"and {type(reopened_optical).__name__} from {DB_FILE.name}"
)

part("Finally: plots, straight from the reopened database")

axes = reopened_ec.plot_measurement()
axes[0].set_title("EC measurement reloaded from SQLite")

ax = reopened_xrd.plot()
ax.set_title("XRD spectrum reloaded from SQLite")

ax = reopened_optical.spectrum_series.heat_plot()
ax.set_title("SpectrumSeries reloaded through an EC-Optical measurement")

print(f"\nDone. The database file is yours to explore: {DB_FILE}")

if not plt.isinteractive():
    plt.show()

DB.backend.close()
