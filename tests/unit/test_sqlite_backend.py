"""Unit tests for the relational schema derivation and the SQLite backend"""

import sqlite3

import numpy as np
import pytest

from ixdat import Measurement, Spectrum
from ixdat.backends import relational
from ixdat.backends.directory_backend import DirBackend
from ixdat.backends.relational import ColumnSchema
from ixdat.backends.sqlite_backend import METADATA_TABLE, SCHEMA_VERSION, SQLiteBackend
from ixdat.data_series import DataSeries, Field, TimeSeries, ValueSeries
from ixdat.db import DB, change_database
from ixdat.exceptions import DataBaseError
from ixdat.measurement_base import Calculator
from ixdat.calculators.ec_calculators import ECCalibration
from ixdat.spectra import MultiSpectrum
from ixdat.techniques.ec_ms import ECMSMeasurement


@pytest.fixture
def sqlite_backend(tmp_path):
    """An SQLite backend on a fresh database file, set as the active backend"""
    original_backend = DB.backend
    change_database("sqlite", db_path=tmp_path / "test.sqlite")
    yield DB.backend
    DB.backend.close()
    DB.set_backend(original_backend)


class TestRelationalSchema:
    """Tests for the derivation of table descriptions from Saveable classes"""

    def test_measurement_main_table(self):
        schema = relational.main_table_schema(Measurement)
        assert schema.name == "measurement"
        assert set(schema.column_names) == {
            "id",
            "name",
            "technique",
            "metadata",
            "aliases",
            "sample_name",
            "tstamp",
        }
        # deterministic order, with id and name first:
        assert schema.column_names[:2] == ["id", "name"]

    def test_extension_table(self):
        (schema,) = relational.extension_table_schemas(ECMSMeasurement)
        assert schema.name == "ecms_measurements"
        assert schema.extends == "measurement"
        # the extension table's id is a foreign key to the main table:
        assert schema.columns[0].name == "id"
        assert schema.columns[0].foreign_key == ("measurement", "id")
        assert {column.name for column in schema.data_columns} == {
            "ec_technique",
            "tspan_bg",
        }

    def test_multiply_inheriting_class_only_gets_its_own_extension_table(self):
        """Known limitation: `extra_column_attrs` is shadowed, not merged, across
        multiple inheritance (see the "Known limitation" section of
        docs/source/diving_deeper/backend.rst). `ECMSMeasurement` inherits from both
        `ECMeasurement` and `MSMeasurement`, but its own `extra_column_attrs`
        replaces `ECMeasurement`'s rather than adding to it, so `ec_measurements` -
        the extension table `ECMeasurement` itself declares - is absent here.
        """
        extension_table_names = {
            schema.name
            for schema in relational.extension_table_schemas(ECMSMeasurement)
        }
        assert extension_table_names == {"ecms_measurements"}
        assert "ec_measurements" not in extension_table_names

    def test_linker_tables(self):
        linkers = {
            schema.name: schema
            for schema in relational.linker_table_schemas(Measurement)
        }
        series_linker = linkers["measurement_series"]
        assert series_linker.owner_column == "measurement_id"
        assert series_linker.linked_column == "data_series_id"
        assert series_linker.id_attr == "s_ids"
        assert series_linker.is_list
        # a self-referencing linker table gets distinguishable column names:
        component_linker = linkers["component_measurements"]
        assert component_linker.owner_column == "measurement_id"
        assert component_linker.linked_column == "linked_measurement_id"

    def test_column_types_come_from_the_declaring_class(self):
        """Each column's type is the one declared by the class which adds it"""
        columns = {
            column.name: column
            for column in relational.main_table_schema(DataSeries).columns
        }
        assert columns["data"].dtype == "NDARRAY"  # declared by DataSeries
        assert columns["name"].dtype == "TEXT"  # inherited from Saveable

    def test_column_types_are_merged_over_multiple_inheritance(self):
        """Unlike extra_column_attrs, column_types of *both* parents carry over"""
        types = ECMSMeasurement.get_column_types()
        assert types["ec_technique"] == "TEXT"  # from ECMeasurement
        assert types["tspan_bg"] == "JSON"  # from ECMSMeasurement itself
        assert types["aliases"] == "JSON"  # from Measurement
        # and a class only has to declare what it adds itself:
        assert ECMSMeasurement.__dict__["column_types"] == {"tspan_bg": "JSON"}

    def test_column_references_become_foreign_keys(self):
        columns = {
            column.name: column
            for column in relational.main_table_schema(Spectrum).columns
        }
        assert columns["field_id"].foreign_key == ("data_series", "id")
        # a column referring to another row always holds that row's integer id:
        assert columns["field_id"].dtype == "INTEGER"

    def test_a_class_outside_ixdat_can_declare_its_own_column_types(self):
        """A plugin class gets its columns typed without ixdat knowing about it"""

        class PluginMeasurement(Measurement):
            extra_column_attrs = {"plugin_measurements": {"instrument_settings"}}
            column_types = {"instrument_settings": "JSON"}

        (schema,) = relational.extension_table_schemas(PluginMeasurement)
        (column,) = schema.data_columns
        assert column.name == "instrument_settings"
        assert column.dtype == "JSON"

    def test_unknown_column_type_is_rejected(self):
        # a stand-in rather than a real Saveable subclass, since a class with a
        # broken column type would stay in Saveable.__subclasses__() for the rest
        # of the session and break the schema of every test after this one:
        class MistypedClass:
            table_name = "mistyped"
            column_attrs = {"whoops"}

            @classmethod
            def get_column_types(cls):
                return {"whoops": "NDARAY"}  # typo for "NDARRAY"

            @classmethod
            def get_column_references(cls):
                return {}

        with pytest.raises(DataBaseError, match="unknown type"):
            relational.main_table_schema(MistypedClass)

    def test_family_table_schemas(self):
        """All classes sharing a main table contribute to the family schema"""
        extension_schemas, linker_schemas = relational.family_table_schemas(
            DataSeries
        )
        # TimeSeries extends the data_series table with the tstamps table:
        assert "tstamps" in {schema.name for schema in extension_schemas}
        # Field links data_series rows to their axes' data_series rows:
        assert "field_axes" in {schema.name for schema in linker_schemas}


class TestSQLiteBackend:
    """Tests saving and loading ixdat objects with the SQLite backend"""

    def test_value_series_round_trip(self, sqlite_backend):
        tseries = TimeSeries(
            name="time / s", unit_name="s", data=np.array([0.0, 1.0, 2.0]), tstamp=1.6e9
        )
        vseries = ValueSeries(
            name="potential / V",
            unit_name="V",
            data=np.array([1.0, 2.0, 4.0]),
            tseries=tseries,
        )
        i = vseries.save()
        loaded = DataSeries.get(i)
        assert isinstance(loaded, ValueSeries)
        assert loaded.unit_name == "V"
        assert np.allclose(loaded.data, vseries.data)
        # the linked TimeSeries is also loaded from the database:
        assert loaded.tseries.tstamp == tseries.tstamp
        assert np.allclose(loaded.t, tseries.data)

    def test_data_is_loaded_lazily(self, sqlite_backend):
        tseries = TimeSeries(
            name="time / s", unit_name="s", data=np.array([0.0, 1.0]), tstamp=1.6e9
        )
        i = tseries.save()
        loaded = DataSeries.get(i)
        assert loaded._data is None  # the array was not fetched by get()...
        assert np.allclose(loaded.data, tseries.data)  # ...but is available on demand
        assert loaded._data is not None

    def test_spectrum_round_trip(self, sqlite_backend):
        xseries = DataSeries(
            name="wavelength / nm", unit_name="nm", data=np.linspace(200, 800, 5)
        )
        field = Field(
            name="intensity",
            unit_name="counts",
            data=np.array([1.0, 2.0, 5.0, 2.0, 1.0]),
            axes_series=[xseries],
        )
        spectrum = Spectrum(
            name="test spectrum",
            technique="spectrum",
            tstamp=1.6e9,
            field=field,
            metadata={"scans": 2},
        )
        i = spectrum.save()
        loaded = Spectrum.get(i)
        assert loaded == spectrum
        assert loaded.metadata == {"scans": 2}  # JSON column round trip
        assert np.allclose(loaded.x, spectrum.x)
        assert np.allclose(loaded.y, spectrum.y)

    def test_calculator_round_trip(self, sqlite_backend):
        calibration = ECCalibration(name="my calibration", RE_vs_RHE=0.72, A_el=0.196)
        i = calibration.save()
        loaded = Calculator.get(i)
        # from_dict dispatches on the calculator_type column:
        assert isinstance(loaded, ECCalibration)
        assert loaded.RE_vs_RHE == 0.72
        assert loaded.A_el == 0.196
        assert loaded.R_Ohm is None

    def test_multiply_inheriting_measurement_round_trips_but_skips_parent_table(
        self, sqlite_backend
    ):
        """Known limitation, documented in docs/source/diving_deeper/backend.rst:
        `ECMSMeasurement` (inherits from `ECMeasurement` and `MSMeasurement`) writes
        only to its own extension table. Round-tripping the object is unaffected -
        `ec_technique` survives via `ecms_measurements` - but a query joining only
        `ec_measurements` would miss this row, since no row is written there.
        """
        tseries = TimeSeries(
            name="t", unit_name="s", data=np.array([0.0, 1.0]), tstamp=1.6e9
        )
        ecms = ECMSMeasurement(
            name="synthetic ecms",
            technique="EC-MS",
            tstamp=1.6e9,
            series_list=[tseries],
            ec_technique="Cyclic Voltammetry Advanced",
        )
        i = ecms.save()

        loaded = Measurement.get(i)
        assert isinstance(loaded, ECMSMeasurement)
        assert loaded.ec_technique == "Cyclic Voltammetry Advanced"  # round trip OK

        existing_tables = sqlite_backend._existing_tables()
        assert "ec_measurements" not in existing_tables  # never created: no writer
        row = sqlite_backend.connection.execute(
            'SELECT 1 FROM "ecms_measurements" WHERE "id" = ?', (i,)
        ).fetchone()
        assert row is not None  # the data lives here instead

    def test_load_returns_newest_with_name(self, sqlite_backend):
        DataSeries(name="twin", unit_name="V", data=np.array([1.0])).save()
        i_newest = DataSeries(name="twin", unit_name="A", data=np.array([2.0])).save()
        loaded = DataSeries.load("twin")
        assert loaded.id == i_newest
        assert loaded.unit_name == "A"
        with pytest.raises(DataBaseError):
            DataSeries.load("no series has this name")

    def test_update_with_force(self, sqlite_backend):
        series = DataSeries(name="before", unit_name="V", data=np.array([1.0]))
        i = series.save()
        series.name = "after"
        assert sqlite_backend.save(series, force=True) == i
        assert DataSeries.get(i).name == "after"

    def test_get_missing_row_raises(self, sqlite_backend):
        with pytest.raises(DataBaseError):
            Measurement.get(999)

    def test_codecs(self, sqlite_backend):
        json_column = ColumnSchema("metadata", dtype="JSON")
        value = {"a": np.int64(1), "b": np.array([1.5, 2.5])}
        encoded = sqlite_backend._encode(json_column, value)
        assert sqlite_backend._decode(json_column, encoded) == {
            "a": 1,
            "b": [1.5, 2.5],
        }

        array_column = ColumnSchema("data", dtype="NDARRAY")
        data = np.linspace(0, 1, 5)
        encoded = sqlite_backend._encode(array_column, data)
        assert np.array_equal(sqlite_backend._decode(array_column, encoded), data)

        # numpy scalars are stored as their python equivalents:
        assert sqlite_backend._encode(ColumnSchema("x"), np.float64(1.5)) == 1.5

    def test_dynamic_column_rejects_containers(self, sqlite_backend):
        """Containers require a logical type declared in the class's column_types"""
        with pytest.raises(DataBaseError):
            sqlite_backend._encode(ColumnSchema("unregistered"), {"nested": "dict"})

    def test_multispectrum_round_trip(self, sqlite_backend):
        xseries = DataSeries(
            name="two theta / deg", unit_name="deg", data=np.linspace(20, 80, 7)
        )
        fields = [
            Field(
                name=name,
                unit_name="counts",
                data=np.linspace(0, 1, 7) * scale,
                axes_series=[xseries],
            )
            for name, scale in (("intensity", 100.0), ("error", 10.0))
        ]
        multi = MultiSpectrum(
            name="test multispectrum", technique="XRD", tstamp=1.6e9, fields=fields
        )
        i = multi.save()
        loaded = MultiSpectrum.get(i)
        assert loaded == multi
        assert len(loaded.fields) == 2
        assert np.allclose(loaded.x, multi.x)

    def test_schema_report(self, sqlite_backend):
        report = sqlite_backend.schema_report()
        assert 'CREATE TABLE IF NOT EXISTS "measurement"' in report
        assert 'CREATE TABLE IF NOT EXISTS "measurement_series"' in report
        assert 'REFERENCES "data_series"("id")' in report

    def test_complete_graph_save_is_atomic(self, sqlite_backend):
        series = TimeSeries(
            name="time / s", unit_name="s", data=np.array([0.0]), tstamp=0.0
        )
        measurement = Measurement(
            name="invalid metadata",
            technique="simple",
            metadata={"not_json": object()},
            series_list=[series],
        )
        original_backend = series._backend
        original_id = series._id

        with pytest.raises(TypeError):
            measurement.save()

        assert sqlite_backend.connection.execute(
            'SELECT COUNT(*) FROM "data_series"'
        ).fetchone()[0] == 0
        assert sqlite_backend.connection.execute(
            'SELECT COUNT(*) FROM "measurement"'
        ).fetchone()[0] == 0
        assert series._backend is original_backend
        assert series._id is original_id

    def test_reopen_and_lazy_load_from_inactive_backend(self, sqlite_backend):
        series = DataSeries(
            name="persistent", unit_name="V", data=np.array([1.0, 2.0])
        )
        measurement = Measurement(
            name="persistent graph", technique="simple", series_list=[series]
        )
        i = measurement.save()
        db_path = sqlite_backend.db_path
        original_active_backend = DB.backend
        sqlite_backend.close()

        with SQLiteBackend(db_path=db_path) as reopened:
            loaded = Measurement.get(i, backend=reopened)
            assert DB.backend is original_active_backend
            loaded_series = loaded.series_list[0]
            assert DB.backend is original_active_backend
            assert loaded_series._data is None
            assert np.array_equal(loaded_series.data, series.data)

    def test_two_connections_share_a_database(self, sqlite_backend):
        with SQLiteBackend(db_path=sqlite_backend.db_path) as second_backend:
            first_id = sqlite_backend.save(
                DataSeries(name="first", unit_name="V", data=np.array([1.0]))
            )
            second_id = second_backend.save(
                DataSeries(name="second", unit_name="A", data=np.array([2.0]))
            )
            assert first_id != second_id
            assert DataSeries.get(second_id, backend=sqlite_backend).name == "second"
            assert sqlite_backend.connection.execute(
                "PRAGMA foreign_key_check"
            ).fetchall() == []

    def test_in_memory_database_uses_sqlite_convention(self):
        with SQLiteBackend(db_path=":memory:") as first_backend:
            i = first_backend.save(
                DataSeries(name="temporary", unit_name="V", data=np.array([1.0]))
            )
            assert DataSeries.get(i, backend=first_backend).name == "temporary"
            with SQLiteBackend(db_path=":memory:") as second_backend:
                assert first_backend != second_backend
                assert not second_backend.contains("data_series", i)

    def test_additive_schema_migration(self, tmp_path):
        db_path = tmp_path / "legacy.sqlite"
        connection = sqlite3.connect(db_path)
        connection.execute(
            'CREATE TABLE "data_series" '
            '("id" INTEGER PRIMARY KEY, "name" TEXT)'
        )
        connection.close()

        with SQLiteBackend(db_path=db_path) as backend:
            series = DataSeries(
                name="after migration", unit_name="A", data=np.array([1.0])
            )
            backend.save(series)
            columns = {
                row[1] for row in backend.connection.execute(
                    'PRAGMA table_info("data_series")'
                )
            }
            assert {"id", "name", "unit_name", "series_type", "data"} <= columns
            assert backend.connection.execute(
                f'SELECT "value" FROM "{METADATA_TABLE}" WHERE "key" = ?',
                ("schema_version",),
            ).fetchone()[0] == str(SCHEMA_VERSION)

    def test_incompatible_schema_requires_manual_migration(self, tmp_path):
        db_path = tmp_path / "incompatible.sqlite"
        connection = sqlite3.connect(db_path)
        connection.execute(
            'CREATE TABLE "data_series" ("id" TEXT PRIMARY KEY, "name" TEXT)'
        )
        connection.close()

        with SQLiteBackend(db_path=db_path) as backend:
            with pytest.raises(DataBaseError, match="manual database migration"):
                backend.save(
                    DataSeries(name="x", unit_name="V", data=np.array([1.0]))
                )

    def test_newer_schema_version_is_rejected(self, tmp_path):
        db_path = tmp_path / "future.sqlite"
        connection = sqlite3.connect(db_path)
        connection.execute(
            f'CREATE TABLE "{METADATA_TABLE}" '
            '("key" TEXT PRIMARY KEY, "value" TEXT NOT NULL)'
        )
        connection.execute(
            f'INSERT INTO "{METADATA_TABLE}" VALUES (?, ?)',
            ("schema_version", str(SCHEMA_VERSION + 1)),
        )
        connection.commit()
        connection.close()

        with pytest.raises(DataBaseError, match="newer ixdat schema"):
            SQLiteBackend(db_path=db_path)

    def test_failed_explicit_backend_load_restores_active_backend(
        self, sqlite_backend, tmp_path
    ):
        original_backend = DB.backend
        with SQLiteBackend(db_path=tmp_path / "other.sqlite") as other_backend:
            with pytest.raises(DataBaseError):
                DataSeries.load("missing", backend=other_backend)
            assert DB.backend is original_backend

    def test_lookup_indexes_are_created(self, sqlite_backend):
        sqlite_backend._ensure_tables(Measurement)
        indexes = {
            row[1]
            for row in sqlite_backend.connection.execute(
                'PRAGMA index_list("measurement")'
            )
        }
        assert "ixdat_measurement_name_id" in indexes
        assert "ixdat_measurement_technique" in indexes


def test_directory_load_uses_exact_unescaped_name(tmp_path):
    original_backend = DB.backend
    backend = DirBackend(directory=tmp_path, project_name="name_collision")
    DB.set_backend(backend)
    try:
        slash_id = DataSeries(
            name="a/b", unit_name="V", data=np.array([1.0])
        ).save()
        DataSeries(name="a_DIV_b", unit_name="A", data=np.array([2.0])).save()
        assert DataSeries.load("a/b").id == slash_id

        measurement = Measurement(
            name="numpy metadata",
            technique="simple",
            metadata={"counts": np.array([1, 2]), "scans": np.int64(2)},
        )
        measurement.save()
        assert Measurement.load(measurement.name).metadata == {
            "counts": [1, 2],
            "scans": 2,
        }
    finally:
        DB.set_backend(original_backend)
