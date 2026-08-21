"""Unit tests for the relational schema derivation and the SQLite backend"""

import sqlite3

import numpy as np
import pytest

from ixdat import Measurement, Spectrum
from ixdat.backends import relational
from ixdat.backends.directory_backend import DirBackend
from ixdat.backends.memory_backend import MemoryBackend
from ixdat.backends.relational import ColumnSchema
from ixdat.backends.sqlite_backend import (
    METADATA_TABLE_NAME,
    SCHEMA_VERSION,
    SQLiteBackend,
)
from ixdat.data_series import DataSeries, Field, TimeSeries, ValueSeries
from ixdat.db import (
    DB,
    PlaceHolderObject,
    Relationship,
    Saveable,
    change_database,
    same_short_identity,
)
from ixdat.exceptions import DataBaseError
from ixdat.measurement_base import Calculator
from ixdat.calculators.ec_calculators import ECCalibration
from ixdat.calculators.ms_calculators import MSCalibration, MSCalResult
from ixdat.spectra import MultiSpectrum
from ixdat.techniques.ec_ms import ECMSMeasurement
from ixdat.techniques.spectroelectrochemistry import ECOpticalMeasurement


def create_schema_metadata(connection, version=SCHEMA_VERSION):
    """Mark a hand-built test database as an ixdat SQLite schema."""
    connection.execute(
        f'CREATE TABLE "{METADATA_TABLE_NAME}" '
        '("key" TEXT PRIMARY KEY, "value" TEXT NOT NULL)'
    )
    connection.execute(
        f'INSERT INTO "{METADATA_TABLE_NAME}" VALUES (?, ?)',
        ("schema_version", str(version)),
    )


@pytest.fixture
def sqlite_backend(tmp_path):
    """An SQLite backend on a fresh database file, set as the active backend"""
    original_backend = DB.backend
    change_database("sqlite", db_path=tmp_path / "test.sqlite")
    yield DB.backend
    DB.backend.close()
    DB.set_backend(original_backend)


class MistypedMeasurement(Measurement):
    """A class whose column type is a typo, used by two tests below

    Defined at module level, since a Saveable subclass stays in
    `Saveable.__subclasses__()` for the rest of the session once it exists. That
    is the point: the tests check that its presence leaves every other
    measurement loadable.
    """

    extra_column_attrs = {"mistyped_measurements": {"whoops"}}
    column_types = {"whoops": "NDARAY"}  # invalid: column types are Python types


class ArrayColumnObject(Saveable):
    """Plugin-like object with a NumPy-array column other than ``data``."""

    table_name = "array_column_objects"
    column_attrs = {"name", "coefficients"}
    column_types = {"coefficients": np.ndarray}

    def __init__(self, name, coefficients):
        super().__init__()
        self.name = name
        self.coefficients = coefficients


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

    def test_extension_tables_merge_over_multiple_inheritance(self):
        """Each ancestor contributes one extension table."""
        schemas = {
            schema.name: schema
            for schema in relational.extension_table_schemas(ECMSMeasurement)
        }
        assert set(schemas) == {"ec_measurements", "ecms_measurements"}
        ecms_schema = schemas["ecms_measurements"]
        assert ecms_schema.base_table == "measurement"
        # the extension table's id is a foreign key to the main table:
        assert ecms_schema.columns[0].name == "id"
        assert ecms_schema.columns[0].foreign_key == ("measurement", "id")
        assert {column.name for column in schemas["ec_measurements"].data_columns} == {
            "ec_technique"
        }
        assert {column.name for column in ecms_schema.data_columns} == {"tspan_bg"}

    def test_collecting_all_attrs_does_not_modify_main_table_metadata(self):
        """Linker attributes stay out of the main table schema."""
        original_column_attrs = Measurement.column_attrs.copy()

        all_attrs = Measurement.get_all_column_attrs()

        assert {"s_ids", "m_ids", "c_ids"} <= all_attrs
        assert Measurement.column_attrs == original_column_attrs
        assert not {"s_ids", "m_ids", "c_ids"} & set(
            relational.main_table_schema(Measurement).column_names
        )

    def test_linker_tables(self):
        linkers = {
            schema.name: schema
            for schema in relational.linker_table_schemas(Measurement)
        }
        series_linker = linkers["measurement_series"]
        assert series_linker.owner_column == "measurement_id"
        assert series_linker.linked_column == "data_series_id"
        assert series_linker.id_attr == "s_ids"
        assert series_linker.many
        # a self-referencing linker table gets distinguishable column names:
        component_linker = linkers["component_measurements"]
        assert component_linker.owner_column == "measurement_id"
        assert component_linker.linked_column == "linked_measurement_id"
        optical_linkers = {
            schema.name: schema
            for schema in relational.linker_table_schemas(ECOpticalMeasurement)
        }
        assert "ec_optical_measurements" not in optical_linkers

    def test_scalar_relationship_in_an_extension_table(self):
        schemas = {
            schema.name: schema
            for schema in relational.extension_table_schemas(ECOpticalMeasurement)
        }
        reference_column = {
            column.name: column
            for column in schemas["ec_optical_measurements"].data_columns
        }["ref_id"]
        assert reference_column.foreign_key == ("spectrums", "id")
        assert "position" not in schemas["ec_optical_measurements"].column_names

    def test_relationship_connects_objects_to_their_storage(self):
        relationship = Measurement.get_relationships()["series_list"]
        assert isinstance(relationship, Relationship)
        assert relationship.id_attr == "s_ids"
        assert relationship.storage_table == "measurement_series"
        assert relationship.many
        assert relationship.save_related

    def test_relationship_many_decides_whether_ids_form_a_list(self):
        """The id attribute's name does not decide the shape of a new relationship."""

        class PluginMeasurement(Measurement):
            relationships = {
                "related_calculators": Relationship(
                    "calculator",
                    "calculator_references",
                    many=True,
                    storage_table="plugin_measurement_calculators",
                )
            }

        linkers = {
            schema.name: schema
            for schema in relational.linker_table_schemas(PluginMeasurement)
        }
        assert linkers["plugin_measurement_calculators"].many

    def test_legacy_linker_names_still_decide_whether_ids_form_a_list(self):
        """Older plugin declarations keep their established ``_ids`` rule."""

        class LegacyPluginMeasurement(Measurement):
            extra_linkers = {
                "legacy_measurement_calculators": ("calculator", "calculator_ids")
            }

        linkers = {
            schema.name: schema
            for schema in relational.linker_table_schemas(LegacyPluginMeasurement)
        }
        assert linkers["legacy_measurement_calculators"].many

    def test_many_relationship_needs_a_storage_table(self):
        with pytest.raises(
            ValueError, match="containing many objects needs a storage_table"
        ):
            Relationship("calculator", "calculator_ids", many=True)

    def test_column_types_merge_from_the_declaring_classes(self):
        """Column types come from every declaring class in the ancestry."""
        columns = {
            column.name: column
            for column in relational.main_table_schema(DataSeries).columns
        }
        assert columns["data"].dtype is np.ndarray  # declared by DataSeries
        assert columns["name"].dtype is str  # inherited from Saveable

        types = ECMSMeasurement.get_column_types()
        assert types["ec_technique"] is str  # from ECMeasurement
        assert types["tspan_bg"] is list  # from ECMSMeasurement itself
        assert types["aliases"] is dict  # from Measurement
        # and a class only has to declare what it adds itself:
        assert ECMSMeasurement.__dict__["column_types"] == {"tspan_bg": list}

    def test_column_references_become_foreign_keys(self):
        columns = {
            column.name: column
            for column in relational.main_table_schema(Spectrum).columns
        }
        assert columns["field_id"].foreign_key == ("data_series", "id")
        # a column referring to another row always holds that row's integer id:
        assert columns["field_id"].dtype is int

    def test_a_class_outside_ixdat_can_declare_its_own_column_types(self):
        """A plugin class gets its columns typed without ixdat knowing about it"""

        class PluginMeasurement(Measurement):
            extra_column_attrs = {"plugin_measurements": {"instrument_settings"}}
            column_types = {"instrument_settings": dict}

        (schema,) = relational.extension_table_schemas(PluginMeasurement)
        (column,) = schema.data_columns
        assert column.name == "instrument_settings"
        assert column.dtype is dict

    def test_unknown_column_type_is_rejected(self):
        with pytest.raises(DataBaseError, match="unknown type"):
            relational.extension_table_schemas(MistypedMeasurement)

    def test_metadata_for_an_unknown_column_is_rejected(self):
        class UnknownColumnMeasurement(Measurement):
            column_references = {"missing_id": "calculator"}

        with pytest.raises(DataBaseError, match="unknown column.*missing_id"):
            relational.main_table_schema(UnknownColumnMeasurement)

    def test_family_table_schemas(self):
        """All classes sharing a main table contribute to the family schema"""
        extension_schemas, linker_schemas = relational.family_table_schemas(DataSeries)
        # TimeSeries extends the data_series table with the tstamps table:
        assert "tstamps" in {schema.name for schema in extension_schemas}
        # Field links data_series rows to their axes' data_series rows:
        assert "field_axes" in {schema.name for schema in linker_schemas}


class TestSQLiteBackend:
    """Tests saving and loading ixdat objects with the SQLite backend"""

    def test_non_data_array_column_loads_with_object(self, sqlite_backend):
        """A differently named NumPy-array column loads with its object."""
        coefficients = np.array([1.0, 0.5, 0.25])
        i = ArrayColumnObject("calibration", coefficients).save()

        loaded = ArrayColumnObject.get(i)

        assert np.array_equal(loaded.coefficients, coefficients)

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
        assert loaded.a_ids == [(sqlite_backend, tseries.id)]
        assert loaded._data is None
        assert np.allclose(loaded.data, vseries.data)
        assert loaded._data is not None
        # the linked TimeSeries is also loaded from the database:
        assert loaded.tseries.tstamp == tseries.tstamp
        assert np.allclose(loaded.t, tseries.data)

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
        assert loaded.field_id == (sqlite_backend, field.id)
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

    def test_multiply_inheriting_measurement_populates_each_extension_table(
        self, sqlite_backend
    ):
        """An EC-MS row carries its EC and EC-MS attributes in separate tables."""
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
        assert loaded.ec_technique == "Cyclic Voltammetry Advanced"

        assert sqlite_backend.connection.execute(
            'SELECT "ec_technique" FROM "ec_measurements" WHERE "id" = ?', (i,)
        ).fetchone() == ("Cyclic Voltammetry Advanced",)
        assert sqlite_backend.connection.execute(
            'SELECT "tspan_bg" FROM "ecms_measurements" WHERE "id" = ?', (i,)
        ).fetchone() == (None,)
        ecms_columns = {
            row[1]
            for row in sqlite_backend.connection.execute(
                'PRAGMA table_info("ecms_measurements")'
            )
        }
        assert "ec_technique" not in ecms_columns

    def test_ms_calibration_results_are_loaded_lazily(self, sqlite_backend):
        """A loaded MS calibration knows its results' id's before loading them"""
        calibration = MSCalibration(
            name="lazy calibration",
            ms_cal_results=[MSCalResult(name="O2 at M32", mol="O2", mass="M32", F=1.5)],
        )
        i = calibration.save()

        loaded = Calculator.get(i)
        assert all(
            isinstance(result, PlaceHolderObject)
            for result in loaded._ms_cal_results  # not fetched by get()...
        )
        assert loaded.ms_cal_result_ids  # ...and not by asking for their id's...
        assert all(
            isinstance(result, PlaceHolderObject) for result in loaded._ms_cal_results
        )
        assert loaded.ms_cal_results[0].F == 1.5  # ...but available on demand
        assert not any(
            isinstance(result, PlaceHolderObject) for result in loaded._ms_cal_results
        )

    def test_a_broken_sibling_class_does_not_block_saving_or_loading(
        self, sqlite_backend
    ):
        """`MistypedMeasurement` shares the "measurement" table with these rows"""
        measurement = Measurement(
            name="healthy",
            technique="simple",
            series_list=[
                TimeSeries(
                    name="t", unit_name="s", data=np.array([0.0, 1.0]), tstamp=1.6e9
                )
            ],
        )
        i = measurement.save()
        assert Measurement.get(i).name == "healthy"
        assert Measurement.load("healthy").name == "healthy"
        # the report covers the healthy classes:
        report = sqlite_backend.schema_report()
        assert 'CREATE TABLE IF NOT EXISTS "measurement"' in report
        assert "mistyped_measurements" not in report

        with pytest.raises(DataBaseError, match="unknown type"):
            relational.family_table_schemas(MistypedMeasurement)
        with pytest.raises(DataBaseError, match="MistypedMeasurement"):
            MistypedMeasurement(name="broken", technique="simple").save()

    def test_missing_and_named_row_lookup(self, sqlite_backend):
        DataSeries(name="twin", unit_name="V", data=np.array([1.0])).save()
        i_newest = DataSeries(name="twin", unit_name="A", data=np.array([2.0])).save()
        loaded = DataSeries.load("twin")
        assert loaded.id == i_newest
        assert loaded.unit_name == "A"
        with pytest.raises(DataBaseError):
            DataSeries.load("no series has this name")
        with pytest.raises(DataBaseError):
            Measurement.get(999)

    def test_update_with_force(self, sqlite_backend):
        series = DataSeries(name="before", unit_name="V", data=np.array([1.0]))
        i = series.save()
        series.name = "after"
        assert sqlite_backend.save(series, force=True) == i
        assert DataSeries.get(i).name == "after"

    def test_codecs(self, sqlite_backend):
        json_column = ColumnSchema("metadata", dtype=dict)
        value = {"a": np.int64(1), "b": np.array([1.5, 2.5])}
        encoded = sqlite_backend._encode(json_column, value)
        assert sqlite_backend._decode(json_column, encoded) == {
            "a": 1,
            "b": [1.5, 2.5],
        }
        tuple_column = ColumnSchema("range", dtype=tuple)
        encoded = sqlite_backend._encode(tuple_column, (1.0, 2.0))
        assert sqlite_backend._decode(tuple_column, encoded) == (1.0, 2.0)

        array_column = ColumnSchema("data", dtype=np.ndarray)
        data = np.linspace(0, 1, 5)
        encoded = sqlite_backend._encode(array_column, data)
        assert np.array_equal(sqlite_backend._decode(array_column, encoded), data)

        # numpy scalars are stored as their python equivalents:
        assert sqlite_backend._encode(ColumnSchema("x"), np.float64(1.5)) == 1.5

    def test_string_object_array_round_trip(self, sqlite_backend):
        """String arrays from pandas are saved safely without pickle."""
        data = np.array(["1 µA", "10 µA", "100 µA"], dtype=object)
        series = DataSeries(name="Current range", unit_name="A", data=data)

        loaded = DataSeries.get(series.save())

        assert loaded.data.dtype.kind == "U"
        assert loaded.data.tolist() == data.tolist()

    def test_non_string_object_array_is_rejected(self, sqlite_backend):
        """SQLite does not enable pickle for arbitrary Python objects."""
        array_column = ColumnSchema("data", dtype=np.ndarray)
        data = np.array(["text", {"some": "object"}], dtype=object)

        with pytest.raises(ValueError, match="Object arrays cannot be saved"):
            sqlite_backend._encode(array_column, data)

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

        assert (
            sqlite_backend.connection.execute(
                'SELECT COUNT(*) FROM "data_series"'
            ).fetchone()[0]
            == 0
        )
        assert (
            sqlite_backend.connection.execute(
                'SELECT COUNT(*) FROM "measurement"'
            ).fetchone()[0]
            == 0
        )
        assert series._backend is original_backend
        assert series._id is original_id

    def test_reopen_and_lazy_load_from_inactive_backend(self, sqlite_backend):
        series = DataSeries(name="persistent", unit_name="V", data=np.array([1.0, 2.0]))
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

    def test_connections_to_one_file_share_storage_and_rows(
        self, sqlite_backend, tmp_path
    ):
        """Separate connections recognize rows stored in the same file."""
        with SQLiteBackend(db_path=sqlite_backend.db_path) as second_backend:
            with SQLiteBackend(db_path=tmp_path / "other.sqlite") as other_backend:
                assert sqlite_backend is not second_backend
                assert sqlite_backend != second_backend
                assert sqlite_backend.shares_storage_with(second_backend)
                assert second_backend.shares_storage_with(sqlite_backend)
                assert not sqlite_backend.shares_storage_with(other_backend)
                assert sqlite_backend.address == second_backend.address
                assert len({sqlite_backend, second_backend, other_backend}) == 3

                first_id = sqlite_backend.save(
                    DataSeries(name="first", unit_name="V", data=np.array([1.0]))
                )
                second_id = second_backend.save(
                    DataSeries(name="second", unit_name="A", data=np.array([2.0]))
                )
                assert first_id != second_id
                loaded = DataSeries.get(second_id, backend=sqlite_backend)
                assert loaded.name == "second"
                assert loaded.short_identity == (sqlite_backend, second_id)
                loaded_again = DataSeries.get(second_id, backend=second_backend)
                assert loaded.short_identity != loaded_again.short_identity
                assert same_short_identity(
                    loaded.short_identity, loaded_again.short_identity
                )
                assert second_backend.save(loaded) is None
                assert sqlite_backend.connection.execute(
                    'SELECT COUNT(*) FROM "data_series" WHERE "id" = ?', (second_id,)
                ).fetchone() == (1,)
                assert (
                    sqlite_backend.connection.execute(
                        "PRAGMA foreign_key_check"
                    ).fetchall()
                    == []
                )

    def test_in_memory_database_uses_sqlite_convention(self):
        with SQLiteBackend(db_path=":memory:") as first_backend:
            first_series = DataSeries(
                name="first temporary", unit_name="V", data=np.array([1.0])
            )
            first_id = first_backend.save(first_series)
            assert (
                DataSeries.get(first_id, backend=first_backend).name == "first temporary"
            )
            with SQLiteBackend(db_path=":memory:") as second_backend:
                assert not second_backend.contains("data_series", first_id)
                second_series = DataSeries(
                    name="second temporary", unit_name="A", data=np.array([2.0])
                )
                second_id = second_backend.save(second_series)

                assert first_backend != second_backend
                assert not first_backend.shares_storage_with(second_backend)
                assert first_backend.address != second_backend.address
                assert first_id == second_id == 1
                assert first_series.full_identity != second_series.full_identity
                assert not same_short_identity(
                    first_series.short_identity, second_series.short_identity
                )
                assert (
                    DataSeries.get(second_id, backend=second_backend).name
                    == "second temporary"
                )

    def test_chimera_save_reuses_local_series_and_copies_memory_series(
        self, sqlite_backend
    ):
        ec_series = DataSeries(name="EC data", unit_name="V", data=np.array([1.0]))
        ec = Measurement(name="EC", technique="EC", series_list=[ec_series])
        ec.save()
        ec_series_id = ec_series.id

        ms_series = DataSeries(name="MS data", unit_name="A", data=np.array([2.0]))
        ms = Measurement(name="MS", technique="MS", series_list=[ms_series])
        ecms = ec + ms

        assert ec_series.short_identity == (sqlite_backend, ec_series_id)
        assert ms_series.short_identity[0].backend_type == "memory"

        ecms.save()

        assert ec_series.short_identity == (sqlite_backend, ec_series_id)
        assert ms_series.short_identity == (sqlite_backend, ms_series.id)
        assert sqlite_backend.connection.execute(
            'SELECT COUNT(*) FROM "data_series" WHERE "id" = ?', (ec_series_id,)
        ).fetchone() == (1,)
        assert sqlite_backend.connection.execute(
            'SELECT COUNT(*) FROM "data_series"'
        ).fetchone() == (2,)

    def test_unversioned_database_is_rejected_without_changes(self, tmp_path):
        db_path = tmp_path / "unversioned.sqlite"
        with sqlite3.connect(db_path) as connection:
            connection.execute(
                'CREATE TABLE "data_series" ' '("id" INTEGER PRIMARY KEY, "name" TEXT)'
            )
        schema_before = db_path.read_bytes()

        with pytest.raises(DataBaseError, match="no ixdat schema version"):
            SQLiteBackend(db_path=db_path)

        assert db_path.read_bytes() == schema_before

    def test_missing_columns_require_manual_migration(self, tmp_path):
        db_path = tmp_path / "old-layout.sqlite"
        with sqlite3.connect(db_path) as connection:
            create_schema_metadata(connection)
            connection.execute(
                'CREATE TABLE "data_series" ' '("id" INTEGER PRIMARY KEY, "name" TEXT)'
            )

        with SQLiteBackend(db_path=db_path) as backend:
            with pytest.raises(DataBaseError, match="manual database migration"):
                backend.save(DataSeries(name="x", unit_name="V", data=np.array([1.0])))
            columns = {
                row[1]
                for row in backend.connection.execute('PRAGMA table_info("data_series")')
            }
            assert columns == {"id", "name"}

    def test_incompatible_schema_requires_manual_migration(self, tmp_path):
        db_path = tmp_path / "incompatible.sqlite"
        with sqlite3.connect(db_path) as connection:
            create_schema_metadata(connection)
            connection.execute(
                'CREATE TABLE "data_series" ("id" TEXT PRIMARY KEY, "name" TEXT)'
            )

        with SQLiteBackend(db_path=db_path) as backend:
            with pytest.raises(DataBaseError, match="manual database migration"):
                backend.save(DataSeries(name="x", unit_name="V", data=np.array([1.0])))

    @pytest.mark.parametrize(
        ("version", "relation"),
        [(SCHEMA_VERSION - 1, "older"), (SCHEMA_VERSION + 1, "newer")],
    )
    def test_other_schema_versions_are_rejected(self, tmp_path, version, relation):
        db_path = tmp_path / f"schema-{version}.sqlite"
        with sqlite3.connect(db_path) as connection:
            create_schema_metadata(connection, version)

        with pytest.raises(DataBaseError, match=f"{relation} ixdat schema"):
            SQLiteBackend(db_path=db_path)

    def test_get_does_not_create_missing_tables(self, tmp_path):
        with SQLiteBackend(db_path=tmp_path / "empty.sqlite") as backend:
            tables_before = backend._existing_tables()
            with pytest.raises(DataBaseError, match="has no table"):
                DataSeries.get(1, backend=backend)
            assert backend._existing_tables() == tables_before

    def test_get_does_not_recreate_missing_indexes(self, tmp_path):
        db_path = tmp_path / "no-read-writes.sqlite"
        with SQLiteBackend(db_path=db_path) as backend:
            series_id = backend.save(
                DataSeries(name="x", unit_name="V", data=np.array([1.0]))
            )
            backend.connection.execute('DROP INDEX "ixdat_data_series_name_id"')

        with SQLiteBackend(db_path=db_path) as backend:
            indexes_before = backend.connection.execute(
                'PRAGMA index_list("data_series")'
            ).fetchall()
            assert DataSeries.get(series_id, backend=backend).name == "x"
            assert (
                backend.connection.execute('PRAGMA index_list("data_series")').fetchall()
                == indexes_before
            )

    def test_failed_explicit_backend_load_restores_active_backend(
        self, sqlite_backend, tmp_path
    ):
        original_backend = DB.backend
        with SQLiteBackend(db_path=tmp_path / "other.sqlite") as other_backend:
            with pytest.raises(DataBaseError):
                DataSeries.load("missing", backend=other_backend)
            assert DB.backend is original_backend

    def test_lookup_indexes_are_created(self, sqlite_backend):
        sqlite_backend._ensure_tables_for_save(Measurement)
        indexes = {
            row[1]
            for row in sqlite_backend.connection.execute(
                'PRAGMA index_list("measurement")'
            )
        }
        assert "ixdat_measurement_name_id" in indexes
        assert "ixdat_measurement_technique" in indexes


def test_memory_backend_load_by_name():
    """`load` is implemented for the memory backend too, not just the real ones"""
    backend = MemoryBackend()
    backend.save(DataSeries(name="twin", unit_name="V", data=np.array([1.0])))
    newest = DataSeries(name="twin", unit_name="A", data=np.array([2.0]))
    backend.save(newest)

    assert backend.load(DataSeries, "twin") is newest
    with pytest.raises(DataBaseError):
        backend.load(DataSeries, "no series has this name")


def test_load_data_takes_a_backend_and_deprecates_db(sqlite_backend):
    """`load_data(db=...)` used to take a DataBase; it now takes a backend"""
    series = DataSeries(name="s", unit_name="V", data=np.array([1.0, 2.0]))
    loaded = DataSeries.get(series.save())
    assert np.array_equal(loaded.load_data(), [1.0, 2.0])
    assert np.array_equal(loaded.load_data(sqlite_backend), [1.0, 2.0])
    with pytest.deprecated_call():
        assert np.array_equal(loaded.load_data(db=sqlite_backend), [1.0, 2.0])


def test_directory_load_uses_exact_unescaped_name(tmp_path):
    original_backend = DB.backend
    backend = DirBackend(directory=tmp_path, project_name="name_collision")
    DB.set_backend(backend)
    try:
        slash_id = DataSeries(name="a/b", unit_name="V", data=np.array([1.0])).save()
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


def test_directory_backend_stores_local_reference_ids(tmp_path):
    original_backend = DB.backend
    backend = DirBackend(directory=tmp_path, project_name="tuple_references")
    DB.set_backend(backend)
    try:
        series = DataSeries(name="signal", unit_name="V", data=np.array([1.0]))
        measurement = Measurement(
            name="measurement", technique="test", series_list=[series]
        )
        measurement_id = measurement.save()

        stored = backend.get_row_as_dict("measurement", measurement_id)
        assert stored["s_ids"] == [series.id]

        loaded = Measurement.get(measurement_id)
        assert loaded.s_ids == [(backend, series.id)]
    finally:
        DB.set_backend(original_backend)
