"""Unit tests for the relational schema derivation and the SQLite backend"""

import numpy as np
import pytest

from ixdat import Measurement, Spectrum
from ixdat.backends import relational
from ixdat.backends.relational import ColumnSchema
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
        value = {"a": 1, "b": [1.5, "two"]}
        encoded = sqlite_backend._encode(json_column, value)
        assert sqlite_backend._decode(json_column, encoded) == value

        array_column = ColumnSchema("data", dtype="NDARRAY")
        data = np.linspace(0, 1, 5)
        encoded = sqlite_backend._encode(array_column, data)
        assert np.array_equal(sqlite_backend._decode(array_column, encoded), data)

        # numpy scalars are stored as their python equivalents:
        assert sqlite_backend._encode(ColumnSchema("x"), np.float64(1.5)) == 1.5

    def test_dynamic_column_rejects_containers(self, sqlite_backend):
        """Containers require a logical type registered in relational.COLUMN_TYPES"""
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
