from pathlib import Path

import numpy as np
import pytest
from pytest import fixture

from ixdat.readers.zilien import ZilienTSVReader, to_mass, parse_metadata_line
from ixdat.data_series import TimeSeries, ValueSeries


DATA_DIR = Path(__file__).parent.parent.parent / "submodules" / "ixdat-large-test-files"


class TestZilienTSVReader:
    """Tests for Zilien TSV Reader."""

    def test_create_series_objects(self):
        """Test creation of Ixdat series objects."""

        reader = ZilienTSVReader()
        reader._timestamp = 123456

        column_headers = ["time/s", "experiment_number", "I/mA", "Ewe/V"]
        names_and_units = [
            ("Biologic time/s", "s", None),
            ("experiment_number", "", None),
            ("I/mA", "mA", None),
            ("Ewe/V", "V", None),
        ]
        data = np.array(
            [
                [0.1, 1, 3, "NaN"],
                [0.2, 1, 4, "NaN"],
            ],
            dtype=float,
        )

        series = reader._create_series_objects(
            column_headers, names_and_units, data
        )

        assert len(series) == 2
        assert isinstance(series[0], TimeSeries) is True
        assert isinstance(series[1], ValueSeries) is True
        assert series[0] is series[1].tseries
        assert series[0].tstamp == 123456
        assert series[0].unit.name == "s"
        assert series[1].unit.name == "mA"
        assert series[0].name == "Biologic time/s"
        assert series[1].name == "I/mA"
        assert list(series[0].data) == [0.1, 0.2]
        assert list(series[1].data) == [3, 4]

        # when time is not first, the dataset will not be parsed
        column_headers = ["experiment_number", "I/mA", "Ewe/V", "time/s"]
        with pytest.raises(AssertionError):
            reader._create_series_objects(
                column_headers, names_and_units, data
            )

    @pytest.mark.parametrize(
        "series_header, column_headers, expected_aliases",
        (
            ("MFC1 value", ["Time [s]", "Flow [ml/min]"], {}),
            ("C4M26", ["Time [s]", "M28-N2/CO [A]"], {"M26": ["M26 [A]"]}),
        ),
    )
    def test_zilien_dataset_part(
        self, series_header, column_headers, expected_aliases
    ):
        """Test forming of Zilien part of a dataset."""

        reader = ZilienTSVReader()
        reader._timestamp = 123456
        reader._metadata = {
            "MFC1 value_MFC1 value_count": 2,
            "C4M26_C4M26_count": 2,
        }

        data = np.array(
            [
                [0.1, 1],
                [0.2, 1],
                ["NaN", "NaN"],
            ],
            dtype=float,
        )

        series, aliases = reader._zilien_dataset_part(
            series_header, column_headers, data
        )

        assert len(series) == 2
        assert series[0].data.size == 2
        assert series[1].data.size == 2
        assert dict(aliases) == expected_aliases

    @pytest.mark.parametrize(
        "data, lines_count, expected",
        (
            (
                np.array(
                    [
                        [0.1, 1, "NaN", 3],
                        [0.2, 1, "NaN", 4],
                        [0.3, 1, "NaN", 5],
                        ["NaN", "NaN", "NaN", "NaN"],
                    ],
                    dtype=float,
                ),
                3,
                (3, 3),
            ),
            (
                np.array(
                    [
                        [0.1, 1, "NaN", 3],
                        [0.2, 2, "NaN", 4],
                        [0.3, 2, "NaN", 5],
                        ["NaN", "NaN", "NaN", "NaN"],
                    ],
                    dtype=float,
                ),
                3,
                (1, 1, 2, 2),
            ),
            (
                np.array(
                    [
                        [0.1, 1, "NaN", 3],
                        [0.2, 2, "NaN", 4],
                        [0.3, 2, "NaN", 5],
                        [0.4, 1, "NaN", 3],
                        [0.5, 1, "NaN", 4],
                        [0.6, 2, "NaN", 5],
                        ["NaN", "NaN", "NaN", "NaN"],
                    ],
                    dtype=float,
                ),
                6,
                (1, 1, 2, 2, 2, 2, 1, 1),
            ),
        ),
    )
    def test_biologic_dataset_part(
        self, data, lines_count, expected
    ):
        """Test forming of Biologic part of a dataset."""

        reader = ZilienTSVReader()
        reader._timestamp = 123456
        reader._metadata = {
            "MFC1 value_MFC1 value_count": 2,
            "C4M26_C4M26_count": 2,
            "EC-lab_EC-lab_count": lines_count,
        }

        column_headers = ["time/s", "technique_number", "Ewe/V", "I/mA"]

        series = reader._biologic_dataset_part(column_headers, data)

        assert len(series) == len(expected)

        for index, count in enumerate(expected):
            assert series[index].data.size == count


class TestZilienTSVReaderUtils:
    """Tests for Zilien TSV Reader utility methods."""

    @pytest.mark.parametrize(
        "line, expected_full_name, expected_value",
        (
            (
                "label_suffix\tLabel_suffix\tC15M56\tstring\tHex\n",
                "C15M56_label_suffix",
                "Hex",
            ),
            (
                'C4M26_count\tNumber of points in the "C4M26" data series\tC4M26\tint\t2\n',
                "C4M26_C4M26_count",
                2,
            ),
            (
                "start_time_unix\tStart time as UNIX time stamp\t\tdouble\t1663579682.70404\n",
                "start_time_unix",
                1663579682.70404,
            ),
            ("subscribe\tSend events\tC15M56\tbool\ttrue\n", "C15M56_subscribe", True),
        ),
    )
    def test_parse_metadata_line(self, line, expected_full_name, expected_value):
        """Test parsing metadata lines."""
        assert parse_metadata_line(line) == (expected_full_name, expected_value)

    def test_parse_metadata_line_error(self):
        """Test raising an error."""
        with pytest.raises(TypeError):
            parse_metadata_line("\t\t\t\t\n")

    @pytest.mark.parametrize(
        "data, expected",
        (
            ("not-a-match", None),
            ("C0M2", "2"),
            ("C2M15", "15"),
            ("C10M44", "44"),
        ),
    )
    def test_to_mass(self, data, expected):
        """Test conversion to mass."""
        assert to_mass(data) == expected

    @pytest.mark.parametrize(
        "data, expected",
        (
            ([], ([], [])),
            (["pot", "", "", ""], ([(0, None)], ["pot"])),
            (
                ["pot", "", "", "", "MFC1 setpoint", ""],
                (
                    [(0, 4), (4, None)],
                    ["pot", "MFC1 setpoint"],
                ),
            ),
            (
                ["pot", "", "", "", "MFC1 setpoint", "", "EC-lab", "", ""],
                (
                    [(0, 4), (4, 6), (6, None)],
                    ["pot", "MFC1 setpoint", "EC-lab"],
                ),
            ),
        ),
    )
    def test_get_series_splits(self, data, expected):
        """Test series splits."""
        reader = ZilienTSVReader()
        assert reader._get_series_splits(data) == expected

    @pytest.mark.parametrize(
        "data, expected",
        (
            ([], []),
            ([1, 1], [(0, 2)]),
            ([1, 1, 2, 2, 2], [(0, 2), (2, 5)]),
            (
                [1, 1, 2, 2, 2, 1, 1, 1, 2, 3, 3],
                [
                    (0, 2),
                    (2, 5),
                    (5, 8),
                    (8, 9),
                    (9, 11),
                ],
            ),
        ),
    )
    def test_get_biologic_splits(self, data, expected):
        """Test biologic splits."""
        reader = ZilienTSVReader()
        assert reader._get_biologic_splits(data) == expected

    @pytest.mark.parametrize(
        "series_header, column_header, expected",
        (
            ("", "", ("", "", None)),
            ("pot", "Time [s]", ("Potential time [s]", "s", None)),
            ("MFC1 value", "Time [s]", ("MFC1 value time [s]", "s", None)),
            ("EC-lab", "time/s", ("Biologic time/s", "s", None)),
            ("pot", "Voltage [V]", ("Voltage [V]", "V", None)),
            (
                "MFC1 setpoint",
                "Flow [ml/min]",
                ("MFC1 setpoint [ml/min]", "ml/min", None),
            ),
            ("Pirani1 value", "Pressure [mbar]", ("Pirani1 value [mbar]", "mbar", None)),
            ("C1M4", "M4-He [A]", ("M4 [A]", "A", "M4")),
            ("EC-lab", "experiment_number", ("experiment_number", "", None)),
            ("EC-lab", "I/mA", ("I/mA", "mA", None)),
        ),
    )
    def test_form_names_and_unit(self, series_header, column_header, expected):
        """Test forming names and units."""
        reader = ZilienTSVReader()
        assert reader._form_names_and_unit(series_header, column_header) == expected
