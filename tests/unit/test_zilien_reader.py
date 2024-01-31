from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.exceptions import TechniqueError
from ixdat.readers.zilien import (
    ZilienTSVReader,
    determine_class,
    to_mass,
    to_snake_case,
    parse_metadata_line,
)
from ixdat.techniques import ECMeasurement, MSMeasurement, ECMSMeasurement


DATA_DIR = Path(__file__).parent.parent.parent / "submodules" / "ixdat-large-test-files"


def patch_biologic_dataset_part(column_headers_cut, data_cut):
    """A patch function to replace the original one inside the _form_series method."""
    return [
        (column_header, data_cut[0, i])
        for i, column_header in enumerate(column_headers_cut)
    ]


def patch_zilien_dataset_part(series_header, column_headers_cut, data_cut):
    """A patch function to replace the original one inside the _form_series method."""
    mass = to_mass(series_header)
    aliases = {f"M{mass}": [f"M{mass} [A]"]} if mass is not None else {}

    return [
        (f"{series_header}_{column_header}", data_cut[0, i])
        for i, column_header in enumerate(column_headers_cut)
    ], aliases


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

        series = reader._create_series_objects(column_headers, names_and_units, data)

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
        with pytest.raises(ValueError):
            reader._create_series_objects(column_headers, names_and_units, data)

    @pytest.mark.parametrize(
        "series_header, column_headers, expected_aliases",
        (
            ("MFC1 value", ["Time [s]", "Flow [ml/min]"], {}),
            ("C4M26", ["Time [s]", "M28-N2/CO [A]"], {"M26": ["M26 [A]"]}),
        ),
    )
    def test_zilien_dataset_part(self, series_header, column_headers, expected_aliases):
        """Test forming of Zilien part of a dataset.

        The headers are parametrized, because some headers do produce an alias
        (C4M26 -> M26) and some do not. The data stays the same, because the
        values in the columns don't have an effect on the parsing.

        """
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
        "data, lines_count, expected_series_length",
        # fmt: off
        (
            (
                np.array([
                    [0.1, 1, 1, "NaN", 3],
                    [0.2, 1, 1, "NaN", 4],
                    [0.3, 1, 1, "NaN", 5],
                    ["NaN", "NaN", "NaN", "NaN", "NaN"],
                ], dtype=float),
                3, (3, 3),
            ),
            (
                np.array([
                    [0.1, 1, 1, "NaN", 3],
                    [0.2, 1, 2, "NaN", 4],
                    [0.3, 1, 2, "NaN", 5],
                    ["NaN", "NaN", "NaN", "NaN", "NaN"],
                ], dtype=float),
                3, (1, 1, 2, 2),
            ),
            (
                np.array([
                    [0.1, 1, 1, "NaN", 3],
                    [0.2, 1, 2, "NaN", 4],
                    [0.3, 1, 2, "NaN", 5],
                    [0.4, 2, 1, "NaN", 3],
                    [0.5, 2, 1, "NaN", 4],
                    [0.6, 2, 2, "NaN", 5],
                    ["NaN", "NaN", "NaN", "NaN", "NaN"],
                ], dtype=float),
                6, (1, 1, 2, 2, 2, 2, 1, 1),
            ),
            (
                np.array([
                    [0.1, 1, 1, "NaN", 3],
                    [0.2, 1, 1, "NaN", 4],
                    [0.3, 2, 1, "NaN", 5],
                    [0.4, 2, 1, "NaN", 3],
                    ["NaN", "NaN", "NaN", "NaN", "NaN"],
                ], dtype=float),
                4, (2, 2, 2, 2),
            ),
        ),
        # fmt: on
    )
    def test_biologic_dataset_part(self, data, lines_count, expected_series_length):
        """Test forming of Biologic part of a dataset.

        The data change, because the "experiment_number" and the "technique_number"
        columns are used to split the data chuck into more chunks. The headers
        stay the same, because nothing else is produces from them.

        """
        reader = ZilienTSVReader()
        reader._timestamp = 123456
        reader._metadata = {
            "MFC1 value_MFC1 value_count": 2,
            "C4M26_C4M26_count": 2,
            "EC-lab_EC-lab_count": lines_count,
        }

        column_headers = [
            "time/s",
            "experiment_number",
            "technique_number",
            "Ewe/V",
            "I/mA",
        ]

        series = reader._biologic_dataset_part(column_headers, data)

        assert len(series) == len(expected_series_length)

        for i, count in enumerate(expected_series_length):
            assert series[i].data.size == count

    @pytest.mark.parametrize(
        "metadata, series_headers, column_headers, expected",
        # fmt: off
        (
            (  # good case with EC-lab data
                {"pot_pot_count": 3},
                ["pot", "", "", "", "Iongauge value", "", "EC-lab", "", ""],
                [
                    "Time [s]", "Voltage [V]", "Current [mA]", "Cycle [n]", "Time [s]", "Pressure [mbar]",  # noqa: E501
                    "time/s", "experiment_number", "technique_number", "Ewe/V"
                ],
                "EC-MS"
            ),
            (  # bug case without EC-lab data
                {"pot_pot_count": 3},
                ["pot", "", "", "", "Iongauge value", "", "EC-lab"],
                [
                    "Time [s]", "Voltage [V]", "Current [mA]", "Cycle [n]", "Time [s]", "Pressure [mbar]",  # noqa: E501
                ],
                "MS"
            ),
            (  # bug case without pot data
                {},
                ["pot", "", "", "", "Iongauge value", ""],
                [
                    "Time [s]", "Voltage [V]", "Current [mA]", "Cycle [n]", "Time [s]", "Pressure [mbar]",  # noqa: E501
                ],
                "MS"
            ),
            (  # good case without EC measurement run
                {},
                ["Iongauge value", ""],
                ["Time [s]", "Pressure [mbar]"],
                "MS"
            ),
            (  # random and unexpected case
                {"pot_pot_count": 3},
                ["Iongauge value", ""],
                ["Time [s]", "Pressure [mbar]"],
                "MS"
            ),
            (  # random and unexpected case
                {},
                ["Iongauge value", "", "EC-lab", ""],
                ["Time [s]", "Pressure [mbar]", "time/s", "Ewe/V"],
                "MS"
            ),
        )
        # fmt: on
    )
    @patch("sys.stderr.write")
    def test_get_technique(self, _, metadata, series_headers, column_headers, expected):
        reader = ZilienTSVReader()
        reader._metadata = metadata
        reader._series_headers = series_headers
        reader._column_headers = column_headers

        assert reader._get_technique() == expected

    @pytest.mark.parametrize(
        "klass, expected_aliases, expected_series",
        # fmt: off
        (
            (
                ECMSMeasurement, {"M40": ["M40 [A]"]}, [
                    ("pot_Time [s]", 0.1),
                    ("pot_Voltage [V]", 11.0),
                    ("pot_Current [mA]", 12.0),
                    ("pot_Cycle [n]", 13.0),
                    ("MFC1 value_Time [s]", 0.2),
                    ("MFC1 value_Flow [ml/min]", 21.0),
                    ("C8M40_Time [s]", 0.3),
                    ("C8M40_M40-Ar [A]", 31.0),
                    ("time/s", 0.4),
                    ("technique_number", 41.0),
                    ("Ewe/V", 42.0),
                    ("I/mA", 43.0),
                ],
            ),
            (
                ECMeasurement, {}, [
                    ("pot_Time [s]", 0.1),
                    ("pot_Voltage [V]", 11.0),
                    ("pot_Current [mA]", 12.0),
                    ("pot_Cycle [n]", 13.0),
                    ("MFC1 value_Time [s]", 0.2),
                    ("MFC1 value_Flow [ml/min]", 21.0),
                    ("time/s", 0.4),
                    ("technique_number", 41.0),
                    ("Ewe/V", 42.0),
                    ("I/mA", 43.0),
                ],
            ),
            (
                MSMeasurement, {"M40": ["M40 [A]"]}, [
                    ("MFC1 value_Time [s]", 0.2),
                    ("MFC1 value_Flow [ml/min]", 21.0),
                    ("C8M40_Time [s]", 0.3),
                    ("C8M40_M40-Ar [A]", 31.0),
                ],
            ),
        ),
        # fmt: on
    )
    def test_form_series(self, klass, expected_aliases, expected_series):
        """Test forming series."""

        reader = ZilienTSVReader()
        reader._cls = klass
        reader._timestamp = 123456
        # fmt: off
        reader._series_headers = [
            "pot", "", "", "", "MFC1 value", "", "C8M40", "", "EC-lab", "", "", "",
        ]
        reader._column_headers = [
            "Time [s]", "Voltage [V]", "Current [mA]", "Cycle [n]",
            "Time [s]", "Flow [ml/min]",
            "Time [s]", "M40-Ar [A]",
            "time/s", "technique_number", "Ewe/V", "I/mA",
        ]
        reader._data = np.array([
            [0.1, 11, 12, 13, 0.2, 21, 0.3, 31, 0.4, 41, 42, 43],
        ])
        # fmt: on

        with patch.multiple(
            reader,
            _zilien_dataset_part=patch_zilien_dataset_part,
            _biologic_dataset_part=patch_biologic_dataset_part,
        ):
            series, aliases = reader._form_series()

        assert series == expected_series
        assert dict(aliases) == expected_aliases


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
                'C4M26_count\tNumber of points in the "C4M26" data series\tC4M26\tint\t2\n',  # noqa: E501
                "C4M26_C4M26_count",
                2,
            ),
            (
                "start_time_unix\tStart time as UNIX time stamp\t\tdouble\t1663579682.70404\n",  # noqa: E501
                "start_time_unix",
                1663579682.70404,
            ),
            ("subscribe\tSend events\tC15M56\tbool\ttrue\n", "C15M56_subscribe", True),
        ),
    )
    def test_parse_metadata_line(self, line, expected_full_name, expected_value):
        """Test parsing metadata lines."""
        full_name, value = parse_metadata_line(line)

        assert full_name == expected_full_name
        if isinstance(value, float):
            assert value == pytest.approx(expected_value)
        else:
            assert value == expected_value

    def test_parse_metadata_line_error(self):
        """Test raising an error."""
        with pytest.raises(TypeError):
            parse_metadata_line("\t\t\t\t\n")

    @pytest.mark.parametrize(
        "data, expected",
        (
            ("", ""),
            ("a", "a"),
            ("A", "a"),
            (" ", "_"),
            ("_", "_"),
            ("A String to-convert", "a_string_to-convert"),
            ("TRAILING  SPACE  ", "trailing__space__"),
        ),
    )
    def test_to_snake_case(self, data, expected):
        """Test converting a string to snakecase."""
        assert to_snake_case(data) == expected

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
            ("EC-MS", ECMSMeasurement),
            ("EC", ECMeasurement),
            ("MS", MSMeasurement),
        ),
    )
    def test_determine_class(self, data, expected):
        """Test determining a class type."""
        assert determine_class(data) == expected

    @pytest.mark.parametrize(
        "data",
        ("ec-ms", "ec", "ms", "read the exception message finally"),
    )
    def test_determine_class_error(self, data):
        """Test an exception while determining a class type."""
        with pytest.raises(TechniqueError):
            determine_class(data)

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
