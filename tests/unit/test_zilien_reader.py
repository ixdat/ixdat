from pathlib import Path

import pytest
from pytest import fixture

from ixdat import Measurement
from ixdat.readers.zilien import ZilienTSVReader, to_mass


DATA_DIR = Path(__file__).parent.parent.parent / "submodules" / "ixdat-large-test-files"


@fixture(scope="class")
def reader():
    return ZilienTSVReader()


class TestZilienTSVReader:
    """Tests for Zilien TSV Reader."""

    def _test_first(self):
        m = Measurement.read(
            DATA_DIR
            / "zilien_version_2/2022-09-19 11_27_59 test-1/2022-09-19 11_27_59 test-21.tsv",  # noqa: E501
            reader="zilien",
        )
        print(m.series_list)


class TestZilienTSVReaderUtils:
    """Tests for Zilien TSV Reader utility methods."""

    @pytest.mark.parametrize(
        "data, expected",
        (
            ("not-a-match", None),
            ("C0M2", "2"),
            ("C2M15", "15"),
            ("C10M44", "44"),
        )
    )
    def test_to_mass(self, data, expected):
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
    def test_get_series_splits(self, reader, data, expected):
        """Test series splits."""
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
    def test_get_biologic_splits(self, reader, data, expected):
        """Test biologic splits."""
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
    def test_form_names_and_unit(self, reader, series_header, column_header, expected):
        """Test forming names and units."""
        assert reader._form_names_and_unit(series_header, column_header) == expected
