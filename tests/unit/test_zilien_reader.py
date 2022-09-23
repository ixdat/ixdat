from pathlib import Path

from pytest import fixture

from ixdat import Measurement
from ixdat.readers.zilien import ZilienTSVReader


DATA_DIR = Path(__file__).parent.parent.parent / "submodules" / "ixdat-large-test-files"


@fixture(scope="class")
def reader():
    return ZilienTSVReader()


class TestZilienTSVReader:
    """Tests for Zilien TSV Reader."""

    def _test_first(self):
        m = Measurement.read(
            DATA_DIR
            / "zilien_version_2/2022-09-19 11_27_59 test-1/2022-09-19 11_27_59 test-21.tsv",
            reader="zilien",
        )
        print(m.series_list)


class TestZilienTSVReaderUtils:
    """Tests for Zilien TSV Reader utility methods."""

    def test_get_series_splits_empty(self, reader):
        """Test series splits for an empty list."""
        data = []
        assert reader._get_series_splits(data) == ([], [])

    def test_get_series_splits_one(self, reader):
        """Test series splits for a list with one group."""
        data = ["pot", "", "", ""]
        assert reader._get_series_splits(data) == ([(0, None)], ["pot"])

    def test_get_series_splits_two(self, reader):
        """Test series splits for a list with two groups."""
        data = ["pot", "", "", "", "MFC1 setpoint", ""]
        assert reader._get_series_splits(data) == (
            [(0, 4), (4, None)],
            ["pot", "MFC1 setpoint"],
        )

    def test_get_series_splits_multi(self, reader):
        """Test series splits for a list with multiple groups."""
        # fmt: off
        data = ["pot", "", "", "", "MFC1 setpoint", "", "MFC1 value", "", "EC-lab", "", ""]
        # fmt: on
        assert reader._get_series_splits(data) == (
            [(0, 4), (4, 6), (6, 8), (8, None)],
            ["pot", "MFC1 setpoint", "MFC1 value", "EC-lab"],
        )

    def test_get_biologic_splits_empty(self, reader):
        """Test biologic splits for an empty list."""
        data = []
        assert reader._get_biologic_splits(data) == []

    def test_get_biologic_splits_one(self, reader):
        """Test biologic splits for a list with one group."""
        data = [1, 1]
        assert reader._get_biologic_splits(data) == [(0, 2)]

    def test_get_biologic_splits_two(self, reader):
        """Test biologic splits for a list with two groups."""
        data = [1, 1, 2, 2, 2]
        assert reader._get_biologic_splits(data) == [(0, 2), (2, 5)]

    def test_get_biologic_splits_multi(self, reader):
        """Test biologic splits for a list with multiple groups."""
        data = [1, 1, 2, 2, 2, 1, 1, 1, 2, 3, 3]
        assert reader._get_biologic_splits(data) == [
            (0, 2),
            (2, 5),
            (5, 8),
            (8, 9),
            (9, 11),
        ]
