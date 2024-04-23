from pathlib import Path
from unittest.mock import patch

import pytest

from ixdat import Measurement
from ixdat.exceptions import SeriesNotFoundError
from ixdat.techniques.ec_ms import (
    ECMSMeasurement,
    MSMeasurement,
)


# fmt: off
DIR_FAILING = Path(__file__).parent.parent.parent / "submodules/ixdat-large-test-files/failing_datasets"  # noqa: E501

OK_ECMS_DATASET = DIR_FAILING / "2024-01-17 15_21_57 ec-included/2024-01-17 15_21_57 ec-included.tsv"  # noqa: E501
OK_MS_DATASET = DIR_FAILING / "2024-01-17 15_06_41 no-ec-nonincluded/2024-01-17 15_06_41 no-ec-nonincluded.tsv"  # noqa: E501
MISSING_POT_SERIES_BUG = DIR_FAILING / "2024-01-17 15_08_29 no-ec-included/2024-01-17 15_08_29 no-ec-included.tsv"  # noqa: E501
MISSING_ECLAB_MPTS_SERIES_BUG = DIR_FAILING / "2024-01-17 15_13_45 ec-included-mpt-missing/2024-01-17 15_13_45 ec-included-mpt-missing.tsv"  # noqa: E501
# fmt: on


@pytest.mark.external
class TestRegressions:
    """Test patched way and the old way of dataset parsing."""

    def test_parse_ecms_measurement(self):
        """Test parsing EC-MS dataset without specifying a technique."""
        m = Measurement.read(OK_ECMS_DATASET)
        assert isinstance(m, ECMSMeasurement)

    def test_parse_ms_measurement(self):
        """Test parsing MS dataset without specifying a technique."""
        m = Measurement.read(OK_MS_DATASET)
        assert isinstance(m, MSMeasurement)

    def test_parse_ms_measurement_read_as_ecms_regression(self):
        """Test parsing MS dataset the old way.

        The dataset used to be created with ECMSMeasurement class by default.
        That was failing, because an essential series 't' was missing.

        """
        with pytest.raises(SeriesNotFoundError):
            Measurement.read(OK_MS_DATASET, technique="EC-MS")

    def test_parse_with_missing_zilien_pot_data(self):
        """Test patched parsing a buggy MS dataset with empty "pot" part."""
        m = Measurement.read(MISSING_POT_SERIES_BUG)
        assert isinstance(m, MSMeasurement)

    def test_parse_with_missing_zilien_pot_data_regression(self):
        """Test parsing a buggy MS dataset with empty "pot" part the old way.

        The dataset used to be created with ECMSMeasurement class by default.
        That was failing, because of the internal bug in Zilien, which caused
        including empty "pot" series.

        """
        with pytest.raises(KeyError):
            Measurement.read(MISSING_POT_SERIES_BUG, technique="EC-MS")

    def test_parse_with_missing_biologic_data(self):
        """Test patched parsing a buggy EC-MS dataset with missing EC-lab part."""

        # silence the warning message for tests
        with patch("sys.stderr.write"):
            m = Measurement.read(MISSING_ECLAB_MPTS_SERIES_BUG)
            assert isinstance(m, MSMeasurement)

    def test_parse_with_missing_biologic_data_regression(self):
        """Test parsing a buggy EC-MS dataset with missing EC-lab part the old way.

        The dataset used to be created with ECMSMeasurement class by default.
        That was failing, because all the "EC-lab" series were missing, when no
        .mpt files were found after Zilien measurement ended.

        """
        # silence the warning message for tests
        with pytest.raises(ValueError):
            Measurement.read(MISSING_ECLAB_MPTS_SERIES_BUG, technique="EC-MS")
