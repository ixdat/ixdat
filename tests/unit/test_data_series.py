"""Unit tests for data series"""

from datetime import timezone, timedelta, datetime
from string import ascii_uppercase
from unittest.mock import patch

import pytest
from dateutil.relativedelta import relativedelta
from numpy import arange, array

from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.config import config


@pytest.fixture
def cfg():
    """Fixture that patches out certain items on config.config"""
    with patch.object(config, "timezone", timezone(timedelta(hours=2), "CEST")):
        with patch.object(config, "timestamp_string_format", "native"):
            yield config


class TestTimeSeries:
    """Test the TimeSeries class"""

    def test_str(self, cfg):
        """Test __str__ method"""
        for n, month_letter in enumerate(ascii_uppercase[:12]):
            # The timestamp is equal to 20A18 15:17:14 in UTC+2
            dt = datetime.fromtimestamp(1642511834.9103568) + relativedelta(months=n)
            tseries = TimeSeries(
                "time series", "s", arange(10.9, 1000.0, 100), dt.timestamp()
            )
        assert (
            str(tseries)
            == f"TimeSeries: 'time series', Span: 11..911s @ 22{month_letter}18 15:17:14"
        )


class TestValueSeries:
    """test the ValueSeries class"""

    def test_str(self):
        """Test __str__ method"""
        value_series = ValueSeries(
            "MFC1 flow",
            "ml/min",
            array([1.21111e-3, 4.58888e-6]),
            tseries=TimeSeries("time", "s", array([0.0, 0.1]), tstamp=100.0),
        )
        assert (
            str(value_series) == "ValueSeries: 'MFC1 flow', Span: 4.6e-06..1.2e-03ml/min"
        )
