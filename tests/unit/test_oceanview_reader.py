"""Tests for OceanView optical spectrum-series integration."""

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pytest

from ixdat import Spectrum
from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.spectra import Spectrum as SingleSpectrum
from ixdat.techniques import ECMeasurement, TECHNIQUE_CLASSES
from ixdat.techniques.spectroelectrochemistry import (
    ECOpticalMeasurement,
    OpticalSpectrumSeries,
)


FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "test_data"
    / "oceanview_sec"
    / "mini_oceanview__0__15-02-35-123.txt"
)
DATED_FIXTURE = FIXTURE.parent / "spectrum_t_demo_QEP100221__0__18-51-43-086.txt"


def make_ec_measurement():
    """Return a minimal EC measurement for combination tests."""
    tseries = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0, 2.0, 3.0]),
        tstamp=1756386155.123,
    )
    potential = ValueSeries(
        name="raw_potential",
        unit_name="V",
        data=np.array([0.1, 0.2, 0.3, 0.4]),
        tseries=tseries,
    )
    current = ValueSeries(
        name="raw_current",
        unit_name="mA",
        data=np.array([1.0, 1.1, 1.2, 1.3]),
        tseries=tseries,
    )
    return ECMeasurement(
        name="test EC",
        technique="EC",
        series_list=[tseries, potential, current],
    )


def test_oceanview_reader_returns_optical_spectrum_series():
    optical = Spectrum.read(FIXTURE, reader="oceanview")

    assert isinstance(optical, OpticalSpectrumSeries)
    assert optical.technique == "Optical"
    assert optical.field.data.shape == (40, 11)
    np.testing.assert_allclose(
        optical.field.axes_series[1].data,
        [400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0],
    )
    np.testing.assert_allclose(optical.t, np.arange(40.0))


def test_optical_technique_mapping_keeps_spectrum_series_indexing():
    optical = Spectrum.read(FIXTURE, reader="oceanview")

    assert TECHNIQUE_CLASSES["Optical"] is OpticalSpectrumSeries
    spectrum = optical[0]
    assert isinstance(spectrum, SingleSpectrum)
    np.testing.assert_allclose(
        spectrum.x,
        [400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0],
    )
    np.testing.assert_allclose(
        spectrum.y,
        [220.0, 240.0, 270.0, 310.0, 360.0, 420.0, 360.0, 310.0, 270.0, 240.0, 220.0],
    )


def test_ec_plus_optical_spectrum_series_returns_ec_optical_measurement():
    ec = make_ec_measurement()
    optical = Spectrum.read(FIXTURE, reader="oceanview")

    ec_optical = ec + optical

    assert isinstance(ec_optical, ECOpticalMeasurement)
    assert ec_optical.technique == "EC-Optical"
    assert ec_optical.spectrum_series is optical


@pytest.mark.parametrize("tstamp_source", ["header", "spectrum"])
def test_series_and_field_timestamps_agree(tstamp_source):
    """ixdat stores the start time twice; the two must not drift apart.

    If they disagree, ixdat re-anchors the field's time axis and shifts the
    relative times to compensate, which silently moves the data.
    """
    optical = Spectrum.read(
        DATED_FIXTURE, reader="oceanview", tstamp_source=tstamp_source
    )

    assert optical.tstamp == optical.field.axes_series[0].tstamp
    assert optical.field.axes_series[0].data[0] == 0.0


def test_placeholder_date_is_rejected_for_the_spectrum_source():
    """1970-01-01 is what OceanView writes when it has no real date."""
    with pytest.raises(ValueError, match="1970-01-01 placeholder"):
        Spectrum.read(FIXTURE, reader="oceanview", tstamp_source="spectrum")


def unix_time(year, month, day, hour, minute, second, microsecond, offset_hours):
    """Return the unix timestamp of a wall-clock time in a fixed-offset zone."""
    return datetime(
        year,
        month,
        day,
        hour,
        minute,
        second,
        microsecond,
        tzinfo=timezone(timedelta(hours=offset_hours)),
    ).timestamp()


def write_oceanview_file(
    directory, name, stamps, header="Wed Aug 05 12:00:00 CEST 2026"
):
    """Write a minimal OceanView export with the given data-row stamps."""
    path = directory / name
    header_line = f"Date: {header}\n" if header else ""
    rows = "\n".join(f"{stamp}\t1.0\t2.0\t3.0" for stamp in stamps)
    path.write_text(
        f"Data from {name} Node\n"
        f"{header_line}"
        "User: ixdat\n"
        ">>>>>Begin Spectral Data<<<<<\n"
        "400.0\t450.0\t500.0\n"
        f"{rows}\n"
    )
    return path


# --- the numeric GMT offset -----------------------------------------------


@pytest.mark.parametrize(
    "header_offset, offset_hours",
    [("GMT+02:00", 2), ("GMT-05:30", -5.5)],
)
def test_numeric_gmt_offset_gives_the_exact_utc_time(
    tmp_path, header_offset, offset_hours
):
    """OceanView is a Java application: 'GMT+02:00' means two hours ahead."""
    path = write_oceanview_file(
        tmp_path,
        "gmt_offset.txt",
        ["2026-08-05 12:00:00.000000"],
        header=f"Wed Aug 05 12:00:00 {header_offset} 2026",
    )

    optical = Spectrum.read(path, reader="oceanview")

    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 0, 0, offset_hours)


# --- assume_timezone -------------------------------------------------------


def test_assume_timezone_used_when_the_header_states_none(tmp_path):
    path = write_oceanview_file(
        tmp_path,
        "no_header_timezone.txt",
        ["2026-08-05 12:00:00.000000"],
        header="Wed Aug 05 12:00:00 2026",
    )

    optical = Spectrum.read(path, reader="oceanview", assume_timezone="+02:00")

    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 0, 0, 2)


def test_assume_timezone_used_when_the_header_timezone_is_unreadable(tmp_path):
    """'CST' is US Central and China Standard, so the reader will not guess."""
    path = write_oceanview_file(
        tmp_path,
        "unreadable_timezone.txt",
        ["2026-08-05 12:00:00.000000"],
        header="Wed Aug 05 12:00:00 CST 2026",
    )

    optical = Spectrum.read(path, reader="oceanview", assume_timezone="+08:00")

    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 0, 0, 8)


def test_error_when_neither_the_file_nor_the_caller_gives_a_timezone(tmp_path):
    path = write_oceanview_file(
        tmp_path,
        "no_timezone_at_all.txt",
        ["2026-08-05 12:00:00.000000"],
        header="Wed Aug 05 12:00:00 2026",
    )

    with pytest.raises(ValueError, match="assume_timezone"):
        Spectrum.read(path, reader="oceanview")


def test_a_timezone_stated_by_the_file_wins(tmp_path):
    """assume_timezone fills a gap; it does not override the file."""
    path = write_oceanview_file(
        tmp_path, "stated_timezone.txt", ["2026-08-05 12:00:00.000000"]
    )

    optical = Spectrum.read(path, reader="oceanview", assume_timezone="+09:00")

    # The header says CEST, so +02:00 applies rather than the +09:00 passed in.
    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 0, 0, 2)


# --- row timestamp formats -------------------------------------------------


def test_rows_with_a_full_date_and_time_are_read(tmp_path):
    path = write_oceanview_file(
        tmp_path,
        "full_datetime.txt",
        ["2026-08-05 12:00:00.000000", "2026-08-05 12:00:02.000000"],
    )

    optical = Spectrum.read(path, reader="oceanview")

    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 0, 0, 2)
    np.testing.assert_allclose(optical.t, [0.0, 2.0])


def test_spectrum_mode_anchors_to_the_first_row(tmp_path):
    """The rows run two seconds behind the header in this file."""
    path = write_oceanview_file(
        tmp_path,
        "spectrum_anchor.txt",
        ["2026-08-05 12:00:02.000000", "2026-08-05 12:00:04.000000"],
    )

    optical = Spectrum.read(path, reader="oceanview", tstamp_source="spectrum")

    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 2, 0, 2)
    np.testing.assert_allclose(optical.t, [0.0, 2.0])


@pytest.mark.parametrize("tstamp_source", ["header", "spectrum"])
@pytest.mark.parametrize(
    "stamps",
    [
        pytest.param(["12:00:00.000000", "12:00:02.000000"], id="clock_time"),
        pytest.param(["12-00-00-000", "12-00-02-000"], id="hyphen_clock_time"),
        pytest.param(["0.0", "2.0"], id="elapsed_seconds"),
    ],
)
def test_rows_without_a_date_are_rejected(tmp_path, stamps, tstamp_source):
    """This reader requires a full date on every row, in both modes.

    A clock time alone cannot anchor an absolute timestamp, and would run
    backwards across midnight.
    """
    path = write_oceanview_file(tmp_path, "no_date.txt", stamps)

    with pytest.raises(ValueError, match="not a full date and time"):
        Spectrum.read(path, reader="oceanview", tstamp_source=tstamp_source)


@pytest.mark.parametrize("tstamp_source", ["header", "spectrum"])
def test_unreadable_row_timestamps_are_rejected(tmp_path, tstamp_source):
    path = write_oceanview_file(
        tmp_path, "unreadable_rows.txt", ["not-a-time", "also-not-a-time"]
    )

    with pytest.raises(ValueError, match="not a full date and time"):
        Spectrum.read(path, reader="oceanview", tstamp_source=tstamp_source)


def test_relative_times_are_correct_across_midnight(tmp_path):
    """Clock-time arithmetic would make the second row jump back a day."""
    path = write_oceanview_file(
        tmp_path,
        "midnight.txt",
        [
            "2026-08-05 23:59:59.000000",
            "2026-08-06 00:00:01.000000",
            "2026-08-06 00:00:03.000000",
        ],
    )

    optical = Spectrum.read(path, reader="oceanview")

    np.testing.assert_allclose(optical.t, [0.0, 2.0, 4.0])


def test_placeholder_rows_are_read_in_header_mode(tmp_path):
    """The fake date cancels out of the row-to-row differences."""
    path = write_oceanview_file(
        tmp_path,
        "placeholder.txt",
        ["1970-01-01 00:00:00.000000", "1970-01-01 00:00:02.000000"],
    )

    optical = Spectrum.read(path, reader="oceanview")

    assert optical.tstamp == unix_time(2026, 8, 5, 12, 0, 0, 0, 2)
    np.testing.assert_allclose(optical.t, [0.0, 2.0])
