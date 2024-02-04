"""This module performs top-level functional tests on the Zilien reader"""

import datetime
import time
from collections import namedtuple
from functools import reduce
from pathlib import Path

import numpy as np
import pytest
from pytest import approx, fixture

from ixdat import Measurement
from ixdat.techniques.ec_ms import (
    ECMSMeasurement,
    MSMeasurement,
    ECMeasurement,
)


DIR_VERSION_1 = (
    Path(__file__).parent.parent.parent
    / "submodules/ixdat-large-test-files/zilien_version_1"
)
DIR_VERSION_2 = (
    Path(__file__).parent.parent.parent
    / "submodules/ixdat-large-test-files/zilien_version_2"
)

PATH_TO_DATAFILE = DIR_VERSION_1 / "2022-04-06 16_17_23 full set.tsv"
TIME_FORMAT = "%Y-%m-%d %H_%M_%S"
EXPECTED_SERIES_ECMS = {
    "M32 [A]",
    "Pirani3 value time [s]",
    "Pirani2 value [mbar]",
    "MFC1 value time [s]",
    "PC value [mbar]",
    "Voltage [V]",
    "PC setpoint time [s]",
    "MFC2 setpoint [ml/min]",
    "C4M32 time [s]",
    "Current [mA]",
    "Pirani2 value time [s]",
    "Iongauge value [mbar]",
    "C3M28 time [s]",
    "PC setpoint [mbar]",
    "MFC2 setpoint time [s]",
    "M28 [A]",
    "M4 [A]",
    "MFC2 value [ml/min]",
    "Pirani1 value time [s]",
    "C2M18 time [s]",
    "MFC1 setpoint [ml/min]",
    "Potential time [s]",
    "M18 [A]",
    "Pirani1 value [mbar]",
    "MFC1 value [ml/min]",
    "MFC1 setpoint time [s]",
    "M2 [A]",
    "C0M2 time [s]",
    "Iongauge value time [s]",
    "Pirani3 value [mbar]",
    "MFC2 value time [s]",
    "Cycle [n]",
    "PC value time [s]",
    "C1M4 time [s]",
}
EXPECTED_SERIES_MS = {
    series
    for series in EXPECTED_SERIES_ECMS
    if series not in {"Potential time [s]", "Voltage [V]", "Current [mA]", "Cycle [n]"}
}
EXPECTED_SERIES_EC = {
    series for series in EXPECTED_SERIES_ECMS if series[2] != "M" and "[A]" not in series
}
TECHNIQUE_NAME_TO_CLASS = {
    "EC-MS": ECMSMeasurement,
    "EC": ECMeasurement,
    "MS": MSMeasurement,
}
series_details = namedtuple("series_details", "t_first,t_last,v_first,v_last,count")
SERIES_DETAILS = {
    "Voltage [V]": series_details(
        9.584931e01, 2.134163e02, 5.920536e-01, 1.014159e-01, 733
    ),
    "Current [mA]": series_details(
        9.584931e01, 2.134163e02, 1.055188e-06, -2.818254e-02, 733
    ),
    "Cycle [n]": series_details(9.584931e01, 2.134163e02, 1, 2, 733),
    "Iongauge value [mbar]": series_details(
        3.479784e-01, 2.171076e02, 6.160000e-07, 6.170000e-07, 457
    ),
    "MFC1 setpoint [ml/min]": series_details(
        5.353265e-01, 2.173076e02, 1.000000e00, 1.000000e00, 457
    ),
    "MFC1 value [ml/min]": series_details(
        5.353265e-01, 2.173076e02, 1.000000e00, 1.000000e00, 457
    ),
    "MFC2 setpoint [ml/min]": series_details(
        5.817895e-01, 2.173550e02, 0.000000e00, 0.000000e00, 457
    ),
    "MFC2 value [ml/min]": series_details(
        5.817895e-01, 2.173550e02, -6.200000e-02, -6.200000e-02, 457
    ),
    "PC setpoint [mbar]": series_details(
        6.286705e-01, 2.174025e02, 9.945849e02, 9.945849e02, 457
    ),
    "PC value [mbar]": series_details(
        6.286705e-01, 2.174025e02, 9.946649e02, 9.945849e02, 457
    ),
    "Pirani1 value [mbar]": series_details(
        3.951142e-01, 2.171688e02, 2.000000e-02, 1.990000e-02, 457
    ),
    "Pirani2 value [mbar]": series_details(
        4.423413e-01, 2.172148e02, 1.950000e-01, 1.980000e-01, 457
    ),
    "Pirani3 value [mbar]": series_details(
        4.888763e-01, 2.172614e02, 1.880000e-02, 1.870000e-02, 457
    ),
    "M2 [A]": series_details(4.798638e-01, 2.163797e02, 2.319726e-09, 2.187911e-09, 404),
    "M4 [A]": series_details(5.548562e-01, 2.164546e02, 3.650363e-08, 3.444133e-08, 404),
    "M18 [A]": series_details(
        6.288487e-01, 2.165286e02, 3.607096e-09, 3.401646e-09, 404
    ),
    "M28 [A]": series_details(
        7.038411e-01, 2.166036e02, 1.864088e-10, 1.828584e-10, 404
    ),
    "M32 [A]": series_details(
        7.778335e-01, 2.166776e02, 1.934676e-10, 1.851902e-10, 404
    ),
}

# the amount of Zilien data columns in the test datasets
# it is used when testing if the amount of created series
# by the "biologic" reader and by the "zilien" reader matches
ZILIEN_COLUMNS_AMOUNT = 56


# fmt: off
@fixture(scope="class", params=(
    (
        DIR_VERSION_2 / "2022-09-19 11_27_59 test-1/2022-09-19 11_27_59 test-21.tsv",
        [
            DIR_VERSION_2 / "2022-09-19 11_27_59 test-1/2022-09-19 11_27_59 test-21_01_01_CV_C01.mpt",  # noqa: E501
            DIR_VERSION_2 / "2022-09-19 11_27_59 test-1/2022-09-19 11_27_59 test-21_01_02_CP_C01.mpt",  # noqa: E501
        ],
        (1.344293e+01, 1.344293e+01)
    ),
    (
        DIR_VERSION_2 / "2022-09-19 11_29_04 test-2/2022-09-19 11_29_04 test-22.tsv",
        [
            DIR_VERSION_2 / "2022-09-19 11_29_04 test-2/2022-09-19 11_29_04 test-22_01_01_CV_C01.mpt",  # noqa: E501
            DIR_VERSION_2 / "2022-09-19 11_29_04 test-2/2022-09-19 11_29_04 test-22_01_02_CP_C01.mpt",  # noqa: E501
            DIR_VERSION_2 / "2022-09-19 11_29_04 test-2/2022-09-19 11_29_04 test-22_02_01_CV_C01.mpt",  # noqa: E501
            DIR_VERSION_2 / "2022-09-19 11_29_04 test-2/2022-09-19 11_29_04 test-22_02_02_CP_C01.mpt",  # noqa: E501
        ],
        (7.223014e+00, 7.223014e+00, 3.286301e+01, 3.286301e+01)
    ),
    (
        DIR_VERSION_2 / "2022-09-19 11_31_19 test-3/2022-09-19 11_31_19 test-23.tsv",
        (
            DIR_VERSION_2 / "2022-09-19 11_31_19 test-3/2022-09-19 11_31_19 test-23_01_C01.mpt",  # noqa: E501
            DIR_VERSION_2 / "2022-09-19 11_31_19 test-3/2022-09-19 11_31_19 test-23_02_C01.mpt",  # noqa: E501
        ),
        (3.294874e+00, 4.138740e+01)
    ),
))
# fmt: on
def datasets(request):
    """Parametrized fixture for reading datasets and return created series.

    One parameter set contains a path to the Zilien TSV file, a list of paths
    to the Biologic MPT files and a tuple with time offsets corresponding to
    the MPT files.

    Args:
        request: A pytest request class with the parameters.
    """
    tsv_file_path = request.param[0]
    mpt_file_paths = request.param[1]
    mpt_time_offsets = request.param[2]

    tsv = Measurement.read(tsv_file_path, reader="zilien")
    mpts = [Measurement.read(path, reader="biologic") for path in mpt_file_paths]

    return tsv, mpts, mpt_time_offsets


@pytest.mark.skip("TEMP SKIP. BROKEN, SEE: https://github.com/ixdat/ixdat/issues/158")
@pytest.mark.external
@pytest.mark.parametrize(
    ["cls_or_technique", "expected_series"],
    (
        (Measurement, EXPECTED_SERIES_ECMS),
        (ECMSMeasurement, EXPECTED_SERIES_ECMS),
        (MSMeasurement, EXPECTED_SERIES_MS),
        (ECMeasurement, EXPECTED_SERIES_EC),
        ("EC-MS", EXPECTED_SERIES_ECMS),
        ("EC", EXPECTED_SERIES_EC),
        ("MS", EXPECTED_SERIES_MS),
    ),
)
def test_read(cls_or_technique, expected_series):
    """Test parsing the file as an ECMSMeasurement. Zilien generated data only."""

    if isinstance(cls_or_technique, str):
        expected_measurement_class = TECHNIQUE_NAME_TO_CLASS[cls_or_technique]
        measurement = Measurement.read(
            PATH_TO_DATAFILE, reader="zilien", technique=cls_or_technique
        )
    else:
        measurement = cls_or_technique.read(PATH_TO_DATAFILE, reader="zilien")
        if cls_or_technique is Measurement:
            expected_measurement_class = ECMSMeasurement
        else:
            expected_measurement_class = cls_or_technique

    assert type(measurement) == expected_measurement_class

    # Check series_names and timestamp
    assert measurement.series_names == expected_series
    timestamp_string = " ".join(measurement.name.split(" ")[:2])
    timestamp = time.mktime(
        datetime.datetime.strptime(timestamp_string, TIME_FORMAT).timetuple()
    )
    assert measurement.tstamp == approx(timestamp)

    # Test series and series data properties
    for series_name in expected_series:
        if series_name.endswith("[s]"):
            continue
        single_series_details = SERIES_DETAILS[series_name]
        time_data = measurement[series_name].tseries.data
        value_data = measurement[series_name].data
        assert len(time_data) == len(value_data) == single_series_details.count
        assert time_data[0] == approx(single_series_details.t_first)
        assert time_data[-1] == approx(single_series_details.t_last)
        assert value_data[0] == approx(single_series_details.v_first)
        assert value_data[-1] == approx(single_series_details.v_last)


class TestZilienIntegrated:
    """Tests for reading the Zilien files with integrated Biologic dataset."""

    @pytest.mark.skip(
        "TEMP SKIP. BROKEN, SEE: https://github.com/ixdat/ixdat/issues/158"
    )
    @pytest.mark.external
    def test_series_match(self, datasets):
        """Test presence of all the series from the mpt files in the tsv file."""

        # coming from parametrized fixture
        tsv, mpts, _ = datasets

        # if one name from a mpt file would match a name from a different mpt
        # in the integrated tsv dataset, then the test for series amount will
        # fail, so it is not necessary to check for the exact belonging
        for mpt in mpts:
            for mpt_name in mpt.series_names:
                if mpt_name == "time/s":
                    assert f"Biologic {mpt_name}" in tsv.series_names
                else:
                    assert mpt_name in tsv.series_names

    @pytest.mark.external
    def test_series_amount(self, datasets):
        """Test the same amount of created series"""

        # coming from parametrized fixture
        tsv, mpts, _ = datasets

        tsv_amount = sum([len(mpt.series_list) for mpt in mpts])
        mpt_amount = len(tsv.series_list) - ZILIEN_COLUMNS_AMOUNT

        assert tsv_amount == mpt_amount

    @pytest.mark.external
    def test_timeseries_match(self, datasets):
        """Test the same timeseries."""

        # coming from parametrized fixture
        tsv, mpts, mpts_time_offsets = datasets

        # same amount
        tsv_time_series = [
            ser for ser in tsv.series_list if ser.name == "Biologic time/s"
        ]
        assert len(tsv_time_series) == len(mpts_time_offsets)

        # same data
        for i, offset in enumerate(mpts_time_offsets):
            tsv_times = tsv_time_series[i].data
            mpt_times = mpts[i]["time/s"].data + offset

            assert np.allclose(tsv_times, mpt_times)

    @pytest.mark.skip(
        "TEMP SKIP. BROKEN, SEE: https://github.com/ixdat/ixdat/issues/158"
    )
    @pytest.mark.external
    def test_data_match(self, datasets):
        """Test the same data."""

        # coming from parametrized fixture
        tsv, mpts, mpts_time_offsets = datasets

        all_mpts = reduce(lambda x, y: x + y, mpts)

        for name in all_mpts.series_names:
            if name == "time/s":
                continue

            assert np.allclose(all_mpts[name].data, tsv[name].data)
