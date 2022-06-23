"""
The module tests correct parsing of mpt files generated by EC-lab with use of Zilien.
The dataset with the files exists as a submodule

NOTE: Test of PEIS technique is skipped, because it contains two columns with
the same name and therefore ixdat parses both of them as a one object that is not
possible to plot later, because the data list length differs from the rest of
the columns.
"""
from pathlib import Path

from ixdat import Measurement

from pytest import fixture, mark


DATA_DIR = (
    Path(__file__).parent.parent.parent / "submodules/ixdat-large-test-files/multiple_techniques_dataset/"
)

"""
filename_with_technique: {
  "column_name": ("unit", (first_column_value, last_column_value)),
  ...
}
"""
TEST_DATA = {
    "multiple_techniques_dataset_01_02_CVA_C01.mpt": {
        "mode": ("", (2, 2)),
        "ox/red": ("red", (0, 1)),
        "error": ("", (0, 0)),
        "control changes": ("", (0, 0)),
        "Ns changes": ("", (0, 0)),
        "counter inc.": ("", (0, 0)),
        "time/s": ("s", (2.868160128949967E+001, 5.592763878855476E+002)),
        "control/V": ("V", (9.1791922E-001, 8.9992851E-001)),
        "Ewe/V": ("V", (9.1770262E-001, 9.0080881E-001)),
        "<I>/mA": ("mA", (-4.376391007099301E-004, 1.126740008349650E-002)),
        "cycle number": ("", (1.000000000000000E+000, 6.000000000000000E+000)),
        "(Q-Qo)/C": ("C", (0.0000000E+000, 1.4638975E-003)),
        "I Range": ("", (41, 41)),
        "Rcmp/Ohm": ("Ohm", (0.0000000E+000, 0.0000000E+000)),
        "P/W": ("W", (-4.0162254E-007, 1.0149774E-005)),
        "Ns": ("", (0, 0))  # this is automatically added by ixdat as a ConstantValue
    },
    "multiple_techniques_dataset_01_03_CP_C01.mpt": {
        "mode": ("", (1, 1)),
        "ox/red": ("red", (1, 0)),
        "error": ("", (0, 0)),
        "control changes": ("", (1, 1)),
        "Ns changes": ("", (0, 0)),
        "counter inc.": ("", (0, 0)),
        "Ns": ("", (0, 2)),
        "time/s": ("s", (5.602763878602855E+002, 5.893783871251071E+002)),
        "control/mA": ("mA", (0.0000000E+000, -2.0000001E-003)),
        "<Ewe>/V": ("V", (8.8483948E-001, 8.2798851E-001)),
        "I/mA": ("mA", (1.2781935E-006, -2.0001377E-003)),
        "dQ/C": ("C", (7.1232675E-010, -1.9990615E-005)),
        "(Q-Qo)/C": ("C", (7.1232675E-010, 5.0315875E-006)),
        "half cycle": ("", (0, 1)),
        "Q charge/discharge/mA.h": ("mA.h", (1.978685425384362E-010, -5.552948449702752E-006)),
        "I Range": ("", (42, 42)),
        "Rcmp/Ohm": ("Ohm", (0.0000000E+000, 0.0000000E+000)),
        "Energy charge/W.h": ("W.h", (1.750818973293703E-013, 6.555031900229182E-009)),
        "Energy discharge/W.h": ("W.h", (0.000000000000000E+000, -4.919033634676710E-009)),
        "Capacitance charge/µF": ("µF", (0.000000000000000E+000, 2.273762151181344E+002)),
        "Capacitance discharge/µF": ("µF", (0.000000000000000E+000, 1.369492186275801E+002)),
        "Q discharge/mA.h": ("mA.h", (0.000000000000000E+000, 5.552948449702752E-006)),
        "Q charge/mA.h": ("mA.h", (1.978685425384362E-010, 0.000000000000000E+000)),
        "Capacity/mA.h": ("mA.h", (1.978685425384362E-010, 5.552948449702752E-006)),
        "Efficiency/%": ("%", (0.0000000E+000, 7.9891510E+001)),
        "cycle number": ("", (0.000000000000000E+000, 0.000000000000000E+000)),
        "P/W": ("W", (1.1309961E-009, -1.6560911E-006)),
    },
    "multiple_techniques_dataset_01_04_CP_C01.mpt": {
        "mode": ("", (1, 1)),
        "ox/red": ("red", (0, 0)),
        "error": ("", (0, 0)),
        "control changes": ("", (0, 0)),
        "Ns changes": ("", (0, 0)),
        "counter inc.": ("", (0, 0)),
        "Ns": ("", (0, 0)),
        "time/s": ("s", (5.894805871225253E+002, 5.994803868699091E+002)),
        "control/mA": ("mA", (-2.0000001E-003, -2.0000001E-003)),
        "Ewe/V": ("V", (8.2613617E-001, 6.2719619E-001)),
        "I/mA": ("mA", (-1.9717729E-003, -1.9393268E-003)),
        "dQ/C": ("C", (0.0000000E+000, -1.9855208E-005)),
        "(Q-Qo)/C": ("C", (0.0000000E+000, -1.9855208E-005)),
        "half cycle": ("", (0, 0)),
        "Q charge/discharge/mA.h": ("mA.h", (0.000000000000000E+000, -5.515335538398682E-006)),
        "I Range": ("", (41, 41)),
        "Rcmp/Ohm": ("Ohm", (0.0000000E+000, 0.0000000E+000)),
        "Energy charge/W.h": ("W.h", (0.000000000000000E+000, 0.000000000000000E+000)),
        "Energy discharge/W.h": ("W.h", (0.000000000000000E+000, -4.130971864776140E-009)),
        "Capacitance charge/µF": ("µF", (0.000000000000000E+000, 0.000000000000000E+000)),
        "Capacitance discharge/µF": ("µF", (0.000000000000000E+000, 1.050930137626443E+002)),
        "Q discharge/mA.h": ("mA.h", (0.000000000000000E+000, 5.515335538398682E-006)),
        "Q charge/mA.h": ("mA.h", (0.000000000000000E+000, 0.000000000000000E+000)),
        "Capacity/mA.h": ("mA.h", (0.000000000000000E+000, 5.515335538398682E-006)),
        "Efficiency/%": ("%", (0.0000000E+000, 0.0000000E+000)),
        "cycle number": ("", (0.000000000000000E+000, 0.000000000000000E+000)),
        "P/W": ("W", (-1.6289529E-006, -1.2163383E-006)),
    },
    "multiple_techniques_dataset_01_05_CA_C01.mpt": {
        "mode": ("", (2, 2)),
        "ox/red": ("red", (0, 0)),
        "error": ("", (1, 0)),
        "control changes": ("", (0, 0)),
        "Ns changes": ("", (0, 0)),
        "counter inc.": ("", (0, 0)),
        "Ns": ("", (0, 2)),
        "time/s": ("s", (6.005809868421056E+002, 6.295807861095091E+002)),
        "control/V": ("V", (5.0035376E-002, 1.9995773E-001)),
        "Ewe/V": ("V", (4.9797375E-002, 2.0026733E-001)),
        "<I>/mA": ("mA", (-3.640213997266856E-002, -1.001744437708906E-003)),
        "dQ/C": ("C", (-3.6402140E-005, -1.4945253E-004)),
        "(Q-Qo)/C": ("C", (-3.6402140E-005, 4.7911013E-005)),
        "I Range": ("", (41, 41)),
        "Rcmp/Ohm": ("Ohm", (0.0000000E+000, 0.0000000E+000)),
        "Q charge/discharge/mA.h": ("mA.h", (-1.011170550352997E-005, -4.151459001554434E-005)),
        "half cycle": ("", (0, 3)),
        "Energy charge/W.h": ("W.h", (0.000000000000000E+000, 1.143620235377176E-007)),
        "Energy discharge/W.h": ("W.h", (-5.035363883729353E-010, -8.224704148866127E-009)),
        "Capacitance charge/µF": ("µF", (0.000000000000000E+000, 1.845039442322728E+002)),
        "Capacitance discharge/µF": ("µF", (0.000000000000000E+000, 2.743783797316349E+004)),
        "Q discharge/mA.h": ("mA.h", (1.011170550352997E-005, 4.151459001554434E-005)),
        "Q charge/mA.h": ("mA.h", (0.000000000000000E+000, 0.000000000000000E+000)),
        "Capacity/mA.h": ("mA.h", (1.011170550352997E-005, 4.151459001554434E-005)),
        "Efficiency/%": ("%", (0.0000000E+000, 5.8165386E+001)),
        "cycle number": ("", (0.000000000000000E+000, 1.000000000000000E+000)),
        "P/W": ("W", (-1.8127311E-006, -2.0061668E-007)),
    },
    "multiple_techniques_dataset_01_06_CA_C01.mpt": {
        "mode": ("", (2, 2)),
        "ox/red": ("red", (0, 0)),
        "error": ("", (0, 0)),
        "control changes": ("", (0, 0)),
        "Ns changes": ("", (0, 0)),
        "counter inc.": ("", (0, 0)),
        "Ns": ("", (0, 0)),
        "time/s": ("s", (6.296813861069677E+002, 6.396811858543515E+002)),
        "control/V": ("V", (5.0035376E-002, 5.0035376E-002)),
        "Ewe/V": ("V", (1.9815466E-001, 5.0038058E-002)),
        "I/mA": ("mA", (-2.5009227E-001, -1.8450420E-003)),
        "dQ/C": ("C", (0.0000000E+000, -4.4130393E-005)),
        "(Q-Qo)/C": ("C", (0.0000000E+000, -4.4130393E-005)),
        "I Range": ("", (41, 41)),
        "Rcmp/Ohm": ("Ohm", (0.0000000E+000, 0.0000000E+000)),
        "Q charge/discharge/mA.h": ("mA.h", (0.000000000000000E+000, -6.561739509278495E-006)),
        "half cycle": ("", (0, 11)),
        "Energy charge/W.h": ("W.h", (0.000000000000000E+000, 3.031160052865908E-013)),
        "Energy discharge/W.h": ("W.h", (0.000000000000000E+000, -3.327533907569938E-010)),
        "Capacitance charge/µF": ("µF", (0.000000000000000E+000, 0.000000000000000E+000)),
        "Capacitance discharge/µF": ("µF", (0.000000000000000E+000, 5.665210811458426E+005)),
        "Q discharge/mA.h": ("mA.h", (0.000000000000000E+000, 6.561739509278495E-006)),
        "Q charge/mA.h": ("mA.h", (0.000000000000000E+000, 0.000000000000000E+000)),
        "Capacity/mA.h": ("mA.h", (0.000000000000000E+000, 6.561739509278495E-006)),
        "Efficiency/%": ("%", (0.0000000E+000, 1.0970160E+005)),
        "cycle number": ("", (0.000000000000000E+000, 5.000000000000000E+000)),
        "P/W": ("W", (-4.9556947E-005, -9.2322317E-008)),
    },
    "multiple_techniques_dataset_01_07_ZIR_C01.mpt": {
        "freq/Hz": ("Hz", None),
        "Re(Z)/Ohm": ("Ohm", None),
        "-Im(Z)/Ohm": ("Ohm", None),
        "|Z|/Ohm": ("Ohm", None),
        "Phase(Z)/deg": ("deg", None),
        "time/s": ("s", None),
        "<Ewe>/V": ("V", None),
        "<I>/mA": ("mA", None),
        "I Range": ("", None),
        "|Ewe|/V": ("V", None),
        "|I|/A": ("A", None),
        "Re(Y)/Ohm-1": ("Ohm-1", None),
        "Im(Y)/Ohm-1": ("Ohm-1", None),
        "|Y|/Ohm-1": ("Ohm-1", None),
        "Phase(Y)/deg": ("deg", None),
    },
    "multiple_techniques_dataset_01_08_CVA_C01.mpt": {
        "mode": ("", (2, 2)),
        "ox/red": ("red", (0, 0)),
        "error": ("", (1, 0)),
        "control changes": ("", (1, 0)),
        "Ns changes": ("", (0, 0)),
        "counter inc.": ("", (0, 0)),
        "time/s": ("s", (6.409993858441158E+002, 1.426767165993951E+003)),
        "control/V": ("V", (-6.1258295E-005, 8.7560779E-001)),
        "Ewe/V": ("V", (-3.8347600E-004, 8.7638128E-001)),
        "<I>/mA": ("mA", (-2.500922679901123E-001, 2.756041360240096E-004)),
        "cycle number": ("", (1.000000000000000E+000, 7.000000000000000E+000)),
        "(Q-Qo)/C": ("C", (0.0000000E+000, -3.1572717E-004)),
        "I Range": ("", (41, 41)),
        "Rcmp/Ohm": ("Ohm", (0.0000000E+000, 0.0000000E+000)),
        "P/W": ("W", (9.5904383E-008, 2.4153431E-007)),
    }
}


@fixture(scope="module", params=tuple(TEST_DATA.items()))
@mark.submodule_multiple_techniques
def measurements_with_data(request):
    filename = request.param[0]
    column_data = request.param[1]
    measurement = Measurement.read(DATA_DIR / filename, reader="biologic")
    return measurement, column_data


@mark.submodule_multiple_techniques
def test_shape(measurements_with_data):
    measurement = measurements_with_data[0]
    column_names = measurements_with_data[1].keys()

    # test amount of columns
    assert len(measurement.series_list) == len(column_names)

    # test names of columns
    for name in column_names:
        if name not in measurement.series_names:
            assert False


@mark.submodule_multiple_techniques
def test_units(measurements_with_data):
    measurement = measurements_with_data[0]
    column_data = measurements_with_data[1]

    for column_name, unit_values in column_data.items():
        expected_unit = unit_values[0]

        # TODO do this easier
        #  Does Ixdat have similar method to `grab()`,
        #  that returns an object from the `series_list`?
        for series in measurement.series_list:
            if column_name == series.name:
                assert expected_unit == series.unit_name
                break


@mark.submodule_multiple_techniques
def test_values(measurements_with_data):
    measurement = measurements_with_data[0]
    column_data = measurements_with_data[1]

    # test first and last values in columns
    for column_name, unit_values in column_data.items():
        first_value = unit_values[1][0]
        last_value = unit_values[1][1]
        parsed_values = measurement.grab(column_name)[1]

        assert first_value == parsed_values[0]
        assert last_value == parsed_values[-1]
