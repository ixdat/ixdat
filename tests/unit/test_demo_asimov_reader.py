import importlib.util
from pathlib import Path

import pytest

DEMO_PATH = (
    Path(__file__).parents[2]
    / "development_scripts"
    / "reader_demonstrators"
    / "demo_asimov_reader.py"
)
spec = importlib.util.spec_from_file_location("demo_asimov_reader", DEMO_PATH)
demo_asimov_reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo_asimov_reader)
prepare_measurements_for_demo = demo_asimov_reader.prepare_measurements_for_demo


class DummyMeasurement:
    U_name = "U"
    J_name = "J"
    E_name = "E"
    I_name = "I"

    def __init__(self, name):
        self.name = name
        self.metadata = {}

    def grab(self, series_name):
        return [0, 1], [0, 1]

    def as_cv(self):
        return self

    def plot(self, *args, **kwargs):
        return None


class IncompatibleMeasurement(DummyMeasurement):
    def __init__(
        self, name, U_name="U", J_name="J", E_name="E", I_name="I", cv_error=None
    ):
        super().__init__(name)
        self.U_name = U_name
        self.J_name = J_name
        self.E_name = E_name
        self.I_name = I_name
        self.cv_error = cv_error

    def as_cv(self):
        if self.cv_error:
            raise self.cv_error
        return self

    def grab(self, series_name):
        if series_name in {"potential", self.U_name, "raw_potential", self.E_name}:
            return [0, 1], [0, 1]
        if series_name in {"current", self.J_name, "raw_current", self.I_name}:
            return [0, 1], [0, 1]
        raise demo_asimov_reader.SeriesNotFoundError(series_name)


class EmptySeriesMeasurement(IncompatibleMeasurement):
    def grab(self, series_name):
        return [], []


class MissingPotentialMeasurement(IncompatibleMeasurement):
    def grab(self, series_name):
        if series_name in {"potential", self.U_name, "raw_potential", self.E_name}:
            raise demo_asimov_reader.SeriesNotFoundError(series_name)
        return super().grab(series_name)


def test_prepare_measurements_for_demo_rejects_incompatible_measurements():
    measurements = [
        IncompatibleMeasurement("good", U_name="E", J_name="j"),
        MissingPotentialMeasurement("bad", U_name="U", J_name="j", I_name="j"),
    ]

    with pytest.raises(ValueError, match="cannot be plotted together"):
        prepare_measurements_for_demo(measurements)


def test_prepare_measurements_for_demo_reports_empty_series_and_cv_failure():
    measurements = [
        IncompatibleMeasurement("good", U_name="E", J_name="j"),
        EmptySeriesMeasurement(
            "empty",
            U_name="E",
            J_name="j",
            cv_error=ValueError("zero-size array to reduction operation maximum"),
        ),
    ]

    with pytest.raises(
        ValueError, match="empty voltage series for key 'potential'"
    ) as exception:
        prepare_measurements_for_demo(measurements)

    assert "cannot convert to CV: ValueError" in str(exception.value)


def test_prepare_measurements_for_demo_accepts_alias_equivalent_series_names():
    measurements = [
        IncompatibleMeasurement("a", U_name="E/V", J_name="<I>/mA"),
        IncompatibleMeasurement("b", U_name="raw_potential", J_name="I/mA"),
    ]

    reference_u_name, reference_j_name, cv_measurements = prepare_measurements_for_demo(
        measurements
    )

    assert reference_u_name == "potential"
    assert reference_j_name == "current"
    assert cv_measurements == measurements
