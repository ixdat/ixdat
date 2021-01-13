
"""Dummy test file, should be deleted as soon as we have the first real test"""


def test_dummy():
    pass


def test_blank_measurement():
    from ixdat.measurements import Measurement

    meas = Measurement(
        name="blank",
    )
    print(meas)
    assert len(meas.value_names) == 0


if __name__ == "__main__":
    test_blank_measurement()
