"""Dummy test file, should be deleted as soon as we have the first real test"""


def test_dummy():
    pass


def dont_test_blank_measurement():  # FIXME: tox can't handle matplotlib

    from ixdat.measurements import Measurement

    meas = Measurement(
        name="blank",
    )
    print(meas)
    assert len(meas.value_names) == 0


if __name__ == "__main__":
    test_dummy()
    # test_blank_measurement()
