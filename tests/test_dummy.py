"""Dummy test file, should be deleted as soon as we have the first real test"""


def test_dummy():
    pass


def test_blank_measurement():  # FIXME: tox can't handle matplotlib

    from ixdat.measurements import Measurement

    import matplotlib

    print(f"mpl version is {matplotlib.__version__}")

    meas = Measurement(
        name="blank",
    )
    print(meas)
    assert len(meas.value_names) == 0


if __name__ == "__main__":
    test_dummy()
    test_blank_measurement()
