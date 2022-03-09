"""For use in development of the ixdat .csv reader (and exporter)."""

from ixdat import Measurement


if False:  # test the version that's online on the tutorials page
    meas = Measurement.read_url(
        "https://raw.githubusercontent.com/ixdat/tutorials/"
        + "main/loading_appending_and_saving/co_strip.csv",
        reader="ixdat",
    )
    meas.plot_measurement()

else:
    meas = Measurement.read(
        "../../test_data/biologic/Pt_poly_cv_CUT.mpt", reader="biologic"
    )
    meas.calibrate_RE(0.715)

    meas.correct_ohmic_drop(R_Ohm=100)

    meas.normalize_current(0.196)

    meas.as_cv().export("test.csv")

    meas_loaded = Measurement.read("test.csv", reader="ixdat")

    meas_loaded.plot()

