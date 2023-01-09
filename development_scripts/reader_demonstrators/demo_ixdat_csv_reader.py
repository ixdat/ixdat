"""For use in development of the ixdat .csv reader (and exporter)."""

from ixdat import Measurement


if True:  # test back-compatability with the ixdat v0p1 online on the tutorials page
    meas = Measurement.read_url(
        "https://raw.githubusercontent.com/ixdat/tutorials/"
        + "ixdat_v0p1/loading_appending_and_saving/co_strip.csv",
        reader="ixdat",
        aliases={
            "t": ["time/s"],
            "raw_current": ["raw current / [mA]"],
            "raw_potential": ["raw potential / [V]"],
        },
    )
    meas.plot_measurement()
elif False:  # test with the ixdat v0p2 version online on the tutorials page
    meas = Measurement.read_url(
        "https://raw.githubusercontent.com/ixdat/tutorials/"
        + "main/electrochemistry/data/co_strip.csv",
        reader="ixdat",
    )
    meas.plot_measurement()

else:
    meas = Measurement.read(
        "../../test_data/biologic/Pt_poly_cv_CUT.mpt", reader="biologic"
    )
    meas.calibrate_RE(0.01)

    meas.correct_ohmic_drop(R_Ohm=100)

    meas.normalize_current(0.196)

    cv = meas.as_cv()

    cv.export("test.csv")

    meas_loaded = Measurement.read("test.csv", reader="ixdat")

    meas_loaded.plot()
