from ixdat import Measurement

meas = Measurement.read_url(
    "https://raw.githubusercontent.com/ixdat/tutorials/"
    + "main/loading_appending_and_saving/co_strip.csv",
    reader="ixdat",
)
meas.plot_measurement()
