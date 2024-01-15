from pathlib import Path
from ixdat import Measurement

for data_dir in [
    Path.home() / "Dropbox/ixdat_resources/test_data/biologic/17J04_Pt_isotope_exchange",
    Path.home() / "Dropbox/ixdat_resources/test_data/biologic/22I27_London",
    Path.home() / "Dropbox/ixdat_resources/test_data/biologic/22K14_Tempo",
]:

    combined_meas = None

    for file in data_dir.iterdir():
        if not file.suffix == ".mpr":
            continue
        meas = Measurement.read(file, reader="biologic")
        print(meas)
        print("... was read successfully!")
        meas.plot()
        if combined_meas:
            combined_meas = combined_meas + meas
        else:
            combined_meas = meas

    combined_meas.plot()
    combined_meas.plot(J_name="selector")
