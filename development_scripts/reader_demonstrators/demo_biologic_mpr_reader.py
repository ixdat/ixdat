from tools_for_demos import DEMO_DATA_DIR
from ixdat import Measurement

combined_measurement_list = []

for folder in [
    "17J04_Pt_isotope_exchange",
    "22I27_London",
    "22K14_Tempo",
]:
    combined_meas = None
    data_dir = DEMO_DATA_DIR / "biologic" / folder

    for file in data_dir.iterdir():
        if not file.suffix == ".mpr":
            continue
        meas = Measurement.read(file, reader="biologic")
        print(meas)
        print("... was read successfully!\n\n")
        # meas.plot()
        if combined_meas:
            combined_meas = combined_meas + meas
        else:
            combined_meas = meas

    combined_meas.plot()
    combined_meas.plot(J_name="selector")

    combined_measurement_list.append(combined_meas)
