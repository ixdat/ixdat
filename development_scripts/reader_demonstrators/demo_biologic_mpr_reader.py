from pathlib import Path
from ixdat import Measurement
import matplotlib.pyplot as plt


# confirm <Ewe/V> gets normalized to <Ewe>/V>
def check_ewe_v_column_normalization(meas, file):
    """Check that <Ewe/V> was normalized to <Ewe>/V> in the given measurement."""
    if "<Ewe/V>" in meas.series_names:
        # This file is relevant for the normalization test
        print(f"\n>>> Testing column normalization in: {file.name}")
        assert (
            "<Ewe>/V" in meas.series_names
        ), "Expected normalized name <Ewe>/V not found"
        assert (
            "<Ewe/V>" not in meas.series_names
        ), "Unnormalized name <Ewe/V> should be normalized"
        assert meas.aliases.get("WE_potential") == [
            "<Ewe>/V"
        ], "Alias for WE_potential incorrect"
        assert len(meas["<Ewe>/V"]) > 0, "<Ewe>/V series should have data"
        print("Column normalization passed!\n")


for data_dir in [
    # Path.home() / "Dropbox/ixdat_resources/test_data/biologic/17J04_Pt_isotope_exchange",
    # Path.home() / "Dropbox/ixdat_resources/test_data/biologic/22I27_London",
    # Path.home() / "Dropbox/ixdat_resources/test_data/biologic/22K14_Tempo",
    Path(__file__).parent
    / "../../test_data/biologic",
]:
    combined_meas = None

    for file in data_dir.iterdir():
        if not file.suffix in [".mpr", ".mpt"]:
            continue
        meas = Measurement.read(file, reader="biologic")
        print(meas)
        print("... was read successfully!\n\n")

        #  run column normalization test if needed
        check_ewe_v_column_normalization(meas, file)

        meas.plot()
        if combined_meas:
            combined_meas = combined_meas + meas
        else:
            combined_meas = meas

    if combined_meas:
        combined_meas.plot()
        combined_meas.plot(J_name="selector")
        plt.show()  # force QT to plot even in WSL
