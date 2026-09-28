"""Demo of the biologic .mpr reader. Requires the demo data."""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir


def read_mpr_folder(folder, data_dir=None):
    """Return an ECMeasurement combining all the .mpr files in a biologic folder."""
    root = get_demo_data_dir(data_dir)
    combined_meas = None
    for file in sorted((root / "biologic" / folder).iterdir()):
        if not file.suffix == ".mpr":
            continue
        meas = Measurement.read(file, reader="biologic")
        print(f"read {file.name}")
        if combined_meas:
            combined_meas = combined_meas + meas
        else:
            combined_meas = meas
    return combined_meas


def load_pt_isotope_exchange(data_dir=None):
    """Return an ECMeasurement from a Pt isotope exchange experiment."""
    return read_mpr_folder("17J04_Pt_isotope_exchange", data_dir=data_dir)


def load_london(data_dir=None):
    """Return an ECMeasurement from the 22I27_London folder."""
    return read_mpr_folder("22I27_London", data_dir=data_dir)


def load_tempo(data_dir=None):
    """Return an ECMeasurement from a TEMPO experiment."""
    return read_mpr_folder("22K14_Tempo", data_dir=data_dir)


EXAMPLES = {
    "pt_isotope_exchange": load_pt_isotope_exchange,
    "london": load_london,
    "tempo": load_tempo,
}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    if show:
        for meas in results.values():
            meas.plot()
            meas.plot(J_name="selector")
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
