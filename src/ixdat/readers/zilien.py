from . import TECHNIQUE_CLASSES
from .ec_ms_pkl import measurement_from_ec_ms_dataset

ECMSMeasurement = TECHNIQUE_CLASSES["EC-MS"]


class ZilienTSVReader:
    """Class for reading files saved by Spectro Inlets' Zilien software"""

    def read(self, path_to_file, cls=None, name=None, **kwargs):
        """Read a zilien file

        TODO: This is a hack using EC_MS to read the .tsv. Will be replaced.
        """
        cls = cls or ECMSMeasurement
        from EC_MS import Zilien_Dataset

        ec_ms_dataset = Zilien_Dataset(path_to_file)
        return measurement_from_ec_ms_dataset(
            ec_ms_dataset.data,
            cls=cls,
            name=name,
            reader=self,
            technique="EC-MS",
            **kwargs
        )


if __name__ == "__main__":
    """Module demo here.

    To run this module in PyCharm, open Run Configuration and set
        Module name = ixdat.readers.zilien,
    and *not*
        Script path = ...
    """

    from pathlib import Path
    from ixdat.measurements import Measurement

    path_to_test_file = Path.home() / (
        "Dropbox/ixdat_resources/test_data/"
        # "zilien_with_spectra/2021-02-01 14_50_40.tsv"
        "zilien_with_ec/2021-02-01 17_44_12.tsv"
    )

    ecms_measurement = Measurement.read(
        reader="zilien",
        path_to_file=path_to_test_file,
    )

    ecms_measurement.plot_measurement()
