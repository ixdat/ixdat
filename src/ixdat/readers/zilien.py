from . import TECHNIQUE_CLASSES
from .ec_ms_pkl import measurement_from_ec_ms_dataset

ECMSMeasurement = TECHNIQUE_CLASSES["EC-MS"]


class ZilienTSVReader:
    def read(self, path_to_file, cls=None, name=None, **kwargs):
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

    from pathlib import Path
    from ixdat.measurements import Measurement

    path_to_test_file = Path.home() / (
        "Dropbox/ixdat_resources/test_data/"
        # "2021-02-01 14_50_40 Chip 2 test water/"
        # "2021-02-01 14_50_40 Chip 2 test water.tsv"
        "2021-02-01 17_44_12 NMC vs Li 0.1C/2021-02-01 17_44_12 NMC vs Li 0.tsv"
    )

    ecms_measurement = Measurement.read(
        reader="zilien",
        path_to_file=path_to_test_file,
    )

    ecms_measurement.plot_measurement()
