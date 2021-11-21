"""A reader for text exports from the RGA Software of Stanford Instruments"""

from EC_MS import Dataset
from .ec_ms_pkl import measurement_from_ec_ms_dataset
from ..techniques import ECMeasurement


class CHInstrumentsTXTReader:
    path_to_file = None

    def read(self, path_to_file, cls=None):
        """Read a .txt file exported by CH Instruments software.

        TODO: Write a new reader that doesn't use the old EC_MS package

        Args:
            path_to_file (Path or str): The file to read
            cls (Measurement subclass): The class to return. Defaults to ECMeasuremnt
        """
        self.path_to_file = path_to_file
        cls = cls if (cls and not issubclass(ECMeasurement, cls)) else ECMeasurement
        ec_ms_dataset = Dataset(path_to_file, data_type="CHI")
        return measurement_from_ec_ms_dataset(
            ec_ms_dataset.data, cls=cls, reader=self, technique="EC"
        )
