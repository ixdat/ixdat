"""A reader for text exports from the potentiostat software of CH Instruments"""

from .ec_ms_pkl import measurement_from_ec_ms_dataset
from ..techniques import MSMeasurement


class StanfordRGASoftReader:
    path_to_file = None

    def read(self, path_to_file, cls=None):
        """Read a .txt file exported by CH Instruments software.

        TODO: Write a new reader that doesn't use the old EC_MS package

        Args:
            path_to_file (Path or str): The file to read
            cls (Measurement subclass): The class to return. Defaults to ECMeasuremnt
        """
        try:
            from EC_MS import Dataset
        except ImportError:
            print(
                "The ixdat StanfordRGASoftReader relies on the EC_MS package.\n"
                "Use `pip install EC_MS`. \n"
                "Alternatively considering writing a new Reader for ixdat!"
            )

        # with open(path_to_file, "r") as f:
        #     timestamp_string = f.readline().strip()
        # tstamp = timestamp_string_to_tstamp(
        #     timestamp_string,
        #     form="%b %d, %Y  %I:%M:%S %p",  # like "Mar 05, 2020  09:50:34 AM"
        # )   # ^ For later. EC_MS actually gets this right.

        self.path_to_file = path_to_file
        cls = cls if (cls and not issubclass(MSMeasurement, cls)) else MSMeasurement
        ec_ms_dataset = Dataset(
            path_to_file,
            data_type="RGA",  # tstamp=tstamp
        )
        return measurement_from_ec_ms_dataset(
            ec_ms_dataset.data, cls=cls, reader=self, technique="MS"
        )
