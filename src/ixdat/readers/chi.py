"""A reader for text exports from the RGA Software of Stanford Instruments"""

from .ec_ms_pkl import measurement_from_ec_ms_dataset
from ..techniques import ECMeasurement


CHI_LEGACY_ALIASES = {
    # TODO: These should change to what Zilien calls them. Right now the alias's
    #   reflect the way the lagacy EC_MS code renames essential series
    "t": ["time/s"],
    "raw_potential": ["Ewe/V", "<Ewe>/V"],
    "raw_current": ["I/mA", "<I>/mA"],
    "cycle": ["cycle number"],
}


class CHInstrumentsTXTReader:
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
                "The ixdat CHInstrumentsTXTReader relies on the EC_MS package.\n"
                "Use `pip install EC_MS`. \n"
                "Alternatively considering writing a new Reader for ixdat!"
            )

        self.path_to_file = path_to_file
        cls = cls if (cls and not issubclass(ECMeasurement, cls)) else ECMeasurement
        ec_ms_dataset = Dataset(path_to_file, data_type="CHI")
        return measurement_from_ec_ms_dataset(
            ec_ms_dataset.data,
            cls=cls,
            reader=self,
            technique="EC",
            aliases=CHI_LEGACY_ALIASES,
        )
