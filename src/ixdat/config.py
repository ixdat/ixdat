"""Configuration variables including standard data directory

"""
from pathlib import Path


class Config:
    """
    Attributes:
        standard_data_directory (Path): the directory in which to save by default with
            the default database backend, which is to save files in a directory.
            ixdat will make the directory if it does not exist.
        standard_metadata_suffix (str): The file ext. for JSON format metadata files
        standard_data_suffix (str): The file extension for numpy.save format data files
    """

    def __init__(self):
        self.standard_metadata_suffix = "ix"
        self.standard_data_suffix = "ixdata"
        self.standard_data_directory = Path.home() / "ixdat"

        if not self.standard_data_directory.exists():
            try:
                self.standard_data_directory.mkdir()
            except Exception:
                raise  # TODO, figure out what gets raised, then except with line below
                # raise ConfigError(f"Cannot make dir '{self.standard_data_directory}'")


CFG = Config()
