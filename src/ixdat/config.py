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
        self.standard_metadata_suffix = ".ix"
        self.standard_data_suffix = ".ix.npy"
        self.standard_data_directory = Path.home() / "ixdat"
        self.default_project_name = "test"


CFG = Config()


def prompt_for_permission(prompt):
    yn = input(prompt + "\nEnter y for yes or anything else for no.")
    return yn in ["y", "yes", "Y", "Yes"]
