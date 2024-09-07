"""Configuration variables including standard data directory

This module also defines user-specific options, such as whether to use plugins

The object `plugins` created here gives direct access to these options.
Example useage:

```
import ixdat

ixdat.config.plugins.use_si_quant = True  # use the spectro_inlets_quantification package
```

See `help(ixdat.options.plugins)` for information.
"""
import datetime
from pathlib import Path

# For back-compatibility until 0.3.1:
from .plugins import plugins  # noqa: F401


class _Config:
    """
    Attributes:
        standard_data_directory (Path): the directory in which to save by default with
            the default database backend, which is to save files in a directory.
            ixdat will make the directory if it does not exist.
        standard_metadata_suffix (str): The file ext. for JSON format metadata files
        standard_data_suffix (str): The file extension for numpy.save format data files
        timestamp_string_format (str): A format string for datetime.datetime.strftime.
            Defaults to ixdats custom datetime format: 22E18 14:34:55
        timezone (datetime.timezone): The timezone timestamps should use when formatted
            to string. Defaults to the current local timestamp.
    """

    def __init__(self):
        self.standard_metadata_suffix = ".ix"
        self.standard_data_suffix = ".ix.npy"
        self.standard_ixdat_directory = Path.home() / ".ixdat"
        self.standard_data_directory = self.standard_ixdat_directory / "projects"
        self.default_project_name = "test"
        self.timestamp_string_format = "native"
        self.timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo

    @property
    def ixdat_temp_dir(self):
        temp_dir = self.standard_data_directory / "temp"
        if not temp_dir.exists():
            temp_dir.mkdir(parents=True)
        return temp_dir


config = _Config()


def prompt_for_permission(prompt):
    yn = input(prompt + "\nEnter y for yes or anything else for no.")
    return yn in ["y", "yes", "Y", "Yes"]
