"""Configuration variables including standard data directory

This module also defines user-specific options, such as whether to use plugins

The object `plugins` created here gives direct access to these options.
Example useage:

```
import ixdat

ixdat.config.plugins.USE_QUANT = True  # use the spectro_inlets_quantification package
```

See `help(ixdat.options.plugins)` for information.
"""

from pathlib import Path


class _Config:
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

    @property
    def ixdat_temp_dir(self):
        temp_dir = self.standard_data_directory / "temp"
        if not temp_dir.exists():
            temp_dir.mkdir(parents=True)
        return temp_dir


CFG = _Config()


def prompt_for_permission(prompt):
    yn = input(prompt + "\nEnter y for yes or anything else for no.")
    return yn in ["y", "yes", "Y", "Yes"]


class _PluginOptions:
    """A class to store preferences and information about plugin packages

    These packages need to be separately installed.

    Packages
    --------

    spectro_inlets_quantification:
        - USE_QUANT (bool): Set this to True to use the `spectro_inlets_quantification`
            package. This changes the behaviour of some methods in `MSMeasurement` and
            inheriting classes. Defaults to False.
        - QUANT_DIRECTORY (Path): Set this to the location of your data directory for
            `spectro_inlets_quantification`. The QUANT_DIRECTORY should be a folder with
            a subfolder called "molecules", into which molecule data files will go.
            Defaults to a small set of molecular data located within the ixdat source
            code: "ixdat/src/ixdat/plugin/data/ms_quant"
    """

    def __init__(self):
        self._USE_QUANT = False
        self._QUANT_DIRECTORY = None

    @property
    def USE_QUANT(self):
        return self._USE_QUANT

    @USE_QUANT.setter
    def USE_QUANT(self, use_quant):
        self._USE_QUANT = use_quant
        if self._USE_QUANT:

            from spectro_inlets_quantification.config import Config

            quant_config = Config()
            quant_config.data_directory = self.QUANT_DIRECTORY

    @property
    def QUANT_DIRECTORY(self):
        if not self._QUANT_DIRECTORY:
            self._QUANT_DIRECTORY = Path(__file__).parent / "plugin_data/ms_quant"
        return self._QUANT_DIRECTORY

    @QUANT_DIRECTORY.setter
    def QUANT_DIRECTORY(self, quant_directory):
        self._QUANT_DIRECTORY = quant_directory
        if self._USE_QUANT:
            self._USE_QUANT = True  # gets the quant directory to be reset


plugins = _PluginOptions()
