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
from .tools import deprecate


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


class _PluginOptions:
    """A class for activating plugins and giving access to items from plugins

    This class has only one instance, initiated on import of ixdat. Use it by
    `from ixdat.config import plugins`

    These packages need to be separately installed.

    Packages
    --------

    spectro_inlets_quantification
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    - use_si_quant (bool): Read-only. If this is True, ixdat uses the
        `spectro_inlets_quantification` package. This changes the behaviour of some
        methods in `MSMeasurement` and inheriting classes. Defaults to False.
    - activate_si_quant(): Sets use_si_quant to True and initializes the quant. package
    - deactivate_si_quant(): Sets use_si_quant to False.
    - si_quant (_QuantDeps). This is where the needed imports from the external
        quantification package go. See the docstring of _QuantDeps
    """

    def __init__(self):
        self._use_siq = False
        self.siq = _SIQ()
        self.cinfdata = _CinfData()

    @property
    def use_siq(self):
        return self._use_siq

    def activate_siq(self):
        """Changes mass spec methods to use spectro_inlets_quantification"""
        self._use_siq = True
        self.siq.populate()

    def deactivate_siq(self):
        """Changes mass spec methods to use ixdat's native MS quantification"""
        self._use_siq = False

    def activate_cinfdata(self):
        self.cinfdata.activate()

    @deprecate(
        last_supported_release="0.2.6",
        update_message="`siq` is the universal abreviation "
        "for `spectro_inlets_quantification`. "
        "Thus, `ixdat.plugins.activate_si_quant` "
        "is now `ixdat.plugins.activate_siq`",
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def activate_si_quant(self):
        return self.activate_siq()

    @property
    @deprecate(
        last_supported_release="0.2.6",
        update_message="`siq` is the universal abreviation "
        "for `spectro_inlets_quantification`. "
        "Thus, `ixdat.plugins.use_si_quant` "
        "is now `ixdat.plugins.use_siq`",
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def use_si_quant(self):
        return self.use_siq

    @property
    @deprecate(
        last_supported_release="0.2.6",
        update_message="`siq` is the universal abreviation "
        "for `spectro_inlets_quantification`. "
        "Thus, `ixdat.plugins.si_quant` is now `ixdat.plugins.siq`",
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def si_quant(self):
        return self.siq


class _SIQ:
    """Class storing items of `spectro_inlets_quantification`.

    This class has one instance, which is an attribute of `plugins`. To print this
    docstring, you would type:
    ```
    from ixdat.config import plugins

    plugins.use_si_quant = True  # Activates plugins.quant.
    help(plugins.si_quant)  #gives information on the si_quant package
    ```

    The attributes of this class are `None` until the property `plugins.use_si_quant` is
    set to True, triggering their population (activating quant).

    Once activated, the attributes of `plugins.quant` are:
    - `Chip`: Class describing the MS inlet. More powerful than ixdat's `MSInlet`
    - `Molecule`: Class with data about molecules relevant to (EC-)MS quantification
    - `CalPoint`: Class with data about an MS calibration experiment
    - `Calibration`: Class for storing, visualizing, and using multiple CalPoints.
    - `Quantifier`: Class for using a Calibration to quantify MS data
    - `quant_config`: The config object of the external quantification package

    `plugins.quant` also has
    - `QUANT_DIRECTORY`: A property for getting and setting the data directory used by
      the external quantification package.
    """

    def __init__(self):
        self.Chip = None
        self.Molecule = None
        self.Calibration = None
        self.CalPoint = None
        self.Quantifier = None
        self.quant_config = None
        self._QUANT_DIRECTORY = None

    def populate(self):
        from spectro_inlets_quantification.config import Config
        from spectro_inlets_quantification.medium import Medium
        from spectro_inlets_quantification.molecule import Molecule
        from spectro_inlets_quantification.chip import Chip
        from spectro_inlets_quantification.calibration import CalPoint, Calibration
        from spectro_inlets_quantification.quantifier import Quantifier

        self.quant_config = Config()
        self.medium = Medium()
        self.Molecule = Molecule
        self.Chip = Chip
        self.CalPoint = CalPoint
        self.Calibration = Calibration
        self.Quantifier = Quantifier

    @property
    def QUANT_DIRECTORY(self):
        if not self._QUANT_DIRECTORY:
            self._QUANT_DIRECTORY = (
                config.standard_ixdat_directory / "plugin_data/ms_quant"
            )
            if not self._QUANT_DIRECTORY.exists():
                self._QUANT_DIRECTORY.mkdir(parents=True)
        return self._QUANT_DIRECTORY

    @QUANT_DIRECTORY.setter
    def QUANT_DIRECTORY(self, quant_directory):
        self._QUANT_DIRECTORY = quant_directory
        self.quant_config.aux_data_directory = self.QUANT_DIRECTORY


class _CinfData:
    """Class implement direct database read access using external module Cinfdata"""

    def __init__(
        self,
    ):
        self.cinfdata = None
        self._managed_cinfdata_object = None
        self._context_manager_kwargs = None

    def activate(self):
        from cinfdata import Cinfdata

        self.DB = Cinfdata

    def connect(self, setup_name=None, grouping_column=None):
        """setup_name (str): The name of the table inside the database
        grouping_column (str): Either the 'timestamp'/'comment' or 'Comment' column"""
        return self.DB(setup_name=setup_name, grouping_column=grouping_column)

    def __call__(self, **kwargs):
        """**kwargs: setup_name (str) and grouping_column (str) (see connect())"""
        self._context_manager_kwargs = kwargs
        return self

    def __enter__(self):
        if self._managed_cinfdata_object:
            self._context_manager_kwargs = None
            raise RuntimeError(
                "Using the cinfdata plugin as a context manager more "
                "than once at the same time is not supported"
            )

        self._managed_cinfdata_object = self.connect(**self._context_manager_kwargs)
        return self._managed_cinfdata_object

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._managed_cinfdata_object.connection.close()
        self._managed_cinfdata_object = None
        self._context_manager_kwargs = None


plugins = _PluginOptions()
