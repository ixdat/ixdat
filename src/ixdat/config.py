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
    """

    def __init__(self):
        self.standard_metadata_suffix = ".ix"
        self.standard_data_suffix = ".ix.npy"
        self.standard_data_directory = Path.home() / ".ixdat"
        self.default_project_name = "test"

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
        self._use_si_quant = False
        self.si_quant = _SIQuantDeps()
        self.cinfdata = _CinfData()

    @property
    def use_si_quant(self):
        return self._use_si_quant

    def activate_si_quant(self):
        """Changes mass spec methods to use spectro_inlets_quantification"""
        self._use_si_quant = True
        self.si_quant.populate()

    def deactivate_si_quant(self):
        """Changes mass spec methods to use ixdat's native MS quantification"""
        self._use_si_quant = False

    @property
    @deprecate(
        last_supported_release="0.2.4",
        update_message="This has changed to lower-case: `use_si_quant`",
        remove_release="0.2.5",
    )
    def USE_QUANT(self):
        return self.use_si_quant

    @USE_QUANT.setter
    @deprecate(
        last_supported_release="0.2.4",
        update_message="Use the methods `activate_si_quant()` and "
        "`deactivate_si_quant()` instead.",
        remove_release="0.2.5",
    )
    def USE_QUANT(self, use_quant):
        if use_quant:
            self.activate_si_quant()
        else:
            self.deactivate_si_quant()

    def activate_cinfdata(self):
        self.cinfdata.activate()


class _SIQuantDeps:
    """Class storing items of the external MS quantification package.

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
        from spectro_inlets_quantification.chip import Chip
        from spectro_inlets_quantification.molecule import Molecule
        from spectro_inlets_quantification.calibration import Calibration, CalPoint
        from spectro_inlets_quantification.quantifier import Quantifier
        from spectro_inlets_quantification.config import Config

        self.Chip = Chip
        self.Molecule = Molecule
        self.Calibration = Calibration
        self.CalPoint = CalPoint
        self.Quantifier = Quantifier

        self.quant_config = Config()
        self.quant_config.data_directory = self.QUANT_DIRECTORY

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


class _CinfData:
    def __init__(
        self,
    ):
        self.cinfdata = None

    def activate(self):
        from cinfdata import Cinfdata

        self.DB = Cinfdata

    def connect(self, setup_name=None, grouping_column=None):
        return self.DB(setup_name=setup_name, grouping_column=grouping_column)


plugins = _PluginOptions()
