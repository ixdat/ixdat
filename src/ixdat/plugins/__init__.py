from .siq_plugin import SIQ_Plugin
from .cinfdata_plugin import CinfData_Plugin
from ..tools import deprecate


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
        classes and the siq `Calculator` are then available as `siq.<class_name>`
    - deactivate_si_quant(): Sets use_si_quant to False.
    - si_quant (_QuantDeps). This is where the needed imports from the external
        quantification package go. See the docstring of _QuantDeps

    cinfdata
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    - activate_cinfdata(): sets up a database as `cinfdata` The database has to be
        connected using `cinfdata.connect(...)`
    """

    def __init__(self):
        self._use_siq = False
        self.siq = SIQ_Plugin()
        self.cinfdata = CinfData_Plugin()

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
        hard_deprecation_release="0.3.1",
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
        hard_deprecation_release="0.3.1",
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
        hard_deprecation_release="0.3.1",
        remove_release="1.0.0",
    )
    def si_quant(self):
        return self.siq


plugins = _PluginOptions()
