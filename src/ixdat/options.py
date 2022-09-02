"""This module defines user-specific options, such as whether to use plugins"""

from pathlib import Path


class PluginOptions:
    def __init__(self):
        self._USE_QUANT = False
        self.QUANT_DIRECTORY = None

    @property
    def USE_QUANT(self):
        return self._USE_QUANT

    @USE_QUANT.setter
    def USE_QUANT(self, use_quant):
        self._USE_QUANT = use_quant
        if self._USE_QUANT:
            if not self.QUANT_DIRECTORY:
                self.QUANT_DIRECTORY = Path(__file__).parent / "plugin_data/ms_quant"

            from spectro_inlets_quantification.config import Config

            quant_config = Config()
            quant_config.data_directory = self.QUANT_DIRECTORY


plugins = PluginOptions()
