"""initialize ixdat, giving top-level access to a few of the important structures
"""
__version__ = "0.1.6"
__title__ = "ixdat"
__description__ = "The in-situ experimental data tool"
__url__ = "https://github.com/ixdat/ixdat"
__author__ = "Soren B. Scott, Kevin Krempl, Kenneth Nielsen"
__email__ = "scott.soren@gmail.com"  # maybe we should get an orgianization email?
# __copyright__ = "Copyright (c) 2020 ixdat"
__license__ = "MIT"

from .measurements import Measurement
from .spectra import Spectrum
from . import db
from . import techniques
from . import plotters
from . import exporters

# I like this to be sure I'm importing from where I think I am:
print(f"importing ixdat v{__version__} from {__file__}")
