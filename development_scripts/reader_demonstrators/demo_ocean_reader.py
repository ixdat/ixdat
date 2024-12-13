from pathlib import Path
from ixdat import Spectrum
import matplotlib.pyplot as plt

data_dir = Path.home() / r"Desktop\Thesis\Data\Sample 45"

my_spec = Spectrum.read(data_dir / "Absorbance__0__11-24-39-251.txt",reader="ocean")

ax = my_spec.plot()

plt.show()