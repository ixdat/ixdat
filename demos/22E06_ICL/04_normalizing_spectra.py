import numpy as np
from matplotlib import pyplot as plt
from ixdat import Spectrum


to_plot = {
    "Ru": ("RuO$_2$", "Rude4_omega0p5.xrdml", "black"),
    "RuTi_1": ("Ru$_{0.9}$Ti$_{0.1}$O$_2$", "RutiZelensky_omega0p5.xrdml", "blue"),
    "RuTi_2": ("Ru$_{0.75}$Ti$_{0.25}$O$_2$", "RutiMacron_omega0p5.xrdml", "green"),
    "Ti": ("TiO$_2$", "Poseidon_omega0p5.xrdml", "red"),
}

spectra = {}

for name, (formula, file_name, color) in to_plot.items():
    xrd = Spectrum.read("../../data/04/" + file_name, reader="xrdml")
    ax = xrd.plot(color=color)
    ax.set_title(formula)
    spectra[name] = xrd

fig, ax = plt.subplots()
for name, (formula, file_name, color) in to_plot.items():
    xrd = spectra[name]
    x, y = xrd.x, xrd.y
    y_fto_peak = max(y[np.logical_and(37 < x, x < 38.5)])
    ax.plot(x, y / y_fto_peak, label=formula, color=color)

ax.set_xlabel("two theta / [deg]")
ax.set_ylabel("norm. intensity")
ax.legend()
