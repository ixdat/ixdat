"""For use in development of the MSRH SEC reader. Requires access to sample data.
MSRH = molecular science research hub, at Imperial College London.
"""

from pathlib import Path
from ixdat import Measurement

from matplotlib import pyplot as plt

plt.close("all")

data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/sec"
sec_meas = Measurement.read(
    data_dir / "test-7SEC.csv",
    path_to_ref_spec_file=data_dir / "WL.csv",
    path_to_U_J_file=data_dir / "test-7_JV.csv",
    scan_rate=1,
    tstamp="now",
    reader="msrh_sec",
)

sec_meas.calibrate_RE(RE_vs_RHE=0.26)  # provide RE potential in [V] vs RHE
sec_meas.normalize_current(A_el=1)  # provide electrode area in [cm^2]

# This plotting makes use of a built-in DeltaODCalculator:
sec_meas.set_reference_spectrum(V_ref=0.4)

axes = sec_meas.plot_measurement(
    cmap_name="jet",
    make_colorbar=True,
)
ax = sec_meas.plot_waterfall()
ax.get_figure().savefig("sec_waterfall.png")

if False:  # This makes use of a built-in DeltaODCalculator

    spec1 = sec_meas.calc_dOD_spectrum(V=1.0, V_ref=0.66)
    spec2 = sec_meas.calc_dOD_spectrum(V=1.4, V_ref=1.0)
    spec3 = sec_meas.calc_dOD_spectrum(V=1.7, V_ref=1.4)

else:  # another way to do the above, explicitly using the calculator instead
    from ixdat.calculators import DeltaODCalculator


    calc = DeltaODCalculator(measurement=sec_meas)

    spec1 = calc.spectrum(V=1.0, V_ref=0.66)
    spec2 = calc.spectrum(V=1.4, V_ref=1.0)
    spec3 = calc.spectrum(V=1.7, V_ref=1.4)

ax = spec1.plot(color="b", label="redox 1")
spec2.plot(color="g", label="redox 2", ax=ax)
spec3.plot(color="r", label="redox 3", ax=ax)
ax.legend()

if True:  # test export and reload
    # Suggestion: command-line switching for development scripts.
    #  https://github.com/ixdat/ixdat/pull/30/files#r810014299
    export_name = "exported_sec.csv"
    sec_meas.export(export_name)
    sec_reloaded = Measurement.read(export_name, reader="ixdat")
    sec_reloaded.set_reference_spectrum(V_ref=0.66)
    sec_reloaded.plot_vs_potential(cmap_name="jet")
    sec_reloaded.continuous = False
    sec_reloaded.plot_vs_potential(cmap_name="jet")



ax = sec_meas.plot_waterfall(
    V_ref=0.4,
    cmap_name="jet",
    make_colorbar=True,
)

axes2 = sec_meas.plot_vs_potential(V_ref=0.66, cmap_name="jet", make_colorbar=False)
axes2 = sec_meas.plot_vs_potential(
    V_ref=0.66, vspan=[1.0, 1.5], wlspan=[500, 700], cmap_name="jet", make_colorbar=False
)


from ixdat.calculators import Surfer
# a Surfer is a calculator that makes a ValueSeries by following a specific value of x
#  in a SpectrumSeries.


sec.add_calculator(Surfer(xs={"w450": 460, "w600": 600, "w850": 850}))
axes = sec_meas.plot_xs_vs_potential(["w460", "w600", "w850"])
axes[0].set_ylabel("intense!")

ax.legend()
