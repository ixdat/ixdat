"""For use in development of the cinfdata reader. Requires access to sample data."""

from pathlib import Path

from ixdat import Measurement

path_to_file = (
    Path.home()
    / "Dropbox/ixdat_resources/test_data/cinfdata/Trimarco2018_fig3/QMS_1.txt"
)
ms_meas = Measurement.read(path_to_file, reader="cinfdata")
ms_meas.plot_measurement()

path_to_ec_file_start = (
    Path.home() / "Dropbox/ixdat_resources/test_data/cinfdata/Trimarco2018_fig3/09_fig4"
)
ec_meas = Measurement.read_set(path_to_ec_file_start, reader="biologic")
ec_meas.calibrate(RE_vs_RHE=0.65, A_el=0.196)
ec_meas.plot_measurement()

ecms_meas = ec_meas + ms_meas
axes = ecms_meas.plot_measurement(
    mass_lists=[["M4", "M28"], ["M44", "M2"]],
    tspan_bg=[None, [30, 40]],
    legend=False,
    unit="pA",
)

axes[2].set_ylim([-7, 70])
axes[0].set_ylim([-1.8e3, 18e3])
fig = axes[0].get_figure()
fig.tight_layout()
# fig.savefig("../../docs/source/figures/ec_ms.svg")

ecms_meas.set_bg(tspan_bg=[0, 10])

cv = ecms_meas.as_cv()
cv.redefine_cycle(start_potential=0.39, redox=False)

axes_cv = cv[2].plot(
    mass_list=["M2", "M44"],
    logplot=False,
)
axes_cv = cv[1].plot(
    mass_list=["M2", "M44"],
    linestyle="--",
    axes=axes_cv,
    logplot=False,
)
axes_cv[0].get_figure().savefig("Trimarco2018_ixdat.png")
