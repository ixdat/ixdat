from ixdat import plugins

plugins.activate_si_quant()


from pathlib import Path


data_dir = Path.home() / "Dropbox/ixdat_resources/presentations_and_workshops/22K16_ICL_quant/2022-09-27 21_15_42 22I27_London"

for f in data_dir.iterdir():
    print(f)


from ixdat import Measurement

ms = Measurement.read_set(data_dir, suffix=".tsv", technique="MS")
# ms.plot()
ec = Measurement.read_set(data_dir / "03", suffix=".mpt")
# ec.plot()


from ixdat.readers.biologic import fix_WE_potential

fix_WE_potential(ec)

# ec.plot()
ecms = ec + ms
ecms.plot(J_name="selector", tspan=[900, 1600])

cal_H2_M2_EC = ecms.ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan_bg=[50, 100],
    selector_list=[2, 5, 8],
    t_steady_pulse=30
)
cal_O2_M32_EC = ecms.ecms_calibration_curve(
    mol="O2",
    mass="M32",
    n_el=4,  # 2 H2O --> O2 + 4 (H+ + e-)
    tspan_bg=[50, 100],
    selector_list=[10, 12, 14],
    t_steady_pulse=30
)
cal_CO2_M44_EC = ecms.ecms_calibration_curve(
    mol="CO2",
    mass="M44",
    n_el=2,
    tspan_bg=[2520, 2540],
    selector_list=[28, 32, 34],
    t_steady_pulse=30
)

cal_He = ecms.gas_flux_calibration(
    mol="He",
    mass="M4",
    tspan=[0, 50]
)
cal_CO = ecms.gas_flux_calibration(
    mol="CO",
    mass="M28",
    tspan=[2800, 2850]
)

cal_air = ms.multicomp_gas_flux_calibration(
    mol_list=["N2", "O2", "Ar"],
    mass_list=["M28", "M32", "M40"],
    gas={"N2": 0.79, "O2": 0.20, "Ar": 0.01},
    tspan=[7100, 7200]
)

from ixdat.si_quant_patch import append_sensitivity_factors


calibration = append_sensitivity_factors(
    cal_H2_M2_EC, cal_O2_M32_EC, cal_CO2_M44_EC, cal_He, cal_CO, cal_air
)
ax = calibration.plot_as_spectrum()
ax.get_figure().savefig("calibration_as_spectrum.png")
ax = calibration.plot_F_vs_f()
ax.get_figure().savefig("calibration_trend.png")


ecms.set_quantifier(
    calibration=calibration, mol_list=["CO2"], mass_list=["M44"]
)
ecms.export(mol_list=["CO2"], mass_list=["M44"])