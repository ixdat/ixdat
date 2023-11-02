from ixdat import Measurement

ms = Measurement.read(
    "../data/03/2022-04-28 19_18_01 Ruti.tsv", reader="zilien", technique="MS"
) + Measurement.read(
    "../data/03/2022-04-28 23_25_58 Ruti.tsv", reader="zilien", technique="MS"
)

ec = Measurement.read_set("../data/03/01", reader="biologic", suffix=".mpt")
ecms = ec + ms

axes = ecms.plot()

O2_M32, ax_cal = ecms.ecms_calibration_curve(
    mol="O2",
    mass="M32",
    n_el=4,
    tspan_list=[(1300, 1350), (1900, 1950), (2500, 2550)],
    tspan_bg=(950, 1000),
    return_ax=True,
    axes_measurement=axes
)
ecms.calibrate(ms_cal_results=[O2_M32], RE_vs_RHE=0, A_el=0.196)
