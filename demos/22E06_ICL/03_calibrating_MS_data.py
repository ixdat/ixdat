from ixdat import Measurement

ms_1 = Measurement.read(
    "../../data/03/2022-04-28 19_18_01 Ruti.tsv", reader="zilien", technique="MS"
)
ms_1.plot()

ms_2 = Measurement.read(
    "../../data/03/2022-04-28 23_25_58 Ruti.tsv", reader="zilien", technique="MS"
)
ms_2.plot()

ms = ms_1 + ms_2
ms.plot()

ec_1 = Measurement.read_set("../../data/03/01", reader="biologic", suffix=".mpt")
ec_1.plot()

ecms_1 = ec_1 + ms
ecms_1.plot()

ec_2 = Measurement.read_set("../../data/03/06", reader="biologic", suffix=".mpt")
ecms_2 = ec_2 + ms
ecms_2.plot()

O2_M32 = ecms_1.ecms_calibration_curve(
    mol="O2",
    mass="M32",
    n_el=4,
    tspan_list=[(1300, 1350), (1900, 1950), (2500, 2550)],
    tspan_bg=(950, 1000),
)

ecms_2.calibrate(ms_cal_results=[O2_M32], RE_vs_RHE=0, A_el=0.196)
ecms_2.plot(mol_list=["O2"], tspan=[0, 300], logplot=False, unit="nmol/s")
