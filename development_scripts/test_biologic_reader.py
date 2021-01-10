from pathlib import Path
from matplotlib import pyplot as plt
from ixdat.measurements import Measurement

test_data_dir = Path(__file__).parent.parent / "test_data/biologic_mpt_and_zilien_tsv"

path_to_test_file = test_data_dir / "2020-07-29 10_30_39 Pt_poly_cv_01_02_CVA_C01.mpt"

ec_measurement = Measurement.read(
    reader="biologic",
    path_to_file=path_to_test_file,
)

t, v = ec_measurement.get_potential(tspan=[0, 100])

ec_measurement.tstamp -= 20
t_shift, v_shift = ec_measurement.get_potential(tspan=[0, 100])

fig, ax = plt.subplots()
ax.plot(t, v, "k", label="original tstamp")
ax.plot(t_shift, v_shift, "r", label="shifted tstamp")
ax.legend()
