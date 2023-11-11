from pathlib import Path
from ixdat import Measurement

data_dir = (
    Path.home() / "Dropbox/ixdat_resources/test_data/biologic/17J04_Pt_isotope_exchange"
)

file = data_dir / "05_O2dose_COox_04_CA_C01.mpt"

lines = []

with open(file, "r", encoding="ISO-8859-1") as f:
    print(f)
    line = True
    while line:
        line = f.readline()
        fixed_line = line.replace("XXX", "0")
        print(fixed_line)
        lines.append(fixed_line)

fixed_file = data_dir / (file.stem + "_fixed.mpt")
with open(fixed_file, "w") as f:
    f.writelines(lines)

ax = Measurement.read_set(fixed_file, suffix=".mpt").plot()

ax.get_figure().savefig()
