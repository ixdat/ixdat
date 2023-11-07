from pathlib import Path
from ixdat import Measurement


data_dir = (
    Path.home() / "Dropbox/ixdat_resources/test_data/biologic/17J04_Pt_isotope_exchange"
)
data_file = data_dir / "05_O2dose_COox_01_LSV_C01.mpr"

import eclabfiles

# data is a list of dictionaries, one for each line of the file:
# data, meta = eclabfiles.process(data_file)

df = eclabfiles.to_df(data_file)


from ixdat import Measurement

meas = Measurement.read_set(
    data_dir / "05_O2dose_COox_03_LSV_C01.mpr", reader="biologic"
)
# unfortunately the CVA and CA files in that same folder do not work.

meas.plot(J_str="loop number")
