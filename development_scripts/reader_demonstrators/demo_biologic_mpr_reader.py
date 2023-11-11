from pathlib import Path
from ixdat import Measurement


data_dir = (
    Path.home() / "Dropbox/ixdat_resources/test_data/biologic/17J04_Pt_isotope_exchange"
)

import eclabfiles


combined_meas = None

for file in data_dir.iterdir():
    if not file.suffix == ".mpr":
        continue
    try:
        meas = Measurement.read(file, reader="biologic")
    except Exception as e:
        print(f"{file.name} did not work. Error: {e}")
        try:
            df = eclabfiles.to_df(file)
        except Exception as e1:
            print(f"\t `eclabfiles.to_df` gives {e1}")
        try:
            data, meta = eclabfiles.process(file)
        except Exception as e2:
            print(f"`\t eclabfiles.process` gives {e2}")
    else:
        meas.plot()
        print(f"{file.name} worked :) ")
        if combined_meas:
            combined_meas = combined_meas + meas
        else:
            combined_meas = meas
    # unfortunately the CVA and CA files in that same folder do not work.

combined_meas.plot(J_name="loop number")
