from pathlib import Path
from ixdat import Spectrum

data_dir = Path(
    r"C:\Users\scott\Dropbox\WORKSPACES\Caiwu\22E23_beamline_data\constant potential"
)

spec = Spectrum.read(data_dir / "540117_IrO2_crys_0.60V_1.dat", reader="qexafs")
print(spec)
