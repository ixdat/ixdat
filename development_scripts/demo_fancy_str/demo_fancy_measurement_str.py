"""For use in development of the zilien reader(s). Requires access to sample data."""

from pathlib import Path

from ixdat import Measurement
from ixdat.techniques import MSMeasurement, ECMeasurement

data_dir = Path(
    r"/home/kenneth/Spectro Inlets Dropbox/Kenneth Nielsen/ixdat_resources/test_data/zilien_with_ec"
)

path_to_file = data_dir / "2021-02-01 17_44_12.tsv"

# This imports it with the EC data
ecms = Measurement.read(path_to_file, reader="zilien")
voltage = ecms["Voltage [V]"]
print()
print("Voltage __repr__:", repr(voltage))
print("Voltage __str__:", voltage)
print()
voltage_time = voltage_time = voltage.tseries
print("Voltage time __repr__:", repr(voltage_time))
print("Voltage time __str__:", voltage_time)
print()
print("ECMS repr", repr(ecms))
print("The __str__ is:")
print(ecms)
