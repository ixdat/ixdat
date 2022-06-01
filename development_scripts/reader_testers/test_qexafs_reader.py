from pathlib import Path
from ixdat import Spectrum, Measurement

data_dir = Path(
    r"C:\Users\scott\Dropbox\WORKSPACES\Caiwu\22E23_beamline_data\constant potential"
)

if True:
    xas = Spectrum.read(
        data_dir / "540117_IrO2_crys_0.60V_1.dat", reader="qexafs", technique="xas"
    )
    xas.plot()

if False:
    multi_spec = Spectrum.read(
        data_dir / "540117_IrO2_crys_0.60V_1.dat",
        reader="qexafs",
    )
    for spectrum in multi_spec.spectrum_list:
        spectrum.plot()

if True:
    xas_series = Spectrum.read_set(
        part=data_dir / "IrO2_crys", suffix=".dat", reader="qexafs", technique="xas"
    )
    ax = xas_series.heat_plot()
    ax.get_figure().savefig("xas_heat_plot.png")

    ec = Measurement.read(data_dir / "IrO2_CA_0p60_C02.mpt", reader="biologic")
    ec.plot()

    ec_xas = ec + xas_series

    ec_xas.plot()
