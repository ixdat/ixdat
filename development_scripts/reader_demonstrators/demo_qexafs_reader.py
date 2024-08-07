"""Sandbox script for use in the development of the qexafs reader"""

from pathlib import Path
from ixdat import Spectrum, Measurement

data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/qexafs/constant potential"

if True:
    xas = Spectrum.read(
        data_dir / "540117_IrO2_crys_0.60V_1.dat", reader="qexafs", technique="XAS"
    )
    xas.plot()

if True:
    multi_spec = Spectrum.read(
        data_dir / "540117_IrO2_crys_0.60V_1.dat",
        reader="qexafs",
    )
    ax = multi_spec["QexafsFFI0"].plot()
    ax2 = ax.twinx()
    multi_spec["I0"].plot(ax=ax2, color="r")
    ax2.set_ylim([-1e6, 2e6])

if True:
    xas_series = Spectrum.read_set(
        part=data_dir / "IrO2_crys", suffix=".dat", reader="qexafs", technique="XAS"
    )
    ax = xas_series.heat_plot()
    ax.get_figure().savefig("xas_heat_plot.png")

    ec = Measurement.read(data_dir / "IrO2_CA_0p60_C02.mpt", reader="biologic")
    ec.plot()

    ec_xas = ec + xas_series
    ec_xas.tstamp += ec_xas.t[0]

    axes = ec_xas.plot(xspan=[11200, 11300])
