"""Compare ixdat's Matplotlib and Plotly backends across supported data types.

The EC, NMR, and XRD examples use repository data. The remaining examples use
dense synthetic signals. Running this file opens a gallery containing every
Matplotlib and Plotly pair.

Install the Plotly extra before running the demo:

    pip install ".[plotly]"
"""

# %%
from base64 import b64encode
from collections import OrderedDict
from html import escape
from io import BytesIO
from pathlib import Path
import webbrowser
import warnings

from matplotlib import MatplotlibDeprecationWarning
import numpy as np
import plotly.io as pio

from ixdat import Measurement, Spectrum
from ixdat.data_series import DataSeries, Field, TimeSeries, ValueSeries
from ixdat.spectra import SpectrumSeries
from ixdat.techniques.ec import ECMeasurement
from ixdat.techniques.ec_ms import ECMSMeasurement
from ixdat.techniques.ms import MSMeasurement, MSSpectroMeasurement
from ixdat.techniques.nmr import NMRSpectrum
from ixdat.techniques.reactor import ReactorMeasurement, ReactorSpectroMeasurement
from ixdat.techniques.spectroelectrochemistry import SpectroECMeasurement
from ixdat.techniques.xrd import XRDSpectrum
from ixdat.techniques.xrf import ECTRXRFMeasurement, TRXRFMeasurement


warnings.filterwarnings(
    "ignore",
    message=r"Temperature is not factorial converted.*",
    category=UserWarning,
)
warnings.filterwarnings(
    "ignore",
    message=r"Can't convert original unit.*",
    category=UserWarning,
)
warnings.filterwarnings(
    "ignore",
    category=MatplotlibDeprecationWarning,
    module=r"ixdat\..*",
)
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    module=r"nmrglue\..*",
)
warnings.filterwarnings(
    "ignore",
    category=ResourceWarning,
    module=r"subprocess",
)
pio.renderers.default = "plotly_mimetype"

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEMO_TIME = np.linspace(0, 240, 241)
PHASE = 2 * np.pi * DEMO_TIME / DEMO_TIME[-1]
EC_ALIASES = {
    "t": ["time"],
    "raw_potential": ["Ewe/V"],
    "raw_current": ["<I>/mA"],
}
MS_VALUES = [
    (
        "M2",
        "A",
        0.3e-9
        + 2.5e-9 * np.exp(-(((DEMO_TIME - 75) / 24) ** 2))
        + 1.1e-9 * np.exp(-(((DEMO_TIME - 175) / 35) ** 2)),
    ),
    (
        "M32",
        "A",
        0.2e-9
        + 1.4e-9 * np.exp(-(((DEMO_TIME - 125) / 30) ** 2))
        + 0.08e-9 * np.sin(PHASE * 3) ** 2,
    ),
]
REACTOR_VALUES = [
    ("temperature", "K", 315 + 85 * DEMO_TIME / DEMO_TIME[-1]),
    ("pressure", "bar", 1.02 + 0.06 * np.sin(PHASE - 0.3)),
]
XRF_VALUES = [
    (
        "FF_over_I0",
        "",
        0.25 + 0.1 * np.exp(-(((DEMO_TIME - 90) / 42) ** 2)) + 0.025 * np.sin(PHASE * 2),
    )
]


def make_measurement(
    measurement_class,
    name,
    values,
    time=None,
    aliases=None,
    spectrum_series=None,
):
    """Build one measurement from aligned arrays."""
    time_data = np.asarray(DEMO_TIME if time is None else time)
    time_series = TimeSeries(
        name="time",
        unit_name="s",
        data=time_data,
        tstamp=0,
    )
    series_list = [time_series]
    for value_name, unit_name, data in values:
        series_list.append(
            ValueSeries(
                name=value_name,
                unit_name=unit_name,
                data=np.asarray(data),
                tseries=time_series,
            )
        )
    kwargs = {
        "name": name,
        "series_list": series_list,
        "aliases": aliases,
        "tstamp": 0,
    }
    if spectrum_series is not None:
        kwargs["spectrum_series"] = spectrum_series
    return measurement_class(**kwargs)


def make_spectrum_series(name="spectrum series", times=None):
    """Build a smooth spectrum series with moving and growing peaks."""
    spectrum_times = np.asarray(
        np.linspace(0, DEMO_TIME[-1], 21) if times is None else times
    )
    wavelength = np.linspace(380, 760, 160)
    spectra = []
    for index, spectrum_time in enumerate(spectrum_times):
        fraction = index / (len(spectrum_times) - 1)
        main_peak = (0.7 + 0.8 * fraction) * np.exp(
            -(((wavelength - 515 - 25 * fraction) / 35) ** 2)
        )
        shoulder = (0.35 + 0.15 * np.sin(fraction * np.pi)) * np.exp(
            -(((wavelength - 650 + 10 * fraction) / 55) ** 2)
        )
        baseline = 0.08 + 0.015 * np.sin(wavelength / 25 + fraction * np.pi)
        spectra.append(
            Spectrum.from_data(
                x=wavelength,
                y=baseline + main_peak + shoulder,
                x_name="wavelength / [nm]",
                y_name="intensity / [a.u.]",
                name=name,
                tstamp=float(spectrum_time),
            )
        )
    return SpectrumSeries.from_spectrum_list(spectra)


def make_nmr_spectrum():
    """Load the repository NMR data and reduce its display size."""
    data_path = REPOSITORY_ROOT / "test_data/bruker/MTBLS1_ADG19007u_162_10"
    try:
        source = Spectrum.read(data_path, reader="bruker")
        step = max(1, len(source.x) // 2000)
        return NMRSpectrum.from_data(
            x=source.x[::step],
            y=source.y[::step],
            x_name="chemical shift / [ppm]",
            y_name="intensity / [a.u.]",
            name=source.name,
            tstamp=source.tstamp,
        )
    except ImportError:
        chemical_shift = np.linspace(-2, 12, 1800)
        intensity = (
            0.04
            + 1.5 * np.exp(-(((chemical_shift - 4.7) / 0.06) ** 2))
            + 0.8 * np.exp(-(((chemical_shift - 3.25) / 0.12) ** 2))
            + 0.45 * np.exp(-(((chemical_shift - 1.3) / 0.09) ** 2))
        )
        return NMRSpectrum.from_data(
            x=chemical_shift,
            y=intensity,
            x_name="chemical shift / [ppm]",
            y_name="intensity / [a.u.]",
            name="NMR spectrum",
            tstamp=0,
        )


def make_xrd_spectrum():
    """Load repository XRD values and interpolate them for display."""
    source_path = REPOSITORY_ROOT / "test_data/xrd/no_header.xye"
    source = Spectrum.read(source_path, reader="xrdxy")
    x = np.linspace(source.x.min(), source.x.max(), 180)
    x_series = DataSeries(
        name=source.x_name,
        unit_name=source.xseries.unit_name,
        data=x,
    )
    fields = [
        Field(
            name=source.y_name,
            unit_name=source.fields[0].unit_name,
            data=np.interp(x, source.x, source.y),
            axes_series=[x_series],
        ),
        Field(
            name="intensity_error",
            unit_name=source.fields[1].unit_name,
            data=np.interp(x, source.x, source.y_err),
            axes_series=[x_series],
        ),
    ]
    return XRDSpectrum(
        name=source.name,
        technique=source.technique,
        tstamp=source.tstamp,
        fields=fields,
    )


def get_matplotlib_figure(plot_result):
    """Return the Matplotlib figure from an axis or axis collection."""
    if hasattr(plot_result, "figure"):
        return plot_result.figure
    for axis in plot_result:
        if axis is not None:
            return axis.figure
    raise TypeError("A Matplotlib plot must return an axis or axis collection.")


def add_comparison(comparisons, title, matplotlib_result, plotly_figure):
    """Add Matplotlib and Plotly versions of one plot."""
    plotly_figure.update_layout(template="plotly_white")
    comparisons[title] = (
        get_matplotlib_figure(matplotlib_result),
        plotly_figure,
    )


def matplotlib_data_uri(figure):
    """Encode one Matplotlib figure for the comparison gallery."""
    image = BytesIO()
    figure.savefig(image, format="png", dpi=125, bbox_inches="tight")
    return "data:image/png;base64," + b64encode(image.getvalue()).decode("ascii")


def open_gallery(comparisons):
    """Write the backend comparisons beside this script and open them."""
    sections = []
    for index, (title, (matplotlib_figure, plotly_figure)) in enumerate(
        comparisons.items()
    ):
        matplotlib_image = matplotlib_data_uri(matplotlib_figure)
        plotly_html = plotly_figure.to_html(
            full_html=False,
            include_plotlyjs=True if index == 0 else False,
            config={"responsive": True},
        )
        sections.append(
            """
<section>
  <h2>{title}</h2>
  <div class="comparison">
    <div class="backend">
      <h3>Matplotlib</h3>
      <img src="{matplotlib_image}" alt="Matplotlib {title}">
    </div>
    <div class="backend">
      <h3>Plotly</h3>
      {plotly_html}
    </div>
  </div>
</section>
""".format(
                title=escape(title),
                matplotlib_image=matplotlib_image,
                plotly_html=plotly_html,
            )
        )
    document = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ixdat plot backends</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem auto; max-width: 1800px; }}
    h1, h2 {{ margin-left: 1rem; }}
    section {{ border-top: 1px solid #ddd; margin-top: 2rem; padding-top: 1rem; }}
    .comparison {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }}
    .backend {{ min-width: 0; padding: 0 1rem; }}
    .backend h3 {{ text-align: center; }}
    .backend img {{
      display: block;
      height: auto;
      margin: auto;
      max-width: 100%;
      width: 640px;
    }}
    .backend .plotly-graph-div {{
      margin: auto;
      max-width: 640px;
      width: 100% !important;
    }}
    @media (max-width: 900px) {{
      .comparison {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <h1>ixdat plot backends</h1>
  <p>Each row uses the same ixdat data and plot method.</p>
  {sections}
</body>
</html>
""".format(
        sections="\n".join(sections)
    )
    output_path = Path(__file__).resolve().with_suffix(".html")
    output_path.write_text(document, encoding="utf-8")
    webbrowser.open(output_path.as_uri())
    return output_path


# %%
comparisons = OrderedDict()

real_ec = Measurement.read(
    REPOSITORY_ROOT / "test_data/biologic/Pt_poly_cv_CUT.mpt",
    reader="biologic",
)
ec_time, ec_potential = real_ec.grab(real_ec.U_name)
ec_current = real_ec.grab_for_t(real_ec.J_name, ec_time)
real_ec_values = [
    ("Ewe/V", "V", ec_potential),
    ("<I>/mA", "mA", ec_current),
]

generic = make_measurement(
    Measurement,
    "generic",
    [("potential / [V]", "V", ec_potential)],
    time=ec_time,
)
add_comparison(
    comparisons,
    "Generic Measurement",
    generic.plot(),
    generic.plot(backend="plotly"),
)

add_comparison(
    comparisons,
    "EC measurement",
    real_ec.plot(),
    real_ec.plot(backend="plotly"),
)
add_comparison(
    comparisons,
    "EC current vs potential",
    real_ec.plot_vs_potential(),
    real_ec.plot_vs_potential(backend="plotly"),
)

ms = make_measurement(
    MSMeasurement,
    "MS",
    MS_VALUES,
)
add_comparison(
    comparisons,
    "MS measurement",
    ms.plot(logplot=False),
    ms.plot(backend="plotly", logplot=False),
)

ecms_time_scale = (ec_time - ec_time.min()) / (ec_time.max() - ec_time.min())
ecms_values = real_ec_values + [
    (
        "M2",
        "A",
        0.25e-9 + 2.2e-9 * np.exp(-(((ecms_time_scale - 0.42) / 0.14) ** 2)),
    ),
    (
        "M32",
        "A",
        0.18e-9 + 1.3e-9 * np.exp(-(((ecms_time_scale - 0.68) / 0.18) ** 2)),
    ),
]
ecms = make_measurement(
    ECMSMeasurement,
    "EC-MS",
    ecms_values,
    time=ec_time,
    aliases=EC_ALIASES,
)
add_comparison(
    comparisons,
    "EC-MS measurement",
    ecms.plot(logplot=False),
    ecms.plot(backend="plotly", logplot=False),
)

# %%
ms_spectro = make_measurement(
    MSSpectroMeasurement,
    "spectro-MS",
    MS_VALUES,
    spectrum_series=make_spectrum_series("MS spectra"),
)
add_comparison(
    comparisons,
    "Spectro-MS measurement",
    ms_spectro.plot(logplot=False),
    ms_spectro.plot(backend="plotly", logplot=False),
)

sec = make_measurement(
    SpectroECMeasurement,
    "spectroelectrochemistry",
    real_ec_values,
    time=ec_time,
    aliases=EC_ALIASES,
    spectrum_series=make_spectrum_series(
        "optical spectra",
        times=np.linspace(ec_time.min(), ec_time.max(), 21),
    ),
)
add_comparison(
    comparisons,
    "Spectroelectrochemistry measurement",
    sec.plot(),
    sec.plot(backend="plotly"),
)

tpms = make_measurement(
    ReactorMeasurement,
    "TP-MS",
    MS_VALUES + REACTOR_VALUES,
)
add_comparison(
    comparisons,
    "TP-MS measurement",
    tpms.plot(logplot=False, x_unit="s"),
    tpms.plot(backend="plotly", logplot=False, x_unit="s"),
)

tpms_spectro = make_measurement(
    ReactorSpectroMeasurement,
    "spectro-TP-MS",
    MS_VALUES + REACTOR_VALUES,
    spectrum_series=make_spectrum_series("reactor spectra"),
)
add_comparison(
    comparisons,
    "Spectro-TP-MS measurement",
    tpms_spectro.plot(logplot=False, x_unit="s"),
    tpms_spectro.plot(backend="plotly", logplot=False, x_unit="s"),
)

# %%
xrf = make_measurement(
    TRXRFMeasurement,
    "TR-XRF",
    XRF_VALUES,
)
add_comparison(
    comparisons,
    "Time-resolved XRF measurement",
    xrf.plot(),
    xrf.plot(backend="plotly"),
)

ec_xrf = make_measurement(
    ECTRXRFMeasurement,
    "EC-TR-XRF",
    real_ec_values
    + [
        (
            "FF_over_I0",
            "",
            0.24 + 0.11 * np.exp(-(((ecms_time_scale - 0.55) / 0.2) ** 2)),
        )
    ],
    time=ec_time,
    aliases=EC_ALIASES,
)
add_comparison(
    comparisons,
    "EC and time-resolved XRF measurement",
    ec_xrf.plot(),
    ec_xrf.plot(backend="plotly"),
)

# %%
nmr = make_nmr_spectrum()
generic_spectrum = Spectrum.from_data(
    x=nmr.x,
    y=nmr.y,
    x_name=nmr.x_name,
    y_name=nmr.y_name,
    name="spectrum",
    tstamp=nmr.tstamp,
)
add_comparison(
    comparisons,
    "Spectrum",
    generic_spectrum.plot(),
    generic_spectrum.plot(backend="plotly"),
)
add_comparison(
    comparisons,
    "NMR spectrum",
    nmr.plot(),
    nmr.plot(backend="plotly"),
)

xrd = make_xrd_spectrum()
add_comparison(
    comparisons,
    "XRD reader data",
    xrd.plot(),
    xrd.plot(backend="plotly"),
)

# %%
spectrum_series = make_spectrum_series()
add_comparison(
    comparisons,
    "SpectrumSeries heatmap",
    spectrum_series.plot(make_colorbar=True),
    spectrum_series.plot(backend="plotly", make_colorbar=True),
)
add_comparison(
    comparisons,
    "SpectrumSeries waterfall",
    spectrum_series.plotter.plot_waterfall(),
    spectrum_series.plotter.plot_waterfall(backend="plotly"),
)
add_comparison(
    comparisons,
    "SpectrumSeries stacked spectra",
    spectrum_series.plotter.plot_stacked_spectra(index_list=[0, 5, 10, 15, 20]),
    spectrum_series.plotter.plot_stacked_spectra(
        backend="plotly",
        index_list=[0, 5, 10, 15, 20],
    ),
)

# %%
comparisons["EC-MS measurement"][1]

# %%
comparisons["Spectro-TP-MS measurement"][1]

# %%
comparisons["SpectrumSeries heatmap"][1]


if __name__ == "__main__":
    gallery_path = open_gallery(comparisons)
    print(f"Opened {gallery_path}")
