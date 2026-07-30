.. _plotters:

Plotters: visualizing ``ixdat`` data
====================================
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/plotters

In ``ixdat``, straightforward plotting of your data is a priority. This page documents
plotters for the supported experimental data types. Use the menu on the right to move
between sections.

Plotter backends
----------------

Measurements and spectra use Matplotlib by default. Install ixdat's Plotly extra to
create interactive figures::

    pip install "ixdat[plotly]"

Select Plotly when calling a plot method::

    from ixdat import Measurement

    measurement = Measurement.read(
        "experiment.mpt",
        reader="biologic",
    )
    figure = measurement.plot(backend="plotly")
    figure.show()

Plot methods with a registered adapter accept the ``backend`` argument at plot time,
so one measurement can produce Matplotlib and Plotly figures. Plotly support covers
the default plots for generic, EC, MS, EC-MS, spectro-MS, spectroelectrochemistry,
TP-MS, spectro-TP-MS, and time-resolved XRF measurements. It also covers spectra,
NMR spectra, XRD spectra, and spectrum-series heatmap, waterfall, and stacked plots.
Series names such as ``<I>`` appear as literal text in Plotly legends and axis titles.
Multi-panel time plots link their x-axes and show the time label on the bottom panel.
Hovering at one time reports values from traces in every panel. Plotly uses a 640 by
480 canvas with the default Matplotlib settings and follows configured Matplotlib
figure dimensions. It frames each panel with a black border, omits grid lines,
uses outward ticks to space tick labels from the border, reserves margins for axis
titles, and places legends above the data. A pair of left and right y-axis traces
uses matching colors for each trace, axis title, ticks, and tick labels. The colored
axes identify those traces, and the paired traces omit legend entries. Spectro-TP-MS
mass traces use the MS palette, and stacked Plotly spectra use a different color for
each trace. TP-MS temperature and pressure axes include the series units. The
``TP_units`` argument on TP-MS plots and ``meta_units`` on spectro-TP-MS plots update
the displayed unit and apply the corresponding scale factor. A requested time unit
applies to every panel in a composite plot. Two-panel figures follow the same panel
proportions as Matplotlib. Spectrum-series heatmaps use ``continuous`` and recorded
spectrum durations to place their cells. Common ``color``, ``linestyle``, and
``linewidth`` arguments carry into Plotly line traces where the plot method accepts
them.

Pass a plain Plotly figure through ``figure`` to add the subplot layout and traces::

    import plotly.graph_objects as go

    figure = go.Figure(layout={"width": 800, "height": 500})
    measurement.plot(backend="plotly", figure=figure)

A figure returned by the same ixdat plot method has a compatible subplot layout and
can receive more traces. ixdat reports an incompatible subplot grid with a
``ValueError``.

For Plotly and registered extension backends, an adapter describes the plotter call
with a :class:`~ixdat.plotters.plot_spec.PlotSpec`. A plot specification contains
panels, axes, line traces, error bands, heatmaps, and continuous color scales. The
selected renderer translates those components into a plotting library's objects.

Readers construct the appropriate ixdat data type and assign its Matplotlib plotter.
A new reader gains the adapters registered for that plotter class. A new plotter works
through Matplotlib without any adapter. If a user requests Plotly for that plotter,
ixdat issues a ``PlotterBackendWarning`` and returns the Matplotlib result. Plotly
support can be added later by registering an adapter outside the Matplotlib plotter
class.

ixdat ships Matplotlib and Plotly renderers. Register another renderer class with
``register_plotter_backend(name, renderer_class)``. The class implements
``render(plot_spec, **kwargs)``. This single interface supports libraries such as
Seaborn wherever the plot specification uses components that the renderer handles.
``available_plotter_backends()`` returns registered names. Register plot coverage
with ``register_plotter_adapter(plotter_class, method_name, adapter)``. The adapter
receives the plotted data object and a mapping of arguments bound to the original
Matplotlib method, then returns a ``PlotSpec``.

The example ``development_scripts/demo_plotly_backend.py`` compares Matplotlib and
Plotly for every supported Plotly path. It uses repository EC, NMR, and XRD data and
dense synthetic data for the remaining techniques. Script execution writes
``development_scripts/demo_plotly_backend.html`` and opens the complete gallery in the
browser. The short ``development_scripts/demo_plotly_backend.ipynb`` notebook shows
the two backends side by side.

Basic
-----

The ``backends`` module
.......................

.. automodule:: ixdat.plotters.backends
    :members:

The ``plot_spec`` module
........................

.. automodule:: ixdat.plotters.plot_spec
    :members:

The ``renderers`` module
........................

.. automodule:: ixdat.plotters.renderers
    :members:

The ``base_mpl_plotter`` module
...............................

.. automodule:: ixdat.plotters.base_mpl_plotter
    :members:

The ``value_plotter`` module
............................

.. automodule:: ixdat.plotters.value_plotter
    :members:

Electrochemistry
----------------

.. _`ec-plotter`:

The ``ec_plotter`` module
.........................

.. automodule:: ixdat.plotters.ec_plotter
    :members:

Mass Spectrometry
-----------------

.. _`ms-plotter`:

The ``ms_plotter`` module
.........................

.. automodule:: ixdat.plotters.ms_plotter
    :members:

EC-MS
-----

The ``ecms_plotter`` module
...........................

.. _ecms-plotter:

.. automodule:: ixdat.plotters.ecms_plotter
    :members:

Spectra
-------

The ``spectrum_plotter`` module
...............................

.. _spectrum-plotter:

.. automodule:: ixdat.plotters.spectrum_plotter
    :members:

Spectroelectrochemistry
-----------------------

.. _sec-plotter:

The ``sec_plotter`` module
..........................

.. automodule:: ixdat.plotters.sec_plotter
    :members:
