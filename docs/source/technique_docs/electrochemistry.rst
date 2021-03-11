.. _electrochemistry:

Electrochemistry
================

The main TechniqueMeasurement class for electrochemistry is the `ECMeasurement`.
Sublcasses of `ECMeasurement` include `CyclicVoltammagram` and `CyclicVoltammagramDiff`.

Direct-current electrochemsitry measurements (`ixdat` does not yet offer specific
functionality for impedance data) are characterized by the essential quantities being
working-electrode current (in loop with the counter electrode) and potential (vs the
reference electrode) as a function of time. Either current or potential can be controlled
as the input variable, so the other acts at the response, and it is common to plot
current vs potential, but in all cases both are tracked or controlled as a function of
time. This results in the essential variables `t` (time), `v` (potential), and `j`
(current). The main job of `ECMeasurement` and subclasses is to give standardized,
convenient, and powerful access to these three variables for data selection, analysis,
and visualization, regardless of which hardware the data was acquired with.

The default plotter, :ref:`ECPlotter <ec-plotter>`, plots these variables.
The default exporter, ECExporter, exports these variables as well as an incrementer for
selecting data, ``cycle``.

The ``ec`` module
-----------------

.. automodule:: ixdat.techniques.ec
    :members:

The ``cv`` module
-----------------

.. automodule:: ixdat.techniques.cv
    :members: