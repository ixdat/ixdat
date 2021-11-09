.. _electrochemistry:

Electrochemistry
================

The main TechniqueMeasurement class for electrochemistry is the ``ECMeasurement``.
Sublcasses of ``ECMeasurement`` include ``CyclicVoltammagram`` and ``CyclicVoltammagramDiff``.

Direct-current electrochemsitry measurements (``ixdat`` does not yet offer specific
functionality for impedance data) are characterized by the essential quantities being
working-electrode current (in loop with the counter electrode) and potential (vs the
reference electrode) as a function of time. Either current or potential can be controlled
as the input variable, so the other acts at the response, and it is common to plot
current vs potential, but in all cases both are tracked or controlled as a function of
time. This results in the essential variables ``t`` (time), ``v`` (potential), and ``j``
(current). The main job of ``ECMeasurement`` and subclasses is to give standardized,
convenient, and powerful access to these three variables for data selection, analysis,
and visualization, regardless of which hardware the data was acquired with.

The default plotter, :ref:`ECPlotter <ec-plotter>`, plots these variables.
The default exporter, ECExporter, exports these variables as well as an incrementer for
selecting data, ``cycle``.

Electrochemistry is the most thoroughly developed technique in ``ixdat``. For in-depth
examples of the functionality in the ``ECMeasurement`` class and its subclasses, see
the following Tutorials:

- `Loading appending and saving <https://github.com/ixdat/tutorials/blob/main/loading_appending_and_saving/export_demo_data_as_csv.ipynb>`_

- `Analyzing cyclic voltammagrams <https://github.com/ixdat/tutorials/blob/main/simple_ec_analysis/difference_between_two_cvs.ipynb>`_

The ``ec`` module
-----------------
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques/ec.py

.. figure:: ../figures/ec_subplots.svg
  :width: 600
  :alt: Example plots. left: ``ECMeasurement.plot_vs_potential()`` right: ``ECMeasurement.plot_measurement()``

  left: ``ECMeasurement.plot_vs_potential()`` right: ``ECMeasurement.plot_measurement()``. `See tutorial <https://github.com/ixdat/tutorials/blob/main/simple_ec_analysis/difference_between_two_cvs.ipynb>`_

.. automodule:: ixdat.techniques.ec
    :members:

.. _`cyclic_voltammetry`:

The ``cv`` module
-----------------
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques/cv.py

.. figure:: ../figures/cv_diff.svg
  :width: 300
  :alt: Example ``CyclicVoltammagramDiff`` plot

  output of ``CyclicVoltammagramDiff.plot()``.  `Tutorial <https://github.com/ixdat/tutorials/blob/main/loading_appending_and_saving/export_demo_data_as_csv.ipynb>`_.

.. automodule:: ixdat.techniques.cv
    :members: