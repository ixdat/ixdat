Instructions
============

Use this as a staging ground for CHANGES.rst. In other words, describe the
changes and additions to ixdat's API associated with your contribution. The idea is
that what you write here informs other developers what is on its way and then will be
copied to CHANGES.rst when the next version of ixdat is distributed. Please include
links to relevant Issues, Discussions, and PR's on github.

Coming changes for ixdat 0.2.0
==============================

API Changes
-----------

ixdat.measurement
^^^^^^^^^^^^^^^^^

- Generalization of multiple calibrations:

  ``Measurement.calibration`` is deprecated in favor of ``Measurement.calibration_list``.
  Use ``add_calibration()``. Setting ``calibration`` is deprecated.


ixdat.techniques
^^^^^^^^^^^^^^^^

- Renamed measurement technique class: ``CyclicVoltammogram`` (replaces ``CyclicVoltammagram``)
  The old name is deprecated.

- Renamed properties in ``ECMeasurement`` and inheriting classes:

  - ``v_name`` (replaces ``V_str``)
  - ``j_name`` (replaces ``J_str``)
  - ``E_name`` (replaces ``E_str``)
  - ``I_name`` (replaces ``I_str``)

  The old property names are deprecated.

- Renamed keyward in ``MSMeasurement.plot_flux()``:

  - ``remove_backround`` (replaces ``removebackground``

  The old argument name is deprecated

ixdat.plotters
^^^^^^^^^^^^^^

- Axes order: Two-panel figures with shared x-axis always return a list of axes in the order
  ``axes=[top_left, bottom_left]`` in the case of one y-axis of each panel, or
  ``[top_left, bottom_left, top_right, bottom_right]`` in the case of an additional y-axis
  in either panel. Note that ``axes[2]`` or ``axes[3]`` might be ``None``. For example,
  ``axes = ECMSPlotter.plot_measurement()`` by default results in ``axes[0]`` being the
  MS data, ``axes[1]`` electrode potential, ``axes[2]=None`` and ``axes[3]`` being the
  electrode current.

- Renamed key-word arguments in EC, EC-MS, and SEC plotting functions:

  - ``v_name`` (replaces ``V_str``)
  - ``j_name`` (replaces ``J_str``)
  - ``v_color`` (replaces ``V_color``)
  - ``j_color`` (replaces ``J_color``)

  The old key-word argument names are deprecated.

- Renamed key-word argument in MS and EC-MS plotting functions:

  - ``remove_backround`` (replaces ``removebackground``

  The old argument name is deprecated

Bug Fixes
---------
