0.2.0 (2022-3-20)
=================

API Changes
-----------

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

  The old arguments names are deprecated.

- Renamed key-word argument in MS and EC-MS plotting functions:

  - ``remove_backround`` (replaces ``removebackground``

  The old argument name is deprecated

Bug Fixes
---------

