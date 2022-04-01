Instructions
============

Use this as a staging ground for CHANGES.rst. In other words, describe the
changes and additions to ixdat's API associated with your contribution. The idea is
that what you write here informs other developers what is on its way and then will be
copied to CHANGES.rst when the next version of ixdat is distributed. Please include
links to relevant Issues, Discussions, and PR's on github with the following format
(replace XX):

`Issue #XX <https://github.com/ixdat/ixdat/issues/XX>`_

`PR #XX <https://github.com/ixdat/ixdat/pulls/XX>`_


Debugging for 0.2.1
-------------------

Will add the following when we finish debugging (will be 0.2.1 then). But some will go
in the 0.2.0 section of CHANGES.rst as indicated:

API changes
...........

techniques
^^^^^^^^^^
- Lookups instead of properties (as of 0.2.0)

  - ``my_ec_meas["raw_potential"]`` replaces ``ECMeasurement.raw_potential``
  - ``my_ec_meas["raw_current"]`` replaces ``ECMeasurement.raw_current``
  - ``my_cv["scan_rate"]`` replaces ``CyclicVoltammogram.scan_rate``

  The old properties are deprecated.

- ``ECMSCyclicVoltammogram.diff_with()`` raises a ``NotImplementedError`` rather than
  an obscure error message. (as of 0.2.1)

- ``ECMSPlotter.plot_vs_potential`` can accept ``color`` as a keyword argument. This
  colors both the U-J curve in the lower panel and all the signals in the top panel,
  so best to use together with a ``mass_list`` or ``mol_list`` of length 1. (as of 0.2.1)

exporters
^^^^^^^^^

- ``export()`` functions now take a ``time_step`` argument. It must be used together with
  a ``tspan``. It is effective in reducing the file size of exported .csv's. `PR #60 <https://github.com/ixdat/ixdat/pull/60>`_

- Renamed keyword argument: ``columns`` replaces ``v_list`` in exporters (``ECExporter`` in v0.2.0, ``ECMSExporter`` in 0.2.1)

  The old name is deprecated.

plotters
^^^^^^^^

- Added interactive range selection functionality to all plotters,
  `PR #61 <https://github.com/ixdat/ixdat/pull/61>`_. Left mouse click will set the
  left marker, right mouse click the right marker, and double clicking with any of the
  buttons will remove that marker. When both left and right markers are in place, the
  selected span will be printed out.

constants
^^^^^^^^^
- Dynamic viscosities are temperature-dependent, `PR #55 <https://github.com/ixdat/ixdat/pull/55>`_ 
  This enables accurate gas flux MS calibration (for the gases included) accross a range of temperatures.

Debugging
.........

- ``ECMSExporter`` works as of 0.2.1 (it had been broken in 0.2.0).

  This solves `Issue #56 <https://github.com/ixdat/ixdat/issues/56>`_

- Other exporters were also sometimes prone to having extra ``aliases`` leading to
  circular lookups in the measurement when they were loaded. This has been fixed (as of 0.2.1)

- ``RE_vs_RHE=0`` counts as calibrated potential!