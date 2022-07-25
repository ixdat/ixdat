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

API changes
-----------
techniques
^^^^^^^^^^^
- Improved docstring for ``ECMSMeasurement.ecms_calibration_curve()`` to include the new additions from previous release.

- Added MSInlet``gas_flux_calibration_curve`` to enable multiple point calibration using calculated gas flux
either with different concentrations in carrier gas or at different inlet pressures. Note, concentration needs to be given in ppm, as the flux calculation uses various constants from the carrier gas molecule instead of a mixture, which will lead to significant inaccuracy for high concentrations.

- ``Measurement.select`` is now even more versatile. A user can specify a ``selector_name``
  for args to work on. This enables selection based on columns with a space in them, like
  "cycle number".
  Resolves `Issue #77 <https://github.com/ixdat/ixdat/issues/77>`_

readers
^^^^^^^
``Mesurement.read_set`` can now take a folder as the first argument, in which case  all
files in that folder (with the specified suffix) are appended.
Resolves `Issue #88 <https://github.com/ixdat/ixdat/issues/88>`_
``Measurement.read_set`` now also raises a ``ReadError` rather than returning ``None`` in
the case of no matching files.

``Measurement.read`` (and by extension ``Measurement.read_set``) can now be called
without a specified reader for certain known file types. To see which file types, use::

  from ixdat.readers.reading_tools import DEFAULT_READER_NAMES
  print(DEFAULT_READER_NAMES)

plotters
^^^^^^^^

``ECMSPlotter.plot_measurment`` and ``plot_potential`` can take string options "all",
"ms", and "ec" to specify the time span of the data to plot.
Resolves `Issue #91 <https://github.com/ixdat/ixdat/issues/91>`_

Debugging
---------

readers
^^^^^^^
The biologic reader now checks for "Ns" and "cycle number" rather than assuming it
knows which EC techniques have which of these selector-defining series.
Resolves `Issue #87 <https://github.com/ixdat/ixdat/issues/87>`_

techniques
^^^^^^^^^^^
- ``Measurement.cut`` now skips empty component measurements rather than crashing.
  Resolves `Issue #93 <https://github.com/ixdat/ixdat/issues/93>`_

- ``MSMeasurement.reset_bg`` works again! It now adds a new calibration with bg=0 for
  masses that had previously had a bg set.
  Resolves `Issue #82 <https://github.com/ixdat/ixdat/issues/82>`_

- ``_get_tspan_list`` in ``ECMSMeasurement`` now defaults ``t_steady_pulse`` to ``None``
  instead of ``0``, which simplifies explanation in docstring and makes it more clear what
  it does (i.e. if now a pulse time of 0 is given it will actually use 0s instead of the
  entire pulse)


Other changes
-------------

- Tests for EC-lab parser using a dataset with multiple techniques and a dataset with looping techniques.

