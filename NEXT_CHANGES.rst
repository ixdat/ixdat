Instructions
============

Use this as a staging ground for CHANGES.rst. In other words, describe the
changes and additions to ixdat's API associated with your contribution. The idea is
that what you write here informs other developers what is on its way and then will be
copied to CHANGES.rst when the next version of ixdat is distributed. Please include
links to relevant Issues, Discussions, and PR's on github with the following format
(replace XX):

`Issue #XX <https://github.com/ixdat/ixdat/issues/XX>`_
`PR #XX <https://github.com/ixdat/ixdat/pull/XX>`_

For ixdat 0.3.1
===============

readers
^^^^^^^

- The ``NordicTDMSReader`` (reader="nordic") has been added for reading the .tdms files
  produced by Nordic potentiostat. It reads current, potential, and impedance, with units
  taken directly from the TDMS channel properties (current is converted from A to mA to
  follow ixdat convention). An ISO 8601 datetime is included in the measurement metadata.
  If a ``*.EC_Macro`` file is present alongside the data file, the experiment sequence
  (hardware control, gas, CV, and hold steps) is also parsed into the metadata.
  `PR #167 <https://github.com/ixdat/ixdat/pull/167>`_

- The ``XRDXYReader`` (reader="xrdxy") has been added for reading generic two- or
  three-column powder diffraction text files (.xy, .xye, or similar). It handles
  both whitespace- and comma-separated data, and scans comment and bare header lines to
  detect whether the x axis is 2-theta or Q-space (with units).
  `PR #203 <https://github.com/ixdat/ixdat/pull/203>`_

- The ``EChemDBReader`` (reader="echemdb") has been added in `PR #194 <https://github.com/ixdat/ixdat/pull/194>`_ for reading CV reference data from echemdb.org, a curated open-access repository for digitized electrochemical datasets.

techniques
^^^^^^^^^^

- ``XRDSpectrum`` (a ``MultiSpectrum`` subclass) has been added in ``techniques/xrd.py``
  as the dedicated spectrum class for XRD data. The ``XRDXYReader`` returns it for all
  .xy and .xye files. For .xye files the per-point intensity error is stored as a second
  field alongside the intensity.
  `PR #203 <https://github.com/ixdat/ixdat/pull/203>`_

plotters
^^^^^^^^

- ``XRDSpectrumPlotter`` has been added in ``plotters/xrd_plotter.py`` as the default
  plotter for ``XRDSpectrum``. It plots intensity vs x and, when error data is present
  (i.e. for .xye files), overlays a shaded y+/-e band.
  `PR #203 <https://github.com/ixdat/ixdat/pull/203>`_
