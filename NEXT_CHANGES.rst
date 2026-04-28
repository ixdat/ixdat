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

- The ``EChemDBReader`` (reader="echemdb") has been added in `PR #194 <https://github.com/ixdat/ixdat/pull/194>`_ for reading CV reference data from echemdb.org, a curated open-access repository for digitized electrochemical datasets.

- The ``SpecDATReader`` (reader="spec") has been added for reading SPEC-format ``.dat``
  files produced by synchrotron beamlines running Certified Scientific Software's SPEC
  instrument-control program, such as ESRF BM31. This is distinct from
  ``QexafsDATReader`` (reader="qexafs"), which reads the Diamond B18-Core format: SPEC
  files use a structured multi-line header (``#F``, ``#E``, ``#S``, ``#L``, ``#O``,
  ``#P``, ``#T``, ``#UMI``, ...), space-separated columns, multiple scans per file
  delimited by ``#S``, and energy in keV. The reader supports transmission XAS
  (``technique="XAS"``), fluorescence XAS (``technique="XAS_fluorescence"``), and a
  full multi-column ``MultiSpectrum`` (no technique). Multiple scans within one file can
  be selected by ``scan_numbers`` and averaged with ``average_scans=True``. The x axis
  (``ZapEnergy`` by default) is converted from keV to eV. Because column names are
  beamline-specific, ``y_name`` (It) and ``ref_name`` (I0) must be supplied explicitly
  for XAS techniques. Metadata stored on the returned object includes the scan number,
  scan command, count time, motor positions (from ``#O``/``#P`` headers), UMI lines,
  original beamline file path (``#F``), and experiment comment (``#C``).
