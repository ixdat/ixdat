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

- The ``BrukerNMRReader`` (reader="bruker") has been added for reading Bruker
  TopSpin 1D NMR experiment folders. It uses the optional ``nmrglue`` package
  to parse the ``acqus`` parameter file and the processed real spectrum from
  ``pdata/<procno>/``, and returns an ``NMRSpectrum`` with a chemical-shift
  axis in ppm. Pass ``processed=False`` to get the raw FID as a ``FIDSpectrum``
  with a time axis in seconds reconstructed from ``SW_h``. Key acquisition
  parameters (``PULPROG``, ``SOLVENT``, ``BF1``, ``NS``, ``TE``, ``DATE`` ...)
  and processing parameters (``SI``, ``SF``, ``OFFSET`` ...) are lifted into
  the spectrum's ``metadata``, and the full ``acqus``/``procs`` dictionaries
  are preserved alongside. The constants ``ACQUS_KEYS`` and ``PROCS_KEYS`` are
  public and can be imported from ``ixdat.readers.bruker``.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_

techniques
^^^^^^^^^^

- New ``NMRSpectrum`` and ``NMRSpectrumSeries`` classes (technique ``"NMR"``
  / ``"NMR_spectra"``) and ``FIDSpectrum`` class (technique ``"FID"``) added
  in ``ixdat.techniques.nmr``, mirroring the ``FTIRSpectrum`` /
  ``FTIRSpectrumSeries`` pattern. ``NMRSpectrum`` uses ``NMRPlotter`` which
  inverts the x-axis to follow NMR convention (high ppm on the left).
  ``FIDSpectrum`` uses a plain plotter for the time-domain signal.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_

plotters
^^^^^^^^

- New ``NMRPlotter`` in ``ixdat.plotters.nmr_plotter``: subclasses
  ``SpectrumPlotter`` and inverts the x-axis so ``spec.plot()`` gives the
  standard NMR view out of the box.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_

tools
^^^^^

- New ``to_jsonable`` function in ``ixdat.tools``: recursively converts numpy
  arrays, numpy scalars, and byte strings into JSON-safe Python primitives.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_

