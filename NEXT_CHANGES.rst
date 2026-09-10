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

- The ``AsimovReader`` (reader="asimov") has been added for retrieving measurements,
  spectra, and spectrum series from asimov.enci.dk via REST API. Authentication uses
  the OAuth 2.0 Device Flow via Keycloak. The reader auto-detects the object type
  from the payload.
  `PR #198 <https://github.com/ixdat/ixdat/pull/198>`_

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
  
- The ``XRDXYReader`` (reader="xrdxy") has been added for reading generic two- or
  three-column powder diffraction text files (.xy, .xye, or similar). It handles
  both whitespace- and comma-separated data, and scans comment and bare header lines to
  detect whether the x axis is 2-theta or Q-space (with units).
  `PR #203 <https://github.com/ixdat/ixdat/pull/203>`_

- The ``OceanViewTimeSeriesReader`` (reader="oceanview") gained an
  ``average_every`` option: rows are averaged in groups of that size as they
  are parsed, reducing the size of the resulting series. Defaults to 1
  (every row kept). A group smaller than ``average_every`` left over at the
  end of the file is still averaged and kept, not dropped.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

techniques
^^^^^^^^^^

- ``CyclicVoltammogram.redefine_cycle`` gained a ``turning_point`` option:
  with ``turning_point=True``, cycles are defined by reversals in the
  direction of the potential sweep (detected from the sign of dU/dt)
  instead of by crossing a fixed ``start_potential``. ``redox`` selects
  which reversals count (``True`` for negative-to-positive, ``False`` for
  positive-to-negative, ``None`` for either); ``N_points`` and ``N_sep``
  control how many consistent points must follow a candidate reversal and
  how far apart two reversals must be to both register.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``"CV-Optical"`` is now a registered technique combination
  (``TECHNIQUE_CLASSES``), so ``cv + optical_spectrum_series`` (where
  ``cv`` is a ``CyclicVoltammogram``) produces an ``ECOpticalMeasurement``
  directly, the same as ``ec + optical_spectrum_series`` already did for a
  plain ``ECMeasurement``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``ECOpticalMeasurement.get_dOD_cycle`` can now split a cycle into its
  anodic/cathodic halves (``direction=0``/``1``) using the same
  turning-point detection as ``redefine_cycle``, instead of only a fixed
  ``start_potential``-based split. ``get_dOD_difference_spectra``,
  ``get_dOD_cycle_noise``, and ``get_convergent_spectra`` accept the same
  ``direction``/``N_points`` and pass them through.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``ECOpticalMeasurement.denoise_spectra`` now takes an explicit
  ``spectra_field`` (defaulting to the measurement's own spectra) instead
  of always denoising ``self``, so it can be called on a measurement with
  any ``Field`` of spectral data, e.g. the output of
  ``get_dOD_difference_spectra``. It also gained a ``normalise`` option.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``get_convergent_spectra``'s ``smooth_distances`` is now a bool (use a
  Savitzky-Golay filter on the adjacent-distance array or not) rather than
  a filter size; the new ``window_length``/``polyorder`` control the
  filter. ``get_convergent_spectra`` also gained a ``normalise`` option.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- New ``NMRSpectrum`` and ``NMRSpectrumSeries`` classes (technique ``"NMR"``
  / ``"NMR_spectra"``) and ``FIDSpectrum`` class (technique ``"FID"``) added
  in ``ixdat.techniques.nmr``, mirroring the ``FTIRSpectrum`` /
  ``FTIRSpectrumSeries`` pattern. ``NMRSpectrum`` uses ``NMRPlotter`` which
  inverts the x-axis to follow NMR convention (high ppm on the left).
  ``FIDSpectrum`` uses a plain plotter for the time-domain signal.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_
  
- ``XRDSpectrum`` (a ``MultiSpectrum`` subclass) has been added in ``techniques/xrd.py``
  as the dedicated spectrum class for XRD data. The ``XRDXYReader`` returns it for all
  .xy and .xye files. For .xye files the per-point intensity error is stored as a second
  field alongside the intensity.
  `PR #203 <https://github.com/ixdat/ixdat/pull/203>`_

plotters
^^^^^^^^

- New ``NMRPlotter`` in ``ixdat.plotters.nmr_plotter``: subclasses
  ``SpectrumPlotter`` and inverts the x-axis so ``spec.plot()`` gives the
  standard NMR view out of the box.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_
  
- ``XRDSpectrumPlotter`` has been added in ``plotters/xrd_plotter.py`` as the default
  plotter for ``XRDSpectrum``. It plots intensity vs x and, when error data is present
  (i.e. for .xye files), overlays a shaded y+/-e band.
  `PR #203 <https://github.com/ixdat/ixdat/pull/203>`_

- New ``ECOpticalPlotter.plot_wavelengths_vs_cv``: plots tracked wavelengths'
  dO.D. against potential on one axis, with the CV current on a twinned
  right-hand axis, alongside the existing two-panel
  ``plot_wavelengths_vs_potential``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``ECOpticalPlotter.plot_dOD_difference_spectra`` gained ``denoise``,
  ``denoise_method``, ``PCA_explained_variance``, ``sg_window``, and
  ``sg_poly_order`` to optionally denoise the difference spectra (via
  ``denoise_spectra``) before plotting them.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``plot_waterfall_cycle``, ``plot_dOD_cycle_diff``, and
  ``plot_dOD_difference_spectra`` accept the new ``N_points`` (passed
  through to the underlying ``get_dOD_cycle``/``get_dOD_difference_spectra``
  calls); ``plot_convergent_spectra`` accepts the new
  ``window_length``/``polyorder`` (see the ``get_convergent_spectra``
  entry above), and its ``min_region_separation`` default changed from 50
  to 2 to match ``get_convergent_spectra``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

tools
^^^^^

- New ``to_jsonable`` function in ``ixdat.tools``: recursively converts numpy
  arrays, numpy scalars, and byte strings into JSON-safe Python primitives.
  `PR #200 <https://github.com/ixdat/ixdat/pull/200>`_
