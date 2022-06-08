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


For ixdat 0.2.3
===============

API Changes
-----------

measurement
^^^^^^^^^^^
- Added example usage to docstring of ``Measurement.select_values``

readers
^^^^^^^
- Added an ``XRDMLReader`` (reader="xrdml") for xml-formatted XRD spectrum files from,
  for example, Empyrion's software. Usage::

    from ixdat import Spectrum
    my_xrd = Spectrum.read("my_file.xrdml", reader="xrdml")
    my_xrd.plot()

- Added an ``AvantageAVGReader`` (reader="avantage") for avantage-exported spectrum files from,
  for example, Thermo Fisher's K-alpha equipment. Usage::

    from ixdat import Spectrum
    my_xps = Spectrum.read("my_file.avg", reader="avantage")
    my_xps.plot()

- Added a ``QEXAFSDATReader`` (reader="qexafs") for .dat files exported by Diamond
  Synchrotron's B18-Core beamline. These data files include X-Ray Absorption across
  multiple detector elements, X-ray fluorescence, and diagnostic information such as primary
  beam intensity. All the data is read together as a new object called a ``MultiSpectrum``.
  If the ``read`` function is called with ``technique="XAS"``, a single spectrum with the
  XAS data (the "FFIO" column) is returned::

    from ixdat import Spectrum
    multi_spectrum = Spectrum.read("my_file.dat", reader="qexafs")
    my_xas = multi_spectrum["QexafsFFI0"]   # index retrieves a Spectrum from a MultiSpectrum
    # OR
    my_xas = Spectrum.read("my_file.dat", reader="qexafs", technique="XAS")
    my_xas.plot()

- A new reading method ``Spectrum.read_set()`` imports and appends multiple spectrum
  files. It is thus similar to ``Measurement.read_set()``. The appended spectra are
  returned as a ``SpectrumSeries``, which has a heat plot as its default ``plot()``.

- Both ``Measurement.read_set()`` and the new ``Spectrum.read_set()`` allow the user
  to specify a string that appears in the middle of the names of the files to be read
  and appended, rather than just the beginning::

    from ixdat import Spectrum
    spectrum_series = Spectrum.read_set(
       part="data/XAS_spectrum", suffix=".dat", technique="XAS"
    )  # reads .dat files in the directory "data" with "XAS_spectrum" in their name
    spectrum_series.plot()  # heat plot with time on x-axis and energy on y-axis


techniques
^^^^^^^^^^

- ``CyclicVoltammogram`` now has a ``plot_cycles`` function. This plots all
  the cycles in the cv color coded by cycle number, with a scale bar.

- Refactoring of ``SpectroECMeasurement`` in order to generalize aspects of combining
  time-resolved measurements with spectra taken simultaneously. This results in the
  following classes:

  - ``SpectroMeasurement(Measurement)`` for any time-resolved measurement data with spectra.
  - ``SpectroECMeasurement(SpectroMeasurement, ECMeasurement)``, ``technique="SEC"``
    for an EC measurement with spectra. By default plots the EC data on a lower panel and
    the spectral data as a heat plot on the upper panel.
  - ``ECXASMeasurement(SpectroECMeasurement)``, ``technique="EC-XAS"`` for an EC measurement
    with XAS spectra.
  - ``ECOpticalMeasurement(SpectroECMeasurement)``, ``technique=EC-Optical`` for an EC
    measurement with optical spectroscopy. A reference spectrum is required, and can also
    be set using a reference time or reference potential after the data is loaded. By
    default, optical density = **-log(y/y0)**, with ``y=spectrum_series.y`` and
    ``y0=reference_spectrum.y``, is plotted instead of raw data.

  Before the refactor, ``ECOpticalMeasurement`` had been called ``SpectroECMeasurement``.

- Addition of a ``Measurement`` and a ``SpectrumSeries`` gives a ``SpectroMeasurement``
  or a subclass thereof determined by hyphenating the technique. For example::

    from ixdat import Measurement, Spectrum

    my_ec = Measurement.read("my_ec_data.mpt", reader="biologic")
    print(my_ec, ",", my_ec.technique)
    # ECMeasurement(...) , EC

    my_xas = Spectrum.read_set("data/my_xas", suffix=".dat", technique="XAS")
    print(my_xas, ",", my_ec.technique)
    # SpectrumSeries(...) , XAS

    my_ec_xas = my_ec + my_xas
    print(my_ec_xas, ",", my_ec_xas.technique)
    # ECXASMeasurement(...) , EC-XAS

    my_ec_xas.plot()  # XAS data heat plot in top panel and EC data in bottom panel.

plotters
^^^^^^^^

- Plotters for spectroelectrochemistry were refactored in a similar way to the technique
  classes (see above, under "techniques"). All ``SpectroECMeasurements`` come with plotters
  with ``plot_measurement()`` (heat plot and EC data vs time) and ``plot_vs_potential()``
  (heat plot and EC data vs potential). For ``ECOpticalMeasurements``, these plot
  optical density rather than raw data on the heat plots.

exporters
^^^^^^^^^
- The main .csv file exported from an ``ECOpticalMeasurement`` refers to its auxilliary
  files as::

    'spectrum_series' in file: 'exported_sec_spectra.csv'
    'reference_spectrum' in file: 'exported_sec_reference.csv'

  Files exported by ``SpectroECMeasurement`` in ixdat<0.2.3 will need these two lines modified, i.e.
  "spectra"->"spectrum_series" and "reference"->"reference_spectrum", before they can
  be imported. After this modification they can be read as before by
  ``Measurement.read(path_to_main_csv_file, reader="ixdat")``.

Debugging
---------

- ``PVMassSpecReader`` (reader="pfeiffer") now identifies columns as masses when they have
  names assigned during the data acquisition (e.g. "32_amu_(Oxygen)" as a column name),
  whereas before it would only have identified as masses columns which ended with "_amu".

- ``ECPlotter.plot_measurement`` and ``ECMSPlotter.plot_measuremnt`` raise warnings instead
  of ``SeriesNotFoundError``\ s if they can't find the requested
  potential (``U_name``) or current (``J_name``).

- ``SpectrumSeriesPlotter.heat_plot()`` now correctly orients its data. Previously it
  had been flipped in the vertical direction.