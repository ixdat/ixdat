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

For ixdat 0.4.1
===============

development_scripts
^^^^^^^^^^^^^^^^^^^

- The demos have moved to the new ``ixdat.demos`` package, grouped like ixdat into
  ``readers``, ``calculators`` and ``techniques``. Each demo has loader functions
  which return its ixdat objects, an ``EXAMPLES`` dictionary of these loaders, and a
  ``main(show=True)`` which plots the results and returns them as a dictionary. Other
  scripts import the loaders, e.g.
  ``from ixdat.demos.readers.demo_autolab_reader import load_single_file``.
  The demo data is an archive on ERDA, given in ``ixdat/demos/demo_data.ini``.
  ``python -m ixdat.demos.download`` downloads it with ``pooch``, checks its sha256
  checksum and unzips it in ixdat's cache folder. ``get_demo_data_dir()`` returns the
  demo data folder: the environment variable ``IXDAT_DEMO_DATA_DIR``, else
  ``demo_data/`` in the repository root, else the downloaded folder.
  ``development_scripts/run_all_demos.py`` runs every demo and lists the ones that
  fail. ``development_scripts/save_all_demo_data_in_sqlite.py`` saves every demo
  example in an SQLite file and shows its tables in a browser.
  `PR #212 <https://github.com/ixdat/ixdat/pull/212>`_

plotters
^^^^^^^^

- The plotters get colormaps with ``plt.get_cmap``. Matplotlib 3.9 removed
  ``mpl.cm.get_cmap``, which broke waterfall plots, colorbars, ``plot_cycles`` and the
  SEC wavelength plots.
  `PR #212 <https://github.com/ixdat/ixdat/pull/212>`_

techniques
^^^^^^^^^^

- ``as_cv()`` of an ``ECMSSpectroMeasurement``, e.g. a zilien file with mass scans,
  returns an ``ECMSCyclicVoltammogram`` with the EC and MS data and leaves out the mass
  spectra. It raised a ``TechniqueError`` before.
  `PR #212 <https://github.com/ixdat/ixdat/pull/212>`_

spectra
^^^^^^^

- Indexing a ``SpectrumSeries`` whose technique is registered to a ``SpectrumSeries``
  class, e.g. ``meas[1]`` on a zilien MS measurement with mass scans, returns a
  ``Spectrum`` again. Since the ``Spectrum.from_dict`` added in 0.4.0 it returned a
  broken ``SpectrumSeries`` which could not be plotted.
  `PR #212 <https://github.com/ixdat/ixdat/pull/212>`_
