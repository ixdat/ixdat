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

Quant integration
=================

general
-------
An object called ``plugins`` can be imported from the ``config`` module. This gives access
to variables controlling plugin behaviour.``plugins.USE_QUANT`` determines whether an
external quantification package is used for quantification by ``MSMeasurement``.
``plugins.QUANT_DIRECTORY`` determines where that external quantification package looks
for reference data.

A module ``quant_patch`` is included with functions that act on objects of the external
quantification package, with the hope that the functionality migrates to that package in
the near future.

techniques
----------
``MSMeasurement`` has new methods which can be used if ``plugins.USE_QUANT`` is True:

- ``MSMeasurement.gas_flux_calibration`` for sensitivity factor determination by flowing
  a pure gas through an inlet chip
- ``MSMeasurement.multicomp_gas_flux_calibration`` for sensitivity matrix determination
  by flowing a gas with multiple analyte components. This uses reference spectra and
  results in sensitivity factors that are as inaccurate as the reference spectra. Its use
  is therefore discouraged (but sometimes necessary).
- ``MSMeasurement.set_quantifier`` sets the *quantifier*, which then determines how
  ``MSMeasurement.grab_flux`` and ``MSMeasurement.grab_fluxes`` calculate fluxes of
  molecules to the vacuum chamber.
- ``MSMeasurement.grab_fluxes`` uses the measurement's quantifier to calculate the fluxes
  of all the molecules in the quantifier's ``mol_list`` with the signals at all the
  masses in the quantifier's ``mass_list`` as inputs. It takes the tspan and background
  arguments familiar in ixdat from other ``grab`` methods.

The workings of ``MSMeasurement.grab_flux`` are changed if ``plugins.USE_QUANT`` is True.
In that case, it invokes the *quantifier* via ``MSMeasurement.grab_fluxes`` and returns
just the flux of the requested molecule.

New Zilien reader
=================

A new Zilien reader for the SpectroInlets' new Zilien dataset file format
version. The new dataset version is able to integrate the Biologic EC-lab
dataset during a measurement. Such dataset contains a new series with a header
name ``EC-lab`` and two meta columns ``experiment_number`` and
``technique_number`` in the series. The new reader is using the columns during
the process of creating Ixdat series objects. The objects, created from the
Zilien dataset, match exactly the objects created from the Biologic MPT files,
like when read one by one, **except for** the timestamps from the Biologic
series.  The Biologic series timestamps are incremented by a time offset when
the Biologic EC-lab measurement was triggered.

E.g. when you start a Zilien measurement and then trigger an EC-lab measurement
after five seconds, the timestamps in the series of the Biologic dataset part
will be higher by five, comparing to the timestamps in the MPT files.


Debugging
=========

- xrdml reader can now import files where the data is labeled "counts" rather than
  "intensities", as the text exports from the Royce Institute XRD