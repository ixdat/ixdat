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
  results in sensitivity factors that are as innacurrate as the reference spectra. Its use
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

