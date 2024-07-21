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


For ixdat 0.2.13
===============

Debugging
---------
- Fixed timestamp form in ``QexafsDATReader`` to correctly parse timezone all year.


API changes
-----------

- Time-resolved x-ray flouresence (``technique = "TRXRF"``) implemented in `PR #168 <https://github.com/ixdat/ixdat/pull/168>`_:

  - ``B18TRXRFReader`` (reader="b18_trxrf") implemented for reading TRXRF data from the Diamond lightsource beamline B18TRXRFReader

  - ``TRXRFMeasurement`` with a series constructor method for the value series of interest, "FF_over_I0", and ``TRXRFPlotter`` for plotting the TRXRF data.

  - Hyphenation of TRXRF with EC (``technique = "EC-TRXRF"``) implemented (syntax: ``ec_txrf = ec + trxrf``) in ``ECTRXRFMeasurement`` and ``ECTRXRFPlotter``
