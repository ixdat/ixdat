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

Debugging
---------

Other changes
-------------

- Tests for EC-lab parser using a dataset with multiple techniques and a dataset with looping techniques.

techniques
^^^^^^^^^^^
- ``_get_tspan_list`` in ``ECMSMeasurement`` now defaults ``t_steady_pulse`` to ``None`` instead of ``0``, which simplifies explanation in docstring and makes it more clear what it does (i.e. if now a pulse time of 0 is given it will actually use 0s instead of the entire pulse)

