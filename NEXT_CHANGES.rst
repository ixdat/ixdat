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

ixdat 0.2.12
============

techniques
^^^^^^^^^^

- Added ``faradaic_efficiency`` to correct the molecule flux as calculated from EC if not 100% FE,
	also works for isotopes

- Added ``integrate_flux()`` method to techniques\ms.py (analogous to existing integrate_signal(),
  allows for integration of the calibrated signal)
 

plotters
^^^^^^^^

- Added default colors for "H2@M2", "H2@M3", "H2@M4" to ms_plotter.py



For ixdat 0.3.0
===============
