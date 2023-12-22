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


ixdat 0.2.9
===========


techniques
^^^^^^^^^^

``MSMeasurement.remove_matrix_interference`` associates a background with the measurement
  which subtracts a ratio (measured when there's only solvent) times the signal of the primary solvent
  mass from the signal of specified masses. To use the background, put ``remove_background=True``
  as argument to ``plot()`` or ``grab``.
  Background classes are implemented in a way parallel to Calibration classes.

ECOpticalMeasurement.get_spectrum() now has an option not to interpolate.
  Passing the argument `interpolate=False` gets it to choose nearest spectrum instead.