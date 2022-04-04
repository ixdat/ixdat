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


reader
^^^^^^^^

- ``cinfdata_db`` is a new direct db reader for cinfdatabase at DTU SurfCat

plotters
^^^^^^^^

- ``MSPlotter.plot_measurment`` now create a new axis if no initial axis is given
  before initiating right axes in the case of dual plotting on left and right axes.
  Fixes `Issue #97 <https://github.com/ixdat/ixdat/issues/97>`_

techniques
^^^^^^^^

- ``ReactorMeasurement`` class, technique="reactor", with a ``TPMSPlotter``. This
  technique is analogous to EC-MS with temperature replacing potential and
  pressure replacing current.
