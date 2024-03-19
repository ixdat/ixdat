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
Added FTIR and EC-FTIR. The latter for now just inherits from SpectroECMeasurement

readers
^^^^^^^
Added "opus_ftir" reader for text-exported files from Opus FTIR spectrumeter.

plotters
^^^^^^^^
Added a plot_stacked_spectra method to SpectrumSeriesPlotter, SpectroMeasurementPlotter,
and SECPlotter. This is the default plotting method for FTIR.