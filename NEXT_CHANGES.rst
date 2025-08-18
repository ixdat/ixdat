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

For ixdat 0.3.1
===============

calculators
^^^^^^^^^^^

- The ``CalculatorPack`` has been added in `PR #196 <https://github.com/ixdat/ixdat/pull/196>`_, adressing `Issue #184 <https:github.com/ixdat/ixdat/issues/184>`_, to allow capturing multiple calculators from a measurement, saving them to JSON, and re-attaching them to other measurements for convenient, reproducible, and portable reuse of calibration workflows.
