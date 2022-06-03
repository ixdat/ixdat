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

API Changes
-----------

techniques
^^^^^^^^^^

``ECMSMeasurement.ecms_calibration_curve`` now supports data specification using a
a selector. To do so, specify the section to use as numbers in the argument ``selector_list``,
the counter defining the sections (if different from the default selector) in ``selector_name``,
and the length of the steady-state period at the end of the pulse in ``t_steady_pulse``.
This can be much more convenient than manually specifying ``tspans``.
Implemented in #74.