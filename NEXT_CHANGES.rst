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

ixdat 0.2.7
===========

API Changes
-----------



Debugging
---------

general
^^^^^^^
- ``EC_MS`` is no longer a dependency
  Resolves `Issue #128 <https://github.com/ixdat/ixdat/issues/124>`_

exporters
^^^^^^^^^
- Fixed exporting and re-importing of ``ECOpticalMeasurment``s for new pandas version
  Resolves `Issue #124 <https://github.com/ixdat/ixdat/issues/124>`_

constants
^^^^^^^^^
- ``BOLTZMAN_CONSTANT`` renamed ``BOLTZMANN_CONSTANT``
  Resolves `Issue #125 <https://github.com/ixdat/ixdat/issues/125>`_