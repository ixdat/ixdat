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

readers
^^^^^^^

- Normalized EC-Lab column names like `<Ewe/V>` to `<Ewe>/V>` in `BiologicReader` to ensure consistent parsing and aliasing.
  This prevents failures in tools like `fix_WE_potential` when non-standard column names are encountered.
  `Issue #193 <https://github.com/ixdat/ixdat/issues/193>`_
  `PR #195 <https://github.com/ixdat/ixdat/pull/195>`_
