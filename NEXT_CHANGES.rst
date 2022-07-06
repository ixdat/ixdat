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

readers
^^^^^^^
``Mesurement.read_set`` can now take a folder as the first argument, in which case  all
files in that folder (with the specified suffix) are appended.
Resolves `Issue #88 <https://github.com/ixdat/ixdat/issues/88>`_


Debugging
---------

readers
^^^^^^^
The biologic reader now checks for "Ns" and "cycle number" rather than assuming it
knows which EC techniques have which of these selector-defining series.
Resolves `Issue #87 <https://github.com/ixdat/ixdat/issues/87>`_
