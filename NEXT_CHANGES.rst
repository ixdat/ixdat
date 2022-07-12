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

``Measurement.read`` (and by extension ``Measurement.read_set``) can now be called
without a specified reader for certain known file types. To see which file types, use::

  from ixdat.readers.reading_tools import DEFAULT_READER_NAMES
  print(DEFAULT_READER_NAMES)

``Measurement.select`` is now even more versatile. A user can specify a ``selector_name``
for args to work on. This enables selection based on columns with a space in them, like
"cycle number".
Resolves `Issue #77 <https://github.com/ixdat/ixdat/issues/77>`_

Debugging
---------

readers
^^^^^^^
The biologic reader now checks for "Ns" and "cycle number" rather than assuming it
knows which EC techniques have which of these selector-defining series.
Resolves `Issue #87 <https://github.com/ixdat/ixdat/issues/87>`_

techniques
^^^^^^^^^^
``MSMeasurement.reset_bg`` works again! It now adds a new calibration with bg=0 for
masses that had previously had a bg set.
Resolves `Issue #82 <https://github.com/ixdat/ixdat/issues/82>`_
