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


For ixdat 0.2.3
===============

API Changes
-----------

readers
^^^^^^^
- Added an ``XRDMLReader`` (reader="xrdml") for xml-formatted XRD spectrum files from,
  for example, Empyrion's software. Usage::

    from ixdat import Spectrum
    my_xrd = Spectrum.read("my_file.xrdml", reader="xrdml")
    my_xrd.plot()

- Added an ``AvantageAVGReader`` (reader="avantage") for avantage-exported spectrum files from,
  for example, Thermo Fisher's K-alpha equipment. Usage::

    from ixdat import Spectrum
    my_xps = Spectrum.read("my_file.avg", reader="avantage")
    my_xps.plot()


techniuqes
^^^^^^^^^^

- CyclicVoltammogram now has a ``plot_cycles`` function. This plots all
  the cycles in the cv color coded by cycle number, with a scale bar.


Debugging
---------

- PVMassSpecReader (reader="pfeiffer") now identifies columns as masses when they have
  names assigned during the data acquisition (e.g. "32_amu_(Oxygen)" as a column name),
  whereas before it would only have identified as masses columns which ended with "_amu".
