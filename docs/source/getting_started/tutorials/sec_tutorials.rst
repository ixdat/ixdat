.. _tutorials:

=========================
Electrochemistry tutorial
=========================

``ixdat`` has a growing number of tutorials and examples available.

Jupyter notebook tutorials
--------------------------
Jupyter notebooks are available in the ixdat Tutorials repository:
https://github.com/ixdat/tutorials/

This repository is a bit of a mess at the moment, apologies, but the tutorials themselves are
not bad, if we may say so ourselves. More are needed. Right now there are two,
both based on electrochemistry data:

Electrochemistry Tutorial 2: Comparing cycles of a CV
*****************************************************

Location (ixdat v0.2+): `electrochemistry/02_comparing_cycles.ipynb <https://github.com/ixdat/tutorials/blob/ixdat_v0p2/electrochemistry/02_comparing_cycles.ipynb>`_

Location (ixdat v0.1.x): `simple_ec_analysis/difference_between_two_cvs.ipynb <https://github.com/ixdat/tutorials/blob/main/simple_ec_analysis/difference_between_two_cvs.ipynb>`_

This tutorial, together with the previous one, shows ``ixdat``'s API for electrochemistry data.
It demonstrates, with CO stripping as an example, the following features:

- Selecting cyclic voltammatry cycles

- Integrating current to get charge passed

- Lining seperate cycles up with respect to potential

It reads ixdat-exported data directly from github.


Spectroelectrochemistry
***********************

.. _sec-tutorial:

Location:`spectroelectrochemistry/ <https://github.com/ixdat/tutorials/blob/main/spectroelectrochemistry/>`_

This tutorial demonstrates importing, plotting, and exporting in-operando UV-Vis data
as an example of spectroelectrochemistry (S-EC).
It shows delta optical density calculation and both calculation and plotting of the full 2-D data field and
cross sections (i.e. spectra and wavelength-vs-time).

The example data is not yet publically available.
