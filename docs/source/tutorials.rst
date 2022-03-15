.. _tutorials:

=========
Tutorials
=========

``ixdat`` has a growing number of tutorials and examples available.

Jupyter notebook tutorials
--------------------------
Jupyter notebooks are available in the ixdat Tutorials repository:
https://github.com/ixdat/tutorials/

This repository is a bit of a mess at the moment, apologies, but the tutorials themselves are
not bad, if we may say so ourselves. More are needed. Right now there are two,
both based on electrochemistry data:

Electrochemistry Tutorial 1: Reading and using data
***************************************************

Loading, selecting, calibrating, and exporting data, with a set of electrochemistry
files from a long complex measurement, taken with Biologic's EC-Lab, as an example.

Location (ixdat v0.2+): `electrochemistry/01_reading_and_using_data.ipynb <https://github.com/ixdat/tutorials/blob/ixdat_v0p2/electrochemistry/01_reading_and_using_data.ipynb>`_

Location (ixdat v0.1.x): `loading_appending_and_saving/export_demo_data_as_csv.ipynb <https://github.com/ixdat/tutorials/blob/main/loading_appending_and_saving/export_demo_data_as_csv.ipynb>`_

This tutorial shows with electrochemistry data how to load, append, and export data.
It shows, among other things, the **appending + operator** and how to use the **backend** (save() and get()).

It requires the data files `here <https://www.dropbox.com/sh/ag3pq7vqwuapd0o/AAB2Vqs6ZLZuFuMGp2ZeeWisa?dl=0>`_.


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

Article repositories
--------------------

Calibrating EC-MS data
**********************
See these two examples, respectively, for making and using an ixdat EC-MS calibration (here with isotope-labeled data):

- https://github.com/ScottSoren/pyCOox_public/blob/main/paper_I_fig_S1/paper_I_fig_S1.py

- https://github.com/ScottSoren/pyCOox_public/blob/main/paper_I_fig_2/paper_I_fig_2.py

EC-MS data analysis
*******************

This article has examples of analyzing and manually plotting data imported by ixdat

https://github.com/ScottSoren/Huang2021


Development scripts
-------------------
The basics of importing and plotting from each reader are demonstrated in
the **development_scripts/reader_testers** folder of the repository:
https://github.com/ixdat/ixdat/tree/user_ready/development_scripts/reader_testers