tutorials
#########
jupyter notebook tutorials for ``ixdat``
========================================
This is a collection of tutorials and worked examples for use of ``ixdat``,
the `In-situ Experimental Data Tool <https://ixdat.readthedocs.io>`_.

The tutorials here are intended to demonstrate general features of ``ixdat`` and
give a sense of the design of the package, as well as key features.

Links to specific examples of ixdat usage showing how to make published figures 
and derive published results can be found in article repositories. For a complete
list of tutorials and article repositories, see the ixdat documentation:

https://ixdat.readthedocs.io/en/latest/tutorials.html

Setup
=====

- Install python with, at a minimum: numpy, scipy, matplotlib, and Jupyter notebook.
  Anaconda python is recommended

- Install the latest version of ixdat with ``pip install --upgrade ixdat``

  This main version of the tutorials repository works with ixdat v0.2
  For a version compatible with ixdat v0.1.x, see the `ixdat_v0p1 branch <https://github.com/ixdat/tutorials/tree/ixdat_v0p1>`_

- Download or clone this repository or the tutorial you are interested in

- If necessary, download the example data it uses (it will be described below).

- Open the tutorial in Jupyter notebook and it should run!

- Help us with your feedback! If something was unclear in the tutorial, it's probably
  because we need to improve it.


The Tutorials
=============

Electrochemistry
----------------


Tutorial 1: Reading and Using Data
..................................

Location: `electrochemistry/01_reading_and_using_data.ipynb <https://github.com/ixdat/tutorials/blob/main//electrochemistry/01_reading_and_using_data.ipynb>`_

This tutorial shows with electrochemistry data how to load, append, and export data.
It shows, among other things, the **appending + operator** and how to use the **backend** (save() and get()).

It requires the data files `here <https://www.dropbox.com/sh/ag3pq7vqwuapd0o/AAB2Vqs6ZLZuFuMGp2ZeeWisa?dl=0>`_.


Tutorial 2: Comparing Cycles
............................

Location: `electrochemistry/02_comparing_cycles.ipynb <https://github.com/ixdat/tutorials/blob/main//electrochemistry/02_comparing_cycles.ipynb>`_

This tutorial, together with the previous one, shows the ``ixdat``'s API for electrochemistry data.
It demonstrates, with CO stripping as an example, the following features:

- Selecting cyclic voltammatry cycles

- Integrating current to get charge passed

- Lining seperate cycles up with respect to potential

It reads ixdat-exported data directly from github.
A worked example based on the methods in this tutorial


Spectroelectrochemistry
-----------------------

Location: `spectroelectrochemistry/spectroelectrochemistry_demo.ipynb <https://github.com/ixdat/tutorials/blob/main/spectroelectrochemistry/spectroelectrochemistry_demo.ipynb>`_

The sample data is not yet publically available.

EC-MS
-----

Location: `ec_ms_quantification/EC-MS_ixdat_tutorial.ipynb <https://github.com/ixdat/tutorials/blob/main//ec_ms_quantification/EC-MS_ixdat_tutorial.ipynb>`_

This tutorial shows how to use ixdat to analyse electrochemical and gas calibration data and how to apply a calibration to EC-MS data, as well as some ways of how to use standard matplotlib functions to modify EC-MS plots generated using ixdat.

In order to run this tutorial, download the data from https://zenodo.org/record/8400063 (DOI:10.5281/zenodo.8400063) and place it in the same folder as the ipython notebook.


If interested, check out these two examples, respectively, for making and using an ixdat EC-MS calibration:

- https://github.com/ScottSoren/pyCOox_public/blob/main/paper_I_fig_S1/paper_I_fig_S1.py

- https://github.com/ScottSoren/pyCOox_public/blob/main/paper_I_fig_2/paper_I_fig_2.py


Demos
=====

22E06_ICL
---------

Four short scripts use to demonstrate ixdat on May 6, 2022 at Imperial College London.

The data is available here:
https://www.dropbox.com/sh/0xy80ytu9oyykkr/AACcXTRozwDbESFUzeV7bOx5a?dl=0

The contents should be downloaded into a folder **tutorials/data** for the
relative paths used in the scripts to be accurate.

Demos 02 and 03 work with ixdat version 0.2.2. Demos 01 and 04 require aspects of version 0.2.3
(``CyclicVoltammogram.plot_cycles()`` and ``XRDMLReader``, respectively) which are available
pre-release as ixdat 0.2.3dev0: https://pypi.org/project/ixdat/0.2.3.dev0/

For developers
==============
Clear output of jupyter notebooks before committing! Add a pre-commit hook that clears all ipython notebook output, as described here:
https://medium.com/somosfit/version-control-on-jupyter-notebooks-6b67a0cf12a3

A pre-commit hook is prepared for you. Just copy it to the git folder::

  cp pre_commit_hook .git/hooks/pre-commit
