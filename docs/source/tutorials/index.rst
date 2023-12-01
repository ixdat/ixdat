.. _tutorials:

Tutorials
=========

``ixdat`` has a growing number of tutorials and examples available, addressing different steps in your ``ixdat`` workflow.

Browse the tutorials in the pages of the documentation here, or download the jupyter notebook files to try yourself.


Overview
--------

Each technique comes with its own methods developed specifically for the types of data handling that are of interest. For example, the electrochemistry technique subclass comes with a subclass ``CyclicVoltammogram`` which allows you to select your data according to cycles, the mass spectrometry technique subclass comes with a method to calibrate the MS signals, and the spectroelectrochemistry technique subclass with a method to set a reference spectrum.
The tutorials for the technique-specific methods are still work in progress - please, bear with us - you can find all the ones already available below.

As an open source project, ``ixdat`` is always happy to get contributions from the users, e.g. in form of new techniques.

.. list-table::
   :widths: 30 30
   :header-rows: 1

   * - Subject
     - Tutorials
   * - **EC**
     - :ref:`ec_tutorials`
   * - **EC-MS**
     - :ref:`ecms_tutorial`
   * - **SEC**
     - :ref:`sec_tutorial`
   * - **Data organization**
     - :ref:`data_tutorial`

..
  * - Biologic
    - :ref:`biologic_tutorials`
  * - Zilien
    - :ref:`zilien_tutorials`
  * - Database
  * - :ref:`database_tutorials`

  This commented-out list can be used for tutorials we want.

.. _readers_tutorial:

Readers
-------


In ``ixdat`` your data is imported with a ``reader`` (specific to your data format) into a ``Measurement`` object and can then be plotted and/or
exported into a format of your choice. You can figure out how in this tutorial:


Jupyter notebook file: https://github.com/ixdat/tutorials/blob/af1b38f2096555904786bb9fbbaf79bc209d4ef9/L1_basic_concepts/ixdat_intro_readers.ipynb

.. toctree::
  :maxdepth: 0

  tutorials_repo/L1_basic_concepts/ixdat_intro_readers


Technqiues
----------

.. _ec_tutorials:

Electrochemistry
^^^^^^^^^^^^^^^^
Jupyter notebook file: https://github.com/ixdat/tutorials/blob/af1b38f2096555904786bb9fbbaf79bc209d4ef9/L2_techniques/electrochemistry/ixdat_demo_cyclic_voltammetry.ipynb

.. toctree::
  :maxdepth: 0

  tutorials_repo/L2_techniques/electrochemistry/ixdat_demo_cyclic_voltammetry


.. _ecms_tutorial:

EC-MS
^^^^^
Jupyter notebook file: https://github.com/ixdat/tutorials/blob/af1b38f2096555904786bb9fbbaf79bc209d4ef9/L2_techniques/ec_ms_quantification/ixdat_demo_EC-MS.ipynb

.. toctree::
  :maxdepth: 0

  tutorials_repo/L2_techniques/ec_ms_quantification/ixdat_demo_EC-MS

.. _sec_tutorial:

spectroelectrochemistry
^^^^^^^^^^^^^^^^^^^^^^^
Jupyter notebook file: https://github.com/ixdat/tutorials/blob/af1b38f2096555904786bb9fbbaf79bc209d4ef9/L2_techniques/spectroelectrochemistry/ixdat_demo_spectroelectrochemistry.ipynb

.. toctree::
  :maxdepth: 0

  tutorials_repo/L2_techniques/spectroelectrochemistry/ixdat_demo_spectroelectrochemistry


.. _data_tutorial:


Data concepts
-------------

While the existing technique subclasses allow for a wide range of general and technique specific data treatment, sometimes this is not enough for your
individual needs. Luckily, ``ixdat`` also allows for direct data handling in array-type form. How to access data in this way is demonstrated in the tutorials on advanced data handling.

Jupyter notebook file: https://github.com/ixdat/tutorials/blob/af1b38f2096555904786bb9fbbaf79bc209d4ef9/L3_data_structure/ixdat_demo_manipulating_electrochemistry_data.ipynb

.. toctree::
  :maxdepth: 0

  tutorials_repo/L3_data_structure/ixdat_demo_manipulating_electrochemistry_data


..
  Advanced data handling
  ^^^^^^^^^^^^^^^^^^^^^^

  In the future there will also be a tutorial on data manipulation and treatment including how to cut and select different sections of your data, as well as
  an introduction to calibration.



The jupyter notebooks
---------------------

All tutorials linked to above are also hosted and version controlled in the ixdat Tutorials repository:
https://github.com/ixdat/tutorials/

This repository is a bit of a mess at the moment, apologies, but the tutorials themselves are
not bad, if we may say so ourselves. More are in progress.


Write your own tutorial!
^^^^^^^^^^^^^^^^^^^^^^^^
Let the world know how to analyze data of your technique. We are very happy
for pull requests to the tutorials repository! Please clear all output of
your jupyter notebook when you make your pull request to that repository,
as jupiter notebook output does not verson control well. Please also make
it clear in first markdown cell where one can download the data needed to
run it (link to another repository, dropbox folder, or zenodo DOI, for
example). More details to come on the tutorials repo README.


Development scripts
-------------------
For developers:

The basics of importing and plotting from each reader are demonstrated in
the **development_scripts/reader_demonstrators** folder of the repository:
https://github.com/ixdat/ixdat/tree/main/development_scripts/reader_demonstrators/