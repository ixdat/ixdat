.. _tutorials:

=========
Tutorials
=========

``ixdat`` has a growing number of tutorials and examples available, addressing different steps in your ``ixdat`` workflow.

Getting data in and out
-----------------------

In ``ixdat`` your data is imported with a ``reader`` (specific to your data format) into a ``Measurement`` object and can then be plotted and/or 
exported into a format of your choice. You can figure out how in this tutorial: 

.. toctree::
    :maxdepth: 1

    data_in_and_out

In the future there will also be a tutorial on data manipulation and treatment including how to cut and select different sections of your data, as well as 
an introduction to calibration.

Techniques
----------

Each technique comes with its own methods developed specifically for the types of data handling that are of interest. For example, the electrochemistry technique subclass comes with a subclass ``CyclicVoltammogram`` which allows you to select your data according to cycles, the mass spectrometry technique subclass comes with a method to calibrate the MS signals, and the spectroelectrochemistry technique subclass with a method to set a reference spectrum. 
The tutorials for the technique-specific methods are still work in progress - please, bear with us - you can find all the ones already available below.

.. list-table:: 
   :widths: 30 30
   :header-rows: 1

   * - Technique
     - Tutorials
   * - EC
     - :ref:`ec_tutorials`
   * - MS
     - 
   * - EC-MS
     - :ref:`ecms_tutorial`
   * - SEC
     - :ref:`sec-tutorial`  
   * - XRD
     -


Advanced data handling
----------------------

While the existing technique subclasses allow for a wide range of general and technique specific data treatment, sometimes this is not enough for your 
individual needs. Luckily, ``ixdat`` also allows for direct data handling in array-type form. How to access data in this way is demonstrated in the tutorials on advanced data handling.

As an open source project, ``ixdat`` is always happy to get contributions from the users, e.g. in form of new techniques.


.. toctree::
    :maxdepth: 1

    advanced_tutorials


Download Jupyter notebook tutorials
-----------------------------------
All tutorials linked to above are also available as fully interactive Jupyter notebooks in the ixdat Tutorials repository:
https://github.com/ixdat/tutorials/

This repository is a bit of a mess at the moment, apologies, but the tutorials themselves are
not bad, if we may say so ourselves. More are in progress.


Development scripts
-------------------
For developers:

The basics of importing and plotting from each reader are demonstrated in
the **development_scripts/reader_testers** folder of the repository:
hhttps://github.com/ixdat/ixdat/tree/main/development_scripts/reader_demonstrators