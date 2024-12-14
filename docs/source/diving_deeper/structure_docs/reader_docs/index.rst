.. _readers:

Readers: getting data into ``ixdat``
====================================
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/readers

On this page, you can find the documentation for the different readers for importing data obtained using different techniques. For easier navigation between the different sections, use the menu on the left. 


Initiating a measurement
------------------------

A typical workflow is to start by reading a file. For convenience, most readers are
accessible directly from ``Measurement``. So, for example, to read a .mpt file exported
by Biologic's EC-Lab, one can type:

>>> from ixdat import Measurement
>>> ec_meas = Measurement.read("my_file.mpt", reader="biologic")

See :ref:`readers <readers>` for a description of the available readers.

The biologic reader (``ixdat.readers.biologic.BiologicMPTReader``) ensures that the
object returned, ``ec_meas``, is of type ``ECMeasurement``.

Another workflow starts with loading a measurement from the active ``ixdat`` backend.
This can also be done straight from ``Measurement``, as follows:

>>> from ixdat import Measurement
>>> ec_meas = Measurement.get(3)

Where the row with id=3 of the measurements table represents an electrochemistry
measurement. Here the column "technique" in the measurements table specifies which
TechniqueMeasurement class is returned. For row three of the measurements
table, the entry "technique" is "EC", ensuring ``ec_meas`` is an object of type
``ECMeasurement``.


A full list of the readers thus accessible and their names can be viewed by typing:

>>> from ixdat.readers import READER_CLASSES
>>> READER_CLASSES

.. _ixdatcsvreader:

Reading .csv files exported by ixdat: The ``IxdatCSVReader``
------------------------------------------------------------

``ixdat`` can export measurement data in a .csv format with necessary information in the
header. See :ref:`exporters`. It can naturally read the data that it exports itself. Exporting and reading,
however, may result in loss of raw data (unlike ``save()``).

The ``ixdat_csv`` module
........................

.. automodule:: ixdat.readers.ixdat_csv
    :members:

Importing from other experimental data platforms
------------------------------------------------

**cinfdata** is a web-based database system for experimental data, developed and used at DTU SurfCat
(formerly CINF) in concert with The ``PyExpLabSys`` suite of experimental data acquisition tools.
Both are available at https://github.com/CINF.

As of yet, ``ixdat`` only has a text-file reader for data exported from **cinfdata**, but
in the future it will also have a reader which downloads from the website given e.g. a
setup and date.

The ``cinfdata`` module
.......................

.. automodule:: ixdat.readers.cinfdata
    :members:

Electrochemistry and sub-techniques
------------------------------------
These are readers which by default return an ``ECMeasurement``.
(See :ref:`electrochemistry`)

The ``biologic`` module
.......................

.. automodule:: ixdat.readers.biologic
    :members:

The ``autolab`` module
......................

.. automodule:: ixdat.readers.autolab
    :members:

The ``ivium`` module
....................

.. automodule:: ixdat.readers.ivium
    :members:

The ``chi`` module
..................

.. automodule:: ixdat.readers.chi
    :members:

Mass Spectrometry and sub-techniques
------------------------------------
These are readers which by default return an ``MSMeasurement``.
(See :ref:`mass-spec`)

The ``pfeiffer`` module
.......................

.. automodule:: ixdat.readers.pfeiffer
    :members:

The ``rgasoft`` module
......................

.. automodule:: ixdat.readers.rgasoft
    :members:

EC-MS and sub-techniques
------------------------
These are readers which by default return an ``ECMSMeasurement``.
(See :ref:`ec-ms`)

The ``zilien`` module
.....................

.. automodule:: ixdat.readers.zilien
    :members:


The ``ec_ms_pkl`` module
........................

.. automodule:: ixdat.readers.ec_ms_pkl
    :members:


SEC and sub-techniques
----------------------
These are readers which by default return a ``SpectroECMeasurement``.
(See :ref:`sec`)

The ``msrh_sec`` module
.......................

.. automodule:: ixdat.readers.msrh_sec
    :members:


Other techniques
----------------

The ``avantage`` module
.......................

.. automodule:: ixdat.readers.avantage
    :members:


The ``xrdml`` module
....................

.. automodule:: ixdat.readers.xrdml
    :members:


The ``qexafs`` module
.....................

.. automodule:: ixdat.readers.qexafs
    :members:

