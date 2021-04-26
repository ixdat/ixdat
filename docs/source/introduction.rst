============
Introduction
============

.. figure:: figures/ixdat_profile_pic.svg
    :width: 200

    The power of combining techniques (fig made with ``EC_Xray``, an ``ixdat`` precursor)


``ixdat`` will provide a powerful **object-oriented** interface to experimental data, especially in-situ experimental data for which it is of interest to combine data obtained simultaneously from multiple techniques.

``ixdat`` will replace the existing electrochemistry - mass spectrometry data tool, `EC_MS <https://github.com/ScottSoren/EC_MS>`_, and will thus become a powerful stand-alone tool for analysis and visualization of data acquired by the equipment of `Spectro Inlets <https://spectroinlets.com>`_ and other EC-MS solutions.
It will also replace the existing electrochemistry - synchrotron GIXRD data tool, `EC_Xray <https://github.com/ScottSoren/EC_Xray>`_ when needed.
Over time, it will acquire functionality for more and more techniques.

In addition to a **pluggable** parser interface for importing your data format, it will include pluggable exporters and plotters, as well as a database interface.

We will update this documentation as features are added.

``ixdat`` is free and open source software and we welcome input and new collaborators.
The source is here: https://github.com/ixdat/ixdat

For a long motivation, see :ref:`concept`.


Installation
------------

To use ``ixdat``, you need to have python installed. We recommend
`Anaconda python <https://www.anaconda.com/products/individual>`_.

To install ``ixdat``, just type in your terminal or Anaconda prompt::

    $ pip install ixdat

And hit enter.

``ixdat`` is under development, and to make use of the newest features,
you may need to upgrade to the latest version. This is also easy. Just type::

    $ pip install --upgrade ixdat

Tutorials
----------
``ixdat`` has a growing number of tutorials available as jupyter notebooks. The tutorials
are here: https://github.com/ixdat/tutorials

.. toctree::
    :maxdepth: 2

    extended-concept