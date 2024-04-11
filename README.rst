.. figure:: docs/source/figures/logo.svg
    :width: 200

=============================================
``ixdat``: The In-situ Experimental Data Tool
=============================================

With ``ixdat``, you can import, combine, and export complex experimental datasets
as simply as::

    ec = Measurement.read_set("awesome_EC_data", reader="biologic")
    ec.plot_measurement()

    ms = Measurement.read("2021-03-30 16_59_35 MS data.tsv", reader="zilien")
    ms.plot_measurement()

    ecms = ec + ms
    ecms.plot_measurement()

    ecms.export("my_combined_data.csv")

Output:

.. figure:: docs/source/figures/ixdat_example_figures.png
    :width: 700

    In-situ experimental data made easy

Or rather than exporting, you can take advantage of ``ixdat``'s powerful analysis
tools and database backends to be a one-stop tool from messy raw data to public
repository accompanying your breakthrough publication and advancing our field.

Version
-------
This is the latest version.

For changes up to this version, see `CHANGES.rst <https://github.com/ixdat/ixdat/blob/main/CHANGES.rst>`_

For ixdat 0.1.9 see the `v0.1.x branch <https://github.com/ixdat/ixdat/tree/v0.1.x>`_.

About
-----

``ixdat`` provides a powerful **object-oriented** interface to experimental data, especially in-situ experimental data for which it is of interest to combine data obtained simultaneously from multiple techniques.

Documentation is at https://ixdat.readthedocs.io

In addition to a **pluggable** parser interface for importing your data format, ``ixdat`` also includes
pluggable exporters and plotters, as well as a database interface. A relational model of experimental data is
designed into every level.

.. list-table:: Techniques and Readers
   :widths: 20 15 50
   :header-rows: 1


   * - Measurement technique
     - Status
     - Readers
   * - Electrochemistry
     - Released
     - - biologic: .mpt files from Biologic's EC-Lab software
       - autolab: ascii files from AutoLab's NOVA software
       - ivium: .txt files from Ivium's IviumSoft software
   * - Mass Spectrometry
     - Released
     - - pfeiffer: .dat files from Pfeiffer Vacuum's PVMassSpec software
       - cinfdata: text export from DTU Physics' cinfdata system
       - zilien: .tsv files from Spectro Inlets' Zilien software
   * - Electrochemistry - Mass Spectrometry
     - Released
     - - zilien: .tsv files from Spectro Inlets' Zilien software
       - EC_MS: .pkl files from the legacy EC_MS python package
   * - Spectroelectrochemistry
     - Released
     - - msrh_sec: .csv file sets from Imperial College London's SEC system
   * - X-ray photoelectron spectroscopy (XPS)
     - Development
     - - avantage: .avg files from Thermo Scientific's Avantage software
   * - X-ray diffraction (XRD)
     - Development
     - - xrdml: .xrdml files from e.g. PanAnalytical's Empyereon
   * - In-situ Electrochemistry - X-ray adsorption spectroscopy
     - Development
     - - qexafs: .dat files from Diamond's B18 beamline
   * - Low-Energy Ion Scattering (LEIS)
     - Future
     -

Tutorials are provided at https://ixdat.readthedocs.io/en/latest/tutorials/index.html

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


Article repositories
--------------------

``ixdat`` is shown in practice in a growing number of open repositories of data and analysis
for academic publications:

- Soren B. Scott, et al.  **Tracking oxygen atoms in electrochemical CO oxidation –Part I: Oxygen exchange via CO2 hydration**. `Electrochimica Acta, 374, 137842 <https://doi.org/10.1016/j.electacta.2021.137842>`_, **2021**.

  Repository: https://github.com/ScottSoren/pyCOox_public

- Soren B. Scott, et al.  **Tracking oxygen atoms in electrochemical CO oxidation –Part II: Lattice oxygen reactivity in oxides of Pt and Ir**. `Electrochimica Acta, 374, 137844 <https://doi.org/10.1016/j.electacta.2021.137844>`_, **2021**.

  Repository: https://github.com/ScottSoren/pyCOox_public

- Kevin Krempl, et al. **Dynamic Interfacial Reaction Rates from Electrochemistry - Mass Spectrometry**. `Journal of Analytical Chemistry. 93, 7022-7028 <https://doi.org/10.1021/acs.analchem.1c00110>`_, **2021**

  Repository: https://github.com/kkrempl/Dynamic-Interfacial-Reaction-Rates

- Junheng Huang, et al. **Online Electrochemistry−Mass Spectrometry Evaluation of the Acidic Oxygen Evolution Reaction at Supported Catalysts**. `ACS Catal. 11, 12745-12753 <https://doi.org/10.1021/acscatal.1c03430>`_, **2021**

  Repository: https://github.com/ScottSoren/Huang2021


Join us
-------

``ixdat`` is free and open source software and we welcome input and new collaborators. Please help us improve ``ixdat``!

Contact us (https://github.com/ixdat/ixdat/discussions or sbscott@ic.ac.uk) or just
`get started developing <https://ixdat.readthedocs.io/en/latest/developing/index.html>`_.
