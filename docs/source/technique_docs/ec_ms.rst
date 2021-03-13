.. _ec-ms:

Electrochemistry - Mass Spectrometry (EC-MS)
============================================

The main class for EC-MS data is the ECMSMeasurement.

It comes with the :ref:`EC-MS plotter <ec_ms-plotter>` which makes EC-MS plots like this one:

.. figure:: ../figures/ec_ms_annotated.svg
    :width: 600

    ``ECMSMeasurement.plot_measurement()``. Data from Trimarco, 2018.

Other than that it doesn't have much but inherits from both ``ECMeasurement`` and ``MSMeasurement``.
An ``ECMSMeasurement`` can be created either by adding an ``ECMeasurement`` and an ``MSMeasurement``
using the ``+`` operator, or by directly importing data using an EC-MS :ref:`reader <readers>`
such as "zilien".

Deconvolution, described in a publication under review, is implemented in the deconvolution module,
in a class inheriting from ``ECMSMeasurement``.

ixdat will soon have all the functionality and more for EC-MS data and analysis as the
legacy `EC_MS <https://github.com/ScottSoren/EC_MS>`_ package. This includes the tools
behind the EC-MS analysis and visualization in the puplications:

- Daniel B. Trimarco and Soren B. Scott, et al. **Enabling real-time detection of electrochemical desorption phenomena with sub-monolayer sensitivity**. `Electrochimica Acta, 2018 <https://doi.org/10.1016/j.electacta.2018.02.060>`_.

- Claudie Roy, Bela Sebok, Soren B. Scott, et al.  **Impact of nanoparticle size and lattice oxygen on water oxidation on NiFeOxHy**. `Nature Catalysis, 2018 <https://doi.org/10.1038/s41929-018-0162-x>`_.

- Anna Winiwarter and Luca Silvioli, et al. **Towards an Atomistic Understanding of Electrocatalytic Partial Hydrocarbon Oxidation: Propene on Palladium**. `Energy and Environmental Science, 2019 <https://doi.org/10.1039/C8EE03426E>`_.

- Soren B. Scott and Albert Engstfeld, et al.  **Anodic molecular hydrogen formation on Ru and Cu electrodes**. `Catalysis Science & Technology, 2020 <https://doi.org/10.1039/d0cy01213k>`_.

- Anna Winiwarter, et al.  **CO as a Probe Molecule to Study Surface Adsorbates during Electrochemical Oxidation of Propene**. `ChemElectroChem, 2021 <https://doi.org/10.1002/celc.202001162>`_.

- Soren B. Scott, et al.  **Tracking oxygen atoms in electrochemical CO oxidation –Part I: Oxygen exchange via CO2 hydration**. `Electrochimica Acta, 2021 <https://doi.org/10.1016/j.electacta.2021.137842>`_.

- Soren B. Scott, et al.  **Tracking oxygen atoms in electrochemical CO oxidation –Part II: Lattice oxygen reactivity in oxides of Pt and Ir**. `Electrochimica Acta, 2021 <https://doi.org/10.1016/j.electacta.2021.137844>`_.


The ``ec_ms`` module
--------------------
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques/ec_ms.py

.. automodule:: ixdat.techniques.ec_ms
    :members:

The ``deconvolution`` module
----------------------------
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques/deconvolution.py

.. automodule:: ixdat.techniques.deconvolution
    :members: