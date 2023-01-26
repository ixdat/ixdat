.. _sec-tutorial:

================================
Spectroelectrochemistry tutorial
================================

This tutorial demonstrates importing, plotting, and exporting in-operando UV-Vis data
as an example of spectroelectrochemistry (S-EC).
It shows delta optical density calculation and both calculation and plotting of the full 2-D data field and
cross sections (i.e. spectra and wavelength-vs-time).

The example data is not yet publically available, thus there is no code execution shown here. 

Spectroelectrochemistry with ixdat
==================================

The data
--------

.. code:: ipython3

    from pathlib import Path
    
    data_dir = Path.home() / "Dropbox/ixdat_resources/test_data/sec"
    
    print("data directory containes files and folders:")
    for path in data_dir.iterdir():
        print(path)

SEC while scanning potential
----------------------------

.. code:: ipython3

    from ixdat import Measurement
    
    sec = Measurement.read(
        data_dir / "test-7SEC.csv",
        path_to_ref_spec_file = data_dir / "WL.csv",
        path_to_U_J_file = data_dir / "test-7_JV.csv",
        scan_rate=1,
        tstamp=1,
        reader="msrh_sec"
    )

.. code:: ipython3

    sec.plot_measurement()

The raw data as numpy arrays is accessible as such:

.. code:: ipython3

    wavelength = sec.wavelength.data
    absorption_matrix = sec.spectra.data
    potential = sec.potential.data
    
    print("wavelength:")
    print(wavelength)
    
    print("potential:")
    print(potential)
    
    print("absorption:")
    print(absorption_matrix)
    

.. code:: ipython3

    sec.calibrate_RE(RE_vs_RHE=0.26)
    sec.normalize_current(A_el=1)
    sec.plot_measurement(make_colorbar=True)

.. code:: ipython3

    sec.plot_waterfall()

.. code:: ipython3

    sec.reference_spectrum.plot()

.. code:: ipython3

    sec.set_reference_spectrum(V_ref=0.66)
    sec.reference_spectrum.plot()

.. code:: ipython3

    help(sec)
    print(sec.calc_dOD().data)
    sec.plot_waterfall()
    

.. code:: ipython3

    axes = sec.plot_measurement(cmap_name="jet")
    axes[0].set_xlabel("time / [s]")

.. code:: ipython3

    sec.plot_vs_potential()

.. code:: ipython3

    sec.plot_vs_potential(cmap_name="jet")
    
    axes = sec.plot_vs_potential(cmap_name="jet", vspan=[1.0, 1.45], make_colorbar=True, color="g")
    axes[1].set_ylim([0, 0.1])

.. code:: ipython3

    sec.plot_wavelengths(wavelengths=["w800", "w580"])

.. code:: ipython3

    sec.plot_wavelengths_vs_potential(wavelengths=["w800", "w580", "w480"])

.. code:: ipython3

    spec_1 = sec.get_spectrum(V=1.0)
    
    ax = spec_1.plot(color="r", label="1.0 V vs RHE")
    sec.reference_spectrum.plot(ax=ax, label="reference")
    ax.legend()

.. code:: ipython3

    spec_1 = sec.get_dOD_spectrum(V=1.0, V_ref=0.66)
    spec_1.plot()

.. code:: ipython3

    ax = sec.get_dOD_spectrum(V=1.0, V_ref=0.66).plot(color="b", label="species 1")
    sec.get_dOD_spectrum(V=1.4, V_ref=1.0).plot(color="g", label="species 2", ax=ax)
    sec.get_dOD_spectrum(V=1.7, V_ref=1.4).plot(color="r", label="species 3", ax=ax)
    ax.legend()

.. code:: ipython3

    sec.export(data_dir / "sec_export.csv")
    
    print(sec.reference_spectrum.name)

.. code:: ipython3

    
    
    sec_reloaded = Measurement.read(data_dir / "sec_export.csv", reader="ixdat")
    
    sec_reloaded.set_reference_spectrum(V_ref=0.66)
    
    sec_reloaded.plot_vs_potential(cmap_name="jet")

.. code:: ipython3

    [(s.name, s.shape) for s in sec_reloaded.series_list]
    
    

Open-circuit potential decay
----------------------------

.. code:: ipython3

    sec_decay = Measurement.read(
        data_dir / "decay/PDtest-1.35-1OSP-SP.csv",
        path_to_ref_spec_file=data_dir / "WL.csv",
        path_to_t_U_file=data_dir / "decay/PDtest-1.35-1OSP-E-t.csv",
        path_to_t_J_file=data_dir / "decay/PDtest-1.35-1OSP-J-t.csv",
        tstamp=1,
        reader="msrh_sec_decay",
    )

.. code:: ipython3

    sec_decay.plot_measurement()

.. code:: ipython3

    sec_decay.calibrate_RE(RE_vs_RHE=0.26)
    sec_decay.set_reference_spectrum(t_ref=5)
    sec_decay.plot_measurement(cmap_name="jet")

.. code:: ipython3

    sec_decay.plot_wavelengths(wavelengths=["w480", "w600", "w850"])

.. code:: ipython3

    sec_decay.export(data_dir / "sec_decay_export.csv")

The End
-------



