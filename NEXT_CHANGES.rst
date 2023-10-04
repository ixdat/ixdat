Instructions
============

Use this as a staging ground for CHANGES.rst. In other words, describe the
changes and additions to ixdat's API associated with your contribution. The idea is
that what you write here informs other developers what is on its way and then will be
copied to CHANGES.rst when the next version of ixdat is distributed. Please include
links to relevant Issues, Discussions, and PR's on github with the following format
(replace XX):

`Issue #XX <https://github.com/ixdat/ixdat/issues/XX>`_

`PR #XX <https://github.com/ixdat/ixdat/pulls/XX>`_

ixdat 0.2.6
===========

API changes
-----------

Readers
^^^^^^^
- A function ``fix_WE_potential`` has been added to the ``biologic`` module. This
  corrects the mistake in some .mpt files that are exported with zeros for "<Ewe>/V".

Techniques
^^^^^^^^^^
- MS measurements now have access to the powerful ``spectro_inlets_quantification`` (siq)
  package as a plugin. See `siq's documentation <https://spectroinlets-spectro-inlets-quantification.readthedocs-hosted.com/en/v1.1/theory/quant_ms.html>`_
  for details.

  To activate the plugin, use::

    import ixdat
    ixdat.plugins.activate_si_quant()

  after activation, the following calibration and quantification methods will use siq:
  - ``MSMeasurement.gas_flux_calibration`` now returns a siq ``CalPoint``. Note that this
    method is only available with siq activated - native ixdat gas flux calibration
    is only available as a method of ``MSInlet``. Otherwise they work basically the same.
  - ``MSMeasurement.multicomp_gas_flux_calibration`` returns a ``Calibration``. Note that it
    solves a matrix equation to deconvolute overlapping peaks in a multi-analyte
    calibration gas.
  - ``ECMSMeasurement.ecms_calibration`` and ``ECMSMeasuremnt.ecms_calibration_curve`` each
    return an object of ``CalPoint``.
  - ``MSMeasurement.set_quantifier`` is used to attache a siq quantifier to the measurement.
    Usage: ``ms.set_quantifier(calibration=my_siq_calibration)``
  - ``MSMeasurement.grab_fluxes`` is a new method which uses the attached quantifier to
    quantify the fluxes of all molecules covered therein. In general it solves a matrix
    equation to deconvolute overlapping signals.
  - ``MSMeasurement.grab_flux`` will, when siq is activated, run ``MSMeasurement.grab_fluxes``
    and return the requested flux vector.

  SIQ comes with data for a small batch of molecules. To supplement this collection of
  yaml-formatted molecule data, place new files in the QUANT_DIRECTORY. This can be set
  as follows (change to the folder where you have your data)::

    ixdat.plugins.si_quant.QUANT_DIRECTORY = "~/projects/batteries/quantification_data"

=======
readers
^^^^^^^

- ``cinfdata_db`` is a new direct db reader for cinfdatabase at DTU SurfCat

plotters
^^^^^^^^

- ``MSPlotter.plot_measurment`` now create a new axis if no initial axis is given
  before initiating right axes in the case of dual plotting on left and right axes.
  Fixes `Issue #97 <https://github.com/ixdat/ixdat/issues/97>`_


- ``SpectrumSeriesPlotter.heat_plot`` now accept max_threshold and min_threshold and 
   scanning_mask to include or exclude specific values from scanning variable
  
- ``SpectroMSPlotter`` new plotter for ``SpectroMSMeasurment`` now create a new axis if no initial axis is given
  before initiating right axes in the case of dual plotting on left and right axes.

techniques
^^^^^^^^^^

- ``ReactorMeasurement`` class, technique="reactor", with a ``TPMSPlotter``. This
  technique is analogous to EC-MS with temperature replacing potential and
  pressure replacing current.

- ``SpectroMSMeasurement`` class set ``SpectroMSPlotter`` as default plotter

dev
^^^
- Renamed development scripts which are not software tests "demo" instead of "test".

- Skip py36 because github is having a problem building it. See, for example, here:
https://github.com/ixdat/ixdat/actions/runs/3876991446/jobs/6611480640#step:3:7

- Do black test before software tests in github CI to save time
