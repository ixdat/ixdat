.. _techniques:

Techniques: ``ixdat``'s measurement subclasses
==============================================
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques

TechniqueMeasurement classes (interchangeable with Techniques or Measurement subclasses)
inherit from the ``Measurement`` class (see :ref:`measurement`)

A full list of the techniques and there names is in the ``TECHNIQUE_CLASSES`` dictionary::

    >>> from ixdat.techniques import TECHNIQUE_CLASSES
    >>> TECHNIQUE_CLASSES  # note, more techniques may have been added since!
    {
	'simple': <ixdat.measurements.Measurement>,
	'EC': <ixdat.techniques.ec.ECMeasurement>,
	'CV': <ixdat.techniques.cv.CyclicVoltammogram>,
	'MS': <ixdat.techniques.ms.MSMeasurement>,
	'EC-MS': <ixdat.techniques.ec_ms.ECMSMeasurement>,
	'XRD': <ixdat.spectra.Spectrum>,
	'XPS': <ixdat.spectra.Spectrum>,
	'XAS': <ixdat.spectra.Spectrum>,
	'EC-Optical': <ixdat.techniques.spectroelectrochemistry.ECOpticalMeasurement>,
	'SEC': <ixdat.techniques.spectroelectrochemistry.SpectroECMeasurement>,
	'EC-XAS': <ixdat.techniques.spectroelectrochemistry.ECXASMeasurement>
    }

.. toctree::
    :maxdepth: 2

    electrochemistry
    mass_spec
    ec_ms
    sec
    spectra
