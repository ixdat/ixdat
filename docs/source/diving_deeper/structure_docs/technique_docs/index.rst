.. _techniques:

Techniques: ``ixdat``'s measurement subclasses
==============================================

TechniqueMeasurement classes (interchangeable with Techniques or Measurement subclasses)
inherit from the ``Measurement`` class (see :ref:`measurement`). Techniques allow for methods
and properties that are specific to the analysis of a certain type of data, providing maximum 
functionality for each of the supported techniques or combinations of techniques. 

Check out the following pages for documentation of the different technique modules published so far:

.. toctree::
    :maxdepth: 2

    electrochemistry
    mass_spec
    ec_ms
    sec
    spectra

An up-to-date, full list of the techniques and their names is given in the ``TECHNIQUE_CLASSES`` dictionary::

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