.. _cheat_sheet:

===========
Cheat sheet
===========

Importing a measurement
-----------------------

.. code:: ipython3

	my_measurement = ixdat.Measurement.read(path_to_file, reader, **kwargs)

 	"""	
   	Return a Measurement object from parsing a file with the specified reader
    
    	Args:
       		path_to_file (Path or str): The path to the file to read
        	reader (str or Reader class): The (name of the) reader to read the file with.
        	kwargs: key-word arguments are passed on to the reader's read() method.
	"""


Check what's in a ``Measurement``
---------------------------------
If you want to know which data your ``Measurement`` object contains, you can ask for:

``my_measurement.series_names`` - returns the names of the data and time series that can be used in ``grab()``

``my_measurement.series_list`` - returns a list of all data and time series in the object


Plotting
--------

To make the standard plot associated with a ``Measurement`` object, simply write:

``my_measurement.plot()`` or ``my_measurement.plot_measurement()``

use ``help(my_measurement.plot)`` to find out which arguments are accepted.

Calling the above functions to make an ixdat plot will return a list of axes. If you want to be able to modify the standard plot, save the axes in a variable.

``axes = my_measurement.plot()``

 This you can then use get a handle on the figure:

``fig_1 = axes[0].get_figure()``

You can use any standard matplotlib commands to make modifications to these axes and figure.


Calibration
-----------
Calibrations are stored in the ``Measurement`` object as a list. Depending on the type of data (and therefore ``Measurement`` object, the options for calibration will vary.

If you have your calibration stored in a ``Calibration`` object, you apply it like this:

``my_measurement.add_calibration(my_calibration_object)``

Other calibrations can be added using ``calibrate``, for example:

``my_ec_measurement.calibrate(RE_vs_RHE=0.67,  A_el=0.197)``
``my_ms_measurement.calibrate(ms_cal_results = [cal_h2, cal_o2, cal_c2h4])``, where cal_h2, cal_o2 and cal_c2h4 are the MSCalResults objects containing the calibration factors for H2, O2 and C2H4, respectively

Built-in constants
------------------
ixdat has a range of built-in constants saved in the sub-module ``constants``, accesible e.g. via:

``from ixdat.constants import FARADAY_CONSTANT``



Use ``help()`` to get the documentation
---------------------------------------

For any object or function, you can use ``help(my_random_object)`` to return the docstring of my_random_object - this can be very useful e.g. for finding methods associated with that object.