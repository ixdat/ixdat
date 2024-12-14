# Calculators 

Calculations involving data from `Measurement`s in ixdat are handled by 
objects called `Calculator`s. Calculators can be attached to measurements
in order to make the data that they calculate seem seamlessly part of the
measurement.

Example workflow:

<img style="float: left;" src="../../../figures/calculator_workflow.svg" width="400px">

<pre>
from ixdat import Measurement
from ixdat.calculators import <b>MSCalculator</b>, MSInlet

# calibration experiment
<p style="color:blue;">
meas1 = Measurement.read(...)
</p>

# calibration

<b>calc = MSCalculator.gas_flux_calibration</b>(
<p style="color:blue;">
    measurement=meas1,
    mol=“H2”,
    mass=“M2”,
    tspan=[100, 200],
    inlet=MSInlet(),
</p>
)

# data to be quantified 
<p style="color:green;">meas2 = Measurement.read(...)</p>
<p style="color:green;">meas2.add_calculator</p>(<b>calc</b>)

# quantification
<p style="color:green;">
t, n_dot_H2 = meas2.grab(“n_dot_H2”, tspan=[400, 500])

meas2.plot(mol_list=[“H2”])
</p>
</pre>