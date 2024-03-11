"""Module for deconvolution of mass transport effects."""

from .ec_ms import ECMSMeasurement
from scipy.optimize import curve_fit  # noqa
from scipy.interpolate import interp1d  # noqa
from scipy import signal  # noqa
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth  # noqa
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft, ifftshift, fftfreq  # noqa
import numpy as np

# FIXME: too much abbreviation in this module.
# TODO: Implement the PR review here: https://github.com/ixdat/ixdat/pull/4
#  Perhaps best to merge [master] into [deconvolution], improve the module on the
#  latter branch, and then reopen the PR.

# changes made by Anna from the original PR by Kevin:
# *) some syntax fixes to adjust for changes made to ixdat overall (up to 0.2.8)
# *) remove the class DecoMeasurement and put the related methods directly into
# the ECMSMeasurement class. I think it's extremely confusing and inconvenient 
# if I have to import my ECMS data in a different way to perform deconvolution.
# It should be possible to perform this on any ECMS data I have.
# *) Change the name "Kernel" to "ECMSImpulseResponse". Might be the correct mathematical
# term, but I found it confusing (especially since Python also has Kernels)
# *) Make the ECMSImpulseResponse class inherit from Measurement (if I figure out how to) or even ECMSMeasurement?
# *) I don't think there should be a separation in the class whether the object is
# constructed from measured data or parameters (since the object will be the same, won't it? It's like
# the sensitivity factor will be the same object whether calculated or measured).
# instead there should be different methods to construct it.

class DecoMeasurement(ECMSMeasurement):
    """Class implementing deconvolution of EC-MS data"""

    def __init__(self, **kwargs):
        """Initialize a deconvolution EC-MS measurement

        Args:
            name (str): The name of the measurement"""
        super().__init__(**kwargs)

    def grab_partial_current(
        self, signal_name, kernel_obj, tspan=None, tspan_bg=None, snr=10
    ):
        """Return the deconvoluted partial current for a given signal

        Args:
            signal_name (str): Name of molecule for which deconvolution is to
                be carried out.
            kernel_obj (Kernel): Kernel object which contains the mass transport
                parameters
            tspan (list): Timespan for which the partial current is returned.
            tspan_bg (list): Timespan that corresponds to the background signal.
            snr (int): signal-to-noise ratio used for Wiener deconvolution.
        """
        # TODO: comments in this method so someone can tell what's going on!

        t_sig, v_sig = self.grab_flux(signal_name, tspan=tspan, tspan_bg=tspan_bg)

        kernel = kernel_obj.calculate_kernel(
            dt=t_sig[1] - t_sig[0], duration=t_sig[-1] - t_sig[0]
        )
        kernel = np.hstack((kernel, np.zeros(len(v_sig) - len(kernel))))
        H = fft(kernel)
        # TODO: store this as well.
        partial_current = np.real(
            ifft(fft(v_sig) * np.conj(H) / (H * np.conj(H) + (1 / snr) ** 2))
        )
        partial_current = partial_current * sum(kernel)
        return t_sig, partial_current

    def extract_kernel(self, signal_name, cutoff_pot=0, tspan=None, tspan_bg=None):
        """Extracts a Kernel object from a measurement.

        Args:
            signal_name (str): Signal name from which the kernel/impulse
                response is to be extracted.
            cutoff_pot (int): Potential which the defines the onset of the
                impulse. Must be larger than the resting potential before the
                impulse.
            tspan(list): Timespan from which the kernel/impulse response is
                extracted.
            tspan_bg (list): Timespan that corresponds to the background signal.
        """
        x_curr, y_curr = self.grab("current", tspan=tspan) # not sure if this still works
        x_pot, y_pot = self.grab("potential", tspan=tspan) # not sure if this still works
        x_sig, y_sig = self.grab_signal(signal_name, tspan=tspan, tspan_bg=tspan_bg)

        # TODO make this more generally applicable than these 3 masses
        if signal_name == "M32":
            t0 = x_curr[np.argmax(y_pot > cutoff_pot)]  # time of impulse
        elif signal_name == "M2" or signal_name == "M17":
            t0 = x_curr[np.argmax(y_pot < cutoff_pot)]
        else:
            print("mass not found")

        x_sig = x_sig - t0

        y_sig = y_sig[x_sig > 0]
        x_sig = x_sig[x_sig > 0]

        y_curr = y_curr[x_curr > t0]
        x_curr = x_curr[x_curr > t0]
        y_pot = y_pot[x_pot > t0]
        x_pot = x_pot[x_pot > t0]

        kernel = Kernel(
            MS_data=np.array([x_sig, y_sig]),
            EC_data=np.array([x_curr, y_curr, x_pot, y_pot]),
        )

        return kernel

class ECMSImpulseResponse:
    """
    Class implementing impulse response deconvolution of ECMS data.
    
    # TODO: Make class inherit from Measurement, add properties to store kernel
    # TODO: Reference equations to paper.
    
    """
    
    def __init__(
        self,
        mol,
        data=None,
        working_distance=None,
        diff_const=None,
        henry_vola=None,
        chip=None
    ):
        """Initializes a Kernel object either in functional form by defining the
        mass transport parameters or in the measured form by passing of EC-MS
        data.

        Args:
            mol: Molecule to calculate the impulse response of
            data: ECMSMeasurement object. Optional. If passed, the impulse response
            will be calculated based on the measured data, overwriting the parameters
            passed for a calculated one.
            work_dist: Working distance between electrode and gas/liq interface. Optional, 
            though necessary if no data is provided.
            diff_const: Diffusion constant in liquid. Optional. Default will check
            diffusion constant in water in Molecule data from si_quant.
            henry_vola: Dimensionless Henry volatility. Optional. Default will check
            Henry volatility constant in water in Molecule data from si_quant.
            chip: Optional. Needed to define capillary flow. Default will use
            SpectroInlets chip from si_quant.
        """

        if data is None and working_distance is None:
            raise DeconvolutionError # TODO need to make sure this error type exists
            
        if data is not None:            
            self.mol = mol
            self.data = data
            print("Generating ECMSImpulseResponse from measured data.")
            self.type = "measured"
            if working_distance is not None:
                raise UserWarning("Data was used to generate ECMSImpulseResponse ignoring the given working_distance.")
        else:
            self.mol = mol
            self.type = "functional"
            
            # find the Molecule parameters 
            # chip properties: molecule flux through capillary, and also the size of the sample volume.
            
            # dont I need the sample area for something as well??
            
            self.params = {working_distance: working_distance, diff_const: diff_const,
                           henry_vola: henry_vola, }       
            
            raise Exception(
                "Kernel can only be initialized with data OR parameters, not both"
            )
      
    
    def calculate_kernel(self, dt=0.1, duration=100, norm=True, matrix=False):
        """Calculates a kernel/impulse response.

        Args:
            dt (int): Timestep for which the kernel/impulse response is calculated.
                Has to match the timestep of the measured data for deconvolution.
            duration(int): Duration in seconds for which the kernel/impulse response is
                calculated. Must be long enough to reach zero.
            norm (bool): If true the kernel/impulse response is normalized to its
                area.
            matrix (bool): If true the circulant matrix constructed from the kernel/
                impulse reponse is returned.
        """
        if self.type == "functional":

            t_kernel = np.arange(0, duration, dt)
            t_kernel[0] = 1e-6

            diff_const = self.params["diff_const"]
            work_dist = self.params["work_dist"]
            vol_gas = self.params["vol_gas"]
            volflow_cap = self.params["volflow_cap"]
            henry_vola = self.params["henry_vola"]
            el_A = self.params["el_A"]

            tdiff = t_kernel * diff_const / (work_dist**2)

            def fs(s):
                # See Krempl et al, 2021. Equation 6.
                #     https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
                return 1 / (
                    sqrt(s) * sinh(sqrt(s))
                    + (vol_gas * henry_vola / (el_A * 1e-4 * work_dist))
                    * (s + volflow_cap / vol_gas * work_dist**2 / diff_const)
                    * cosh(sqrt(s))
                )

            kernel = np.zeros(len(t_kernel))
            for i in range(len(t_kernel)):
                kernel[i] = invertlaplace(fs, tdiff[i], method="talbot")
                print(tdiff[i])
                print(kernel[i])

        elif self.type == "measured":
            kernel = self.MS_data[1]
            t_kernel = self.MS_data[0]

        if norm:
            area = np.trapz(kernel, t_kernel)
            kernel = kernel / area

        if matrix:
            kernel = np.tile(kernel, (len(kernel), 1))
            i = 1
            while i < len(t_kernel):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i = i + 1

        return kernel
    
    @property
    def sig_area(self):
        """Integrates a measured impulse response and returns the area."""
        delta_sig = self.MS_data[1] - self.MS_data[1][-1]
        sig_area = np.trapz(delta_sig, self.MS_data[0])

        return sig_area

    @property
    def charge(self):
        """Integrates the measured current over the time."""
        y_curr = self.EC_data[1]

        mask = np.isclose(y_curr, y_curr[0], rtol=1e-1)

        Q = np.trapz(y_curr[mask], self.EC_data[0][mask])

        return Q
        
        

class Kernel:
    """Kernel class implementing datatreatment of kernel/impulse response data."""

    # TODO: Make class inherit from Measurement, add properties to store kernel
    # TODO: Reference equations to paper.
    def __init__(
        self,
        parameters={},  # FIXME: no mutable default arguments!
        MS_data=None,
        EC_data=None,
    ):
        """Initializes a Kernel object either in functional form by defining the
        mass transport parameters or in the measured form by passing of EC-MS
        data.

        Args:
            parameters (dict): Dictionary containing the mass transport
                parameters with the following keys:
                    diff_const: Diffusion constant in liquid
                    work_dist: Working distance between electrode and gas/liq interface
                    vol_gas: Gas sampling volume of the chip
                    volflow_cap: Volumetric capillary flow
                    henry_vola: Dimensionless Henry volatility
                MS_data (list): List of numpy arrays containing the MS signal
                    data.
                EC_data (list): List of numpy arrays containing the EC (time,
                    current, potential).
        """

        if MS_data and parameters:  # TODO: Make two different classes
            raise Exception(
                "Kernel can only be initialized with data OR parameters, not both"
            )
        if EC_data and MS_data:
            print("Generating kernel from measured data")
            self.type = "measured"
        elif parameters:
            print("Generating kernel from parameters")
            self.type = "functional"
        else:
            print("Generating blank kernel")
            self.type = None

        self.params = parameters
        self.MS_data = MS_data
        self.EC_data = EC_data  # x_curr, y_curr, x_pot, y_pot

    @property
    def sig_area(self):
        """Integrates a measured impulse response and returns the area."""
        delta_sig = self.MS_data[1] - self.MS_data[1][-1]
        sig_area = np.trapz(delta_sig, self.MS_data[0])

        return sig_area

    @property
    def charge(self):
        """Integrates the measured current over the time."""
        y_curr = self.EC_data[1]

        mask = np.isclose(y_curr, y_curr[0], rtol=1e-1)

        Q = np.trapz(y_curr[mask], self.EC_data[0][mask])

        return Q

    def plot(self, dt=0.1, duration=100, ax=None, norm=True, **kwargs): 
        
        """Returns a plot of the kernel/impulse response.
        
        TODO: this shouldn't be it's separate function but rather use the ECMS plotter!?
        """
        if ax is None:
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)

        if self.type == "functional":
            t_kernel = np.arange(0, duration, dt)
            ax.plot(
                t_kernel,
                self.calculate_kernel(dt=dt, duration=duration, norm=norm),
                **kwargs,
            )

        elif self.type == "measured":
            ax.plot(
                self.MS_data[0],
                self.calculate_kernel(dt=dt, duration=duration, norm=norm),
                **kwargs,
            )

        else:
            raise Exception("Nothing to plot with blank kernel")

        return ax

    def calculate_kernel(self, dt=0.1, duration=100, norm=True, matrix=False):
        """Calculates a kernel/impulse response.

        Args:
            dt (int): Timestep for which the kernel/impulse response is calculated.
                Has to match the timestep of the measured data for deconvolution.
            duration(int): Duration in seconds for which the kernel/impulse response is
                calculated. Must be long enough to reach zero.
            norm (bool): If true the kernel/impulse response is normalized to its
                area.
            matrix (bool): If true the circulant matrix constructed from the kernel/
                impulse reponse is returned.
        """
        if self.type == "functional":

            t_kernel = np.arange(0, duration, dt)
            t_kernel[0] = 1e-6

            diff_const = self.params["diff_const"]
            work_dist = self.params["work_dist"]
            vol_gas = self.params["vol_gas"]
            volflow_cap = self.params["volflow_cap"]
            henry_vola = self.params["henry_vola"]
            el_A = self.params["el_A"]

            tdiff = t_kernel * diff_const / (work_dist**2)

            def fs(s):
                # See Krempl et al, 2021. Equation 6.
                #     https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
                return 1 / (
                    sqrt(s) * sinh(sqrt(s))
                    + (vol_gas * henry_vola / (el_A * 1e-4 * work_dist))
                    * (s + volflow_cap / vol_gas * work_dist**2 / diff_const)
                    * cosh(sqrt(s))
                )

            kernel = np.zeros(len(t_kernel))
            for i in range(len(t_kernel)):
                kernel[i] = invertlaplace(fs, tdiff[i], method="talbot")
                print(tdiff[i])
                print(kernel[i])

        elif self.type == "measured":
            kernel = self.MS_data[1]
            t_kernel = self.MS_data[0]

        if norm:
            area = np.trapz(kernel, t_kernel)
            kernel = kernel / area

        if matrix:
            kernel = np.tile(kernel, (len(kernel), 1))
            i = 1
            while i < len(t_kernel):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i = i + 1

        return kernel
