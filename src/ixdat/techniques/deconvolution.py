"""Module for deconvolution of mass transport effects."""

import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  # noqa
from scipy.interpolate import interp1d  # noqa
from scipy import signal  # noqa
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth  # noqa
from numpy.fft import fft, ifft, ifftshift, fftfreq  # noqa

from ..exceptions import TechniqueError
from ..config import plugins
from ..constants import R, STANDARD_TEMPERATURE, STANDARD_PRESSURE
from ..plotters.plotting_tools import calc_linear_background

# TODO: Would be very useful to have a method to calculate the condition number for
# a certain analyte

class ECMSImpulseResponse:
    """
    Class implementing impulse response deconvolution of ECMS measurement.
    
    Extracts an ECMSImpulseResponse object from a measurement using the
    algorithm developed by Krempl et al. 2021
    https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
    
    some args for initializing are optional if a 
    
    TODO: add something about the requirements of the ECMSMeasurement. eg t_zero needs to 
    be where the current starts (does it?), it needs to be quantified data! only one mol at a time

    # TODO: Make class inherit from Calculator, add properties to store kernel
    # TODO: Reference equations to paper.
    """

    def __init__(
        self,
        mol,
        measurement=None,
        working_distance=None,
        A_el=None,
        D=None,
        H_v_cc=None,
        n_dot=None,
        T = STANDARD_TEMPERATURE,
        p = STANDARD_PRESSURE,
        carrier_gas=None,
        gas_volume=1e-10,
    ):
        """
        Initializes a ECMSImpulseResponse object either in functional form by defining
        the mass transport parameters or in the measured form by passing of EC-MS measurement.
        
        Args:
            mol: Molecule to calculate the impulse response of
            measurement: ECMSMeasurement object. Optional. If passed, the impulse response
                will be calculated based on the measurement, overwriting the parameters
                passed for a calculated one.
            working_distance: Working distance between electrode and gas/liq interface in
                m. Optional, though necessary if no data is provided.
            A_el: Geometric electrode area in cm2. Optional, though necessary if no data
                is provided.
            D: Diffusion constant in liquid. Optional. Default will check
            diffusion constant in water in Molecule data from siq.
            H_v_cc: Dimensionless Henry volatility. Optional. Default will check
            Henry volatility constant in water in Molecule data from siq.
            n_dot: Capillary flux in mol/s. Optional if siq is activated.
            T: Temperature in K. Defaults to STANDARD_TEMPERATURE (298.15 K)
            p: Pressure (on the high pressure side) in Pa. Defaults to STANDARD_PRESSURE (100000 Pa)
            carrier_gas: The carrier gas used to calculate capillary flow when using siq. Defaults to He
            gas_volume: the volume of the headspace volume in the chip. Default is the
                volume of the SpectroInlets chip.
        """
        self.mol = mol # Johannes said that might not be needed anymore, assigning all the selfs here?
        self.measurement = measurement
        self.working_distance = working_distance
        self.A_el = A_el
        self.D = D
        self.H_v_cc = H_v_cc
        self.n_dot = n_dot
        self.T = T
        self.p = p
        self.carrier_gas = carrier_gas
        self.gas_volume = gas_volume

    @classmethod
    def from_measurement(
            cls,
            mol,
            measurement=None,
            A_el=None,
            T=STANDARD_TEMPERATURE,
            p=STANDARD_PRESSURE,
            tspan=None,
            tspan_bg=None,
            norm=True,
            matrix=False
            ):
        """Generate an ECMSImpulseResponse from measurement.
        TODO: add option to plot the implulse response from measurement/ return an axis to
            co-plot with the measurement

        Args:
            tspan (list): tspan over which to calculate the impulse response
            tspan_bg (list): tspan of background to subtract. if list of tspans is
                passed, will interpolate between the points using
                calc_linear_background()
            norm (bool): If true the impulse response is normalized to its
                area. Default is True.
            matrix (bool): If true the circulant matrix constructed from the
                impulse reponse is returned. Default is False.
        Return t_kernel (np.array), kernel (np/.array): tuple of time and value arrays of the measured impulse response 
        """
        cls.type = "measured"
        print("Generating `ECMSImpulseResponse` from measurement.")
        if type(tspan_bg[0]) is not list:
            t_kernel, kernel = measurement.grab_flux(
                mol=mol, tspan=tspan, tspan_bg=tspan_bg
            )
        else:
            t_kernel, kernel_raw = measurement.grab_flux(mol=mol, tspan=tspan)
            bg = calc_linear_background(t_kernel, kernel_raw, tspans=tspan_bg)
            kernel = kernel_raw - bg
        if norm:
            area = np.trapezoid(kernel, t_kernel)
            kernel = kernel / area
        if matrix:
            kernel = np.tile(kernel, (len(kernel), 1))
            i = 1
            while i < len(t_kernel):  # TODO: pythonize this?
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i = i + 1
        return t_kernel, kernel # can this return more than just this tuple of arrays? I.e. some meta info as defined in the class?
    # I want this to return something like the quantification does. Maybe wait for calculators.
            
    @classmethod
    def from_parameters(
            cls,
            mol,
            measurement=None,
            working_distance=None,
            A_el=None,
            D=None,
            H_v_cc=None,
            n_dot=None,
            T = STANDARD_TEMPERATURE,
            p = STANDARD_PRESSURE,
            carrier_gas=None,
            gas_volume=1e-10,
            dt=0.1, # TODO need to think about how to do this when extracting deconvoluted data
            # as that requires recalculation of this to fit exactly to the length of data. 
            # some kind of shorthand for an object containing all the info needed except
            # these two would be really helpful, but this has gotten lost with this classmethod
            # stuff
            duration=100,
            norm=True,
            matrix=False
            ):
        """Generate an ECMSImpulseResponse from parameters.
        Args:
            dt (int): Timestep for which the impulse response is calculated.
                Has to match the timestep of the measurement for deconvolution.

            duration(int): Duration in seconds for which the kernel/impulse response is
                calculated. Must be long enough to reach zero.
            norm (bool): If true the impulse response is normalized to its
                area.
            matrix (bool): If true the circulant matrix constructed from the
                impulse reponse is returned.
        Return t_kernel (np.array), kernel (np/.array): tuple of time and value arrays of the modelled impulse response 
        
        """
        cls.type = "functional"
        if plugins.use_siq:
            # calculate the capillary flow for the specified gas & chip
            # initiate the chip
            Chip = plugins.siq.Chip
            chip = Chip()
            # set the T and p of the chip according to the values given
            chip.T = T
            chip.p = p
            n_dot = chip.calc_n_dot_0(gas=carrier_gas)
            # find the other parameters from the siq Molecule files
            Molecule = plugins.siq.Molecule
            molecule = Molecule.load(mol)
            if D is None:
                D = molecule.D
                # TODO double check units and understand why siq integration is not
                # working - maybe because I initiate the molecule wrongly! 
                print(D)
            if H_v_cc is None:
                H_v_cc = 1 / (R * chip.T * molecule.H_0) * 100
                print(H_v_cc)
        else:
            if D is None:
               raise TechniqueError(
                   "D is required to initialize ECMSImpluseResponse. Alternatively" 
                   "activate siq as plugin to automatically load D for molecule."
                   ) 
            if H_v_cc is None:
                raise TechniqueError(
                    "H_v_cc is required to initialize ECMSImpluseResponse. Alternatively" 
                    "activate siq as plugin to automatically load H_v_cc for molecule."
                    ) 
        # convert to volumetric flux using T and p given by chip.
        V_dot = n_dot * R * T / p
        # calculate the residence time in the chip headspace
        residence_t = V_dot / gas_volume                
        
        t_kernel = np.arange(0, duration, dt)
        t_kernel[0] = 1e-6
        
        tdiff = t_kernel * D / (working_distance**2)

        def fs(s):
            # See Krempl et al, 2021. Equation 6.
            #     https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
            return 1 / (
                sqrt(s) * sinh(sqrt(s))
                + (gas_volume * H_v_cc / (A_el * 1e-4 *  working_distance))
                * (s + residence_t *  working_distance**2 / D)
                * cosh(sqrt(s))
            )

        kernel = np.zeros(len(t_kernel))
        for i in range(len(t_kernel)):
            kernel[i] = invertlaplace(fs, tdiff[i], method="talbot")

        if norm:
            area = np.trapz(kernel, t_kernel)
            kernel = kernel / area
        if matrix:
            kernel = np.tile(kernel, (len(kernel), 1))
            i = 1
            while i < len(t_kernel):  # TODO: pythonize this?
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i = i + 1
            
        return t_kernel, kernel

    @np.vectorize
    def get_cond_number(sampling_freq, working_dist):
        """
        Function to calculate condition number for a given sampling frequency
        and working distance.
        TODO: finish this method
        """

        # Define mass transport parameters.
        params = {
            "D": 5.05e-9,
            "work_dist": working_dist * 1e-6,
            "vol_gas": 1.37e-10,
            "volflow_cap": 1e-10,
            "H_v_cc": 52,
        }

        model_kernel = ECMSImpulseResponse(parameters=params)

        kernel_matrix = model_kernel.calculate_kernel(
            dt=1 / sampling_freq, duration=40, matrix=True
        )

        cond = np.linalg.cond(kernel_matrix)
        print(str(cond))
        return cond

    def get_limit_frequency():
        """Calculate the limiting frequency for deconvolution for a given
        ECMSImpulseResponse.
        TODO: finish this method"""

    def get_limit_time_res():
        """Calculate the limiting time resolution for deconvolution for a given
        ImpulseResponse.
        TODO: finish this method"""
