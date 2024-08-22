"""Module for deconvolution of mass transport effects."""

import numpy as np
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
# a certain analyte. see Krempl et al. 2021
# https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
# for definition of condition number

# TODO: rename kernel to something else


class ECMSImpulseResponse:
    """
    Class for signal impulse response for ECMS data for modelling diffusion.
    see Krempl et al. 2021
    https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110

    This class currently handles only one molecule at a time. If deconvolution of
    multiple molecules is required this needs to be handled by separate objects.

    # TODO: Make class inherit from Calculator
    """

    def __init__(
        self,
        mol,
        t_kernel,
        kernel,
        measurement=None,
        working_distance=None,
        A_el=None,
        D=None,
        H_v_cc=None,
        n_dot=None,
        T=None,
        p=None,
        carrier_gas=None,
        gas_volume=None,
        dt=None,
        duration=None,
        ir_type=None,
    ):
        """
        Initializes a ECMSImpulseResponse object either in functional form by defining
        the mass transport parameters or in the measured form by passing of EC-MS
        measurement.

        Args:
            mol (str): Molecule to calculate the impulse response of.
            t_kernel (array): Impulse response time series
            kernel (array): Impulse response value series
            measurement (ECMSMeasurement): Measurement including impulse response
                measurment. Optional. Will be included if generated from measurement.
            working_distance (float): Working distance between electrode and gas/liq
                interface in [m]. Optional.  Will be included if generated from
                parameters.
            A_el (float): Geometric electrode area in [cm^2]. Optional. Will
                automatically be included if generated from parameters.
            D (float): Diffusion constant in liquid in [m2/s]. Optional. Will
                automatically be included if generated from parameters.
            H_v_cc (float): Dimensionless Henry volatility. Optional. Will automatically
                be included if generated from parameters.
            n_dot (float): Capillary flux in [mol/s]. Optional. Will automatically be
                included if generated from parameters.
            T (float): Temperature in K. Optional. Will automatically be included if
                generated from parameters.
            p (float): Pressure (on the high pressure side) in Pa. Optional. Will
                automatically be included if generated from parameters.
            carrier_gas: The carrier gas used to calculate capillary flow. Optional.
                Will automatically be included if generated from parameters AND using
                siq.
            gas_volume (float): the volume of the headspace volume in the chip. Default
                is the volume of the SpectroInlets chip.
            dt (float): Timestep for which the impulse response is calculated in [s].
                Has to match the timestep of the measurement for deconvolution.
            duration (float): Duration in [s] for which the kernel/impulse response is
                calculated. Must be long enough to reach zero.
        """
        self.mol = mol
        self.t_kernel = t_kernel
        self.kernel = kernel
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
        self.ir_type = ir_type
        self.dt = dt
        self.duration = duration

    @classmethod
    def from_measurement(
        cls,
        mol,
        measurement,
        tspan=None,
        tspan_bg=None,
        norm=True,
        matrix=False,
        **kwargs,
    ):
        """
        Generate an ECMSImpulseResponse from measurement.

        Args:
            mol (str): Molecule to calculate the impulse response of
            measurement (ECMSMeasurement): Measurement including impulse response
                measurment. Measurement needs to contain calibration data for mol
                (either using ixdat native ECMSMeasurement.calibrate() or external
                 package (eg siq's "quantifier"))
            tspan (list): tspan over which to calculate the impulse response. Needs to
                include zero.
            tspan_bg (list): tspan of background to subtract. If list of tspans
                (list of lists) is passed, will interpolate between the points using
                calc_linear_background()
            norm (bool): If true the impulse response is normalized to its
                area. Default is True.

            Additional keyword arguments are passed on to ECMSImpulseResponse.__init__
        Return ECMSImpulseResponse
        """
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
            area = np.trapz(kernel, t_kernel)
            kernel = kernel / area
        return cls(mol, t_kernel, kernel, ir_type="measured", **kwargs)

    @classmethod
    def from_parameters(
        cls,
        mol,
        working_distance,
        A_el=1,
        D=None,
        H_v_cc=None,
        n_dot=None,
        T=STANDARD_TEMPERATURE,
        p=STANDARD_PRESSURE,
        carrier_gas="He",
        gas_volume=1e-10,
        dt=0.1,
        duration=100,
        norm=True,
        matrix=False,
        **kwargs,
    ):
        """Generate an ECMSImpulseResponse from parameters.

            All references to equations refer to the Supporting Information of Krempl
            et al 2021: https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110

            Requires activation of siq for optimal performance.

        Args:
            mol (str): Molecule to calculate the impulse response of.
            working_distance (float): Working distance between electrode and gas/liq
                interface in [m].
            A_el (float): Geometric electrode area in [cm^2]. Defaults to 1.
                #TODO is that a good idea?
            D (float): Diffusion constant in liquid in [m2/s]. Optional if using siq,
                will then be looked up in molecule.yml file
            H_v_cc (float): Dimensionless Henry volatility. Optional if using siq, will
                then be looked up in molecule.yml file
            n_dot (float): Capillary flux in [mol/s]. Optional if using siq, will then
                be calculated based on the carrier_gas selected.
            T (float): Temperature in K. Optional. Defaults to STANDARD_TEMPERATURE
                (298.15K)
            p (float): Pressure (on the high pressure side) in Pa. Optional. Defaults to
                STANDARD_PRESSURE (1e5 Pa).
            carrier_gas (str): The carrier gas used to calculate capillary flow. Only
                used when using siq. Defaults to He.
            gas_volume (float): The volume of the headspace volume in the chip in [m3].
                Defaults to the volume of the SpectroInlets chip.
            dt (float): Timestep for which the impulse response is calculated in [s].
                Has to match the timestep of the measurement for deconvolution.
            duration (float): Duration in [s] for which the kernel/impulse response is
                calculated. Must be long enough to reach zero.
            norm (bool): If true the impulse response is normalized to its
                area. Default is True.

            Additional keyword arguments are passed on to ECMSImpulseResponse.__init__
         Return ECMSImpulseResponse
        """

        if plugins.use_siq:
            if n_dot is None:
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
        # convert the mol flux to volumetric flux using T and p given by chip.
        V_dot = n_dot * R * T / p

        # make a time array
        t_kernel = np.arange(0, duration, dt)
        t_kernel[0] = 1e-6

        tdiff = (
            t_kernel * D / (working_distance**2)
        )  # See Krempl et al, 2021. SI: Equ. (6)

        def calc_H(s):
            # See Krempl et al, 2021. SI: Equ 24 + definitions for kappa and lamda
            # (between Equ 11 and 12)
            kappa = (A_el * 1e-4 * working_distance) / gas_volume

            lambda_ = V_dot / gas_volume * working_distance**2 / D

            H = 1 / (
                H_v_cc * (s + lambda_) / kappa * cosh(sqrt(s)) + sqrt(s) * sinh(sqrt(s))
            )
            return H

        # now calculate the modelled impulse response: See Krempl et al, 2021. SI, page
        # S4: "calculating the discrete values of the impulse response at a sampling
        # frequency equal to the sampling frequency of the measured mass spectrometer
        # signal by numerical inverse Laplace transformation with the Talbot method"
        kernel = np.zeros(len(t_kernel))
        for i in range(len(t_kernel)):
            kernel[i] = invertlaplace(calc_H, tdiff[i], method="talbot")

        if (
            norm
        ):  # normalize the kernel intensity to the total area under the ImpulseResponse
            area = np.trapz(kernel, t_kernel)
            kernel = kernel / area

        # return cls(mol, t_kernel, kernel, ir_type="calculated", **kwargs)
        # the kwargs thing doesnt work somehow
        return cls(
            mol,
            t_kernel,
            kernel,
            working_distance=working_distance,
            A_el=A_el,
            D=D,
            H_v_cc=H_v_cc,
            n_dot=n_dot,
            T=T,
            p=p,
            carrier_gas=carrier_gas,
            gas_volume=gas_volume,
            dt=dt,
            duration=duration,
            ir_type="calculated",
            **kwargs,
        )


# TODO: other potentially useful methods:
# https://github.com/ixdat/ixdat/blob/f577a434a966e486cf4cb66253677f7839fe117a/src/
# ixdat/techniques/deconvolution.py#L242
