"""Module for analysis of EC-MS measurements with impulse response and
deconvolution of mass transport effects to obtain partial current densities."""

from .ec_ms import ECMSMeasurement
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import signal
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft, ifftshift, fftfreq
import numpy as np


class DecoMeasurement(ECMSMeasurement):
    """Class implementing deconvolution of EC-MS data"""

    def __intit__(self, name, **kwargs):
        """initialize a deconvolution EC-MS measurement

        Args:
            name (str): The name of the measurement"""
        super().__init__(name, **kwargs)

    def get_partial_current(
        self, signal_name, kernel_obj, tspan=None, t_bg=None, snr=10
    ):
        """Return the deconvoluted partial current for a given signal

        Args:
            signal_name (str): Name of signal for which deconvolution is to
                be carried out.
            kernel_obj (Kernel): Kernel object which contains the mass transport
                parameters
            tspan (list): Timespan for which the partial current is returned.
            t_bg (list): Timespan that corresponds to the background signal.
            snr (int): signal-to-noise ratio used for Wiener deconvolution.
        """

        t_sig, v_sig = self.get_calib_signal(signal_name, tspan=tspan, t_bg=t_bg)

        kernel = kernel_obj.calculate_kernel(dt=t_sig[1] - t_sig[0])
        kernel = np.hstack((kernel, np.zeros(len(v_sig) - len(kernel))))
        H = fft(kernel)
        # TODO: store this as well.
        partial_current = np.real(
            ifft(fft(v_sig) * np.conj(H) / (H * np.conj(H) + (1 / snr) ** 2))
        )
        partial_current = partial_current * sum(kernel)
        return t_sig, partial_current

    def extract_kernel(self, signal_name, cutoff_pot=0, tspan=None, t_bg=None):
        """Extracts a Kernel object from a measurement.

        Args:
            signal_name (str): Signal name from which the kernel/impule
                response is to be extracted.
            cutoff_pot (int): Potential which the defines the onset of the
                impulse. Must be larger than the resting potential before the
                impulse.
            tspan(list): Timespan from which the kernel/impulse response is
                extracted.
            t_bg (list): Timespan that corresponds to the background signal.
        """
        x_curr, y_curr = self.get_current(tspan=tspan)
        x_pot, y_pot = self.get_potential(tspan=tspan)
        x_sig, y_sig = self.get_signal(signal_name, tspan=tspan, t_bg=t_bg)

        if signal_name == "M32":
            t0 = x_curr[np.argmax(y_pot > cutoff_pot)]  # time of impulse
        elif signal_name == "M2":
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


class Kernel:
    """Kernel class implementing datatreatment of kernel/impulse response data."""

    # TODO: Make class inherit from Measurement, add properties to store kernel
    # TODO: Reference equations to paper.
    def __init__(
        self,
        parameters={},
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

        if MS_data is not None and parameters:  # TODO: Make two different classes
            raise Exception(
                "Kernel can only be initialized with data OR parameters, not both"
            )
        if EC_data is not None and MS_data is not None:
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
        """Returns a plot of the kernel/impulse response."""
        if ax is None:
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)

        if self.type is "functional":
            t_kernel = np.arange(0, duration, dt)
            ax.plot(
                t_kernel,
                self.calculate_kernel(dt=dt, duration=duration, norm=norm),
                **kwargs,
            )

        elif self.type is "measured":
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
        if self.type is "functional":

            t_kernel = np.arange(0, duration, dt)
            t_kernel[0] = 1e-6

            diff_const = self.params["diff_const"]
            work_dist = self.params["work_dist"]
            vol_gas = self.params["vol_gas"]
            volflow_cap = self.params["volflow_cap"]
            henry_vola = self.params["henry_vola"]

            tdiff = t_kernel * diff_const / (work_dist ** 2)
            fs = lambda s: 1 / (
                sqrt(s) * sinh(sqrt(s))
                + (vol_gas * henry_vola / 0.196e-4 / work_dist)
                * (s + volflow_cap / vol_gas * work_dist ** 2 / diff_const)
                * cosh(sqrt(s))
            )

            kernel = np.zeros(len(t_kernel))
            for i in range(len(t_kernel)):
                kernel[i] = invertlaplace(fs, tdiff[i], method="talbot")

        elif self.type is "measured":
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
