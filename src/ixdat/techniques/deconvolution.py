"""Module for deconvolution of mass transport effects."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  # noqa
from scipy.interpolate import interp1d  # noqa
from scipy import signal  # noqa
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth  # noqa
from numpy.fft import fft, ifft, ifftshift, fftfreq  # noqa

from ..exceptions import TechniqueError
from ..config import plugins
from ..constants import R
from ..plotters.plotting_tools import calc_linear_background

# FIXME: too much abbreviation in this module.
# TODO: Implement the PR review here: https://github.com/ixdat/ixdat/pull/4
# TODO: Would be very useful to have a method to calculate the condition number for
# a certain analyte


class ECMSImpulseResponse:
    """
    Class implementing impulse response deconvolution of ECMS data.

    # TODO: Make class inherit from Calculator, add properties to store kernel
    # TODO: Reference equations to paper.

    """

    def __init__(
        self,
        mol,
        data=None,
        working_distance=None,
        A_el=None,
        diff_const=None,
        henry_vola=None,
        chip=None,
        carrier_gas=None,
        gas_volume=1e-10,
    ):
        """
        Initializes a ECMSImpulseResponse object either in functional form by defining
        the mass transport parameters or in the measured form by passing of EC-MS data.
        Args:
            mol: Molecule to calculate the impulse response of
            data: ECMSMeasurement object. Optional. If passed, the impulse response
                will be calculated based on the measured data, overwriting the parameters
                passed for a calculated one.
            working_distance: Working distance between electrode and gas/liq interface in
                m. Optional, though necessary if no data is provided.
            A_el: Geometric electrode area in cm2. Optional, though necessary if no data
                is provided.
            diff_const: Diffusion constant in liquid. Optional. Default will check
            diffusion constant in water in Molecule data from siq.
            henry_vola: Dimensionless Henry volatility. Optional. Default will check
            Henry volatility constant in water in Molecule data from siq.-CURRENTLY WRONG
            chip: Optional. Needed to define capillary flow. Default will use
                SpectroInlets chip from siq.
            carrier_gas: The carrier gas used to calculate capillary flow. Defaults to He
            gas_volume: the volume of the headspace volume in the chip. Default is the
                volume of the SpectroInlets chip.
        """
        if data is not None:
            self.mol = mol
            self.data = data
            self.type = "measured"
            print("Generating `ECMSImpulseResponse` from measured data.")
            if working_distance is not None or A_el is not None:
                raise UserWarning(
                    "Data was used to generate `ECMSImpulseResponse` ignoring the given"
                    "working_distance/electrode area."
                )

        elif working_distance is not None and A_el is not None:
            self.mol = mol
            self.type = "functional"
            if not plugins.use_siq:
                raise TechniqueError(
                    "`ECMSImpulseResponse` will only work properly when using "
                    # TODO this should be improved.- doesnt need to fully depend on siq
                    # integration
                    "`spectro_inlets_quantification` "
                    "(`ixdat.plugins.activate_siq()`). "
                )
            # calculate the capillary flow for the specified gas & chip
            Chip = plugins.siq.Chip
            chip = chip or Chip()
            n_dot = chip.calc_n_dot_0(gas=carrier_gas)
            # convert to volumetric flux using T and p given by chip.
            V_dot = n_dot * R * chip.T / chip.p
            # find the other parameters from the siq Molecule files
            Molecule = plugins.siq.Molecule
            molecule = Molecule(mol)
            # print(molecule)
            if diff_const is None:
                raise TechniqueError("Default diffusion constant not yet implemented")
                diff_const = molecule.D
                # TODO double check units and understand why siq integration is not
                # working
                print(diff_const)
            if henry_vola is None:
                raise TechniqueError("Default henry volatility not yet implemented")
                # henry_vola = molecule.H_0
                # TODO convert this to the right unit! (dimensionless henry volatility)
            self.params = {
                "working_distance": working_distance,
                "A_el": A_el,
                "diff_const": diff_const,
                "henry_vola": henry_vola,
                "V_dot": V_dot,
                "gas_volume": gas_volume,
            }
        else:
            raise TechniqueError(
                "Cannot initialize ECMSImpluseResponse without either data or working"
                "distance + electrode area being provided."
            )

    def model_impulse_response_from_params(
        self, dt=0.1, duration=100, norm=True, matrix=False
    ):
        """Calculates an impulse response from parameters used to initialize the
        ECMSImpulseResponse object.
        TODO: might make more sense to pass parameters to the method instead?
        Args:
            dt (int): Timestep for which the impulse response is calculated.
                Has to match the timestep of the measured data for deconvolution.

            duration(int): Duration in seconds for which the kernel/impulse response is
                calculated. Must be long enough to reach zero.
            norm (bool): If true the impulse response is normalized to its
                area.
            matrix (bool): If true the circulant matrix constructed from the
                impulse reponse is returned.
        """
        if self.type == "functional":

            t_kernel = np.arange(0, duration, dt)
            t_kernel[0] = 1e-6
            diff_const = self.params["diff_const"]
            work_dist = self.params["working_distance"]
            vol_gas = self.params["gas_volume"]
            volflow_cap = self.params["V_dot"]
            henry_vola = self.params["henry_vola"]
            A_el = self.params["A_el"]

            tdiff = t_kernel * diff_const / (work_dist**2)

            def fs(s):
                # See Krempl et al, 2021. Equation 6.
                #     https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
                return 1 / (
                    sqrt(s) * sinh(sqrt(s))
                    + (vol_gas * henry_vola / (A_el * 1e-4 * work_dist))
                    * (s + volflow_cap / vol_gas * work_dist**2 / diff_const)
                    * cosh(sqrt(s))
                )

            kernel = np.zeros(len(t_kernel))
            for i in range(len(t_kernel)):
                kernel[i] = invertlaplace(fs, tdiff[i], method="talbot")
                # print(tdiff[i])
                # print(kernel[i])
            if norm:
                area = np.trapz(kernel, t_kernel)
                kernel = kernel / area
            if matrix:
                kernel = np.tile(kernel, (len(kernel), 1))
                i = 1
                while i < len(t_kernel):  # TODO: pythonize this?
                    kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                    i = i + 1
        else:
            raise TechniqueError(
                "Cannot model impulse response if not initialized with parameters."
            )
        return t_kernel, kernel

    def calc_impulse_response_from_data(
        self, dt=0.1, duration=100, tspan=None, tspan_bg=None, norm=True, matrix=False
    ):
        """Calculates impulse response from data.
        TODO: add possibility of subtracting a linear background!!!
        TODO: add option to plot the implulse response from data/ return an axis to
            co-plot with the data
        TODO: figure out if it's ok to just get rid of dt and duration (not used here)
        TODO: change the name kernel to something else

        Args:
            tspan (list): tspan over which to calculate the impulse response
            tspan_bg (list): tspan of background to subtract. if list of tspans is
                passed, will interpolate between the points using
                calc_linear_background()
            norm (bool): If true the impulse response is normalized to its
                area. Default is True.
            matrix (bool): If true the circulant matrix constructed from the
                impulse reponse is returned. Default is False.
        """
        if self.type == "measured":
            if type(tspan_bg[0]) is not list:
                t_kernel, kernel = self.data.grab_flux(
                    mol=self.mol, tspan=tspan, tspan_bg=tspan_bg
                )
            else:
                t_kernel, kernel_raw = self.data.grab_flux(mol=self.mol, tspan=tspan)
                bg = calc_linear_background(t_kernel, kernel_raw, tspans=tspan_bg)
                kernel = kernel_raw - bg
            if norm:
                area = np.trapz(kernel, t_kernel)
                kernel = kernel / area
            if matrix:
                kernel = np.tile(kernel, (len(kernel), 1))
                i = 1
                while i < len(t_kernel):  # TODO: pythonize this?
                    kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                    i = i + 1
        else:
            raise TechniqueError("Cannot calculate impulse response without data.")
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
            "diff_const": 5.05e-9,
            "work_dist": working_dist * 1e-6,
            "vol_gas": 1.37e-10,
            "volflow_cap": 1e-10,
            "henry_vola": 52,
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

    def plot_cond_number_heatmap(cond_number, Z, cc):
        """
        TODO: finish this method

        Parameters
        ----------
        cond_number : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        # Create figure
        fig = plt.figure(figsize=(7.5 / 2, 3.5))
        ax = fig.add_subplot(111)

        # Construct meshgrid for heatmap
        x_sampling_freq = np.logspace(-1, 1, 100)  # Sampling frequency
        y_working_dist = np.linspace(100, 200, 50)  # Working distance
        X, Y = np.meshgrid(x_sampling_freq, y_working_dist)

        # Create and format heatmap plot
        Z[Z > 100] = 100  # can't get `vmax=100` colorbar to work right otherwise.
        cp = plt.contourf(X, Y, Z, 100, cmap=cc.cm.rainbow)
        # for c in cp.collections:
        #     c.set_edgecolor("face")
        cbar = plt.colorbar(cp, ticks=[1, 20, 40, 60, 80, 100])
        cbar.ax.set_ylim([0, 100])
        cbar.ax.set_yticklabels([r"$1$", r"$20$", r"$40$", r"$60$", r"$80$", r"$>100$"])
        plt.xscale("log")
        ax.set_title("Condition number")
        ax.set_ylabel(r"Working distance $L$ / [$\mu$m]")
        ax.set_xlabel(r"Sampling frequency $f$ / [$Hz$]")
        plt.tight_layout()
        plt.savefig("Plots/heatmap.png")
        # plt.savefig("Plots/heatmap.png", dpi=1000, format="png")
