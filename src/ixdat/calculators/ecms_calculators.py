# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 20:48:06 2024

@author: Søren
"""
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize  # noqa
from scipy.interpolate import interp1d  # noqa
from scipy import signal  # noqa
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth  # noqa
from numpy.fft import fft, ifft, ifftshift, fftfreq  # noqa
from ..plotters.ms_plotter import STANDARD_COLORS
from ..constants import FARADAY_CONSTANT, R, STANDARD_TEMPERATURE, STANDARD_PRESSURE
from .ms_calculators import MSCalibration, MSCalResult
from ..measurements import Calculator
from ..data_series import ValueSeries, TimeSeries
from ..exceptions import TechniqueError
from ..config import plugins
from ..plotters.plotting_tools import calc_linear_background


class ECMSCalibration(Calculator):
    """Class for calibrations done using ECMSMeasurements

    Its classmethods return MSCalibration objects.
    """

    @classmethod
    def ecms_calibration(cls, measurement, mol, mass, n_el, tspan, tspan_bg=None):
        """Calibrate for mol and mass based on one period of steady electrolysis

        Args:
            measurement (ECMSMeasurement): measurement with calibration data
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            n_el (str): Number of electrons passed per molecule produced (remember the
                sign! e.g. +4 for O2 by OER and -2 for H2 by HER)
            tspan (tspan): The timespan of steady electrolysis
            tspan_bg (tspan): The time to use as a background

        Return MSCalResult: The result of the ecms_calibration
        """
        if plugins.use_siq:
            warnings.warn(
                "spectro_inlets_quantification is active but you are using the native "
                "ixdat version of `ECMSMeasurement.ecms_calibration`"
            )
        Y = measurement.integrate_signal(mass, tspan=tspan, tspan_bg=tspan_bg)
        Q = measurement.integrate("raw_current", tspan=tspan) * 1e-3
        n = Q / (n_el * FARADAY_CONSTANT)
        F = Y / n
        cal = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type="ecms_calibration",
            F=F,
        )
        return MSCalibration(ms_cal_results=[cal], measurement=measurement)

    @classmethod
    def ecms_calibration_curve(
        cls,
        measurement,
        mol,
        mass,
        n_el,
        tspan_list=None,
        selector_name=None,
        selector_list=None,
        t_steady_pulse=None,
        tspan_bg=None,
        force_through_zero=False,
        ax="new",
        axes_measurement=None,
        axes_measurement_J_name="raw_current",
        return_ax=False,
    ):
        """Fit mol's sensitivity at mass based on steady periods of EC production.

        Args:
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            n_el (str): Number of electrons passed per molecule produced (remember the
                sign! e.g. +4 for O2 by OER and -2 for H2 by HER)
            tspan_list (list of tspan): The timespans of steady electrolysis
            selector_name (str): Name of selector which identifies the periods
                of steady electrolysis for automatic selection of timespans of steady
                electrolysis. E.g. "selector" or "Ns" for biologic EC data
            selector_list (list): List of values for selector_name for automatic
                selection of timespans of steady electrolysis
            t_steady_pulse (float): Length of steady electrolysis for each segment
                given by selector_list. Defaults to None = entire length of segment
            tspan_bg (tspan): The time to use as a background
            force_through_zero (boolean): Whether to force the calibration curve through
                zero. This can be done when confident in the background subtraction.
            ax (Axis): The axis on which to plot the ms_calibration curve result.
                Defaults to a new axis.
            axes_measurement (list of Axes): The EC-MS plot axes to highlight the
                ms_calibration on. Defaults to None. These axes are not returned.
            axes_measurement_J_name (str): The J_name used in the axis passed
                to axes_measurement. Must be passed manually as the axis does not "know"
                its J_name. Defaults to "raw_current". IMPORTANT: the method still uses
                "raw_current" to calculate the sensitivity factor, this J_name is only
                used for plotting.
            return_ax (bool): Whether to return the axis on which the calibration curve
                is plotted together with the MSCalResult. Defaults to False.

        Return MSCalResult(, Axis): The result of the ms_calibration (and calibration
            curve axis if requested) based on integration of selected time periods.
        """

        axis_ms = axes_measurement[0] if axes_measurement else None
        axis_current = axes_measurement[3] if axes_measurement else None
        Y_list = []
        n_list = []
        if not tspan_list:
            tspan_list = measurement._get_tspan_list(
                selector_list, selector_name, t_steady_pulse
            )
        for tspan in tspan_list:
            Y = measurement.integrate_signal(
                mass, tspan=tspan, tspan_bg=tspan_bg, ax=axis_ms
            )
            # FIXME: plotting current by giving integrate() an axis doesn't work great.
            if (
                axes_measurement
            ):  # FIXME: need to run twice, once to plot, once to calculate Q
                measurement.integrate(
                    axes_measurement_J_name, tspan=tspan, ax=axis_current
                )
            Q = measurement.integrate("raw_current", tspan=tspan)
            Q *= 1e-3  # mC --> [C]
            n = Q / (n_el * FARADAY_CONSTANT)
            Y_list.append(Y)
            n_list.append(n)
        n_vec = np.array(n_list)
        Y_vec = np.array(Y_list)
        n_fit = np.array([0, max(n_vec)])
        if force_through_zero:

            def rms_error(F_guess):
                return np.mean((Y_vec - F_guess * n_vec) ** 2)

            F_guess_0 = np.sum(Y_vec) / np.sum(n_vec)
            res = minimize(rms_error, F_guess_0)
            F = res.x[0]
            Y_fit = n_fit * F
        else:
            pfit = np.polyfit(n_vec, Y_vec, deg=1)
            F = pfit[0]
            Y_fit = n_fit * pfit[0] + pfit[1]

        if ax:
            color = STANDARD_COLORS[mass]
            if ax == "new":
                ax = measurement.plotter.new_ax()
                ax.set_xlabel("amount produced / [nmol]")
                ax.set_ylabel("integrated signal / [nC]")
            ax.plot(n_vec * 1e9, Y_vec * 1e9, "o", color=color)
            ax.plot(n_fit * 1e9, Y_fit * 1e9, "--", color=color)

        cal = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type="ecms_calibration_curve",
            F=F,
        )
        return MSCalibration(ms_cal_results=[cal], measurement=measurement)

    from_siq = MSCalibration.from_siq


class ECMSImpulseResponse(Calculator):
    """
    Class for signal impulse response for ECMS data for modelling diffusion.
    see Krempl et al. 2021
    https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110

    This class currently handles only one molecule at a time. If deconvolution of
    multiple molecules is required this needs to be handled by separate objects.

    # TODO: Make class inherit from Calculator
    """

    calculator_type = "ecms_deconvolution"

    def __init__(
        self,
        mol,
        t_kernel,
        kernel,
        name=None,
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
        name = name or f"ECMSImpulseResponse(mol={mol})"
        super().__init__(name=name, technique="EC-MS", measurement=measurement)
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

    def grab_deconvoluted_signal(
        self, mol, measurement=None, tspan=None, tspan_bg=None, snr=10
    ):
        """Return the mass transport deconvoluted MS signal for a given MS signal using
        the algorithm developed by Krempl et al.
        https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110

        Note, this actually doesnt need the EC data - it is calculated
        from calibrated MS data. It is only meaningful for ECMS measurements
        at the moment, justifying placement here.

        Args:
            mol (str): Name of molecule for which deconvolution is to be carried out.
            impulse_response (ECMSImpulseResponse): The impulse response must contain all
                the attribute necessary to calculate a new ECMSImpulseResponse with
                adjusted time and sample frequency. Therefore currently only CALCULATED
                ECMSImpulseResponse objects are allowed.
            tspan (list): Timespan for which the deconvolued signal is returned. Needs
                          to contain time zero (for alignment with the model).
            tspan_bg (list): Timespan that corresponds to the background signal.
            snr (int): signal-to-noise ratio used for Wiener deconvolution.
        Returns tuple of (time [s], deconvoluted MS signal [mol/s])
        """
        # grab the calibrated data
        measurement = measurement or self.measurement
        t_sig, v_sig = measurement.grab_flux(mol, tspan=tspan, tspan_bg=tspan_bg)

        # first check that the impulse reponse passed is from parameters, because
        # otherwise this might not work for now. Should make this work in the future!
        # There is no reason why a measured impulse response shouldnt be just as useful
        # for deconvolution, but only possible if dt and duration are the same as in
        # the measurement that is to be deconvoluted.
        if self.ir_type != "calculated":
            raise TechniqueError(
                "You need to pass an ECMSImpulseResponse object calculated"
                "from parameters to calculate deconvoluted currents. "
                "(for now)"
            )  # TODO: check instead whether the object contains all the necessary
            # attributes

        # Check if the one passed already has the correct dt and duration
        dt = t_sig[1] - t_sig[0]
        duration = t_sig[-1] - t_sig[0]

        if dt == self.dt and duration == self.duration:
            # no need to recalculate if these parameters fit
            signal_response = self
        else:
            # re-calculate the impulse response
            # FIXME: There must be a better way!
            signal_response = ECMSImpulseResponse.from_parameters(
                mol=mol,
                measurement=self,
                working_distance=self.working_distance,
                A_el=self.A_el,
                D=self.D,
                H_v_cc=self.H_v_cc,
                n_dot=self.n_dot,
                T=self.T,
                p=self.p,
                carrier_gas=self.carrier_gas,
                gas_volume=self.gas_volume,
                dt=dt,  # as defined from measurement above
                duration=duration,  # as defined from measurement above
            )
        # make sure signal_response and v_sig are same length for the steps below by
        # adding zeros to the end of whichever array is shorter
        if len(v_sig) >= len(signal_response.kernel):
            kernel = np.hstack(
                (
                    signal_response.kernel,
                    np.zeros(len(v_sig) - len(signal_response.kernel)),
                )
            )
            v_sig_corr = v_sig
        else:
            kernel = signal_response.kernel
            v_sig_corr = np.hstack(
                (v_sig, np.zeros(len(signal_response.kernel) - len(v_sig)))
            )
        # calculate the convolution function from the calculated kernel
        H = fft(kernel)  # (see Krempl et al 2019, SI, page S4 bottom)
        # TODO: cache this somehow
        decon_signal = np.real(
            ifft(fft(v_sig_corr) * np.conj(H) / (H * np.conj(H) + (1 / snr) ** 2))
        )
        # see Krempl et al 2019, SI, eq. 26 and paragraph below) -
        # SNR in equ = (1 / snr) ** 2 here?
        decon_signal = decon_signal * sum(kernel)  # what does this do????
        # Now finally make sure t_sig and the calculated deconvoluted signal are the
        # same length (for plotting etc later)
        if len(t_sig) < len(decon_signal):
            delta = len(decon_signal) - len(t_sig)
            decon_signal = decon_signal[:-delta]
        return t_sig, decon_signal

    def deconvolute_for_tspans(
        self,
        tspan_list,
        mol,
        measurement=None,
        t_zero_list=None,
        plot=True,
        name=None,
        t_bg=[-10, -1],
        snr=7,
        return_t_v_list=False,
        export_data=False,
    ):
        """
        Loops though list of tspans and associated t_zero list to deconvolute using
        the given impulse response (from model, but could also be from data) for
        the molecule the impulse resp
        Args:
            tspan_list (list): list of tspans to devonvolute data over. if no t_zero
                is given needs to include zero.
            t_zero_list (list): Optional. zero point where the deconvolution start
            measurement (ECMSMeasurement): A measurement object. Has to be calibrated
                so that the flux of mol is available.
            mol (str): molecule
            plot (bool): will return a plot of each tspan and save using name.The plots
                will only be meaningful if the current is calibrated.
            name (str): str to use for title in figure and saving
            t_bg (tspan): tspan to be used as background IN RELATION TO t_zero.
                Will be the same for each tspan. # TODO: add list option
            snr (int): signal-to-noise ratio used for Wiener deconvolution
                (see grab_deconvoluted_signal()). Defaults to 7 (works with test data).
            return_t_v_list (bool): Whether to return list of (time, deconvoluted
                MS signal), Defaults to False
            export_data (bool): save raw and deconvoluted data as csv using name

        Return t_v_list (list): list of tuple of (time [s], deconvoluted MS signal
                                [mol/s]) as returned from grab_deconvoluted_signal()
        """
        measurement = measurement or self.measurement
        t_v_list = []
        if t_zero_list is None:
            t_zero_list = [None for tspan in tspan_list]
        for tspan, t_zero in zip(tspan_list, t_zero_list):
            print(
                "Now working on tspan {}, starting sequence at t_zero {}s".format(
                    tspan, t_zero
                )
            )
            data_snippet = measurement.cut(tspan=tspan, t_zero=t_zero)
            t_decon_signal, v_decon_signal = self.grab_deconvoluted_signal(
                mol=mol,
                measurement=data_snippet,
                tspan=None,
                tspan_bg=t_bg,
                snr=snr,
            )
            t_v_list.append((t_decon_signal, v_decon_signal))

            if plot:
                axes = data_snippet.plot(
                    mol_list=[mol], logplot=False, tspan_bg=t_bg, alpha=0.5
                )
                # TODO: find a way to have different hue on the lines of mol and EC data
                axes[0].plot(t_decon_signal, v_decon_signal, color=STANDARD_COLORS[mol])
                axes[0].get_figure().savefig(name + "tstart_" + str(t_zero) + ".png")

            if export_data:
                # TODO use the ixdat csv exporter (after converted to Calculator)
                raw_I_t, raw_I = data_snippet.grab("raw_current")
                i_t, i = data_snippet.grab("J / [mA cm$^{-2}$]")
                raw_U_t, raw_U = data_snippet.grab("raw_potential")
                raw_sig_t, raw_sig = data_snippet.grab_flux(mol)
                export_dict = {
                    "time raw current / s": raw_I_t,
                    "raw current / mA": raw_I,
                    "time J / s": i_t,
                    "J / [mA cm$^{-2}$]": i,
                    "time potential / s": raw_U_t,
                    "raw_potential / V": raw_U,
                    "time measured " + mol + " flux/ s": raw_sig_t,
                    "measured " + mol + " flux / mol/s": raw_sig,
                    " time deconvoluted " + mol + " flux / s": t_decon_signal,
                    "deconvoluted " + mol + " flux / mol/s": v_decon_signal,
                }
                export_df = pd.DataFrame(
                    {key: pd.Series(value) for key, value in export_dict.items()}
                )
                export_df.to_csv(name + "tstart_" + str(t_zero) + ".csv", index=False)
        if return_t_v_list:
            return t_v_list

    # TODO: Would be very useful to have a method to calculate the condition number for
    #   a certain analyte. see Krempl et al. 2021
    #   https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
    #   for definition of condition number

    @property
    def available_series_names(self):
        return set([f"n_dot_{self.mol}-deconvoluted"])

    def calculate_series(self, key, measurement=None):
        mol = key.removesuffix("-deconvoluted").removeprefix("n_dot_")
        if mol != self.mol:
            print(
                "Warning: tried (and failed) to look up "
                f"{key} in a deconvolution for {self.mol}"
            )
            return

        t, n_dot_decon = self.grab_deconvoluted_signal(measurement=measurement, mol=mol)
        return ValueSeries(
            name=key,
            data=n_dot_decon,
            unit_name="mol/s",
            tseries=TimeSeries(
                name="time for {key}", data=t, unit_name="s", tstamp=measurement.tstamp
            ),
        )


# TODO: other potentially useful methods:
# https://github.com/ixdat/ixdat/blob/f577a434a966e486cf4cb66253677f7839fe117a/src/
