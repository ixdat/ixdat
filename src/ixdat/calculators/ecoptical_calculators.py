# -*- coding: utf-8 -*-
"""
Created on Dec  14 15:07 2025

@author: Jiaze
"""
from ..measurement_base import Calculator
from ..data_series import ValueSeries
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import copy

class OpticalCVFitting(Calculator):
    """
    A Fitting for CV-type ECOptical data. The raw data is input as the ec-optical obj, 
    it will be ECCalibrated based on RE_vs_RHE, A_el, and/or R_Ohm values provided when the fiiting is initialized
    
    
    """
    calculator_type = "ECOptical fitting"
    def __init__(
        self,
        name=None,
        technique="EC-Optical",
        tstamp=None,
        measurement=None,
        RE_vs_RHE=None,
        A_el=None,
        R_Ohm=0,
        ref_spectra_t=0,
    ):
        """Initiate a Fitting

        Args:
            name (str): The name of the fitting
            technique (str): The technique of the fitting
            tstamp (float): The time at which the fitting took place or is valid
            measurement (EC-optical measurement):Mandatory. A measurement to calibrate by default
            RE_vs_RHE (float): The reference electrode potential on the RHE scale in [V]
            A_el (float): The electrode area in [cm^2]
            R_Ohm (float): The ohmic drop resistance in [Ohm]
        """
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.ecoptical = measurement
        self.RE_vs_RHE = RE_vs_RHE
        self.A_el = A_el
        self.R_Ohm = R_Ohm
        optical_t_mask = (self.ecoptical.spectrum_series.t>=self.ecoptical.t[0])&(self.ecoptical.spectrum_series.t<=self.ecoptical.t[-1])
        Uir_array = self.ecoptical.grab("WE(1).Potential (V)")[1]-self.ecoptical.grab('WE(1).Current (A)')[1]*self.R_Ohm # to be replaced by using ECCalibration
        Uir_f = interp1d(x = self.ecoptical.t,y = Uir_array)
        self.optical_U = Uir_f(self.ecoptical.spectrum_series.t[optical_t_mask])
        self.truncated_spectra_t = self.ecoptical.spectrum_series.t[optical_t_mask]
        self.truncated_spectra_data = self.ecoptical.spectra.data[optical_t_mask,:]
        self.truncated_spectra_data_OD = -np.log10(self.truncated_spectra_data/self.truncated_spectra_data[np.argmin(np.abs(self.ecoptical.spectrum_series.t[optical_t_mask]-ref_spectra_t)),:])
    def generate_assigned_spectra(
        self,
        t_span,
        V_ranges,
        assigned_spectra_smooth_param = None,
    ):
        t_mask = (self.truncated_spectra_t>=t_span[0])&(self.truncated_spectra_t<=t_span[1])
        self.assigned_spectra_dict = {}
        for key, value in V_ranges.items():
            if value is not None:
                spectra = self.truncated_spectra_data_OD[t_mask,:][self.optical_U[t_mask]<=value[1],:][-1,:] - self.truncated_spectra_data_OD[t_mask,:][self.optical_U[t_mask]>=value[0],:][0,:]
                print(self.optical_U[self.optical_U<=value[1]])
                if assigned_spectra_smooth_param is not None:
                    spectra = savgol_filter(
                            spectra,
                            assigned_spectra_smooth_param[0],
                            assigned_spectra_smooth_param[1],
                        )

                spectra = spectra / max(abs(spectra))
                self.assigned_spectra_dict[key] = spectra


    def growth_fit(
        self,
        t_span,
        potential_range_constraints = None,
    ):
        t_mask = (self.truncated_spectra_t>=t_span[0])&(self.truncated_spectra_t<=t_span[1])
        self.fit_results_dict = {}
        for key,value in potential_range_constraints.items():
            self.fit_results_dict[key] = [0]
        self.U_array = self.optical_U[t_mask]
        for i in range(self.optical_U[t_mask].shape[0]):
            self.fit_single_spectrum(
                spectrum=self.truncated_spectra_data_OD[t_mask][i],
                U_value=self.optical_U[t_mask][i],
                fit_results_dict=self.fit_results_dict,
                assigned_spectra_dict = self.assigned_spectra_dict,
                potential_range_constraints=potential_range_constraints
            )

    def objective_function(self,weights,components_array, target_vector):
        weighted_sum = np.dot(weights, components_array)
        diff = weighted_sum - target_vector
        return np.sum(np.square(diff))
    def fit_single_spectrum(self,spectrum,U_value,fit_results_dict,assigned_spectra_dict,potential_range_constraints):
        varying_components_list = []
        varying_components_key_list = []
        varying_weights_initial_list = []
        fixed_components_list = []
        fixed_weights_list = []
        for key,value in potential_range_constraints.items():
            if U_value>=value[0] and U_value<= value[1]:
                varying_components_list.append(assigned_spectra_dict[key])
                varying_weights_initial_list.append(fit_results_dict[key][-1])
                varying_components_key_list.append(key)
            else:
                fixed_components_list.append(assigned_spectra_dict[key])
                fixed_weights_list.append(fit_results_dict[key][-1])
        if fixed_components_list :
            target_vector = spectrum-np.dot(np.array(fixed_weights_list), np.array(fixed_components_list))
        else:
            target_vector = spectrum
        
        initial_weights = np.array(varying_weights_initial_list)
        COBYLA_constraints = []
        for i in range(len(initial_weights)):
            COBYLA_constraints.append(
                {"type": "ineq", "fun": lambda x, lb=-np.inf, i=i: x[i] - lb}
            )  # Weight must be greater than 0
            COBYLA_constraints.append(
                {"type": "ineq", "fun": lambda x, ub=np.inf, i=i: ub - x[i]}
            )  # Weight must be less than 10
        # Optimization options
        COBYLA_options = {
            "maxiter": 5000,  # Increase the maximum number of iterations
            "tol": 1e-6,  # Lower the tolerance for convergence
        }
        result = minimize(
            self.objective_function,
            initial_weights,
            args=(np.array(varying_components_list),target_vector,),
            method="COBYLA",
            constraints=COBYLA_constraints,
            options=COBYLA_options,
        )
        
        for key,value in fit_results_dict.items():
            if key in varying_components_key_list:
                value.append(result.x[varying_components_key_list.index(key)])
            else:
                value.append(value[-1])