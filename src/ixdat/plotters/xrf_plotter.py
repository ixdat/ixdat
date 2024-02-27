# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 17:14:59 2024

@author: cl2120
"""

from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter


class TRXRFPlotter(MPLPlotter):
    
    
    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement
    
    def plot_measurement(
            self, measurement=None, ax=None, tspan=None, y_name="FF_over_I0", **kwargs
        ):
        measurement = measurement or self.measurement 
        
        t, y = measurement.grab(y_name, tspan=tspan)
        
        time_name = measurement["FF_over_I0"].tseries.name 
        
        if not ax:
            ax = self.new_ax(xlabel=time_name, ylabel=y_name)
        
        ax.plot(t, y, **kwargs)
        
        return ax




class ECTRXRFPlotter(MPLPlotter):
    
    
    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=self.measurement)
    
    def plot_measurement(
            self, 
            measurement=None, 
            axes=None, 
            tspan=None,
            y_name="FF_over_I0", 
            U_name=None,
            J_name=None,
            U_color=None,
            J_color=None,
            **kwargs
        ):
        
        if not axes:
            axes = self.new_two_panel_axes(n_bottom=2, n_top=1, emphasis=None)
            
        measurement = measurement or self.measurement 
        t, y = measurement.grab(y_name, tspan=tspan)
        time_name = measurement["FF_over_I0"].tseries.name 
        
        axes[0].set_ylabel(y_name)
        axes[0].set_xlabel(time_name)
        axes[0].plot(t, y, **kwargs)
        
        self.ec_plotter.plot_measurement(
            measurement=measurement,
            axes=[axes[1], axes[3]],
            tspan=tspan,
            U_name=U_name,
            J_name=J_name,
            U_color=U_color,
            J_color=J_color,
            **kwargs,
        )
        axes[1].set_xlim(axes[0].get_xlim())
        
        return axes
 
    plot_vs_potential = ECPlotter.plot_vs_potential
        