# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 17:14:59 2024

@author: cl2120
"""

from .base_mpl_plotter import MPLPlotter


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
        