from .ec import ECMeasurement
from ..spectra import Spectrum
from ..data_series import Field
import numpy as np
from scipy.interpolate import interp1d


class SpectroECMeasurement(ECMeasurement):
    @property
    def reference_spectrum(self):
        return Spectrum.from_field(self["reference"])

    @property
    def spectra(self):
        return self["spectra"]

    @property
    def wavelength(self):
        return self.spectra.axes_series[1]

    @property
    def wl(self):
        return self.wavelength.data

    @property
    def plotter(self):
        """The default plotter for ECMeasurement is ECPlotter"""
        if not self._plotter:
            from ..plotters.sec_plotter import SECPlotter

            self._plotter = SECPlotter(measurement=self)

            self.plot_waterfall = self._plotter.plot_waterfall

        return self._plotter

    def calc_dOD(self, V_ref=None):
        counts = self.spectra.data
        if V_ref:
            counts_interpolater = interp1d(self.v, counts, axis=0)
            ref_counts = counts_interpolater(V_ref)
        else:
            ref_counts = self.reference_spectrum.y
        dOD = np.log10(counts / ref_counts)
        dOD_series = Field(
            name="$\Delta$ O.D.",
            unit_name="",
            axes_series=self.spectra.axes_series,
            data=dOD,
        )
        return dOD_series
