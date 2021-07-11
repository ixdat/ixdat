from .ec import ECMeasurement
from ..spectra import Spectrum
from ..data_series import Field
import numpy as np
from scipy.interpolate import interp1d
from ..spectra import SpectrumSeries


class SpectroECMeasurement(ECMeasurement):
    @property
    def reference_spectrum(self):
        return Spectrum.from_field(self["reference"])

    @property
    def spectra(self):
        """The Field that is the spectra of the SEC Measurement"""
        return self["spectra"]

    @property
    def spectrum_series(self):
        """The SpectrumSeries that is the spectra of the SEC Measurement"""
        return SpectrumSeries.from_field(
            self.spectra, tstamp=self.tstamp, name=self.name + " spectra"
        )

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
            # FIXME: The above line is in __init__ in other classes.

        return self._plotter

    def calc_dOD(self, V_ref=None):
        """Calculate the optical density with respect to a given reference potential"""
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

    def get_spectrum(self, V=None, t=None, index=None):
        if V and V in self.v:
            index = int(np.argmax(self.v == V))
        elif t and t in self.t:
            index = int(np.argmax(self.t == t))
        if index:
            return self.spectrum_series[index]
        counts = self.spectra.data
        if V:
            counts_interpolater = interp1d(self.v, counts, axis=0)
            # FIXME: This requires that potential and spectra have same tseries!
            y = counts_interpolater(V)
        elif t:
            t_spec = self.spectra.axes_series[0].t
            counts_interpolater = interp1d(t_spec, counts, axis=0)
            y = counts_interpolater(t)
        else:
            raise TypeError(f"Need t or v or index to select a spectrum!")

        field = Field(
            data=y,
            name=self.spectra.name,
            unit_name=self.spectra.unit_name,
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field, tstamp=self.tstamp)

    def get_dOD_spectrum(self, V=None, t=None, index=None, V_ref=None):
        if V_ref:
            spectrum_ref = self.get_spectrum(V=V_ref)
        else:
            spectrum_ref = self.reference_spectrum
        spectrum = self.get_spectrum(V=V, t=t, index=index)
        field = Field(
            data=np.log10(spectrum.y / spectrum_ref.y),
            name="$\Delta$ OD",
            unit_name="",
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field)
