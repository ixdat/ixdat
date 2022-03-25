import numpy as np
from scipy.interpolate import interp1d

from .ec import ECMeasurement
from ..spectra import Spectrum
from ..data_series import Field, ValueSeries
from ..spectra import SpectrumSeries
from ..exporters.sec_exporter import SECExporter
from ..plotters.sec_plotter import SECPlotter


class SpectroECMeasurement(ECMeasurement):

    default_plotter = SECPlotter
    default_exporter = SECExporter

    def __init__(self, *args, **kwargs):
        """Initialize an SEC measurement. All args and kwargs go to ECMeasurement."""
        ECMeasurement.__init__(self, *args, **kwargs)
        self._reference_spectrum = None
        self.tracked_wavelengths = []
        self.plot_waterfall = self.plotter.plot_waterfall
        self.plot_wavelengths = self.plotter.plot_wavelengths
        self.plot_wavelengths_vs_potential = self.plotter.plot_wavelengths_vs_potential
        self.technique = "S-EC"

    @property
    def reference_spectrum(self):
        """The spectrum which will by default be used to calculate dOD"""
        if not self._reference_spectrum or self._reference_spectrum == "reference":
            self._reference_spectrum = Spectrum.from_field(self["reference"])
        return self._reference_spectrum

    def set_reference_spectrum(
        self,
        spectrum=None,
        t_ref=None,
        V_ref=None,
    ):
        """Set the spectrum used as the reference when calculating dOD.

        Args:
            spectrum (Spectrum or str): If a Spectrum is given, it becomes the reference
                spectrum. The string "reference" can be given to make the reference
                spectrum become (via the reference_spectrum property) one that the
                measurement was loaded with (evt. for definition of wavelengths).
            t_ref (float): The time (with respect to self.tstamp) to use as the
                reference spectrum
            V_ref (float): The potential to use as the reference spectrum. This will
                only work if the potential is monotonically increasing.
        """
        if t_ref and not spectrum:
            spectrum = self.get_spectrum(t=t_ref)
        if V_ref and not spectrum:
            spectrum = self.get_spectrum(V=V_ref)
        if not spectrum:
            raise ValueError("must provide a spectrum, t_ref, or V_ref!")
        self._reference_spectrum = spectrum

    @property
    def spectra(self):
        """The Field that is the spectra of the SEC Measurement"""
        return self["spectra"]

    @property
    def spectrum_series(self):
        """The SpectrumSeries that is the spectra of the SEC Measurement"""
        return SpectrumSeries.from_field(
            self.spectra,
            tstamp=self.tstamp,
            name=self.name + " spectra",
        )

    @property
    def wavelength(self):
        """A DataSeries with the wavelengths for the SEC spectra"""
        return self.spectra.axes_series[1]

    @property
    def wl(self):
        """A numpy array with the wavelengths in [nm] for the SEC spectra"""
        return self.wavelength.data

    def calc_dOD(self, V_ref=None, t_ref=None, index_ref=None):
        """Calculate the optical density with respect to a reference

        Provide at most one of V_ref, t_ref, or index. If none are provided the default
        reference spectrum (self.reference_spectrum) will be used.

        Args:
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Return Field: the delta optical density spanning time and wavelength
        """
        counts = self.spectra.data
        if V_ref or t_ref:
            ref_spec = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
        else:
            ref_spec = self.reference_spectrum
        dOD = -np.log10(counts / ref_spec.y)
        dOD_series = Field(
            name=r"$\Delta$ O.D.",
            unit_name="",
            axes_series=self.spectra.axes_series,
            data=dOD,
        )
        return dOD_series

    def get_spectrum(self, V=None, t=None, index=None, name=None):
        """Return the Spectrum at a given potential V, time t, or index

        Exactly one of V, t, and index should be given. If V (t) is out of the range of
        self.U (self.t), then first or last spectrum will be returned.

        Args:
            V (float): The potential at which to get the spectrum. Measurement.U must
                be monotonically increasing for this to work.
            t (float): The time at which to get the spectrum
            index (int): The index of the spectrum
            name (str): Optional. name to give the new spectrum if interpolated

        Return Spectrum: The spectrum. The data is (spectrum.x, spectrum.y)
        """
        if V and V in self.U:  # woohoo, can skip interpolation!
            index = int(np.argmax(self.U == V))
        elif t and t in self.t:  # woohoo, can skip interpolation!
            index = int(np.argmax(self.t == t))
        if index:  # then we're done:
            return self.spectrum_series[index]
        # otherwise, we have to interpolate:
        counts = self.spectra.data
        end_spectra = (self.spectrum_series[0].y, self.spectrum_series[-1].y)
        if V:
            counts_interpolater = interp1d(
                self.U, counts, axis=0, fill_value=end_spectra, bounds_error=False
            )
            # FIXME: This requires that potential and spectra have same tseries!
            y = counts_interpolater(V)
            name = name or f"{self.spectra.name}_{V}V"
        elif t:
            t_spec = self.spectra.axes_series[0].t
            counts_interpolater = interp1d(
                t_spec, counts, axis=0, fill_value=end_spectra, bounds_error=False
            )
            y = counts_interpolater(t)
            name = name or f"{self.spectra.name}_{t}s"
        else:
            raise ValueError("Need t or V or index to select a spectrum!")

        field = Field(
            data=y,
            name=name,
            unit_name=self.spectra.unit_name,
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field, tstamp=self.tstamp)

    def get_dOD_spectrum(
        self,
        V=None,
        t=None,
        index=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
    ):
        """Return the delta optical density Spectrum given a point and reference point.

        Provide exactly one of V, t, and index, and at most one of V_ref, t_ref, and
        index_ref. For V and V_ref to work, the potential in the measurement must be
        monotonically increasing.

        Args:
            V (float): The potential at which to get the spectrum.
            t (float): The time at which to get the spectrum
            index (int): The index of the spectrum
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Return:
             Spectrum: The dOD spectrum. The data is (spectrum.x, spectrum.y)
        """
        if V_ref or t_ref or index_ref:
            spectrum_ref = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
        else:
            spectrum_ref = self.reference_spectrum
        spectrum = self.get_spectrum(V=V, t=t, index=index)
        field = Field(
            data=-np.log10(spectrum.y / spectrum_ref.y),
            name=r"$\Delta$ OD",
            unit_name="",
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field)

    def track_wavelength(self, wl, width=10, V_ref=None, t_ref=None, index_ref=None):
        """Return and cache a ValueSeries for the dOD for a specific wavelength.

        The caching adds wl_str to the SECMeasurement's data series, where
            wl_str = "w" + int(wl)
            This is dOD. The raw is also added as wl_str + "_raw".
        So, to get the raw counts for a specific wavelength, call this function and
            then use __getitem__, as in: sec_meas[wl_str + "_raw"]
        If V_ref, t_ref, or index_ref are provided, they specify what to reference dOD
            to. Otherwise, dOD is referenced to the SECMeasurement's reference_spectrum.

        Args:
            wl (float): The wavelength to track in [nm]
            width (float): The width around wl to average. For example, if wl=400 and
                width = 20, the spectra will be averaged between 390 and 410 nm to get
                the values. Defaults to 10. To interpolate at the exact wavelength
                rather than averaging, specify `width=0`.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Returns ValueSeries: The dOD value of the spectrum at wl.
        """
        if V_ref or t_ref or index_ref:
            spectrum_ref = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
        else:
            spectrum_ref = self.reference_spectrum
        x = self.wl
        if width:  # averaging
            wl_mask = np.logical_and(wl - width / 2 < x, x < wl + width / 2)
            counts_ref = np.mean(spectrum_ref.y[wl_mask])
            counts_wl = np.mean(self.spectra.data[:, wl_mask], axis=1)
        else:  # interpolation
            counts_ref = np.interp(wl, spectrum_ref.x, spectrum_ref.y)
            counts_wl = []
            for counts_i in self.spectra.data:
                c = np.interp(wl, x, counts_i)
                counts_wl.append(c)
            counts_wl = np.array(counts_wl)
        dOD_wl = -np.log10(counts_wl / counts_ref)
        raw_name = f"w{int(wl)} raw"
        dOD_name = f"w{int(wl)}"
        tseries = self.spectra.axes_series[0]
        raw_vseries = ValueSeries(
            name=raw_name, unit_name="counts", data=counts_wl, tseries=tseries
        )
        dOD_vseries = ValueSeries(
            name=dOD_name, unit_name="", data=dOD_wl, tseries=tseries
        )
        self.replace_series(raw_name, raw_vseries)
        # FIXME: better caching. See https://github.com/ixdat/ixdat/pull/11
        self.replace_series(dOD_name, dOD_vseries)
        # FIXME: better caching. See https://github.com/ixdat/ixdat/pull/11
        self.tracked_wavelengths.append(dOD_name)  # For the exporter.

        return dOD_vseries
