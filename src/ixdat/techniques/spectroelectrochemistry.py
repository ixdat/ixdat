import numpy as np
from scipy.interpolate import interp1d

from .ec import ECMeasurement
from ..db import PlaceHolderObject, Relationship
from ..spectra import Spectrum, SpectroMeasurement, SpectrumSeries
from ..data_series import Field, ValueSeries
from ..exporters import SECExporter
from ..plotters import SECPlotter, ECOpticalPlotter
from ..tools import deprecate


class SpectroECMeasurement(SpectroMeasurement, ECMeasurement):
    """Electrochemistry with spectrometry."""

    default_exporter = SECExporter
    default_plotter = SECPlotter

    def __init__(self, **kwargs):
        """FIXME: Passing the right key-word arguments on is a mess"""
        ec_kwargs = {
            k: v for k, v in kwargs.items() if k in ECMeasurement.get_all_column_attrs()
        }
        spec_kwargs = {
            k: v
            for k, v in kwargs.items()
            if k in SpectroMeasurement.get_all_column_attrs()
        }
        # FIXME: I think the lines below could be avoided with a PlaceHolderObject that
        #  works together with MemoryBackend

        if "series_list" in kwargs:
            ec_kwargs.update(series_list=kwargs["series_list"])
            spec_kwargs.update(series_list=kwargs["series_list"])
        if "component_measurements" in kwargs:
            ec_kwargs.update(component_measurements=kwargs["component_measurements"])
            spec_kwargs.update(component_measurements=kwargs["component_measurements"])
        if "calibration_list" in kwargs:
            ec_kwargs.update(calibration_list=kwargs["calibration_list"])
            spec_kwargs.update(calibration_list=kwargs["calibration_list"])
        if "spectrum_series" in kwargs:
            spec_kwargs.update(spectrum_series=kwargs["spectrum_series"])
        SpectroMeasurement.__init__(self, **spec_kwargs)
        ECMeasurement.__init__(self, **ec_kwargs)


class ECXASMeasurement(SpectroECMeasurement):
    """Electrochemistry with X-ray Absorption Spectroscopy"""

    pass


class OpticalSpectrumSeries(SpectrumSeries):
    """Optical Spectrum Series"""

    extra_column_attrs = {"optical_spectrum_series": {"spectra_type"}}
    column_types = {"spectra_type": str}

    def __init__(self, *args, **kwargs):
        """Initiate an Optical Spectrum Series

        Args:
            spectra_type (str): The form of the raw data: "Intensity" (default),
                "Transmission" (in %), or "Absorption"
            Other args and kwargs are passed on to SpectrumSeries
        """
        self.spectra_type = kwargs.pop("spectra_type", "Intensity")
        super().__init__(*args, **kwargs)


class ECOpticalMeasurement(SpectroECMeasurement):
    """Electrochemistry with optical Spectroscopy

    This adds, to the SpectroElectrochemistry base class, methods for normalizing to a
    reference spectrum to get optical density, and for tracking intensity at specific
    wavelengths.
    """

    default_plotter = ECOpticalPlotter
    built_in_calculator_types = SpectroECMeasurement.built_in_calculator_types + [
        "dOD_calculator"
    ]

    relationships = {
        "reference_spectrum": Relationship(
            "spectrums", "ref_id", storage_table="ec_optical_measurements"
        )
    }

    def __init__(
        self,
        reference_spectrum=None,
        ref_id=None,
        spectra_type=None,
        **kwargs,
    ):
        """Initialize an SEC measurement. All args and kwargs go to ECMeasurement."""
        SpectroECMeasurement.__init__(self, **kwargs)

        if spectra_type is None:
            spectra_type = getattr(self.spectrum_series, "spectra_type", "Intensity")
        self.spectra_type = spectra_type

        if reference_spectrum:
            self._reference_spectrum = reference_spectrum
        elif ref_id:
            self._reference_spectrum = PlaceHolderObject(ref_id, cls=Spectrum)
        else:
            # A reference spectrum can be set later with set_reference_spectrum, or
            # left out for transmission or absorption data.
            self._reference_spectrum = None
        self.tracked_wavelengths = []
        # Other plots are available via self.plotter, e.g.
        # self.plotter.plot_dOD_difference_spectra()
        self.plot_waterfall = self.plotter.plot_waterfall
        self.technique = "EC-Optical"

    @deprecate(
        "0.4.0",
        "Use `plotter.plot_wavelengths()` instead.",
        remove_release="0.5.0",
    )
    def plot_wavelengths(self, **kwargs):
        """Plot dOD at specific wavelengths. See ECOpticalPlotter.plot_wavelengths"""
        return self.plotter.plot_wavelengths(**kwargs)

    @deprecate(
        "0.4.0",
        "Use `plotter.plot_wavelengths_vs_potential()` instead.",
        remove_release="0.5.0",
    )
    def plot_wavelengths_vs_potential(self, **kwargs):
        """Plot dOD at specific wavelengths vs potential. See
        ECOpticalPlotter.plot_wavelengths_vs_potential"""
        return self.plotter.plot_wavelengths_vs_potential(**kwargs)

    def cut(self, tspan, t_zero=None):
        """Select the portion of the data in a given tspan, keeping the reference.

        See :func:`~spectra.SpectroMeasurement.cut`
        """
        cut_measurement = super().cut(tspan, t_zero=t_zero)
        cut_measurement._reference_spectrum = self.reference_spectrum
        return cut_measurement

    @property
    def dOD_calculator(self):
        """The DeltaODCalculator providing "dOD" and ">diff" lookups"""
        return self.calculators["dOD_calculator"]

    @property
    def reference_spectrum(self):
        """The reference spectrum which will by default be used to calculate dOD"""
        if isinstance(self._reference_spectrum, PlaceHolderObject):
            self._reference_spectrum = self._reference_spectrum.get_object()
        return self._reference_spectrum

    @property
    def ref_id(self):
        """The id of the reference spectrum, or None if there is none"""
        if not self._reference_spectrum:
            return None
        return self.reference_spectrum.short_identity

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
        if t_ref is not None and spectrum is None:
            spectrum = self.get_spectrum(t=t_ref)
        if V_ref is not None and spectrum is None:
            spectrum = self.get_spectrum(V=V_ref)
        if spectrum is None:
            raise ValueError("must provide a spectrum, t_ref, or V_ref!")
        self._reference_spectrum = spectrum
        self.clear_cache()  # cached dOD is calculated with the old reference

    def set_reference_spectrum_through_t_range(
        self,
        spectrum=None,
        t_range=None,
        V_ref=None,
    ):
        """Set the reference spectrum, averaging spectra within a time range.

        Args:
            spectrum (Spectrum): If given, it becomes the reference spectrum.
            t_range (iterable): [t_start, t_end] in seconds relative to self.tstamp.
                All recorded spectra within the range, including its endpoints, are
                averaged. Used when spectrum is omitted. The range must contain at
                least one spectrum.
            V_ref (float): The potential at which to get the reference spectrum when
                spectrum and t_range are omitted. The potential must be monotonically
                increasing.
        """
        if t_range is not None and spectrum is None:
            spectra = self.spectra
            tseries = spectra.axes_series[0]
            t = tseries.t + tseries.tstamp - self.tstamp
            mask = (t >= t_range[0]) & (t <= t_range[1])
            if not np.any(mask):
                raise ValueError(f"No spectra in t_range={t_range}.")
            field = Field(
                name=f"{spectra.name} averaged over {t_range}",
                unit_name=spectra.unit_name,
                data=spectra.data[mask].mean(axis=0),
                axes_series=[self.wavelength],
            )
            spectrum = Spectrum.from_field(field, tstamp=self.tstamp + t[mask].mean())
        if V_ref is not None and spectrum is None:
            spectrum = self.get_spectrum(V=V_ref)
        if spectrum is None:
            raise ValueError("must provide a spectrum, t_range, or V_ref!")
        self.set_reference_spectrum(spectrum=spectrum)

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

        See :func:`~calculators.sec_calculators.DeltaODCalculator.calc_dOD`. With no
        reference arguments, this is the same as `self["dOD"]`.

        Args:
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Return Field: the delta optical density spanning time and wavelength
        """
        if V_ref is None and t_ref is None and index_ref is None:
            return self["dOD"]
        return self.dOD_calculator.calc_dOD(
            self, V_ref=V_ref, t_ref=t_ref, index_ref=index_ref
        )

    def get_spectrum(self, V=None, t=None, index=None, name=None, interpolate=True):
        """Return the Spectrum at a given potential V, time t, or index

        Exactly one of V, t, and index should be given. If V (t) is out of the range of
        self.U (self.t), then first or last spectrum will be returned.

        Args:
            V (float): The potential at which to get the spectrum. Measurement.U must
                be monotonically increasing for this to work.
            t (float): The time at which to get the spectrum
            index (int): The index of the spectrum
            name (str): Optional. name to give the new spectrum if interpolated
            interpolate (bool): Optional. Set to false to grab closest spectrum rather
                than interpolating.

        Return Spectrum: The spectrum. The data is (spectrum.x, spectrum.y)
        """
        if V is not None and V in self.U:  # woohoo, can skip interpolation!
            index = int(np.argmax(self.U == V))
        elif t is not None and t in self.t:  # woohoo, can skip interpolation!
            index = int(np.argmax(self.t == t))
        if index is not None:  # then we're done:
            return self.spectrum_series[index]
        # otherwise, we have to interpolate:
        counts = self.spectra.data
        end_spectra = (self.spectrum_series[0].y, self.spectrum_series[-1].y)
        if V is not None:
            if interpolate:
                V_spec = self.grab_for_t(
                    self.U_name,
                    t=self.spectra.axes_series[0].t,
                )
                counts_interpolater = interp1d(
                    V_spec, counts, axis=0, fill_value=end_spectra, bounds_error=False
                )
                y = counts_interpolater(V)
            else:
                U_diff = np.abs(self.U - V)
                index = np.argmin(U_diff)
                y = counts[index]
            name = name or f"{self.spectra.name}_{V}V"
        elif t is not None:
            t_spec = self.spectra.axes_series[0].t
            if interpolate:
                counts_interpolater = interp1d(
                    t_spec, counts, axis=0, fill_value=end_spectra, bounds_error=False
                )
                y = counts_interpolater(t)

            else:
                t_diff = np.abs(t_spec - t)
                index = np.argmin(t_diff)
                y = counts[index]

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
        interpolate=True,
    ):
        """Return the delta optical density Spectrum given a point and reference point.

        Provide exactly one of V, t, and index, and at most one of V_ref, t_ref, and
        index_ref. For V and V_ref to work, the potential in the measurement must be
        monotonically increasing. If none of V_ref, t_ref, or index_ref are provided,
        the default reference spectrum (self.reference_spectrum) is used. Without any
        reference, intensity data raises a ValueError, while absorption and
        transmission data are taken as referenced during data collection.

        Args:
            V (float): The potential at which to get the spectrum.
            t (float): The time at which to get the spectrum
            index (int): The index of the spectrum
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
            interpolate (bool): Optional. Set to false to grab closest spectrum rather
                than interpolating.
        Return:
             Spectrum: The dOD spectrum. The data is (spectrum.x, spectrum.y)
        """

        from ..calculators.sec_calculators import calc_dOD_values

        spectrum = self.get_spectrum(V=V, t=t, index=index, interpolate=interpolate)
        if V_ref is not None or t_ref is not None or index_ref is not None:
            spectrum_ref = self.get_spectrum(
                V=V_ref, t=t_ref, index=index_ref, interpolate=interpolate
            )
        else:
            spectrum_ref = self.reference_spectrum
        y_ref = spectrum_ref.y if spectrum_ref is not None else None
        dOD = calc_dOD_values(spectrum.y, y_ref, self.spectra_type)

        field = Field(
            data=dOD,
            name=r"$\Delta$ OD",
            unit_name="",
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field, tstamp=self.tstamp)

    def track_wavelength(self, wl, width=10, V_ref=None, t_ref=None, index_ref=None):
        """Return and cache a ValueSeries for the dOD for a specific wavelength.

        The caching adds wl_str to the SECMeasurement's data series, where
            wl_str = "w" + int(wl)
            This is dOD. The raw is also added as wl_str + " raw".
        So, to get the raw counts for a specific wavelength, call this function and
            then use __getitem__, as in: sec_meas[wl_str + " raw"]
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
        dOD = self.calc_dOD(V_ref=V_ref, t_ref=t_ref, index_ref=index_ref).data
        counts = self.spectra.data
        x = self.wl
        if width:  # averaging
            wl_mask = np.logical_and(wl - width / 2 < x, x < wl + width / 2)
            counts_wl = np.mean(counts[:, wl_mask], axis=1)
            dOD_wl = np.ma.filled(np.mean(dOD[:, wl_mask], axis=1), np.nan)
        else:  # interpolation
            counts_wl = np.array([np.interp(wl, x, c) for c in counts])
            dOD_wl = np.array([np.interp(wl, x, d) for d in np.ma.filled(dOD, np.nan)])
        raw_name = f"w{int(wl)} raw"
        dOD_name = f"w{int(wl)}"
        tseries = self.spectra.axes_series[0]
        raw_vseries = ValueSeries(
            name=raw_name,
            unit_name=self.spectra.unit_name,
            data=counts_wl,
            tseries=tseries,
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

    def denoise_spectra(
        self,
        spectra_field=None,
        denoise_method=None,
        PCA_explained_variance=None,
        sg_window=None,
        sg_poly_order=None,
        normalise=None,
        wlmin=None,
        wlmax=None,
    ):
        """Return a field of denoised spectra, `self.spectra` by default.

        See :func:`~calculators.sec_calculators.Denoiser.denoise_spectra`. Arguments
        left as None take the `Denoiser` defaults. To look up denoised spectra
        instead, add a `Denoiser` with `self.add_calculator` and use e.g.
        `self["dOD>smooth"]`.

        Args:
            spectra_field (Field): The spectra to denoise. Defaults to self.spectra.
            denoise_method (str): "PCA" (default) or "Savitzky-Golay"
            PCA_explained_variance (float): The fraction of the variance that the
                kept principal components must explain. 0.99 by default.
            sg_window (int): The window length of the Savitzky-Golay filter. 10 by
                default.
            sg_poly_order (int): The polynomial order of the Savitzky-Golay filter.
                3 by default.
            normalise (bool): Whether to scale each spectrum to [0, 1]. False by
                default.
            wlmin (float): The minimum wavelength in [nm] used for normalisation
            wlmax (float): The maximum wavelength in [nm] used for normalisation
        Return Field: The denoised spectra
        """
        from ..calculators.sec_calculators import Denoiser

        if spectra_field is None:
            spectra_field = self.spectra
        xspan = None
        if wlmin is not None or wlmax is not None:
            wl = spectra_field.axes_series[1].data
            xspan = [
                wlmin if wlmin is not None else np.min(wl),
                wlmax if wlmax is not None else np.max(wl),
            ]
        return Denoiser().denoise_spectra(
            spectra_field,
            denoise_method=denoise_method,
            PCA_explained_variance=PCA_explained_variance,
            sg_window=sg_window,
            sg_poly_order=sg_poly_order,
            normalise=normalise,
            xspan=xspan,
        )
