import numpy as np
from scipy.interpolate import interp1d

from .ec import ECMeasurement
from ..db import PlaceHolderObject
from ..spectra import Spectrum, SpectroMeasurement, SpectrumSeries
from ..data_series import Field, ValueSeries, TimeSeries
from ..exporters import SECExporter
from ..plotters import SECPlotter, ECOpticalPlotter
import warnings


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

    def __init__(self, *args, **kwargs):
       """Initiate an Optical Spectrum Series, as an addition to the spectrum series 
       spectra_type: The form of the raw data (Intensity, Transmission, or Absorbance. Defaults to Intensity)
       """
       self.spectra_type = kwargs.pop("spectra_type", "Intensity")
       super().__init__(*args, **kwargs)
        

    def calc_average_spectra(
        self,
        average_number=10,
        field=None,
    ):
        """Return a ValueSeries of averaged spectra for either an OpticalSpectrumSeries
        or specific field containing spectral data.

        Args:
            average_number (int) : The number of spectra to average. 10 by default.
            field (field) : The spectra data to average, if not included in the
                measurement
        Returns ValueSeries: An optical spectrum series containing averaged spectra
        """
        spectrum = self
        averaged_spectra = []
        averaged_time = []
        if field is not None:
            spectra_field = field
            spectra = field.data
            time = field.axes_series[0].data
        else:
            spectra_field = spectrum.field
            spectra = spectra_field.data
            time = spectra_field.axes_series[0].data

        for index in range(0, len(spectra) - average_number, average_number):
            window_end = index + average_number
            averaged_spectra.append(spectra[index:window_end].mean(axis=0))
            averaged_time.append(time[index:window_end].mean(axis=0))

        averaged_spectra = np.array(averaged_spectra)
        averaged_time = np.array(averaged_time)

        tseries = TimeSeries(
            name="time",
            unit_name=spectra_field.axes_series[0].unit_name,
            data=averaged_time,
            tstamp=spectra_field.axes_series[0].tstamp,
        )

        averaged_spectra_field = Field(
            data=averaged_spectra,
            name=f"Averaged {spectra_field.name}",
            unit_name="",
            axes_series=[tseries, spectra_field.axes_series[1]],
        )

        cls = OpticalSpectrumSeries
        averaged_uvvis_series = cls(
            name=averaged_spectra_field.name,
            reader=Spectrum,
            technique="Optical",
            tstamp=tseries.tstamp,
            field=averaged_spectra_field,
            continuous=True,
            spectra_type=self.spectra_type,
        )

        return averaged_uvvis_series


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

    extra_linkers = SpectroECMeasurement.extra_linkers.copy()
    extra_linkers.update({"ec_optical_measurements": ("spectra", "ref_id")})

    def __init__(
        self,
        reference_spectrum=None,
        ref_id=None,
        spectra_type=None,
        **kwargs,
    ):  # Added in fix for reference spectrum issue
        """Initialize an SEC measurement. All args and kwargs go to ECMeasurement."""
        SpectroECMeasurement.__init__(self, **kwargs)

        if spectra_type is None:
            spectra_type = self.spectrum_series.spectra_type
        self.spectra_type = spectra_type
        self.ref_id = ref_id

        if reference_spectrum:
            self._reference_spectrum = reference_spectrum
        elif ref_id:
            self._reference_spectrum = PlaceHolderObject(ref_id, cls=Spectrum)
        else:
            self._reference_spectrum = None #Sometimes the user may want to not provide a reference spectrum, e.g for transmission or absorption data
        self.tracked_wavelengths = []
        self.plot_waterfall = self.plotter.plot_waterfall
        self.plot_wavelengths = self.plotter.plot_wavelengths
        self.plot_dOD_difference_spectra = self.plotter.plot_dOD_difference_spectra
        self.plot_converged_spectra = self.plotter.plot_converged_spectra
        self.plot_fit_and_residuals = self.plotter.plot_fit_and_residuals
        self.plot_fit_reconstruction = self.plotter.plot_fit_reconstruction
        self.plot_wavelengths_vs_cv = self.plotter.plot_wavelengths_vs_cv
        self.plot_wavelengths_vs_potential = self.plotter.plot_wavelengths_vs_potential
        self.technique = "EC-Optical"

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
        if t_range and not spectrum:
            spectrum = self.get_spectrum(t=t_range)
        if V_ref and not spectrum:
            spectrum = self.get_spectrum(V=V_ref)
        if spectrum is None:
            raise ValueError("must provide a spectrum, t_ref, or V_ref!")
        self._reference_spectrum = spectrum

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
        interpolate=True
    ):
        """Return the delta optical density Spectrum given a point and reference point.

        Provide exactly one of V, t, and index, and at most one of V_ref, t_ref, and
        index_ref. For V and V_ref to work, the potential in the measurement must be
        monotonically increasing. If none of V_ref, t_ref, or index_ref are provided the default
        reference spectrum (self.reference_spectrum) will be used. If the default reference
        spectrum is missing, the function wil raise a Value error for intensity data, or 
        assume an unknown pre-set reference for absorption or transmission data.

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

        spectrum = self.get_spectrum(V=V, t=t, index=index, interpolate=interpolate)
        spectrum_ref=None
        #---- Check form of data and reference spectrum. Note that:
            #(1) If no reference spectrum is given and data is in intensity 
            #    form, the dOD cannot be calculated
            #(2) Absorption and Transmission data are already referenced to some 
            #    unknown reference point, presumably set during data collection. 
            #    If no reference spectrum is given, it is assumed that the user 
            #    wants to maintain this unknown reference point.
            #(3) If a reference spectrum is given, the dOD is given relative to 
            #    this reference point.
        if V_ref or t_ref or index_ref:
            spectrum_ref = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref, interpolate=interpolate)
        elif self.reference_spectrum:
            spectrum_ref = self.reference_spectrum
        if not spectrum_ref:
            if self.spectra_type=="Intensity":
                raise ValueError("For intensity data, a reference spectrum must be provided to find the optical density")
            if self.spectra_type == "Absorption":
                warnings.warn("Absorption set as optical data type and no reference spectrum is given. Optical density referred to an unknown reference.")
                dOD = spectrum.y
            if self.spectra_type == "Transmission":
                warnings.warn("Transmission set as optical data type and no reference spectrum is given. Optical density referred to an unknown reference.")
                spec_decimal = spectrum.y / 100
                spec_decimal[spec_decimal == 0] = np.nan
                dOD = -np.log10(spec_decimal)
        else:
            if self.spectra_type=="Intensity" or self.spectra_type=="Transmission":
                ratio = np.divide(
                    spectrum.y,
                    spectrum_ref.y,
                    out=np.full_like(spectrum.y, np.nan),
                    where=(abs(spectrum.y) > 0) & (abs(spectrum_ref.y) > 0),
                )
                dOD = -np.log10(ratio)
            if self.spectra_type=="Absorption":
                dOD = np.divide(
                    spectrum.y,
                    spectrum_ref.y,
                    out = np.full_like(spectrum.y, np.nan),
                    where = abs(spectrum_ref.y)>0
                    )

        dOD = np.ma.masked_invalid(dOD)

        field = Field(
            data=dOD,
            name=r"$\Delta$ OD",
            unit_name="",
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field, tstamp=self.tstamp)
    
    def get_dOD_cycle_noise(
        self,
        c_name="cycle",
        cycle_1=1,
        cycle_2=2,
        V_ref=None,
        t_ref=None,
        index_ref=None,
    ):
        """Takes an ECOpticalMeasurement and returns the irreversible changes in
        absorbance between the two cycles.
        Used to check for system stability.
        Exactly one of V_ref, t_ref, or index_ref should be supplied. If none are supplied,
        the default reference spectrum is used.

        Args:
            J_name (str): Name of the cycle series
            cycle_1 (int) : the cycle number for the first cycle of interest. Defaults to
                cycle 1.
            cycle_2 (int) : the cycle number for the second cycle of interest. Defaults
                to cycle 2.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Returns ValueSeries: The difference spectrum of the dOD of the two cycles , and
            the difference CVs of the two cycles (CyclicVoltammogramDiff).
        """

        measurement = self

        # Get optical data for the two cycles:

        cycle1=measurement.select(cycle=cycle_1)
        cycle2=measurement.select(cycle=cycle_2)
        if not V_ref and not t_ref and not index_ref:
            cycle1.set_reference_spectrum(measurement.reference_spectrum)
            cycle2.set_reference_spectrum(measurement.reference_spectrum)
        dODcycle1=cycle1.calc_dOD(V_ref=V_ref, t_ref=t_ref, index_ref=index_ref)
        dODcycle2=cycle2.calc_dOD(V_ref=V_ref, t_ref=t_ref, index_ref=index_ref)


        # Make sure fields contain same number of points

        t1 = cycle1.spectra.axes_series[0].data
        t2 = cycle2.spectra.axes_series[0].data

        # convert to relative cycle time and normalise
        t1_norm = (t1 - t1[0]) / (t1[-1] - t1[0])
        t2_norm = (t2 - t2[0]) / (t2[-1] - t2[0])

        interp = interp1d(
            t2_norm,
            dODcycle2.data,
            axis=0,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

        cycle2_dOD_interp = interp(t1_norm)

        difference = (cycle2_dOD_interp - dODcycle1.data) * 1000

        diff_field = Field(
            name="Irreversible ΔA",
            unit_name="mOD",
            data=difference,
            axes_series=cycle1.spectra.axes_series,
        )

        return diff_field

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
        if self.spectra_type is None or self.spectra_type == "Intensity":
            if V_ref is not None or t_ref is not None or index_ref is not None:
                spectrum_ref = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
            else:
                spectrum_ref = self.reference_spectrum
        x = self.wl
        if width:  # averaging
            wl_mask = np.logical_and(wl - width / 2 < x, x < wl + width / 2)
            counts_wl = np.mean(self.spectra.data[:, wl_mask], axis=1)
            if self.spectra_type is None or self.spectra_type == "Intensity":
                counts_ref = np.mean(spectrum_ref.y[wl_mask])

        else:  # interpolation
            if self.spectra_type is None or self.spectra_type == "Intensity":
                counts_ref = np.interp(wl, spectrum_ref.x, spectrum_ref.y)
            counts_wl = []
            for counts_i in self.spectra.data:
                c = np.interp(wl, x, counts_i)
                counts_wl.append(c)
            counts_wl = np.array(counts_wl)
        if self.spectra_type == "Absorption":
            dOD_wl = counts_wl
        elif self.spectra_type == "Transmission":
            counts_wl_dec = counts_wl / 100
            dOD_wl = -np.log10(counts_wl_dec)
        else:
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
