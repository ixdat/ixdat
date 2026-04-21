import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

from .ec import ECMeasurement
from ..db import PlaceHolderObject
from ..spectra import Spectrum, SpectroMeasurement,SpectrumSeries
from ..data_series import Field, ValueSeries
from ..exporters import SECExporter
from ..plotters import SECPlotter, ECOpticalPlotter,StaircaseSECPlotter



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
    pass

class ECOpticalMeasurement(SpectroECMeasurement):
    """Electrochemistry with optical Spectroscopy

    This adds, to the SpectroElectrochemistry base class, methods for normalizing to a
    reference spectrum to get optical density, and for tracking intensity at specific
    wavelengths.
    """

    default_plotter = ECOpticalPlotter

    extra_linkers = SpectroECMeasurement.extra_linkers.copy()
    extra_linkers.update({"ec_optical_measurements": ("spectra", "ref_id")})

    def __init__(self, reference_spectrum=None, ref_id=None, **kwargs):
        """Initialize an SEC measurement. All args and kwargs go to ECMeasurement."""
        SpectroECMeasurement.__init__(self, **kwargs)
        if reference_spectrum:
            self._reference_spectrum = reference_spectrum
        elif ref_id:
            self._reference_spectrum = PlaceHolderObject(ref_id, cls=Spectrum)
        self.tracked_wavelengths = []
        self.plot_waterfall = self.plotter.plot_waterfall
        self.plot_wavelengths = self.plotter.plot_wavelengths
        self.plot_wavelengths_vs_potential = self.plotter.plot_wavelengths_vs_potential
        self.technique = "EC-Optical"

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
        if t_ref and not spectrum:
            spectrum = self.get_spectrum(t=t_ref)
        if V_ref and not spectrum:
            spectrum = self.get_spectrum(V=V_ref)
        if not spectrum:
            raise ValueError("must provide a spectrum, t_ref, or V_ref!")
        self._reference_spectrum = spectrum
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
        if not spectrum:
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
            if interpolate:
                counts_interpolater = interp1d(
                    self.U, counts, axis=0, fill_value=end_spectra, bounds_error=False
                )
                # FIXME: This requires that potential and spectra have same tseries!
                y = counts_interpolater(V)
            else:
                U_diff = np.abs(self.U - V)
                index = np.argmin(U_diff)
                y = counts[index]
            name = name or f"{self.spectra.name}_{V}V"
        elif t:
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


class StaircaseSEC:
    default_plotter = StaircaseSECPlotter

    def __init__(
        self,
        SEC_intensity_dataFrame,
        reverse=False,
        ohmic_corrected = False,
        auto_reverse_correction = True,
        global_WL_interval=None,
        global_V_interval=None,
        global_smooth_param=None,
        V_correction=0,
    ):
        """
        Args:
            SEC_intensity_dataFrame: a dataframe of the raw intensities 
            reverse: whether its a cathodic scan (automatic detection is in use in the current version)
            global_WL_interval: wavelength interval
            global_V_interval: voltage interval
            global_smooth_param: smooth parameters [windon length, polyorder]
            V_correction: apperent voltage = raw voltage + V_correction
        """

        self.SEC_intensity_dataFrame = SEC_intensity_dataFrame.copy()
        # for reverse scan, simply reverse the order of columns, then the columns become ascending order, making reverse scan not different from forward scan
        if self.SEC_intensity_dataFrame.columns[-1] < self.SEC_intensity_dataFrame.columns[0] and auto_reverse_correction:
            reverse = True
        if reverse:
            column_reverse = []
            df_new = self.SEC_intensity_dataFrame.copy()
            for i in range(len(self.SEC_intensity_dataFrame.columns)):
                df_new.iloc[:, i] = self.SEC_intensity_dataFrame.iloc[:, -(i + 1)]
                column_reverse.append(self.SEC_intensity_dataFrame.columns[-(i + 1)])
            df_new.columns = column_reverse
            self.SEC_intensity_dataFrame = df_new

        self.SEC_intensity_dataFrame.columns = (
            np.round(pd.to_numeric(self.SEC_intensity_dataFrame.columns), 3) + V_correction
        )

        self.SEC_raw = self.SEC_intensity_dataFrame.copy()

        reference = self.SEC_intensity_dataFrame[self.SEC_intensity_dataFrame.columns[0]].copy()

        self.SEC_dataFrame = self.SEC_intensity_dataFrame.copy()
        self.ohmic_corrected = ohmic_corrected

        for i in self.SEC_dataFrame.columns:
            self.SEC_dataFrame[i] = -np.log10(self.SEC_intensity_dataFrame[i] / reference)
        self.SEC_original = self.SEC_dataFrame.copy()

        # the SEC_dataFrame after initialization is the unsmoothed data with its first voltage as reference with no applied wl and v truncate
        if global_WL_interval is not None:
            self.set_global_wavelength_interval(WL_interval=global_WL_interval)
        if global_V_interval is not None:
            self.set_global_voltage_interval(V_interval=global_V_interval)
        if global_smooth_param is not None:
            self.SEC_dataFrame = self.smooth(
                self.SEC_dataFrame,
                window_length=global_smooth_param[0],
                polyorder=global_smooth_param[1],
            )
        self.init_plotter()
    def init_plotter(self):
        """
        initiate the plotter
        """
        self.plotter = self.__class__.default_plotter(measurement=self)
        self.plot_waterfall = self.plotter.plot_waterfall
        self.plot_differential = self.plotter.plot_differential
        self.plot_spectrum = self.plotter.plot_spectrum
        self.plot_stack_differential = self.plotter.plot_stack_differential
        

    def smooth(self, input_dataframe, window_length=100, polyorder=1):
        """smooth spectra

        Args:
            input_dataframe : the input dataframe
            window_length : window length for savgol_filter
            polyorder : polyorder for savgol_filter

        Returns:
            smoothed dataframe

        """
        dataframe = input_dataframe.copy()
        for i in dataframe.columns:
            dataframe[i] = savgol_filter(dataframe[i], window_length, polyorder)
        return dataframe
    def spectroLSV_smooth(self, input_dataframe, window_length=100, polyorder=1,loop_times = 1):
        """smooth spetro LSV

        Args:
            input_dataframe : the input dataframe
            window_length : window length for savgol_filter
            polyorder : polyorder for savgol_filter

        Returns:
            smoothed dataframe

        """
        dataframe = input_dataframe.copy()

        for i in dataframe.index:
            for j in range(loop_times):
                dataframe.loc[i,:] = savgol_filter(dataframe.loc[i,:], window_length, polyorder)
        return dataframe

    def with_wavelength_interval(self, WL_interval, input_dataframe):
        """truncate the dataframe based a wavelength interval

        Args:
            WL_interval : [min_wavelength,max_wavelength]
            input_dataframe :  the input dataframe

        Returns:
            output dataframe
        """
        dataframe = input_dataframe.loc[
            (input_dataframe.index >= WL_interval[0])
            & (input_dataframe.index <= WL_interval[1]),
            :,
        ]
        return dataframe

    def with_voltage_interval(self, V_interval, input_dataframe):
        """truncate the dataframe based a voltage interval

        Args:
            V_interval : [min_voltage,max_voltage]
            input_dataframe :  the input dataframe

        Returns:
            output dataframe
        """
        dataframe = input_dataframe.loc[
            :,
            (input_dataframe.columns >= V_interval[0])
            & (input_dataframe.columns <= V_interval[1]),
        ]
        return dataframe

    def set_global_wavelength_interval(self, WL_interval):
        """
        set the global wavelenth interval

        Args:
            WL_interval: wavelength interval
    
        """
        # Apply the wavelength interval
        self.SEC_dataFrame = self.SEC_dataFrame.loc[
            (self.SEC_dataFrame.index >= WL_interval[0])
            & (self.SEC_dataFrame.index <= WL_interval[1]),
            :,
        ]

    def set_global_voltage_interval(self, V_interval):
        """
        set the global voltage interval

        Args:
            V_interval: voltage interval
    
        """
        # Apply the voltage interval
        self.SEC_dataFrame = self.SEC_dataFrame.loc[
            :,
            (self.SEC_dataFrame.columns >= V_interval[0])
            & (self.SEC_dataFrame.columns <= V_interval[1]),
        ]

    def set_global_reference_spectrum(self, V_ref):
        """
        set the global ref spectrum (baseline) to that of a given voltage

        Args:
            V_ref: voltage in volts
    
        """
        if V_ref is not False:
            ref_spectra = self.SEC_dataFrame[self.around_V(V_ref)]
            for i in self.SEC_dataFrame.columns:
                self.SEC_dataFrame[i] = self.SEC_dataFrame[i] - ref_spectra
        else:
            self.SEC_dataFrame = self.SEC_raw.copy()

    def data_config(self, dataFrame, WL_interval=None, V_interval=None):
        """a function both truncate the dataframe based a voltage interval and wavelength interval

        Args:
            WL_interval : [min_wavelength,max_wavelength]
            V_interval : [min_voltage,max_voltage]
            input_dataframe :  the input dataframe

        Returns:
            output dataframe
        """

        if WL_interval is not None:
            dataFrame = self.with_wavelength_interval(
                input_dataframe=dataFrame, WL_interval=WL_interval
            )
        if V_interval is not None:
            dataFrame = self.with_voltage_interval(
                input_dataframe=dataFrame, V_interval=V_interval
            )
        return dataFrame

    def around_WL(self, wavelength):
        """get the wavelength that exist for real data that is closest to a given wavelength

        Args:
            wavelength : a given wavelength

        Returns:
            the closest real wavelength present in the data
        """
        diff = np.abs(np.array(self.SEC_dataFrame.index) - wavelength)
        return self.SEC_dataFrame.index[np.argmin(diff)]

    def around_V(self, voltage):
        """get the voltage that exist for real data that is closest to a given voltage

        Args:
            voltage : a given voltage

        Returns:
            the closest real voltage present in the data
        """
        diff = np.abs(np.array(self.SEC_dataFrame.columns) - voltage)
        return self.SEC_dataFrame.columns[np.argmin(diff)]

    @property
    def SEC_dataFrame_column_axis_label(self):
        if self.ohmic_corrected:
            return "U-IR vs RHE (V)"
        else:
            return "U vs RHE (V)"