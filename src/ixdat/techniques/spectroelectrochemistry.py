import numpy as np
from scipy.interpolate import interp1d

from .ec import ECMeasurement
from ..db import PlaceHolderObject
from ..spectra import Spectrum, SpectroMeasurement, SpectrumSeries
from ..data_series import Field, ValueSeries, TimeSeries
from ..exporters import SECExporter
from ..plotters import SECPlotter, ECOpticalPlotter
from scipy.spatial.distance import pdist, squareform
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d
import copy


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

    def average_spectra(
        self,
        average_number=10,
        field=None,
    ):
        """Return a ValueSeries of averaged spectra for either an OpticalSpectrumSeries or specific field containing spectral data.

        Args:
            average_number (int) : The number of spectra to average. 10 by default.
            field (field) : The spectra data to average, if not included in the measurement
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
            averaged_spectra.append(spectra[index : index + average_number].mean(axis=0))
            averaged_time.append(time[index : index + average_number].mean(axis=0))

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
            spectra_type=spectrum.spectra_type,
        )

        return averaged_uvvis_series
    


class ECOpticalMeasurement(SpectroECMeasurement):
    """Electrochemistry with optical Spectroscopy

    This adds, to the SpectroElectrochemistry base class, methods for normalizing to a
    reference spectrum to get optical density, and for tracking intensity at specific
    wavelengths.
    """

    default_plotter = ECOpticalPlotter

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

        if spectra_type == None:
            spectra_type = self.spectrum_series.spectra_type
        self.spectra_type = spectra_type

        if reference_spectrum:
            self._reference_spectrum = reference_spectrum
        elif ref_id:
            self._reference_spectrum = PlaceHolderObject(ref_id, cls=Spectrum)
        self.tracked_wavelengths = []
        self.plot_waterfall = self.plotter.plot_waterfall
        self.plot_wavelengths = self.plotter.plot_wavelengths
        self.plot_waterfall_cycle = self.plotter.plot_waterfall_cycle
        self.plot_dOD_cycle_diff = self.plotter.plot_dOD_cycle_diff
        self.plot_dOD_difference_spectra = self.plotter.plot_dOD_difference_spectra
        self.plot_convergent_spectra = self.plotter.plot_convergent_spectra
        self.plot_fit_and_residuals = self.plotter.plot_fit_and_residuals
        self.plot_fit_reconstruction = self.plotter.plot_fit_reconstruction
        self.plot_wavelengths_vs_cv = self.plotter.plot_wavelengths_vs_cv
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
        if t_ref is not None and spectrum is None:
            spectrum = self.get_spectrum(t=t_ref)
        if V_ref is not None and spectrum is None:
            spectrum = self.get_spectrum(V=V_ref)
        if spectrum is None:
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
        if t_range is not None and spectrum is None:
            spectrum = self.get_spectrum(t=t_range)
        if V_ref is not None and spectrum is None:
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

        Provide at most one of V_ref, t_ref, or index. If none are provided the default
        reference spectrum (self.reference_spectrum) will be used.

        Args:
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Return Field: the delta optical density spanning time and wavelength
        """
        counts = self.spectra.data

        # Check what form data is in. By default it is Intensity
        if self.spectra_type == "Absorption":
            dOD = counts
        elif self.spectra_type == "Transmission":
            counts_decimal = counts / 100
            counts_decimal[counts_decimal == 0] = np.nan
            dOD = -np.log10(counts_decimal)
        else:
            if V_ref is not None or t_ref is not None:
                ref_spec = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
            else:
                ref_spec = self.reference_spectrum
            ref = ref_spec.y
            counts[counts == 0] = np.nan
            ref[ref == 0] = np.nan
            ratio = counts / ref
            dOD = -np.log10(ratio)

        dOD = np.ma.masked_invalid(dOD)
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

        spectrum = self.get_spectrum(V=V, t=t, index=index)
        # Check what form data is in. By default it is Intensity
        if self.spectra_type is not None and self.spectra_type != "Intensity":
            if self.spectra_type == "Absorption":
                dOD = spectrum.y
            elif self.spectra_type == "Transmission":
                spec_decimal = spectrum.y / 100
                spec_decimal[spec_decimal == 0] = np.nan
                dOD = -np.log10(spec_decimal)
        else:
            if V_ref is not None or t_ref is not None or index_ref is not None:
                spectrum_ref = self.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
            else:
                spectrum_ref = self.reference_spectrum
            ratio = np.divide(
                spectrum.y,
                spectrum_ref.y,
                out=np.full_like(spectrum.y, np.nan),
                where=(spectrum.y > 0) & (spectrum_ref.y > 0),
            )
            dOD = -np.log10(ratio)

        dOD = np.ma.masked_invalid(dOD)

        field = Field(
            data=dOD,
            name=r"$\Delta$ OD",
            unit_name="",
            axes_series=[self.wavelength],
        )
        return Spectrum.from_field(field)

    def get_dOD_spectrum_diff(
        self,
        V=None,
        t=None,
        index=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        normalise=True,
        step=1,
        wlmin=None,
        wlmax=None,
    ):
        """Return the difference of two delta optical density Spectra given a point and reference point (the one at the given point, and the next one measured).

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
            normalise (bool): Whether or not to normalise spectra. Defaults to true
        Return:
             Spectrum: The difference dOD spectrum as a field.
        """
        measurement = self
        t_spec = measurement.spectra.axes_series[0].t
        if index is None:
            if t is not None:
                index = int(np.argmin(np.abs(t_spec - t)))
            elif V is not None:
                index = int(np.argmin(np.abs(self.U - V)))
            else:
                raise ValueError("Need one of t, V, or index.")

        wl = measurement.wavelength.data
        if wlmin == None:
            wlmin = np.min(wl)
        if wlmax == None:
            wlmax = np.max(wl)
        mask = (wl >= wlmin) & (wl <= wlmax)

        dOD1 = measurement.get_dOD_spectrum(
            index=index, V_ref=V_ref, t_ref=t_ref, index_ref=index_ref
        )
        dOD2 = measurement.get_dOD_spectrum(
            index=index + step, V_ref=V_ref, t_ref=t_ref, index_ref=index_ref
        )

        dOD_diff = dOD2.y - dOD1.y

        if normalise == True:
            dOD_diff_max = np.max((np.abs(dOD_diff[mask])))
            print(dOD_diff_max)
            dOD_diff = dOD_diff / dOD_diff_max

        dOD_diff_field = Field(
            data=dOD_diff,
            name="dOD difference spectrum",
            unit_name="",
            axes_series=[measurement.wavelength],
        )

        return Spectrum.from_field(dOD_diff_field)

    @staticmethod
    def _check_direction(direction):
        """Raise ValueError unless direction is 0 (anodic), 1 (cathodic), or None."""
        if direction is not None and direction not in (0, 1):
            raise ValueError(
                "direction must be 0 (anodic), 1 (cathodic), or None (full cycle)"
            )

    def get_dOD_cycle(
        self,
        J_name="cycle",
        cycle_number=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        direction=None,
        N_points=10
    ):
        """Return a ValueSeries for the dOD for a specific cycle.
        If V_ref, t_ref, or index_ref are provided, they specify what to reference dOD
            to. Otherwise, dOD is referenced to the SECMeasurement's reference_spectrum.

        Args:
            J_name (str): Name of the cycle series
            cycle_number (int) : the cycle of interest
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
            direction (0, 1): the scan direction. 0 = anodic, 1 = cathodic, nothing is full cycle. Full cycle by default.
        Returns ValueSeries: The dOD value of the spectrum at wl.
        """
        self._check_direction(direction)
        measurement = self
        # get cycle values
        cycle_series = measurement[J_name]

        cycle = np.interp(
            measurement.spectra.axes_series[0].t,
            cycle_series.tseries.t,
            cycle_series.data,
        )
        # if in-between cycles, round up
        cycle = np.ceil(cycle)

        # Mask spectra belonging to this cycle. Get cycle 1 by default.
        if cycle_number is not None:
            mask = cycle == cycle_number
        else:
            mask = cycle == 1
        dOD = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref, index_ref=index_ref)
        dOD_data = dOD.data[mask]
        t_spec=measurement.spectra.axes_series[0].data[mask]
        
        if direction == 0 or direction == 1:
            U_interp = np.interp(t_spec, measurement.t, measurement.U)
            dUdt = np.gradient(U_interp)
            sign = np.sign(dUdt)
            
            turning_indices = np.where(
                sign[:-1] != sign[1:]
            )[0] + 1
                   
            valid_indices = []
            if len(turning_indices) >= 1:
                for idx in turning_indices:
                    if idx<N_points or len(t_spec)-idx<N_points:
                        continue
                    next_points = sign[idx:idx + N_points]
                
                    if len(next_points) > 0:
                        if len(next_points) > 0:
                            same_sign = (np.all(next_points==sign[idx]))
                            if same_sign and sign[idx]!=0:
                                valid_indices.append(idx)
                            
            #For SEC, you should only have 1 turning point per cycle
            if len(valid_indices)>1:
                raise ValueError(
                    "Multiple turning points detected. Set cycle limits so that there is only one turning point per cycle"
                )
            if len(valid_indices) == 0:
                raise ValueError(
                    "No valid turning point detected in this cycle."
                )
            
            valid_index=valid_indices[0]
            
            if sign[valid_index+1]>0 and sign[valid_index-1]<0:
                anodic_mask=np.arange(len(t_spec))>valid_index
                cathodic_mask=np.arange(len(t_spec))<valid_index
                
            elif sign[valid_index+1]<0 and sign[valid_index-1]>0:
                anodic_mask=np.arange(len(t_spec))<valid_index
                cathodic_mask=np.arange(len(t_spec))>valid_index
            
            if direction==0:
                t_spec=t_spec[anodic_mask]
                dOD_data=dOD_data[anodic_mask]
            elif direction==1:
                t_spec=t_spec[cathodic_mask]
                dOD_data=dOD_data[cathodic_mask]

        tseries = TimeSeries(
            name=measurement.spectra.axes_series[0].name,
            unit_name=measurement.spectra.axes_series[0].unit_name,
            data=t_spec,
            tstamp=measurement.tstamp,
        )

        dOD_cycle = Field(
            name=dOD.name,
            unit_name=dOD.unit_name,
            data=dOD_data,
            axes_series=[
                tseries,
                measurement.spectra.axes_series[1],
            ],
        )
        return dOD_cycle

    def get_dOD_difference_spectra(
        self,
        J_name="cycle",
        cycle_number=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        direction=None,
        normalise=True,
        wlmin=None,
        wlmax=None,
        step=1,
        cycle_field=None,
        N_points=None
    ):
        """Return a ValueSeries of dOD difference spectra for a specific cycle.
        If V_ref, t_ref, or index_ref are provided, they specify what to reference dOD
            to. Otherwise, dOD is referenced to the SECMeasurement's reference_spectrum.

        Args:
            J_name (str): Name of the cycle series. Defaults to "cycle"
            cycle_number (int) : the cycle of interest. Defaults to 1.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
            direction (0, 1): the scan direction. 0 = anodic, 1 = cathodic, nothing is full cycle. Full cycle by default.
            normalise (bool): Whether or not to normalise spectra. Defaults to true
            wlmin (float): minimum wavelength to consider
            wlmax (float): maximum wavelength to consider
            step (int): the step size for difference spectra. 1 by default
        Returns ValueSeries: The difference dOD value of the spectra in the cycle.
        """
        self._check_direction(direction)

        if step < 1 and step is not None:
            raise ValueError("Step must be a positive integer")

        measurement = self
        if cycle_field is not None:
            dOD_cycle = cycle_field
        else:
            dOD_cycle = measurement.get_dOD_cycle(
                J_name=J_name,
                cycle_number=cycle_number,
                V_ref=V_ref,
                t_ref=t_ref,
                index_ref=index_ref,
                direction=direction,
                N_points=N_points
            )

        wl = dOD_cycle.axes_series[1].data

        # trim data to ignore regions outside of detection range before normalising
        if wlmin == None:
            wlmin = np.min(wl)
        if wlmax == None:
            wlmax = np.max(wl)

        dOD_diff = dOD_cycle.data[step::step] - dOD_cycle.data[0:-step:step]

        dOD_diff = np.array(dOD_diff)

        if normalise == True:
            # Select wavelength range
            wl_range = (wl >= wlmin) & (wl <= wlmax)
            max_vals = np.max(dOD_diff[:, wl_range], axis=1)
            dOD_diff = dOD_diff/ max_vals[:, np.newaxis]
            dOD_diff = np.ma.masked_invalid(dOD_diff)
            
            # dOD_diff=dOD_diff[:,wl_range] values outside of wavelength range are non-physical

        indices = np.arange(step, len(dOD_cycle.data), step)
        
        tseries = TimeSeries(
        name=dOD_cycle.axes_series[0].name,
        unit_name=dOD_cycle.axes_series[0].unit_name,
        data=dOD_cycle.axes_series[0].data[indices],
        tstamp=dOD_cycle.axes_series[0].tstamp,
        )

        field = Field(
            data=dOD_diff,
            name=r"$\Delta$ A (U-$\delta$U)(normalised)",
            unit_name="",
            axes_series=[
                tseries,
                dOD_cycle.axes_series[1],
            ],
        )

        return field

    def get_convergent_spectra(
        self,
        J_name="cycle",
        cycle_number=None,
        diff_spectra_field=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        direction=None,
        normalise=True,
        wlmin=400,
        wlmax=900,
        step=1,
        conv_limit=0.01,
        convergence_metric="correlation",
        smooth_distances=False,
        window_length=4,
        polyorder=3,
        prominence=None,
        threshold=None,
        distance=None,
        height=None,
        min_region_width=2,
        min_region_separation=2,
        converged_spectrum_form="average",
        n_start = 3
    ):
        """Return a ValueSeries of convergent dOD spectra.
        If V_ref, t_ref, or index_ref are provided, they specify what to reference dOD
            to. Otherwise, dOD is referenced to the SECMeasurement's reference_spectrum.
            Note: If data is too noisy, it may not meet convergence criterion.

        Args:
            J_name (str): Name of the cycle series. Defaults to "cycle"
            cycle_number (int) : the cycle of interest. Defaults to 1.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
            direction (0, 1): the scan direction. 0 = anodic, 1 = cathodic, nothing is full cycle. Full cycle by default.
            normalise (bool): Whether or not to normalise spectra. Defaults to true
            step (int): the step size for difference spectra. 1 by default
            wlmin (float): minimum wavelength to consider. 400nm by default.
            wlmax (float): maximum wavelength to consider. 900 nm by default
            conv_limt (float): the limit for convergence of spectra. Set to 0.01 by default.
            convergence_metric (str) : the choice of distance algorithm to pass to pdist. Correlation by default.
            smooth_distances (bool) : whether to smooth the distance array (used to find
                similar spectra) with a Savitzky-Golay filter before peak-finding.
            window_length (int) : the window length to pass to the Savitzky-Golay
                filter when smooth_distances is True. 4 by default.
            polyorder (int) : the polynomial order to pass to the Savitzky-Golay
                filter when smooth_distances is True. 3 by default.
            prominence (float): To be passed to find_peaks. How far the peak protrudes from the rest
            threshold (float) : To be passed to find_peaks. The required threshold of peaks, the vertical distance to its neighboring samples
            distance (float) : To be passed to find_peaks. The required minimum horizontal distance in samples between neighbouring peaks. Smaller peaks are removed first until the condition is met.
            height (float) : To be passed to find_peaks. The required height of peaks.
            min_region_width (int) : Minimum number of points in a convergent region. 2 by default.
            min_region_separation (int) : The minimum number of spectra that need to be between two convergent regions.
            converged_spectrum_form (str) : The way the converged spectra are handled. "average" returns an average spectrum over the convergent region.
        Returns: The set of convergent spectra for the specific cycle, the potentials at which they occur, and the plot of the distances used for the convergence criterion
        """
        # FIX: if data is noisy, does not meet convergence criterion. Smoothing doesn't seem to help much either.

        measurement = self
        if diff_spectra_field is None:
            spectra_field = measurement.get_dOD_difference_spectra(
                J_name=J_name,
                cycle_number=cycle_number,
                V_ref=V_ref,
                t_ref=t_ref,
                index_ref=index_ref,
                direction=direction,
                normalise=normalise,
                wlmin=wlmin,
                wlmax=wlmax,
                step=step,
            )
        else:
            spectra_field = diff_spectra_field

        wl = spectra_field.axes_series[1].data
        spectra = spectra_field.data

        # The spectra will be nonsense outside of the region in which we can actually detect signals. For the CHEAC spectrometer, this is between 400 and 900 nm.
        # To make the convergence maths work, we must neglect spectra outside of this region
        wl_range = (wl >= wlmin) & (wl <= wlmax)
        spectra = spectra[:, wl_range]

        # find the stastistical difference between each spectrum
        distance_matrix = squareform(pdist(spectra, convergence_metric))
        # Adjacent distances - the distance between each successive spectrum is on the offset diagonal of the above matrix
        adjacent_distances = np.diag(distance_matrix, k=1)

        # spectra=np.array(spectra)
        # diff_of_spectra=np.diff(spectra, axis=0)
        # adjacent_distances=np.sqrt(np.mean(diff_of_spectra**2, axis=1))

        # Apply smoothing if necessary. Aids peak detection.
        if smooth_distances == True:
            from scipy.signal import savgol_filter
            adjacent_distances = savgol_filter(
                adjacent_distances,
                window_length=window_length,
                polyorder=polyorder
            )
    
        # use find_peaks to identify the minima of the adjacent_distances
        mins, _ = find_peaks(
            -adjacent_distances,
            prominence=prominence,
            threshold=threshold,
            distance=distance,
            height=height,
        )
        # Ensure the minima are at indices where the adjacent_distance is below the convergence limit
        mins = mins[adjacent_distances[mins] < conv_limit]
        
        
        #check to see if first point is a miminum

        if len(adjacent_distances) >= n_start:
            start_gradient = np.gradient(adjacent_distances[:n_start])
        
            if (
                adjacent_distances[0] < conv_limit
                and np.all(start_gradient > 0)
            ):
                mins = np.insert(mins, 0, 0)

        # The convergent spectra will likely persist over a timeframe. Noisy data may incraese this timespan. Identify the regions where spectra are convergent.
        regions = np.split(mins, np.where(np.diff(mins) > min_region_separation)[0] + 1)
        
            
        # Discard regions of smaller width than min_region_width (noise)
        regions = [r for r in regions if len(r) >= min_region_width]

        # Find the timespans for the regions of convergence
        time = spectra_field.axes_series[0].data
        timespans = []
        potentialspans = []
        time_convergence = []
        potentials = []
        convergent_spectra = []
        for r in regions:
            start_time = time[r[0]]
            end_time = time[r[-1] + 1]
            times, potentials_array = measurement.grab(
                "potential", tspan=(start_time, end_time)
            )
            if converged_spectrum_form == "average":
                convergent_spectra.append(spectra[r[0] : r[-1] + 1].mean(axis=0))
                potentials.append(potentials_array.mean())
                time_convergence.append(times.mean())
            if converged_spectrum_form == "raw":
                convergent_spectra.append(spectra[r])
                potentials.append(potentials_array)
                time_convergence.append(times)
            elif (
                converged_spectrum_form is not None
                and converged_spectrum_form != "average"
            ):
                raise ValueError("converged_spectrum_form must be average or raw")
            timespans.append((start_time, end_time))
            potentialspans.append((potentials_array[0], potentials_array[-1]))
        
        if normalise == True:
            convergent_spectra = np.asarray(convergent_spectra)
            # Select wavelength range
            min_vals = np.min(convergent_spectra, axis=1)
            max_vals = np.max(convergent_spectra, axis=1)
            ranges = max_vals - min_vals
            ranges[ranges == 0] = np.nan
            dOD_diff = (convergent_spectra - min_vals[:, np.newaxis]) / ranges[:, np.newaxis]
            dOD_diff = np.ma.masked_invalid(dOD_diff)
        
        
        wl_series = copy.copy(spectra_field.axes_series[1])
        wl_series._data = wl_series.data[wl_range]

        tseries = TimeSeries(
            name="time",
            unit_name=spectra_field.axes_series[0].unit_name,
            data=time_convergence,
            tstamp=spectra_field.axes_series[0].tstamp,
        )

        convergent_spectra = Field(
            data=convergent_spectra,
            name=r"$\Delta$ A (U-$\delta$U)(normalised)",
            unit_name="",
            axes_series=[tseries, wl_series],
        )

        return (
            convergent_spectra,
            potentials,
            adjacent_distances,
            timespans,
            potentialspans,
        )

    def get_dOD_cycle_noise(
        self,
        J_name="cycle",
        cycle_1=None,
        cycle_2=None,
        V_ref1=None,
        t_ref1=None,
        index_ref1=None,
        V_ref2=None,
        t_ref2=None,
        index_ref2=None,
        direction=None,
        N_points=None
    ):
        """Takes an ECOpticalMeasurement and returns the irreversible changes in absorbance between the two cycles.
        Used to check for system stability.

        Args:
            J_name (str): Name of the cycle series
            cycle_1 (int) : the cycle number for the first cycle of interest. Defaults to cycle 1.
            cycle_2 (int) : the cycle number for the second cycle of interest. Defaults to cycle 2.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Returns ValueSeries: The difference spectrum of the dOD of the two cycles , and the difference CVs of the two cycles (CyclicVoltammogramDiff).
        """

        measurement = self

        if cycle_1 == None:
            cycle_1 = 1
        if cycle_2 == None:
            cycle_2 = 2

        # Get optical data for the two cycles:

        dOD_cycle_1 = measurement.get_dOD_cycle(
            J_name=J_name,
            V_ref=V_ref1,
            t_ref=t_ref1,
            index_ref=index_ref1,
            cycle_number=cycle_1,
            direction=direction,
            N_points=N_points
        )
        dOD_cycle_2 = measurement.get_dOD_cycle(
            J_name=J_name,
            V_ref=V_ref2,
            t_ref=t_ref2,
            index_ref=index_ref2,
            cycle_number=cycle_2,
            direction=direction,
            N_points=N_points
        )

        # Make sure fields contain same number of points

        t1 = dOD_cycle_1.axes_series[0].t
        t2 = dOD_cycle_2.axes_series[0].t

        # convert to relative cycle time and normalise
        t1_norm = (t1 - t1[0]) / (t1[-1] - t1[0])
        t2_norm = (t2 - t2[0]) / (t2[-1] - t2[0])

        interp = interp1d(
            t2_norm,
            dOD_cycle_2.data,
            axis=0,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

        dOD_cycle_2_interp = interp(t1_norm)

        difference = (dOD_cycle_2_interp - dOD_cycle_1.data) * 1000

        diff_field = Field(
            name="Irreversible ΔA (mOD)",
            unit_name="",
            data=difference,
            axes_series=dOD_cycle_1.axes_series,
        )

        return diff_field

    def fit_convergent_spectra(
        self,
        cycle_data=None,
        converged_data=None,
        wlmin=400,
        wlmax=900,
        bounds_bool=False,
        noise_bound=0,
    ):
        """Return a ValueSeries of fitting coefficients for a set of convergent spectra found for a measurement.
        cycle_data and converged_data must be provided
            Note: If data is too noisy, it may not meet convergence criterion.
            Note 2: Ensure wlmin and wlmax are consistent between fit_convergent_spectra and get_convergent_spectra

        Args:
            cycle_data=output of get_dOD_cycle
            converged_data=output of get_convergent_spectra
            wlmin (float): minimum wavelength to consider. 400nm by default.
            wlmax (float): maximum wavelength to consider. 900 nm by default
            bounds_bool (bool) : switch to apply a physicality constraint such that the components must increase as the spectral intensity is always increasing
                                The strictness of this constraint removes some flexibility so in real data one might relax it to say that the next component can only be the noise threshold below the previous component.
                                False by default
            noise_bound (float) : the noise bound of the system. 0 by default
        Returns: The set of fitting coefficients, the reconstructed spectra from the fit, and the residuals for the spectra
        """
        from scipy.optimize import lsq_linear

        measurement = self

        if cycle_data == None or converged_data == None:
            raise ValueError("Cycle_data and converged_data must be provided")

        # make sure axes match
        cycle = cycle_data.data
        cycle = np.ma.filled(cycle_data.data, 0)  # NaN will break lsq
        wl = cycle_data.axes_series[1].data
        wl_range = (wl >= wlmin) & (wl <= wlmax)
        cycle = cycle[:, wl_range]
        wl = wl[wl_range]

        converged_spectra_field = converged_data[0]
        converged_spectra = np.array(converged_spectra_field.data)

        converged_spectra = np.ma.filled(converged_spectra, 0)  # NaN will break lsq
        coefficients = []

        for i, spectrum in enumerate(cycle):
            # if i>0 the lower bound of each component is the fit results of the previous iteration
            if bounds_bool == True:
                if i > 0:
                    bounds = (
                        np.asarray(coefficients[i - 1]) - noise_bound / 1000,
                        np.inf,
                    )
                else:
                    bounds = (0, np.inf)
                # print("bounds=", bounds)
                result = lsq_linear(converged_spectra.T, spectrum, bounds=bounds)

            else:
                result = lsq_linear(converged_spectra.T, spectrum)
            coefficients.append(np.asarray(result.x))

        coefficients = np.array(coefficients)

        coefficients_valueseries = ValueSeries(
            data=coefficients,
            name=r"Fitting coefficients",
            unit_name="(arb. units)",
            tseries=cycle_data.axes_series[0],
        )

        # get reconstructed spectra

        reconstructed = coefficients @ converged_spectra

        wl_series = copy.copy(cycle_data.axes_series[1])
        wl_series._data = wl_series.data[wl_range]

        reconstructed_field = Field(
            name="Reconstructed ΔA (OD)",
            unit_name="",
            data=reconstructed,
            axes_series=(cycle_data.axes_series[0], wl_series),
        )

        # get residuals

        residuals = cycle - reconstructed
        residuals = residuals * 1000

        residuals_field = Field(
            name="ΔA residuals (mOD)",
            unit_name="",
            data=residuals,
            axes_series=(cycle_data.axes_series[0], wl_series),
        )

        return coefficients_valueseries, reconstructed_field, residuals_field

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
        denoise_method="PCA",
        PCA_explained_variance=0.99,
        sg_window=10,
        sg_poly_order=3,
        normalise=False,
        wlmin=None,
        wlmax=None
    ):
        """Return a ValueSeries of denoised spectra for either a specific field containing spectral data or for the spectra in a measurement.

        Args:
            spectra_field (str) : the field containing the data for denoising
            denoise_method (str): The method to use for denoising spectra. PCA (principal component analysis) by default. See https://www.nature.com/articles/s41598-018-32713-7
            PCA_n_components (float): The value of explained variance to pass to PCA. 0.99 by default (use the number of components such that the amount of variance that needs to be explained is greater than 99%)
            sg_window (int) : the window size to use for the Savitzky-Golay filter. 10 by default.
            sg_poly_order (int) : the order of the polynomial used to fit the data in the window. 3 by default.
        Returns ValueSeries: The field containing denoised spectra.
        """
        
        if spectra_field==None:
            spectra_field=self.spectra
        
        spectra=spectra_field.data

        if denoise_method == "PCA":
            from sklearn.decomposition import PCA

            pca = PCA()
            pca.fit(spectra)
            explained_variance = np.cumsum(pca.explained_variance_ratio_)
            components = np.argmax(explained_variance >= PCA_explained_variance) + 1

            analysis = PCA(n_components=components)
            scores = analysis.fit_transform(spectra)
            spectra_denoised = analysis.inverse_transform(scores)

        if denoise_method == "Savitzky-Golay":
            from scipy.signal import savgol_filter

            spectra_denoised = savgol_filter(
                spectra, window_length=sg_window, polyorder=sg_poly_order, axis=1
            )
        if normalise == True:
            # Select wavelength range
            wl = spectra_field.axes_series[1].data
            # trim data to ignore regions outside of detection range before normalising
            if wlmin == None:
                wlmin = np.min(wl)
            if wlmax == None:
                wlmax = np.max(wl)
            wl_range = (wl >= wlmin) & (wl <= wlmax)
            
            min_vals = np.min(spectra_denoised[:, wl_range], axis=1)
            max_vals = np.max(spectra_denoised[:, wl_range], axis=1)
            ranges = max_vals - min_vals
            ranges[ranges == 0] = np.nan
            spectra_denoised = (spectra_denoised - min_vals[:, np.newaxis]) / ranges[:, np.newaxis]
            spectra_denoised = np.ma.masked_invalid(spectra_denoised)

        spectra_denoised = Field(
            data=spectra_denoised,
            name=r"$\Delta$ A (U-$\delta$U)(normalised, denoised)",
            unit_name="",
            axes_series=spectra_field.axes_series,
        )

        return spectra_denoised

    
