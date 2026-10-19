"""Calculators for spectroelectrochemistry (SEC) measurements"""

import warnings

import numpy as np
from scipy.interpolate import interp1d

from ..measurement_base import Calculator
from ..data_series import DataSeries, Field, TimeSeries, ValueSeries
from ..exceptions import QuantificationError


def get_x_mask(x, xspan=None):
    """Return a boolean mask of the values of `x` within `xspan` (inclusive).

    Args:
        x (np.array): The values to mask, typically wavelengths
        xspan (iterable): [x_min, x_max]. Defaults to the full range of `x`.
    """
    if xspan is None:
        return np.ones_like(x, dtype=bool)
    return (x >= xspan[0]) & (x <= xspan[1])


def normalise_spectra(spectra, x, xspan=None):
    """Return each spectrum divided by its maximum absolute value within xspan.

    Args:
        spectra (np.array): 2-D array with one spectrum per row
        x (np.array): The x values (e.g. wavelengths) of the spectra
        xspan (iterable): [x_min, x_max] in which to find the maximum.
    """
    mask = get_x_mask(x, xspan)
    max_vals = np.nanmax(np.abs(spectra[:, mask]), axis=1)
    max_vals[max_vals == 0] = np.nan
    return spectra / max_vals[:, np.newaxis]


def calc_cycle_difference(measurement, cycle_1=1, cycle_2=2, c_name="cycle"):
    """Return the change in dOD from one cycle to another, to check for stability.

    Both cycles are selected from the measurement and their dOD is calculated with
    the measurement's reference spectrum. The dOD of `cycle_2` is interpolated onto
    the times of `cycle_1`, with each cycle's time scaled from 0 at its start to 1 at
    its end, and the dOD of `cycle_1` is subtracted. A reversible system gives zero.

    Args:
        measurement (ECOpticalMeasurement): The measurement with both cycles
        cycle_1 (int): The first cycle. Defaults to 1.
        cycle_2 (int): The second cycle. Defaults to 2.
        c_name (str): The name of the cycle series. Defaults to "cycle".
    Return Field: The dOD of cycle_2 minus that of cycle_1 in [mOD], on the time and
        wavelength axes of cycle_1
    """
    cycle1 = measurement.select(**{c_name: cycle_1})
    cycle2 = measurement.select(**{c_name: cycle_2})
    dOD1 = cycle1["dOD"]
    dOD2 = cycle2["dOD"]

    # Scale each cycle's time to [0, 1], so the cycles can be compared point by point
    t1 = dOD1.axes_series[0].t
    t2 = dOD2.axes_series[0].t
    t1_norm = (t1 - t1[0]) / (t1[-1] - t1[0])
    t2_norm = (t2 - t2[0]) / (t2[-1] - t2[0])
    dOD2_interp = interp1d(
        t2_norm,
        np.ma.filled(np.ma.asarray(dOD2.data, dtype=float), np.nan),
        axis=0,
        bounds_error=False,
        fill_value=np.nan,
    )(t1_norm)
    difference = (dOD2_interp - dOD1.data) * 1000

    return Field(
        name="Irreversible ΔA",
        unit_name="mOD",
        data=np.ma.masked_invalid(difference),
        axes_series=dOD1.axes_series,
    )


class DeltaODCalculator(Calculator):
    """Calculator for change in optical density (dOD) and its finite difference.

    Attached by default to every `ECOpticalMeasurement` as `dOD_calculator`. It
    provides:

    - "dOD": The change in optical density with respect to the measurement's
      `reference_spectrum`, spanning time and wavelength. See `calc_dOD`.
    - Keys ending in ">diff", e.g. "dOD>diff": The difference between spectra of
      the field named before ">diff", taken at steps in index, time, or potential.
      See `calc_diff`.

    The step and normalisation used by ">diff" lookups are attributes of the
    calculator. A cached lookup is not recalculated when they change, so call
    `measurement.clear_cache()` after changing them, or call `calc_diff` directly.
    """

    calculator_type = "dOD_calculator"
    available_series_names = {"dOD"}
    available_suffixes = {">diff"}
    extra_column_attrs = {
        "dOD_calculators": {
            "normalise",
            "xspan",
            "idx_step",
            "tstep",
            "Vstep",
            "interpolate",
        }
    }
    column_types = {"xspan": list, "idx_step": int, "tstep": float, "Vstep": float}

    def __init__(
        self,
        name=None,
        technique="EC-Optical",
        tstamp=None,
        measurement=None,
        normalise=True,
        xspan=None,
        idx_step=None,
        tstep=None,
        Vstep=None,
        interpolate=True,
    ):
        """Initiate a DeltaODCalculator

        Args:
            name (str): The name of the calculator
            technique (str): The technique of the calculator
            tstamp (float): The time at which the calculator is valid
            measurement (ECOpticalMeasurement): Optional. The measurement to use by
                default
            normalise (bool): Whether to normalise each difference spectrum to its
                maximum absolute value within `xspan`. Defaults to True.
            xspan (iterable): [wl_min, wl_max], the wavelength range in [nm] used
                for normalisation. Defaults to the full range.
            idx_step (int): The step in spectrum index between difference spectra
            tstep (float): The step in time in [s] between difference spectra
            Vstep (float): The step in potential in [V] between difference spectra.
                Requires the potential to be monotonic, i.e. a single sweep.
            interpolate (bool): Whether to interpolate spectra at the steps in time
                or potential. If False, the closest spectrum is used. Defaults to
                True.
        """
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        # bool() because the sqlite backend returns booleans as 0 or 1
        self.normalise = bool(normalise)
        self.xspan = xspan
        self.idx_step = idx_step
        self.tstep = tstep
        self.Vstep = Vstep
        self.interpolate = bool(interpolate)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(normalise={self.normalise}, "
            f"xspan={self.xspan}, idx_step={self.idx_step}, tstep={self.tstep}, "
            f"Vstep={self.Vstep}, interpolate={self.interpolate})"
        )

    def __add__(self, other):
        """Adding a DeltaODCalculator replaces this one with the newer one."""
        return other

    def calculate_series(self, key, measurement=None):
        """Return the field for "dOD" or for a key ending in ">diff"."""
        measurement = measurement or self.measurement
        if key == "dOD":
            return self.calc_dOD(measurement)
        if key.endswith(">diff"):
            return self.calc_diff(measurement, field_name=key.removesuffix(">diff"))
        raise QuantificationError(f"{self} can not calculate {key}")

    def calc_dOD(self, measurement=None, V_ref=None, t_ref=None, index_ref=None):
        """Calculate the optical density with respect to a reference

        Provide at most one of V_ref, t_ref, or index_ref. If none are provided, the
        measurement's `reference_spectrum` is used. If that is missing too, a
        ValueError is raised for intensity data, while absorption and transmission
        data are taken as referenced to the reference used during data collection.

        Args:
            measurement (ECOpticalMeasurement): The measurement. Defaults to
                self.measurement.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
        Return Field: the delta optical density spanning time and wavelength
        """
        measurement = measurement or self.measurement
        spectra = measurement.spectra
        counts = spectra.data
        spectra_type = measurement.spectra_type

        ref = None
        if V_ref is not None or t_ref is not None or index_ref is not None:
            ref = measurement.get_spectrum(V=V_ref, t=t_ref, index=index_ref)
        elif measurement.reference_spectrum is not None:
            ref = measurement.reference_spectrum

        if ref is None:
            if spectra_type == "Absorption":
                warnings.warn(
                    "Absorption data and no reference spectrum given. Optical density"
                    " is referred to the reference used during data collection."
                )
                dOD = counts
            elif spectra_type == "Transmission":
                warnings.warn(
                    "Transmission data and no reference spectrum given. Optical"
                    " density is referred to the reference used during data"
                    " collection."
                )
                # Transmission is in percent
                transmission = counts / 100
                transmission[transmission == 0] = np.nan
                dOD = -np.log10(transmission)
            else:
                raise ValueError(
                    "For intensity data, a reference must be provided to find the"
                    " optical density"
                )
        elif spectra_type == "Absorption":
            dOD = counts - ref.y
        else:  # Intensity or Transmission
            ratio = np.divide(
                counts,
                ref.y,
                out=np.full_like(counts, np.nan, dtype=float),
                where=(abs(counts) > 0) & (abs(ref.y) > 0),
            )
            dOD = -np.log10(ratio)

        return Field(
            name="dOD",
            unit_name="",
            axes_series=spectra.axes_series,
            data=np.ma.masked_invalid(dOD),
        )

    def calc_diff(
        self,
        measurement=None,
        field_name="dOD",
        normalise=None,
        xspan=None,
        idx_step=None,
        tstep=None,
        Vstep=None,
        interpolate=None,
    ):
        """Return the differences between spectra at regular steps

        The spectra of the field named `field_name` are sampled at regular steps in
        index, time, or potential, starting from the first spectrum. Each difference
        spectrum is the sampled spectrum minus the one sampled a step earlier, and is
        placed at the time of the later of the two.

        Arguments left as None are taken from the calculator's attributes. Arguments
        that are given update the calculator's attributes. Exactly one of
        `idx_step`, `tstep`, and `Vstep` must be set after this update.

        Args:
            measurement (ECOpticalMeasurement): The measurement. Defaults to
                self.measurement.
            field_name (str): The name of the field to difference. Defaults to "dOD".
            normalise (bool): Whether to normalise each difference spectrum to its
                maximum absolute value within `xspan`.
            xspan (iterable): [wl_min, wl_max], the wavelength range in [nm] used
                for normalisation.
            idx_step (int): The step in spectrum index
            tstep (float): The step in time in [s]
            Vstep (float): The step in potential in [V]. The potential must be
                monotonic, i.e. the measurement a single sweep.
            interpolate (bool): Whether to interpolate spectra at the steps in time
                or potential. If False, the closest spectrum is used.
        Return Field: The difference spectra spanning time and wavelength
        """
        measurement = measurement or self.measurement
        given = dict(
            normalise=normalise,
            xspan=xspan,
            idx_step=idx_step,
            tstep=tstep,
            Vstep=Vstep,
            interpolate=interpolate,
        )
        for attr, value in given.items():
            if value is not None:
                setattr(self, attr, value)
        if any(given[k] is not None for k in ("idx_step", "tstep", "Vstep")):
            # A newly given step replaces the steps set before.
            for k in ("idx_step", "tstep", "Vstep"):
                setattr(self, k, given[k])

        steps = {
            k: getattr(self, k)
            for k in ("idx_step", "tstep", "Vstep")
            if getattr(self, k) is not None
        }
        if len(steps) != 1:
            raise ValueError(
                f"{self} needs exactly one of idx_step, tstep, or Vstep, got {steps}"
            )

        field = measurement[field_name]
        tseries, xseries = field.axes_series
        data = np.ma.filled(np.ma.asarray(field.data, dtype=float), np.nan)
        t = tseries.t

        if self.idx_step is not None:
            if self.idx_step < 1 or self.idx_step >= len(t):
                raise ValueError(
                    f"idx_step={self.idx_step} must be between 1 and the number of"
                    f" spectra ({len(t)})."
                )
            indices = np.arange(0, len(t), self.idx_step)
            sampled = data[indices]
            sample_t = t[indices]
        else:
            if self.tstep is not None:
                x = t
                step = self.tstep
            else:
                # The potential at the time of each spectrum
                x = measurement.grab_for_t(
                    "potential", t=t + tseries.tstamp - measurement.tstamp
                )
                if not (np.all(np.diff(x) > 0) or np.all(np.diff(x) < 0)):
                    raise ValueError(
                        "Vstep requires a monotonic potential. Select a single sweep"
                        " first, e.g. with `redefine_cycle(turning_point=True)` and"
                        " `select(cycle=...)`."
                    )
                # Step in the direction of the sweep:
                step = abs(self.Vstep) * np.sign(x[-1] - x[0])
            if step == 0 or abs(step) > abs(x[-1] - x[0]):
                raise ValueError(
                    f"The step {step} is zero or larger than the range of the data."
                )
            grid = np.arange(x[0], x[-1] + step / 1e6, step)
            kind = "linear" if self.interpolate else "nearest"
            sampled = interp1d(x, data, axis=0, kind=kind)(grid)
            sample_t = interp1d(x, t, kind=kind)(grid)

        diff = sampled[1:] - sampled[:-1]
        if self.normalise:
            diff = normalise_spectra(diff, xseries.data, self.xspan)

        diff_tseries = TimeSeries(
            name=tseries.name,
            unit_name=tseries.unit_name,
            data=sample_t[1:],
            tstamp=tseries.tstamp,
        )
        return Field(
            name=field_name + ">diff",
            unit_name="" if self.normalise else field.unit_name,
            data=np.ma.masked_invalid(diff),
            axes_series=[diff_tseries, xseries],
        )


class Denoiser(Calculator):
    """Calculator for denoised spectra.

    Provides keys ending in ">smooth", e.g. "dOD>smooth" or "dOD>diff>smooth": the
    field named before ">smooth", denoised by principal component analysis (PCA) or a
    Savitzky-Golay filter. See `denoise_spectra`.

    Add it to a measurement with `measurement.add_calculator(Denoiser(...))`. A cached
    lookup is not recalculated when the parameters change, so call
    `measurement.clear_cache()` after changing them, or call `denoise_spectra`
    directly.
    """

    calculator_type = "denoiser"
    available_series_names = set()
    available_suffixes = {">smooth"}
    denoise_methods = ("PCA", "Savitzky-Golay")
    extra_column_attrs = {
        "denoisers": {
            "denoise_method",
            "PCA_explained_variance",
            "sg_window",
            "sg_poly_order",
            "normalise",
            "xspan",
        }
    }
    column_types = {
        "denoise_method": str,
        "PCA_explained_variance": float,
        "sg_window": int,
        "sg_poly_order": int,
        "xspan": list,
    }

    def __init__(
        self,
        name=None,
        technique=None,
        tstamp=None,
        measurement=None,
        denoise_method="PCA",
        PCA_explained_variance=0.99,
        sg_window=10,
        sg_poly_order=3,
        normalise=False,
        xspan=None,
    ):
        """Initiate a Denoiser

        Args:
            name (str): The name of the calculator
            technique (str): The technique of the calculator
            tstamp (float): The time at which the calculator is valid
            measurement (SpectroMeasurement): Optional. The measurement to use by
                default
            denoise_method (str): "PCA" (principal component analysis, default) or
                "Savitzky-Golay". For PCA, see
                https://www.nature.com/articles/s41598-018-32713-7
            PCA_explained_variance (float): The fraction of the variance that the
                kept principal components must explain. 0.99 by default.
            sg_window (int): The window length of the Savitzky-Golay filter, in
                number of wavelengths. 10 by default.
            sg_poly_order (int): The order of the polynomial fitted in each window of
                the Savitzky-Golay filter. 3 by default.
            normalise (bool): Whether to scale each denoised spectrum to the range
                [0, 1] within `xspan`. False by default.
            xspan (iterable): [wl_min, wl_max], the wavelength range in [nm] used
                for normalisation. Defaults to the full range.
        """
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.denoise_method = denoise_method
        self.PCA_explained_variance = PCA_explained_variance
        self.sg_window = sg_window
        self.sg_poly_order = sg_poly_order
        self.normalise = bool(normalise)
        self.xspan = xspan

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(denoise_method={self.denoise_method!r}, "
            f"PCA_explained_variance={self.PCA_explained_variance}, "
            f"sg_window={self.sg_window}, sg_poly_order={self.sg_poly_order}, "
            f"normalise={self.normalise}, xspan={self.xspan})"
        )

    def __add__(self, other):
        """Adding a Denoiser replaces this one with the newer one."""
        return other

    def calculate_series(self, key, measurement=None):
        """Return the denoised field for a key ending in ">smooth"."""
        measurement = measurement or self.measurement
        if key.endswith(">smooth"):
            spectra_field = measurement[key.removesuffix(">smooth")]
            return self.denoise_spectra(spectra_field, name=key)
        raise QuantificationError(f"{self} can not calculate {key}")

    def denoise_spectra(
        self,
        spectra_field,
        denoise_method=None,
        PCA_explained_variance=None,
        sg_window=None,
        sg_poly_order=None,
        normalise=None,
        xspan=None,
        name=None,
    ):
        """Return a field of denoised spectra

        Arguments left as None are taken from the calculator's attributes. Arguments
        that are given update the calculator's attributes.

        Args:
            spectra_field (Field): The spectra to denoise, spanning time and
                wavelength
            denoise_method (str): "PCA" or "Savitzky-Golay"
            PCA_explained_variance (float): The fraction of the variance that the
                kept principal components must explain
            sg_window (int): The window length of the Savitzky-Golay filter
            sg_poly_order (int): The polynomial order of the Savitzky-Golay filter
            normalise (bool): Whether to scale each spectrum to [0, 1] within xspan
            xspan (iterable): [wl_min, wl_max] used for normalisation
            name (str): The name of the returned field. Defaults to the name of
                `spectra_field` with ">smooth" appended.
        Return Field: The denoised spectra, with the axes of `spectra_field`
        """
        given = dict(
            denoise_method=denoise_method,
            PCA_explained_variance=PCA_explained_variance,
            sg_window=sg_window,
            sg_poly_order=sg_poly_order,
            normalise=normalise,
            xspan=xspan,
        )
        for attr, value in given.items():
            if value is not None:
                setattr(self, attr, value)

        spectra = spectra_field.data

        if self.denoise_method == "PCA":
            try:
                from sklearn.decomposition import PCA
            except ImportError as e:
                raise ImportError(
                    "Denoising with PCA uses the scikit-learn package. Install it"
                    " with `pip install scikit-learn` and try again, or use"
                    ' denoise_method="Savitzky-Golay".'
                    f"\noriginal error:\n{e}"
                )
            pca = PCA()
            pca.fit(spectra)
            explained_variance = np.cumsum(pca.explained_variance_ratio_)
            components = (
                np.argmax(explained_variance >= self.PCA_explained_variance) + 1
            )
            analysis = PCA(n_components=components)
            scores = analysis.fit_transform(spectra)
            spectra_denoised = analysis.inverse_transform(scores)
        elif self.denoise_method == "Savitzky-Golay":
            from scipy.signal import savgol_filter

            spectra_denoised = savgol_filter(
                spectra,
                window_length=self.sg_window,
                polyorder=self.sg_poly_order,
                axis=1,
            )
        else:
            raise ValueError(
                f"denoise_method must be one of {self.denoise_methods},"
                f" not {self.denoise_method!r}"
            )

        if self.normalise:
            mask = get_x_mask(spectra_field.axes_series[1].data, self.xspan)
            min_vals = np.min(spectra_denoised[:, mask], axis=1)
            max_vals = np.max(spectra_denoised[:, mask], axis=1)
            ranges = max_vals - min_vals
            ranges[ranges == 0] = np.nan
            spectra_denoised = (spectra_denoised - min_vals[:, np.newaxis]) / ranges[
                :, np.newaxis
            ]

        return Field(
            name=name or spectra_field.name + ">smooth",
            unit_name="" if self.normalise else spectra_field.unit_name,
            data=np.ma.masked_invalid(spectra_denoised),
            axes_series=spectra_field.axes_series,
        )


class SpectrumDeconvolutor(Calculator):
    """Calculator that fits spectra as linear combinations of converged spectra.

    Converged spectra are found where adjacent difference spectra stop changing,
    e.g. where one redox transition is underway. Build the calculator from a
    measurement with `from_spectromeasurement`, which finds the converged spectra and
    stores them with the diagnostics of the search. Then add it to the same or
    another measurement with `add_calculator` and look up:

    - "coefficient_<i>": The fit coefficient of converged spectrum i vs time
    - "reconstructed": The fitted spectra, spanning time and wavelength
    - "residuals": The measured minus the fitted spectra in [mOD]

    The fitted field is `field_name` ("dOD" by default) of the measurement.
    """

    calculator_type = "spectrum_deconvolutor"
    extra_column_attrs = {
        "spectrum_deconvolutors": {
            "converged_data",
            "wavelength",
            "converged_times",
            "converged_tstamp",
            "field_name",
            "xspan",
            "bounds_bool",
            "noise_bound",
            "adjacent_distances",
            "timespans",
            "potentialspans",
        }
    }
    column_types = {
        "converged_data": np.ndarray,
        "wavelength": np.ndarray,
        "converged_times": np.ndarray,
        "converged_tstamp": float,
        "field_name": str,
        "xspan": list,
        "noise_bound": float,
        "adjacent_distances": np.ndarray,
        "timespans": list,
        "potentialspans": list,
    }

    def __init__(
        self,
        converged_data,
        wavelength,
        converged_times=None,
        converged_tstamp=None,
        name=None,
        technique="EC-Optical",
        tstamp=None,
        measurement=None,
        field_name="dOD",
        xspan=None,
        bounds_bool=False,
        noise_bound=0,
        adjacent_distances=None,
        timespans=None,
        potentialspans=None,
    ):
        """Initiate a SpectrumDeconvolutor

        Args:
            converged_data (np.array): 2-D array with one converged spectrum per row
            wavelength (np.array): The wavelengths in [nm] of the converged spectra
            converged_times (np.array): The times at which the converged spectra
                occur, relative to `converged_tstamp`
            converged_tstamp (float): The unix time corresponding to t=0 in
                `converged_times`
            name (str): The name of the calculator
            technique (str): The technique of the calculator
            tstamp (float): The time at which the calculator is valid
            measurement (ECOpticalMeasurement): Optional. The measurement to use by
                default
            field_name (str): The field of the measurement to fit. Defaults to
                "dOD".
            xspan (iterable): [wl_min, wl_max], the wavelength range in [nm] of the
                fit. Defaults to the range of `wavelength`.
            bounds_bool (bool): Whether to constrain each coefficient to be at least
                its value for the previous spectrum minus `noise_bound`, i.e. to
                only increase in time. False by default.
            noise_bound (float): The noise level in [mOD] allowed by the bounds. 0
                by default.
            adjacent_distances (np.array): Diagnostic from `from_spectromeasurement`:
                the distance between each pair of adjacent difference spectra
            timespans (list): Diagnostic: the timespan of each convergent region
            potentialspans (list): Diagnostic: the potential span of each
                convergent region
        """
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.converged_data = np.asarray(converged_data)
        self.wavelength = np.asarray(wavelength)
        self.converged_times = converged_times
        self.converged_tstamp = converged_tstamp
        self.field_name = field_name
        self.xspan = xspan
        self.bounds_bool = bool(bounds_bool)
        self.noise_bound = noise_bound
        self.adjacent_distances = adjacent_distances
        self.timespans = timespans
        self.potentialspans = potentialspans
        # Set by `deconvolute`:
        self.coefficients = None
        self.reconstructed_field = None
        self.residuals_field = None
        self._fitted_measurement = None

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"{len(self.converged_data)} converged spectra, "
            f"field_name={self.field_name!r}, xspan={self.xspan})"
        )

    def __add__(self, other):
        """Adding a SpectrumDeconvolutor replaces this one with the newer one."""
        return other

    @property
    def available_series_names(self):
        names = {f"coefficient_{i}" for i in range(len(self.converged_data))}
        return names | {"reconstructed", "residuals"}

    @property
    def wavelength_series(self):
        return DataSeries(name="wavelength", unit_name="nm", data=self.wavelength)

    @property
    def converged_spectra(self):
        """The converged spectra as a Field spanning time and wavelength"""
        times = self.converged_times
        if times is None:
            times = np.arange(len(self.converged_data), dtype=float)
        tseries = TimeSeries(
            name="time", unit_name="s", data=times, tstamp=self.converged_tstamp
        )
        return Field(
            name="converged spectra",
            unit_name="",
            data=self.converged_data,
            axes_series=[tseries, self.wavelength_series],
        )

    @property
    def spectrum_list(self):
        """The converged spectra as a list of Spectrum objects"""
        from ..spectra import Spectrum

        return [
            Spectrum.from_field(
                Field(
                    name=f"converged spectrum {i}",
                    unit_name="",
                    data=y,
                    axes_series=[self.wavelength_series],
                ),
                tstamp=self.converged_tstamp,
            )
            for i, y in enumerate(self.converged_data)
        ]

    @classmethod
    def from_spectromeasurement(
        cls,
        measurement,
        field_name="dOD>diff",
        normalise=True,
        xspan=None,
        idx_step=None,
        tstep=None,
        Vstep=None,
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
        n_start=3,
        fit_field_name="dOD",
        bounds_bool=False,
        noise_bound=0,
    ):
        """Find the converged spectra of a measurement and return a deconvolutor

        Converged spectra are found in the difference spectra `field_name`
        ("dOD>diff" by default). If any of `idx_step`, `tstep`, `Vstep`, or
        `normalise` and `xspan` are given, they are set on the measurement's
        `dOD_calculator` first.

        The distance between each pair of adjacent difference spectra is
        calculated with scipy's pdist
        (https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html).
        With the default "correlation" metric, a distance of 0 means the two
        spectra have the same shape, and 1 means they are uncorrelated. A
        convergent region is a run of local minima of this distance below
        `conv_limit`. Its converged spectrum is the average of the difference
        spectra in the region. Noisy data may not meet the convergence criterion.

        Args:
            measurement (ECOpticalMeasurement): The measurement
            field_name (str): The difference spectra to search. Defaults to
                "dOD>diff".
            normalise (bool): Whether to normalise the difference spectra and the
                converged spectra to their maximum absolute value within `xspan`
            xspan (iterable): [wl_min, wl_max], the wavelength range in [nm] to
                consider. Spectra outside the detectable range would otherwise
                break the convergence maths.
            idx_step (int): The step in spectrum index of the difference spectra
            tstep (float): The step in time in [s] of the difference spectra
            Vstep (float): The step in potential in [V] of the difference spectra
            conv_limit (float): The maximum distance between adjacent spectra for
                convergence. Dimensionless. 0.01 by default.
            convergence_metric (str): The distance metric passed to pdist.
                "correlation" by default.
            smooth_distances (bool): Whether to smooth the distances with a
                Savitzky-Golay filter before finding their minima
            window_length (int): The Savitzky-Golay window length. 4 by default.
            polyorder (int): The Savitzky-Golay polynomial order. 3 by default.
            prominence (float): Passed to scipy.signal.find_peaks: how far a
                minimum must stand out from its surroundings
            threshold (float): Passed to find_peaks: the required vertical distance
                of a minimum to its neighbours
            distance (float): Passed to find_peaks: the minimum number of samples
                between neighbouring minima
            height (float): Passed to find_peaks: the required depth of minima
            min_region_width (int): The minimum number of minima in a convergent
                region. 2 by default.
            min_region_separation (int): The number of spectra separating two
                convergent regions. 2 by default.
            n_start (int): The number of first distances checked to decide whether
                the first spectrum is a minimum. 3 by default.
            fit_field_name (str): The field the deconvolutor fits. "dOD" by default.
            bounds_bool (bool): See `__init__`
            noise_bound (float): See `__init__`
        Return SpectrumDeconvolutor: with the converged spectra and the diagnostics
            `adjacent_distances`, `timespans`, and `potentialspans`
        """
        from scipy.spatial.distance import pdist, squareform
        from scipy.signal import find_peaks

        diff_settings = dict(idx_step=idx_step, tstep=tstep, Vstep=Vstep)
        if any(v is not None for v in diff_settings.values()):
            spectra_field = measurement.dOD_calculator.calc_diff(
                measurement,
                field_name=field_name.removesuffix(">diff"),
                normalise=normalise,
                xspan=xspan,
                **diff_settings,
            )
        else:
            spectra_field = measurement[field_name]

        wl = spectra_field.axes_series[1].data
        mask = get_x_mask(wl, xspan)
        spectra = np.ma.filled(np.ma.asarray(spectra_field.data, dtype=float), 0)
        spectra = spectra[:, mask]

        # The distance between each pair of spectra. For the correlation metric:
        #   1 - (u - mean(u)) . (v - mean(v)) / (|u - mean(u)| |v - mean(v)|)
        distance_matrix = squareform(pdist(spectra, convergence_metric))
        # Adjacent spectra are on the first off-diagonal
        adjacent_distances = np.diag(distance_matrix, k=1)

        if smooth_distances:
            from scipy.signal import savgol_filter

            adjacent_distances = savgol_filter(
                adjacent_distances, window_length=window_length, polyorder=polyorder
            )

        # The minima of the distances, below the convergence limit
        mins, _ = find_peaks(
            -adjacent_distances,
            prominence=prominence,
            threshold=threshold,
            distance=distance,
            height=height,
        )
        mins = mins[adjacent_distances[mins] < conv_limit]

        # find_peaks never returns the first point; check it separately
        if len(adjacent_distances) >= n_start:
            start_gradient = np.gradient(adjacent_distances[:n_start])
            if adjacent_distances[0] < conv_limit and np.all(start_gradient > 0):
                mins = np.insert(mins, 0, 0)

        # Convergence persists over several spectra, more so for noisy data. Group
        # the minima into regions and drop regions too narrow to be real.
        regions = np.split(mins, np.where(np.diff(mins) > min_region_separation)[0] + 1)
        regions = [r for r in regions if len(r) >= min_region_width]
        if not regions:
            raise ValueError(
                f"No convergent regions found in {field_name} with"
                f" conv_limit={conv_limit}. Try a larger conv_limit, a narrower"
                " xspan, or a larger step between difference spectra."
            )

        tseries = spectra_field.axes_series[0]
        time = tseries.t
        t_offset = tseries.tstamp - measurement.tstamp
        timespans = []
        potentialspans = []
        converged_times = []
        converged_spectra = []
        for r in regions:
            start, end = r[0], r[-1] + 1
            tspan = (time[start] + t_offset, time[end] + t_offset)
            potentials = measurement.grab_for_t("potential", t=np.array(tspan))
            converged_spectra.append(spectra[start:end].mean(axis=0))
            converged_times.append(time[start:end + 1].mean())
            timespans.append(tspan)
            potentialspans.append((potentials[0], potentials[-1]))

        converged_spectra = np.array(converged_spectra)
        if normalise:
            converged_spectra = normalise_spectra(converged_spectra, wl[mask])

        return cls(
            converged_data=converged_spectra,
            wavelength=wl[mask],
            converged_times=np.array(converged_times),
            converged_tstamp=tseries.tstamp,
            field_name=fit_field_name,
            xspan=xspan,
            bounds_bool=bounds_bool,
            noise_bound=noise_bound,
            adjacent_distances=adjacent_distances,
            timespans=timespans,
            potentialspans=potentialspans,
        )

    def calculate_series(self, key, measurement=None):
        """Return a coefficient series, or the reconstructed or residual field."""
        measurement = measurement or self.measurement
        if measurement is self._fitted_measurement:
            coefficients = self.coefficients
        else:
            coefficients = self.deconvolute(measurement)
        if key == "reconstructed":
            return self.reconstructed_field
        if key == "residuals":
            return self.residuals_field
        if key.startswith("coefficient_"):
            return coefficients[int(key.removeprefix("coefficient_"))]
        raise QuantificationError(f"{self} can not calculate {key}")

    def deconvolute(self, measurement=None):
        """Fit each spectrum of the measurement as a sum of the converged spectra

        The spectra of `self.field_name` are interpolated onto the wavelengths of
        the converged spectra, within `self.xspan`, and fitted by linear least
        squares (scipy.optimize.lsq_linear). The reconstructed spectra and the
        residuals are stored as `self.reconstructed_field` and
        `self.residuals_field`.

        Args:
            measurement (ECOpticalMeasurement): The measurement. Defaults to
                self.measurement.
        Return list of ValueSeries: The coefficient of each converged spectrum
        """
        from scipy.optimize import lsq_linear

        measurement = measurement or self.measurement
        field = measurement[self.field_name]
        tseries, xseries = field.axes_series

        wl_mask = get_x_mask(self.wavelength, self.xspan)
        wl = self.wavelength[wl_mask]
        converged = np.ma.filled(self.converged_data[:, wl_mask], 0)  # NaN breaks lsq
        data = np.ma.filled(np.ma.asarray(field.data, dtype=float), 0)
        # Put the measured spectra on the wavelengths of the converged spectra
        spectra = interp1d(xseries.data, data, axis=1)(wl)

        coefficients = []
        for i, spectrum in enumerate(spectra):
            if self.bounds_bool:
                # The lower bound of each component is its previous value, minus
                # the noise
                lower = coefficients[i - 1] - self.noise_bound / 1000 if i else 0
                result = lsq_linear(converged.T, spectrum, bounds=(lower, np.inf))
            else:
                result = lsq_linear(converged.T, spectrum)
            coefficients.append(np.asarray(result.x))
        coefficients = np.array(coefficients)

        wl_series = DataSeries(name="wavelength", unit_name="nm", data=wl)
        reconstructed = coefficients @ converged
        self.coefficients = [
            ValueSeries(
                name=f"coefficient_{i}",
                unit_name="",
                data=c,
                tseries=tseries,
            )
            for i, c in enumerate(coefficients.T)
        ]
        self._fitted_measurement = measurement
        self.reconstructed_field = Field(
            name="reconstructed",
            unit_name=field.unit_name,
            data=reconstructed,
            axes_series=[tseries, wl_series],
        )
        self.residuals_field = Field(
            name="residuals",
            unit_name="mOD",
            data=(spectra - reconstructed) * 1000,
            axes_series=[tseries, wl_series],
        )
        return self.coefficients
