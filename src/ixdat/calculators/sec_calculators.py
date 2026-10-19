"""Calculators for spectroelectrochemistry (SEC) measurements"""

import warnings

import numpy as np
from scipy.interpolate import interp1d

from ..measurement_base import Calculator
from ..data_series import Field, TimeSeries
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
        self.normalise = normalise
        self.xspan = xspan
        self.idx_step = idx_step
        self.tstep = tstep
        self.Vstep = Vstep
        self.interpolate = interpolate

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
