"""Representation and analysis of thermal catalysis (TP) with MS measurements"""
from .ms import MSMeasurement, MSSpectroMeasurement
from ..measurements import Calibration
from ..plotters.tpms_plotter import TPMSPlotter, SpectroTPMSPlotter
from ..data_series import ValueSeries
import warnings
import numpy as np


class ReactorMeasurement(MSMeasurement):
    """Class implementing thermal catalysis measurement.

    Now the class has implemented temperature and pressure, and mass spectrometry
    through the MSMeasurement. Future updates should include GC measurements.

    The main job of this class is making sure that the ValueSeries most essential for
    visualizing thermal catalysis measurements are always available in the
    correct form as the measurement is added with others, reduced to a selection,
    calibrated and normalized, etc. These most important ValueSeries are:

    - `temperature`: The temperature inside the reactor typically in [K].

    - `pressure`: The pressure inside the reactor typically in [bar].

    These ValueSeries are directly accesable as properties of the class

    - `measurement.T` yields the data series for temperature
    - `measurement.P` yields the data series for pressure
    - `measurement.t` yield the TimeSeries for the measurement

    The inverse of the temperature series is often used in themal catalysis when fitting
    an arrhenius equation to the data.
    The inverse of a `ValueSeries` is accesable by adding "inverse_" in front of the
    ValueSeries.name (i.e measurement["inverse_temperature"]).

    For temperature this has can be accessed directly as a property

    - `measurement.inverse_T`

    `ReactorMeasurement` inherits from `MSMeasurement` to acces and calibrate MS data

    `ReactorMeasurement` has a default plotter `TPMSPlotter` which can plot a two panel
    plot with MS data in one panel and temperature and pressure in the other panel
    against time (`measurement.plot_measurement()`) or
    plot in one panel with MS data on left y-axis and temperature on right y-axis against
    time (`measurement.plotter.plot_measurement_in_one_panel()`) or
    plot arrhenius like plot with MS data plotted against inverse of temperature and fit
    arrhenius equation in specific tspan_list (`measurement.plotter.plot_arrhenius()`)

    """

    default_plotter = TPMSPlotter
    essential_series_names = ("temperature", "pressure")

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
        self.add_calibration(ReactorCalibration(name="inverse_calibration"))
        self.activation_energy = {}

    @property
    def T_name(self):
        return self["temperature"].name

    @property
    def P_name(self):
        return self["pressure"].name

    @property
    def t_name(self):
        return self["temperature"].tseries.name

    @property
    def T(self):
        return self["temperature"].data

    @property
    def inverse_T(self):
        return self["inverse_temperature"].data

    @property
    def P(self):
        return self["pressure"].data

    @property
    def t(self):
        return self["temperature"].t

    @property
    def meta_list(self):
        """List of the ValueSeries contained in the measurement that is not a mass"""
        return [
            col
            for col in self.series_names
            if col not in self.mass_list and col[-2:] != "-x"
        ]

    def unit_converter(self, v_name, new_unit):
        """Convert dataseries from one unit to another using self.correct_data from super
        in the case of temperature the unit_factor should be added not subtracted.
        Typical usage:
            measurement.unit_converter('temperature', 'K')

        Args:
            v_name (str): name of a DataSeries to convert unit from
            new_unit (str): Name of the new unit to convert to
        Return:
           DataSeries

        """
        self.clear_cache()  # To achieve correct behavior with calibrated series
        original_unit = self[v_name].unit.name
        if not original_unit:
            warnings.warn(
                f"Series {v_name} does not contain a unit."
                "Set unit before using this function",
                stacklevel=2,
            )
            return
        unit_factor, unit_name = self._get_unit_factor(  # self,
            original_unit=original_unit, new_unit=new_unit
        )

        if original_unit == unit_name:
            warnings.warn(f"Did not convert {original_unit} to {new_unit}", stacklevel=2)
            return
        data = self[v_name].data

        if (
            "temperatur" in v_name.lower()
        ):  # converting from kelvin to celsius is addition
            new_data = data + unit_factor
        else:
            new_data = data * unit_factor

        self.correct_data(v_name, new_data)

        self[v_name].unit.name = new_unit

    def fit_to_arrhenius_equation(self, inverse_T, k, logdata=False):
        """Method to fit data to an arrhenius expression
        k = A exp(Ea/R T) # per mole
        k = A exp(Ea/kB T) # per molecule
        kB_J = 1.380649e-23 J K-1, or kB_eV =  8.617333e-5 eV K-1

        args:
            inverse_T (list): list of inverse temperatures
            k (list): list of rates at T (often mol/s)
            logdata (bool): Used if k is np.log(k). Default False.
        return (either or):
            A, Ea (if logdata)
            popt, popv (if not logdata)
        """
        if len(inverse_T) != len(k):
            warnings.warn("Length of T and k has to be equal")
            return

        # if data is ln(k) the regression is linear
        coef = np.polyfit(inverse_T, k, 1)
        kB = 8.617333e-5
        R = 8.31446261815324
        Ea = -R * coef[0]
        A = np.exp(coef[1])

        if logdata:
            print(
                f"pre exponential factor A = {A},\n and activity energy Ea = {Ea}\n"
                "universal gas constant R = 8.314.. J mol-1 K-1\n"
            )
            print(
                f"activity energy Ea/R = {coef[0]} => Ea[J K-1] = {-R * coef[0]}\n"
                f"activity energy Ea/kB = {coef[0]} => Ea [eV K-1] = {-kB * coef[0]}\n"
            )
            return coef  # A, Ea
        from scipy.optimize import curve_fit

        popt, pcov = curve_fit(self._func, inverse_T, k)
        print(
            f"pre factor A = {popt[0]}, Ea/R = {popt[1]}, Ea = {popt[1]*R} [J/K mol]"
            f"activity energy Ea/kB = {popt[1]} => Ea [eV K-1] = {-kB * popt[1]}\n"
        )  # kB = 8.617e-5 eV

        return popt, pcov

    def _func(self, x, a, b):
        """helper function for fitting arrhenius equation to quantified mass data"""
        return a * np.exp(-b * x)

    def _get_unit_factor(
        self,
        original_unit,
        new_unit,
    ):
        """Return a unit_factor to convert from one unit to another unit.
        helper function to convert different units used in a TPMS measurement
        Args:
            original_unit (str): The original unit of a DataSeries
            new_unit (str): The new unit for DataSerie
        Return:
            unit_factor, new_unit_name

        """

        try:
            if "1/" in (original_unit or new_unit):
                warnings.warn(
                    f"Recieved an inverted serie with unit '{original_unit}'."
                    f"Cannot convert units of inverted series'",
                    stacklevel=2,
                )
            if original_unit == "celsius" or original_unit == "C":
                unit_factor = {
                    "C": 0,
                    "celsius": 0,
                    "K": 273.15,
                    "kelvin": 273.15,
                }[new_unit]

            elif original_unit == "kelvin" or original_unit == "K":
                unit_factor = {
                    "K": 0,
                    "kelvin": 0,
                    "celsius": -273.15,
                    "C": -273.15,
                }[new_unit]

            elif original_unit == "mbar":
                unit_factor = {
                    "mbar": 1,
                    "bar": 1000,
                    "hPa": 1,
                    "kPa": 0.1,
                }[new_unit]

            elif original_unit == "bar":
                unit_factor = {
                    "mbar": 1e-3,
                    "bar": 1,
                    "hPa": 1e-3,
                    "kPa": 0.1e-3,
                }[new_unit]

            else:
                unit_factor = {
                    # Time conversion
                    "s": 1,
                    "min": 1 / 60,
                    "minutes": 1 / 60,
                    "h": 1 / 3600,
                    "hr": 1 / 3600,
                    "hour": 1 / 3600,
                    "hours": 1 / 3600,
                    "d": 1 / (3600 * 24),
                    "days": 1 / (3600 * 24),
                    # Pressure conversion
                    "mbar": 1,
                    "bar": 1000,
                    "hPa": 1,
                    "kPa": 0.1,
                    # Temperature conversion
                    "K": 273.15,
                    "kelvin": 273.15,
                }[new_unit]
        except KeyError:
            warnings.warn(
                f"Can't convert original unit '{original_unit}' to new unit"
                f"'{new_unit}'. You can try to clear cache with "
                "measurement.clear_cache(). Using original unit {original_unit} !",
                stacklevel=2,
            )
            unit_factor = 1
            new_unit = original_unit
        return unit_factor, new_unit


class ReactorSpectroMeasurement(ReactorMeasurement, MSSpectroMeasurement):
    default_plotter = SpectroTPMSPlotter


class ReactorCalibration(Calibration):
    """A reactor calibration to calibrate inverse of meta_series"""

    def __repr__(self):
        # TODO: make __repr__'s consistent throught ixdat
        return (
            f"{self.__class__.__name__}"
            f"(Calibration={self.name} for setup={self.setup} on date={self.date})"
        )

    def calibrate_series(self, key, measurement=None):
        """Return a calibrated series for key based on the raw data in the measurement.

        Key should start with "inverse". Anything else will return None.

        - "temperature": the calibration looks up "temperature" in the measurement,
        translate it to 1/T and then returns a calibrated temperature
        """

        measurement = measurement or self.measurement
        if key.startswith("inverse_"):
            meta = key.split("_")[-1]
            signal_series = measurement[meta]
            unit_name = measurement[meta].unit.name
            y = signal_series.data
            y_inverse = 1 / y
            return ValueSeries(
                name=key,
                unit_name=f"1/{unit_name}",
                data=y_inverse,
                tseries=signal_series.tseries,
            )
        elif key.startswith("log_") or key.startswith("ln_"):
            meta = key.split("_")[-1]
            signal_series = measurement[meta]
            unit_name = measurement[meta].unit.name
            y = signal_series.data
            log_y = np.log(y)
            return ValueSeries(
                name=key,
                unit_name=f"ln({unit_name})",
                data=log_y,
                tseries=signal_series.tseries,
            )
