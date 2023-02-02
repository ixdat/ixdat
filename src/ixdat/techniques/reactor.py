from .ms import MSMeasurement, SpectroMSMeasurement
from ..measurements import Calibration
from ..plotters.tpms_plotter import TPMSPlotter, SpectroTPMSPlotter
from ..db import Saveable
from ..data_series import ValueSeries
import warnings
import numpy as np

class ReactorMeasurement(MSMeasurement):

    default_plotter = TPMSPlotter
    essential_series_names = ("temperature", "pressure")

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
    def T_inverse(self):
        return self["inverse_temperature"].data

    @property
    def P(self):
        return self["pressure"].data

    @property
    def t(self):
        return self["temperature"].t

    @property
    def meta_list(self):
        """List of the ValueSeries contained in the measurement that is not a mass """
        return [col for col in self.series_names if col not in self.mass_list and col[-2:] != '-x']

    def unit_converter(self, v_name, new_unit):
        """Convert dataseries from one unit to another using self.correct_data from super

        Typical usage:
            measurement.unit_converter('temperature', 'K')

        Args:
            v_name (str): name of a DataSeries to convert unit from
            new_unit (str): Name of the nw unit to convert to
        Return

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

        if "temperatur" in v_name.lower():  # converting from kelvin to celsius is addition
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
        print(f"pre exponential factor A = {A},'\n', and activity energy Ea = {Ea}"
              "universal gas constant R = 8.314.. J mol-1 K-1"
              )
        print(f"activity energy Ea/R = {coef[0]} => Ea[J K-1] = {-R * coef[0]}"
              f"activity energy Ea/kB = {coef[0]} => Ea [eV K-1] = {-kB * coef[0]}"
              )

        if logdata:
            return coef#A, Ea
        from scipy.optimize import curve_fit

        popt, pcov = curve_fit(self._func, inverse_T, k)
        print(f"pre exponential factor A = {popt[0]}, Ea/R = {popt[1]}, Ea = {popt[1]*R}"
               )# kB = 8.617e-5 eV

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
                warning.warn(
                        f"Recieved an inverted serie with unit '{original_unit}'."
                        f"Cannot convert units of inverted series'",
                        stacklevel=2
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
                    "K": 1,
                    "kelvin": 1,
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


class SpectroReactorMeasurement(ReactorMeasurement, SpectroMSMeasurement):
    default_plotter = SpectroTPMSPlotter


class ReactorUnitConverterResult(Saveable):
    """A class for a reactor reactor_calibration result.

    FIXME: I think that something inheriting directly from Saveable does not belong in
        a technique module.
    """

    table_name = "reactor_cal_results"
    column_attrs = {"name", "mass", "new_unit"}

    def __init__(
        self,
        name=None,
        mass=None,
        new_unit=None,
    ):
        super().__init__()
        self.name = name or f"{mass}@{new_unit}"
        self.mass = mass
        self.new_unit = new_unit

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name={self.name}, mol={self.mass}, "
            f"calibration type={self.cal_type})"
        )

    # @property
    # def color(self):
    #    return STANDARD_COLORS[self.mass]


class ReactorCalibration(Calibration):
    """A reactor calibration to calibrate inverse of meta_series"""

    extra_column_attrs = {
        "reactor_calibration_results": ("reactor_cal_results", "reactor_cal_result_id")
    }
    # TODO: https://github.com/ixdat/ixdat/pull/11#discussion_r677552828

    def __init__(
        self,
        name=None,
        date=None,
        tstamp=None,
        setup=None,
        reactor_cal_results=None,
        technique="reactor",
        new_unit=None,
        measurement=None,
    ):
        """Initiate a Calibration

        Args:
            name (str): The name of the reactor_calibration
            date (str): Date of the reactor_calibration
            setup (str): Setup where the measurmenet was done
            technique (str): The technique of the calibration
            tstamp (float): The time at which the calibration took place or is valid
            measurement (TPMSMeasurement): Optional.A measurement to calibrate by default
        """
        super().__init__(
            name=name or f"TPMS reactor_calibration for {setup} on {date}",
            technique=technique,
            tstamp=tstamp,
            measurement=measurement,
        )

        self.name = name
        self.date = date
        self.setup = setup
        self.new_unit = new_unit

        # self.reactor_unit_convert = reactor_unit_results or []

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
            y_inverse = 1 / y
            return ValueSeries(
                name=key,
                unit_name=f"ln({unit_name})",
                data=y_inverse,
                tseries=signal_series.tseries,
            )
