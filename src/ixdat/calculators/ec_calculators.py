# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:29:09 2024

@author: Søren
"""
from ..measurement_base import Calculator
from ..data_series import ValueSeries
from ..plotters.ec_plotter import EC_FANCY_NAMES
from ..exceptions import QuantificationError
from .scan_rate_tools import calc_sharp_v_scan


class ECCalibration(Calculator):
    """An electrochemical calibration with RE_vs_RHE, A_el, and/or R_Ohm"""

    calculator_type = "EC calibration"
    extra_column_attrs = {"ec_calibration": {"RE_vs_RHE", "A_el", "R_Ohm"}}
    # TODO: https://github.com/ixdat/ixdat/pull/11#discussion_r677552828

    def __init__(
        self,
        name=None,
        technique="EC",
        tstamp=None,
        measurement=None,
        RE_vs_RHE=None,
        A_el=None,
        R_Ohm=None,
    ):
        """Initiate a Calibration

        Args:
            name (str): The name of the calibration
            technique (str): The technique of the calibration
            tstamp (float): The time at which the calibration took place or is valid
            measurement (ECMeasurement): Optional. A measurement to calibrate by default
            RE_vs_RHE (float): The reference electrode potential on the RHE scale in [V]
            A_el (float): The electrode area in [cm^2]
            R_Ohm (float): The ohmic drop resistance in [Ohm]
        """
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.RE_vs_RHE = RE_vs_RHE
        self.A_el = A_el
        self.R_Ohm = R_Ohm

    def __repr__(self):
        # TODO: make __repr__'s consistent throught ixdat
        return (
            f"{self.__class__.__name__}"
            f"(RE_vs_RHE={self.RE_vs_RHE}, A_el={self.A_el}, R_Ohm={self.R_Ohm})"
        )

    available_series_names = {"potential", "current"}
    # This is a constant class attribute rather than a property because an
    # ECCalibration always returns potential and current, as correctly as possible.

    def calculate_series(self, key, measurement=None):
        """Return a calibrated series for key based on the raw data in the measurement.

        Key should be "potential" or "current". Anything else will return None.

        - "potential": the calibration looks up "raw_potential" in the measurement,
        shifts it to the RHE potential if RE_vs_RHE is available, corrects it for
        Ohmic drop if R_Ohm is available, and then returns a calibrated potential
        series with a name indicative of the corrections done.
        - "current": The calibration looks up "raw_current" in the measurement,
        normalizes it to the electrode area if A_el is available, and returns a
        calibrated current series with a name indicative of whether the normalization
        was done.
        """
        measurement = measurement or self.measurement
        if key == "potential":
            raw_potential = measurement["raw_potential"]
            name = raw_potential.name
            U = raw_potential.data
            if self.RE_vs_RHE is not None:
                U = U + self.RE_vs_RHE
                name = EC_FANCY_NAMES["potential"]
            if self.R_Ohm is not None:
                I_mA = measurement.grab_for_t("raw_current", t=raw_potential.t)
                U = U - self.R_Ohm * I_mA * 1e-3  # [V] = [Ohm*mA*(A/mA)]
                name = name + " $_{ohm. corr.}$"
            return ValueSeries(
                name=name,
                unit_name=raw_potential.unit_name,
                data=U,
                tseries=raw_potential.tseries,
            )

        if key == "current":
            raw_current = measurement["raw_current"]
            name = raw_current.name
            J = raw_current.data
            if self.A_el is not None:
                J = J / self.A_el
                name = EC_FANCY_NAMES["current"]
            return ValueSeries(
                name=name,
                unit_name=raw_current.unit_name,
                data=J,
                tseries=raw_current.tseries,
            )

        raise QuantificationError(f"{self} can not calculate {key}")

    def __add__(self, other):
        """Add two ECCalibrations. The second one (`other`) takes priority.

        In other words, if both `ec_cal_1` and `ec_cal_2` have their own value of
        `RE_vs_RHE`, then `(ec_cal_1 + ec_cal_2).RE_vs_RHE` will be the same as
        `ec_cal_2.RE_vs_RHE`. This ensures that the values in the`ECCalibration` most
        recently added to a measurement (e.g. to update a value) is used over older ones.
        """
        RE_vs_RHE = other.RE_vs_RHE if other.RE_vs_RHE is not None else self.RE_vs_RHE
        A_el = other.A_el if other.A_el is not None else self.A_el
        R_Ohm = other.R_Ohm if other.R_Ohm is not None else self.R_Ohm
        name = self.name + " + " + other.name
        return ECCalibration(
            RE_vs_RHE=RE_vs_RHE,
            A_el=A_el,
            R_Ohm=R_Ohm,
            name=name,
        )


class ScanRateCalculator:
    calculator_type = "scan_rate_calculator"
    available_series_names = {"scan_rate"}

    def calculate_series(self, key, measurement, res_points=10):
        """The scan rate as a ValueSeries"""
        if not key == "scan_rate":
            raise QuantificationError(f"{self} can not calculate {key}")
        t, v = measurement.grab("potential")
        scan_rate_vec = calc_sharp_v_scan(t, v, res_points=res_points)
        scan_rate_series = ValueSeries(
            name="scan_rate",
            unit_name="V/s",  # TODO: unit = potential.unit / potential.tseries.unit
            data=scan_rate_vec,
            tseries=measurement.potential.tseries,
        )
        return scan_rate_series
