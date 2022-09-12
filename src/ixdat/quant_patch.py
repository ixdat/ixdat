"""This module has been renamed si_quant_patch. """
from .si_quant_patch import (
    append_sensitivity_factors as real_append_sensitivity_factors,
    append_two_sensitivity_factors as real_append_two_sensitivity_factors,
    add_isotopes as real_add_isotopes,
)
from .tools import deprecate


@deprecate(
    last_supported_release="0.2.4",
    update_message="The module name has changed to `si_quant_patch`",
    remove_release="0.3",
)
def append_sensitivity_factors(*args):
    return real_append_sensitivity_factors(*args)


@deprecate(
    last_supported_release="0.2.4",
    update_message="The module name has changed to `si_quant_patch`",
    remove_release="0.3",
)
def append_two_sensitivity_factors(sf_1, sf_2):
    return real_append_two_sensitivity_factors(sf_1, sf_2)


@deprecate(
    last_supported_release="0.2.4",
    update_message="The module name has changed to `si_quant_patch`",
    remove_release="0.3",
)
def add_isotopes(calibration, isotope_spec):
    return real_add_isotopes(calibration, isotope_spec)
