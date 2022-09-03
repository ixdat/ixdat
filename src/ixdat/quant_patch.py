"""Patches corresponding to code that should be added to Spectro Inlets quantification"""


def append_sensitivity_factors(*args):
    """Recursively append elements in a list"""
    sf = append_two_sensitivity_factors(args[0], args[1])
    if len(args) == 2:
        return sf
    remaining_args = [sf] + list(args[2:])
    return append_sensitivity_factors(*remaining_args)


def append_two_sensitivity_factors(sf_1, sf_2):
    """This should be incorporated into __add__ of SensitivityFactor and SensitivityList

    Appends all the sensitivity factors in sf_1 and sf_2, which can each be a single
    SensitivityFactor (or CalPoint) or a SensitivityList (or Calibration).

    Args:
        sf_1 (SensitivityFactor or SensitivityList): The first sensitivity factor or list
            thereof.
        sf_2 (CalPoint or Calibration): The second sensitivity factor or list thereof
    """
    from spectro_inlets_quantification.sensitivity import (
        SensitivityList, SensitivityFactor
    )
    from spectro_inlets_quantification.calibration import Calibration

    if isinstance(sf_1, SensitivityList):
        if isinstance(sf_2, SensitivityList):
            return sf_1 + sf_2   # this is implemented in spectro_inlets_quantification
        elif isinstance(sf_2, SensitivityFactor):
            sf_list = sf_1.sf_list + [sf_2]
            obj_as_dict = sf_1.as_dict()
            if isinstance(sf_1, Calibration):
                sf_1.pop("cal_dicts")
                obj_as_dict["cal_list"] = sf_list
                return Calibration(**obj_as_dict)
            else:
                sf_1.pop("sf_dicts")
                obj_as_dict["sf_list"] = sf_list
                return SensitivityList(**obj_as_dict)
        else:
            raise TypeError(f"Can't add {sf_1} and {sf_2}")
    elif isinstance(sf_1, SensitivityFactor):
        if isinstance(sf_2, SensitivityFactor):
            sf_list = [sf_1, sf_2]
            return Calibration(cal_list=sf_list)
        else:
            return append_sensitivity_factors(sf_2, sf_1)
    else:
        raise TypeError(f"Can't add {sf_1} and {sf_2}")


def add_isotopes(calibration, isotope_spec):
    pass
