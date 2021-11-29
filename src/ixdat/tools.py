"""This module contains general purpose tools"""

import numpy as np


def thing_is_close(thing_one, thing_two):
    """Return whether two things are (nearly) equal, looking recursively if necessary"""
    if type(thing_one) != type(thing_two):
        return False

    if isinstance(thing_one, list):
        if not list_is_close(thing_one, thing_two):
            return False
    elif isinstance(thing_one, dict):
        if not dict_is_close(thing_one, thing_two):
            return False
    else:
        if not value_is_close(thing_one, thing_two):
            return False

    return True


def value_is_close(value_one, value_two):
    """Return whether `value_one` and `value_two` are equal (or close for floats)"""
    if isinstance(value_one, float) or isinstance(value_two, float):
        return np.isclose(value_one, value_two)
    elif isinstance(value_one, np.ndarray) and isinstance(value_two, np.ndarray):
        return np.allclose(value_one, value_two)

    return value_one == value_two


def dict_is_close(dict_one, dict_two):
    """Return True if the dicts are equal, except floats are allowed to just be close

    Args:
        dict_one (dict): The first dictionary to compare
        dict_two (dict): The second dictionary to compare

    .. warning::
       This function is recursive (also together with :ref:func:`list_is_close` and
       **will** result in an infinite loop if the values reference back to itself
    """
    if len(dict_one) != len(dict_two):
        return False
    if dict_one.keys() != dict_two.keys():
        return False

    for key, value_one in dict_one.items():
        value_two = dict_two[key]
        if not thing_is_close(value_one, value_two):
            return False

    return True


def list_is_close(list_one, list_two):
    """Return True if the lists are equal, except floats are allowed to just be close

    Args:
        list_one (list): The first list to compare
        list_two (list): The second list to compare

    .. warning::
       This function is recursive (also together with :ref:func:`list_is_close` and
       **will** result in an infinite loop if the values reference back to itself
    """
    if len(list_one) != len(list_two):
        return False

    for value_one, value_two in zip(list_one, list_two):
        if type(value_one) != type(value_two):
            return False
        if not thing_is_close(value_one, value_two):
            return False

    return True
