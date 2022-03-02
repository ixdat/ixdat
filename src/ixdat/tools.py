"""This module contains general purpose tools"""
import inspect
import warnings
from functools import wraps

import numpy as np

from ixdat.exceptions import DeprecationError

warnings.simplefilter("default")


def thing_is_close(thing_one, thing_two):
    """Return whether two things are (nearly) equal, looking recursively if necessary"""
    if type(thing_one) != type(thing_two):
        return False

    if isinstance(thing_one, list):
        return list_is_close(thing_one, thing_two)
    elif isinstance(thing_one, dict):
        return dict_is_close(thing_one, thing_two)
    else:
        return value_is_close(thing_one, thing_two)


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


def _construct_deprecation_message(
    identity,
    last_release_when_supported,
    soft_deprecation_period,
    hard_deprecation_period,
    message,
    kwarg_name,
):
    """Return a deprecation message, see :func:`deprecate` for details on the arguments

    Args:
        identity (str): The identity of the callable being deprecated e.g.
            "fuction 'myfunction'"

    All other arguments are as in :func:`deprecate`

    """
    property_part = f"property named '{kwarg_name}' in " if kwarg_name else ""
    return (
        f"The {property_part}{identity} is deprecated, its last supported version being "
        f"{last_release_when_supported}.\n"
        "This means that:\n"
        "* It is soft deprecated (issuing warnings) in all releases for "
        f"{soft_deprecation_period} after the last supported release\n"
        "* It is hard deprecated (raising exceptions) in all releases for an additional "
        f"{hard_deprecation_period} after that\n"
        "* After the hard deprecation period ends, it will be removed\n\n"
        "See instruction below on how to update your code to avoid this message:\n"
        f"{message}"
    )


def deprecate(
    last_release_when_supported,
    soft_deprecation_period,
    hard_deprecation_period,
    message,
    kwarg_name=None,
    hard_deprecate=False,
):
    """Mark a function, method or class for deprecation

    The deprecator supports soft and hard deprecation which will either issue warnings or
    raise exceptions. The deprecation periods are calculated from the date of
    `last_release_when_supported`.

    Args:
        last_release_when_supported (str): The name of the last version when the deprecated
            functionality was fully supported e.g. "1.2.3"
        soft_deprecation_period (str): The period after `last_release_when_supported` when
            the functionality will be soft deprecated, meaning using warnings e.g. "3 months"
        hard_deprecation_period (str): The period after `last_release_when_supported` +
            `soft_deprecation_period` when the functionality will raise Exception, but still
            be present. E.g. "1 month". After the period the functionality will be removed
            entirely.
        message (str): A message to the user instructing with instruction on how to upgrade
            to avoid the deprecated functionality

    Keyword Args:
        kwarg_name (str): If given, is the name of the keyword argument which is deprecated.
            If not given it is assumed that the entire callable is deprecated.
        hard_deprecate (bool): Whether the functionality is now hard deprecated. Defaults to
            False.

    Examples:

    Used to deprecate a class::
     @deprecate("1.2.3", "1 month", "1 month", "Please use `MyNewClass` instead")
     class MyClass:
         ...

    Used to deprecate a method or an argument in a method::
        class MyClass:

            @deprecate("1.2.3", "1 month", "1 month", "Please use `mynewmethod` instead")
            def mymethod(self):
                ...

            @classmethod
            @deprecate(
                "1.2.3",
                "1 month",
                "1 month",
                "Please use `kwargb` instead",
                kwarg_name="kwarga"
            )
            def mymethod(cls, arg, kwarga=None, kwargb=None):
                ...

    Note: In the example above, when this decorator is applied to class or static methods,
    it must be applied as the first decorator (closest to the def).

    """

    def decorator(callable_):
        """Decorate a callable with"""
        # Form an identity string for the object which is being decorated, which is used
        # in the message to the user
        if inspect.isclass(callable_):
            identity = f"class '{callable_.__qualname__}'"
        else:
            if "." in callable_.__qualname__:
                identity = f"method '{callable_.__qualname__}'"
            else:
                identity = f"function '{callable_.__qualname__}'"

        compound_message = _construct_deprecation_message(
            identity,
            last_release_when_supported,
            soft_deprecation_period,
            hard_deprecation_period,
            message,
            kwarg_name,
        )

        # Get the argument signature of the callable
        callable_signature = inspect.signature(callable_)

        @wraps(callable_)
        def inner_function(*args, **kwargs):
            # Form bound arguments, so both args and kwargs gets mapped to their name
            bound_arguments = {}
            if kwarg_name:
                bound_arguments = callable_signature.bind(*args, **kwargs).arguments

            # If something deprecated is being used, issue the warning or exception
            if kwarg_name is None or kwarg_name in bound_arguments:
                if hard_deprecate:
                    raise DeprecationError(compound_message)
                else:
                    warnings.warn(compound_message, DeprecationWarning, stacklevel=2)

            # Calculate the return value of the original object and return
            return_value = callable_(*args, **kwargs)

            return return_value

        return inner_function

    return decorator
