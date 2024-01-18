"""This module contains general purpose tools"""
import datetime
import inspect
import time
import warnings
from functools import wraps
from string import ascii_uppercase
from typing import Optional

import numpy as np
from packaging import version

from ixdat.exceptions import DeprecationError
import ixdat.config

# from ixdat.config import CFG
from ixdat import __version__

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
    callable_,
    last_supported_release,
    update_message,
    hard_deprecation_release,
    remove_release,
    kwarg_name,
):
    """Return a deprecation message

    Args:
        callable_ (Callable): The callable to form a deprecation message for

    All other arguments are as in :func:`deprecate`

    """
    # Form an identity string for the object which is being deprecated, which is used
    # in the message to the user
    identity = f"argument named '{kwarg_name}' in " if kwarg_name else ""
    if inspect.isclass(callable_):
        identity += f"class '{callable_.__qualname__}'"
    else:
        if "." in callable_.__qualname__:
            identity += f"method '{callable_.__qualname__}'"
        else:
            identity += f"function '{callable_.__qualname__}'"

    message = (
        f"The {identity} is deprecated, its last supported version being "
        f"{last_supported_release}:\n"
    )
    # Add information on potential hard deprecation
    if hard_deprecation_release is not None:
        message += (
            "* It will become hard deprecated, raising exceptions rather than "
            f"issuing warnings, from version {hard_deprecation_release}\n"
        )
    else:
        message += (
            "* It will continue to be soft deprecated, issuing warnings, for the "
            "foreseeable future\n"
        )

    # Add information of potential removal
    if remove_release is not None:
        message += f"* It is planned for complete removal in version {remove_release}\n"

    message += (
        "\n"
        "See instructions below on how to update your code to avoid this message:\n"
        f"{update_message}"
    )
    return message


def deprecate(
    last_supported_release,
    update_message,
    hard_deprecation_release=None,
    remove_release=None,
    kwarg_name=None,
):
    """Mark a function, method, programmed property or class for deprecation

    The deprecator supports soft and hard deprecation, which will either issue warnings
    or raise exceptions, as well as providing information about an update path and the
    potential time the functionality will be completely removed.

    Args:
        last_supported_release (str): The name of the last version when the deprecated
            functionality was fully supported e.g. "1.2.3"
        update_message (str): A message to the user with instructions on how to upgrade
            to avoid the deprecated functionality

    Keyword Args:
        hard_deprecation_release (str): The release with which the deprecation will raise
            exceptions rather than issue warnings e.g. "1.4.0". The default None means
            that this deprecation will remain soft for the foreseeable future.
        remove_release (str): The release with which the deprecated functionality is
            planned for deletion
        kwarg_name (str): If given, is the name of the keyword argument which is
            deprecated. If not given it is assumed that the entire callable is
            deprecated.

    Examples:

    Used to deprecate a class::
     @deprecate("1.2.3", "Please use `MyNewClass` instead", "1.4.0", "2.0.0")
     class MyClass:
         ...

    Used to deprecate a method or an argument in a method::
        class MyClass:

            @deprecate("1.2.3", "Please use `mynewmethod` instead", "1.4.0", "2.0.0")
            def mymethod(self):
                ...

            @classmethod
            @deprecate(
                "1.2.3",
                "Please use `kwargb` instead",
                "1.4.0",
                "2.0.0",
                kwarg_name="kwarga"
            )
            def mymethod(cls, arg, kwarga=None, kwargb=None):
                ...

    Used to decorate a property::
        class MyClass:
            def __init__(self):
                self._internal = 8

            @property
            @deprecate("1.2.3", "Please use `new_external` instead")
            def external(self):
                return self._internal

            @external.setter
            @deprecate("1.2.3", "Please use `new_external` instead")
            def external(self, value):
                self._internal = value

    .. note::
       In the examples above, when this decorator is applied to a class method, static
       method or programmed property, it must be applied as the first decorator (closest
       to the def).

    """

    def decorator(callable_):
        """Decorate a callable with"""

        compound_message = _construct_deprecation_message(
            callable_,
            last_supported_release,
            update_message,
            hard_deprecation_release,
            remove_release,
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
                if hard_deprecation_release is not None and version.parse(
                    __version__
                ) >= version.parse(hard_deprecation_release):
                    raise DeprecationError(compound_message)
                else:
                    # The stacklevel argument here is used to make the warning reference
                    # the line at which the decorated callable was called, as the place
                    # where the warning is raised, instead of here in `innner_function`,
                    # which would be useless
                    warnings.warn(compound_message, DeprecationWarning, stacklevel=2)

            # Calculate the return value of the original object and return
            return_value = callable_(*args, **kwargs)

            return return_value

        return inner_function

    return decorator


@deprecate(
    last_supported_release="0.2.8",
    update_message=(
        "Please use `ixdat.tools.tstamp_to_string` with the keyword argument "
        "`string_format='native_date'` instead."
    ),
    hard_deprecation_release="0.2.12",
    remove_release="0.2.14",
)
def tstamp_to_yyMdd(tstamp: float) -> str:
    """Return the date in compact form "yyMdd" format given the unix time (float).
    In this format the month is given as a capital letter, starting with A for January.
    E.g. June 4th, 2022 will become 22F04.
    """
    a = time.localtime(tstamp)
    year = a.tm_year
    month = a.tm_mon
    day = a.tm_mday
    date_string = "{0:02d}{1:1s}{2:02d}".format(
        year % 100, chr(ord("A") + month - 1), day
    )
    return date_string


def tstamp_to_string(tstamp: float, string_format: Optional[str] = None) -> str:
    """Return a string representation of unix timestamps `tstamp`

    Args:
        tstamp (float): The unix time stamp to convert
        string_format (str): Optional. The datetime string format to use. If not given,
            the value of ``ixdat.config.config.timestamp_string_format`` will be used.
            Accepted values are format strings accepted by `datetime.datetime.strftime`
            or "native" or "native_date", which will produce ixdat native datetime
            strings or date strings respectively.

    """
    dt = datetime.datetime.fromtimestamp(tstamp, ixdat.config.config.timezone)
    if string_format is None:
        string_format = ixdat.config.config.timestamp_string_format

    if string_format in ("native", "native_date"):
        # ixdat shows months as capital letters, where Jan->A, Feb->B etc.
        month_letter = ascii_uppercase[dt.month - 1]
        if string_format == "native":
            # Brings to the total format to: 22E18 14:34:55
            string_format = f"%y{month_letter}%d %H:%M:%S"
        else:
            string_format = f"%y{month_letter}%d"

    return dt.strftime(string_format)


if __name__ == "__main__":
    t0 = time.time()
    print(tstamp_to_string(t0))
