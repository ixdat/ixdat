"""Tests for tools.py"""
import warnings
from unittest.mock import patch

import pytest
from packaging import version

from ixdat.exceptions import DeprecationError
from ixdat.tools import deprecate

# Standard arguments for deprecate
DEPRECATE_STANDARD_ARGS = {
    "last_supported_release": "1.2.3",
    "update_message": "Do SOMETHING ELSE now.",
}


def generate_function_and_classes(decorator, *args):
    """Generate decorated function and classes"""

    class MyClass:
        def __init__(self, intarg, intkwarg=2, one_more_kwarg=None):
            pass

        @decorator
        def mymethod(self, intarg, intkwarg=2, one_more_kwarg=None):
            return 47 + intarg + intkwarg

        @classmethod
        @decorator
        def myclassmethod(cls, intarg, intkwarg=2, one_more_kwarg=None):
            return 47 + intarg + intkwarg

        @staticmethod
        @decorator
        def mystaticmethod(intarg, intkwarg=2, one_more_kwarg=None):
            return 47 + intarg + intkwarg

    DecoratedMyClass = decorator(MyClass)

    @decorator
    def myfunction(intarg, intkwarg=2, one_more_kwarg=None):
        return 47 + intarg + intkwarg

    return MyClass, DecoratedMyClass, myfunction


class TestDeprecate:
    """Test the `deprecate` decorator"""

    @pytest.mark.parametrize("decorate_kwarg", [None, "intkwarg"])
    @pytest.mark.parametrize("current_release", ["1.3.0", "1.9.0"])
    @pytest.mark.parametrize("remove_release", ["2.0.0", None])
    @pytest.mark.parametrize("hard_deprecate_release", ["1.4.0", None])
    @pytest.mark.parametrize(
        "callable_name",
        ["myfunction", "mymethod", "myclassmethod", "mystaticmethod", "MyClass"],
    )
    def test_function_decorator(
        self,
        callable_name,
        hard_deprecate_release,
        remove_release,
        current_release,
        decorate_kwarg,
    ):
        """Test all permutations of the callable to decorate, decorate kwarg or not and hard
        deprecate or not

        Args:
            callable_name (str): The name of the kind of callable to test
            decorate_kwarg (bool): Whether to test deprecating a kwarg in the callable
            hard_deprecate (bool): Whether to test hard deprecation

        """
        # First form the decorator
        decorator = deprecate(
            **DEPRECATE_STANDARD_ARGS,
            hard_deprecation_release=hard_deprecate_release,
            remove_release=remove_release,
            kwarg_name=decorate_kwarg,
        )

        # Then get the class and functions we may need
        MyClass, DecoratedMyClass, myfunction = generate_function_and_classes(decorator)

        # Extract or set the appropriate callable, according to what we are testing
        if callable_name == "myclassmethod":
            decorated_callable = getattr(MyClass, callable_name)
        elif callable_name == "MyClass":
            decorated_callable = DecoratedMyClass
        elif callable_name == "myfunction":
            decorated_callable = myfunction
        else:
            my_obj = MyClass(1)
            decorated_callable = getattr(my_obj, callable_name)

        # Determine whether we are in hard deprecation
        hard_deprecate = hard_deprecate_release is not None and version.parse(
            current_release
        ) >= version.parse(hard_deprecate_release)

        # Patch out the imported version of __version__ with test value
        with patch("ixdat.tools.__version__", current_release):
            # If we are testing deprecating a kwarg
            if decorate_kwarg:
                message = self._test_decorate_kwarg(
                    decorated_callable, MyClass, hard_deprecate
                )
            else:
                message = self._test_decorate_callable(
                    decorated_callable, MyClass, hard_deprecate
                )

        # Make sure all the `deprecate` text inputs are in the exception or warning
        # messages
        required_text_snippets = list(DEPRECATE_STANDARD_ARGS.values()) + [
            callable_name,
            "SOMETHING ELSE",
        ]

        if hard_deprecate_release is not None:
            required_text_snippets.append(hard_deprecate_release)
        else:
            required_text_snippets.append(
                "soft deprecated, issuing warnings, for the foreseeable future"
            )

        if remove_release:
            required_text_snippets.append(remove_release)

        for required_text_snippet in required_text_snippets:
            assert required_text_snippet in message

    def _test_decorate_kwarg(self, decorated_callable, MyClass, hard_deprecate):
        """Test deprecation of kwarg"""
        # Make sure no warning is issues if we do not use the kwarg
        # Instriuctions for how to test no warnings is from here:
        # https://docs.pytest.org/en/latest/how-to/capture-warnings.html#additional-use-
        # cases-of-warnings-in-tests
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            return_value = decorated_callable(3)
        if decorated_callable.__name__ == "MyClass":
            assert isinstance(return_value, MyClass)
        else:
            assert return_value == 52

        # In the case of hard deprecation, test that it raises an exception, otherwise it
        # should raise a warning
        if hard_deprecate:
            with pytest.raises(DeprecationError) as exception:
                decorated_callable(3, intkwarg=3)
            message = str(exception.value)
        else:
            with pytest.warns(DeprecationWarning) as captured_warnings:
                return_value = decorated_callable(3, intkwarg=3)
            if decorated_callable.__name__ == "MyClass":
                assert isinstance(return_value, MyClass)
            else:
                assert return_value == 53
            warning = captured_warnings.pop()
            message = str(warning.message)

        # Make sure deprecation is triggered both when using the kwarg as a kwarg and
        # when reaching it as an arg
        if hard_deprecate:
            with pytest.raises(DeprecationError) as exception:
                decorated_callable(3, 3)
            message2 = str(exception.value)
        else:
            with pytest.warns(DeprecationWarning) as captured_warnings:
                return_value = decorated_callable(3, 3)
            if decorated_callable.__name__ == "MyClass":
                assert isinstance(return_value, MyClass)
            else:
                assert return_value == 53
            warning = captured_warnings.pop()
            message2 = str(warning.message)

        # No-matter how we reached the deprecation, the message should be the same
        assert message == message2
        return message

    def _test_decorate_callable(self, decorated_callable, MyClass, hard_deprecate):
        """Test ordinary deprecation of callable"""
        # In the case of hard deprecation, test that it raises an exception, otherwise it
        # should raise a warning
        if hard_deprecate:
            with pytest.raises(DeprecationError) as exception:
                decorated_callable(3, 3)
            message = str(exception.value)
        else:
            with pytest.warns(DeprecationWarning) as captured_warnings:
                return_value = decorated_callable(3, 3)
            if decorated_callable.__name__ == "MyClass":
                assert isinstance(return_value, MyClass)
            else:
                assert return_value == 53
            warning = captured_warnings.pop()
            message = str(warning.message)

        return message
