"""Tests for tools.py"""
import pytest

from ixdat.exceptions import DeprecationError
from ixdat.tools import deprecate

# Standard arguments for deprecate
DEPRECATE_STANDARD_ARGS = {
    "last_release_when_supported": "1.2.3",
    "soft_deprecation_period": "3 months",
    "hard_deprecation_period": "6 months",
}


class TestDeprecate:
    """Test the `deprecate` decorator"""


    @pytest.mark.parametrize(
        "callable_name",
        ["myfunction", "mymethod", "myclassmethod", "mystaticmethod", "MyClass"],
    )
    @pytest.mark.parametrize("decorate_kwarg", [None, "intkwarg"])
    @pytest.mark.parametrize("hard_deprecate", [False, True])
    def test_function_decorator(self, callable_name, hard_deprecate, decorate_kwarg):
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
            message="Do SOMETHING ELSE now",
            hard_deprecate=hard_deprecate,
            kwarg_name=decorate_kwarg,
        )
        if "method" in callable_name or "MyClass" in callable_name:
            # If we need a class, ee build it inside the test function to avoid decorating
            # the methods more than once
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

            # Get the appropriate callable from the class
            if "classmethod" in callable_name:
                decorated_callable = getattr(MyClass, callable_name)
            elif callable_name == "MyClass":
                decorated_callable = decorator(MyClass)
            else:
                my_obj = MyClass(1)
                decorated_callable = getattr(my_obj, callable_name)
        else:
            # We just need a normal function
            def myfunction(intarg, intkwarg=2, one_more_kwarg=None):
                return 47 + intarg + intkwarg

            decorated_callable = decorator(myfunction)

        # If we are testing decorating a kwarg
        if decorate_kwarg:
            # Make sure nothing happens if we do not use the kwarg
            with pytest.warns(None) as captured_warnings:
                return_value = decorated_callable(3)
            if decorated_callable.__name__ == "MyClass":
                assert isinstance(return_value, MyClass)
            else:
                assert return_value == 52
            assert len(captured_warnings) == 0

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

            # Make sure deprecation is triggered both when using the kwarg as a kwarg as a
            # kwarg and when reaching it as a arg
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

        else:  # We are decorating a callable, not a kwarg within it
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

        # Make sure all the `deprecate` text inputs are in the exception or warning message
        for snippet in list(DEPRECATE_STANDARD_ARGS.values()) + [
            callable_name,
            "SOMETHING ELSE",
        ]:
            assert snippet in message
