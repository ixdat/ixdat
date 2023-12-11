"""Module dealing with Units

TODO: look into using an available unit's package like astropy's
"""
from pint import UnitRegistry, UndefinedUnitError


ureg = UnitRegistry()


class Unit:
    """TODO: flesh out this class or find an appropriate 3rd-party to use instead"""

    def __init__(self, name):
        self.name = name or ""
        try:
            self.u = ureg(self.name)
        except UndefinedUnitError:
            self.u = None

    def __repr__(self):
        return f"Unit('{self.name}')"

    def __eq__(self, other):
        return self.name == other.name
