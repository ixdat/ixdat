"""Module dealing with Units

Important: The only import from pint in ixdat is in this module. This is essential
for ureg, since units of seperately initiated UnitRegistries cannot be compared.
"""
from pint import UnitRegistry, UndefinedUnitError, DimensionalityError, Quantity


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
