"""Module dealing with Units

Important: The only import from pint in ixdat is in this module. This is essential
for ureg, since units of seperately initiated UnitRegistries cannot be compared.
"""

import warnings
from pint import UnitRegistry, UndefinedUnitError, DimensionalityError, Quantity



ureg = UnitRegistry(autoconvert_offset_to_baseunit = False)


class Unit:
    """TODO: flesh out this class or find an appropriate 3rd-party to use instead"""

    def __init__(self, name):
        self.name = name or ""
        try:
            self.u = ureg(self.name)
        except UndefinedUnitError:
            self.u = None
            
    def set_unit(self, new_unit_name):
        if ureg(self.name) == ureg(""):
            self.name = new_unit_name
            self.u = ureg(self.name)
        else:
            warnings.warn(
                f"DataSeries already has assigned a unit {self.u} with the name" 
                f" {self.name}. Consider using 'reset_unit()' if the unit is wrong",
                stacklevel=2,
            )
                             
    def reset_unit(self,new_unit_name):
        warnings.warn(
            f"resetting unit of dataseries to {new_unit_name} from unit {self.u} with the name {self.name}.",
            stacklevel=2,
        )
        self.name = new_unit_name
        self.u = ureg(self.name)


    def __repr__(self):
        return f"Unit('{self.name}')"

    def __eq__(self, other):
        return self.name == other.name
