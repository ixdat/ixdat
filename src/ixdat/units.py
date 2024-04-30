"""Module dealing with Units

Important: The only import from pint in ixdat is in this module. This is essential
for ureg, since units of seperately initiated UnitRegistries cannot be compared.
"""
from pint import UnitRegistry, UndefinedUnitError, DimensionalityError, Quantity
import warnings

ureg = UnitRegistry()

ureg.setup_matplotlib(True)  # I think this should be closer to the measurement or specific plotter
ureg.autoconvert_offset_to_baseunit = True  # This is tmo help with logarithmic plotting. I dont know if it is better to simply set the dimensions to "dimensionless" prior to plotting with units



class Unit:
    """TODO: flesh out this class or find an appropriate 3rd-party to use instead"""

    def __init__(self, name):
        self.name = name or ""
        try:
            self.u = ureg(self.name)
        except UndefinedUnitError:
            self.u = None
            
                
    def set_unit(self,new_unit_name):
        warnings.warn(
            f"setting unit of dataseries to {new_unit_name} from unit {self.u}"
            " with the name {self.name}. This does NOT alter the values in "
            "this dataseries. "
            f"If conversion to this unit {new_unit_name} is expected this "
            "can be done with inplace_to_unit() from the measurement class",
            stacklevel=2,
        )
        self.name = new_unit_name

        self.u = ureg(self.name)
                     

    def __repr__(self):
        return f"Unit('{self.name}')"

    def __eq__(self, other):
        return self.name == other.name


STANDARD_UNITS = {
    "current":"A",
    "time":"s",
    "charge":"C",
    "electrical_potential":"V",
    "resistance":"Ω",
    "temperature":"K",
    "pressure":"mbar"
}