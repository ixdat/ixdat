"""Module dealing with Units

Important: The only import from pint in ixdat is in this module. This is essential
for ureg, since units of seperately initiated UnitRegistries cannot be compared.
"""
from pint import UnitRegistry, UndefinedUnitError, DimensionalityError, Quantity
from pint.util import UnitsContainer
from pint import formatting
import warnings

ureg = UnitRegistry()


try:
    @formatting.register_unit_format("ixdat")
    def format_ixdat_ms(unit: UnitsContainer, registry: UnitRegistry, technique='MS', **options) -> str:
        formatted_unit = formatting.formatter(
            unit.items(),
            as_ratio=False,
            single_denominator=False,
            product_fmt=" ",
            division_fmt="/",
            power_fmt="{}$^{}$",
            parentheses_fmt="({})",
            exp_call=_ixdat_fmt_exponent,
            **options,
        )
        

        dimension = _get_dimension(unit=unit, technique=technique)
        print("dimension: ", dimension)

        return f"{dimension} / [" + formatted_unit + "]"
except ValueError:
    print("Warning")


def _ixdat_fmt_exponent(num) -> str:
    return "{"+f"{num}"+"}"
    
    

ureg.setup_matplotlib(True)  # I think this should be closer to the measurement or specific plotter
ureg.mpl_formatter = "{:~ixdat}"
ureg.autoconvert_offset_to_baseunit = True  # This is to help with logarithmic plotting. I dont know if it is better to simply set the dimensions to "dimensionless" prior to plotting with units


class Unit:
    """TODO: flesh out this class or find an appropriate 3rd-party to use instead"""

    def __init__(self, name):
        self.name = name or ""
        try:
            self.u = ureg(self.name).u
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
        self.u = ureg(self.name).u

    def __repr__(self):
        return f"Unit('{self.name}')"

    def __eq__(self, other):
        return self.name == other.name




def _get_dimension(technique, unit):
    u = ureg.dimensionless
    for key, value in unit.items():
        u *= ureg(key).u ** int(value)
    dimensionality = u.dimensionality
    print(u, dimensionality)
    print("cal. sig.", (ureg.mol / ureg.s).dimensionality)
    for label, dimension in STANDARD_LABELS_BY_TECHNIQUE[technique].items():
        if dimension.__eq__(dimensionality):
            return label
    return "signal"

STANDARD_LABELS_BY_TECHNIQUE = {
    'MS':{
        "signal":ureg.A.dimensionality,
        "time":ureg.s.dimensionality,
        "temperature":ureg.kelvin.dimensionality,
        "pressure":ureg.bar.dimensionality,
        "cal. sig.":(ureg.mol / ureg.s).dimensionality,
    },
    'reactor':{


    }
}
