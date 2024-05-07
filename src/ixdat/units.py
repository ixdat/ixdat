"""Module dealing with Units

Important: The only import from pint in ixdat is in this module. This is essential
for ureg, since units of seperately initiated UnitRegistries cannot be compared.
"""
from pint import UnitRegistry, UndefinedUnitError, DimensionalityError, Quantity, Context
from pint.util import UnitsContainer
from pint import formatting
import warnings

ureg = UnitRegistry()

  
ureg.setup_matplotlib(True)  # I think this should be closer to the measurement or specific plotter
ureg.mpl_formatter = "{:~ixdat_default}"
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


# Adding some context to go back and forth from mol/s to mol/s/cm^2 
## NEED more work to figure out where and when the A_el constant should be added to the registry and defining this context. 

# def add_context_to_measurement(measurement):
#     ureg.define('_A_el = measurement.A_el')
#     c = Context('area_calibration')
#     c.add_transformation('[substance] / [lenght] **2 / [time]', '[substance] / [time]', lambda ureg, x: x * _A_el)
#     c.add_transformation('[substance] / [time]', '[substance] / [lenght] ** 2 / [time]', lambda ureg, x: x / _A_el)

def register_format_function(plotter_type):
    try:
        @formatting.register_unit_format(f"ixdat_{plotter_type.lower()}")
        def format_function(unit: UnitsContainer, registry: UnitRegistry, **options) -> str:
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
            
            dimension = _get_dimension(unit=unit, plotter=plotter_type)
            if plotter_type == "LOG":
                return f"{dimension} / [ln(" + formatted_unit + ")]"
            return f"{dimension} / [" + formatted_unit + "]"
    except ValueError:
        pass

# Register the format functions for each plotter type
for plotter in ["DEFAULT", "EC", "LOG"]:
    register_format_function(plotter)



def _ixdat_fmt_exponent(num) -> str:
    return "{"+f"{num}"+"}"
    


def _get_dimension(plotter, unit):
    u = ureg.dimensionless
    for key, value in unit.items():
        u *= ureg(key).u ** int(value)
    dimensionality = u.dimensionality
    #print(u, dimensionality)
    #print("cal. sig.", (ureg.mol / ureg.s).dimensionality)
    for label, dimension in STANDARD_LABELS_BY_PLOTTER[plotter].items():
        if dimension.__eq__(dimensionality):
            return label[0] if isinstance(label, tuple) else label
    return "signal"

STANDARD_LABELS_BY_PLOTTER = {
    'DEFAULT':{
        "signal":ureg.A.dimensionality,
        ("cal. sig.","mol/s"):(ureg.mol / ureg.s).dimensionality,
        ("cal. sig.","mol/s/cm"):(ureg.mol / ureg.s / ureg.cm ** 2).dimensionality,
        "time":ureg.s.dimensionality,
        "temperature":ureg.kelvin.dimensionality,
        "pressure":ureg.bar.dimensionality,
    },
    'LOG':{
        "ln(signal)":ureg.A.dimensionality,
        "ln(cal. sig.)":(ureg.mol / ureg.s).dimensionality,
    },
    
    'EC':{
        "E$_{we}$":ureg.V.dimensionality,
        "raw current":ureg.A.dimensionality,
        "time":ureg.s.dimensionality,
}
}
