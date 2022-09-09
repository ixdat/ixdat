"""Patches corresponding to code that should be added to Spectro Inlets quantification"""


def append_sensitivity_factors(*args):
    """Recursively append elements in a list"""
    sf = append_two_sensitivity_factors(args[0], args[1])
    if len(args) == 2:
        return sf
    remaining_args = [sf] + list(args[2:])
    return append_sensitivity_factors(*remaining_args)


def append_two_sensitivity_factors(sf_1, sf_2):
    """This should be incorporated into __add__ of SensitivityFactor and SensitivityList

    Appends all the sensitivity factors in sf_1 and sf_2, which can each be a single
    SensitivityFactor (or CalPoint) or a SensitivityList (or Calibration).

    Args:
        sf_1 (SensitivityFactor or SensitivityList): The first sensitivity factor or list
            thereof.
        sf_2 (CalPoint or Calibration): The second sensitivity factor or list thereof

    Returns SensitivityList: a `SensitivityList` containing all of the sensitivity
        factors in `sf_1` and `sf_2`. This will be a `Calibration` unless both `sf_1` and
        `sf_2` are both themselves plain `SensitivityList`s
    """
    from spectro_inlets_quantification.sensitivity import (
        SensitivityList,
        SensitivityFactor,
    )
    from spectro_inlets_quantification.calibration import Calibration

    if isinstance(sf_1, SensitivityList):
        if isinstance(sf_2, SensitivityList):
            return sf_1 + sf_2  # this is implemented in spectro_inlets_quantification
        elif isinstance(sf_2, SensitivityFactor):
            sf_list = sf_1.sf_list + [sf_2]
            obj_as_dict = sf_1.as_dict()
            if isinstance(sf_1, Calibration):
                sf_1.pop("cal_dicts")
                obj_as_dict["cal_list"] = sf_list
                return Calibration(**obj_as_dict)
            else:
                sf_1.pop("sf_dicts")
                obj_as_dict["sf_list"] = sf_list
                return SensitivityList(**obj_as_dict)
        else:
            raise TypeError(f"Can't add {sf_1} and {sf_2}")
    elif isinstance(sf_1, SensitivityFactor):
        if isinstance(sf_2, SensitivityFactor):
            sf_list = [sf_1, sf_2]
            return Calibration(cal_list=sf_list)
        else:
            return append_sensitivity_factors(sf_2, sf_1)
    else:
        raise TypeError(f"Can't add {sf_1} and {sf_2}")


def add_isotopes(calibration, isotope_spec):
    """Duplicate sensitivity factor(s) in the calibration to cover different isotopes

    This method adds CalPoints in the calibration for each new tracked isotope. It
    assumes that the sensitivity factor for each isotope at its respective mass is
    the same.

    The isotopes have to be treated as different molecules, or they will end up on the
      same row of a SensitivityMatrix, convoluting their quantification. This means
      that the `mol` attribute of their CalPoints must indicate the isotope, here
      done with "@" and the mass.
    Because SI quant's Quantifier object makes sure each of the molecules in the
      calibration are in its MoleculeDict, Molecules of the new name must be added to
      the MoleculeDict

    Args:
        calibration (Calibration): The calibration to expand
        isotope_spec (dict): A specification of the isotopes to expand the calibration
            with. The keys are molecules and the values are a tuple with the base mass,
            which already exists in the calibration, followed by a list of masses to
            add. An example for CO2 and O2 in 18-O labeling experiments is:
             {"CO2": ("M44", ["M46", "M48"]), "O2": ("M32", ["M34", "M36"])}
    """
    from spectro_inlets_quantification.calibration import CalPoint
    from spectro_inlets_quantification.molecule import MoleculeDict, Molecule

    mdict = MoleculeDict()

    for mol, (mass, new_masses) in isotope_spec.items():
        cal_point = calibration.get(mol, mass)
        molecule_as_dict = mdict.get(mol).as_dict()
        for new_mass in new_masses:
            new_mol = mol + "@" + new_mass
            new_cal_point = CalPoint(
                mol=new_mol, mass=new_mass, F=cal_point.F, F_type=cal_point.F_type
            )
            calibration.append(new_cal_point)
            new_molecule_as_dict = molecule_as_dict.copy()
            new_molecule_as_dict["mol"] = new_mol
            # To avoid quant incorrectly predicting sensitivity factors at other masses,
            #    we set a spectrum predicting intensity only at the specified mass.
            new_molecule_as_dict["spectrum"] = {new_mass: 1}
            new_molecule = Molecule(**new_molecule_as_dict)
            mdict[new_mol] = new_molecule
