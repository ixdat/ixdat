from scipy import constants as scipy_constants

# short-form aliases for a few scipy constants
c = scipy_constants.c  # speed of light / (m/s)
qe = scipy_constants.e  # fundamental charge / (C)
h = scipy_constants.h  # planck's constant / (J*s)
hbar = scipy_constants.hbar  # reduced planck's constant / (J*s)
NA = scipy_constants.N_A  # Avogadro's number /(mol) or dimensionless
me = scipy_constants.m_e  # mass of electron / (kg)
kB = scipy_constants.k  # Boltzman constant / (J/K)
u0 = scipy_constants.mu_0  # permeability of free space / (J*s^2/(m*C^2))
e0 = scipy_constants.epsilon_0  # permittivity of free space / (C^2/(J*m))
R = scipy_constants.R  # gas constant / (J/(mol*K))

# a few extra derived constants
amu = 1e-3 / NA  # atomic mass unit / (kg)    # amu=(1g/mol)/NA
Far = NA * qe  # Faraday's constant, C/mol

# long-form aliases
FARADAY_CONSTANT = Far
AVOGADROS_CONSTANT = NA
BOLTZMAN_CONSTANT = kB

# standard conditions
STANDARD_TEMPERATURE = 298.15  # Standard temperature of 25 C in [K]
STANDARD_PRESSURE = 1e5  # Standard pressure of 1 bar in [Pa]

# molecule properties (should probably come from elsewhere).
DYNAMIC_VISCOSITIES = {
    "O2": 2.07e-05,
    "N2": 1.79e-05,
    "Ar": 2.27e-05,
    "He": 1.99e-05,
    "CO": 1.78e-05,
    "CO2": 1.503e-05,
    "H2": 8.90e-06,
}  # in [Pa*s]
MOLECULAR_DIAMETERS = {
    "O2": 3.55e-10,
    "N2": 3.7e-10,
    "Ar": 3.58e-10,
    "He": 2.15e-10,
    "CO": 3.76e-10,
    "CO2": 3.34 - 10,
    "H2": 2.71e-10,
}  # in [m]
MOLAR_MASSES = {
    "O2": 31.998,
    "N2": 28.014,
    "Ar": 39.948,
    "He": 4.002,
    "CO": 28.010,
    "CO2": 44.010,
    "H2": 2.016,
}  # in [g/mol]
