import numpy as np
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
AVOGADRO_CONSTANT = NA
BOLTZMANN_CONSTANT = kB

# standard conditions
STANDARD_TEMPERATURE = 298.15  # Standard temperature of 25 C in [K]
STANDARD_PRESSURE = 1e5  # Standard pressure of 1 bar in [Pa]

# molecule properties (should probably come from elsewhere).
# TODO: Move these viscosities to external MS quantification package
#    They are only used by MSInlet, which is deprecated.
#    See https://github.com/ixdat/ixdat/issues/122

DYNAMIC_VISCOSITIES = {  # Values is found on engineeringtoolbox.com [22C30]
    "O2": np.array(
        [
            [200, 14.72e-06],  # Temperature [K], Dynamic Viscosit [Pa*s]
            [220, 15.98e-06],
            [240, 17.20e-06],
            [260, 18.38e-06],
            [280, 19.53e-06],
            [300, 20.65e-06],
            [320, 21.74e-06],
            [340, 22.80e-06],
            [360, 23.84e-06],
            [400, 25.84e-06],
            [500, 30.49e-06],
            [600, 34.73e-06],
            [700, 38.65e-06],
            [800, 42.33e-06],
        ]
    ),
    "N2": np.array(
        [
            [200, 12.91e-06],
            [220, 13.97e-06],
            [240, 15.00e-06],
            [260, 15.99e-06],
            [280, 16.96e-06],
            [300, 17.89e-06],
            [320, 18.80e-06],
            [340, 19.68e-06],
            [360, 20.55e-06],
            [400, 22.21e-06],
            [500, 26.06e-06],
            [600, 29.58e-06],
            [700, 32.83e-06],
            [800, 35.89e-06],
        ]
    ),
    "Ar": np.array(
        [
            [273.15, 2.10e-05],
            [293.15, 2.23e-05],
            [298.15, 2.27e-05],
            [323.15, 2.42e-05],
            [373.15, 2.73e-05],
            [473.15, 3.28e-05],
            [573.15, 3.77e-05],
            [673.15, 4.22e-05],
            [773.15, 4.64e-05],
            [873.15, 5.04e-05],
        ]
    ),
    "He": np.array(
        [
            [273.15, 1.87e-05],
            [293.15, 1.96e-05],
            [323.15, 2.10e-05],
            [373.15, 2.32e-05],
            [473.15, 2.73e-05],
            [573.15, 3.12e-05],
            [673.15, 3.48e-05],
            [773.15, 3.84e-05],
            [873.15, 4.18e-05],
        ]
    ),
    "CO": np.array(
        [
            [273.15, 1.66e-05],
            [293.15, 1.74e-05],
            [323.15, 1.88e-05],
            [373.15, 2.10e-05],
            [473.15, 2.52e-05],
            [573.15, 2.90e-05],
            [673.15, 3.25e-05],
            [773.15, 3.56e-05],
            [873.15, 3.86e-05],
        ]
    ),
    "CO2": np.array(
        [
            [220, 11.06e-06],
            [240, 12.07e-06],
            [250, 12.57e-06],
            [260, 13.06e-06],
            [280, 14.05e-06],
            [300, 15.02e-06],
            [320, 15.98e-06],
            [340, 16.93e-06],
            [360, 17.87e-06],
            [400, 19.70e-06],
            [450, 21.90e-06],
            [500, 24.02e-06],
            [600, 28.00e-06],
            [650, 29.87e-06],
            [850, 36.71e-06],
            [1050, 42.69e-06],
        ]
    ),
    "H2": np.array(
        [
            [273.15, 0.84e-05],
            [293.15, 0.88e-05],
            [323.15, 0.94e-05],
            [373.15, 1.04e-05],
            [473.15, 1.21e-05],
            [573.15, 1.37e-05],
            [673.15, 1.53e-05],
            [773.15, 1.69e-05],
            [873.15, 1.84e-05],
        ]
    ),
}
MOLECULAR_DIAMETERS = {
    "O2": 3.55e-10,
    "N2": 3.7e-10,
    "Ar": 3.58e-10,
    "He": 2.15e-10,
    "CO": 3.76e-10,
    "CO2": 3.34e-10,
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
