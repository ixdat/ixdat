import numpy as np

c = 2.997925e8  # speed of light / (m/s)
qe = 1.60219e-19  # fundamental charge / (C)
h = 6.62620e-34  # planck's constant / (J*s)
hbar = h / (2 * np.pi)  # reduced planck's constant / (J*s)
NA = 6.02217e23  # Avogadro's number /(mol) or dimensionless
me = 9.10956e-31  # mass of electron / (kg)
kB = 1.38062e-23  # Boltzman constant / (J/K)

u0 = 4 * np.pi * 1e-7  # permeability of free space / (J*s^2/(m*C^2))
e0 = 1 / (u0 * c ** 2)  # permittivity of free space / (C^2/(J*m))

R = NA * kB  # gas constant / (J/(mol*K))                 #NA in /mol
amu = 1e-3 / NA  # atomic mass unit / (kg)    # amu=1g/NA     #NA dimensionless
Far = NA * qe  # Faraday's constant, C/mol

FARADAY_CONSTANT = Far
AVOGADROS_CONSTANT = NA
BOLTZMAN_CONSTANT = kB

STANDARD_TEMPERATURE = 298.15  # Standard temperature of 25 C in [K]
STANDARD_PRESSURE = 1e5  # Standard pressure of 1 bar in [Pa]

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
    "CO2": 3.34-10,
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
