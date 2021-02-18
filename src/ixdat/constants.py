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
Faraday = NA * qe  # Faraday's constant, C/mol
amu = 1e-3 / NA  # atomic mass unit / (kg)    # amu=1g/NA     #NA dimensionless
