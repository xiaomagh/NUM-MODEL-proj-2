# -*- coding: utf-8 -*-
"""
student number: 25806676
"""

from __future__ import division
import numpy as np


dict = {
    # The grid-length, d, unit is m
    'd' : 25000.
        }

type(dict)

# The gravitational acceleration, g, unit is m/s^2
g = 10.

# The linear drag coefficient, gama, unit is s^-1
gama = 1e-06

# The uniform density, ro, unit is kg/m^3
ro = 1000.

# The resting depth of the fluid which is assumed constant, H, unit is m
H = 1000.

# The range of the numerical square, L, unit is m
L = 1e06

# The time step length, dt, unit is s
dt = 120.

# The number of grids, dimensionless
nd = int(L/dict['d'])

# The parameter for calculating Coriolis parameter in Beta-plane approximation
beta = 1e-11

# The initial Coriolis parameter in Beta-plane approximation
f0 = 1e-04

# Function for calculating Coriolis parameter based on beta-plane approximation
def f(y):
    f = f0 + beta*y
    
    return float(f)
    
# Function for the wind stress acting on the surface of the fluid, ta
def ta(y):
    tau0 = 0.2
    taux = -np.cos(np.pi*y/L)*tau0
    
    return float(taux)