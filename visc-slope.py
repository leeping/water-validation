#!/usr/bin/env python

import numpy as np
from simtk.unit import *
import sys

# First argument: Temperature in kelvin
# Second argument: Slope of the diffusion constant vs. reciprocal box length
# in units of 1e-5 cm^2 s^-1 Angstrom
# Third argument: Uncertainty.

X = -1.0e-5 * (centimeter ** 2) * angstrom * float(sys.argv[2]) / second
eta = 2.837 * BOLTZMANN_CONSTANT_kB * float(sys.argv[1]) * kelvin / ((6.0 * np.pi) * X)

X1 = -1.0e-5 * (centimeter ** 2) * angstrom * (float(sys.argv[2]) - float(sys.argv[3])) / second

deta = eta - (2.837 * BOLTZMANN_CONSTANT_kB * float(sys.argv[1]) * kelvin / ((6.0 * np.pi) * X1))

# Viscosity is in millipascal seconds, so
# m L^-1 / T
print float(sys.argv[1]), eta / millipascal / second, deta / millipascal / second
