"""door-closure-optimize.py: Use GEKKO optimization library to find the optimal
design for waterjet door closure mechanism.

Created: 7 Dec 2021
Author: Ryan Larson

"""

from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

m = GEKKO(remote=True)

# Design variables
tdead = m.Var(15.0, lb=0.0, ub=90.0)
torque = m.Var(2000.0, lb=100.0, ub=10000.0)
nsprng_init = m.Var(2, lb=0, ub=20, integer=True)
nsprng_dead = m.Var(4, lb=0, ub=20, integer=True)

# # Initial values
# tdead.value = 15.0
# torque.value = 2000.0
# nsprng_init.value = 2
# nsprng_dead.value = 4
#
# # Set nsprng_init and nsprng_dead as integer variables
# nsprng_init.integer = True
# nsprng_dead.integer = True
#
# # Upper bounds
# tdead.upper = 90.0
# torque.upper = 10000.0
# nsprng_init.upper = 20
# nsprng_dead.upper = 20
#
# # Lower bounds
# tdead.lower = 0.0
# torque.lower = 100.0
# nsprng_init.lower = 0
# nsprng_dead.lower = 0

# Other analysis variables
min_open = m.Const(value=1.0)
max_open = m.Const(value=90.0)
npts = m.Const(value=90)
theta = m.Param(value=np.linspace(min_open, max_open, npts, endpoint=True))
force = np.zeros(theta.shape)

theta_max = m.Const(value=90.0)    # Deflection at which springs achieve maximum torque, deg
theta_max_rad = m.Const(value=theta_max*np.pi/180.0)
kappa = m.Const(value=torque_rating/theta_max_rad)  # Torsional spring constant (lbf*in/rad)
weight = m.Const(value=600.0)
height = m.Const(value=45.0)

# Analysis functions
for i in range(len(theta)):
    angle = theta[i]
    if angle < tdead:
        force[i] = (weight/2.0)*np.sin(angle*np.pi/180.0) - nsprng_init*(kappa*(angle-tdead)*np.pi)/(180.0*height)
    else:
        force[i] = (weight/2.0)*np.sin(angle*np.pi/180.0) - (nsprng_init+nsprng_dead)*(kappa*(angle-tdead)*np.pi)/(180.0*height)

maxforce = np.abs(max(force))   # Calculate the maximum force experienced

# Constraints
m.Equation(maxforce<=120.0)

# Objective function
m.Minimize(maxforce)

# Set global options
m.options.IMODE = 3 # steady state optimization

# Solve simulation
m.solve()

# Results
print("\nDead angle: {} degrees".format(tdead.value))
print("Spring torque rating: {} lbf*in".format(torque.value))
print("Number of initial springs: {}".format(nsprng_init.value))
print("Number of springs added after dead angle: {}".format(nsprng_dead.value))
