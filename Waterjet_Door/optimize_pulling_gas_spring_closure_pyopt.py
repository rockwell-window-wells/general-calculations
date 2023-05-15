# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:37:15 2022

@author: Ryan.Larson
"""

import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

# Create GEKKO model
m = GEKKO()

# Fixed parameters:
Ldoor = m.Const(value=45.0)         # Length of long leg of L-shaped door
cm_factor = m.Const(value=2/3)      # factor for placing the center of mass of the door
Lcm = m.Const(value=cm_factor*Ldoor)    # Distance along door to center of mass
W = m.Const(value=250.0)            # Door weight (lbs)
ns = m.Const(value=2) # Number of gas springs in the design
wdoor = m.Const(value=5.75)         # Length of short leg of door

min_angle = 1.0
max_angle = 90.0
npts = 90
theta = np.linspace(min_angle, max_angle, num=npts, endpoint=True)

# Define thetarad as a model parameter array
thetarad = m.Array(m.Param, npts, lb=min_angle, ub=max_angle)
for i,ti in enumerate(theta):
    thetarad[i].value = ti*np.pi/180.0

# Design parameters:
ds = m.Var(28.0, lb=1.0, ub=Ldoor)
S = m.Var(200.0, lb=0.0, ub=250.0)
hp = m.Var(17.5, lb=15.0, ub=20.0)
wp = m.Var(2.0, lb=-5.0, ub=5.0)
sumabsforce = m.Var()

# Derived parameters
# S:
xs = [m.Intermediate(wdoor*m.cos(thetarad[i]) + ds*m.sin(thetarad[i])) for i in range(npts)]
ys = [m.Intermediate(ds*m.cos(thetarad[i]) - wdoor*m.sin(thetarad[i])) for i in range(npts)]
Ls = m.Intermediate(m.sqrt(wdoor**2 + ds**2))

# W:
xw = [m.Intermediate(wdoor*m.cos(thetarad[i]) + Lcm*m.sin(thetarad[i])) for i in range(npts)]
yw = [m.Intermediate(Lcm*m.cos(thetarad[i]) - wdoor*m.sin(thetarad[i])) for i in range(npts)]
Lw = m.Intermediate(m.sqrt(wdoor**2 + Lcm**2))

# F:
xf = [m.Intermediate(wdoor*m.cos(thetarad[i]) + Ldoor*m.sin(thetarad[i])) for i in range(npts)]
yf = [m.Intermediate(Ldoor*m.cos(thetarad[i]) - wdoor*m.sin(thetarad[i])) for i in range(npts)]
Lf = m.Intermediate(m.sqrt(wdoor**2 + Ldoor**2))

# Gas spring length
ls = [m.Intermediate(m.sqrt((xs[i] - wp)**2 + (hp - ys[i])**2)) for i in range(npts)]

# Angles
phi = [m.Intermediate(m.atan((hp-ys[i])/(xs[i]-wp))) for i in range(npts)]
beta = [m.Intermediate(m.atan(ys[i]/xs[i])) for i in range(npts)]
gamma = [m.Intermediate(m.atan(yw[i]/xw[i])) for i in range(npts)]
kappa = [m.Intermediate(m.atan(yf[i]/xf[i])) for i in range(npts)]
alpha = [m.Intermediate(beta[i] + phi[i]) for i in range(npts)]
psi = [m.Intermediate(np.pi/2 - gamma[i]) for i in range(npts)]

# Calculate the force required at each angle:
force = [m.Intermediate((W*Lw*m.sin(psi[i]) - ns*S*Ls*m.sin(alpha[i])) / Lf) for i in range(npts)]

# Spring length constraint
maxspring = m.Intermediate(ls[0].value)
minspring = m.Intermediate(100.0)
for i,lsi in enumerate(ls):
    if lsi.value > maxspring.value:
        maxspring.value = lsi.value
    if lsi.value < minspring.value:
        minspring.value = lsi.value

m.Equation(maxspring<(2*minspring))     # Reject unrealistic gas spring extension

absforce = [m.Intermediate(m.abs3(force[i])) for i in range(npts)]
m.Equation(sumabsforce==m.sum(absforce))

m.Minimize(sumabsforce)

m.options.SOLVER=1
m.solve()
print(ds.value[0], S.value[0], hp.value[0], wp.value[0])
