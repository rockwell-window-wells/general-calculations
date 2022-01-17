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
ds = m.Var(8.5, lb=1.0, ub=23.5)
S = m.Var(200.0, lb=5.0, ub=500.0)
xp = m.Var(-17.5, lb=-40.0, ub=-1.0)
yp = m.Var(-10.0, lb=-23.5, ub=-1.0)
sumabsforce = m.Var()

# Derived parameters
# S:
xs = [m.Intermediate(ds*m.cos(thetarad[i])) for i in range(npts)]
ys = [m.Intermediate(ds*m.sin(thetarad[i])) for i in range(npts)]
# Ls = m.Intermediate(m.sqrt(wdoor**2 + ds**2))

# W:
xw = [m.Intermediate(wdoor*m.cos(thetarad[i]) + Lcm*m.sin(thetarad[i])) for i in range(npts)]
yw = [m.Intermediate(Lcm*m.cos(thetarad[i]) - wdoor*m.sin(thetarad[i])) for i in range(npts)]
Lw = m.Intermediate(m.sqrt(wdoor**2 + Lcm**2))

# F:
xf = [m.Intermediate(wdoor*m.cos(thetarad[i]) + Ldoor*m.sin(thetarad[i])) for i in range(npts)]
yf = [m.Intermediate(Ldoor*m.cos(thetarad[i]) - wdoor*m.sin(thetarad[i])) for i in range(npts)]
Lf = m.Intermediate(m.sqrt(wdoor**2 + Ldoor**2))

# Gas spring length
ls = [m.Intermediate(m.sqrt((m.abs3(xp) + xs[i])**2 + (yp - ys[i])**2)) for i in range(npts)]

# Angles
gamma = [m.Intermediate(m.atan(yw[i]/xw[i])) for i in range(npts)]
kappa = [m.Intermediate(m.atan(yf[i]/xf[i])) for i in range(npts)]
psi = [m.Intermediate(np.pi/2 - gamma[i]) for i in range(npts)]
phi = [m.Intermediate(m.atan((yp-ys[i])/(m.abs3(xp)+xs[i]))) for i in range(npts)]
beta = [m.Intermediate(thetarad[i] - phi[i]) for i in range(npts)]
# alpha = [m.Intermediate(beta[i] + phi[i]) for i in range(npts)]

# Calculate the force required at each angle:
force = [m.Intermediate((W*Lw*m.sin(psi[i]) - ns*S*ds*m.sin(beta[i])) / Lf) for i in range(npts)]

# Spring length constraint
maxspring = m.Intermediate(0.0) # Maxspring starts from 0 and grows until the max value is found
minspring = m.Intermediate(1000.0) # Minspring starts from 1000 and is diminished until the min is found
for i,lsi in enumerate(ls):
    # print(lsi.value)
    if lsi.value > maxspring.value:
        maxspring.value = lsi.value
    if lsi.value < minspring.value:
        minspring.value = lsi.value

# print("minspring: {}".format(minspring.value))
# print("maxspring: {}".format(maxspring.value))
# maxzero = m.Param(value=0.0)
# maxspring = m.Param(m.max2(ls, maxzero))
# print(maxspring.value)
# min1000 = m.Param(value=1000.0)
# minspring = m.Param(m.min2(ls, min1000))
# print(minspring.value)

m.Equation(maxspring<(2*minspring))     # Reject unrealistic gas spring extension

absforce = [m.Intermediate(m.abs3(force[i])) for i in range(npts)]
m.Equation(sumabsforce==m.sum(absforce))

m.Minimize(sumabsforce)

m.options.SOLVER=1
# m.open_folder()
m.solve()
print(ds.value[0], S.value[0], xp.value[0], yp.value[0], minspring.value, maxspring.value)
