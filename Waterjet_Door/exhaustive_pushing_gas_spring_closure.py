"""Optimization by exhaustive grid search of pusher gas spring waterjet door
closure.
"""

import numpy as np
import matplotlib.pyplot as plt

optforce = None
ds = None
S = None
xp = None
yp = None

# Fixed parameters
Ldoor = 45.0            # Length of door
cm_factor = 2/3
Lcm = cm_factor*Ldoor   # Length along door to center of mass
W = 250.0               # Weight of door (lbf)

wdoor = 5.75    # Length from the pivot point to the line of the door (the short part of the L shape made by the door)
ns = 2          # Number of gas springs in the design

min_angle = 1.0
max_angle = 90.0
npts = 90
theta = np.linspace(min_angle, max_angle, num=npts, endpoint=True)
thetarad = theta*np.pi/180.0

### Design parameters ###
npts_var = 10

ds_ub = 23.5
ds_lb = 1.0
ds_vec = np.linspace(ds_lb, ds_ub, npts_var, endpoint=True)
S_ub = 500.0
S_lb = 5.0
S_vec = np.linspace(S_lb, S_ub, npts_var, endpoint=True)       # Force rating (lbf) of the gas spring

xp_ub = -1.0
xp_lb = -40.0
xp_vec= np.linspace(xp_lb, xp_ub, npts_var, endpoint=True)        # Horizontal position of the fixed end of the gas spring, relative to the door pivot point
yp_ub = -1.0
yp_lb = -23.5
yp_vec = np.linspace(yp_lb, yp_ub, npts_var, endpoint=True)         # Vertical position of the fixed end of the gas spring, relative to the door pivot point (max value 23.5 inches, negative)

### Derived parameters ###
for j in range(npts_var):
    for k in range(npts_var):
        for m in range(npts_var):
            for n in range(npts_var):
                # S:
                xs = ds_vec[j]*np.cos(thetarad)
                ys = ds_vec[j]*np.sin(thetarad)

                # W:
                xw = wdoor*np.cos(thetarad) + Lcm*np.sin(thetarad)
                yw = Lcm*np.cos(thetarad) - wdoor*np.sin(thetarad)
                Lw = np.sqrt(wdoor**2 + Lcm**2)

                # F:
                xf = wdoor*np.cos(thetarad) + Ldoor*np.sin(thetarad)
                yf = Ldoor*np.cos(thetarad) - wdoor*np.sin(thetarad)
                Lf = np.sqrt(wdoor**2 + Ldoor**2)

                # Gas spring length
                ls = np.sqrt((np.abs(xp_vec[m]) + xs)**2 + (yp_vec[n] - ys)**2)

                # Angles
                gamma = np.arctan(yw/xw)
                kappa = np.arctan(yf/xf)
                psi = np.pi/2 - gamma
                phi = np.arctan((yp_vec[n]-ys)/(np.abs(xp_vec[m])+xs))
                beta = thetarad - phi


                ### Force and extension ###
                force = np.zeros(thetarad.shape)
                for i in range(len(thetarad)):
                    force[i] = (W*Lw*np.sin(psi[i]) - ns*S_vec[k]*ds_vec[j]*np.sin(beta[i])) / Lf

                absforce = np.abs(force)
                sumforce = sum(absforce)

                # Get maximum spring extension
                maxspring = np.around(max(ls), 1)
                minspring = np.around(min(ls), 1)
                extension = np.around(maxspring - minspring, 1)
                if maxspring > 2*minspring:
                    feasible = False
                else:
                    feasible = True

                if feasible is True and optforce is None:
                    optforce = sumforce
                    ds = ds_vec[j]
                    S = S_vec[k]
                    xp = xp_vec[m]
                    yp = yp_vec[n]
                elif feasible is True and sumforce < optforce:
                    optforce = sumforce
                    ds = ds_vec[j]
                    S = S_vec[k]
                    xp = xp_vec[m]
                    yp = yp_vec[n]


print("Optimal design:")
print("ds: {}".format(ds))
print("S: {}".format(S))
print("xp: {}".format(xp))
print("yp: {}".format(yp))
