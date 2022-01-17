"""Optimization by exhaustive grid search of pusher gas spring waterjet door
closure.
"""

import numpy as np
import matplotlib.pyplot as plt

optforce = None
ds = None
S = None
wp = None
hp = None

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

ds_ub = 45.0
ds_lb = 1.0
ds_vec = np.linspace(ds_lb, ds_ub, npts_var, endpoint=True)
S_ub = 500.0
S_lb = 5.0
S_vec = np.linspace(S_lb, S_ub, npts_var, endpoint=True)       # Force rating (lbf) of the gas spring

wp_ub = 5.0
wp_lb = -5.0
wp_vec= np.linspace(wp_lb, wp_ub, npts_var, endpoint=True)        # Horizontal position of the fixed end of the gas spring, relative to the door pivot point
hp_ub = 15.0
hp_lb = 17.5
hp_vec = np.linspace(hp_lb, hp_ub, npts_var, endpoint=True)         # Vertical position of the fixed end of the gas spring, relative to the door pivot point (max value 23.5 inches, negative)

### Derived parameters ###
for j in range(npts_var): # ds_vec
    for k in range(npts_var): # S_vec
        for m in range(npts_var): # wp_vec
            for n in range(npts_var): # hp_vec
                # S:
                xs = wdoor*np.cos(thetarad) + ds_vec[j]*np.sin(thetarad)
                ys = ds_vec[j]*np.cos(thetarad) - wdoor*np.sin(thetarad)
                Ls = np.sqrt(wdoor**2 + ds_vec[j]**2)

                # W:
                xw = wdoor*np.cos(thetarad) + Lcm*np.sin(thetarad)
                yw = Lcm*np.cos(thetarad) - wdoor*np.sin(thetarad)
                Lw = np.sqrt(wdoor**2 + Lcm**2)

                # F:
                xf = wdoor*np.cos(thetarad) + Ldoor*np.sin(thetarad)
                yf = Ldoor*np.cos(thetarad) - wdoor*np.sin(thetarad)
                Lf = np.sqrt(wdoor**2 + Ldoor**2)

                # Gas spring length
                ls = np.sqrt((xs - wp_vec[m])**2 + (hp_vec[n] - ys)**2)

                # Angles
                phi = np.arctan((hp_vec[n]-ys)/(xs-wp_vec[m]))
                beta = np.arctan(ys/xs)
                gamma = np.arctan(yw/xw)
                kappa = np.arctan(yf/xf)
                alpha = beta + phi
                psi = np.pi/2 - gamma


                ### Force and extension ###
                force = np.zeros(thetarad.shape)
                for i in range(len(thetarad)):
                    force[i] = (W*Lw*np.sin(psi[i]) - ns*S_vec[k]*Ls*np.sin(alpha[i])) / Lf

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
                    wp = wp_vec[m]
                    hp = hp_vec[n]
                elif feasible is True and sumforce < optforce:
                    optforce = sumforce
                    ds = ds_vec[j]
                    S = S_vec[k]
                    wp = wp_vec[m]
                    hp = hp_vec[n]


print("Optimal design:")
print("ds: {}".format(ds))
print("S: {}".format(S))
print("wp: {}".format(wp))
print("hp: {}".format(hp))
