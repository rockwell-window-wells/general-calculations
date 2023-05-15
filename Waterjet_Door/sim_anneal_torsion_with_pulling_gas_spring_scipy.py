# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 07:51:37 2022

@author: Ryan.Larson
"""

"""Optimization by simulated annealing of the combined torsion spring and
pushing gas spring concept for the waterjet door.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import dual_annealing

def objective(parameters):
    """Given parameter values, calculate the force curve of the design.
    """
    # Read in parameter values from parameter list
    ds = parameters[0]
    S = parameters[1]
    wp = parameters[2]
    hp = parameters[3]
    

    ### Fixed parameters ###
    Ldoor = 45.0            # Length of door
    cm_factor = 2/3
    Lcm = cm_factor*Ldoor   # Length along door to center of mass
    W = 250.0               # Weight of door (lbf)

    wdoor = 5.75    # Length from the pivot point to the line of the door (the short part of the L shape made by the door)
    ns = 2          # Number of gas springs in the design

    min_angle = 0.0
    max_angle = 90.0
    npts = 90
    theta = np.linspace(min_angle, max_angle, num=npts, endpoint=True)
    thetarad = theta*np.pi/180.0

    ### Torsion spring parameters ###
    theta_max_rad = np.pi/2.0
    nsprings = 2                # Number of torsion springs
    torque_rating = 2000.0      # Torsion spring torque rating (in*lbs) (4 1877.0 rated springs balances the door at 90 degrees)
    kappa_tor = torque_rating/theta_max_rad

    ### Derived parameters ###
    # S:
    xs = wdoor*np.cos(thetarad) + ds*np.sin(thetarad)
    ys = ds*np.cos(thetarad) - wdoor*np.sin(thetarad)
    Ls = np.sqrt(wdoor**2 + ds**2)

    # W:
    xw = wdoor*np.cos(thetarad) + Lcm*np.sin(thetarad)
    yw = Lcm*np.cos(thetarad) - wdoor*np.sin(thetarad)
    Lw = np.sqrt(wdoor**2 + Lcm**2)

    # F:
    xf = wdoor*np.cos(thetarad) + Ldoor*np.sin(thetarad)
    yf = Ldoor*np.cos(thetarad) - wdoor*np.sin(thetarad)
    Lf = np.sqrt(wdoor**2 + Ldoor**2)

    # Gas spring length
    ls = np.sqrt((xs - wp)**2 + (hp - ys)**2)
    lsmin = np.min(ls)
    lsmax = np.max(ls)

    # Angles
    phi = np.arctan((hp-ys)/(xs-wp))
    beta = np.arctan(ys/xs)
    gamma = np.arctan(yw/xw)
    kappa = np.arctan(yf/xf)
    alpha = beta + phi
    psi = np.pi/2 - gamma


    ### Force and extension ###
    force = np.zeros(thetarad.shape)
    for i in range(len(thetarad)):
        angle = thetarad[i]
        force[i] = (W*Lw*np.sin(psi[i]) - ns*S*Ls*np.sin(alpha[i]) - kappa_tor*angle*nsprings) / Lf

    # Minimize the range between the min and max forces
    minabsforce = np.abs(min(force))
    maxabsforce = np.abs(max(force))
    sumforce = minabsforce + maxabsforce

    return sumforce
    

# Define range for inputs
ds_min, ds_max = 1.0, 30.0
S_min, S_max = 5.0, 500.0
wp_min, wp_max = -15.0, 5.0
hp_min, hp_max = 10.0, 16.5

# Define bounds on the search
bounds = [[ds_min, ds_max], [S_min, S_max], [wp_min, wp_max], [hp_min, hp_max]]

# Perform the simulated annealing search
result = dual_annealing(objective, bounds)

# Summarize the result
print("Status: {}".format(result["message"]))
print("Total Evaluations: {}".format(result["nfev"]))

# Evaluate solution
solution = result['x']
evaluation = objective(solution)
print("Solution: f({}) = {}".format(solution, evaluation))
