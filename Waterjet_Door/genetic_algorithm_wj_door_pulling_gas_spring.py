# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:17:54 2022

@author: Ryan.Larson
"""

import numpy as np
import matplotlib.pyplot as plt


def objective(parameters):
    """Given parameter values, calculate the force curve of the design.
    """
    # Read in parameter values from parameter list
    ds = parameters[0]
    S = parameters[1]
    wp = parameters[2]
    hp = parameters[3]
    l_ext = parameters[4]
    l_comp = parameters[5]
    
    # Assume the design is feasible
    feasible = True

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
    
    if lsmax > l_ext:
        feasible = False
    if lsmin < l_comp:
        feasible = False

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

    return sumforce, feasible


def generate_feasible_individual():
    # Define range for inputs
    ds_min, ds_max = 1.0, 30.0
    S_min, S_max = 5.0, 500.0
    wp_min, wp_max = -15.0, 5.0
    hp_min, hp_max = 10.0, 16.5
    l_ext_min, l_ext_max = 15.0, 35.0
    l_comp_min, l_comp_max = 5.0, 15.0
    
    l_ratio_min = 0.55
    
    feasible = False
    
    while feasible == False:
        ds = np.random.uniform(ds_min, ds_max)
        S = np.random.uniform(S_min, S_max)
        wp = np.random.uniform(wp_min, wp_max)
        hp = np.random.uniform(hp_min, hp_max)
        l_ext = np.random.uniform(l_ext_min, l_ext_max)
        l_comp = np.random.uniform(l_comp_min, l_comp_max)
        
        feasible = check_design_feasibility(ds, S, wp, hp, l_ext, l_comp, l_ratio_min)
            
    parameters = [ds, S, wp, hp, l_ext, l_comp]
    
    return parameters


def check_design_feasibility(ds, S, wp, hp, l_ext, l_comp, l_ratio_min):
    feasible = False
    
    if l_comp/l_ext >= l_ratio_min:
        feasible = True
    else:
        print("\nInfeasible compression ratio")
        feasible = False
        
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
    
    if lsmax > l_ext:
        feasible = False
        print("Infeasible extension")
    if lsmin < l_comp:
        feasible = False
        print("Infeasible compression")
        
    return feasible
    

def genetic_algorithm(n_pop):
    # Initial population
    pop = [generate_feasible_individual() for x in range(n_pop)]
    
    return pop
    

if __name__ == '__main__':
    pop = genetic_algorithm(20)