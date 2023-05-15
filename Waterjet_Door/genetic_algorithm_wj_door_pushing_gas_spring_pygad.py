# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:24:55 2023

@author: Ryan.Larson
"""

import numpy as np
import pygad

def fitness_func(ga_instance, solution, solution_idx):
    # Read in parameter values from parameter list
    ds = solution[0]
    hs = solution[1]
    S = solution[2]
    xp = solution[3]
    yp = solution[4]
    l_ext = solution[5]
    l_comp = solution[6]
    
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
    nsprings = 4                # Number of torsion springs
    torque_rating = 2000.0      # Torsion spring torque rating (in*lbs) (4 1877.0 rated springs balances the door at 90 degrees)
    kappa_tor = torque_rating/theta_max_rad

    ### Derived parameters ###
    # S:
    nu = np.arctan(hs/ds)
    sigma = thetarad - nu
    Ls = np.sqrt(ds**2 + hs**2)
    xs = Ls*np.cos(sigma)
    ys = Ls*np.sin(sigma)
    
    # Line equation for gas spring
    m_closed = (ys[0] - yp) / (xs[0] - xp)
    m_open = (ys[-1] - yp) / (xs[-1] - xp)
    
    b_closed = ys[0] - m_closed*xs[0]
    b_open = ys[-1] - m_open*xs[-1]
    
    # Corner points that must not be intersected
    x_corner_closed = 4.75
    y_corner_closed = 2.25 - 0.5
    x_corner_open = 2.25
    y_corner_open = -4.75 - 0.5
    
    # Check if the gas spring line is below the corner point
    spring_closed_ht = m_closed*x_corner_closed + b_closed
    spring_open_ht = m_open*x_corner_open + b_open

    # W:
    xw = wdoor*np.cos(thetarad) + Lcm*np.sin(thetarad)
    yw = Lcm*np.cos(thetarad) - wdoor*np.sin(thetarad)
    Lw = np.sqrt(wdoor**2 + Lcm**2)

    # F:
    xf = wdoor*np.cos(thetarad) + Ldoor*np.sin(thetarad)
    yf = Ldoor*np.cos(thetarad) - wdoor*np.sin(thetarad)
    Lf = np.sqrt(wdoor**2 + Ldoor**2)

    # Gas spring length
    ls = np.sqrt((np.abs(xp) + xs)**2 + (yp - ys)**2)
    maxspring = np.around(max(ls), 1)
    minspring = np.around(min(ls), 1)
    
    if maxspring > l_ext:
        feasible = False
    if minspring < l_comp:
        feasible = False
    if spring_closed_ht > y_corner_closed:
        feasible = False
    if spring_open_ht > y_corner_open:
        feasible = False

    # Angles
    gamma = np.arctan(yw/xw)
    kappa = np.arctan(yf/xf)
    psi = np.pi/2 - gamma
    phi = np.arctan((yp-ys)/(np.abs(xp)+xs))
    beta = sigma + phi
    
    # Actual closing force curve
    close_angles = np.array([1, 30, 45, 60, 90])
    close_angles_rad = (np.pi/180)*close_angles
    close_force = np.array([22, 70, 80, 89, 0])
    # Develop function for new angles
    z_close = np.polyfit(close_angles_rad, close_force, 3)
    f_close = np.poly1d(z_close)

    ### Force and extension ###
    force = np.zeros(thetarad.shape)
    for i in range(len(thetarad)):
        angle = thetarad[i]
        # force[i] = (W*Lw*np.sin(psi[i]) - ns*S*Ls*np.sin(beta[i])  - kappa_tor*angle*nsprings) / Lf
        door_component = f_close(angle)
        gas_spring_component = ns*S*Ls*np.sin(beta[i]) / Lf

        force[i] = door_component - gas_spring_component

    # Minimize the range between the min and max forces
    minabsforce = np.abs(min(force))
    maxabsforce = np.abs(max(force))
    sumforce = minabsforce + maxabsforce
    
    fitness = 1.0 / sumforce

    return fitness


num_generations = 200
num_parents_mating = 2

fitness_function = fitness_func

sol_per_pop = 8
num_genes = 7
# num_genes = len(function_inputs)

init_range_low = -2
init_range_high = 5

parent_selection_type = "sss"
keep_parents = 1

crossover_type = "single_point"

mutation_type = "random"
mutation_percent_genes = 10

gene_space = [{'low': 12.0, 'high': 17.0},
              {'low': -0.9, 'high': 0.75},
              {'low': 240.0, 'high':250.0},
              {'low': -7.5, 'high': -4.75},
              {'low': -0.85, 'high': 2.0},
              {'low': 16.0, 'high': 35.5},
              {'low': 17.63, 'high': 22.0}]

ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       fitness_func=fitness_function,
                       sol_per_pop=sol_per_pop,
                       gene_space=gene_space,
                       num_genes=num_genes,
                       init_range_low=init_range_low,
                       init_range_high=init_range_high,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       crossover_type=crossover_type,
                       mutation_type=mutation_type,
                       mutation_percent_genes=mutation_percent_genes)

ga_instance.run()

ga_instance.plot_fitness()

solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("Parameters of the best solution : {solution}".format(solution=solution))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
print("Index of the best solution : {solution_idx}".format(solution_idx=solution_idx))

if ga_instance.best_solution_generation != -1:
    print("Best fitness value reached after {best_solution_generation} generations.".format(best_solution_generation=ga_instance.best_solution_generation))