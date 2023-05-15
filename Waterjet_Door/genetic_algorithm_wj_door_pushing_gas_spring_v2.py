# -*- coding: utf-8 -*-
"""
Created on Thur Jun 30 11:36:00 2022

@author: Ryan.Larson
"""

import numpy as np
import matplotlib.pyplot as plt


def objective_function(parameters):
    """Given parameter values, calculate the force curve of the design.
    """
    # Read in parameter values from parameter list
    ds = parameters[0]
    hs = parameters[1]
    S = parameters[2]
    xp = parameters[3]
    yp = parameters[4]
    l_ext = parameters[5]
    l_comp = parameters[6]

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

    return sumforce


def generate_feasible_individual(l_ratio_min):
    # Define range for inputs
    ds_min, ds_max = 12.0, 17.0
    hs_min, hs_max = -0.9, 0.75
    S_min, S_max = 240.0, 250.0
    xp_min, xp_max = -7.5, -4.75
    yp_min, yp_max = -0.85, 2.0
    l_ext_min, l_ext_max = 16.0, 35.5
    l_comp_min, l_comp_max = 17.63, 22.0
    
    # l_ratio_min = 0.55
    
    feasible = False
    
    while feasible == False:
        ds = np.random.uniform(ds_min, ds_max)
        hs = np.random.uniform(hs_min, hs_max)
        S = np.random.uniform(S_min, S_max)
        xp = np.random.uniform(xp_min, xp_max)
        yp = np.random.uniform(yp_min, yp_max)
        l_comp = np.random.uniform(l_comp_min, l_comp_max)
        if l_comp/l_ratio_min <= l_ext_max:
            l_ext_alt = l_comp/l_ratio_min
        else:
            l_ext_alt = l_ext_max
            
        l_ext = max(np.random.uniform(l_comp, l_ext_max), l_ext_alt)
        
        feasible = check_design_feasibility(ds, hs, S, xp, yp, l_ext, l_comp, l_ratio_min, False)
            
    parameters = [ds, hs, S, xp, yp, l_ext, l_comp]
    # print("Parameters: {}".format(parameters))
    
    return parameters


def check_design_feasibility(ds, hs, S, xp, yp, l_ext, l_comp, l_ratio_min, explain):
    feasible = True
    
    while True:
        # 1st check: l_ext > l_comp
        if l_ext <= l_comp:
            feasible = False
            break
        
        # 2nd check: Adequate l_ratio above l_ratio_min
        if l_comp/l_ext >= l_ratio_min:
            feasible = True
        else:
            # print("\nInfeasible compression ratio")
            feasible = False
            if explain == True:
                print("\nInfeasible compression ratio: {}".format(l_comp/l_ext))
            break
         
        # 3rd check: Minimum and maximum extensions based on position    
        wdoor = 5.75    # Length from the pivot point to the line of the door (the short part of the L shape made by the door)
    
        min_angle = 0.0
        max_angle = 90.0
        npts = 90
        theta = np.linspace(min_angle, max_angle, num=npts, endpoint=True)
        thetarad = theta*np.pi/180.0
        
        nu = np.arctan(hs/ds)
        sigma = thetarad - nu
        xs = np.sqrt(ds**2 + hs**2)*np.cos(sigma)
        ys = np.sqrt(ds**2 + hs**2)*np.sin(sigma)
    
        # Gas spring length
        ls = np.sqrt((np.abs(xp) + xs)**2 + (yp - ys)**2)
        maxspring = np.around(np.max(ls), 1)    # Maximum spring extension
        minspring = np.around(np.min(ls), 1)    # Minimum spring extension
        
        if maxspring > l_ext:
            feasible = False
            if explain == True:
                print("\nMax spring length exceeded by {}".format(maxspring-l_ext))
            break
        if minspring < l_comp:
            feasible = False
            if explain == True:
                print("\nMin spring length violated by {}".format(l_comp-minspring))
            break
        
        # # 4th check: Clearance under the frame corner
        # # Line equation for gas spring
        # m_closed = (ys[0] - yp) / (xs[0] - xp)
        # m_open = (ys[-1] - yp) / (xs[-1] - xp)
        
        # b_closed = ys[0] - m_closed*xs[0]
        # b_open = ys[-1] - m_open*xs[-1]
        
        # # Corner points that must not be intersected
        # x_corner_closed = 4.75
        # y_corner_closed = 2.25 - 0.25
        # x_corner_open = 2.25
        # y_corner_open = -4.75 - 0.25
        
        # # Check if the gas spring line is below the corner point
        # spring_closed_ht = m_closed*x_corner_closed + b_closed
        # spring_open_ht = m_open*x_corner_open + b_open
        
        # if spring_closed_ht > y_corner_closed:
        #     feasible = False
        # if spring_open_ht > y_corner_open:
        #     feasible = False
        
        # If the design has passed all feasibility checks, exit the while loop
        if feasible == True:
            break
        
    return feasible

# Tournament selection
def selection(pop, scores, k=3):
    # First random selection
    selection_ix = np.random.randint(len(pop))
    for ix in np.random.randint(0, len(pop), k-1):
        # Check if better (e.g. perform a tournament)
        if scores[ix] < scores[selection_ix]:
            selection_ix = ix
    return pop[selection_ix]


# Crossover two parents to create two children
def crossover(p1, p2, r_cross, r_mut):
    # Children are copies of parents by default
    c1, c2 = p1.copy(), p2.copy()
    
    # Check for recombination
    if np.random.rand() < r_cross:
        # select crossover point that is not on the end of the parameter list
        pt = np.random.randint(1, len(p1)-2)
        # perform crossover
        c1 = p1[:pt] + p2[pt:]
        c2 = p2[:pt] + p1[pt:]
        
    mutation(c1, r_mut)
    mutation(c2, r_mut)
    return [c1, c2]


# Mutation
def mutation(individual, r_mut):
    # Define range for inputs
    ds_min, ds_max = 12.0, 17.0
    hs_min, hs_max = -0.9, 0.75
    S_min, S_max = 240.0, 250.0
    xp_min, xp_max = -7.5, -4.75
    yp_min, yp_max = -0.85, 2.0
    l_ext_min, l_ext_max = 16.0, 35.5
    l_comp_min, l_comp_max = 17.63, 22.0
        
    for i in range(len(individual)):
        # Check for a mutation
        if np.random.rand() < r_mut:
            if i == 0:
                individual[i] = np.random.uniform(ds_min, ds_max)
            elif i == 1:
                individual[i] = np.random.uniform(hs_min, hs_max)
            elif i == 2:
                individual[i] = np.random.uniform(S_min, S_max)
            elif i == 3:
                individual[i] = np.random.uniform(xp_min, xp_max)
            elif i == 4:
                individual[i] = np.random.uniform(yp_min, yp_max)
            elif i == 5:
                individual[i] = np.random.uniform(l_ext_min, l_ext_max)
            elif i == 6:
                individual[i] = np.random.uniform(l_comp_min, l_comp_max)
    

def genetic_algorithm(objective, n_iter, n_pop, r_cross, r_mut, l_ratio_min):
    # Initial population
    pop = [generate_feasible_individual(l_ratio_min) for x in range(n_pop)]
    
    # Keep track of best solution
    best, best_eval = 0, objective(pop[0])
    
    best_evals = []
    # Enumerate generations
    for gen in range(n_iter):
        print("Generation {}".format(gen))
        # Evaluate all candidates in the population
        scores = [objective(c) for c in pop]
        # Check for new best solution
        for i in range(n_pop):
            if scores[i] < best_eval:
                best, best_eval = pop[i], scores[i]
                print("\nGeneration {}: new best {} = {}".format(gen, pop[i], scores[i]))
        best_evals.append(best_eval)
        # select parents
        selected = [selection(pop, scores) for x in range(n_pop)]
        
        # Create the next generation
        children = []
        for i in range(0, n_pop, 2):
            # Get selected parents in pairs
            p1, p2 = selected[i], selected[i+1]
            # crossover and mutation
            feasible = False
            loop_counter = 0
            while feasible == False:
                c1, c2 = crossover(p1, p2, r_cross, r_mut)
                c1_feasible = check_design_feasibility(c1[0], c1[1], c1[2], c1[3], c1[4], c1[5], c1[6], l_ratio_min, False)
                c2_feasible = check_design_feasibility(c2[0], c2[1], c2[2], c2[3], c2[4], c2[5], c1[6], l_ratio_min, False)
                if (c1_feasible or c2_feasible) is False:
                    feasible = False
                    loop_counter += 1
                    if loop_counter > 5000:
                        # print("Maximum checks reached")
                        c1 = generate_feasible_individual(l_ratio_min)
                        c2 = generate_feasible_individual(l_ratio_min)
                        break
                else:
                    feasible = True
                # print("Feasibility check: {}".format(loop_counter))
            
            children.append(c1)
            children.append(c2)
        
        # replace population
        pop = children
    return [best, best_eval], best_evals
    

if __name__ == '__main__':
    n_iter = 200
    n_pop = 200
    r_cross = 0.9
    r_mut = 0.3
    l_ratio_min = 0.5
    
    [best, score], best_evals = genetic_algorithm(objective_function, n_iter, n_pop, r_cross, r_mut, l_ratio_min)
    
    generations = list(range(n_iter))
    plt.plot(generations, best_evals)