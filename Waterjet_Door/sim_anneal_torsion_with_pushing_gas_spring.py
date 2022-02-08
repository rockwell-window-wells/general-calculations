"""Optimization by simulated annealing of the combined torsion spring and
pushing gas spring concept for the waterjet door.
"""
import numpy as np
import matplotlib.pyplot as plt

def objective(parameters):
    """Given parameter values, calculate the force curve of the design.
    """
    # Read in parameter values from parameter list
    ds = parameters[0]
    S = parameters[1]
    xp = parameters[2]
    yp = parameters[3]

    ### Fixed parameters ###
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

    ### Torsion spring parameters ###
    theta_max_rad = np.pi/2.0
    nsprings = 2                # Number of torsion springs
    torque_rating = 2000.0      # Torsion spring torque rating (in*lbs) (4 1877.0 rated springs balances the door at 90 degrees)
    kappa_tor = torque_rating/theta_max_rad

    ### Derived parameters ###
    # S:
    xs = ds*np.cos(thetarad)
    ys = ds*np.sin(thetarad)

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

    # Angles
    gamma = np.arctan(yw/xw)
    kappa = np.arctan(yf/xf)
    psi = np.pi/2 - gamma
    phi = np.arctan((yp-ys)/(np.abs(xp)+xs))
    beta = thetarad - phi

    ### Force and extension ###
    force = np.zeros(thetarad.shape)
    for i in range(len(thetarad)):
        angle = thetarad[i]
        force[i] = (W*Lw*np.sin(psi[i]) - ns*S*ds*np.sin(beta[i])  - kappa_tor*angle*nsprings) / Lf

    # Minimize the range between the min and max forces
    minabsforce = np.abs(min(force))
    maxabsforce = np.abs(max(force))
    sumforce = minabsforce + maxabsforce

    return sumforce

def Boltzmann(dE, dEavg, T):
    P = np.exp(-dE/(dEavg*T))
    return P

def SimAnneal(ds0, S0, xp0, yp0):
    # Starting design
    xs = np.array([ds0, S0, xp0, yp0])
    xub = np.array([5.75, 500.0, -1.0, -1.0])
    xlb = np.array([1.0, 5.0, -40.0, -23.5])
    if xs[0] > xub[0]:
        raise ValueError("ds0 out of bounds")
    elif xs[1] > xub[1]:
        raise ValueError("S0 out of bounds")
    elif xs[2] > xub[2]:
        raise ValueError("xp0 out of bounds")
    elif xs[3] > xub[3]:
        raise ValueError("yp0 out of bounds")
    elif xs[0] < xlb[0]:
        raise ValueError("ds0 out of bounds")
    elif xs[1] < xlb[1]:
        raise ValueError("S0 out of bounds")
    elif xs[2] < xlb[2]:
        raise ValueError("xp0 out of bounds")
    elif xs[3] < xlb[3]:
        raise ValueError("yp0 out of bounds")
    fs = objective(xs)
    xsearch = [xs]

    # Select Ps, Pf, N, and calculate Ts, Tf, and F
    Ps = 0.3                # Probability of acceptance at start
    Pf = 0.00001             # Probability of acceptance at finish
    N = 200                 # Number of cycles

    Ts = -1/np.log(Ps)      # Temperature at start
    Tf = -1/np.log(Pf)      # Temperature at finish
    F = (Tf/Ts)**(1/(N-1))  # Temperature reduction factor each cycle

    # Perturbation information
    delta = 2.0               # Max perturbation
    n = 2                   # Starting number of perturbations per cycle

    # Holding variables
    dE = 0.0
    dEavg = 0.0
    perturbations = list(range(N))
    objvals = [0] * N

    # Set starting values
    xc = xs
    fc = fs
    T = Ts

    # Step through the cycles
    for i, perturb in enumerate(perturbations):
        # Add the current objective value to the objective vector for plotting
        objvals[i] = objective(xc)

        # Step through the perturbations
        for j in range(n):
            # Perturb xc by some random value within delta. If any perturbed
            # value falls outside the specified bounds, retry the perturbation.
            while True:
                dsp = np.random.uniform(-delta, delta)
                Sp = np.random.uniform(-delta, delta)
                xpp = np.random.uniform(-delta, delta)
                ypp = np.random.uniform(-delta, delta)
                perturb = np.array([dsp, Sp, xpp, ypp])
                xp = xc + perturb

                if xp[0] > xub[0]:
                    continue
                elif xp[1] > xub[1]:
                    continue
                elif xp[2] > xub[2]:
                    continue
                elif xp[3] > xub[3]:
                    continue
                elif xp[0] < xlb[0]:
                    continue
                elif xp[1] < xlb[1]:
                    continue
                elif xp[2] < xlb[2]:
                    continue
                elif xp[3] < xlb[3]:
                    continue
                else:
                    break

            # print(xp)

            # Get the objective value at the perturbed point
            fp = objective(xp)

            # Calculate values for Boltzmann function in case they're needed
            dE = np.abs(fp - fc)
            if i == 1 and j == 1:
                dEavg = dE
            else:
                dEavg = (dEavg + dE)/2

            P = Boltzmann(dE, dEavg, T)

            # Check if the new design is better than the old design
            if fp < fc:
                xc = xp     # Accept as current design if better
                fc = objective(xc)
            else:
                # If the new design is worse, generate a random number and
                # compare to the Boltzmann probability. If the random number is
                # lower than the Boltzmann probability, accept the worse design
                # as the current design
                randnum = np.random.uniform(0,1)

                if randnum < P:
                    xc = xp
                    fc = objective(xc)

        # Decrease the temperature by factor F
        T = F*T

        # Increase the number of perturbations every few cycles
        if (i % 3) == 1:
            n += 1

        # Save the new search position at the end of each cycle
        xsearch.append(xc)

    return perturbations, objvals, xsearch

if __name__ == '__main__':
    perturbations, objvals, xsearch = SimAnneal(2.0, 100.2, -1.3, -8.3)

    ds_end = xsearch[-1][0]
    S_end = xsearch[-1][1]
    xp_end = xsearch[-1][2]
    yp_end = xsearch[-1][3]

    best_ind = objvals.index(min(objvals))
    ds_best = xsearch[best_ind][0]
    S_best = xsearch[best_ind][1]
    xp_best = xsearch[best_ind][2]
    yp_best = xsearch[best_ind][3]

    # Plot the cooling curve
    fig1 = plt.figure(1, figsize=(12,8))
    plt.plot(perturbations, objvals)
    plt.xlabel("Cycles")
    plt.ylabel("Objective")
    plt.title("Cooling Curve - Torsion Spring and Pushing Gas Spring")
    end_annotation = "Final design:\nObj: {}\nds: {}\nS: {}\nxp: {}\nyp: {}".format(objvals[-1],ds_end, S_end, xp_end, yp_end)
    plt.annotate(end_annotation, (perturbations[-1], objvals[-1]), (0.8*perturbations[-1], 0.8*np.max(objvals)), arrowprops=dict(facecolor='black', shrink = 0.01, width=0.5))
    best_annotation = "Best design:\nObj: {}\nds: {}\nS: {}\nxp: {}\nyp: {}".format(objvals[best_ind],ds_best, S_best, xp_best, yp_best)
    mid_ind = int(len(perturbations)/3)
    plt.annotate(best_annotation, (perturbations[best_ind], objvals[best_ind]), (0, 0.25*(np.max(objvals)-np.min(objvals))+np.min(objvals)), arrowprops=dict(facecolor='black', shrink = 0.01, width=0.5))
    plt.show()
