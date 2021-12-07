"""door-closure-calc.py: Evaluate force profiles for waterjet door closure
options.

Date Created: 6 Dec 2021
Author: Ryan Larson

"""

import numpy as np
import matplotlib.pyplot as plt

########## Torsion spring option ##########
# Initialize angles in degrees for open position
min_open = 0.0
max_open = 90.0
npts = 91
theta_open = np.linspace(min_open, max_open, num=npts, endpoint=True)

# print(theta_open)

# Set parameters
theta_dead = 0.0   # Angle at which the springs start to engage
theta_max = 90.0    # Assumed angle at which maximum torque is achieved
theta_max_rad = theta_max*np.pi/180.0
torque_rating = 2000.0    # Spring torque rating, lbf*in
kappa = torque_rating/theta_max_rad     # Torsional spring constant, (lbf*in/deg)
weight = 600.0      # Door weight, lbs
height = 45.0       # Height from pivot point to top of door
nsprings = 4

force = np.zeros(theta_open.shape)

for i in range(len(theta_open)):
    angle = theta_open[i]
    if angle < theta_dead:
        force[i] = (weight/2.0)*np.sin(angle*np.pi/180.0)
    else:
        force[i] = (weight/2.0)*np.sin(angle*np.pi/180.0) - nsprings*(kappa*(angle-theta_dead)*np.pi)/(180.0*height)

print("\nMaximum force, Concept 1: {} lbs".format(max(force)))

fig1 = plt.figure(1)
plt.plot(theta_open,force)
plt.xlabel("Angle (deg)")
plt.ylabel("Force (lbf)")
titlestring = "{} Torsion Springs, Dead Angle {} Degrees".format(nsprings, theta_dead)
plt.title(titlestring)
fig1.suptitle("Concept 1: Torsion Springs Only")
# plt.show()




########## Pulley and cable option ##########
P0 = 1.0    # Preload exerted by rotational spring on the door, lbf
delta_p = 5.0                   # Pulley height above the door
h_pulley = height + delta_p     # Pulley actual height
L = height                      # Renaming door height to make equations short
weight_rating = 150.0           # Weight rating of garage door spring
nrotations = 8                  # Expected number of rotations (approx. # feet of garage door)
theta_max_rad = 2*np.pi*nrotations
R = 2.0     # Radius of rotational spring drum, in
kappa = weight_rating*R / theta_max_rad
# torque_rating = 9.5   # Spring torque rating, lbf*in
# kappa = torque_rating/theta_max_rad     # lbf*in/rad
nsprings = 2

force = np.zeros(theta_open.shape)

for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    c = np.sqrt(h_pulley**2 + L**2 -2*h_pulley*L*np.cos(angle))
    omega = c/R
    P = P0 + nsprings*kappa*c/R
    force[i] = (weight/2.0)*np.sin(angle) - P*np.sin(np.arcsin(h_pulley*np.sin(angle)/c))

print("\nMaximum force, Concept 2: {} lbs".format(max(force)))

fig2 = plt.figure(2)
plt.plot(theta_open,force)
plt.xlabel("Angle (deg)")
plt.ylabel("Force (lbf)")
titlestring = "{} Torsion Springs".format(nsprings)
plt.title(titlestring)
fig2.suptitle("Concept 2: Cables and Pulley")
plt.show()
