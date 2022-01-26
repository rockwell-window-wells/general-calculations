"""door-closure-calc.py: Evaluate force profiles for waterjet door closure
options.

Date Created: 6 Dec 2021
Author: Ryan Larson

"""

import numpy as np
import matplotlib.pyplot as plt

########## Torsion spring option ##########
# Initialize angles in degrees for open position
min_open = 1.0
max_open = 90.0
npts = 90
theta_open = np.linspace(min_open, max_open, num=npts, endpoint=True)

# print(theta_open)

# Set parameters
theta_dead = 15.0   # Angle at which the springs start to engage
theta_max = 90.0    # Assumed angle at which maximum torque is achieved
theta_max_rad = theta_max*np.pi/180.0
torque_rating = 2000.0    # Spring torque rating, lbf*in
kappa = torque_rating/theta_max_rad     # Torsional spring constant, (lbf*in/rad)
weight = 250.0      # Door weight, lbs (estimate 250)
height = 45.0       # Height from pivot point to top of door
nsprings_initial = 2
nsprings_dead = 3

force = np.zeros(theta_open.shape)

for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    if angle < theta_dead:
        force[i] = (weight/2.0)*np.sin(angle*np.pi/180.0) - nsprings_initial*(kappa*(angle-theta_dead)*np.pi)/(180.0*height)
    else:
        force[i] = (weight/2.0)*np.sin(angle*np.pi/180.0) - (nsprings_initial+nsprings_dead)*(kappa*(angle-theta_dead)*np.pi)/(180.0*height)

# print("\nMaximum force, Concept 1: {} lbs".format(max(force)))

# fig1 = plt.figure(1)
# plt.plot(theta_open,force)
# plt.xlabel("Angle (deg)")
# plt.ylabel("Force (lbf)")
# titlestring = "{} Total Springs, {} Initial, {} Added, Dead Angle {} Degrees".format(nsprings_initial + nsprings_dead, nsprings_initial, nsprings_dead, theta_dead)
# plt.title(titlestring)
# fig1.suptitle("Concept 1: Torsion Springs Only")
# plt.show()

########## Pulley and cable option ##########
P0 = 1.0    # Preload exerted by rotational spring on the door, lbf
delta_p = -22.5                   # Pulley height relative to top of the door
h_pulley = height + delta_p     # Pulley actual height
cable_offset = 11.25
L_cable = height - cable_offset     # Cable moment arm
L_door = height     # Door height
cm_factor = 2/3     # Factor for where center of mass of door is relative to pivot
l_c_door = cm_factor*L_door # Moment arm for center of mass of door
R = 2.625     # Radius of rotational spring drum, inches
nsprings = 2


# Adjust the original concept using some fixed values:
force = np.zeros(theta_open.shape)

for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    force[i] = (1/L_door)*(weight*l_c_door*np.sin(angle) - kappa*angle*nsprings)
    
fig1 = plt.figure(1)
plt.plot(theta_open,force)
plt.xlabel("Angle (deg)")
plt.ylabel("Force (lbf)")
titlestring = "{} Springs at Pivot Point".format(nsprings)
plt.title(titlestring)
fig1.suptitle("Concept 1: Torsion Springs Only")

###### Redoing the dead angle concept ######
force = np.zeros(theta_open.shape)

for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    if angle < theta_dead*np.pi/180.0:
        force[i] = (1/L_door)*(weight*l_c_door*np.sin(angle) - kappa*angle*nsprings_initial)
    else:
        force[i] = (1/L_door)*(weight*l_c_door*np.sin(angle) - kappa*angle*nsprings_initial - kappa*(angle-theta_dead*np.pi/180.0)*nsprings_dead)
        
fig1 = plt.figure(1)
plt.plot(theta_open,force)
plt.xlabel("Angle (deg)")
plt.ylabel("Force (lbf)")
titlestring = "{} Total Springs, {} Initial, {} Added, Dead Angle {} Degrees".format(nsprings_initial + nsprings_dead, nsprings_initial, nsprings_dead, theta_dead)
plt.title(titlestring)
fig1.suptitle("Concept 1: Torsion Springs Only")
plt.show()



# Determine the required spring constant to balance the door using two springs and the original concept
kappa_necessary = np.zeros(theta_open.shape)
for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    kappa_necessary[i] = (weight * l_c_door * np.sin(angle))/(nsprings*angle)

# fig2 = plt.figure(2)
# plt.plot(theta_open, kappa_necessary)   # Not sure this is calculated correctly, since the plot is not what I would expect
# plt.xlabel("Angle (deg)")
# plt.ylabel("Spring Constant (lbf*in)")
# titlestring = "Concept 1: Required Spring Constant, {} Springs".format(nsprings)
# plt.title(titlestring)


# Plot the required spring constant per spring over the range of opening angles, cable concept
kappa_necessary = np.zeros(theta_open.shape)
for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    c = np.sqrt(h_pulley**2 + l_c_door**2 - 2*h_pulley*l_c_door*np.cos(angle))
    psi = np.arcsin((h_pulley*np.sin(angle)) / c)
    omega = c/R

    P = (weight * l_c_door * np.sin(angle))/(L_cable * np.sin(psi))
    kappa_necessary[i] = (P * R**2)/(nsprings * c)      # Necessary spring constant per spring, lbf*in



# Plot kappa_necessary over the range of angles for the cable concept
# fig3 = plt.figure(3)
# plt.plot(theta_open,kappa_necessary)
# plt.xlabel("Angle (deg)")
# plt.ylabel("Spring Constant (lbf*in)")
# titlestring = "Cable Concept Required Spring Constant, {} Springs".format(nsprings)
# plt.title(titlestring)


# Calculate the force experienced by the operator using kappa_necessary springs
force = np.zeros(theta_open.shape)
for i in range(len(theta_open)):
    angle = theta_open[i]*np.pi/180.0
    c = np.sqrt(h_pulley**2 + l_c_door**2 - 2*h_pulley*l_c_door*np.cos(angle))
    psi = np.arcsin((h_pulley*np.sin(angle)) / c)
    # kappa = max(kappa_necessary)
    kappa = 40.0
    P = P0 + nsprings*c*kappa / R**2
    force[i] = (1/L_door) * (weight*l_c_door*np.sin(angle) - P*L_cable*np.sin(psi))

print(force)
print("\nMaximum force, Concept 2: {} lbs".format(np.abs(max(force))))

# fig4 = plt.figure(4)
# plt.plot(theta_open,force)
# plt.xlabel("Angle (deg)")
# plt.ylabel("Operator Force (lbf)")
# titlestring = "{} Torsion Springs, {} in-lbs each".format(nsprings, kappa)
# plt.title(titlestring)
# fig4.suptitle("Concept 2: Cables and Pulley")
# plt.show()
