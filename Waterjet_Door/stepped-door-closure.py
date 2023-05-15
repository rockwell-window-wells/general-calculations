# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 14:26:15 2021

@author: Ryan.Larson
"""

import numpy as np
import matplotlib.pyplot as plt

########## Torsion spring option ##########
# Initialize angles in degrees for open position
min_open = 1.0
max_open = 90.0
npts = 90
theta = np.linspace(min_open, max_open, num=npts, endpoint=True)

# Door parameters
weight = 250.0                          # Door weight, lbs (estimate 250)
L_door = 45.0                           # Height from pivot point to top of door
cm_factor = 2/3     # Factor for where center of mass of door is relative to pivot
l_c_door = cm_factor*L_door # Moment arm for center of mass of door

# Calculate parameters of existing spring, using a few assumptions
theta_max = 90.0                        # Assumed angle at which maximum torque is achieved
theta_max_rad = theta_max*np.pi/180.0   # Convert the maximum angle
tr0 = 2000.0                  # Spring torque rating, lbf*in
k0 = tr0/theta_max_rad        # Torsional spring constant, (lbf*in/rad)
print("\nExisting spring constant: {} lbf*in/rad\n".format(k0))

##### Design parameters #####
theta_engage  = [0.0]          # Angles at which a new number of springs start to engage (deg)
nsprings      = [1]                  # Number of springs at each angle that will engage
torque_rating = [7800.0]     # Spring torque rating (lbf*in)

# Make checks and make intermediate calculations
assert (len(theta_engage) == len(nsprings)) and (len(theta_engage) == len(torque_rating))
tot_springs = sum(nsprings)
theta_engage = np.array(theta_engage)
nsprings = np.array(nsprings)
torque_rating = np.array(torque_rating)
kappa = torque_rating/theta_max_rad

###### Redoing the dead angle concept ######
force = np.zeros(theta.shape)
theta_engage_rad = theta_engage*np.pi/180.0

for i in range(len(theta)):
    angle = theta[i]
    anglerad = theta[i]*np.pi/180.0
    force[i] = (1/L_door)*(weight*l_c_door*np.sin(anglerad))
    
    for j in range(len(theta_engage)):
        if angle >= theta_engage[j]:
            force[i] -= (1/L_door)*(kappa[j]*(anglerad-theta_engage_rad[j])*nsprings[j])
    

# Plot the results 
fig1 = plt.figure(1)
plt.plot(theta,force)
plt.xlabel("Angle (deg)")
plt.ylabel("Force (lbf)")

titlestring = "{} Total Springs: ".format(tot_springs)
for i in range(len(theta_engage)):
    addstring = "{} {}s at {} deg, ".format(nsprings[i], int(torque_rating[i]), theta_engage[i])
    titlestring += addstring
titlestring = titlestring[:-2]
plt.title(titlestring)

maxstring = "Max push force: {}".format(np.around(max(force)),1)
maxind = np.argmax(force)
maxpt = (30.0, force[maxind]/2)
print(maxpt)

minstring = "Max pull force: {}".format(np.around(min(force)),1)
minind = np.argmin(force)
minpt = (30.0, force[maxind]/3)

plt.annotate(maxstring, maxpt)
if min(force) < 0:
    plt.annotate(minstring, minpt)
    
fig1.suptitle("Concept 1: Torsion Springs Only")
plt.show()

print("Max push force: {}".format(max(force)))
if min(force) < 0:
    print("Max pull force: {}".format(min(force)))