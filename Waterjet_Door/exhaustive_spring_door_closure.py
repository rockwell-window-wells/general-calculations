import numpy as np
import matplotlib.pyplot as plt

########## Torsion spring option ##########
npts_var = 1000
optforce = None

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

##### Design parameters #####
nsprings      = 4                  # Number of springs at each angle that will engage
torque_rating_ub = 10000.0     # Spring torque rating (lbf*in)
torque_rating_lb = 1.0
torque_rating_vec = np.linspace(torque_rating_lb, torque_rating_ub, npts_var, endpoint=True)

###### Redoing the dead angle concept ######
force = np.zeros(theta.shape)

for k in range(npts_var):
    kappa = torque_rating_vec[k]/theta_max_rad
    for i in range(len(theta)):
        angle = theta[i]
        anglerad = theta[i]*np.pi/180.0
        force[i] = (1/L_door)*(weight*l_c_door*np.sin(anglerad) - kappa*anglerad*nsprings)

    # Minimize the range between the min and max forces
    minabsforce = np.abs(min(force))
    maxabsforce = np.abs(max(force))
    sumforce = minabsforce + maxabsforce

    # Minimize the total absolute value of force experienced (less useful in this case)
    # absforce = np.abs(force)
    # sumforce = sum(absforce)
    #
    if optforce is None:
        optforce = sumforce
    elif sumforce < optforce:
        optforce = sumforce
        torque_rating = torque_rating_vec[k]

print("Optimal design:")
print("Torque rating: {} lbf*in".format(torque_rating))
