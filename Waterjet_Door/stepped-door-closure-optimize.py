import numpy as np

from gekko import GEKKO

m = GEKKO()

# Declare model parameters STILL WORKING ON FIXING THIS!!!!!!!!
theta = m.Array(m.Param, (90))
for i in range(len(theta)):
    if i == 0:
        theta[i].value = 1.0
    else:
        theta[i].value  = theta[i-1].value + 1.0   # Not sure if this results in endpoint inclusive

force = m.Array(m.Param, (90))
for f in force:
    f.value = 0.0

weight = m.Const(value=300.0)
L_door = m.Const(value=45.0)
cm_factor = m.Const(value=2/3)
l_c_door = cm_factor*L_door

theta_max = m.Const(value=90.0)
theta_max_rad = theta_max*np.pi/180.0


# Declare variables and initial values
n_engage = m.Var(2, lb=1, ub=10, integer=True)

theta_engage = m.Array(m.Var, (n_engage.value))
theta_engage[0].value = 0.0
theta_engage[0].lower = 0.0
theta_engage[0].upper = 85.0
theta_engage[1].value = 15.0
theta_engage[1].lower = 0.0
theta_engage[1].upper = 85.0

theta_engage_rad = m.Array(m.Var, (n_engage,1))
for i in range(n_engage):
    theta_engage_rad[i].value = theta_engage[i].value*np.pi/180.0
    theta_engage_rad[i].lower = theta_engage[i].lower*np.pi/180.0
    theta_engage_rad[i].upper = theta_engage[i].upper*np.pi/180.0

nsprings = m.Array(m.Var, (n_engage))
nsprings[0].value = 2
nsprings[0].lower = 1
nsprings[0].upper = 10
nsprings[1].value = 2
nsprings[1].lower = 1
nsprings[1].upper = 10

totsprings = m.Param(value=sum(nsprings))

torque_rating = m.Array(m.Var, (n_engage))
torque_rating[0].value = 2000.0
torque_rating[0].lower = 10.0
torque_rating[0].upper = 5000.0
torque_rating[1].value = 2000.0
torque_rating[1].lower = 10.0
torque_rating[1].upper = 5000.0

kappa = m.Array(m.Var, (n_engage))
for i in range(n_engage):
    kappa[i].value = torque_rating[i].value/theta_max_rad
    kappa[i].lower = torque_rating[i].lower/theta_max_rad
    kappa[i].upper = torque_rating[i].upper/theta_max_rad

# Calculate force
for i in range(len(theta)):
    angle = theta[i]
    anglerad = theta[i]*np.pi/180.0
    force[i] = (1/L_door)*(weight*l_c_door*np.sin(anglerad))

    for j in range(len(theta_engage)):
        if angle >= theta_engage[j]:
            force[i] -= (1/L_door)*(kappa[j]*(anglerad-theta_engage_rad[j])*nsprings[j])

# Objective
sumforce = sum(np.abs(force))

# Minimize objective
m.Minimize(sumforce)

# Solve optimization
m.solve()

# Plot the results
fig1 = plt.figure(1)
plt.plot(theta,force)
plt.xlabel("Angle (deg)")
plt.ylabel("Force (lbf)")

titlestring = "{} Total Springs: ".format(tot_springs)
for i in range(len(theta_engage)):
    addstring = "{} at {} deg, ".format(nsprings[i], theta_engage[i])
    titlestring += addstring
titlestring = titlestring[:-2]
plt.title(titlestring)
fig1.suptitle("Concept 1: Torsion Springs Only")
plt.show()
