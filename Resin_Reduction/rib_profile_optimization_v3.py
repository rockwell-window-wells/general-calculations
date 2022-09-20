# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 08:12:02 2022

@author: Ryan.Larson
"""

"""
Goal: Find a rib profile that gives equivalent or greater stiffness with less
area (meaning less resin is required to make it)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import minimize


#### REDOING THIS TO USE THE ACTUAL RIB PROFILE WITH TWO FEET AND A TRAPEZOID

##### Geometric Calculation Functions for Hollow Trapezoid #####
def I_rectangle(b,h):
    Ix = (b*h**3)/12
    Iy = (h*b**3)/12
    A = b*h
    return Ix, Iy, A

def I_triangle(b,h,a):
    Ix = (b*h**3)/36
    Iy = (1/36)*(h*b**3 - h*a*b**2 + b*h*a**2)
    A = b*h/2
    return Ix, Iy, A

def parallel_axis_theorem(Ix, A, d):
    Ix_prime = Ix + A*d**2
    return Ix_prime

def triangle_centroid_height(h):
    return h/3

def rectangle_centroid_height(h):
    return h/2

def rect_shift(t,b_t,h):
    theta = np.arctan(b_t/h)
    B = np.pi/2 - theta
    rect_shift = t*(np.tan(B/2))
    return rect_shift

def I_open_trapezoid_about_base(ht,b1,b2,b3,t):
    # Angle offset c
    theta = np.arctan(2*ht/(b1-b2))
    phi = (np.pi - theta)/2
    c = t/np.tan(phi)
    
    h = ht + t
    
    # Rectangle 1
    Ix_rect1, _, A_rect1 = I_rectangle(b2, ht)
    d_rect1 = rectangle_centroid_height(ht)
    
    # Rectangle 2
    br2 = b3 - c
    Ix_rect2, _, A_rect2 = I_rectangle(br2, ht)
    d_rect2 = rectangle_centroid_height(ht) + t
    
    # Triangle
    b_t = (b1-b2)/2
    a = b_t     # Right triangle (a is the distance along the x axis from the left base point to the upper point)
    Ix_tri, _, A_tri = I_triangle(b_t, ht, a)
    d_tri_low = triangle_centroid_height(ht)
    d_tri_high = h - d_tri_low
    
    # Parallel axis theorem and combining
    Ix_rect1_prime = parallel_axis_theorem(Ix_rect1, A_rect1, d_rect1)
    Ix_rect2_prime = parallel_axis_theorem(Ix_rect2, A_rect2, d_rect2)
    Ix_tri1_prime = parallel_axis_theorem(Ix_tri, A_tri, d_tri_low)
    Ix_tri2_prime = parallel_axis_theorem(Ix_tri, A_tri, d_tri_high)
    
    w = 2*b3 + b1    
    Ix_solid_rectangle, _, A_solid_rectangle = I_rectangle(w, h)
    Ix_prime_solid_rectangle = parallel_axis_theorem(Ix_solid_rectangle, A_solid_rectangle, h/2)
    
    # Subtract smaller rectangles and triangles from larger rectangle
    Ix_prime_trapezoid = Ix_prime_solid_rectangle - Ix_rect1_prime - 2*(Ix_rect2_prime + Ix_tri1_prime + Ix_tri2_prime)
    
    A_hollow_trapezoid = A_solid_rectangle - A_rect1 - 2*(A_rect2 + 2*A_tri)
    
    return Ix_prime_trapezoid, A_hollow_trapezoid
    

##### Genetic Algorithm Functions #####
def objfunc(x):
    """
    Given parameter values, calculate the fitness of the design, with the goal
    to minimize A/Ix

    Parameters
    ----------
    parameters : list
        [h, b1, b2, t]

    Returns
    -------
    None.

    """
    # Read in parameter values from parameter list
    ht = x[0]
    b1 = x[1]
    b2 = x[2]
    b3 = x[3]
    t = x[4]
    
    Ix, A = I_open_trapezoid_about_base(ht, b1, b2, b3, t)

    return A
    
def constraint2(x):
    Ix, _ = I_open_trapezoid_about_base(x[0],x[1],x[2],x[3],x[4])
    return Ix - 0.34

if __name__ == "__main__":
    # Initial values
    ht0 = 1.075
    b10 = 4.248
    b20 = 2.56994615
    b30 = 0.911
    t0 = 0.075
    x0 = [ht0, b10, b20, b30, t0]
    Ix_ref, A_ref = I_open_trapezoid_about_base(x0[0], x0[1], x0[2], x0[3], x0[4])
    
    # Lower and upper bounds
    ht_lb = 0.1
    b1_lb = 2.0
    b2_lb = 0.01
    b3_lb = 0.911
    t_lb = 0.075
    
    ht_ub = 1.5
    b1_ub = 6.0
    b2_ub = 5.0
    b3_ub = 5.0
    t_ub = 0.1    
    
    res = minimize(objfunc,
                   x0,
                   constraints = (
                       {'type': 'ineq', 'fun': lambda x: x[1] - 2*x[2]},
                       {'type': 'ineq', 'fun': constraint2}
                           ),
                   bounds = ((ht_lb,ht_ub), (b1_lb,b1_ub), (b2_lb,b2_ub), (b3_lb,b3_ub), (t_lb,t_ub)),
                   method='SLSQP')
    
    print(res)
    Ix_res, A_res = I_open_trapezoid_about_base(res.x[0], res.x[1], res.x[2], res.x[3], res.x[4])
    
    A_pct = A_res/A_ref
    Ix_pct = Ix_res/Ix_ref
    
    print("\nArea % of original:\t{}".format(A_pct))
    print("Ix % of original:\t{}".format(Ix_pct))
    
    
    ### Plot the design ###
    # Helpful variables
    w = 2*res.x[3] + res.x[1]
    w_ref = 2*b30 + b10
    theta = np.arctan(2*res.x[0]/(res.x[1]-res.x[2]))
    phi = (np.pi - theta)/2
    c = res.x[4]/np.tan(phi)
    h = res.x[0] + res.x[4]
    
    theta = np.arctan(2*ht0/(b10-b20))
    phi = (np.pi - theta)/2
    c_ref = t0/np.tan(phi)
    h_ref = ht0 + t0
    
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(211,aspect='equal',autoscale_on=True)
    plt.xlim([-2*w/3,2*w/3])
    plt.ylim([0, h+0.25])
    x = [-w/2, -w/2, -res.x[1]/2 - c, -res.x[2]/2 - c, res.x[2]/2 + c, res.x[1]/2 + c, w/2, w/2, res.x[1]/2, res.x[2]/2, -res.x[2]/2, -res.x[1]/2]
    y = [0, res.x[4], res.x[4], h, h, res.x[4], res.x[4], 0, 0, res.x[0], res.x[0], 0]
    ax.add_patch(patches.Polygon(xy=list(zip(x,y)), fill=True))
    ax.title.set_text("Optimized Design: {} Area, {} Stiffness".format(np.around(A_pct,2), np.around(Ix_pct,2)))    
    
    ax = fig.add_subplot(212,aspect='equal',autoscale_on=True)
    plt.xlim([-2*w_ref/3,2*w_ref/3])
    plt.ylim([0, h_ref+0.25])
    x = [-w_ref/2, -w_ref/2, -b10/2 - c_ref, -b20/2 - c_ref, b20/2 + c_ref, b10/2 + c_ref, w_ref/2, w_ref/2, b10/2, b20/2, -b20/2, -b10/2]
    y = [0, t0, t0, h_ref, h_ref, t0, t0, 0, 0, ht0, ht0, 0]
    ax.add_patch(patches.Polygon(xy=list(zip(x,y)), fill=True))
    ax.title.set_text("Current Design")
    
    plt.show()