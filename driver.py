#!/usr/bin/python3

# from functions import initialize_slab_geometry, evolve_slab_geometry
# from functions import plot_slab_shape, solve_slab_velocities

import numpy as np
import time, sys, os

#################################
# set geometric variables
slab_thickness      = 100.    # km
init_slab_dip       = 70.     # deg
init_slab_depth     = 500     # km
init_rad_curvature  = 250.    # km
slab_halfwidth      = 1000.   # km
op_thickness        = 50.     # km
lower_mantle_velocity  =  0   # m/s
#################################
# set mechanical variables
slab_density        = 3400.    # kg/m3    
crust_density       = 2800.    # kg/m3
prism_density       = 2800.    # kg/m3
asthen_density      = 3300.    # kg/m3
asthen_visc         = 2.5e20   # Pa s
lith_rigidity       = 10**19   # N m
g                   = 9.81     # m/s2
################################
# time-stepping variables
nt                  = 100      # num. timesteps
dt_myrs             = 0.25     # Myrs
#################################

##### set up initial slab shape
# slab_coords = initialize_slab_geometry(slab_thickness,init_rad_curvature,init_slab_depth,init_slab_dip)
plot_filename = 'plots/slab-shape.t%s.pdf' % 0
# plot_slab_shape(slab_coords, plot_filename)

for i in range(nt):

    time_myrs = dt_myrs + (i * dt_myrs)
    print("---------Solving at %.2f Myrs---------" % time_myrs)

    ##### solve
    # slab_vels = solve_slab_velocities(... lots of variables ...)
    
    ##### evolve slab geometry using the velocity solution
    # slab_coords = evolve_slab_geometry(slab_coords, slab_vels, dt_myrs)

    ##### plot new slab geometry
    plot_filename = 'plots/slab-shape.t%s.pdf' % i
    # plot_slab_shape(slab_coords, plot_filename)
