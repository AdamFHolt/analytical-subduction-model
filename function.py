#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt


# from functions import initialize_slab_geometry, evolve_slab_geometry
# from functions import plot_slab_shape, solve_slab_velocities
###

def create_slab_geometry(slab_thickness, radius_curvature, slab_depth, slab_dip, x_range, flat_range=500):
	"""
	Create an array with the geometry of a curved slab with specified dip, radius of curvature, slab depth, and thickness.
	X is horizontal distance, y is vertical depth.
	"""
	# Convert slab dip to radians
	slab_dip_rad = np.deg2rad(slab_dip)

	# Create x values for the flat and curved parts
	x_values_flat = np.linspace(0, flat_range, num=500)
	x_values_curved = np.linspace(flat_range, x_range + flat_range, num=1000)

	# Calculate y values based on the geometry of the slab
	y_values_flat = np.zeros_like(x_values_flat)
	y_values_curved = (radius_curvature - np.sqrt(radius_curvature**2 - ((x_values_curved - flat_range) - slab_thickness / np.tan(slab_dip_rad))**2))

	# Calculate the normal to the curve
	dx = np.gradient(x_values_curved)
	dy = np.gradient(y_values_curved)
	normals = np.vstack((-dy, dx))
	normals /= np.linalg.norm(normals, axis=0)

	# Create the additional slab profile at a constant profile-normal distance
	x_values_curved_deeper = x_values_curved + slab_thickness * normals[0]
	y_values_curved_deeper = y_values_curved + slab_thickness * normals[1]

	# For the flat part, the normal points vertically downwards
	x_values_flat_deeper = x_values_flat
	y_values_flat_deeper = y_values_flat + slab_thickness

	# Combine x and y values into a 2D array
	x_values = np.concatenate((x_values_flat, x_values_curved, x_values_curved_deeper[::-1], x_values_flat_deeper[::-1], x_values_flat[0:1]))
	y_values = np.concatenate((y_values_flat, y_values_curved, y_values_curved_deeper[::-1], y_values_flat_deeper[::-1], y_values_flat[0:1]))
	slab_coords = np.vstack((x_values, y_values))

	return slab_coords


slab_coords = create_slab_geometry(100., 250., 500., 70., 1000.)

def plot_slab(slab_coords):

	plt.figure()
	plt.plot(slab_coords[0], slab_coords[1])
	plt.xlabel('Horizontal Distance')
	plt.ylabel('Depth')
	plt.title('Slab Geometry')
	plt.gca().invert_yaxis()  # Invert y axis as depth increases downwards
	plt.show()

plot_slab(slab_coords)





# def solve_slab_velocities(...):
# 	...
# 	return(???)

# def evolve_slab_geometry(slab_coords, slab_vels, dt_myrs):
# 	slab_coords=1. # temp
# 	return(slab_coords)


# def plot_slab_shape(slab_coords, plot_filename):
# 	pass




