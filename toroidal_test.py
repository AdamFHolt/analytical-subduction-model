import numpy as np
import time, sys, os
from scipy.linalg import solve
import time
from numba import jit

T1=time.time()
#################################
# set geometric variables
slab_thickness      = 100.    # km
init_slab_dip       = 70.     # deg
init_slab_depth     = 500     # km
init_rad_curvature  = 250.    # km
slab_halfwidth      = 1000.   # km
op_thickness        = 50.     # km
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

# Toroidal flow -Tao

# set-up basic parameters
up_mantle_thick = 660e3 # km--> m
litho_thick = 100e3 # km--> m
lam = up_mantle_thick - litho_thick # Hele-shaw flow channel thickness -- sublithosphere thickness
Vt = 3.e-2  # cm/yr --> m/yr # upper plate velocity
Vm = 0.e-2  # cm/yr --> m/yr # lower mantle velocity
Vr = 3.e-2  # cm/yr--> m/yr  # trench velocity

conversion_factor = 1 / (365.25 * 24 * 3600)  # years to seconds

# Convert m/yr to m/s
Vt  = Vt  * conversion_factor
Vr  = Vr  * conversion_factor
Vm  = Vm  * conversion_factor



viscosity = 1.e20  # N*s / m2
a = 350e3 # km--> m  half-width of slab

# set-up grids along the slab 
xi = 0
yi = np.linspace(0, a, 351, endpoint=True)
#print("yi-shape",yi.shape)


# set-up grids for whole domain
x = np.linspace(-4*a, 4*a, int(8*a/1000)+1)
y = np.linspace(0, 4*a, int(8*a/1000)+1)        #).reshape(-1,1)
#print("x-shape",x.shape,"y-shape", y.shape)



# solve matrix equation [B][A]=[C]
matrix_C = np.zeros_like(yi)
matrix_C[:] = (Vr - ( Vt + Vm )/2 )  * 12 * viscosity / (lam **2)
matrix_C = matrix_C.reshape(-1,1)
#print(matrix_C.shape, matrix_C)


matrix_B = np.zeros((len(yi), len(yi)))
for i in range(len(yi)):
    for j in range(len(yi)):
        matrix_B[i, j] =   - ((yi[i] - yi[j])**2) / (yi[i] - yi[j])**2
#print(matrix_B, matrix_B.shape)
matrix_B = np.nan_to_num(matrix_B)



matrix_A = solve(matrix_B,matrix_C)
#print(matrix_A.shape,matrix_A)

array_A = matrix_A.reshape(len(matrix_A[:,0]))  # transform from list to array, to avoid format warning from python

# calculate the Pressure
Pressure = np.zeros((len(x), len(y)))

#print(Pressure.shape,Pressure)



#print(array_A)
@jit(nopython=True)
def calculate_pressure(x, y, yi, array_A,xi):
    Pressure = np.zeros((len(x), len(y)))
    for p in range(len(x)):
        for q in range(len(y)):
            for i in range(len(yi)):
                denominator = (x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2
                if denominator != 0:
                    Pressure[p, q] += array_A[i] * (x[p] - xi) / ((x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2)
    return Pressure

Pressure = calculate_pressure(x, y, yi, array_A,xi)

#print(Pressure.shape, Pressure)
"""

for p in range(len(x)):
    for q in range(len(y)):
        for i in range(len(yi)):
           # matrix_A_i = float(array_A[i])
           
            Pressure[p, q] += array_A[i] * (x[p] - xi) / ((x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2)

#Pressure = calculate_pressure(x, y, yi, array_A,xi)



"""

T2=time.time()
#print(T2-T1)

#
# plot the pressure field
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x, y)
plt.figure()
plt.contourf(X, Y, Pressure, cmap='viridis')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Pressure')
plt.show()

plt.savefig("toroidal.png")
print(Pressure)
