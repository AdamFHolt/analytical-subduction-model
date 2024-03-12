import numpy as np
import time, sys, os
from scipy.linalg import solve
from numba import jit

T1=time.time()

# Toroidal flow -Tao

# set-up basic parameters
up_mantle_thick = 660e3 # km--> m
litho_thick = 100e3 # km--> m
lam = up_mantle_thick - litho_thick # Hele-shaw flow channel thickness -- sublithosphere thickness
Vt = 0 # cm/yr
Vm = 0  # cm/yr 
Vr = 3. # cm/yr
cmyr_to_ms = 1. / (365.25 * 24 * 60 * 60)  # Conversion factor from cm/yr to m/s
Vt, Vm, Vr = Vt * cmyr_to_ms, Vm * cmyr_to_ms, Vr * cmyr_to_ms
#print(Vr)
viscosity = 1.e20  # N*s / m2
a = 500e3 # km--> m  half-width of slab

# set-up grids along the slab 
xi = 0
yi = np.linspace(-a, a,  int(2*a/2e3)+1, endpoint=True)

# set-up grids for whole domain
x = np.linspace(-4*a, 4*a, int(8*a/10e3)+1)
y = np.linspace(-4*a, 4*a, int(8*a/10e3)+1)      

import matplotlib.pyplot as plt
plt.figure(figsize=(6,6))
plt.scatter([xi]*len(yi), yi, color='r', label='Slab points')
plt.xlim([x.min(), x.max()])
plt.ylim([y.min(), y.max()])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Slab Points on the X-Y Domain')
#plt.legend()
#plt.grid(True)
plt.savefig("slab.png")


# solve matrix equation [B][A]=[C]
matrix_C = np.zeros_like(yi)
matrix_C[:] = (Vr - ( Vt + Vm )/2 )  * 12 * viscosity / (lam **2)
matrix_C = matrix_C.reshape(-1,1)  # transform from list to array, to avoid format warning from python
print(matrix_C.shape)



matrix_B = np.zeros((len(yi), len(yi)))   # consider the points on the slab, xi=0, yi= (-a, a), x = 0, y=(-a,a)
matrix_B[:] = np.nan
for i in range(len(yi)):
    for j in range(len(yi)):
        if i==j: 
            matrix_B[i, j] = 0
        else:
            matrix_B[i, j] =   - 1 / (yi[i] - yi[j])**2

#print(matrix_B.shape, matrix_B)

matrix_A = solve(matrix_B,matrix_C)
print(matrix_A.shape)
array_A = matrix_A.reshape(len(matrix_A[:,0]))  # transform from list to array, to avoid format warning from python
print(array_A.shape)

# calculate the Pressure
@jit(nopython=True)
def calculate_pressure(x, y, yi, array_A,xi):
    Pressure = np.zeros((len(y), len(x)))   # len(y) is the number of rows (height of the array, corresponding to the y-axis), and len(x) is the number of columns (width of the array, corresponding to the x-axis).
    for p in range(len(y)):
        for q in range(len(x)):
            for i in range(len(yi)):
                denominator = (x[q] - xi) ** 2 + (y[p] - yi[i]) ** 2
                if denominator != 0:
                    Pressure[p, q] += array_A[i] * (x[q] - xi) / denominator 
    return Pressure

Pressure = calculate_pressure(x, y, yi, array_A,xi)
print(Pressure.shape)

#normalized pressure
max_pressure = np.max(Pressure)
normalized_pressure = Pressure / max_pressure

# Save the normalized pressure field as an image
plt.imshow(normalized_pressure, cmap='viridis', extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar(label='Normalized Pressure')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Normalized Pressure Field')
plt.savefig("normalized_pressure.png")
T2=time.time()





# calculate the Vx and Vy
@jit(nopython=True)
def calculate_vx(x, y, yi, array_A,xi):
    vx = np.zeros((len(y), len(x)))   # len(y) is the number of rows (height of the array, corresponding to the y-axis), and len(x) is the number of columns (width of the array, corresponding to the x-axis).
    for p in range(len(y)):
        for q in range(len(x)):
            for i in range(len(yi)):
                denominator =   (x[q] - xi) ** 2 + (y[p] - yi[i]) ** 2
                coefficient = (lam **2)/ (12 * viscosity)
                if x[q] ==xi and -a<= y[p] <= a :
                    vx[p, q] = 0
                else:
                    vx[p, q] += ( (Vt+Vm)/2 + array_A[i] * coefficient * ( (x[q] - xi) ** 2 - (y[p] - yi[i]) ** 2) / (denominator **2) ) /Vr
                
                    
                
    return vx

vx = calculate_vx(x, y, yi, array_A,xi)
print(np.max(vx),np.min(vx))

print("vx",vx.shape)

@jit(nopython=True)
def calculate_vy(x, y, yi, array_A,xi):
    vy = np.zeros((len(y), len(x)))
    for p in range(len(y)):
        for q in range(len(x)):
            for i in range(len(yi)):
                coefficient = (lam **2)/ (6 * viscosity)
                denominator =  (x[q] - xi) ** 2 + (y[p] - yi[i]) ** 2
                if x[q] ==xi and -a<= y[p] <= a :
                    vy[p, q] = 0
                else:
                    vy[p, q] += ( coefficient   *  ( array_A[i] * (x[q] - xi) * (y[p] - yi[i]) )  /   ( denominator **2 )   )/Vr
    return vy

vy = calculate_vy(x, y, yi, array_A,xi)
print(vy.shape)
print(np.max(vy),np.min(vy))
# Save the normalized pressure field as an image



# Plot velocity vectors
plt.quiver(x[::10], y[::10], vx[::10,::10], vy[::10,::10],scale=0.1)
plt.xlabel('X')
plt.ylabel('Y')

plt.title('Velocity Vectors')
plt.savefig("velocity_vectors.png")
