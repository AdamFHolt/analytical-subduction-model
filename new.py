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
Vt = 0  # cm/yr
Vm = 0  # cm/yr 
Vr = 3. # cm/yr
cmyr_to_ms = 1. / (365.25 * 24 * 60 * 60)  # Conversion factor from cm/yr to m/s
Vt, Vm, Vr = Vt * cmyr_to_ms, Vm * cmyr_to_ms, Vr * cmyr_to_ms
viscosity = 1.e20  # N*s / m2
a = 500e3 # km--> m  half-width of slab

# set-up grids along the slab 
xi = 0
yi = np.linspace(-a, a,  int(2*a/2e3)+1, endpoint=True)

# set-up grids for whole domain
x = np.linspace(-4*a, 4*a, int(8*a/10e3)+1)
y = np.linspace(-4*a, 4*a, int(8*a/10e3)+1)      


# solve matrix equation [B][A]=[C]
matrix_C = np.zeros_like(yi)
matrix_C[:] = (Vr - ( Vt + Vm )/2 )  * 12 * viscosity / (lam **2)
matrix_C = matrix_C.reshape(-1,1)  # transform from list to array, to avoid format warning from python


matrix_B = np.zeros((len(yi), len(yi)))
for i in range(len(yi)):
    for j in range(len(yi)):
        matrix_B[i, j] =   - ((yi[i] - yi[j])**2) / (yi[i] - yi[j])**2

matrix_B = np.nan_to_num(matrix_B)  #  there is situation that denominator is 0, replace nan with 0 to avoid error message from following solve function


matrix_A = solve(matrix_B,matrix_C)
array_A = matrix_A.reshape(len(matrix_A[:,0]))  # transform from list to array, to avoid format warning from python
print(array_A.shape)

# calculate the Pressure
@jit(nopython=True)
def calculate_pressure(x, y, yi, array_A,xi):
    Pressure = np.zeros((len(x), len(y)))
    for p in range(len(x)):
       # print(p)
        for q in range(len(y)):
            for i in range(len(yi)):
                denominator = (x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2
                if denominator != 0:
                    Pressure[p, q] += array_A[i] * (x[p] - xi) / denominator 
    return Pressure

Pressure = calculate_pressure(x, y, yi, array_A,xi)

T2=time.time()

# calculate the Vx and Vy
@jit(nopython=True)
def calculate_vx(x, y, yi, array_A,xi):
    vx = np.zeros((len(x), len(y)))
    for p in range(len(x)):
        for q in range(len(y)):
            for i in range(len(yi)):
                denominator = (x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2
                if denominator != 0:
                    vx[p, q] +=   array_A[i] * (lam **2)/ (12 * viscosity) * ((x[p] - xi) ** 2 - (y[q] - yi[i]) ** 2) / ( (x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2)
    return vx

vx = calculate_vx(x, y, yi, array_A,xi)

#print("vx",vx)

@jit(nopython=True)
def calculate_vy(x, y, yi, array_A,xi):
    vy = np.zeros((len(x), len(y)))
    for p in range(len(x)):
        for q in range(len(y)):
            for i in range(len(yi)):
                denominator = (x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2
                if denominator != 0:
                    vy[p, q] += (Vt+Vm)/2  - (lam **2) / (12 * viscosity)   *  (-2 * array_A[i] * (x[p] - xi) * (y[q] - yi[i]) )  /   ( (x[p] - xi) ** 2 + (y[q] - yi[i]) ** 2) 
    return vy

vy = calculate_vy(x, y, yi, array_A,xi)
#print(vy)


import matplotlib.pyplot as plt

#normalized pressure
max_pressure = np.max(Pressure)
normalized_pressure = Pressure / max_pressure

# Save the normalized pressure field as an image
plt.imshow(normalized_pressure, cmap='viridis', extent=[x.min(), x.max(), y.min(), y.max()])


# Plot velocity vectors
plt.colorbar(label='Normalized Pressure')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Normalized Pressure Field (Maximum Pressure = 1)')
plt.savefig("normalized_pressure.png")



plt.quiver(x[::50], y[::50], vy[::50,::50] ,vx[::50,::50])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Velocity Vectors')
plt.savefig("velocity_vectors.png")

plt.close()



normalized_pressure_along_y0 = normalized_pressure[:, np.abs(y).argmin()]
# Plot normalized pressure data along the line y=0
plt.plot(x, normalized_pressure_along_y0)
plt.xlabel('X')
plt.ylabel('Normalized Pressure')
plt.title('Normalized Pressure along y=0')
plt.savefig("Normalized_Pressure_along_y=0.png")
plt.close()

Pressure_in_paper = np.zeros_like(x)
Pressure_in_paper[:] =  12 * viscosity * a /(lam**2) * (Vr -(Vt+Vm)/2)  * np.exp(-0.9 * (x[:]/a)**0.8)

Pressure_in_paper = np.nan_to_num(Pressure_in_paper) 
max_pressure = np.max(Pressure_in_paper)
normalized_Pressure_in_paper= Pressure_in_paper / max_pressure
plt.plot(x, normalized_Pressure_in_paper)
#print(normalized_Pressure_in_paper)
plt.xlabel('X')
plt.ylabel('normalized_Pressure_in_paper')
plt.title('normalized_Pressure_in_paper')
plt.savefig("normalized_Pressure_in_paper.png")
plt.close()


Pressure = Pressure[:, np.abs(y).argmin()]
plt.plot(x, Pressure_in_paper, label='Pressure in Paper')
plt.plot(x, Pressure, label='Pressure')
plt.xlabel('X')
plt.ylabel('Normalized Pressure')
plt.title('Comparison of Normalized Pressure')
plt.legend()
plt.savefig("normalized_Pressure_comparison.png")
plt.close()
