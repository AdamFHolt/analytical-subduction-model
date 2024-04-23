import numpy as np
import time, sys, os
from scipy.linalg import solve
from numba import jit

import matplotlib.pyplot as plt
T1=time.time()

# Toroidal flow 

# set-up basic parameters
up_mantle_thick = 660.e3 # km--> m
litho_thick = 100.e3 # km--> m
lam = up_mantle_thick - litho_thick # Hele-shaw flow channel thickness -- sublithosphere thickness
Vt = 0 # cm/yr
Vm = 0  # cm/yr 
Vr = 3. # cm/yr
cmyr_to_ms = 0.01 / (365 * 24 * 60 * 60) # Conversion factor from cm/yr to m/s
Vt, Vm, Vr = Vt * cmyr_to_ms, Vm * cmyr_to_ms, Vr * cmyr_to_ms
#print(Vr)
viscosity = 1.e22  # N*s / m2
a = 500.e3 # km--> m  half-width of slab
print(up_mantle_thick,litho_thick,lam,Vt, Vm, Vr,viscosity,a)


# set-up grids along the slab 
xi = 0
#yi = np.linspace(-a, a,  int(2*a/2e3)+1, endpoint=True)
yi = np.linspace(0, a,  int(a/2e3)+1, endpoint=True)


# set-up grids for whole domain
x = np.linspace(-4*a, 4*a, int(8*a/10e3)+1)
y = np.linspace(-4*a, 4*a, int(8*a/10e3)+1)      




# solve matrix equation [B][A]=[C]
matrix_C = np.zeros((len(yi),1))
matrix_C[:] = (Vr - ( Vt + Vm )/2 )  * 12 * viscosity / (lam**2)


matrix_B = np.zeros((len(yi), len(yi)))  #(251,251) # consider the points on the slab, xi=0, yi= (0, a), x = 0, y=(0,a)

for i in range(matrix_B.shape[0]): 
    for j in range(matrix_B.shape[1]): 
        if i==j: 
            matrix_B[i, j] = 0
        else:
            matrix_B[i, j] =     -  1 /  (yi[j] - yi[i])**2 
       

matrix_A = solve(matrix_B,matrix_C) 

extended_matrix_A = np.concatenate([matrix_A[::-1][:-1], matrix_A], axis=0) # Ai from 0-a to -a to a 

array_A = extended_matrix_A.reshape(len(extended_matrix_A[:,0]))  # transform from list to array, to avoid format warning from python

yi = np.linspace(-a, a,  int(2*a/2e3)+1, endpoint=True)
print(array_A.shape,yi.shape)
# calculate the Pressure
@jit(nopython=True)
def calculate_pressure(x, y, yi, array_A,xi):
    Pressure = np.zeros((len(y), len(x)))   # len(y) is the number of rows (height of the array, corresponding to the y-axis), and len(x) is the number of columns (width of the array, corresponding to the x-axis).
    
    for p in range(len(y)):
        for q in range(len(x)):
            pre_p= np.zeros(len(array_A))

            for i in range(len(array_A)):
                denominator =  (x[q] - xi) ** 2 + (y[p] - yi[i]) ** 2
                if denominator == 0:
                    pass
                else:
                    pre_p[i] = array_A[i] * (x[q] - xi) / denominator
            Pressure[p, q] = np.sum(pre_p[:])
            
    return Pressure

Pressure = calculate_pressure(x, y, yi, array_A,xi)
print(Pressure.shape)

for p in range(len(y)):
        for q in range(len(x)):
            if x[q] == a and y[p] == 0:
                    test=Pressure[p, q]
print(test,"p")

#normalized pressure
max_pressure = np.max(Pressure)
normalized_pressure = Pressure /1e6 # / max_pressure

import math
#expo= (-0.9*300e3/a)**0.8

test_p = 12*viscosity*a/(lam**2) * (Vr-(Vt+Vm)/2) * math.exp(-0.9)
print(test_p,"test_p")





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



plt.figure(figsize=(6,6))
plt.scatter([xi]*len(yi), yi, color='r', label='Slab points')
plt.xlim([x.min(), x.max()])
plt.ylim([y.min(), y.max()])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Slab Points on the X-Y Domain')
plt.savefig("slab.png")
