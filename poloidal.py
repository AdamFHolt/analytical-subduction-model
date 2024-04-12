### Script that calculates stresses along the slab top and velocity field in the mantle wedge ###
############# We start by considering a non-curved slab and constant trench velocity ##############

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#################################
# set geometric variables
slab_thickness      = 100.e3    # m
init_slab_dip       = 40.     # deg
init_slab_depth     = 500.e3     # m
op_thickness        = 50.e3    # m
ymax                = 660.e3    # m
xmax                = 3000.e3   # m
xsp                 = 1000.e3    # m, where does the slab starts
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
# trench velocity and far_field pressure
vt_yr                  = -5.e-2    # m/yr
yr                     = 3.15e7    # s
vt                     = vt_yr/yr # m/s
P                      = 70.e6      # Pa`
#################################

##### set up arrays
ys = np.linspace(0, init_slab_depth, 100) # slab chunks y coordinates
xs = ys/np.tan(np.radians(init_slab_dip))+xsp # slab chunks x coordinates
ro = abs(ys/np.sin(np.radians(init_slab_dip))) # distance from trench
theta = np.linspace(0, init_slab_dip, 100) # angle from trench
vn = np.zeros(len(ys)) # normal velocity
vs = np.zeros(len(ys)) # shear velocity
s = np.zeros(len(ro)) # slab chunks depth
dsn = np.zeros(len(ys)) # normal stress gradient
ts = np.zeros(len(ys)) # shear stress on slab top
sigma_n = np.zeros(len(ys)) # normal stress

s = ro.max() - ro
ds = ro.max()/len(ro)
# plt.plot(xs +xsp, ro*np.sin(np.radians(init_slab_dip)), 'k', label = "ro")
# plt.plot(xs+xsp, ys, 'r', linestyle = "--",  label = "ys")
# plt.plot(xs+xsp, s*np.sin(np.radians(init_slab_dip)), 'b', linestyle = "-.", label = "s")
# plt.ylim(ymax, 0)
# plt.gca().set_aspect('equal')
# plt.legend()
# plt.show()

vn[:] = vt * np.sin(np.radians(init_slab_dip))        # m/yr
vs[:] = -vt * np.cos(np.radians(init_slab_dip))        # m/yr

## initialize summation variable
sum = 0

for i in range(len(ys)):
        sum += vn[i]*s[i]
        dsn[i] = -((12*asthen_visc)/(ro[i]**3 * init_slab_dip**3))*(sum) + ((6*asthen_visc*(vs[i]+vt))/(ro[i]**2 * init_slab_dip**2))
        ts[i] = -(ro[i]*init_slab_dip*dsn[i]) + ((asthen_visc*(vt-vs[i]))/(ro[i]*init_slab_dip))

dsn = np.nan_to_num(dsn)
integ = 0


for j in range(len(dsn)):
        integ = integ + dsn[j]*s[j]
        sigma_n[j] = P + integ

# plt.plot(dsn/1e6, ys, 'k', label = 'Normal stress gradient')
plt.plot(ts/1e6, ys/1e3, 'r', label = 'Shear stress on slab top')
plt.plot(sigma_n/1e6, ys/1e3, 'b', label = 'Normal stress')
plt.ylim(ymax/1e3, 0)
# plt.xlim(0, 6)
plt.ylabel('Depth (km)')
plt.legend()
plt.xlabel('Stress (MPa)')
plt.savefig('slab_top_shear_stress.png', dpi = 1000)
plt.close()



##### plot slab position
plt.plot(xs/1e3, ys/1e3, 'k')
plt.xlim(0, xmax/1e3)
plt.ylim(ymax/1e3, 0)
plt.gca().set_aspect('equal')
plt.savefig('poloidal_slab_position.png', dpi = 1000)
plt.close()



