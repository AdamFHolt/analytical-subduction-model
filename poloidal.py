### Script that calculates stresses along the slab top and velocity field in the mantle wedge ###
############# We start by considering a non-curved slab and constant trench velocity ##############

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#################################
# set geometric variables
slab_thickness      = 100.    # km
init_slab_dip       = 40.     # deg
init_slab_depth     = 500     # km
op_thickness        = 50.     # km
ymax                = 660.    # km
xmax                = 3000.   # km
xsp                 = 1000.    # km, where does the slab starts
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
vt                  = -5.e-2    # m/yr
#################################

##### set up arrays
ys = np.linspace(0, init_slab_depth, 100) # slab chunks y coordinates
xs = ys/np.tan(np.radians(init_slab_dip)) # slab chunks x coordinates
ro = abs(ys/np.sin(np.radians(init_slab_dip))) # distance from trench
theta = np.linspace(0, init_slab_dip, 100) # angle from trench
vn = np.zeros(len(ys)) # normal velocity
vs = np.zeros(len(ys)) # shear velocity
s = np.zeros(len(ro)) # slab chunks depth
dsn = np.zeros(len(ys)) # normal stress gradient
ts = np.zeros(len(ys)) # shear stress on slab top
tt = np.zeros(len(ys)) # shear stress at theta = 0
vr = np.zeros((len(ys), len(theta))) # radial velocity
vtheta = np.zeros((len(ys), len(theta))) # theta velocity

s = ro.max() - ro
ds = ro.max()/len(ro)
# plt.plot(xs +xsp, ro*np.sin(np.radians(init_slab_dip)), 'k', label = "ro")
# plt.plot(xs+xsp, ys, 'r', linestyle = "--",  label = "ys")
# plt.plot(xs+xsp, s*np.sin(np.radians(init_slab_dip)), 'b', linestyle = "-.", label = "s")
# plt.ylim(ymax, 0)
# plt.gca().set_aspect('equal')
# plt.legend()
# plt.show()

# vn[:] = vt * np.sin(np.radians(init_slab_dip))        # m/yr
# vs[:] = -vt * np.cos(np.radians(init_slab_dip))        # m/yr

## initialize summation variable
sum = 0

for i in range(len(ys)):
        sum += vn[i]*s[i]
        dsn[i] = -((12*asthen_visc)/(ro[i]**3 * init_slab_dip**3))*(sum) + ((6*asthen_visc*(vs[i]+vt))/(ro[i]**2 * init_slab_dip**2))
        ts[i] = -(ro[i]*init_slab_dip*dsn[i]) + ((asthen_visc*(vt-vs[i]))/(ro[i]*init_slab_dip))
        # tt[i] = ((asthen_visc*(vt-vs[i]))/(ro[i]*init_slab_dip)) + ((ro[i]*init_slab_dip)/2)*dsn[i]
        # for j in range(len(theta)):
        #     vr[i] = -(ro[i]**2 * init_slab_dip**2)/2 * dsn[i] + tt[i]*ro[i]*theta[j] - vt

plt.plot(dsn/1e9, ys, 'k', label = 'Normal stress gradient')
plt.plot(ts/1e9, ys, 'r', label = 'Shear stress on slab top')
plt.ylim(ymax, 0)
# plt.xlim(-0.1, 1e7)
plt.ylabel('Depth (km)')
plt.legend()
plt.show()


##### plot slab position
plt.plot(xs, ys, 'k')
plt.xlim(0, xmax)
plt.ylim(ymax, 0)
plt.gca().set_aspect('equal')
plt.savefig('poloidal_slab_position.png')
plt.close()



