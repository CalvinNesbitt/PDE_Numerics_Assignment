"This is a workspace to write the semi Lagranian scheme"

import matplotlib.pyplot as plt
import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *
from TimeStep_Error_Plotting import *
from scipy.interpolate import lagrange

# Parameters
xmin = 0
xmax = 1
nx = 100
nt = 20
c = 0.2

# Derived parameters
dx = (xmax - xmin)/nx

# Spatial points for plotting and for defining initial conditions
x = np.arange(xmin, xmax, dx)
#print(x)

phiOld = cosBell(x, 0, 0.75)
plt.plot(x, phiOld, label='Initial')

#def sem_lag(phiOld, c, nt):
"Linear advection of profile in phiOld using Semi Lagranian, Courant number c"
"for nt time-steps"

nx = len(phiOld)

# new time-step array for phi
phi = phiOld.copy()

# Semi Lagrangian for each time-step
for it in range(nt):

    for j in range(nx):
        # Find base_points of interpolating polynomial
        k = int(np.floor(j - c)) # Index below advected point
        #print(k)
        x_base_points = np.array([x[k-1],x[k],x[(k+1)%nx],x[(k+2)%nx]])
        y_base_points = np.array([phiOld[k-1],phiOld[k],phiOld[(k+1)%nx],phiOld[(k+2)%nx]])

        # Build interpolating polynomial
        poly = lagrange(x_base_points, y_base_points)
        # Find down-wind point
        beta = j - c - k
        x_jd = beta * dx + x[k]
        #print(x_jd)
        # Calculate phi upwind
        phi[j] = poly(x_jd)
    phiOld = phi.copy() # Update for next time step

plt.plot(x, phi, label='SL', color='red')
plt.show()

    #return phi
