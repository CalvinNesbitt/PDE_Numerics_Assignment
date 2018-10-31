"This is a workspace to write the semi Lagranian scheme"

### The matplotlib package contains plotting functions              ###
import matplotlib.pyplot as plt
import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *
from TimeStep_Error_Plotting import *

# Parameters
xmin = 0
xmax = 1
nx = 5
nt = 3
c = 0.2

# Derived parameters
dx = (xmax - xmin)/nx

# Spatial points for plotting and for defining initial conditions
x = np.arange(xmin, xmax, dx)
print(x)

phiOld = cosBell(x, 0, 0.75)

#def sem_lag(phiOld, c, nt):
"Linear advection of profile in phiOld using Semi Lagranian, Courant number c"
"for nt time-steps"

nx = len(phiOld)

# new time-step array for phi
#phi = phiOld.copy()

# FTCS for each time-step
for it in range(nt):
    # Loop through all space using remainder after division (%)
    # to cope with periodic boundary conditions
    for j in range(nx):
        k = int(np.floor(abs(j + 1 - c))) # Index defining interpolating polynomial
        # abs term is to deal with j = 0 index
        print(k)
        base_points = x[k-1:len(x)]# Array of the spatial values we will interpolate from
        print(base_points)
        #phi[j] = phiOld[j] - 0.5*c*\
                 #(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])

        # update arrays for next time-step
        #phiOld = phi.copy()

    #return phi
