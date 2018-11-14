# Analysis of the order of convergence of each scheme. Each scheme is run for
# a variety of spatial discretisation. The l2 error is calculated for each
# discretisation.

import numpy as np
from advectionSchemes import *
from initialConditions import *
from diagnostics import *

# Parameters
xmin = 0
xmax = 1
c = 0.5

# Values to test over
nt_values = np.arange(20, 80, 4) # Array to story different number of time steps
nx_values = np.arange(10, 40, 2) # Adjust nx to fix courant number
dx_values = (xmax - xmin)/nx_values # Calculate spatial step

# Error calcualtion for each scheme

error_Matrix = np.zeros(len(dx_values)) # Matrix to store errors

# FTBS scheme
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = nt_values[i] # nt must be changed that we run the scheme for same amount of time t
    phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
    error_Matrix[i] = l2ErrorNorm(phiFTBS, phiAnalytic) # Store L2 errors
plt.plot(dx_values, error_Matrix,'--g^', label = 'FTBS')

# BTCS scheme
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = nt_values[i] # nt must be changed that we run the scheme for same amount of time t
    phiBTCS = BTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
    error_Matrix[i] = l2ErrorNorm(phiBTCS, phiAnalytic) # Store L2 errors
plt.plot(dx_values, error_Matrix,'--r^', label = 'BTCS')

# CTCS scheme
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = nt_values[i] # nt must be changed that we run the scheme for same amount of time t
    phiCTCS = CTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
    error_Matrix[i] = l2ErrorNorm(phiCTCS, phiAnalytic) # Store L2 errors
plt.loglog(dx_values, error_Matrix,'--m^', label = 'CTCS')

# SL scheme
for i,dx in enumerate(dx_values):
    if nt < 60: # Stopping scheme running for too many timesteps
        x = np.arange(xmin, xmax, dx)
        phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
        nt = nt_values[i] # nt must be changed that we run the scheme for same amount of time t
        phi_sem_lag = sem_lag(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[i] = l2ErrorNorm(phi_sem_lag, phiAnalytic) # Store L2 errors
plt.loglog(dx_values, error_Matrix,'--k^', label = 'SL')

#Plot Details
plt.plot(dx_values, dx_values, label = 'x')
plt.loglog(dx_values, dx_values**2, label = 'x^2')
plt.xlabel('Spatial step (dx)')
plt.ylabel('L2 Error (Logarithmic)')
plt.title('L2 Error as function of dx')
plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=5)
plt.savefig('plots/L2Error_vs_dx', bbox_inches='tight')
plt.show()
