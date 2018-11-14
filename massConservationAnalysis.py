# Mass Conservation Analysis

from initialConditions import *
from advectionSchemes import *
from diagnostics import *
from TimeStep_Error_Plotting import *
import numpy as np

def nmc_error(phi, phi_initial):
    "Calculates Normalised Mass Conservation error of phi given initial condition"
    "phi_initial"
    error = (phi.mean() - phi_initial.mean())/phi_initial.mean()
    return error

# Parameters
xmin = 0
xmax = 1
c = 0.4

# Calculate NMC error for different nt

nt_values = np.arange(1, 100, 1) # Array to story different number of time steps
error_Matrix = np.zeros((len(nt_values), 2)) # Matrix to store errors

# FTBS scheme
for i,nt in enumerate(nt_values):
    error_Matrix[i, 0] = nt # Store nt value
    nx = int(nt/c) # nx must be changed that we run the scheme for same amount of time t
    dx = (xmax - xmin)/nx
    x = np.arange(xmin, xmax, dx) # Spatial grid points
    phiOld = cosBell(x, 0, 0.75)
    phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 1] = nmc_error(phiFTBS, phiOld) # Store NMC errors
np.savetxt('results/mass_conservation/FTBS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Mass conservation error of FTBS scheme for different nt')
plt.plot(nt_values, error_Matrix[:, 1],'--g^', label = 'FTBS')

# FTCS scheme
for i,nt in enumerate(nt_values):
    error_Matrix[i, 0] = nt # Store nt value
    nx = int(nt/c) # nx must be changed that we run the scheme for same amount of time t
    dx = (xmax - xmin)/nx
    x = np.arange(xmin, xmax, dx) # Spatial grid points
    phiOld = cosBell(x, 0, 0.75)
    phiFTCS = FTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 1] = nmc_error(phiFTCS, phiOld) # Store NMC errors
np.savetxt('results/mass_conservation/FTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Mass conservation error of FTCS scheme for different nt')
plt.plot(nt_values, error_Matrix[:, 1],'--r^', label = 'FTCS')

# BTCS scheme
for i,nt in enumerate(nt_values):
    error_Matrix[i, 0] = nt # Store nt value
    nx = int(nt/c) # nx must be changed that we run the scheme for same amount of time t
    dx = (xmax - xmin)/nx
    x = np.arange(xmin, xmax, dx) # Spatial grid points
    phiOld = cosBell(x, 0, 0.75)
    phiBTCS = BTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 1] = nmc_error(phiBTCS, phiOld) # Store NMC errors
np.savetxt('results/mass_conservation/BTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Mass conservation error of BTCS scheme for different nt')
plt.plot(nt_values, error_Matrix[:, 1],'--c^', label = 'BTCS')

# CTCS scheme
for i,nt in enumerate(nt_values):
    error_Matrix[i, 0] = nt # Store nt value
    nx = int(nt/c) # nx must be changed that we run the scheme for same amount of time t
    dx = (xmax - xmin)/nx
    x = np.arange(xmin, xmax, dx) # Spatial grid points
    phiOld = cosBell(x, 0, 0.75)
    phiCTCS = CTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 1] = nmc_error(phiCTCS, phiOld) # Store NMC errors
np.savetxt('results/mass_conservation/CTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Mass conservation error of CTCS scheme for different nt')
plt.plot(nt_values, error_Matrix[:, 1],'--m^', label = 'CTCS')


plt.show()
