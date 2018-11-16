# Mass Conservation Analysis

from initialConditions import *
from advectionSchemes import *
from diagnostics import *
import numpy as np

def main(xmin, xmax, c):

    # Range of nt and nx values ensuring c = 0.5

    nt_values = np.arange(20, 80, 4) # Array to story different number of time steps
    nx_values = np.arange(10, 40, 2) # Adjust nx to fix courant number

    error_Matrix = np.zeros(len(nt_values)) # Matrix to store errors

    # Normalised Mass Conservation error calculations

    # FTBS scheme
    for i,nt in enumerate(nt_values):
        dx = (xmax - xmin)/nx_values[i] # Calculate new dx
        x = np.arange(xmin, xmax, dx) # Spatial grid points
        phiOld = cosBell(x, 0, 0.75) # Initial condition
        phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[i] = nmc_error(phiFTBS, phiOld) # Store NMC errors
    plt.plot(nt_values, error_Matrix,'--g.', label = 'FTBS')

    # FTCS scheme
    for i,nt in enumerate(nt_values):
        dx = (xmax - xmin)/nx_values[i] # Calculate new dx
        x = np.arange(xmin, xmax, dx) # Spatial grid points
        phiOld = cosBell(x, 0, 0.75) # Initial condition
        phiFTCS = FTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[i] = nmc_error(phiFTCS, phiOld) # Store NMC errors
    plt.plot(nt_values, error_Matrix,'--r.', label = 'FTCS')

    # BTCS scheme
    for i,nt in enumerate(nt_values):
        dx = (xmax - xmin)/nx_values[i] # Calculate new dx
        x = np.arange(xmin, xmax, dx) # Spatial grid points
        phiOld = cosBell(x, 0, 0.75) # Initial condition
        phiBTCS = BTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[i] = nmc_error(phiBTCS, phiOld) # Store NMC errors
    plt.plot(nt_values, error_Matrix,'--c.', label = 'BTCS')

    # CTCS scheme
    for i,nt in enumerate(nt_values):
        dx = (xmax - xmin)/nx_values[i] # Calculate new dx
        x = np.arange(xmin, xmax, dx) # Spatial grid points
        phiOld = cosBell(x, 0, 0.75) # Initial condition
        phiCTCS = CTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[i] = nmc_error(phiCTCS, phiOld) # Store NMC errors
    plt.plot(nt_values, error_Matrix,'--m.', label = 'CTCS')

    # Saving Results for Conservative Schemes
    plt.xlabel('Number of time steps nt')
    plt.ylabel('Mass Conservation Error')
    plt.title('Mass Conservation Error for c = %s' % c)
    plt.legend()
    save_Location = 'plots/analysis/conservative_schemes_%s.png' % c
    plt.savefig(save_Location, dpi=1000)

    # Semi Lagrangian scheme
    nt_values = nt_values[nt_values <= 100] # Use less nt_values as SL scheme is slow
    error_Matrix = np.zeros(len(nt_values)) # Matrix to store errors

    for i,nt in enumerate(nt_values):
        dx = (xmax - xmin)/nx_values[i] # Calculate new dx
        x = np.arange(xmin, xmax, dx) # Spatial grid points
        phiOld = cosBell(x, 0, 0.75) # Initial condition
        phi_sem_lag = sem_lag(phiOld.copy(), c, nt, x, dx) # Run the scheme for each (c, nt) pair
        error_Matrix[i] = nmc_error(phi_sem_lag, phiOld) # Store NMC errors
    plt.plot(nt_values, error_Matrix,'--k.', label = 'SL')

    #Plot Details
    plt.xlabel('Number of time steps nt')
    plt.ylabel('Mass Conservation Error')
    plt.title('Mass Conservation Error for c = %s' % c)
    plt.legend()
    save_Location = 'plots/analysis/mass_conservation_%s.png' % c
    plt.savefig(save_Location, dpi=1000)
    plt.close()
    return
