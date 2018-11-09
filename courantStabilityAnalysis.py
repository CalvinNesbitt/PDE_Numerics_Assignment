# Stability analysis for a varying courant number c. Requires directories called
# 'plots' and 'results'. For a variety of courant numbers the L1, L2 and Linf
# errors from each scheme in advectionSchemes is calculated.
# The results are saved in a .csv file. Plots are saved as .png files in plots.
from initialConditions import *
from advectionSchemes import *
from diagnostics import *
from TimeStep_Error_Plotting import *
import numpy as np

# Parameters
def main(xmin, xmax, nx, dx, x, phiOld):

    # Error calculation for FTBS scheme
    c_values = np.arange(0.1, 1.2, 0.1)
    error_Matrix = np.zeros((len(c_values), 4)) # Matrix to store errors
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phiFTBS, phiAnalytic) # Store L2 errors
    np.savetxt('results/FTBS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1],color = 'g', label = 'FTBS')

    # Error calculation for FTCS scheme
    c_values = np.arange(0.1, 1.2, 0.1)
    error_Matrix = np.zeros((len(c_values), 4)) # Matrix to store errors
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phiFTCS = FTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phiFTCS, phiAnalytic) # Store L2 errors
    np.savetxt('results/FTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTCS scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1],color = 'r', label = 'FTCS')

    # Error calculation for FTBS scheme
    c_values = [0.2, 1.2]
    error_Matrix = np.zeros((len(c_values), 4)) # Matrix to store errors
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phi_sem_lag = sem_lag(phiOld.copy(), c, nt, x, dx)
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phi_sem_lag, phiAnalytic) # Store L2 errors
    np.savetxt('results/Sem_Lag_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running Semi Lagrangian scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1], color = 'b', label = 'Semi Lagrangian')
    plt.legend(bbox_to_anchor=(0.5, 0.5))
    plt.show()
        #phiBTCS = BTCS(phiOld.copy(), c, nt)
        #phiCTCS = CTCS(phiOld.copy(), c, nt)
xmin = 0
xmax = 1
nx = 100
# Derived parameters
dx = (xmax - xmin)/nx
x = np.arange(xmin, xmax, dx)
phiOld = cosBell(x, 0, 0.75)

main(xmin, xmax, nx, dx, x, phiOld)
