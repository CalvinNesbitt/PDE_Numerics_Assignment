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
    c_values = np.arange(0.1, 2.2, 0.1) # Courant Values
    error_Matrix = np.zeros((len(c_values), 4)) # Matrix to store errors
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phiFTBS, phiAnalytic) # Store L2 errors
    np.savetxt('results/FTBS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1],'--g^', label = 'FTBS')

    # Error calculation for FTCS scheme
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phiFTCS = FTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phiFTCS, phiAnalytic) # Store L2 errors
    np.savetxt('results/FTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTCS scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1],'--r^', label = 'FTCS')

    # Error calculation for BTCS scheme
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phiBTCS = BTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phiBTCS, phiAnalytic) # Store L2 errors
    np.savetxt('results/BTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running BTCS scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1],'--c^', label = 'BTCS')

    # Error calculation for CTCS scheme
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phiCTCS = CTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phiCTCS, phiAnalytic) # Store L2 errors
    np.savetxt('results/CTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running CTCS scheme for different courant numbers')
    plt.plot(c_values, error_Matrix[:, 1],'--m^', label = 'CTCS')

    # Error calculation for Semi lagranian scheme
    c_values = np.arange(0.2, 2.4, 0.2)
    error_Matrix = np.zeros((len(c_values), 4)) # Matrix to store errors
    for j,c in enumerate(c_values):
        nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
        phi_sem_lag = sem_lag(phiOld.copy(), c, nt, x, dx)
        error_Matrix[j, 0] = c # Store Courant value
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
        error_Matrix[j, 1] = l2ErrorNorm(phi_sem_lag, phiAnalytic) # Store L2 errors
    np.savetxt('results/Sem_Lag_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running Semi Lagrangian scheme for different courant numbers')
    plt.semilogy(c_values, error_Matrix[:, 1],'--b^', label = 'Semi Lagrangian') # Logarithmic Axis

    #Plot Details
    plt.xlabel('Courant number c')
    plt.ylabel('L2 Error (Logarithmic)')
    plt.title('L2 Error as function of c')
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=5)
    plt.savefig('plots/L2Error_vs_Courant_Number', bbox_inches='tight')
    plt.show()

    return

xmin = 0
xmax = 1
nx = 100
# Derived parameters
dx = (xmax - xmin)/nx
x = np.arange(xmin, xmax, dx)
phiOld = cosBell(x, 0, 0.75)

main(xmin, xmax, nx, dx, x, phiOld)
