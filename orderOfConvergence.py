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

nx_values = np.arange(2, 100, 1)
dx_values = (xmax - xmin)/nx_values # Calculate spatial step
c = 0.5

# Error calculation for FTBS scheme
error_Matrix = np.zeros((len(dx_values), 4)) # Matrix to store errors
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = int(nx_values[i] / c ) # nt must be changed that we run the scheme for same amount of time t
    phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 0] = dx_values[i] # Store dx value
    phiAnalytic = cosBell((x - c*nt*dx_values[i])%(xmax - xmin), 0, 0.75)
    error_Matrix[i, 1] = l2ErrorNorm(phiFTBS, phiAnalytic) # Store L2 errors
np.savetxt('results/convergence_analysis/FTBS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different spatial discretisation')
plt.plot(dx_values, error_Matrix[:, 1],'--g^', label = 'FTBS')

# Error calculation for FTCS scheme
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = int(nx_values[i] / c ) # nt must be changed that we run the scheme for same amount of time t
    phiFTCS = FTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 0] = dx_values[i] # Store dx value
    phiAnalytic = cosBell((x - c*nt*dx_values[i])%(xmax - xmin), 0, 0.75)
    error_Matrix[i, 1] = l2ErrorNorm(phiFTCS, phiAnalytic) # Store L2 errors
np.savetxt('results/convergence_analysis/FTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different spatial discretisation')
plt.plot(dx_values, error_Matrix[:, 1],'--r^', label = 'FTCS')

# Error calculation for BTCS scheme
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = int(nx_values[i] / c ) # nt must be changed that we run the scheme for same amount of time t
    phiBTCS = BTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 0] = dx_values[i] # Store dx value
    phiAnalytic = cosBell((x - c*nt*dx_values[i])%(xmax - xmin), 0, 0.75)
    error_Matrix[i, 1] = l2ErrorNorm(phiBTCS, phiAnalytic) # Store L2 errors
np.savetxt('results/convergence_analysis/BTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different spatial discretisation')
plt.plot(dx_values, error_Matrix[:, 1],'--c^', label = 'BTCS')

# Error calculation for CTCS scheme
for i,dx in enumerate(dx_values):
    x = np.arange(xmin, xmax, dx)
    phiOld = cosBell(x, 0, 0.75) # Initial condition for given discretisation
    nt = int(nx_values[i] / c ) # nt must be changed that we run the scheme for same amount of time t
    phiCTCS = CTCS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[i, 0] = dx_values[i] # Store dx value
    phiAnalytic = cosBell((x - c*nt*dx_values[i])%(xmax - xmin), 0, 0.75)
    error_Matrix[i, 1] = l2ErrorNorm(phiCTCS, phiAnalytic) # Store L2 errors
np.savetxt('results/convergence_analysis/CTCS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different spatial discretisation')
plt.plot(dx_values, error_Matrix[:, 1],'--m^', label = 'CTCS')

#Plot Details
plt.plot(dx_values, dx_values, label = 'x')
plt.loglog(dx_values, dx_values**2, label = 'x^2')
plt.xlabel('Spatial step (dx)')
plt.ylabel('L2 Error (Logarithmic)')
plt.title('L2 Error as function of dx')
plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=5)
plt.savefig('plots/L2Error_vs_dx', bbox_inches='tight')
plt.show()
