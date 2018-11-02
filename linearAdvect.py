
# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.

### Copy out most of this code. Code commented with 3#s (like this) ###
### is here to help you to learn python and need not be copied      ###


### ./linearAdvect.py                                              ###

### Note that blocks are defined by indentation in Python. You     ###
### should never mix tabs and spaces for indentation - use 4 spaces.###
### Setup your text editor to insert 4 spaces when you press tab    ###

### If you are using Python 2.7 rather than Python 3, import various###
### functions from Python 3 such as to use real number division     ###
### rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1###
#from __future__ import absolute_import, division, print_function

### The matplotlib package contains plotting functions              ###
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *
from TimeStep_Error_Plotting import *

### The main code is inside a function to avoid global variables    ###
def main():
    "Advect the initial conditions using various advection schemes and"
    "compare results"

    # Parameters
    xmin = 0
    xmax = 1
    nx = 40
    nt = 30
    c = 0.2

    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial conditions
    phiOld = cosBell(x, 0, 0.75)
    #phiOld = squareWave(x, 0, 0.75)
    #phiOld = mixed(x, 0, 0.3,0.5,0.75)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)

    # Advect the profile using finite difference for all the time steps
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiBTCS = BTCS(phiOld.copy(), c, nt)
    phi_sem_lag = sem_lag(phiOld.copy(), c, nt, x, dx)

    # Calculate and print out error norms
    print("FTCS l2 error norm = ", l2ErrorNorm(phiFTCS, phiAnalytic))
    print("FTCS linf error norm = ", lInfErrorNorm(phiFTCS, phiAnalytic))
    print("BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic))
    print("BTCS linf error norm = ", lInfErrorNorm(phiBTCS, phiAnalytic))
    print("Semi Lagrangian l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic))
    print("Semi Lagrangian linferror norm = ", lInfErrorNorm(phiBTCS, phiAnalytic))

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS, label='FTCS', color='red')
    plt.plot(x, phiBTCS, label='BTCS', color='blue')
    plt.plot(x, phi_sem_lag, label='SL', color='green')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend(bbox_to_anchor=(0.5, 0.5))
    plt.xlabel('$x$')
    plt.show()
    input('press return to save file and see timestep error comparison')
    plt.savefig('plots/SchemeComparisons.png')

    # Plotting errors as a function of number of time steps
    TimeStepErrors(xmin, xmax, nx, nt, c)

### Run the function main defined in this file                      ###
main()
