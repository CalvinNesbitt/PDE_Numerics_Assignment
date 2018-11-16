# Code for advecting a cosine wave by each scheme. Plots of results are saved in
# directory plots/illustrations.

import matplotlib.pyplot as plt

# Read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def main(xmin, xmax, nx, nt, c):
    "Advect the initial conditions using various advection schemes and"
    "illustrate results in a plot."

    # Derived parameters
    dx = (xmax - xmin)/nx

    # Spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial condition
    phiOld = cosBell(x, 0, 0.75)

    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)

    # Advect the profile using finite difference for all the time steps
    phiFTBS = FTBS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiBTCS = BTCS(phiOld.copy(), c, nt)
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phi_sem_lag = sem_lag(phiOld.copy(), c, nt, x, dx)

    # Plotting the solutions

    # FTCS plotted seperately as it is unstable
    plt.figure(1)
    plt.plot(x, phiOld,'k', label='Initial')
    plt.plot(x, phiAnalytic,'--k', label='Analytic')
    plt.plot(x, phiFTCS, label='FTCS', color='red')
    title = 'nt = %s, nx = %s, c = %s'%(nt, nx, c)
    plt.title(title)
    plt.legend()
    axes = plt.gca()
    axes.set_ylim([-0.1,1.1])
    plt.xlabel('x')
    save_Location = 'plots/illustrations/FTCS_%s.png' % c
    plt.savefig(save_Location, dpi=1000)
    plt.clf()

    plt.figure(2)
    plt.plot(x, phiOld,'k', label='Initial')
    plt.plot(x, phiAnalytic,'--k', label='Analytic')
    plt.plot(x, phiFTBS,'g', label='FTBS', alpha = 0.6)
    plt.plot(x, phiBTCS,'c', label='BTCS', alpha = 0.6)
    plt.plot(x, phiCTCS,'m', label='CTCS', alpha = 0.6)
    plt.plot(x, phi_sem_lag,'y', label='SL', alpha = 0.6)
    plt.title(title)
    plt.legend()
    axes = plt.gca()
    axes.set_ylim([-0.1,1.1])
    plt.xlabel('x')
    save_Location = 'plots/illustrations/Stable_schemes_%s.png' % c
    plt.savefig(save_Location, dpi=1000)
    plt.close('all')

    return
