# This script plots the l2 and LInf Errors of the FTCS and BTCS schems as a function of time steps
import matplotlib.pyplot as plt

import numpy as np
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def TimeStepErrors(xmin, xmax, nx, nt, c):
    plt.clf() # Clear any current figures
    plt.close()

    # Derived parameters
    dx = (xmax - xmin)/nx
    x = np.arange(xmin, xmax, dx)

    # Vectors storing data to be plotted

    nt_values = np.arange(1, 201, 1) # Number of time steps
    ftcs_l2_errors = np.zeros(nt_values.size) # L2 Errors
    btcs_l2_errors = np.zeros(nt_values.size)
    ftcs_lInf_errors = np.zeros(nt_values.size) # LInf Errorrs
    btcs_lInf_errors = np.zeros(nt_values.size)

    phiOld = cosBell(x, 0, 0.75) # Initial condition

    # Calculating Errors for differnt time steps
    for it in range(nt_values.size):
        phiAnalytic = cosBell((x - c* it *dx)%(xmax - xmin), 0, 0.75)

        # Running schemes for different number of time steps
        phiFTCS = FTCS(phiOld.copy(), c, it)
        phiBTCS = BTCS(phiOld.copy(), c, it)

        # Store Errors
        ftcs_l2_errors[it] = l2ErrorNorm(phiFTCS, phiAnalytic)
        ftcs_lInf_errors[it] = lInfErrorNorm(phiFTCS, phiAnalytic)
        btcs_l2_errors[it] = l2ErrorNorm(phiBTCS, phiAnalytic)
        btcs_lInf_errors[it] = lInfErrorNorm(phiBTCS, phiAnalytic)

    # Plotting Errors

    fig, ax = plt.subplots(2, 1) # Create two subplots
    ax[0].set(ylabel='L2 Error',
    title='L2 Error Comparison')
    ax[0].plot(nt_values, ftcs_l2_errors, 'r', label = 'FTCS')
    ax[0].plot(nt_values, btcs_l2_errors, 'b', label = 'BTCS')
    fig.legend()

    ax[1].set(xlabel='Number of timesteps (nt)', ylabel='LInf Error',
    title='LInf error comparison')
    ax[1].plot(nt_values, ftcs_lInf_errors, 'r', nt_values, btcs_lInf_errors, 'b')

    plt.show()
    input('press return to save file and continue')
    fig.savefig('plots/TimeStepErrors.png')

    return
