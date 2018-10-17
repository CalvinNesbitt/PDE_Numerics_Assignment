# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py

# If you are using Python 2.7 rather than Python 3, import various
# functions from Python 3 such as to use real number division
# rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1
#from __future__ import absolute_import, division, print_function

# The numpy package for numerical functions and pi
import numpy as np
from Matrix_test import *

def FTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # FTCS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - 0.5*c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])

        # update arrays for next time-step
        phiOld = phi.copy()

    return phi

def BTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using BTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # Matrix for calculating phi at the next time step
    timeStepMatrix = btcsTimeStep(nx, c)

    # BTCS for each time-step
    for it in range(nt):
        # Solve Matrix equation at each time step to calculate new phi
        phi =  np.linalg.solve(timeStepMatrix, phiOld)

        # update arrays for next time-step
        phiOld = phi.copy()

    return phi
