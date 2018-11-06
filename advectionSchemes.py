# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py

# If you are using Python 2.7 rather than Python 3, import various
# functions from Python 3 such as to use real number division
# rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1
#from __future__ import absolute_import, division, print_function

# The numpy package for numerical functions and pi
import numpy as np
from Matrix_test import *
from scipy.interpolate import lagrange

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
        phi =  np.linalg.solve(timeStepMatrix, phi)

    return phi

def CTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phiOlder = phiOld.copy()
    phiOld = FTCS(phiOlder, c, 1) #Use FTCS for first time step
    phi = phiOld.copy()
    # CTCS for each time-step after the first
    for it in range(1, nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOlder[j] - c*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])

        # update arrays for next time-step
        phiOlder = phiOld.copy()
        phiOld = phi.copy()

    return phi

def sem_lag(phiOld, c, nt, x, dx):
    "Linear advection of profile in phiOld using Semi Lagranian scheme, Courant number c"
    "for nt time-steps, spatial points x, spatial grid steps of length dx"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    for it in range(nt):

        for j in range(nx):

            # Building down wind interpolating polynomial
                # Find base_points of interpolating polynomial
            k = int(np.floor(j - c)) # Index below advected point
            x_base_points = np.array([x[k-1],x[k],x[(k+1)%nx],x[(k+2)%nx]])
            y_base_points = np.array([phiOld[k-1],phiOld[k],phiOld[(k+1)%nx],phiOld[(k+2)%nx]])
            poly = lagrange(x_base_points, y_base_points)

            # Finding down-wind point
            beta = j - c - k
            x_jd = beta * dx + x[k]

            # Calculate phi upwind
            phi[j] = poly(x_jd)
        phiOld = phi.copy() # Update for next time step

    return phi
