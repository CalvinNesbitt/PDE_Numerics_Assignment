# Various function for plotting results and for calculating error measures

# Read in required packages and scripts
import matplotlib.pyplot as plt
import numpy as np
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def l2ErrorNorm(phi, phiExact):
    "Calculates the l2 error norm (RMS error) of phi in comparison to"
    "phiExact"

    # calculate the error and the RMS error norm
    phiError = phi - phiExact
    l2 = np.sqrt(sum(phiError**2)/sum(phiExact**2))

    return l2

def lInfErrorNorm(phi, phiExact):
    "Calculates the linf error norm (maximum normalised error) in comparison"
    "to phiExact"
    phiError = phi - phiExact
    return np.max(np.abs(phiError))/np.max(np.abs(phiExact))

def nmc_error(phi, phi_initial):
    "Calculates Normalised Mass Conservation error of phi given initial condition"
    "phi_initial"
    error = (phi.mean() - phi_initial.mean())/phi_initial.mean()
    return error
