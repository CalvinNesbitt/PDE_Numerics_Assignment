# Outer script for running advection and scheme anaylsis for Numerics assginment.
import linearAdvect
import courantStabilityAnalysis
import orderOfAccuracy
import massConservationAnalysis
import time

# Parameters
xmin = 0 # Minimum spatial grid value
xmax = 1 # Maxmimum spatial grid value
nx = 80 # Number of spatial grid points
nt = 100 # Number of time steps

# Illustration of scheme results for different courant numbers
for c in [0.5, 1.2]:
    linearAdvect.main(xmin, xmax, nx, nt, c)

print('Illustrations of advected cosine wave saved.')

# Stability Analysis plots
courantStabilityAnalysis.main(0, 1, 50, 0.1, 1.8, 100)
print('Stability analysis plots saved.')

# Order of Accuracy plots
orderOfAccuracy.main(0, 1, 0.5)
print('Order of accuracy analysis plots saved.')

# Mass Conservation plots
massConservationAnalysis.main(0, 1, 0.5)

print('Mass conservation plots saved.')
print('Done.')
