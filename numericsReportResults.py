# Outer script for running advection and scheme anaylsis for Numerics assginment.
import linearAdvect
import courantStabilityAnalysis
import orderOfAccuracy
import massConservationAnalysis
import time

start_time = time.time()


# Parameters
xmin = 0 # Minimum spatial grid value
xmax = 1 # Maxmimum spatial grid value
nx = 80 # Number of spatial grid points
nt = 100

# Illustration of scheme results for different courant numbers
for c in [0.9]:
    linearAdvect.main(xmin, xmax, nx, nt, c)

# Stability Analysis plots
courantStabilityAnalysis.main(0, 1, 50, 0.1, 1.8, 100)

# Order of Accuracy plots

orderOfAccuracy.main(0, 1, 0.5)

# Mass Conservation plots

massConservationAnalysis.main(0, 1, 0.5)

print("--- %s seconds ---" % (time.time() - start_time))
