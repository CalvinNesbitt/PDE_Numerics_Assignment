# Mass Conservation Analysis

def nmc_error(phi, phi_initial):
    "Calculates Normalise Mass Conservation error of phi given initial condition"
    "phi_initial"
    error = (phi.mean - phi_initial.mean)/phi_initial.mean
    return error

# Parameters
xmin = 0
xmax = 1
nx = 100
nt = 100
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

# Calculate NMC error for different nt

# FTBS scheme
nt_values = np.arange(0.1, 2.2, 0.1) # Courant Values
error_Matrix = np.zeros((len(c_values), 2)) # Matrix to store errors
for j,c in enumerate(c_values):
    nt = int(1 / (c * dx)) # nt must be changed that we run the scheme for same amount of time t
    phiFTBS = FTBS(phiOld.copy(), c, nt) # Run the scheme for each (c, nt) pair
    error_Matrix[j, 0] = c # Store Courant value
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
    error_Matrix[j, 1] = nmc_error(phiFTBS, phiOld) # Store NMC errors
np.savetxt('results/FTBS_Errors.csv',error_Matrix, delimiter=',', newline='\n', header='Results of running FTBS scheme for different courant numbers')
plt.plot(c_values, error_Matrix[:, 1],'--g^', label = 'FTBS')
