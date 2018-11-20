# Matrix for calculating the next time step in BTCS scheme

import numpy as np

def btcsTimeStep(nx, c):

    M = np.zeros((nx, nx))

    #Bounday conditions
    M[0,0], M[0, 1], M[0, nx-1] = 1, 0.5*c, -0.5*c
    M[nx -1,0], M[nx-1, nx-2], M[nx-1, nx-1] = 0.5*c, -0.5*c, 1

    for j in range(1, nx-1):
        M[j, j] = 1
        M[j, j-1] = -0.5*c
        M[j, j+1] = 0.5*c

    return M
