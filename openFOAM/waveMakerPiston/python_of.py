#velocity profile mean between given x values
from fluidfoam import readscalar, readmesh
from fluidfoam.readof import readvector
import matplotlib.pyplot as plt
from numpy import sqrt, linspace
import numpy as np
from stokes import *
sol = '../waveMakerPiston'
timename = '9.9'
alphawater = readscalar(sol, timename, 'alpha.water', structured=True)
X, Y, Z = readmesh(sol, structured=True)
U = readvector(sol, timename, 'U', structured=True)
nx, nz = X.shape[0], X.shape[2]
a = 0.0093
freq = 2.25
ak = 0.13
k = ak/a

omega = 2*pi*freq
g = 9.81

x_min = 2
x_max = 3 
idx_min = np.nanargmin(np.abs(X[:,0,0]-x_min))
idx_max = np.nanargmin(np.abs(X[:,0,0]-x_max))

print("shapes ", X.shape, U.shape) 
X = X[idx_min:idx_max,:,:]
Z = Z[idx_min:idx_max,:,:]
U = U[:,idx_min:idx_max,:,:]
alphawater = alphawater[idx_min:idx_max,:,:]
print("shapes ", X.shape, U.shape)

alpha = np.zeros((nz))
for i in range(nz):
    mask = (alphawater[:,:,i]>0.99).flatten()
    
    alpha[i] = np.mean(omega/(ak*g)*sqrt(U[0][mask,:,i]**2+U[2][mask,:,i]**2))

z = (Z[0,0,:]-0.335)/0.335
plt.figure()
plt.plot(alpha, z)


z0 = linspace(-1, a, 1000)
plt.plot(alpha_kz_stokes(k, z0), z0)
plt.show()
