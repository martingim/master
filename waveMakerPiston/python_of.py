
from fluidfoam import readscalar, readmesh
from fluidfoam.readof import readvector
import matplotlib.pyplot as plt
from numpy import sqrt, linspace
import numpy as np
from stokes import *
sol = '../waveMakerPiston'
timename = '9.1'
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

x_Start = 1.
x_end = 4.

start_idx = np.nanargmin(abs(X[:,0,0]-x_Start))
end_idx = np.nanargmin(abs(X[:,0,0]-x_end))

print("start end x ", start_idx, end_idx)
print('shapes ', X.shape, U.shape, alphawater.shape)


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
