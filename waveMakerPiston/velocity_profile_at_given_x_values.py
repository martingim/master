from fluidfoam import readscalar, readmesh
from fluidfoam.readof import readvector
import matplotlib.pyplot as plt
from numpy import sqrt, linspace
import numpy as np
from stokes import *
sol = '../waveMakerPiston'
timename = '9.9'
alphawater = readscalar(sol, timename, 'alpha.water')
X, Y, Z = readmesh(sol)
U = readvector(sol, timename, 'U')

a = 0.0093
freq = 2.25
ak = 0.13
k = ak/a

omega = 2*pi*freq
g = 9.81

plt.figure()
###plots at certain x values
X0 = [1, 1.1, 1.2, 1.3, 1.4]


for x0 in X0:
    idx = np.nanargmin(np.abs(X-x0))
    x0 = X[idx]
    indexes = (abs(X-x0))<0.000001
    x = X[indexes]
    z = Z[indexes]
    u = U[:,indexes]
    print(u.shape)
    n_water = alphawater[indexes]
    n = sum(n_water>0.5)
    z = (z[0:n]-0.335)/0.335
    alpha = omega/(ak*g)*sqrt(u[0][0:n]**2+u[2][0:n]**2)
    print('alphashpae ', alpha.shape)
    plt.plot(alpha,z)
    print(z[-1])

z0 = linspace(-1, a, 1000)
plt.plot(alpha_kz_stokes(k, z0), z0)
plt.show()
