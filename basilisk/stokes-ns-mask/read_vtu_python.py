from meshio import read
import matplotlib.pyplot as plt
from numpy import mean
import numpy as np
import scipy
   
filename = "vtu/ascii/TIME-000*_000017/TIME-000*_000017_0.vtu"

mesh = read(filename) 


cell_points = mesh.cells[0].data    
points = mesh.points

def get_cell_centers():
    return mean(points[cell_points],axis=1)
    

    

cellData = mesh.cell_data
mask = cellData['f'][0]>0.5

Ux = cellData['u.x']
U = Ux[0][:,0]
V = Ux[0][:,1]

X = get_cell_centers()[:,0]
Y = get_cell_centers()[:,1]
center_points = np.stack([X, Y], axis=1)
center_U = np.stack([U,V], axis=1)

plt.quiver(X[mask], Y[mask], U[mask], V[mask])
plt.show()

Y0 = np.linspace(-0.53, 0.0093, 1000)
X0 = -0.0755*np.ones_like(Y0)
interpolation_points = np.stack([X0, Y0], axis=1)

interpolated_U = scipy.interpolate.griddata(center_points, center_U, interpolation_points) 

plt.plot(interpolated_U[:,0], interpolation_points[:,1])
plt.xlabel("u[m/s]")
plt.ylabel("y")
plt.show()


matlab_dict  = {}
matlab_dict['X'] = interpolation_points
matlab_dict['U'] = interpolated_U
matlab_filename  = "basilisk_velocity_profile.mat"

scipy.io.savemat(matlab_filename, matlab_dict)