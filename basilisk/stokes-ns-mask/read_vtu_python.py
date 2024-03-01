from meshio import read
import matplotlib.pyplot as plt
from numpy import mean
import numpy as np
import scipy
import glob

matlab_filename  = "basilisk_velocity_profile.mat"
matlab_dict  = {}


nx = 100 #number of points to save in x direction in the matlab matrix
ny = 100 # y dir
filenames = glob.glob('vtu/ascii/*.vtu')
filenames.sort()
 
i = 0
for filename in filenames:
    print('reading file :', filename)
    mesh = read(filename)
    cell_points = mesh.cells[0].data
    points = mesh.points
    #calculate the cell centers
    X = mean(points[cell_points],axis=1)[:,0]
    Y = mean(points[cell_points],axis=1)[:,1]
    center_points = np.stack([X, Y], axis=1)

    
    # get the data for the center of the cells
    cellData = mesh.cell_data
    f = cellData['f'][0]
    mask = f>0.5

    Ux = cellData['u.x']
    U = Ux[0][:,0]
    V = Ux[0][:,1]
    center_U = np.stack([U,V], axis=1)
    x_interpolation_points = np.linspace(min(X), max(X), nx)
    y_interpolation_points = np.linspace(min(Y), max(Y), ny)
    x,y= np.meshgrid(x_interpolation_points, y_interpolation_points)
    interpolation_points = np.stack([x,y], axis=2)
    interpolated_U = scipy.interpolate.griddata(center_points, center_U, interpolation_points)
    interpolated_f = scipy.interpolate.griddata(center_points, f, interpolation_points)
    res = {}
    res['U'] = interpolated_U
    res['X'] = interpolation_points
    res['mask'] = interpolated_f>0.5 
    matlab_dict[f"timestep_{i}"] = res
    i += 1
scipy.io.savemat(matlab_filename, matlab_dict)

"""
### FIND THE WAVE CRESTS ###
#multiply the y coordinates with the water fraction to find the highest point of the wave
Y_F = Y*f
ind = Y_F.argmax()

Y0 = np.linspace(min(Y), max(Y[mask]), 1000)
X0 = X[ind]*np.ones_like(Y0)
interpolation_points = np.stack([X0, Y0], axis=1)

interpolated_U = scipy.interpolate.griddata(center_points, center_U, interpolation_points) 

plt.plot(interpolated_U[:,0], interpolation_points[:,1])
plt.xlabel("u[m/s]")
plt.ylabel("y")
plt.show()
i += 1

matlab_dict['X'] = interpolation_points
matlab_dict['U'] = interpolated_U


plt.quiver(X[mask], Y[mask], U[mask], V[mask])
plt.show()
"""