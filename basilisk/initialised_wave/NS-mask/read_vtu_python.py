from meshio import read
import matplotlib.pyplot as plt
from numpy import mean
import numpy as np
import scipy
import glob

save_dir  = "/home/martin/Documents/master/matlab/PIV_basilisk/basilisk_results/"


nx = 100 #number of points to save in x direction in the matlab matrix
ny = 1000 # y dir
filenames = glob.glob('vtu/ascii/*.vtu')
assert(len(filenames)>0)
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
    mask = f>0.9

    Ux = cellData['u.x']
    U = Ux[0][:,0]
    V = Ux[0][:,1]
    center_U = np.stack([U,V], axis=1)
    x_interpolation_points = np.linspace(min(X), max(X), nx)
    y_interpolation_points = np.linspace(max(Y), min(Y), ny)
    x,y= np.meshgrid(x_interpolation_points, y_interpolation_points)
    interpolation_points = np.stack([x,y], axis=2)
    interpolated_U = scipy.interpolate.griddata(center_points, center_U, interpolation_points)
    interpolated_f = scipy.interpolate.griddata(center_points, f, interpolation_points)
    matlab_dict = {}
    matlab_dict['U'] = interpolated_U
    matlab_dict['X'] = interpolation_points
    matlab_dict['mask'] = interpolated_f>0.5
    
    scipy.io.savemat(save_dir + f"timestep_{i}.mat", matlab_dict)
    i += 1
