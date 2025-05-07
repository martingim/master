from meshio import read
import matplotlib.pyplot as plt
import pyvista as pv
from numpy import mean
import numpy as np
import scipy
import glob
import os
import sys

basilisk_dir = sys.argv[1]
timestep = int(sys.argv[2])
save_dir  = basilisk_dir+ "/vts/matlab/"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

nx = 10000 # number of points to save in x direction in the matlab matrix
ny = 1000 # y dir
filenames = glob.glob(basilisk_dir + '/vts/*.vts')

if(type(filenames)==str):
    filenames = [filenames]
assert(len(filenames)>0)
filenames.sort()

i = timestep
for filename in filenames:
    reader = pv.get_reader(filenames[i])
    mesh = reader.read()
    centers = mesh.cell_centers().points
    centers = centers.reshape((mesh.dimensions[1]-1,mesh.dimensions[0]-1, 3))
    print(centers.shape)
    print(mesh['eta'].shape)

    V = mesh['Velocity'].reshape((mesh.dimensions[1]-1,mesh.dimensions[0]-1, 3))    
    eta = mesh['eta'].reshape((mesh.dimensions[1]-1,mesh.dimensions[0]-1))

    matlab_dict = {}
    matlab_dict['U'] = V
    matlab_dict['X'] = centers
    
    
    try:
        os.remove(save_dir + f"timestep_{i}.mat")
    except:
        pass
    scipy.io.savemat(save_dir + f"timestep_{i}.mat", matlab_dict)
    
    i += 1
