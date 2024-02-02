from meshio import read
import matplotlib.pyplot as plt
from numpy import mean
import numpy as np
   
filename = "vtu/ascii/ascii_0.vtu"

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


plt.quiver(X[mask], Y[mask], U[mask], V[mask])
plt.show()