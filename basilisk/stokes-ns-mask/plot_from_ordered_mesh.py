from meshio import read
import matplotlib.pyplot as plt


filename = "vtk/TIME-000000.vtk"

mesh = read(filename)


cell_points = mesh.cells[0]


points = mesh.points
pointData = mesh.point_data
mask = pointData['f']>0.5
X = points[:,0]
Y = points[:,1]
U = pointData['u.x']
V = pointData['u.y']
plt.quiver(X[mask[:,0]], Y[mask[:,0]], U[mask[:,0]], V[mask[:,0]])
plt.show()