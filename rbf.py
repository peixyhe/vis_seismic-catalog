#################################################################################################
####  The program is used to visualize seismic catalog data based on RBF Kernel resampling.  ####
####  Author: He Pei; 2024.02.25                                                             ####
#################################################################################################

up = 33.4      # geographic region:
down = 21.8    #                           up(-90.0, 90.0)
left = 97.8    #     left(-180.0, 180.0)                     right(-180.0, 180.0)
right = 107    #                          down(-90.0, 90.0)

step = 0.05
sigma_square2 = 2.0 * 5    # square(sigma) * 2  is the parameter of RBF Kernel
a = 0.001                  # a is the parameter of RBF Kernel
scale = 10                 # Scaling factor



import vtk
import math
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

half_setp = round(  step * 0.5, 5  )
ANGLE2METERS = 2 * math.pi * 6371393.0 / 360.0    # EARTH_RADIUS = 6371393.0
z0 = scale * 6000 / ANGLE2METERS

def gaussian_kernel(  x0, y0, x1, y1  ):
    r = np.power(x1 - x0, 2) + np.power(y1 - y0, 2)
    return a * np.exp(  -1.0 * r / sigma_square2   )

colors = vtk.vtkNamedColors()
ugrid = vtk.vtkUnstructuredGrid()

points = vtk.vtkPoints()

freq = vtk.vtkFloatArray()
freq.SetNumberOfComponents(1)
freq.SetName('frequency')

max_mag = vtk.vtkFloatArray()
max_mag.SetNumberOfComponents(1)
max_mag.SetName('MAX magnitude')

rbf = vtk.vtkFloatArray()
rbf.SetNumberOfComponents(1)
rbf.SetName('RBF value')

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")

data = np.array(    df[  ['lon', 'lat', 'mag']  ].values.tolist()    )
data_kdTree = cKDTree(  data[:, :2]  )

x = [  round(x0, 5) for x0 in np.arange(left, right + half_setp, step)  ]
y = [  round(y0, 5) for y0 in np.arange(down, up + half_setp, step)  ]

# for yj in y:
#     for xi in x:
#         indices = data_kdTree.query_ball_point(  [xi, yj], r = step  )
#         if indices is not None:
#             hotmap_list = []
#             mag0 = 0.0
#             rbf0 = 0.0
            
#             for index in indices:
#                 p = data[index]
#                 rbf0 += gaussian_kernel(xi, yj, p[0], p[1])
#                 if p[2] > mag0:
#                     mag0 = p[2]

#             points.InsertNextPoint(xi, yj, z0 + rbf0)
#             freq.InsertNextValue(  len(indices)  )
#             max_mag.InsertNextValue(mag0)
#             rbf.InsertNextValue(rbf0)
#         else:
#             points.InsertNextPoint(xi, yj, z0)
#             freq.InsertNextValue(-1)
#             max_mag.InsertNextValue(-1)
#             rbf.InsertNextValue(-1)

for yj in y:
    for xi in x:
        hotmap_list = []
        mag0 = 0.0
        rbf0 = 0.0
        indices = data_kdTree.query_ball_point(  [xi, yj], r = step  )
        
        for index in indices:
            p = data[index]
            rbf0 += gaussian_kernel(xi, yj, p[0], p[1])
            if p[2] > mag0:
                mag0 = p[2]

        points.InsertNextPoint(xi, yj, z0 + rbf0)
        freq.InsertNextValue(  len(indices)  )
        max_mag.InsertNextValue(mag0)
        rbf.InsertNextValue(rbf0)

x_num = len(x)
y_num = len(y)

for j in range(y_num - 1):
    for i in range(x_num - 1):
        id0 = i + j * x_num
        id1 = id0 + 1
        id2 = id1 + x_num
        id3 = id0 + x_num
        """
           id3------id2
            |        |
            |        |
           id0------id1
        """
        ugrid.InsertNextCell(    vtk.VTK_TRIANGLE, 3, [id3, id0, id1]    )
        ugrid.InsertNextCell(    vtk.VTK_TRIANGLE, 3, [id3, id1, id2]    )    # ugrid.InsertNextCell(    vtk.VTK_QUAD, 4, [id0, id1, id2, id3]    )



ugrid.SetPoints(points)
ugrid.GetPointData().AddArray(freq)
ugrid.GetPointData().AddArray(max_mag)
ugrid.GetPointData().AddArray(rbf)

writer = vtk.vtkUnstructuredGridWriter()
writer.SetInputData(ugrid)
writer.SetFileName(  "rbf_step-" + str(step) + ".vtk"  )
# writer.SetDataModeToAscii()
writer.Update()

# Create a mapper and actor
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(ugrid)

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
actor.GetProperty().SetPointSize(2)

# Visualize
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.SetWindowName('Polyhedron')
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(actor)
renderer.SetBackground(colors.GetColor3d('Salmon'))
renderer.ResetCamera()
renderer.GetActiveCamera().Azimuth(30)
renderer.GetActiveCamera().Elevation(30)
renderWindow.Render()
renderWindowInteractor.Start()