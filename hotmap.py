################################################################################
####  The program is used to polt hotmap diagrams of seismic catalog data.  ####
####  Author: He Pei; 2024.02.25                                            ####
################################################################################

up = 33.4      # geographic region:
down = 21.8    #                           up(-90.0, 90.0)
left = 97.8    #     left(-180.0, 180.0)                     right(-180.0, 180.0)
right = 107    #                          down(-90.0, 90.0)

step = 0.01
scale = 10
window_Width = 5



import vtk
import math
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

half_setp = round(  step * 0.5, 5  )
ANGLE2METERS = 2 * math.pi * 6371393.0 / 360.0    # EARTH_RADIUS = 6371393.0
z0 = scale * 6000 / ANGLE2METERS
search_radius = window_Width * step / 2.0
core_radius = search_radius * 1.5

colors = vtk.vtkNamedColors()
points = vtk.vtkPoints()
ugrid = vtk.vtkUnstructuredGrid()

freq = vtk.vtkFloatArray()
freq.SetNumberOfComponents(1)
freq.SetName('frequency')

max_mag = vtk.vtkFloatArray()
max_mag.SetNumberOfComponents(1)
max_mag.SetName('MAX magnitude')

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")

data = np.array(    df[  ['lon', 'lat', 'mag']  ].values.tolist()    )
data_kdTree = cKDTree(  data[:, :2]  )

x = [  round(x0, 5) for x0 in np.arange(left, right + half_setp, step)  ]
y = [  round(y0, 5) for y0 in np.arange(down, up + half_setp, step)  ]

for yj in y:
    print(yj)
    for xi in x:
        points.InsertNextPoint(xi, yj, z0)
        freq0 = 0
        mag0 = 0.0
        
        indices = data_kdTree.query_ball_point(  [xi, yj], r = core_radius  )
        if indices is not None:
            hotmap_list = []
            
            for index in indices:
                p = data[index]
                if (p[0] >= xi - search_radius) and (p[0] < xi + search_radius) and (p[1] >= yj - search_radius) and (p[1] < yj + search_radius):
                    hotmap_list.append(index)
                    if p[2] > mag0:
                        mag0 = p[2]
            freq0 = len(hotmap_list)
            
        freq.InsertNextValue(freq0)
        max_mag.InsertNextValue(mag0)

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

writer = vtk.vtkUnstructuredGridWriter()
writer.SetInputData(ugrid)
writer.SetFileName(  "hotmap_step-" + str(step) + ".vtk"  )
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
