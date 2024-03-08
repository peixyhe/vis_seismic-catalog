######################################################################################
####  The program is used to polt M-T-F_hotmap diagrams of seismic catalog data.  ####
####  Author: He Pei; 2024.02.25                                                  ####
######################################################################################



step_mouths = 5
min_mag = 1.0



step_mag = 0.1

import vtk
# import math
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from datetime import datetime

def seconds_between_times(time1, time2):
    time1_dt = datetime.strptime(time1, '%Y.%m.%d %H:%M')    # time1_dt = datetime.strptime(time1, '%Y-%m-%d %H:%M:%S')
    time2_dt = datetime.strptime(time2, '%Y.%m.%d %H:%M')
    
    return abs(  (time2_dt - time1_dt).total_seconds()  ) / (3600.0 * 24 * 30 * step_mouths)

colors = vtk.vtkNamedColors()
ugrid = vtk.vtkUnstructuredGrid()
points = vtk.vtkPoints()

max_mag = vtk.vtkFloatArray()
max_mag.SetNumberOfComponents(1)
max_mag.SetName('MAX magnitude')

freq = vtk.vtkFloatArray()
freq.SetNumberOfComponents(1)
freq.SetName('frequency')

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")
time = df['time'].to_list()
mag = df['mag']

x0 = 0.0
x1 = seconds_between_times("2009.1.1 0:0", "2021.12.31 23:59")
# y0 = min(mag)
y0 = min_mag
y1 = max(mag)

year0 = 2009
year1 = 2022
k = (year1 - year0) / (x1 - x0)
print(1 / k, -1.0 / k * 2009)

xy = []
for i in range(  len(mag)  ):
    if mag[i] >= min_mag:
        xi = seconds_between_times(  "2009.1.1 0:0", time[i]  )
        xy.append(    [  xi, mag[i] ]    )
xy_kdTree = cKDTree(xy)

half_mag_step = step_mag * 0.5
x_step = step_mag
half_x_step = x_step * 0.5
kdTree_radius = half_x_step + half_mag_step

x = [    round(xx, 5) for xx in np.arange(x0, x1 + half_x_step, x_step)    ]
y = [    round(yy, 5) for yy in np.arange(y0, y1 + half_mag_step, step_mag)    ]

for yj in y:
    for xi in x:
        p_list = []
        mag0 = 0.0
        indices = xy_kdTree.query_ball_point(  [xi, yj], r = kdTree_radius  )
        
        for index in indices:
            p = xy[index]
            if (p[0] >= xi - half_x_step) and (p[0] < xi + half_x_step) and (p[1] >= yj - half_mag_step) and (p[1] < yj + half_mag_step):
                p_list.append(index)
                if p[1] > mag0:
                    mag0 = p[1]

        points.InsertNextPoint(    xi, yj, 0.0    )
        freq.InsertNextValue(    len(p_list) + 1    )
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
writer.SetFileName(  "hotmap_M-T_minMAG" + str(min_mag) + "_" + str(step_mouths) + "mouths.vtk"  )
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
