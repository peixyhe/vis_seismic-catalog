###############################################################################
####  The program is used to polt M-T-F diagrams of seismic catalog data.  ####
####  Author: He Pei; 2024.02.25                                           ####
###############################################################################



num_mouths = 1

year0 = 2009
year1 = 2022
mag_step = 0.05

if num_mouths == 1:
    minification = 20.0
    t_step = 0.05 * 6.0
elif num_mouths == 3:
    minification = 30.0
    t_step = 0.05 * 2.0
elif num_mouths == 6:
    minification = 50.0
    t_step = 0.05
else:
    minification = 10.0



import vtk
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from datetime import datetime

def seconds_between_times(time1, time2):
    time1_dt = datetime.strptime(time1, '%Y.%m.%d %H:%M')    # time1_dt = datetime.strptime(time1, '%Y-%m-%d %H:%M:%S')
    time2_dt = datetime.strptime(time2, '%Y.%m.%d %H:%M')
    
    seconds_diff = (time2_dt - time1_dt).total_seconds()
    
    return abs(seconds_diff) / (3600.0 * 24 * 30 * num_mouths)    # 3600.0 * 24 * 30 * 6 = seconds for 6 mouths

colors = vtk.vtkNamedColors()
ugrid = vtk.vtkUnstructuredGrid()

points = vtk.vtkPoints()

freq = vtk.vtkFloatArray()
freq.SetNumberOfComponents(1)
freq.SetName('frequency')

max_mag = vtk.vtkFloatArray()
max_mag.SetNumberOfComponents(1)
max_mag.SetName('MAX magnitude')

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")
time = df['time'].to_list()
mag = df['mag']

t_min = time[-1]
t_max = time[0]

x0 = 0.0
x1 = seconds_between_times(t_min, t_max)
y0 = min(mag)
y1 = max(mag)

k = (year1 - year0) / (x1 - x0)
b = year0 - k * x0

half_mag_setp = round(  mag_step * 0.5, 5  )
half_t_step = round(    t_step * 0.5, 5    )
r_step = half_mag_setp + half_t_step

xy = []
for i in range(len(mag)):
    xi = seconds_between_times(  t_min, time[i]  )
    xy.append(  [xi, mag[i] ]  )

xy_kdTree = cKDTree(xy)

x = [    round(xx, 5) for xx in np.arange(x0, x1 + half_t_step,   t_step)    ]
y = [    round(yy, 5) for yy in np.arange(y0, y1 + half_mag_setp, mag_step)    ]

for yj in y:
    for xi in x:
        p_list = []
        mag0 = 0.0
        indices = xy_kdTree.query_ball_point(  [xi, yj], r = r_step  )
        
        for index in indices:
            p = xy[index]
            if (p[0] >= xi - half_t_step) and (p[0] < xi + half_t_step) and (p[1] >= yj - half_mag_setp) and (p[1] < yj + half_mag_setp):
                p_list.append(index)
                if p[1] > mag0:
                    mag0 = p[1]
                    
        points.InsertNextPoint(    k * xi + b, yj, len(p_list) / minification    )
        freq.InsertNextValue(  len(p_list)  )
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
writer.SetFileName(  "m_t_f-" + str(num_mouths) + "_mouths.vtu"  )
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
