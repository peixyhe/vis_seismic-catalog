#############################################################################
####  The program is used to polt M-T-F diagrams of seismic catalog data.  ####
####  Author: He Pei; 2024.02.25                                         ####
#############################################################################

mag_step = 0.05
t_step = 0.05
r_step = (mag_step**2.0 + t_step**2.0) ** 0.5

sigma_square2 = 2.0 * 5
a = 0.005
year0 = 2009
year1 = 2022



import vtk
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from datetime import datetime

half_mag_setp = round(  mag_step * 0.5, 5  )
half_t_step = round(  t_step * 0.5, 5  )

def gaussian_kernel(  x0, y0, x1, y1  ):
    r = np.power(x1 - x0, 2) + np.power(y1 - y0, 2)
    return a * np.exp(  -1.0 * r / sigma_square2   )

def seconds_between_times(time1, time2):
    time1_dt = datetime.strptime(time1, '%Y.%m.%d %H:%M')    # time1_dt = datetime.strptime(time1, '%Y-%m-%d %H:%M:%S')
    time2_dt = datetime.strptime(time2, '%Y.%m.%d %H:%M')
    
    seconds_diff = (time2_dt - time1_dt).total_seconds()
    
    return abs(seconds_diff) / (3600.0 * 24 * 30 * 6)    # 3600.0 * 24 * 30 * 6 = seconds for 6 mouths

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

xy = []
for i in range(len(mag)):
    xi = seconds_between_times(  t_min, time[i]  )
    xy.append(  [ k * xi + b, mag[i] ]  )
    
xy_kdTree = cKDTree(xy)

x = [round(xx, 5) for xx in np.arange(x0 * k + b, x1 * k + b + half_t_step, t_step)]
y = [round(yy, 5) for yy in np.arange(y0, y1 + half_mag_setp, mag_step)]

for yj in y:
    for xi in x:
        p_list = []
        mag0 = 0.0
        rbf0 = 0.0
        indices = xy_kdTree.query_ball_point(  [xi, yj], r = r_step  )
        
        for index in indices:
            p = xy[index]
            rbf0 += gaussian_kernel(xi, yj, p[0], p[1])
            if p[1] > mag0:
                mag0 = p[1]

        # points.InsertNextPoint(xi, yj, rbf0)
        points.InsertNextPoint(xi, yj, len(indices) / 200.0)
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
writer.SetFileName(  "m_t_f.vtu"  )
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