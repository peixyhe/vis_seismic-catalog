#############################################################################
####  The program is used to polt M-T diagrams of seismic catalog data.  ####
####  Author: He Pei; 2024.02.25                                         ####
#############################################################################

year0 = 2009
year1 = 2022
num_mouths = 1


import vtk
import pandas as pd
from datetime import datetime

# Create the polydata where we will store all the geometric data
linesPolyData = vtk.vtkPolyData()
namedColors = vtk.vtkNamedColors()

points = vtk.vtkPoints()
cell_lines = vtk.vtkCellArray()

magnitude = vtk.vtkFloatArray()
magnitude.SetNumberOfComponents(1)
magnitude.SetName('magnitude')

def seconds_between_times(time1, time2):
    time1_dt = datetime.strptime(time1, '%Y.%m.%d %H:%M')    # time1_dt = datetime.strptime(time1, '%Y-%m-%d %H:%M:%S')
    time2_dt = datetime.strptime(time2, '%Y.%m.%d %H:%M')
    
    seconds_diff = (time2_dt - time1_dt).total_seconds()
    
    return abs(seconds_diff) / (3600.0 * 24 * 30 * num_mouths)    # 3600.0 * 24 * 30 * 6 = seconds for 6 mouths

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")
time = df['time'].to_list()
mag = df['mag']

t_min = time[-1]
t_max = time[0]
x0 = 0.0
x1 = seconds_between_times(t_min, t_max)
y0 = min(mag)

k = (year1 - year0) / (x1 - x0)
b = year0 - k * x0

for i in range(  len(mag)  ):
    poly_line = vtk.vtkLine()
    
    x = seconds_between_times(t_min, time[i])
    points.InsertNextPoint(x * k + b, y0, 0.0)
    points.InsertNextPoint(x * k + b, mag[i], 0.0)
    
    magnitude.InsertNextValue(mag[i])
    
    poly_line.GetPointIds().SetId(  0, 2*i      )
    poly_line.GetPointIds().SetId(  1, 2*i + 1  )
    
    cell_lines.InsertNextCell(poly_line)

linesPolyData.SetPoints(points)           # Add the points to the polydata container
linesPolyData.SetLines(cell_lines)        # Add the lines to the polydata container

# linesPolyData.GetCellData().SetScalars(magnitude)
linesPolyData.GetCellData().AddArray(magnitude)

# writer = vtk.vtkUnstructuredGridWriter()
writer = vtk.vtkXMLDataSetWriter()
writer.SetInputData(linesPolyData)
writer.SetFileName(    "M_T-" + str(num_mouths) + "_mouths.vtu"    )
# writer.SetDataModeToAscii()
writer.Update()

# Setup the visualization pipeline
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(linesPolyData)

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetLineWidth(5)

renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(namedColors.GetColor3d("SlateGray"))

window = vtk.vtkRenderWindow()
window.SetWindowName("ColoredLines")
window.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)

# Visualize
window.Render()
interactor.Start()
