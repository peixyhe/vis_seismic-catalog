#############################################################################
####  The program is used to polt M-T diagrams of seismic catalog data.  ####
####  Author: He Pei; 2024.02.25                                         ####
#############################################################################



step_mouths = 5
min_mag = 1.0

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
    
    return abs(seconds_diff) / (3600.0 * 24 * 30 * step_mouths)    # 3600.0 * 24 * 30 * 6 = seconds for 6 mouths

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")
time = df['time'].to_list()
mag = df['mag']

x0 = 0.0
x1 = seconds_between_times("2009.1.1 0:0", "2021.12.31 23:59")
# y0 = min(mag)
y0 = min_mag

year0 = 2009
year1 = 2022
k = (year1 - year0) / (x1 - x0)
print(1 / k, -1.0 / k * 2009)

j = 0
for i in range(  len(mag)  ):
    if mag[i] >= min_mag:
        poly_line = vtk.vtkLine()
        
        x = seconds_between_times(  "2009.1.1 0:0", time[i]  )
        points.InsertNextPoint(x, y0, 0.0)
        points.InsertNextPoint(x, mag[i], 0.0)
        
        magnitude.InsertNextValue(mag[i])
        
        poly_line.GetPointIds().SetId(  0, 2*j      )
        poly_line.GetPointIds().SetId(  1, 2*j + 1  )
        
        cell_lines.InsertNextCell(poly_line)
        
        j+=1

linesPolyData.SetPoints(points)           # Add the points to the polydata container
linesPolyData.SetLines(cell_lines)        # Add the lines to the polydata container

# linesPolyData.GetCellData().SetScalars(magnitude)
linesPolyData.GetCellData().AddArray(magnitude)

# writer = vtk.vtkUnstructuredGridWriter()
writer = vtk.vtkXMLDataSetWriter()
writer.SetInputData(linesPolyData)
writer.SetFileName(    "M-T_minMAG" + str(min_mag) + "_" + str(step_mouths) + "mouths.vtu"    )
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
