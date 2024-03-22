###################################################################################
####  The program is used to visualize seismic catalog to sphere in 2d space.  ####
####  Author: He Pei; 2024.02.25                                               ####
###################################################################################

scale = 10
min_mag = 1.0


import vtk
import pandas as pd
import math

df = pd.read_csv("data\\mag(1_7.2)-2009_2021.csv")
lon = df['lon']
lat = df['lat']
dep = df['dep']
mag = df['mag']

ANGLE2KILOMETERS = 2 * math.pi * 6371.393 / 360.0    # EARTH_RADIUS = 6371393.0(m)
z0 = 0.05
# z0 = scale * 6.0 / ANGLE2KILOMETERS

def mag2radius(mag0):
    # return math.exp(  (mag - 4.56) / 1.96  )
    return math.pow(  2.0, (mag0 - 4.56) / 1.96  )
max_radius = mag2radius(  max(mag)  )

polydata = vtk.vtkPolyData()

vtk_points = vtk.vtkPoints()

diameter_array = vtk.vtkFloatArray()
diameter_array.SetName("diameter")

dep_array = vtk.vtkFloatArray()
dep_array.SetName("depth(km)")

mag_array = vtk.vtkFloatArray()
mag_array.SetName("magnitude")

for i in range(len(mag)):
    if mag[i] >= CONST.min_mag:
        diameter = (  mag2radius(mag[i]) / max_radius  )
        
        vtk_points.InsertNextPoint(  lon[i], lat[i], -1.0 * dep[i] / CONST.ANGLE2KILOMETERS  )
        diameter_array.InsertNextValue(diameter)
        dep_array.InsertNextValue(dep[i])
        mag_array.InsertNextValue(mag[i])

polydata.SetPoints(vtk_points)
polydata.GetPointData().SetScalars(diameter_array)    # SetScalars; not AddArray
polydata.GetPointData().AddArray(dep_array)
polydata.GetPointData().AddArray(mag_array)

# 创建vtkSphereSource
sphere_source = vtk.vtkSphereSource()
sphere_source.SetPhiResolution(10)
sphere_source.SetThetaResolution(10)

# 使用vtkGlyph3D将球体放置在点上
glyph = vtk.vtkGlyph3D()
glyph.SetInputData(polydata)
glyph.SetSourceConnection(sphere_source.GetOutputPort())
glyph.SetScaleModeToScaleByScalar()

# 创建vtkPolyDataWriter将数据写入VTK文件
writer = vtk.vtkPolyDataWriter()
writer.SetFileName(".\\result_Data\\catalog_3d.vtk")
writer.SetInputConnection(glyph.GetOutputPort())
writer.Write()

# # Create a mapper and actor
# mapper = vtk.vtkDataSetMapper()
# mapper.SetInputData(glyph.GetOutputPort())

# actor = vtk.vtkActor()
# actor.SetMapper(mapper)
# actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
# actor.GetProperty().SetPointSize(2)

# # Visualize
# renderer = vtk.vtkRenderer()
# renderWindow = vtk.vtkRenderWindow()
# renderWindow.SetWindowName('Polyhedron')
# renderWindow.AddRenderer(renderer)
# renderWindowInteractor = vtk.vtkRenderWindowInteractor()
# renderWindowInteractor.SetRenderWindow(renderWindow)

# renderer.AddActor(actor)
# renderer.SetBackground(colors.GetColor3d('Salmon'))
# renderer.ResetCamera()
# renderer.GetActiveCamera().Azimuth(30)
# renderer.GetActiveCamera().Elevation(30)
# renderWindow.Render()
# renderWindowInteractor.Start()
