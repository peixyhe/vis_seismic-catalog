####################################################################################################################
####  The program is used to visualize seismic catalog data based on RBF Kernel resampling by CPU parallelism.  ####
####  Author: He Pei; 2024.02.25                                                                                ####
####################################################################################################################



from joblib import Parallel, delayed
import math
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import vtk

# compute
## geographic region:
##                                         up(-90.0, 90.0)
##                    left(-180.0, 180.0)                     right(-180.0, 180.0)
##                                         down(-90.0, 90.0)
up = 33.4      
down = 21.8
left = 97.8   
right = 107 

step = 0.01
scale = 10                

half_setp = round(step * 0.5, 5)
ANGLE2METERS = 2 * math.pi * 6371.393 / 360.0   
z0 = scale * 0.00001 / ANGLE2METERS

e2 = (6378137.0**2 - 6356752.31414**2) / (6378137.0**2)

def BLH2XYZ(L, B):
    sinB = math.sin(math.radians(B))
    cosB = math.cos(math.radians(B))
    sinL = math.sin(math.radians(L))
    cosL = math.cos(math.radians(L))
    H = 0.0
    
    N = 6378137.0 / math.sqrt(1 - e2 * (sinB**2))    
    X = (N + H) * cosB * cosL
    Y = (N + H) * cosB * sinL
    Z = (N * (1 - e2) + H) * sinB
    
    return X / 1000.0, Y / 1000.0, Z / 1000.0

def gaussian_kernel(d, a, sigma):
    return a * np.exp(    -1.0 * (d * d) / (  2.0 * (sigma ** 2)  )    )

def process_point(yj, x, data, data_kdTree):
    print(yj)
    result = []
    for xi in x:
        xyz0 = BLH2XYZ(xi, yj)
        mag0 = 0.0
        rbf0 = 0.0
        indices = data_kdTree.query_ball_point([xi, yj], r=4.0)
        
        for index in indices:
            p1 = data[index]
            xyz1 = BLH2XYZ(p1[0], p1[1])
            sig = math.pow(2.0, p1[-1])
            distance = math.sqrt(    (xyz1[0] - xyz0[0]) ** 2 + (xyz1[1] - xyz0[1]) ** 2 + (xyz1[2] - xyz0[2]) ** 2    )
            
            if distance < 4.0 * sig:
                A = p1[-1] / 500.0
                gaussian = round(  gaussian_kernel(distance, A, sig), 5  )
            else:
                gaussian = 0.0
            rbf0 += gaussian
            
            if p1[2] > mag0:
                mag0 = p1[2]
        
        rbf0 = round(rbf0, 8)
        freq0 = len(indices)
        result.append(  [xi, yj, rbf0, freq0, mag0]  )
    return result

if __name__ == '__main__':
    df = pd.read_csv(".\\data\\mag(1_7.2)-2009_2021.csv")
    data = np.array(df[['lon', 'lat', 'mag']].values.tolist())
    data_kdTree = cKDTree(data[:, :2])

    x = [round(x0, 5) for x0 in np.arange(left, right + half_setp, step)]
    y = [round(y0, 5) for y0 in np.arange(down, up + half_setp, step)]

    results = Parallel(n_jobs=-1)(delayed(process_point)(yj, x, data, data_kdTree) for yj in y)

    with open(".\\rbf_Mcpu.csv", "w") as w:
        w.write("lon,lat,rbf_value,freq,max_mag\n")
        for result in results:
            for r in result:
                w.write(','.join(map(str, r)) + '\n')
    w.close()



    # visualize
    colors = vtk.vtkNamedColors()
    ugrid = vtk.vtkUnstructuredGrid()

    points = vtk.vtkPoints()

    rbf = vtk.vtkFloatArray()
    rbf.SetNumberOfComponents(1)
    rbf.SetName('RBF value')
    
    freq = vtk.vtkFloatArray()
    freq.SetNumberOfComponents(1)
    freq.SetName('frequency')

    max_mag = vtk.vtkFloatArray()
    max_mag.SetNumberOfComponents(1)
    max_mag.SetName('MAX magnitude')
    
    for lin in results:
        for line in lin:
            points.InsertNextPoint(    line[0], line[1], line[2]+z0    )
            rbf.InsertNextValue(line[2])
            freq.InsertNextValue(line[3])
            max_mag.InsertNextValue(line[4])
    
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
            ugrid.InsertNextCell(    vtk.VTK_TRIANGLE, 3, [id3, id0, id1]    )    # ugrid.InsertNextCell(    vtk.VTK_QUAD, 4, [id0, id1, id2, id3]    )
            ugrid.InsertNextCell(    vtk.VTK_TRIANGLE, 3, [id3, id1, id2]    )    

    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(freq)
    ugrid.GetPointData().AddArray(max_mag)
    ugrid.GetPointData().AddArray(rbf)

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(  ".\\rbf_Mcpu_step-" + str(step) + ".vtk"  )
    # writer.SetDataModeToAscii()
    writer.Update()

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
    actor.GetProperty().SetPointSize(2)

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
    
