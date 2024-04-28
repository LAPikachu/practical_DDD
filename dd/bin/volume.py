#!/usr/bin/env python

def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

import numpy as np
from scipy.spatial import Delaunay
import vtk

list=np.loadtxt('../out/vertexlist')

Nbcvxdom=list[0,0]

limit=list[1,0]
limitold=0

pointsvtk = vtk.vtkPoints()
tetra = vtk.vtkTetra()
cellArray = vtk.vtkCellArray()
pointstot=0
volume=0

for k in range(0,int(Nbcvxdom)):


  limitold=limitold+1
  limit=list[limitold,0]

  points=np.abs(np.array(list[int(limitold+1):int(limitold+1+limit),:]))
  limitold=limitold+limit


  tri = Delaunay(points)

  if hasattr(tri, 'simplices'):
    simplices=tri.simplices
  else:
    simplices=tri.vertices

  tets = tri.points[simplices]

  vol = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                  tets[:, 2], tets[:, 3]))


  volume=volume+vol



  spnt=np.shape(points)


  for i in range(0,spnt[0]):
    pointsvtk.InsertNextPoint(points[i,:])



  unstructuredGrid = vtk.vtkUnstructuredGrid()
  unstructuredGrid.SetPoints(pointsvtk)

  #tetra = vtk.vtkTetra()
  #cellArray = vtk.vtkCellArray()

  for i in range(0,tri.nsimplex):
    tetra.GetPointIds().SetId(0, simplices[i,0]+int(pointstot))
    tetra.GetPointIds().SetId(1, simplices[i,1]+int(pointstot))
    tetra.GetPointIds().SetId(2, simplices[i,2]+int(pointstot))
    tetra.GetPointIds().SetId(3, simplices[i,3]+int(pointstot))
    cellArray.InsertNextCell(tetra)

  pointstot=pointstot+limit

unstructuredGrid.SetCells(vtk.VTK_TETRA, cellArray)

writer=vtk.vtkUnstructuredGridWriter()
writer.SetFileName('../out/simulationvolume.vtk')

if vtk.VTK_MAJOR_VERSION < 6 :
  writer.SetInput(unstructuredGrid)
else :
  writer.SetInputData(unstructuredGrid)

writer.Write()

f=open('../out/volumepy.txt','w')

f.write("%10.2f" %volume)

f.close

# mapper1 = vtk.vtkDataSetMapper()
# if vtk.VTK_MAJOR_VERSION <= 5:
#     mapper1.SetInputConnection(unstructuredGrid.GetProducerPort())
# else:
#     mapper1.SetInputData(unstructuredGrid)
#
# actor1 = vtk.vtkActor()
# actor1.SetMapper(mapper1)
#
# # Create a renderer, render window, and interactor
# renderer = vtk.vtkRenderer()
# renderWindow = vtk.vtkRenderWindow()
# renderWindow.AddRenderer(renderer)
# renderWindowInteractor = vtk.vtkRenderWindowInteractor()
# renderWindowInteractor.SetRenderWindow(renderWindow)
#
# # Add the actor to the scene
# renderer.AddActor(actor1)
# renderer.SetBackground(.3, .6, .3) # Background color green

# Render and interact
#renderWindow.Render()
#renderWindowInteractor.Start()







# pts = np.random.rand(10, 3)
# dt = Delaunay(pts)
# tets = dt.points[dt.simplices]
# vol = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
#                                 tets[:, 2], tets[:, 3]))
