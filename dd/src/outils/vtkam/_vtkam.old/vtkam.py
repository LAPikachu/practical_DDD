#########################################################################################
##### This a script to generate a vkt vizualization of a microMegas microstructure  #####
##### The vtk wrapper for python must be installed                                  #####
##### vtkam.input to set-up the vizualisation                                       #####
#########################################################################################


#######################  IMPORT PACKAGES  ###############################################
import vtk
import math
import numpy as np
from pylab import *
import FortranFile as fofi
import struct
import sys

print ""
print "***********************************************************************************"
print "************************************ NEW VTKAM ************************************"
print "***********************************************************************************"
print ""

if vtk.VTK_MAJOR_VERSION <= 5:
	print "This programm has been designed for vtk >= 6. Your version is currently : ", real(vtk.VTK_MAJOR_VERSION)
	print "Please update vtk to a newer version.\n"
	quit()

try:
	stepIni= int(raw_input("Visualisation start at step:  ? "))
	perio= int(raw_input("Snapshot frequency:  ? "))
except ValueError:
	print "Cannot understand starting step or snapshot frequency"

#######################  COLORS  ########################################################
## Setup two colors - one for each line
red = [255, 0, 0]
green = [0, 255, 0]
blue = [0, 0, 255]

col={}
with open ("../in/couleur.micmeg", "r") as mycolors:
	tata = mycolors.readlines()
	nbsys=int(tata[24].strip())
	col['0']=[255,255,0]
	col['1']=[255,128,0]
	col['2']=[255,0,255]
	col['3']=[255,0,127]
	col['4']=[0,0,255]
	col['5']=[0,128,255]
	col['6']=[0,204,0]
	col['7']=[51,255,153]
	col['8']=[51,102,0]
	col['9']=[0,76,153]
	col['10']=[0,255,255]
	col['11']=[255,153,153]
	col['-1']=[255,255,255] #Junc
	col['-2']=[128,128,128] #GD

## Setup the colors array
colors = vtk.vtkUnsignedCharArray()
colors.SetNumberOfComponents(3)
colors.SetName("Colors")

#######################  CAMERA AND RENDERER  ########################################################
## Create a Camera
camera =vtk.vtkCamera()

## Create a Renderer
ren1 = vtk.vtkRenderer()
ren1.SetActiveCamera(camera); ##Activate camera

renWin = vtk.vtkRenderWindow()

#######################  READ INPUT  #################################################
## VTKAM config
myConfDic={} ##Initialize Dictionary
TrueList=['true','True','TRUE', '1', 't', 'T', 'YES', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']

with open ("../in/vtkam/vtkam.conf", "r") as myconffile:
    myconffile.readline()
    tmp=myconffile.readline().rstrip()
    cline=tmp.split('=')
    while(cline[0]!='[eof]'):
        myConfDic[cline[0]]=cline[1]
        tmp=myconffile.readline().rstrip()
        cline=tmp.split('=')

## Load input.dd
with open ("../in/input.dd", "r") as myinputdd:
    line=myinputdd.readline()
    while (line[0]!='-'):
        line=myinputdd.readline()
    MaterialFileName = myinputdd.readline().rstrip().split()[0] #rstrip used to remove "\n" from string
    ControlFileName  = myinputdd.readline().rstrip().split()[0]
    SegmentFileName  = myinputdd.readline().rstrip().split()[0]

    print "Material     : ", MaterialFileName
    print "Control file : ", ControlFileName
    print "Segment file : ", SegmentFileName

    MaterialFilePath="../in/"+MaterialFileName
    ControlFilePath="../in/"+ControlFileName
    SegmentFilePath="../in/"+SegmentFileName

with open (SegmentFilePath, "r") as mysegfile:
    mysegfile.readline()
    mysegfile.readline()
    boxSize=mysegfile.readline().rstrip().split()
    print "\nSimulation box size : ", boxSize[0],boxSize[1],boxSize[2],"\n"


#######################  BOXES OUTLINE  ################################################
## Main box
box = vtk.vtkCubeSource()
box.SetBounds(0,float(boxSize[0]),0,float(boxSize[1]),0,float(boxSize[2]))

## outline
outline = vtk.vtkOutlineFilter()

#VTK OUPUT
outline.SetInputConnection(box.GetOutputPort())

BoxMapper = vtk.vtkPolyDataMapper()

#VTK OUPUT
BoxMapper.SetInputConnection(outline.GetOutputPort())

BoxActor = vtk.vtkActor()
BoxActor.SetMapper(BoxMapper)
BoxActor.GetProperty().SetColor(red)
ren1.AddActor(BoxActor)

if (myConfDic['EXTRABOX'] in TrueList ):
    ## Secondary box
    extraBox = vtk.vtkCubeSource()
    extraBox.SetBounds(0,200000,102500,307500,95000,285000)
    
    ## outline
    extraOutline = vtk.vtkOutlineFilter()
    
    #VTK OUPUT
    extraOutline.SetInputConnection(extraBox.GetOutputPort())
    
    ExtraBoxMapper = vtk.vtkPolyDataMapper()
    
    #VTK OUPUT
    ExtraBoxMapper.SetInputConnection(extraOutline.GetOutputPort())
    
    ExtraBoxActor = vtk.vtkActor()
    ExtraBoxActor.SetMapper(ExtraBoxMapper)
    #ExtraBoxActor.GetProperty().SetColor(red)
    ren1.AddActor(ExtraBoxActor)

#######################  PLANES  ########################################################
if (myConfDic['PLANES'] in TrueList ):
    with open ("../in/vtkam/vtkam.planes", "r") as myplanes:
        firstline=myplanes.readline()
        if (firstline==""):
            raise NameError('Hi There! vtk.planes seems to be empty !!')
            quit()
        tmp=myplanes.readline()
        pline=tmp.split()
        while(pline[0]!='[eof]'):
            ## Create Plane
            myplane = vtk.vtkPlaneSource()
            #myplane.SetNormal(0,0,10)
            tmp=myplanes.readline()
            pline=tmp.split()
            
            ## Origin
            myplane.SetOrigin(float(pline[1]),float(pline[2]),float(pline[3]))
            tmp=myplanes.readline()
            pline=tmp.split()
            
            ## Point1
            myplane.SetPoint1(float(pline[1]),float(pline[2]),float(pline[3]))
            tmp=myplanes.readline()
            pline=tmp.split()
            
            ## Point2
            myplane.SetPoint2(float(pline[1]),float(pline[2]),float(pline[3]))
            tmp=myplanes.readline()
            pline=tmp.split()
            
            ## Mapper
            PlaneMapper = vtk.vtkPolyDataMapper()
            
            #VTK OUPUT
            PlaneMapper.SetInputConnection(myplane.GetOutputPort())
            
            ## Actor
            PlaneActor = vtk.vtkActor()
            PlaneActor.SetMapper(PlaneMapper)
            PlaneActor.GetProperty().SetOpacity(1)
            ren1.AddActor(PlaneActor)

#######################  READ SEG FILE  #################################################

if (myConfDic['SEGCONF'] in TrueList ):
  ## Read the segment configuration except first line
  print "Reading segment configuration from vtkam.segconf"
  with open ("../in/vtkam/vtkam.segconf", "r") as myseg:
    firstline=myseg.readline()
    segdata=myseg.readlines()

  nblines=len(segdata) ## number of line of interrest in the seg file
  ## vtkPoints represents 3D points. The data model for vtkPoints is an array of
  ## vx-vy-vz triplets accessible by (point or cell) id.

  points = vtk.vtkPoints()
  points.SetNumberOfPoints(2*nblines)

  ## vtkCellArray is a supporting object that explicitly represents cell connectivity.
  ## The cell array structure is a raw integer list of the form:
  ## (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
  ## the cell, and id is a zero-offset index into an associated point list.
  ## Create an array for line storage
  lines = vtk.vtkCellArray()

  ## Render Segments and Kneecaps
  for i in range(0,nblines,1):
	  tmp=segdata[i].split() ##Create a list of strings for the considered line
	  ## Kneecaps
	  if (float(tmp[5])==0):
		if (myConfDic['KNEECAPS'] in TrueList ):
			sphere = vtk.vtkSphereSource() ## Create a sphere for a kneecap
			sphere.SetCenter(float(tmp[2]),float(tmp[3]),float(tmp[4])) ## Set kneecap position
			sphere.SetRadius(400)
			#sphere.SetRadius(int(myConfDic['KNEECAP_RADIUS'])) ## Set kneecap radius
			
			## Create a mapper for the kneecap
			kneecapMapper = vtk.vtkPolyDataMapper()

			## Add kneecap to mapper
			
			#VTK OUPUT
			kneecapMapper.SetInputConnection(sphere.GetOutputPort())

			## Create kneecap Actor
			kneecapActor = vtk.vtkActor()

			## Set color for all kneecap
			#kneecapActor.GetProperty().SetColor(0.0, 0.0, 0.0)

			## Add Mapper to actor and actor to renderer
			kneecapActor.SetMapper(kneecapMapper)
			ren1.AddActor(kneecapActor)
      ## Segments
	  else:
		points.SetPoint(2*i, float(tmp[2]),float(tmp[3]),float(tmp[4]))
		points.SetPoint(2*i+1, float(tmp[6]),float(tmp[7]),float(tmp[8]))

		## Create Segment
		segment = vtk.vtkLine()
		segment.GetPointIds().SetId(0,2*i)      ## Point 1
		segment.GetPointIds().SetId(1,2*i+1)    ## Point 2

		## Add segment in the array

		lines.InsertNextCell(segment)
		
		## Set-up segment color
		if str(tmp[14]) is 'T':
			colors.InsertNextTupleValue(col['-1'])
		elif str(tmp[15]) is 'T':
			colors.InsertNextTupleValue(col['-2'])
		else:
			for i in range(0,nbsys):
				if (int(tmp[1]) in range(8*i+1,8*i+9)):
					colors.InsertNextTupleValue(col[str(i)])
  ## vtkPolyDataMapper is a class that maps polygonal data (i.e., vtkPolyData)
  ## to graphics primitives

  ## vtkPolyData is a data object that is a concrete implementation of vtkDataSet.
  ## vtkPolyData represents a geometric structure consisting of vertices, lines,
  ## polygons, and/or triangle strips

  ## Create PolyData object and add points, segments and color
  polyseg = vtk.vtkPolyData()
  polyseg.SetPoints(points)
  polyseg.SetLines(lines)
  polyseg.GetCellData().SetScalars(colors)
  ## Create a mapper for the segments
  polysegMapper = vtk.vtkPolyDataMapper()
	  
  ## Add segments to mapper
  #Vtk OUTPUT
  polysegMapper.SetInputData(polyseg)
  polysegMapper.Update()

  if (myConfDic['TUBE_FILTER'] in TrueList ):
	  ## Tube filter for rendering
	  tubes = vtk.vtkTubeFilter()
	  tubes.SetInputData(polyseg)
	  tubes.SetRadius(int(myConfDic['TUBE_RADIUS']))
	  tubes.SetNumberOfSides(int(myConfDic['TUBE_NSIDES']))
	  polysegMapper.SetInputConnection(tubes.GetOutputPort())
	  polysegMapper.Update()
  # Create an actor to represent the polygon. The actor orchestrates rendering of
  # the mapper's graphics primitives. An actor also refers to properties via a
  # vtkProperty instance, and includes an internal transformation matrix.

  ## Actor
  polysegActor = vtk.vtkActor()
  #polygonActor.GetProperty().SetLineWidth(100) ## Line width if no tube filter (flat line)

  ## Mapper
  polysegActor.SetMapper(polysegMapper)

  # Create the Renderer and assign actors to it. A renderer is like a
  # viewport. It is part or all of a window on the screen and it is
  # responsible for drawing the actors it has.  We also set the
  # background color here.

  ## Add actor to renderer
  ren1.AddActor(polysegActor)

  #######################  RENDER  ########################################################

  # Automatically set up the camera based on the visible actors.
  # The camera will reposition itself to view the center point of the actors,
  # and move along its initial view plane normal
  # (i.e., vector defined from camera position to focal point) so that all of the
  # actors can be seen.

  ren1.SetBackground(0.0, 0.0, 0.0)
  ren1.ResetCamera(0,float(boxSize[0]),0,float(boxSize[1]),0,float(boxSize[2]))
  #camera.SetFocalPoint(int(boxSize[0])/2,int(boxSize[1])/2,int(boxSize[2])/2)
  camera.GetViewUp()
  #camera.SetViewUp(0.0,0.0,0.0)
  camera.SetFocalPoint(44796, 153728, 191980)

  # Finally we create the render window which will show up on the screen
  # We put our renderer into the render window using AddRenderer. We
  # also set the size to be 500 pixels by 500.

  renWin.AddRenderer(ren1)
  renWin.SetSize(int(myConfDic['WINDOW_SIZE']),int(myConfDic['WINDOW_SIZE']))
  renWin.Render()

  # The vtkRenderWindowInteractor class watches for events (e.g., keypress,
  # mouse) in the vtkRenderWindow. These events are translated into
  # event invocations that VTK understands (see VTK/Common/vtkCommand.h
  # for all events that VTK processes). Then observers of these VTK
  # events can process them as appropriate.
  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)
  iren.Initialize()

  raw_input("")

else:
	## Read film.bin
	f = fofi.FortranFile('../out/film.bin')
	crystalbase={'1':'CS','2':'BCC','3':'CFC','4':'HCP','5':'ORT','6':'MGO','7':'DC'}
	#byte = myfilm.read()
	cristallo   =f.readInts()
	nsegmaxf    =f.readInts()
	avalue      =f.readReals('d')
	nbase       =f.readInts()
	bvec        =f.readInts()
	modur       =f.readInts()

	if ( cristallo>7 or cristallo<1):
		raise NameError('The crystal symmetry is unknown')
	else:
		print "The crystal symmetry is of type",crystalbase[str(cristallo)[1]]

	#print "nsegmaxf",int(nsegmaxf)
	print "avalue",float(avalue)
	#print "nbase",int(nbase)
        print " "

    ## Initialize segment line and normal vector as dictionaries
	bveclin={}
	bvecnor={}
		
    ## Read segment line and normal vector from bvec
	for i in range(0,nbase):
		bveclin[str(i)]=bvec[6*i  :6*i+3]
		bvecnor[str(i)]=bvec[6*i+3:6*i+6]
    
	kk=1
	while (kk>0):
		try:
			data_str = f.readString()
            #Barbarism
			kk=int(np.array(struct.unpack(str(len(data_str[0:4])/struct.calcsize('i'))+'i',data_str[0:4])))
			deltat=float(np.array(struct.unpack(str(len(data_str[4:12])/struct.calcsize('d'))+'d',data_str[4:12])))
			NSEGM=int(np.array(struct.unpack(str(len(data_str[12:16])/struct.calcsize('i'))+'i',data_str[12:16])))
			XSTRESS=float(np.array(struct.unpack(str(len(data_str[16:24])/struct.calcsize('d'))+'d',data_str[16:24])))
			EPSO=float(np.array(struct.unpack(str(len(data_str[24:32])/struct.calcsize('d'))+'d',data_str[24:32])))
            ##
			print "{0:<6}, time {1:6.2e} Eps: {2:2.3f} Sigma: {3:6.2f} MPa Nseg: {4:<8}".format(kk, deltat*kk , EPSO , XSTRESS , NSEGM )
			segdata= f.readInts('i')
			f.readString()

		except IOError as e:
			print ""
			print "***********************************************************************************"
			print "**************************** End of the film.bin file *****************************"
			print "***********************************************************************************"
			print ""
			quit()
		if (kk>=int(stepIni) and kk%int(perio)==0):
			## Render Segments and Kneecaps
			points = vtk.vtkPoints()
			points.SetNumberOfPoints(2*NSEGM)
	 
			## vtkCellArray is a supporting object that explicitly represents cell connectivity.
			## The cell array structure is a raw integer list of the form:
			## (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
			## the cell, and id is a zero-offset index into an associated point list.
		  
			## Create an array for line storage
			lines = vtk.vtkCellArray()
			for i in range(0,NSEGM,1):
				tmp=segdata[9*i:9*i+9]#.split() ##Create a list of strings for the considered line
			
				## Kneecaps
				if (int(tmp[3])==0):
					if (myConfDic['KNEECAPS'] in TrueList ):
						sphere = vtk.vtkSphereSource() ## Create a sphere for a kneecap
						sphere.SetCenter(int(tmp[0]),int(tmp[1]),int(tmp[2])) ## Set kneecap position
						sphere.SetRadius(int(myConfDic['KNEECAP_RADIUS'])) ## Set kneecap radius
					
						## Create a mapper for the kneecap
						kneecapMapper = vtk.vtkPolyDataMapper()
					
						## Add kneecap to mapper
						kneecapMapper.SetInputConnection(sphere.GetOutputPort())
						
						## Create kneecap Actor
						kneecapActor = vtk.vtkActor()
						
						## Set color for all kneecap
						#kneecapActor.GetProperty().SetColor(0.0, 0.0, 0.0)
						
						## Add Mapper to actor and actor to renderer
						kneecapActor.SetMapper(kneecapMapper)
						ren1.AddActor(kneecapActor)
				  ## Segments
				else:
					#print tmp
					O1=int(tmp[0])
					O2=int(tmp[1])
					O3=int(tmp[2])
					E1=int(tmp[0])+int(tmp[3])*int(bveclin[str(tmp[4]-1)][0])
					E2=int(tmp[1])+int(tmp[3])*int(bveclin[str(tmp[4]-1)][1])
					E3=int(tmp[2])+int(tmp[3])*int(bveclin[str(tmp[4]-1)][2])
					
					points.SetPoint(2*i, O1,O2,O3)
					points.SetPoint(2*i+1, E1,E2,E3)
					
					## Create Segment
					segment = vtk.vtkLine()
					segment.GetPointIds().SetId(0,2*i)      ## Point 1
					segment.GetPointIds().SetId(1,2*i+1)    ## Point 2
						
					## Add segment in the array
					lines.InsertNextCell(segment)
					
					## Set-up segment color
					if int(tmp[5])>0:
						colors.InsertNextTupleValue(col['-1'])
					elif int(tmp[5])<0:
						colors.InsertNextTupleValue(col['-2'])
					else:
						for i in range(0,nbsys):
							if (int(tmp[4]) in range(8*i+1,8*i+9)):
								colors.InsertNextTupleValue(col[str(i)])

			## vtkPolyData is a data object that is a concrete implementation of vtkDataSet.
			## vtkPolyData represents a geometric structure consisting of vertices, lines,
			## polygons, and/or triangle strips

			## Create PolyData object and add points, segments and color
			polyseg = vtk.vtkPolyData()
			polyseg.SetPoints(points)
			polyseg.SetLines(lines)
			polyseg.GetCellData().SetScalars(colors)

			## vtkPolyDataMapper is a class that maps polygonal data (i.e., vtkPolyData)
			## to graphics primitives

			## Create a mapper for the segments
			polysegMapper = vtk.vtkPolyDataMapper()

			## Add segments to mapper
			polysegMapper.SetInputData(polyseg)
			polysegMapper.Update()
			if (myConfDic['TUBE_FILTER'] in TrueList ):
				## Tube filter for rendering
				tubes = vtk.vtkTubeFilter()
				tubes.SetInputData(polyseg)
				tubes.SetRadius(int(myConfDic['TUBE_RADIUS']))
				tubes.SetNumberOfSides(int(myConfDic['TUBE_NSIDES']))
				polysegMapper.SetInputConnection(tubes.GetOutputPort())
				polysegMapper.Update()

			# Create an actor to represent the polygon. The actor orchestrates rendering of
			# the mapper's graphics primitives. An actor also refers to properties via a
			# vtkProperty instance, and includes an internal transformation matrix.

			## Actor
			polysegActor = vtk.vtkActor()
			#polygonActor.GetProperty().SetLineWidth(100) ## Line width if no tube filter (flat line)

			## Mapper
			polysegActor.SetMapper(polysegMapper)

			# Create the Renderer and assign actors to it. A renderer is like a
			# viewport. It is part or all of a window on the screen and it is
			# responsible for drawing the actors it has.  We also set the
			# background color here.

			## Add actor to renderer
			ren1.AddActor(polysegActor)

			#######################  RENDER  ########################################################

			# Automatically set up the camera based on the visible actors.
			# The camera will reposition itself to view the center point of the actors,
			# and move along its initial view plane normal
			# (i.e., vector defined from camera position to focal point) so that all of the
			# actors can be seen.
			
			ren1.SetBackground(0.0, 0.0, 0.0)
			#camera.SetFocalPoint(int(boxSize[0])/2,int(boxSize[1])/2,int(boxSize[2])/2)
			camera.GetViewUp()
			#camera.SetViewUp(0.0,0.0,0.0)
			ren1.ResetCamera(0,float(boxSize[0]),0,float(boxSize[1]),0,float(boxSize[2]))
			#camera.SetFocalPoint(44796, 153728, 191980)

			# Finally we create the render window which will show up on the screen
			# We put our renderer into the render window using AddRenderer. We
			# also set the size to be 500 pixels by 500.
			
			renWin.AddRenderer(ren1)
			renWin.SetSize(int(myConfDic['WINDOW_SIZE']),int(myConfDic['WINDOW_SIZE']))
			renWin.Render()

			# The vtkRenderWindowInteractor class watches for events (e.g., keypress,
			# mouse) in the vtkRenderWindow. These events are translated into
			# event invocations that VTK understands (see VTK/Common/vtkCommand.h
			# for all events that VTK processes). Then observers of these VTK
			# events can process them as appropriate.
			iren = vtk.vtkRenderWindowInteractor()
			iren.SetRenderWindow(renWin)
			iren.Initialize()
			#iren.Start()
			raw_input("")
			ren1.RemoveActor(polysegActor)
			points.Initialize()
			lines.Initialize()
			segment.Initialize()
			polyseg.Initialize()
			colors.Initialize()
			polysegMapper.RemoveAllVertexAttributeMappings()
			if (myConfDic['KNEECAPS'] in TrueList ):
				kneecapMapper.RemoveAllVertexAttributeMappings()
				ren1.RemoveActor(kneecapActor)
			# There is no explicit need to free any objects at this point.
			# Once Python exits, memory is automatically freed.
