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

import sys
import datetime

import menu
import utilities as util

import fsrc

class Vtkam(object):
  """
  """
  print ""

  def __init__(self):
    """
    """
    print "***********************************************************************************"
    print "************************************ NEW VTKAM ************************************"
    print "***********************************************************************************"
    print ""

    # Check if installed vtk is compatible
    self.check_version()
    # Define vizualisation options
    menu.define_options(self)
    # Read input files
    self.read_input()
    # Colors for rendering
    self.set_colors()
    # Ask for user inputs
    self.user_input()
    # Print informations
    self.info()

    self.wait = True # Wait for user input

  ### CAMERA AND RENDERER
    ## Create a Camera
    self.camera =vtk.vtkCamera()

    ## Create a Renderer
    self.ren1 = vtk.vtkRenderer()
    self.ren1.SetActiveCamera(self.camera) ##Activate camera

  ### ORIGIN
    sphere0 = vtk.vtkSphereSource()

    sphere0.SetRadius(np.max([ float(self.boxSize[i]) for i in range(3)] )/200) ## Set kneecap radius for the origin
    sphere0.SetCenter(0,0,0) ## Set origin position

    map0 = vtk.vtkPolyDataMapper()
    map0.SetInputConnection(sphere0.GetOutputPort())

    self.actor0 = vtk.vtkActor()
    self.actor0.SetMapper(map0)

    self.ren1.AddActor(self.actor0)

  ### AXES
    self.AxesActor = vtk.vtkAxesActor()
    axelength = np.min([ float(self.boxSize[i]) for i in range(3)] )/2
    self.AxesActor.SetTotalLength(axelength,axelength,axelength)

    self.ren1.AddActor(self.AxesActor)

  ### BOXES OUTLINE
    ## Main box
    self.box = vtk.vtkCubeSource()
    self.box.SetBounds(0,float(self.boxSize[0]),0,float(self.boxSize[1]),0,float(self.boxSize[2]))

    ## outline
    outline = vtk.vtkOutlineFilter()

    #VTK OUPUT
    outline.SetInputConnection(self.box.GetOutputPort())

    self.BoxMapper = vtk.vtkPolyDataMapper()

    #VTK OUPUT
    self.BoxMapper.SetInputConnection(outline.GetOutputPort())

    self.BoxActor = vtk.vtkActor()
    self.BoxActor.SetMapper(self.BoxMapper)
    self.BoxActor.GetProperty().SetColor(self.red)

    self.ren1.AddActor(self.BoxActor)

#    if (self.myConfDic['PLANES'] in self.TrueList ):
#      "The 'self.' has not be written"
#      self.load_vtkplanes()

  ### WINDOW
    ## Create a Window
    self.renWin = vtk.vtkRenderWindow()

    self.ren1.SetBackground(0.0, 0.0, 0.0)

    self.ren1.ResetCamera()

    self.dist = int(self.camera.GetDistance())

    self.camera.ParallelProjectionOn()

    self.camera.SetFocalPoint(0.5*int(self.boxSize[0]),0.5*int(self.boxSize[1]),0.5*int(self.boxSize[2]))
    self.camera.SetPosition(self.dist,0.5*int(self.boxSize[1]),0.5*int(self.boxSize[2]))
    self.camera.SetViewUp(0,0,1)

    print "   Focal point          : ", self.camera.GetFocalPoint()
    print "   Camera position      : ", self.camera.GetPosition()
    print "   Camera distance      : ", self.camera.GetDistance()
    print "   Projection direction : ", self.camera.GetDirectionOfProjection()
    print "   View up vector       : ", self.camera.GetViewUp()

    self.renWin.AddRenderer(self.ren1)

    self.renWin.SetSize(int(self.myConfDic['WINDOW_SIZE']),int(self.myConfDic['WINDOW_SIZE']))

    self.renWin.Render()

  ### INTERACTOR
    iren = vtk.vtkRenderWindowInteractor()

    iren.SetRenderWindow(self.renWin)
    iren.Initialize()


  ### READ AND UPDATE
    if (self.myConfDic['SEGCONF'] in self.TrueList ):
      print ">> NOT WORKING YET << "
      exit
      self.load_filmheader()
      self.load_conf()

      self.generate_vtkdata()
      self.fill_vtkdata()

      self.ren1.ResetCamera()
      self.render()

      Loop = True
      while (Loop): # While choice in menu
        Loop = False
        self.render()
        self.userinput = raw_input(" !> ")
        Loop = menu.apply_option(self)

    else :

      self.load_filmheader()

      #while (not stop)
      self.userinput = 0
      self.kk=self.stepIni

      self.generate_vtkdata()

      print "\n> Type 'm' to display the menu"

      while (self.userinput != 'z') : # while not exit

        self.load_filmstep(self.kk)        # Load Film step

        Loop = True
        kkold = self.kk
        while (Loop): # While choice in menu
          Loop = False
          self.generate_vtkdata()

          self.fill_vtkdata()

          self.render()

          if (self.wait) : self.userinput = raw_input(" !> ")

          self.initialize_vtkobjects() # Clean all

          Loop = menu.apply_option(self)

        if (self.kk == kkold) :
          self.kk = self.kk + self.perio

  def __str__(self):
    string = "VTKAM"
    return string

  def __repr__(self):
    return

#  def __setattr__(self, name, value):
#    if name not in self.__dict__ :
#      print " "*6+"\n!> You can not add new attributes in this class <!"
#    else :
#      super(Vtkam, self).__setattr__(name, value)
#
#  def __delattr__(self, name):
#    print " "*6+"\n!> You can not delete attributes in this class <!"
#    pass

  def info(self):
    """
    """
    print "\nInformations : "
    print "---------------"
    print " Material     : ", self.MaterialFileName
    print " Control file : ", self.ControlFileName
    print " Segment file : ", self.SegmentFileName
    print " Simulation box size : ", self.boxSize[0],self.boxSize[1],self.boxSize[2], "\n"

  def check_version(self):
    """
    """
    #Check updates
    if vtk.VTK_MAJOR_VERSION <= 5:
      print "This programm has been designed for vtk >= 6. Your version is currently : ", real(vtk.VTK_MAJOR_VERSION)
      print "Please update vtk to a newer version.\n"
      quit()

  def ReadVTKSurfaceFile(self, filepath):
    if (filepath == ''):
        print '!> Error: no filepath.'
    print('> Reading VTK surface file.')
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filepath)
    reader.Update()
    surface = reader.GetOutput()

    # create the mapper that corresponds the objects of the vtk file into graphics elements
    surface_mapper = vtk.vtkDataSetMapper()
    surface_mapper.SetInputData(surface)

    # actor
    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)

    return surface_actor


  def user_input(self):
    """
    """
    while True:
      try:
        self.stepIni = util.get_non_negative_int("\n> Visualisation start at step:  ? ")
        break
      except :
        print "\n !> Cannot understand starting step - Exiting"
        quit()

    while True:
      try:
        self.perio   = util.get_non_negative_int("> Snapshot frequency:  ? ")
        break
      except :
        print "\n !> Cannot understand snapshot frequency - Exiting"
        quit()

  def read_input(self):
    """
    """
    #######################  READ INPUT  #################################################
    ## VTKAM config
    self.myConfDic={} ##Initialize Dictionary
    self.TrueList=['true','True','TRUE', '1', 't', 'T', 'YES', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']

    print "> Reading Vtkam configuration file...",
    with open ("../in/vtkam/vtkam.conf", "r") as myconffile:
        myconffile.readline()
        tmp=myconffile.readline().rstrip()
        cline=tmp.split('=')
        while(cline[0]!='[eof]'):
            self.myConfDic[cline[0]]=cline[1]
            tmp=myconffile.readline().rstrip()
            cline=tmp.split('=')
    print "[ok]"

    ## Load input.dd
    print "> Reading mM input file...",
    with open ("../in/input.dd", "r") as myinputdd:
        line=myinputdd.readline()
        while (line[0]!='-'):
            line=myinputdd.readline()
        self.MaterialFileName = myinputdd.readline().rstrip().split()[0] #rstrip used to remove "\n" from string
        self.ControlFileName  = myinputdd.readline().rstrip().split()[0]
        self.SegmentFileName  = myinputdd.readline().rstrip().split()[0]

        MaterialFilePath="../in/"+self.MaterialFileName
        ControlFilePath="../in/"+self.ControlFileName
        SegmentFilePath="../in/"+self.SegmentFileName

        print "[ok]"


    print "> Reading segment file...",
    with open (SegmentFilePath, "r") as mysegfile:
        mysegfile.readline()
        mysegfile.readline()
        self.boxSize=mysegfile.readline().rstrip().split()
        print "[ok]"

    #Reading Mode
    self.mode         = self.myConfDic['READING_MODE']

    #Number of Replica
    self.nrep=[np.int(self.myConfDic['NREP1']),np.int(self.myConfDic['NREP2']),np.int(self.myConfDic['NREP3'])]

    #No CLipping at start
    self.myConfDic['CLIPPING'] = 'F'
    self.origin = [float(self.boxSize[0])/2.,float(self.boxSize[1])/2.,float(self.boxSize[2])/2.]
    self.dx = np.float(self.boxSize[0])/2.
    self.dy = np.float(self.boxSize[1])/2.
    self.dz = np.float(self.boxSize[2])/2.

    self.normal = [None,None,None]
    self.surfaceActor = None

  def set_colors(self):
    """
    """
    print "> Defining colors...",
    #######################  COLORS  ########################################################
    ## Setup two colors - one for each line
    self.red = [255, 0, 0]
    self.green = [0, 255, 0]
    self.blue = [0, 0, 255]
    self.gray = [64,64,64]
    self.white = [255,255,255]

    self.col={}

    if (np.int(self.myConfDic['COLORS'])== 0) :

      self.col['0']=[255,222,0]
      self.col['1']=[255,140,0]
      self.col['2']=[0,153,153]
      self.col['3']=[0,0,255]
      self.col['4']=[51,209,51]
      self.col['5']=[0,102,0]
      self.col['6']=[0,204,204]
      self.col['7']=[0,0,127]
      self.col['8']=[224,114,150]
      self.col['9']=[127,0,127]
      self.col['10']=[255,0,0]
      self.col['11']=[102,0,0]

      self.col['12']=self.white #unknown ?

      self.col['13']=self.gray # [0,204,204] # edge Lomer-Cotrell lock
      self.col['14']=self.gray #[0,0,255] # 60 degree Lomer-Cotrell lock
      self.col['15']=self.gray #[51,209,51] # edge glissile
      self.col['16']=self.gray #[0,102,0] # 60 degree glissile
      self.col['17']=self.gray #[255,102,0] # Hirth lock
      self.col['18']=self.gray #[212,0,0] # Collinear junction
      self.col['19']=self.gray #[127,0,127] # GD segments

    elif (np.int(self.myConfDic['COLORS'])== 1) :

      self.col['0']=self.gray
      self.col['1']=self.gray
      self.col['2']=self.gray
      self.col['3']=self.gray
      self.col['4']=self.gray
      self.col['5']=self.gray
      self.col['6']=self.gray
      self.col['7']=self.gray
      self.col['8']=self.gray
      self.col['9']=self.gray
      self.col['10']=self.gray
      self.col['11']=self.gray

      self.col['12']=self.white #unknown ?

      self.col['13']=[0,204,204] # edge Lomer-Cotrell lock
      self.col['14']=[0,0,255] # 60 degree Lomer-Cotrell lock
      self.col['15']=[51,209,51] # edge glissile
      self.col['16']=[0,102,0] # 60 degree glissile
      self.col['17']=[255,102,0] # Hirth lock
      self.col['18']=[212,0,0] # Collinear junction
      self.col['19']=[127,0,127] # GD segments

    print "[ok]"

  def set_extrabox(self):
    """
    """
    ## Secondary box
    extraBox = vtk.vtkCubeSource()
    xmin = util.get_float(" > xmin : ")
    xmax = util.get_float(" > xmax : ")
    ymin = util.get_float(" > ymin : ")
    ymax = util.get_float(" > ymax : ")
    zmin = util.get_float(" > zmin : ")
    zmax = util.get_float(" > zmax : ")

    extraBox.SetBounds(xmin,xmax,ymin,ymax,zmin,zmax)

    ## outline
    extraOutline = vtk.vtkOutlineFilter()

    #VTK OUPUT
    extraOutline.SetInputConnection(extraBox.GetOutputPort())

    ExtraBoxMapper = vtk.vtkPolyDataMapper()

    #VTK OUPUT
    ExtraBoxMapper.SetInputConnection(extraOutline.GetOutputPort())

    self.ExtraBoxActor = vtk.vtkActor()
    self.ExtraBoxActor.SetMapper(ExtraBoxMapper)
    #ExtraBoxActor.GetProperty().SetColor(red)
    self.ren1.AddActor(self.ExtraBoxActor)

  def del_extrabox(self):
    """
    """
    self.ren1.RemoveActor(self.ExtraBoxActor)

#  def load_vtkplanes(self):
#    """
#    """
#    ###  PLANES
#    with open ("../in/vtkam/vtkam.planes", "r") as myplanes:
#        firstline=myplanes.readline()
#        if (firstline==""):
#            raise NameError('Hi There! vtk.planes seems to be empty !!')
#            quit()
#        tmp=myplanes.readline()
#        pline=tmp.split()
#        while(pline[0]!='[eof]'):
#            ## Create Plane
#            myplane = vtk.vtkPlaneSource()
#            #myplane.SetNormal(0,0,10)
#            tmp=myplanes.readline()
#            pline=tmp.split()
#
#            ## Origin
#            myplane.SetOrigin(float(pline[1]),float(pline[2]),float(pline[3]))
#            tmp=myplanes.readline()
#            pline=tmp.split()
#
#            ## Point1
#            myplane.SetPoint1(float(pline[1]),float(pline[2]),float(pline[3]))
#            tmp=myplanes.readline()
#            pline=tmp.split()
#
#            ## Point2
#            myplane.SetPoint2(float(pline[1]),float(pline[2]),float(pline[3]))
#            tmp=myplanes.readline()
#            pline=tmp.split()
#
#            ## Mapper
#            PlaneMapper = vtk.vtkPolyDataMapper()
#
#            #VTK OUPUT
#            PlaneMapper.SetInputConnection(myplane.GetOutputPort())
#
#            ## Actor
#            PlaneActor = vtk.vtkActor()
#            PlaneActor.SetMapper(PlaneMapper)
#            PlaneActor.GetProperty().SetOpacity(0.1)
#            self.ren1.AddActor(PlaneActor)

  def load_conf(self):
    """
    """

    with open ("../in/vtkam/vtkam.segconf", "r") as myseg:
      firstline=myseg.readline()
      segdata=myseg.readlines()

    nblines=len(segdata) ## number of line of interrest in the seg file

    self.nsegm    =   nblines
    self.veclin   =   np.zeros((nblines))
    self.O1       =   np.zeros((nblines))
    self.O2       =   np.zeros((nblines))
    self.O3       =   np.zeros((nblines))
    self.norme    =   np.zeros((nblines))
    self.junctt   =   np.zeros((nblines))
    self.surface  =   np.zeros((nblines))

    for i in range(0,nblines,1):
      tmp=segdata[i].split() ##Create a list of strings for the considered line

      self.veclin[i]   =   int(tmp[1])
      self.O1[i]       =   int(tmp[2])
      self.O2[i]       =   int(tmp[3])
      self.O3[i]       =   int(tmp[4])
      self.norme[i]    =   int(tmp[5])

      self.surface[i]  =   int(tmp[21])

      if str(tmp[14]) is 'T':
        self.junctt[i]   =   1
      elif str(tmp[15]) is 'T':
        self.junctt[i]   =   -1

  def load_filmheader(self):
    """
    Read the header of film.bin
    """
    cristallo, self.avalue , self.bveclin, self.bvecnor = fsrc.film.readfilmheader()

    if (cristallo == 1): fichbase = "CS"
    if (cristallo == 2): fichbase = "BCC"
    if (cristallo == 3): fichbase = "CFC"
    if (cristallo == 4): fichbase = "HCP"
    if (cristallo == 5): fichbase = "ORT"
    if (cristallo == 6): fichbase = "MGO"
    if (cristallo == 7): fichbase = "DC"

    print ""
    if (cristallo > 7 or cristallo < 1):
      print " Cristallo : Unknown value"          # No way
    else :
      print " The crystal symmetry is of type",fichbase

    print " Avalue ",self.avalue

  def load_filmstep(self,kk):
    """
    Read the a particular step of film.bin
    """

    self.nsegm = fsrc.film.readfilmstep(kk)

    slipsys  = self.myConfDic['SLIPSYS']
    curvact      = self.myConfDic['CURV_ACT']
    ldis_act     = self.myConfDic['CURV_LDIS'] #Play the role of scaling factor in curvature vector calculation

    self.kneepoints, self.looppoints, self.slipsystem, self.loop,self.ntotpoints,self.ntotlines,self.nptsctr,self.seg_center,self.curv_vect = fsrc.film.readfilmstepdata(self.mode, self.bveclin,self.bvecnor, self.nsegm,self.avalue,ldis_act, slipsys,curvact)

  def generate_vtkdata(self):
    """
    """

    ## Setup the colors array
    self.colors = vtk.vtkUnsignedCharArray()
    self.colors.SetNumberOfComponents(3)
    self.colors.SetName("Colors")

    #Points for segment coordinates
    self.points = vtk.vtkPoints()

    ## Create an array for line storage
    ## vtkCellArray is a supporting object that explicitly represents cell connectivity.
    ## The cell array structure is a raw integer list of the form:
    ## (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
    ## the cell, and id is a zero-offset index into an associated point list.
    self.lines = vtk.vtkCellArray()

    ## vtkPolyData is a data object that is a concrete implementation of vtkDataSet.
    ## vtkPolyData represents a geometric structure consisting of vertices, lines,
    ## polygons, and/or triangle strips
    self.polyseg = vtk.vtkPolyData()

    ## Actor
    # The actor orchestrates rendering of
    # the mapper's graphics primitives. An actor also refers to properties via a
    # vtkProperty instance, and includes an internal transformation matrix.
    self.polysegActor = vtk.vtkActor()

    ## Create a mapper for the segments
    ## vtkPolyDataMapper is a class that maps polygonal data (i.e., vtkPolyData)
    ## to graphics primitives
    self.polysegMapper = vtk.vtkPolyDataMapper()

    ## Kneecaps
    if (self.myConfDic['KNEECAPS'] in self.TrueList ):

      self.pointsknee = vtk.vtkPoints()             ## Kneecap coordinate
      self.polyknee = vtk.vtkPolyData()             ## Array of kneecap points
      self.sphere = vtk.vtkSphereSource()           ## Create a sphere for a kneecap
      self.sphere.SetRadius(int(self.myConfDic['KNEECAP_RADIUS'])) ## Set kneecap radius
      self.glyph_knee = vtk.vtkGlyph3D()               ## In order to apply a sphere to every points in Glyph
      self.kneecapActor = vtk.vtkActor()            ## Create kneecap Actor
      self.kneecapMapper = vtk.vtkPolyDataMapper()  ## Create a mapper for the kneecap

    ## Tubes filtering
    if (self.myConfDic['TUBE_FILTER'] in self.TrueList ):
      ## Tube filter for rendering
      self.tubes = vtk.vtkTubeFilter()
      self.tubes.SetRadius(int(self.myConfDic['TUBE_RADIUS']))
      self.tubes.SetNumberOfSides(int(self.myConfDic['TUBE_NSIDES']))

    if (self.myConfDic['CURV_ACT'] in self.TrueList ):
      self.pointscurv = vtk.vtkPoints()             ## Center coordinates
      self.polycurv = vtk.vtkPolyData()             ## Array of center points

      ##  Vectors for curvatureay
      self.vectcurv = vtk.vtkFloatArray()
      self.vectcurv.SetNumberOfComponents(3)
      self.vectcurv.SetName("Vectors")

      self.sphere = vtk.vtkSphereSource()           ## Create a sphere for a kneecap
      self.sphere.SetRadius(int(self.myConfDic['KNEECAP_RADIUS'])) ## Set kneecap radius

      # Source for the glyph filter
      self.arrow_curv = vtk.vtkArrowSource()
      self.arrow_curv.SetTipResolution(10)
      self.arrow_curv.SetTipLength(0.1)
      self.arrow_curv.SetTipRadius(0.05)
      self.arrow_curv.SetShaftRadius(0.01)


      self.glyph_curv = vtk.vtkGlyph3D()

      self.curvActor = vtk.vtkActor()            ## Create kneecap Actor
      self.curvMapper = vtk.vtkPolyDataMapper()  ## Create a mapper for the kneecap

  def fill_vtkdata(self):
    """
    """
    self.points.SetNumberOfPoints(self.ntotpoints)
    self.lines.SetNumberOfCells(0)


    for p in range(self.ntotpoints):
      self.points.SetPoint(p, self.looppoints[:,p])
      self.colors.InsertNextTupleValue(self.col[str(self.slipsystem[p]-1)])

    for line in range(self.ntotlines):
    #for line in range(1):
      nbpoints = self.loop[line,0]
      self.lines.InsertNextCell(nbpoints)

      for p in range(1,nbpoints+1) :
        self.lines.InsertCellPoint(self.loop[line,p])

    ## Create PolyData object and add points, segments and color
    self.polyseg.SetPoints(self.points)
    self.polyseg.SetLines(self.lines)
    self.polyseg.GetPointData().SetScalars(self.colors)

    if (self.myConfDic['KNEECAPS'] in self.TrueList ):

      ntotknee = np.int(self.kneepoints[0,0])

      for iknee in range(1,ntotknee):
        kpt = self.kneepoints[:,iknee]
        if  ( self.myConfDic['CLIPPING'] in self.TrueList and
              (kpt[0] < self.origin[0] - self.dx or kpt[0] > self.origin[0] + self.dx) or \
              (kpt[1] < self.origin[1] - self.dy or kpt[1] > self.origin[1] + self.dy) or \
              (kpt[2] < self.origin[2] - self.dz or kpt[2] > self.origin[2] + self.dz) ):
          pass
        else :
          self.pointsknee.InsertNextPoint(kpt)
          self.polyknee.SetPoints(self.pointsknee)

      self.glyph_knee.SetInputData(self.polyknee)
      self.glyph_knee.SetSourceConnection(self.sphere.GetOutputPort())
      self.glyph_knee.Update()

      self.kneecapMapper.SetInputConnection(self.glyph_knee.GetOutputPort())
      self.kneecapMapper.Update()
      self.kneecapActor.SetMapper(self.kneecapMapper)
      #self.kneecapActor.GetProperty().SetInterpolationToFlat()

    if (self.myConfDic['CURV_ACT'] in self.TrueList ):

      for ictr in range(1,self.nptsctr):
        cpt = self.seg_center[:,ictr]
        if  ( self.myConfDic['CLIPPING'] in self.TrueList and
              (cpt[0] < self.origin[0] - self.dx or cpt[0] > self.origin[0] + self.dx) or \
              (cpt[1] < self.origin[1] - self.dy or cpt[1] > self.origin[1] + self.dy) or \
              (cpt[2] < self.origin[2] - self.dz or cpt[2] > self.origin[2] + self.dz) ):
          pass
        else :
          self.pointscurv.InsertNextPoint(cpt)
          self.vectcurv.InsertNextTupleValue(self.curv_vect[:,ictr])

      self.polycurv.SetPoints(self.pointscurv)

      #self.polycurv.GetPointData().SetScalars(self.vectnorm)
      self.polycurv.GetPointData().SetVectors(self.vectcurv)

      self.glyph_curv.SetInputData(self.polycurv)
      self.glyph_curv.SetSourceConnection(self.arrow_curv.GetOutputPort())
      self.glyph_curv.SetScaleFactor(np.float(self.myConfDic['CURV_SCAL']))
      self.glyph_curv.SetScaleModeToScaleByVector()
      #self.glyph_curv.SetVectorModeToUseNormal()
      #self.glyph_curv.SetColorModeToColorByVector()
      #self.glyph_curv.SetScaleModeToScaleByScalar()
      #self.glyph_curv.OrientOn()
      self.glyph_curv.Update()

      self.curvMapper.SetInputConnection(self.glyph_curv.GetOutputPort())
      self.curvMapper.Update()
      self.curvActor.SetMapper(self.curvMapper)


  def add_replicas(self):
    """
    """

    ntrans= (2*self.nrep[0]+1)*(2*self.nrep[1]+1)*(2*self.nrep[2]+1)-1

    repmat=-666* np.ones((ntrans,3))

    n=-1
    for i in range(-1*self.nrep[0],self.nrep[0]+1):
      for j in range(-1*self.nrep[1],self.nrep[1]+1):
        for k in range(-1*self.nrep[2],self.nrep[2]+1):
          if not (i == 0 and j ==0 and k ==0 ) :
                  n=n+1
                  repmat[n]=[i,j,k]

    actorCollection = self.ren1.GetActors()
    actorCollection.InitTraversal()
    numactors = actorCollection.GetNumberOfItems()

### Method 1
    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputData(self.polyseg)

    for item in range(numactors):
      actortmp = actorCollection.GetNextActor()
      if (item>7) :
        for nr in range(ntrans):
          newactor = vtk.vtkActor()
          newactor.ShallowCopy(actortmp)

          pos = newactor.GetPosition()
          newpos = [pos[i] + repmat[nr][i] * np.int(self.boxSize[i]) for i in range(3)]

          actmap= newactor.GetMapper()
          seg = actmap.GetInput()

          translate = vtk.vtkTransform()
          translate.Translate(newpos)

          pgt=vtk.vtkTransformPolyDataFilter()
          pgt.SetTransform(translate)
          pgt.SetInputData(seg)
          pgt.Update()

          seg = pgt.GetOutput()

          appendFilter.AddInputData(seg)

          appendFilter.Update()

    self.polyseg.Initialize()
    self.polyseg = appendFilter.GetOutput()
    self.polysegMapper.SetInputData(self.polyseg)
    self.polysegMapper.Update()

### Old Method
#    actorCollection = self.ren1.GetActors()
#    actorCollection.InitTraversal()
#    numactors = actorCollection.GetNumberOfItems()
#
#    self.assembly = vtk.vtkPropAssembly()
#
#    for item in range(numactors):
#      actortmp = actorCollection.GetNextActor()
#      if (item>7) :
#        for nr in range(ntrans):
#          newactor = vtk.vtkActor()
#          newactor.ShallowCopy(actortmp)
#          pos = newactor.GetPosition()
#          newpos = [pos[i] + repmat[nr][i] * np.int(self.boxSize[i]) for i in range(3)]
#          newactor.SetPosition(newpos)
#          self.assembly.AddPart(newactor)
#    self.ren1.AddActor(self.assembly)


  def render(self):

    #######################  RENDER  ########################################################

    if (self.myConfDic['TUBE_FILTER'] in self.TrueList ):
      ## Tube filter for rendering
      self.tubes.SetInputData(self.polyseg)
      self.tubes.CappingOn()
      self.polysegMapper.SetInputConnection(self.tubes.GetOutputPort())
      self.polysegMapper.Update()

    else :
      ## Add segments to mapper
      self.polysegMapper.SetInputData(self.polyseg)
      self.polysegMapper.Update()

    self.polysegActor.SetMapper(self.polysegMapper)
    ## Add actor to renderer
    self.ren1.AddActor(self.polysegActor)

    if (self.myConfDic['REPLICA'] in self.TrueList ):
      self.add_replicas()

    # Assign actors to Renderer. A renderer is like a
    # viewport. It is part or all of a window on the screen and it is
    # responsible for drawing the actors it has.  We also set the
    # background color here.


    if (self.myConfDic['KNEECAPS'] in self.TrueList ):
      self.ren1.AddActor(self.kneecapActor)

    if (self.myConfDic['CURV_ACT'] in self.TrueList ):
      self.ren1.AddActor(self.curvActor)

    if (self.surfaceActor is not None ):
      self.ren1.AddActor(self.surfaceActor)

    if ( self.myConfDic['CLIPPING'] in self.TrueList ) :
      self.polyseg = self.clipper(self.polyseg,self.dx,self.dy,self.dz)
      self.polysegMapper.SetInputData(self.polyseg)
      self.polysegMapper.Update()

#    self.ren1.AddActor(self.assembly)
    self.renWin.Render()


  def initialize_vtkobjects(self):
    """
    """
    # Segments
    self.points.Initialize()
    self.lines.Initialize()
    self.polyseg.Initialize()
    self.colors.Initialize()
    self.ren1.RemoveActor(self.polysegActor)

    # Kneecaps
    if (self.myConfDic['KNEECAPS'] in self.TrueList ):
      self.pointsknee.Initialize()
      self.polyknee.Initialize()
      self.ren1.RemoveActor(self.kneecapActor)

    # Cruvature vectors
    if (self.myConfDic['CURV_ACT'] in self.TrueList ):
      self.pointscurv.Initialize()
      self.vectcurv.Initialize()
      self.polycurv.Initialize()
      self.ren1.RemoveActor(self.curvActor)

  def clipper(self, src, dx, dy, dz):
      '''
      Clip a vtkPolyData source.
      A cube is made whose size corresponds the the bounds of the source.
      Then each side is shrunk by the appropriate dx, dy or dz. After
      this operation the source is clipped by the cube.
      :param: src - the vtkPolyData source
      :param: dx - the amount to clip in the x-direction
      :param: dy - the amount to clip in the y-direction
      :param: dz - the amount to clip in the z-direction
      :return: vtkPolyData.
      '''

      clipFunction = vtk.vtkImplicitBoolean()
      clipFunction.SetOperationTypeToUnion()

      if (self.normal[0] != None) :

        norm = np.linalg.norm(self.normal)
        unitnorm = [self.normal[i]/norm for i in range(3)]
        oppunitnorm = [-unitnorm[i] for i in range(3)]

        Origin1 = [self.origin[0] + unitnorm[0]*dx/2., self.origin[1] + unitnorm[1]*dx/2., self.origin[2] + unitnorm[2]*dx/2.]
        Origin2 = [self.origin[0] - unitnorm[0]*dx/2., self.origin[1] - unitnorm[1]*dx/2., self.origin[2] - unitnorm[2]*dx/2.]

        plane1 = vtk.vtkPlane()
        plane1.SetOrigin(Origin1)
        plane1.SetNormal(oppunitnorm)

        plane2 = vtk.vtkPlane()
        plane2.SetOrigin(Origin2)
        plane2.SetNormal(unitnorm)

        clipFunction.AddFunction(plane1)
        clipFunction.AddFunction(plane2)

      else :

        plane1 = vtk.vtkPlane()
        plane1.SetOrigin(self.origin[0] + dx, 0, 0)
        plane1.SetNormal(-1, 0, 0)


        plane2 = vtk.vtkPlane()
        plane2.SetOrigin(self.origin[0] - dx, 0, 0)
        plane2.SetNormal(1, 0, 0)

        plane3 = vtk.vtkPlane()
        plane3.SetOrigin(0,self.origin[1] + dy, 0)
        plane3.SetNormal(0, -1, 0)

        plane4 = vtk.vtkPlane()
        plane4.SetOrigin(0, self.origin[1] - dy, 0)
        plane4.SetNormal(0, 1, 0)

        plane5 = vtk.vtkPlane()
        plane5.SetOrigin(0, 0, self.origin[2] + dz)
        plane5.SetNormal(0, 0, -1)

        plane6 = vtk.vtkPlane()
        plane6.SetOrigin(0, 0, self.origin[2] - dz)
        plane6.SetNormal(0, 0, 1)

        clipFunction.AddFunction(plane1)
        clipFunction.AddFunction(plane2)
        clipFunction.AddFunction(plane3)
        clipFunction.AddFunction(plane4)
        clipFunction.AddFunction(plane5)
        clipFunction.AddFunction(plane6)

      # Clip it.
      clipper =vtk.vtkClipPolyData()
      clipper.SetClipFunction(clipFunction)

      try :
        clipper.SetInputData(src)
      except :
        clipper.SetInputConnection(src)

      clipper.GenerateClipScalarsOff() # 'On' if we want to see what was cutted
      clipper.GenerateClippedOutputOff()
      #clipper.GenerateClippedOutputOn()
      clipper.Update()
      return clipper.GetOutput()

try :
  vtkam=Vtkam()
except KeyboardInterrupt: quit()