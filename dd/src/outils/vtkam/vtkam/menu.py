import numpy as np
import vtk
import utilities as util

def define_options(vtkam):
  vtkam.options = {
#    "a": "...",
    "b": "Boxes",
    "c": "Camera",
    "d": "Colors",
    "e": "Curvature vectors",
    "f": "Load VTK unstructured data files",
    "g": "Slip system selection",
    "h": "Clipping",
#    "i": "...",
#    "j": "...",
    "k": "Kneecaps",
#    "l": "...",
    "m": "Show options (menu)",
#    "n": "...",
#    "o": "...",
#    "p": "...",
#    "q": "...",
    "r": "Reading",
    "s": "Screenshot",
    "t": "Tubes",
#    "u": "...",
    "v": "Objects Visibility",
#    "w": "...",
     "x": "Replica",
#    "y": "...",
    "z": "Exit",
  }

def apply_option(vtkam):
  """
  """
  input = vtkam.userinput

  def print_menu(mymenu):
    for key in sorted(mymenu) :
      print "  ", key,": ",mymenu[key]

  Loop = True

  ###########
  ##>  A  <##
  ###########
  if (input == 'a'):
    pass
  ###########
  ##>  B  <##
  ###########
  elif (input == 'b'):

    suboptions = {
      "1": "Add a Box ",
      "2": "Remove last Box",
      "r": "Return",
    }
    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :
      vtkam.set_extrabox()

    elif (subinput == 2 ) :
      vtkam.del_extrabox()

  ###########
  ##>  C  <##
  ###########
  elif (input == 'c'):

    suboptions = {
      "1": "View direction",
      "2": "Focal point",
      "3": "View up vector",
      "4": "Reset view",
      "r": "Return",
    }

    print_menu(suboptions)

    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      try :
        pos = util.get_int3(" !> Values (x3) ? ")

        norm = np.linalg.norm(pos)
        pos = [pos[i]/norm for i in range(3)]

        Focal = vtkam.camera.GetFocalPoint()
        dist = int(vtkam.camera.GetDistance())

        vtkam.camera.SetPosition(int(dist*pos[0])+Focal[0],int(dist*pos[1])+Focal[1],int(dist*pos[2])+Focal[2])
        vtkam.camera.SetFocalPoint(Focal[0],Focal[1],Focal[2])

        vtkam.renWin.Render()

        print ""
        print "   Focal point          : ", vtkam.camera.GetFocalPoint()
        print "   Camera position      : ", vtkam.camera.GetPosition()
        print "   Camera distance      : ", vtkam.camera.GetDistance()
        print "   Projection direction : ", vtkam.camera.GetDirectionOfProjection()
        print "   View up vector       : ", vtkam.camera.GetViewUp()

      except IndexError :
        print " !/ Try again"


    elif (subinput == 2 ) :

      try :
        pos = util.get_int3(" !> Values (x3) ? ")

        vtkam.camera.SetFocalPoint(pos[0],pos[1],pos[2])

        vtkam.renWin.Render()

        print ""
        print "   Focal point          : ", vtkam.camera.GetFocalPoint()
        print "   Camera position      : ", vtkam.camera.GetPosition()
        print "   Camera distance      : ", vtkam.camera.GetDistance()
        print "   Projection direction : ", vtkam.camera.GetDirectionOfProjection()
        print "   View up vector       : ", vtkam.camera.GetViewUp()

      except IndexError :
        print " !/ Try again"

    elif (subinput == 3 ) :

      try :
        pos = util.get_int3(" !> Values (x3) ? ")

        vtkam.camera.SetViewUp(pos[0],pos[1],pos[2])

        vtkam.renWin.Render()

        print ""
        print "   Focal point          : ", vtkam.camera.GetFocalPoint()
        print "   Camera position      : ", vtkam.camera.GetPosition()
        print "   Camera distance      : ", vtkam.camera.GetDistance()
        print "   Projection direction : ", vtkam.camera.GetDirectionOfProjection()
        print "   View up vector       : ", vtkam.camera.GetViewUp()

      except IndexError :
        print " !/ Try again"

    elif (subinput == 4 ) :
      vtkam.ren1.ResetCamera()
      vtkam.camera.SetFocalPoint(0.5*int(vtkam.boxSize[0]),0.5*int(vtkam.boxSize[1]),0.5*int(vtkam.boxSize[2]))
      vtkam.camera.SetPosition(vtkam.dist,0.5*int(vtkam.boxSize[1]),0.5*int(vtkam.boxSize[2]))
      vtkam.camera.SetViewUp(0,0,1)

  ###########
  ##>  D  <##
  ###########
  elif (input == 'd'):
    suboptions = {
      "1": "Show segments",
      "2": "Show Junctions and GD",
      "r": "Return",
    }

    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      vtkam.myConfDic['COLORS']= 0
      vtkam.set_colors()

    elif (subinput == 2 ) :

      vtkam.myConfDic['COLORS']= 1
      vtkam.set_colors()

  ###########
  ##>  E  <##
  ###########
  elif (input == 'e'):
    suboptions = {
      "1": "Show Curvature vector [On/Off] - Only in reading mode 2",
      "2": "Discretization length          - Only in reading mode 2",
      "3": "Scaling factor                 - Only in reading mode 2",
      "r": "Return",
    }

    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      if (vtkam.myConfDic['CURV_ACT'] in vtkam.TrueList ):
        vtkam.myConfDic['CURV_ACT'] = 'F'
      else :
        vtkam.myConfDic['CURV_ACT'] = 'T'
      print 'Curvature vector activation is set to ', vtkam.myConfDic['CURV_ACT']

    elif (subinput == 2 ) :

      try :
        vtkam.myConfDic['CURV_LDIS'] = util.get_non_negative_float(" !> Discretization length ? ")
        print " > Change Discretization length value will be effective for the next steps"
      except ValueError:
        print "Not an real"

    elif (subinput == 3 ) :

      try :
        vtkam.myConfDic['CURV_SCAL'] = util.get_non_negative_float(" !> Scaling factor ? ")
      except ValueError:
        print "Not an real"

  ###########
  ##>  F  <##
  ###########
  elif (input == 'f'):

    suboptions = {
      "1": "Choose a file",
      "2": "Opacity",
      "r": "Return",
    }

    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      myfilepath = raw_input(' !> Path to vtk unstructured data file [Default : ../out/simulationvolume.vtk]: ') #example : '../out/simulationvolume.vtk'
      if (myfilepath=='') :
          myfilepath = '../out/simulationvolume.vtk'
      vtkam.surfaceActor = vtkam.ReadVTKSurfaceFile(myfilepath)

    elif (subinput == 2 ) :
      if (vtkam.surfaceActor is not None):
        opacity = util.get_non_negative_float(" !> Opacity [0-1]: ")
        vtkam.surfaceActor.GetProperty().SetOpacity(opacity)
      else :
        print ' !> Please load a file first '


  ###########
  ##>  G  <##
  ###########
  elif (input == 'g'):

      try :
        vtkam.myConfDic['SLIPSYS'] = util.get_non_negative_int(" !> Which slip system (> 12 = ALL) ? - Only in reading mode 2 : ")
        print " > Change slip system value will be effective for the next steps"
      except ValueError:
        print "Not an integer"

  ###########
  ##>  H  <##
  ###########
  elif (input == 'h'):
    suboptions = {
      "1": "ON/OFF",
      "2": "Cube",
      "3": "Slice",
      "r": "Return",
    }
    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      if (vtkam.myConfDic['CLIPPING'] in vtkam.TrueList ):
        vtkam.myConfDic['CLIPPING']= 'F'
      else :
        vtkam.myConfDic['CLIPPING'] = 'T'
      print 'Clipping is set to ', vtkam.myConfDic['CLIPPING']

    if (subinput == 2 ) :

        print "Main box center ", np.float(vtkam.boxSize[0])/2.,np.float(vtkam.boxSize[1])/2.,np.float(vtkam.boxSize[2])/2.
        vtkam.origin = util.get_float3("Origin (X,Y,Z) !> ")
        vtkam.dx = util.get_non_negative_float("dx !> ")
        vtkam.dy = util.get_non_negative_float("dy !> ")
        vtkam.dz = util.get_non_negative_float("dz !> ")
        vtkam.normal = [None,None,None]

        vtkam.clipping = True

    elif (subinput == 3 ) :

        vtkam.origin = util.get_float3("Origin (X,Y,Z) !> ")
        vtkam.normal = util.get_int3("Normal (a,b,c) !> ")
        vtkam.dx = util.get_non_negative_float("Distance between planes !> ")

        vtkam.clipping = True

  ###########
  ##>  I  <##
  ###########
#  if (input == 'i'):

  ###########
  ##>  J  <##
  ###########
#  if (input == 'j'):

  ###########
  ##>  K  <##
  ###########
  elif (input == 'k'):
    suboptions = {
      "1": "Show Kneecaps [On/Off] - Reading mode 1 only",
      "2": "Define Radius [Integer]",
      "r": "Return",
    }
    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :
      print " > Current reading mode = ",  vtkam.myConfDic['READING_MODE']
      if (vtkam.myConfDic['KNEECAPS'] in vtkam.TrueList ):
        vtkam.myConfDic['KNEECAPS']='F'
      else :
        vtkam.myConfDic['KNEECAPS']='T'
      print ' > Kneecaps are set to ', vtkam.myConfDic['KNEECAPS']

    elif (subinput == 2 ) :

      try :
        vtkam.myConfDic['KNEECAP_RADIUS']= util.get_non_negative_float(" !> Value ? ")
      except ValueError:
        print "Not an integer"

  ###########
  ##>  L  <##
  ###########
#  if (input == 'l'):

  ###########
  ##>  M  <##
  ###########
  elif (input == 'm'):
    print_menu(vtkam.options)

  ###########
  ##>  N  <##
  ###########
#  if (input == 'n'):

  ###########
  ##>  O  <##
  ###########
#  if (input == 'o'):

  ###########
  ##>  P  <##
  ###########
#  if (input == 'p'):

  ###########
  ##>  Q  <##
  ###########
#  if (input == 'q'):

  ###########
  ##>  R  <##
  ###########
  elif (input == 'r'):
    suboptions = {
       "1": "MODE 1 - All segments and kneecaps are vtk object (slower) ",
       "2": "MODE 2 - Reduced set of vtk object - No kneecaps ",
       "3": "Change Snapshot frequency ",
       "4": "Go to step #",
       "r": "Return",
    }

    print_menu(suboptions)

    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :
      vtkam.mode = 1
      print " > Change for mode 1 reading will be effective for the next steps"

    elif (subinput == 2 ) :
      vtkam.mode = 2
      print " > Change for mode 2 reading will be effective for the next steps"

    elif (subinput == 3 ) :
      vtkam.perio = util.get_non_negative_int(" !> Period ? ")

    elif (subinput == 4 ) :
      vtkam.kk = util.get_non_negative_int(" !> Step ? ")

  ###########
  ##>  S  <##
  ###########
  elif (input == 's'):

    # screenshot code:
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(vtkam.renWin)
    #w2if.SetMagnification(3)
    #w2if.SetInputBufferTypeToRGBA()
    w2if.ReadFrontBufferOff()
    w2if.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName("screenshot.png")
    writer.SetInputConnection(w2if.GetOutputPort())
    writer.Write()
    print "  !> An image 'screenshot.png' has been created in your current folder"

  ###########
  ##>  T  <##
  ###########
  elif (input == 't'):

    suboptions = {
      "1": "Tube Filtering [On/Off]",
      "2": "Radius [Integer]",
      "3": "Tube facets [Integer]",
      "r": "Return",
    }

    print_menu(suboptions)

    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :
      if (vtkam.myConfDic['TUBE_FILTER'] in vtkam.TrueList ):
        vtkam.myConfDic['TUBE_FILTER']='F'
      else :
        vtkam.myConfDic['TUBE_FILTER']='T'
      print 'Tube filter is set to ', vtkam.myConfDic['TUBE_FILTER']

    elif (subinput == 2 ) :
      try :
        vtkam.myConfDic['TUBE_RADIUS']= util.get_non_negative_float(" !> Value ? ")
      except ValueError:
        print "Not an integer"

    elif (subinput == 3 ) :
      try :
        vtkam.myConfDic['TUBE_NSIDES']= util.get_non_negative_int(" !> Value ? ")
      except ValueError:
        print "Not an integer"

  ###########
  ##>  U  <##
  ###########
#  if (input == 'u'):

  ###########
  ##>  V  <##
  ###########
  elif (input == 'v'):
    suboptions = {
      "1": "Main box outline [On/Off]",
      "2": "Coordinate system [On/Off] ",
      "r": "Return",
    }
    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      if (vtkam.BoxActor.GetVisibility()):
        vtkam.BoxActor.VisibilityOff()
      else :
        vtkam.BoxActor.VisibilityOn()

    elif (subinput == 2 ) :

      if (vtkam.AxesActor.GetVisibility()):
        vtkam.actor0.VisibilityOff()
        vtkam.AxesActor.VisibilityOff()
      else :
        vtkam.actor0.VisibilityOn()
        vtkam.AxesActor.VisibilityOn()

  ###########
  ##>  W  <##
  ###########
#  if (input == 'w'):

  ###########
  ##>  X  <##
  ###########
  elif (input == 'x'):

    suboptions = {
      "1": "Show Replica [On/Off]",
      "2": "Number of replica (nx,ny,nz) ",
      "r": "Return",
    }
    print_menu(suboptions)
    subinput = util.get_non_negative_int(" !> ")

    if (subinput == 1 ) :

      if (vtkam.myConfDic['REPLICA'] in vtkam.TrueList ):
        vtkam.myConfDic['REPLICA']='F'
      else :
        vtkam.myConfDic['REPLICA']='T'
      print 'REPLICA are set to ', vtkam.myConfDic['REPLICA']

    elif (subinput == 2 ) :

      try :
        vtkam.nrep[0] = util.get_non_negative_int("nx !> ")
        vtkam.nrep[1] = util.get_non_negative_int("ny !> ")
        vtkam.nrep[2] = util.get_non_negative_int("nz !> ")
      except ValueError:
        print "Not an integer"

      if (vtkam.nrep[0] < 0  or vtkam.nrep[1] < 0  or vtkam.nrep[2] < 0 ) :
        print "Number of Replica must be a positive integer "
        vtkam.nrep[0] =0
        vtkam.nrep[1] =0
        vtkam.nrep[2] =0
      elif ((2*vtkam.nrep[0]+1)*(2*vtkam.nrep[1]+1)*(2*vtkam.nrep[2]+1) > 350 ) :
          print "Are you sure to use such high number for replica ? It may slow down the vizualisation a lot... "
          res = util.get_string(" !> ")
          if (res not in vtkam.TrueList ):
            vtkam.nrep[0] =0
            vtkam.nrep[1] =0
            vtkam.nrep[2] =0
          else :
            print "So be it ... "
  ###########
  ##>  Y  <##
  ###########
#  if (input == 'y'):

  ###########
  ##>  Z  <##
  ###########
  # No Z : program exit

  else :
     Loop = False

  return Loop