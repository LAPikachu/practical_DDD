import numpy as np
import sys

def calculate_faw(a,W,H,type):

  # Sigma : tensile loading constant stress
  # a     : crack length
  # W     : plate weigth
  # H     : plate half height


  def f_Gross(alpha):
    """
    The Stress Analysis of Cracks Handbook
    Hiroshi Tada, Paul Croce Paris, George Rankin Irwin
    """
    return 1.122-0.231*alpha+10.550*alpha**2-21.710*alpha**3+30.382*alpha**4

  def f_Tada1(alpha):
    """
    The Stress Analysis of Cracks Handbook
    Hiroshi Tada, Paul Croce Paris, George Rankin Irwin
    """
    a1 = np.sqrt (2./(np.pi*alpha) *np.tan(np.pi/2.*alpha))
    a2 = (0.752+2.02*alpha+0.37*(1-np.sin(np.pi/2.*alpha))**3)/np.cos(np.pi/2.*alpha)
    return a1*a2

  def f_Tada2(alpha):
    """
    The Stress Analysis of Cracks Handbook
    Hiroshi Tada, Paul Croce Paris, George Rankin Irwin
    """
    a1 = 0.265*(1-alpha)**4
    a2 = (0.857+0.265*alpha)/(1-alpha)**1.5
    return a1+a2

  def f_TheoFett(atest, HdWtest):

    HdW = [1.5, 1.25, 1.00, 0.75, 0.5, 0.4, 0.3, 0.25]
    alpha = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0]

    index_row = 1
    index_col = 1

    while True:
      if (alpha[index_row]>atest) :
        index_row = index_row - 1
        break
      if (index_row >= len(alpha)-1) :
        index_row = len(alpha)-1
        break
      index_row=index_row+1

    while True:
      if (HdW[index_col] <= HdWtest) :
        break
      if (index_col >= len(HdW)-1) :
        index_col = len(HdW)-1
        break
      index_col=index_col+1


    #index_col = min(range(len(HdW)), key=lambda i: abs(HdW[i]-HdWtest))

    if (atest > 1. or atest < 0.) :
      print "!> Alpha should be between 0 and 1"

    data =[
    [1.1215, 1.1215, 1.1215, 1.1215, 1.1215, 1.1215, 1.1215, 1.1215],
    [1.0170, 1.0172, 1.0174, 1.0182, 1.0352, 1.0649, 1.1455, 1.2431],
    [0.9800, 0.9799, 0.9798, 0.9877, 1.0649, 1.1625, 1.3619, 1.5358],
    [0.9722, 0.9723, 0.9729, 0.9840, 1.0821, 1.2134, 1.4892, 1.7225],
    [0.9813, 0.9813, 0.9819, 0.9915, 1.0819, 1.2106, 1.5061, 1.7819],
    [0.9985, 0.9986, 0.9989, 1.0055, 1.0649, 1.1667, 1.4298, 1.7013],
    [1.0203, 1.0203, 1.0204, 1.0221, 1.0496, 1.1073, 1.2898, 1.5061],
    [1.0440, 1.0441, 1.0441, 1.0442, 1.0522, 1.0691, 1.1498, 1.2685],
    [1.0683, 1.0683, 1.0683, 1.0690, 1.0691, 1.0734, 1.0861, 1.1201],
    [1.1215, 1.1215, 1.1215, 1.1215, 1.1215, 1.1215, 1.1215, 1.1215]
    ]

    x1  = alpha[index_row]
    x2  = alpha[index_row+1]
    y1  = HdW[index_col]
    y2  = HdW[index_col-1]

    Q11 = data[index_row][index_col]
    Q21 = data[index_row+1][index_col]
    Q12 = data[index_row][index_col-1]
    Q22 = data[index_row+1][index_col-1]

    a11=(x2-atest)*(y2-HdWtest)/((x2-x1)*(y2-y1))
    a21=(atest-x1)*(y2-HdWtest)/((x2-x1)*(y2-y1))
    a12=(x2-atest)*(HdWtest-y1)/((x2-x1)*(y2-y1))
    a22=(atest-x1)*(HdWtest-y1)/((x2-x1)*(y2-y1))

    resul =  (a11*Q11 + a12*Q12 + a21*Q21 + a22*Q22)/(1-atest)**1.5

    #print "(",atest,HdWtest,")"," - ", "(",x1,x2,y1,y2,")", " - ","(",Q11,Q21,Q12,Q22,")", resul

    return resul


  alpha = a/W

  if (type == 1) :
    #FETT
    print "   > FETT Method"
    HdW   = H/W
    faw = f_TheoFett(alpha,HdW)

  elif (type == 2) :
    print "   > Gross Method"
    #Gross
    faw = f_Gross(alpha)

  elif (type == 3) :
    print "   > Tada1 Method"
    #Tada1
    faw = f_Tada1(alpha)

  elif (type == 4) :
    print "   > Tada2 Method"
    #Tada2
    faw = f_Tada2(alpha)

  #KI = Sigma*np.sqrt(a*np.pi)*faw

  return faw

print " > Calculating geometric factor for SIF"

type  = np.float(sys.argv[1])
a     = np.float(sys.argv[2])
H     = np.float(sys.argv[3])
W     = np.float(sys.argv[4])

faw = calculate_faw(a,W,H,type)
print "   > Gives :",faw
np.savetxt("../out/siffactresult.dat",[faw])