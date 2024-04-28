import numpy as np
import matplotlib.pyplot as plt
from subprocess import call


def f_Gross(faW):
  return 1.122-0.231*faW+10.550*faW**2-21.710*faW**3+30.382*faW**4

def f_Tada1(faW):
  a1 = np.sqrt (2./(np.pi*faW) *np.tan(np.pi/2.*faW))
  a2 = (0.752+2.02*faW+0.37*(1-np.sin(np.pi/2.*faW))**3)/np.cos(np.pi/2.*faW)
  return a1*a2

def f_Tada2(faW):
  a1 = 0.265*(1-faW)**4
  a2 = (0.857+0.265*faW)/(1-faW)**1.5
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

ans=True
while ans:
    print("""
    1.Plot geometric functions
    2.Calculate SIF
    x.Exit/Quit
    """)
    ans=raw_input("What would you like to do? ")
    if ans=="1":

        faW = np.arange(0.01, .99, 0.01)
        faWreduced= np.arange(0.01, .6, 0.01)

        # Plots
        plt.plot(faW, (1-faW)**1.5*f_Tada1(faW),label='Tada_1')
        plt.plot(faW, (1-faW)**1.5*f_Tada2(faW),label='Tada_2')
        plt.plot(faWreduced, (1-faWreduced)**1.5*f_Gross(faWreduced),'r',label = 'Gross-Brown')
        plt.plot(faW, (1-faW)**1.5*f_Gross(faW),'r--')

        toto=[]
        for i in range(np.size(faW)):
          toto = np.append(toto, f_TheoFett(faW[i],0.46))

        plt.plot(faW, (1-faW)**1.5*toto,'+k',label='Theo_Fett(H/W=0.46)')

	plt.xlabel('a/W')
	plt.title('(1-a/W)**3/2*F_r)')
	plt.ylabel('(1-a/W)**3/2*F_r)')
	plt.grid(True)
	plt.ylim((0.96,1.2))
	plt.legend(loc='upper left')
	plt.show()

        print("\n")
    elif ans=="2":
      call(["open", "Single_edge_notch.png"])

      sig = np.float(raw_input("> Sigma = "))
      a = np.float(raw_input("> a = "))
      W = np.float(raw_input("> W = "))
      H = np.float(raw_input("> H = "))
      alpha = a/W
      HdW = H/W

      print " -- SIF -- "
      print "Tada1    :", sig*np.sqrt(np.pi*a)*f_Tada1(alpha)
      print "Tada2    :", sig*np.sqrt(np.pi*a)*f_Tada2(alpha)
      print "Gross    :", sig*np.sqrt(np.pi*a)*f_Gross(alpha)
      print "TheoFett :", sig*np.sqrt(np.pi*a)*f_TheoFett(alpha, HdW)

      print("\n")
    elif ans=="x":
      print("\n Goodbye")
      ans = None
    else:
       print("\n Not Valid Choice Try again")

