import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

def f_Tada1(faW):
  a1 = 1-0.5*faW+0.370*faW**2-0.044*faW**3
  a2 = np.sqrt(1-faW)
  return a1/a2

def f_Tada2(faW):
  a1 = 1-0.025*faW**2+ 0.06*faW**4
  a2 = np.sqrt(1./np.cos(np.pi/2*faW))
  return a1*a2

def f_TheoFett(atest, HdWtest):

  HdW = [1.5, 1.25, 1.00, 0.75, 0.5, 0.35]
  alpha = [0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0]

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
	[1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
	[0.916, 0.924, 0.940, 0.977, 1.051, 1.182],
	[0.888, 0.905, 0.940, 1.008, 1.147, 1.373],
	[0.869, 0.890, 0.942, 1.053, 1.262, 1.562],
	[0.851, 0.877, 0.943, 1.099, 1.391, 1.742],
	[0.827, 0.856, 0.937, 1.130, 1.533, 1.938],
	[0.816, 0.826, 0.914, 1.125, 1.668, 2.197],
	[0.814, 0.818, 0.840, 1.088, 1.689, 2.41],
	[0.826, 0.826, 0.826, 0.826, 0.826, 0.826]
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

  resul =  (a11*Q11 + a12*Q12 + a21*Q21 + a22*Q22)

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
        plt.plot(faW, (1-faW)**0.5*f_Tada1(faW),label='Tada_1')
        plt.plot(faW, (1-faW)**0.5*f_Tada2(faW),label='Tada_2')

        toto=[]
        for i in range(np.size(faW)):
          toto = np.append(toto, f_TheoFett(faW[i],0.46))

        plt.plot(faW, (1-faW)**1.5*toto,'+k',label='Theo_Fett(H/W=0.46)')

	plt.xlabel('a/W')
	plt.title('(1-a/W)**1/2*F_r)')
	plt.ylabel('(1-a/W)**1/2*F_r)')
	plt.grid(True)
#	plt.ylim((0.96,1.2))
	plt.legend(loc='lower left')
	plt.show()

        print("\n")
    elif ans=="2":
      call(["open", "Rectangular_plate_with_an_internal_crack.png"])

      sig = np.float(raw_input("> Sigma = "))
      a = np.float(raw_input("> a = "))
      W = np.float(raw_input("> W = "))
      H = np.float(raw_input("> H = "))
      alpha = a/W
      HdW = H/W

      print " -- SIF -- "
      print "Tada1    :", sig*np.sqrt(np.pi*a)*f_Tada1(alpha)*np.sqrt(1-alpha)
      print "Tada2    :", sig*np.sqrt(np.pi*a)*f_Tada2(alpha)*np.sqrt(1-alpha)
      print "TheoFett :", sig*np.sqrt(np.pi*a)*f_TheoFett(alpha, HdW)*np.sqrt(1-alpha)

      print("\n")
    elif ans=="x":
      print("\n Goodbye")
      ans = None
    else:
       print("\n Not Valid Choice Try again")

