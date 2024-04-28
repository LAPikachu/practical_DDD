import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from subprocess import call


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
    print ("!> Alpha should be between 0 and 1")
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

faW = np.linspace(0.01, .99, 100)
HdW = np.linspace(0.351, 1.5, 100)

toto=np.ones((100,100))
x,y = np.meshgrid(np.linspace(0.01, .99, 100),np.linspace(0.25, 1.5, 100))

for i in range(np.size(faW)):
  for j in range(np.size(HdW)):
    toto[i][j] = f_TheoFett(faW[i],HdW[j])

#
# Plots
plt.pcolormesh(x,y,toto)
cbar=plt.colorbar()
plt.scatter(0.25, 0.46, s=50,color="w", alpha=1.)
plt.text(0.27, 0.48, '(0.25,0.46) = '+ str(f_TheoFett(0.25,0.46)) ,color="w", fontsize=10)

plt.xlabel('a/W')
plt.title(r'(1-a/W)**1/2*F_r')
plt.ylabel('H/W')
plt.grid(True)
plt.xlim((0.01,0.99))
plt.ylim((0.25,1.25))
plt.legend(loc='upper left')
plt.show()

