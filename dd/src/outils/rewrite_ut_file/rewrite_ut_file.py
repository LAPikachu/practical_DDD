import numpy as np
from decimal import *
from sys import argv
import shutil
import os

#Get arguments, scriptname and filename
script, filename = argv

#Backup the .ut file
shutil.copyfile(filename, filename+'.backup')

#Store the whole file in order to get the header
header=open(filename).readlines()

# Store the lines after the header (start at the fourth line)
time=np.loadtxt(filename,skiprows=4)

#Initial Step time
initialstep=time[0,4]

#Remove old ut file
os.remove(filename)

#Open new ut file
f=open(filename,"w")

#Write header
f.writelines( header[0:4] )

#Write time lines with appropriate formatting
np.savetxt(f,time,fmt='%d %d %d %d %.5e')

#close ut file
f.close()
