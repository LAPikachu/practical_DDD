import numpy as np

#xloc=148016./2.
sx=np.array([0. , 778752.])
yloc=389376./2.
#sy=np.array([0. , 481980.])
#zloc=148016./2.
sz=np.array([0. , 121680.])

npointsx=100
#npointsy=100
npointsz=50

#stations=np.array([[-666. for i in range(3)] for j in range(npointsx*npointsy)])
#stations=np.array([[-666. for i in range(3)] for j in range(npointsy*npointsz)])
stations=np.array([[-666. for i in range(3)] for j in range(npointsx*npointsz)])

stepx= (sx[1]-sx[0])/float((npointsx+1))
#stepy= (sy[1]-sy[0])/float((npointsy+1))
stepz= (sz[1]-sz[0])/float((npointsz+1))

row=0

for j in range(npointsx):
#for j in range(npointsy):
    xloc=sx[0]+(float(j)+1)*stepx
    #yloc=sy[0]+(float(j)+1)*stepy

#    for k in range(npointsy):
    for k in range(npointsz):

         #yloc=sy[0]+(float(k)+1)*stepy
         zloc=sz[0]+(float(k)+1)*stepz

         stations[row] = np.array([float(xloc),float(yloc),float(zloc)])
         row=row+1

np.savetxt('locations_dat', stations, fmt='%20.3f', delimiter=' ', newline='\n')
