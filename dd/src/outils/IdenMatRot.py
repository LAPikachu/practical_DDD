from numpy import *
import sys

# program to identify a transformation (rotation) matrix with imposed boundary conditions on two vector vref and vref2
# for instance in what follows, we want a matrix rotating a direction vref as close as possible to [100],
# and vref2 as close as possible to [010]
# Remember that:
# -the matrix used is mm is going to be ([v1],[v2],[v3])
# -the homothetie factor is the inv factor normalizing the matrix
# -the facteur_boite imposed by the final matrix is usually defined as m=h1*h2*0, if [h1,h2,0] is a vector
# existing the BVD.xxx file generated with the transformation matrix
# facteur_boite can also be determined directly with mM, by activating a special procedure existing in 07init.f90
# Such special procedure must be activated directly in the source file with the variable "calculboite"

# Important initialization:
# The doamin of iterations considered to identify the matrix
if (len(sys.argv)>1): limit = int(sys.argv[1])
else : limit = 30
# The tolerance for the scalar product between vectors vref and [100], as well as vref2 and [010]
if (len(sys.argv)>2): error = float(sys.argv[2])
else : error = 0.98

if (len(sys.argv)>3): print 'Too many Arguments';quit()

# The reference vector we want to approach as [100]:
vref=array([1.,-2.,1.])
nvref=vref/sqrt(dot(vref,vref))
# The reference vector we want to approach as [010]:
vref2=array([-1.,0.,1.])
nvref2=vref2/sqrt(dot(vref2,vref2))

# -----------------------------------------
with open('./IdenMatRot_database'+'/Rotmat_('+str(vref[0])+'_'+str(vref[1])+'_'+str(vref[2])+')_('+str(vref2[0])+'_'+str(vref2[1])+'_'+str(vref2[2])+').dat', 'w') as resulfile :
  resulfile.write('>> IdenMatRot \nResul obtained when looking for a rotation matrice with integers with :' )
  resulfile.write('\nv1=' + str(vref) + ' and v2=' + str(vref2) )
  resulfile.write('\nInspection limit was ' + str(limit) +  ' and tolerance was ' + str(error) )

  x=-limit

  nresul=0

  while (x<limit):

    while (sign(x)!=sign(vref[0])) : x=x+1
    else : print 'Looking with x=',x

    y=-limit
    while (y<limit):

      z=-limit
      while (z<limit):

        # The approximate shape of v1 we want to try
        v1=array([x, y, z])
        norme1=sqrt(dot(v1,v1))
        if (norme1<=0) : break
        v1n=v1/norme1
        # The tolerate error between v1 and nvref
        if dot(v1n,nvref) > error:

          # We wants only vector with an integer norm
          if norme1 == int(norme1):

            # Loops to find a vector v2 compatible with v1
            x2=-limit
            while (x2<limit):

              y2=-limit
              while (y2<limit):

                z2=-limit
                while (z2<limit):

                  #the vector v2
                  v2=array([x2, y2 ,z2])
                  norme2=sqrt(dot(v2,v2))
                  if (norme2<=0) : break
                  v2n=v2/norme2

                  # The tolerate error between v2n and vref2
                  if dot(v2n,nvref2) > error:

                    # v1 and v2 must have the same norm
                    if norme1 == norme2:

                      # v3 is normal to v1 and v2
                      v3=cross(v1,v2)/norme2
                      norme3=sqrt(dot(v3,v3))

                      if (norme3<=0) : break

                      if norme1 == norme3:
                        nresul=nresul+1
                        print 'Found a matrice with homothety factor (matrix norme) : ',norme3,' and quality  ',int((sqrt(dot(v1n,nvref)**2+dot(v2n,nvref2)**2)-sqrt(2*error**2))*10000), ' total : ', nresul

                        resulfile.write('\n\nHomothety factor (matrix norme) : '+str(norme3) )
                        resulfile.write('\nQuality -> ' + str(dot(v1n,nvref)) +' '+ str(dot(v2n,nvref2)) )
                        resulfile.write('\nRotation matrice :' )
                        resulfile.write('\n'+str(v1))
                        resulfile.write('\n'+str(v2))
                        resulfile.write('\n'+str(v3))

                  z2=z2+1

                y2=y2+1

              x2=x2+1

        z=z+1

      y=y+1

    x=x+1

if nresul!=0 :
  print '\n >>',nresul, ' Rotation matrices were found with criterias limit = ',limit,' and tolerance = ',error,'<< \n'
else :
  print "You may consider increasing the searching range ('limit') or decrease your tolerance expectations ('error') "