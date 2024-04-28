/*
!==================================================================================================
! mM (for microMegas) is an open source program of DD (Dislocation Dynamics) simulation
! originally developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
! Copyright (C) 1993, 1996, 2000, 2002, 2004   Benoit Devincre, Ladislas Kubin, Marc Condat,
! Christophe Lemarchand, Ronan Madec, Sebastien Groh, Ghiath Monnet, Julien Durinck,
! Phillipe Carrez, Christophe de Sansal, Mohamed-Gazi Tagorti.
!
! This file is part of mM.
!
! mM is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! any later version.
!
! mM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MM; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
!
! For more information see http://zig.onera.fr/mm_home_page/
!===================================================================================================
*/

#include <stdio.h>
#include <math.h>
#include "external.h"

extern int nseg[1];

/**********************************************************************/

void segout(	int *iseg,
		int *junct)

{
int i,id1,id2,id3,id4,id5,id8;

for (i=0;i<*nseg-12*nbube*nbube*nbube;i++)
   {

         id1= *(iseg+5*i);
         id2= *(iseg+1+5*i);
         id3= *(iseg+2+5*i);
         id4= *(iseg+3+5*i);
         id5= *(iseg+4+5*i);
         id8= *(junct+i);

         fprintf(stderr,"%d : X= %d Y= %d Z= %d L= %d Vlin= %d J= %d\n",i+1,id1,id2,id3,id4,id5,id8);

   }
      }
/**********************************************************************/
void unsegout(	int *iseg,
		int *junct,
		int i)


{
int id1,id2,id3,id4,id5,id8;

         id1= *(iseg+5*(i-1));
         id2= *(iseg+1+5*(i-1));
         id3= *(iseg+2+5*(i-1));
         id4= *(iseg+3+5*(i-1));
         id5= *(iseg+4+5*(i-1));
         id8= *(junct+(i-1));

         fprintf(stderr,"%d : X= %d Y= %d Z= %d L= %d Vlin= %d J= %d\n",i,id1,id2,id3,id4,id5,id8);

}
/**********************************************************************/
