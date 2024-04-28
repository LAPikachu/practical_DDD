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

void look_seg(	int *iseg,
		int *ivec,
		int x_ptr,
		int y_ptr,
		int *XO,
		int *YO,
		int *X1,
		int *Y1,
		char *tr,
		int *junct)

/***********************************************************************
*******
******* Ce programme recherche les caracteristiques du segment de dislo-
******* cation pointe par la souris.Pour ce faire,on calcule le cosinus
******* de l'angle dont le sommet est le pointeur de souris et les cotes
******* formes par le pointeur et les extremites du segment vise.
*******
***********************************************************************/

{
int i,j,id,*l,id1,id2,id3;
double x0,x1,y0,y1,ps,mod0,mod1,r;

j = 0;
for (i=0;i<*nseg-12*nbube*nbube*nbube;i++)
   {    /* boucle sur l'ensemble des */
      if (*(tr+i) != 0) continue;  /* segments examines un par un */
      x0=XO[i]-x_ptr;
      x1=X1[i]-x_ptr;
      y0=YO[i]-y_ptr;
      y1=Y1[i]-y_ptr;;
      ps=x0*x1+y0*y1;
      mod0=x0*x0+y0*y0;
      mod0=sqrt(mod0);
      mod1=x1*x1+y1*y1;
      mod1=sqrt((double)mod1);
      r=ps/mod0/mod1;

      if ( r< -.995)
      { /* angle voisin de 180 degres ? -> oui */
         /* un segment est proche */
         /*le parametre -.995 a ete choisi comme bon compromis entre une pre-
            cision raisonnable et un faible nombre d'echecs a la souris */

         id= *(iseg+4+5*i);    /* !!! modif !!! */
         id3= *(iseg+2+5*i);  /* !!! modif !!! */

         l = ivec+id*3 ;
         if (junct[i])
            fprintf(stderr," JUNCTION OF DIRECTION ");
         else
            {
            if ( (id-1)%3 )
               fprintf(stderr," EDGE OF  DIRECTION ");
            else
               fprintf(stderr," SCREW OF DIRECTION ");
            }
         fprintf(stderr,"%d %d %d\n",*(l-3),*(l-2),*(l-1));

         id1= *(iseg+3+5*i);  /* !!! modif !!! */
         id2= *(iseg+4+5*i);  /* !!! modif !!! */

         fprintf(stderr,"Norme = %d Veclin = %d\n",id1,id2);

         id1= *(iseg+5*i);    /* !!! modif !!! */
         id2= *(iseg+1+5*i);  /* !!! modif !!! */

         fprintf(stderr," X = %d Y = %d Z  %d\n",id1,id2,id3);

         fprintf(stderr,"segment number : %d \n",i+1);
         j= 1 ;
      }
   }
if (j == 0 ) fprintf(stderr,"NOTHING FOUND \n");
/* echec de la recherche */
      }
/**********************************************************************/
