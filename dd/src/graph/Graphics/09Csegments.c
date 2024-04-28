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

extern int nseg[1];

/****************************/
/* Equations diophantiennes */
extern	int dioph;
extern int diophySEG;
extern	int *dstabSEG,*dbtabSEG;
/* Equations diophantiennes */

void segments(	float *x,float *y,
		int *XO,int *YO,int *X1,int *Y1)

/***********************************************************************
*******
******* Ce programme transforme les coordonnees des points projetes,x,y
******* en un ensemble de coordonnees segment d'origine XO,YO et d'ex-
******* tremite X1,Y1.
******* nseg represente le nombre de segments a tracer.
*******
***********************************************************************/

{

      int i,j;
      int idio,*dtabSEG;

/* fprintf(stderr,"SEGMENTS\n"); */

if(dioph==0)
{
      for (i=0;i<*nseg;i++) {
      j=i+*nseg;
      XO[i]=x[i];
      YO[i]=y[i];
      X1[i]=x[j];
      Y1[i]=y[j];		}
}
else
{
      dtabSEG=dbtabSEG;
      for (idio=0;idio<diophySEG;idio++) {
      i=abs(dtabSEG[idio]);
      j=i+*nseg;
      XO[i]=x[i];
      YO[i]=y[i];
      X1[i]=x[j];
      Y1[i]=y[j];		}
}

}
/**********************************************************************/
