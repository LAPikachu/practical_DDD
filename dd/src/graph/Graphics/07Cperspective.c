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

#include <math.h>
#include "external.h"

extern	int nseg[1];
extern	unsigned int taille_fenetre; /* tout est ajuste a la taille de la
                              fenetre */
extern	float a_v,a_h,pdv_x,pdv_y,pdv_z,echelle;/* pdv,vecteur point de
                                                   vue */
/*a_v angle de rotation verticale et a_h angle de rotation horizontal*/
/* echelle,profondeur de perspective permettant de faire varier la posi-
tion relative des points de fuite par rapport a l'observateur */

/**********************************************************************/

/****************************/
/* Equations diophantiennes */
extern	int dioph;
extern	float *dsx_in,*dsy_in,*dsz_in;
extern	float *dbx_in,*dby_in,*dbz_in;
extern	float *ddsx_in,*ddsy_in,*ddsz_in;
extern	float *ddbx_in,*ddby_in,*ddbz_in;
extern int *dstabSEG,*dbtabSEG;
extern int diophyPLUS,diophySEG;
extern int *dsxp,*dsyp;
extern int *dbxp,*dbyp;
extern int *ddsxp,*ddsyp;
extern int *ddbxp,*ddbyp;
/* Equations diophantiennes */
/****************************/

void perspective(float *x,float *y,float *z,float *xp,float *yp)

/***********************************************************************
*******
******* Ce programme permet d'effectuer la projection plane d'un ensem-
******* ble de points en y incluant un effet de perspective.Les points
******* sont entres a l'aide de leurs coordonnees trirectangulaires x,y
******* z,mesurees par rapport a un repere orthonorme Oxyz.Le plan Oxz
******* est considere comme un plan vertical dont lanormale passe par un
******* observateur situe au point de vue de coordonnees pvx,pvy,pvz
******* egalement mesure a l'aide de Oxyz.L'axe Oz est un axe vertical
******* autour duquel on peut faire tourner l'ensemble des points d'un
******* angle horizontal a_h (en degres).L'axe Ox est un axe horizontal
******* autour duquel on peut faire tourner les points d'un angle verti
******* cal a_v (en degres).
******* L'observateur regarde le long de Oy.
******* Une forte valeur du parametre echelle permet d'obtenir une pers-
******* pective cavaliere (points de fuite a l'infini)
******* nseg*2 est le nombre de points de la figure.
*******
***********************************************************************/

{
      float pi,vc,vs,hc,hs,pf,a,dim;
      float *dx,*dy,*dz;
      int *dxp,*dyp;
      float *ddx,*ddy,*ddz;
      int *ddxp,*ddyp;
      double h,v;
      int i,idio,*dtabSEG;

/* fprintf(stderr,"PERSPECTIVE\n"); */

      pi=3.141592654/180.;
      dim=taille_fenetre/2.;

      h=a_h*pi;
      v=a_v*pi;

      vc=cos(v);
      vs=sin(v);
      hc=cos(h);
      hs=sin(h);

if(dioph==0)
{
      for (i=0;i<(*nseg)*2;i++) {
      a= x[i]*hs+y[i]*hc;
      yp[i]= a*vs+z[i]*vc-pdv_z;
      xp[i]=x[i]*hc-y[i]*hs-pdv_x;
      pf=echelle/(vc*a-z[i]*vs-pdv_y+echelle);
      xp[i] = dim+xp[i]*pf*Rzoom;
      yp[i] = dim-yp[i]*pf*Rzoom;
      }
}
else
{

dtabSEG=dbtabSEG;

dx=dbx_in;
dy=dby_in;
dz=dbz_in;

dxp=dbxp;
dyp=dbyp;

ddx=ddbx_in;
ddy=ddby_in;
ddz=ddbz_in;

ddxp=ddbxp;
ddyp=ddbyp;

/* fprintf(stderr,"PERSPECTIVE diophSEG  %d\n",diophySEG); */

      for (idio=0;idio<(diophySEG);idio++) {

      i=abs(dtabSEG[idio]);

/* fprintf(stderr,"idio i %d %d\n",idio,i); */

      a= x[i]*hs+y[i]*hc;
      yp[i]= a*vs+z[i]*vc-pdv_z;
      xp[i]=x[i]*hc-y[i]*hs-pdv_x;
      pf=echelle/(vc*a-z[i]*vs-pdv_y+echelle);
      xp[i] = dim+xp[i]*pf*Rzoom;
      yp[i] = dim-yp[i]*pf*Rzoom;

      i=i+*nseg;

      a= x[i]*hs+y[i]*hc;
      yp[i]= a*vs+z[i]*vc-pdv_z;
      xp[i]=x[i]*hc-y[i]*hs-pdv_x;
      pf=echelle/(vc*a-z[i]*vs-pdv_y+echelle);
      xp[i] = dim+xp[i]*pf*Rzoom;
      yp[i] = dim-yp[i]*pf*Rzoom;
      }

/* fprintf(stderr,"PERSPECTIVE diophPLUS  %d\n",diophyPLUS); */

      for (i=0;i<(diophyPLUS);i++) {
      a= dx[i]*hs+dy[i]*hc;
      dyp[i]= a*vs+dz[i]*vc-pdv_z;
      dxp[i]=dx[i]*hc-dy[i]*hs-pdv_x;
      pf=echelle/(vc*a-dz[i]*vs-pdv_y+echelle);
      dxp[i] = dim+dxp[i]*pf*Rzoom;
      dyp[i] = dim-dyp[i]*pf*Rzoom;

      a= ddx[i]*hs+ddy[i]*hc;
      ddyp[i]= a*vs+ddz[i]*vc-pdv_z;
      ddxp[i]=ddx[i]*hc-ddy[i]*hs-pdv_x;
      pf=echelle/(vc*a-ddz[i]*vs-pdv_y+echelle);
      ddxp[i] = dim+ddxp[i]*pf*Rzoom;
      ddyp[i] = dim-ddyp[i]*pf*Rzoom;

      }
}

}
/**********************************************************************/
