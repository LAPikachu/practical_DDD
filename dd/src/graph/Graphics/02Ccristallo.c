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

extern int h,k,l;

/**********************************************************************/

void cristallo(ah,av)

/***********************************************************************
*******
*******	Ce programme permet de passer des indices de Miller h,k,l d'une
******* direction cristallographique aux angles de rotation horizontale
******* et de rotation verticale ah et av qui servent a la projection
******* des points de l'objet a representer.Lorsque les angles de
******* rotation sont nuls,l'axe Ox d'un repere Oxyz est horizontal par
******* rapport a l'observateur et oriente de gauche a droite.L'axe
******* Oy est egalement horizontal et est oriente dans le sens trigono-
******* metrique.L'axe Oz est alors vertical et pointe vers le haut par
******* rapport a l'observateur.L'angle horizontal est pris autour de Oz
******* dans le sens trigonometrique.L'angle vertical est compte autour
******* de Ox,sens positif vers l'observateur.
*******
***********************************************************************/

      float *ah,*av;  /* les angles de rotation horizontale
                           et de rotation verticale */

{
      double mod,pi;
      int  i;

      fprintf(stderr,"\n h k l --> %d %d %d , after dioph \n",h,k,l);

      pi = 180./3.1415926536;

      mod=(double)h*(double)h+(double)k*(double)k;   /* conversion classique sans commentaires */
      mod=sqrt(mod);

      i=0;
      if(h == 0) i=2;
      if((h == 0) && (k == 0)) i=1;

      switch(i) {

                  case 1:
                        *ah = -90.;
                        *av = 90.;
                        break;

                  case 2:
                        if(k>0) *ah = 180.;
                        else *ah= 0.;
                        *av = atan((double)l/mod)*pi;
                        break;

                  default:
                        *ah = -(atan((double)k/(double)h)*pi+90.);
                        *av = atan((double)l/mod)*pi;
                        if(h<0) *ah = *ah+180.;
                        break;
                  }

      fprintf(stderr,"\n ah av %f %f \n",*ah,*av);

}

/**********************************************************************/
