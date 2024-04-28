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
extern int h,k,l;
extern unsigned int taille_fenetre;
extern int bufnseg,buf2nmax,enmic;
extern int lame_matic;
float	d,e,dref,eref;

/**********************************************************************/

void lame(	char *tr,			/* parametre commandant le trace des segments */
		float *x,float *y,float *z,	/* tableau du nuage de points representatif des segments de dislocation */
		int *NMAX)

/***********************************************************************
*******
******* ce programme decoupe l'echantillon en une lame mince.Il ne trace
******* que les segments qui se trouvent dans la lame mince ou qui la
******* traverse.A chaque element est attache un parametre tr qui pro-
******* voque ou non le trace.On entre le nuage de points representatifs
******* de l'origine et de l'extremite des segments de dislocation.
******* La normale a la lame possede des indices de Miller h,k,l qui
******* sont demandes au moment de la fabrication de la lame.
******* Est egalement demandee l'epaisseur de la lame  ainsi que se dis-
******* tance au centre.
*******
***********************************************************************/

{
      float mod,uh,uk,ul,xt,yt,zt,ps1,ps2,d1,d2,k1,k2;
      char p;
      int i,j;

      if (!lame_matic) { /* n'est demande qu'une fois au moment de
                           la fabrication de la lame */
      fprintf(stderr,"\n ENTER THE SPECIMEN THICKNESS (in Microns) : ");
      while (scanf("%f",&eref) != 1 ) getchar();
      fprintf(stderr,"\n ENTER THE DISTANCE FROM THE SPECIMEN CENTER (in Microns) : ");
      while (scanf("%f",&dref) != 1 ) getchar();
      }

      e=(eref*taille_fenetre*enmic)/(buf2nmax*1.0);
      /* (buf2nmax*1.0) remplace (float)*NMAX */
      d=(dref*taille_fenetre*enmic)/(buf2nmax*1.0);
      /* d'ou sort ce foutu 385 ??? */

      mod = sqrt ((1.0*h*h)+(1.0*k*k)+(1.0*l*l));
            /*(float)sqrt((double)(h*h+k*k+l*l));*/
/*printf("eref %f dref %f MOD %f \n",eref,dref,mod);	*/
      uh = (h*1.0)/mod;
      uk = (k*1.0)/mod;
      ul = (l*1.0)/mod;
/*printf("h %f k %f l %f \n",uh,uk,ul);	*/
/* printf("1 d %f e %f \n",d,e);		 */
      d1 = d-(e*0.5);
/* printf("d1 %f\n",d1);
printf("2 d %f e %f \n",d,e);	*/
      d2 = d+(e*0.5);
/* printf("d2 %f\n",d2);
printf("bufnseg %d \n",bufnseg);
printf("buf2nmax %d \n",buf2nmax);
printf("enmic %d \n",enmic); */

      for (i=0;i<bufnseg;i++) { /* on va determiner si un segment se
situe dans ou hors de la lame ou bien encore traverse la lame,ce qui
permet de donner une valeur en consequence au parametre tr */

      j=i+bufnseg;

/*	printf("nb %d %d %d -> %d\n",i,bufnseg,j,*(tr+i));
      printf("i %d x %f y %f z %f\n",i,x[i],y[i],z[i]);
      printf("j %d x %f y %f z %f\n",j,x[j],y[j],z[j]); */

      p=1;
      ps1 = x[i]*uh+y[i]*uk+z[i]*ul;/*lorsqu'un segment traverse la*/
      if (ps1 > d1) p=2;	   /*lame,on calcule le ou les points*/
      if (ps1 > d2) p=3;	   /*d'emergence et un nouveau segment*/


      ps2 = x[j]*uh+y[j]*uk+z[j]*ul;/*incluant ce ou ces points est */
      if (ps2 > d1) p += 3;       /*cree qui remplace le segment de */
      if (ps2 > d2) p += 3;	    /* depart */


      k1 = (d1-ps1)/(ps2-ps1+0.000001); /*le parametre p traduit tous les cas*/
      k2 = (d2-ps1)/(ps2-ps1+0.000001); /*de figure */

      if ((p == 1)||(p == 9)) *(tr+i)=1;/*on ne trace pas le segment*/
      if (p == 4) {
      x[i]=x[i]+k1*(x[j]-x[i]);/* on a fait appel a de la geometrie*/
      y[i]=y[i]+k1*(y[j]-y[i]);/*du plan et de la droite de la plus */
      z[i]=z[i]+k1*(z[j]-z[i]);}/* banale */
      if (p == 2) {
      x[j]=x[i]+k1*(x[j]-x[i]);
      y[j]=y[i]+k1*(y[j]-y[i]);
      z[j]=z[i]+k1*(z[j]-z[i]);}
      if (p == 8) {
      x[j]=x[i]+k2*(x[j]-x[i]);
      y[j]=y[i]+k2*(y[j]-y[i]);
      z[j]=z[i]+k2*(z[j]-z[i]);}
      if (p == 6) {
      x[i]=x[i]+k2*(x[j]-x[i]);
      y[i]=y[i]+k2*(y[j]-y[i]);
      z[i]=z[i]+k2*(z[j]-z[i]);}
      if (p == 7) {
      xt=x[i]+k1*(x[j]-x[i]);
      yt=y[i]+k1*(y[j]-y[i]);
      zt=z[i]+k1*(z[j]-z[i]);
      x[j]=x[i]+k2*(x[j]-x[i]);
      y[j]=y[i]+k2*(y[j]-y[i]);
      z[j]=z[i]+k2*(z[j]-z[i]);
      x[i]=xt;
      y[i]=yt;
      z[i]=zt;}
      if (p == 3) {
      xt=x[i]+k1*(x[j]-x[i]);
      yt=y[i]+k1*(y[j]-y[i]);
      zt=z[i]+k1*(z[j]-z[i]);
      x[i]=x[i]+k2*(x[j]-x[i]);
      y[i]=y[i]+k2*(y[j]-y[i]);
      z[i]=z[i]+k2*(z[j]-z[i]);
      x[j]=xt;
      y[j]=yt;
      z[j]=zt;}
/* 	printf("-> i x %f y %f z %f\n",x[i],y[i],z[i]);
      printf("-> j x %f y %f z %f\n",x[j],y[j],z[j]);
      printf("nb %d %d %d -> %d %d\n",i,bufnseg,j,p,*(tr+i)); */
}
}
/**********************************************************************/
