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
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>
#include <stdlib.h>
#include "cooleur.h"
#include "extercool.h"
#include "external.h"

extern  Display	*dpy;
extern  Window	win;
extern  XEvent	event;
extern  unsigned int taille_fenetre;
extern  float e,d,eref,dref;
extern  float *bx_in,*by_in,*bz_in,*bx,*by;
extern  char *btr;
extern  int *bX0,*bY0,*bX1,*bY1;

/*float a_v=12.604383 ,a_h=-116.565048,echelle=20000.,Rzoom=1.0; */
// float a_v= 0.000,a_h= -135.00000,echelle=20000.,Rzoom=1.1;
float a_v= 10.000,a_h= 10.00000,echelle=20000.,Rzoom=1.1;
float pdv_x=0.,pdv_y= -10500.,pdv_z=1.0;
int h=1,k=1,l=1,extinct=1,delai=1,sysok=1;

/****************************/
/* Equations diophantiennes */
/* int dh=1,dk=-1,dl=0; */
int dh=1,dk=1,dl=1;
int x_centre=0,y_centre=0,z_centre=0; /* pour preciser le plan a visualier */
int dioph=0;
int dixmax=17,diymax=17,dizmax=22;
int dixmyx=0,diymyx=0,dizmyx=0;
int dixmix=0,diymix=0,dizmix=0;
int ddixmyx=0,ddiymyx=0,ddizmyx=0;
int ddixmix=0,ddiymix=0,ddizmix=0;
int TDx[4],TDy[4],TDz[4],nbT=0;
float dD[4];
float nvec1=0,nvec2=0,dPlan=0;
float diodistmin,diodist,diox,dioy,dioz;
float diodystmin;
long long int dioprod=0;

int Lx=0,Ly=0,Lz=0,IndiceMax=0,dIndx=0,dIndy=0,dIndz=0;
int dxy=0,dxz=0,dyx=0,dyz=0,dzx=0,dzy=0,Plan=0;

int tdiplx,tdiply,tdiplz;
int diplx[20000],diply[20000],diplz[20000];

float dmax=0;

long long int hmax,hmin,deltah;
int hmaxint, hminint, deltahint;
int dInd;
int O1,O2,O3; /* pour definir proprement les numeros de plan */

int TailDom;

float diopZ;

/* Equations diophantiennes */
/****************************/

/**********/
/* Images */
/* int imx=0,imy=0,imz=0; */
/* images */
/**********/

int nseg[1],tnseg[1],tnmax[1000],nsmax[1],bufnmax[1],bufnseg,buf2nmax,enmic,matic=1000000,lame_matic=0,new=1,segnum=-1;
/*------------variables statics communes aux deux fonctions-----------*/

      static int i,win_x,win_y,*XO,*YO,*X1,*Y1;

/* win_i et win_y coordonnees du pointeur souris sur l'ecran */
/* XO,YO et X1,Y1 sont respectivement les coordonnees origine et extremite des segments projetes a tracer */

      static int co=2,numfilm=0,makefilm=0,prsegnum=0;
      static float *x,*y,*x_in,*y_in,*z_in; /* x_in,y_in et z_in sont des
coordonnees qui doublent xp,yp et zp et qui servent a restaurer le
cristal entier apres fabrication d'une lame mince sans avoir recours
a un fichier */
      static char *tr; /* parametre qui autorise le trace (les portions de
segments invisibles dans le cas d'un lame mince ne sont evidemment pas
traces */


/*-----------------fin des variables statics--------------------------*/

/**********************************************************************/

int draw(	float *xp,
		float *yp,
		float *zp,
		int *iseg,
		int *ivec,
		int *quit,
		int *junct,
		int *nmax,
		int *fini)

/***********************************************************************
*******
*******	draw est le programme qui pilote le l'affichage graphique.Il est
******* appele par le programme predraw lui-meme appele depuis le pro-
******* gramme FORTRAN.Il a deux modes de fonctionnement:un mode automa-
******* tique qui permet un retour au FORTRAN apres le trace graphique
******* et un mode en boucle infinie qui permet un examen du trace.On
******* oeut cependant coupler les deux modes.
*******
***********************************************************************/

/* iseg et ivec sont des tableaux bi-dimen-
                           sionnes issus du programme FORTRAN au
                           travers du pregramme predraw et qui con-
                           tiennent les informations relatives a
                           l'ensemble des dislocations */
/* junct : tableau de boolean issu du Fortran
               specifiant les jonctions */
/* *xp,*yp,*zp; coordonnees de l'origine des segments de
                              dislocation */
{

int tst;
/* on effectue une allocation */
/* dynamique de memoire a la  */
/* main.Cela est necessite par*/
/* le fait que le programme   */
/* FORTRAN est tres gourmand  */
/* et interdit la coexistence */
/* du calcul et du graphisme. */
/* C'est la raison pour laquel*/
/* le la memoire est restituee*/
/* des que l'on quitte draw   */
/*	j = *nseg*8;
      x_in=(float *)malloc(j);
      y_in=(float *)malloc(j);
      z_in=(float *)malloc(j);
      x=(float *)malloc(j);
      y=(float *)malloc(j);
      tr=(char *)malloc(*nseg);
      j= *nseg*4;
      XO=(int *)malloc(j);
      YO=(int *)malloc(j);
      X1=(int *)malloc(j);
      Y1=(int *)malloc(j);
*/
x_in=bx_in;
y_in=by_in;
z_in=bz_in;
x=bx;
y=by;
tr=btr;
XO=bX0;
YO=bY0;
X1=bX1;
Y1=bY1;

for (i=0;i<*nseg;i++) *(tr+i)=0;
for (i=0;i<*nseg*2;i++) {
         *(x_in+i) = *(xp+i);
         /* xp pointe sur sauvegarde*/
         *(y_in+i) = *(yp+i);/* x  pointe sur affichage*/
         *(z_in+i) = *(zp+i);
         /*printf("i %d & %i x %f y %f z %f\n",i,&x_in[i],x_in[i],y_in[i],z_in[i]);
         printf("i %d & %i x %f y %f z %f\n",i,&xp[i],xp[i],yp[i],zp[i]);*/
} /* recopie */

XClearWindow(dpy,win);	/* effacement de la fenetre graphique */


if (lame_matic) lame(tr,x_in,y_in,z_in,nmax); /* test de la demande
                                                d'un trace en lame mince: si le parametre lame_matic est vrai,une lame
                                                mince est fabriquee a chaque passage dans draw,ce qui permet de suivre
                                                une evolution en lame mince */

perspective(x_in,y_in,z_in,x,y); /* projection avec perspective
                                    d'un brouillard de points representant les origines et les extremites
                                    des segments a projeter.On entre en tridimensionnel (x_in,y_in,z_in et
                                    on ressort en bidimensionnel projete (x,y) */

segments(x,y,XO,YO,X1,Y1); /* regroupement des points en seg-*/
/* ments                          */

if(recool==0){
         recool=1;
         getcool(&modecool,&fondcool,&boitecool,trcool,trR,trG,trB,tron,&cooljonc,&cooldev,&nbube);
}

trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);   /* trace des segments proprement dit*/
if(makefilm==1)traceps(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);   /* trace epsf */

if (matic != 0)
         { /* selon la valeur de matic,on repart au FOR
                              TRAN moyennant un delai permettant d'effectuer la totalite du trace ou
                              bien on passe en boucle infinie,ce qui est toujours le cas au premier
                              pas (ISTEP du FORTRAN = 1 ) */
               for (i=0;i<delai;i++)
                        trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);
               while (XPending(dpy))
                        testevent(xp,yp,zp,iseg,ivec,quit,junct,nmax,fini);
               if (matic != *fini)
                        {
                                 /* restitution de la memoire avant de quitter */
                                 /* le programme draw			      */
                                 /*		    free(x_in);
                                       free(y_in);
                                       free(z_in);
                                       free(x);
                                       free(y);
                                       free(tr);
                                       free(XO);
                                       free(YO);
                                       free(X1);
                                       free(Y1);
                                 */
                  return (0) ;
                  }
               }

      for (;;) {
      /* mise du programme en boucle infinie */
               if ((tst=testevent(xp,yp,zp,iseg,ivec,quit,junct,nmax,fini)))
                        {
                          if (tst==1)
                              return (0);
                          else
                              return(1);
                                    /* sortie de la boucle infinie   */
                                    /* retour au FORTRAN via predraw */
                        }
               }
}

/**********************************************************************/

int testevent(	float *xp,
		float *yp,
		float *zp,
		int *iseg,
		int *ivec,
		int *quit,
		int *junct,
		int *nmax,
		int *fini)

/* iseg et ivec sont des tableaux bi-dimen-
                     sionnes issus du programme FORTRAN au
                     travers du pregramme predraw et qui con-
                     tiennent les informations relatives a
                     l'ensemble des dislocations */
/* junct : tableau de boolean issu du Fortran
               specifiant les jonctions */
/* coordonnees de l'origine des segments de
                     dislocation */
/* procedure qui teste les events : clavier, fenetre... et qui reagit */
{

/****************************/
/* Equations diophantiennes */
int ix_test=0,iy_test=0,iz_test=0;
int ix=0,iy=0,iz=0,itest=0,numDuPlan;
int posx=0,posy=0,posz=0;
long long int vect_1=0,vect_2=0,vect_3=0;
long long int lprodsca=0;
float dmin=0,dtest=0;
/* Equations diophantiennes */
/****************************/

char cas;


      XNextEvent(dpy,&event); /* recherche des evenements pouvant
                              intervenir dans la fenetre graphique:pointage souris,deplacement de
                              la fenetre (auquel cas il faut ordonner un nouveau trace),appui sur
                              une touche camclavier */
      switch (event.type)
               {
               case ConfigureNotify : /* changement taille fenetre */
               if (new == 0)
                        {
                        new =1;
                        return (2);
                        }

if ( event.xconfigure.width > taille_fenetre &&
     event.xconfigure.height > taille_fenetre ) {
        if (event.xconfigure.width < event.xconfigure.height)
          taille_fenetre = event.xconfigure.width;
        else
          taille_fenetre = event.xconfigure.height;
 }
               XResizeWindow(dpy,win,taille_fenetre,taille_fenetre);
               new = 0;
               return (2);
               break;

               case Expose : /* fenetre recouverte puis decouverte:
                             on effectue un nouveau trace */
               if(event.xexpose.count == 0 )
                        {
                        while(XCheckTypedEvent(dpy,Expose,&event));
                        /*trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);*/
                        /* It seems that there is a problem with this case
                        for the moment the corresponding case is commented */
                        }
               return (0);
               break;

               case ButtonPress : /* appui sur un bouton de la
                                    souris */
               if (event.xbutton.button == Button1) /* plus
               exactement,le bouton gauche */
                        {
                        win_x = event.xbutton.x; /* determination */
                        win_y = event.xbutton.y; /* de la position*/
                                                /* de la souris  */

                        fprintf(stderr,"\n MicroMegas> ");
                        scanf("%s",&cas);

                        switch (cas) /*examen de toutes les
                                       possibites d'appui de touche */
                        {

                        case 'x' :      /* on determine s'il n'existe pas un segment de dislocation situe
                                        pres de l'endroit pointe par la souris:si ce segment existe,ses
                                        caracteristiques sont affichees sur la sortie erreur standard */
                        look_seg(iseg,ivec,win_x,win_y,XO,YO,X1,Y1,tr,junct);
                        return(0);
                        break;

                        /* Infos sur un segment*/
                        case 'i' : /* appui sur la touche i */
                        if (segnum==-1)
                        {
                           fprintf(stderr,"\n Segment Number :");
                           while (scanf("%d",&segnum) != 1) getchar();
                        }
                        else
                        {
                           segnum=-1;
                        }
                        if (segnum!=-1) unsegout(iseg,junct,segnum);
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);
                        return(0);
                        break;

/****************************************************************************************/
/* Equations diophantiennes */
                        case 'j' : /* appui sur la touche J */

/*fprintf(stderr,"\n  3  5 %d",modulo( 3, 5));
fprintf(stderr,"\n  6  5 %d",modulo( 6, 5));
fprintf(stderr,"\n -3  5 %d",modulo(-3, 5));
fprintf(stderr,"\n -6  5 %d",modulo(-6, 5));
fprintf(stderr,"\n  3 -5 %d",modulo( 3,-5));
fprintf(stderr,"\n  6 -5 %d",modulo( 6,-5));
fprintf(stderr,"\n -3 -5 %d",modulo(-3,-5));
fprintf(stderr,"\n -6 -5 %d",modulo(-6,-5));*/


                        fprintf(stderr,"\n h k l (0 0 0 pour cube) ? : ");
                        while (scanf("%d%d%d",&dh,&dk,&dl) != 3 ) getchar();
                        if(dh==0 && dk==0 && dl==0)
                        {
                        dioph=0;
                        Rzoom=1.0;
                        extinct = 1;
                        fprintf(stderr,"\n LE MODE GRAPHIQUE BASCULERA EN VOLUME AU PROCHAIN TRACE...\n");
                        }
                        else
                        {
                        dioph=1;
                        fprintf(stderr,"\n CENTRE (%d %d %d) ? : ",x_centre,y_centre,z_centre);
                        while (scanf("%d%d%d",&x_centre,&y_centre,&z_centre) != 3 ) getchar();
                        fprintf(stderr,"\n Lx %d Ly %d Lz %d",Lx,Ly,Lz);
                        fprintf(stderr,"\n dxy %d dxz %d",dxy,dxz);
                        fprintf(stderr,"\n dyx %d dyz %d",dyx,dyz);
                        fprintf(stderr,"\n dzx %d dzy %d",dzx,dzy);
                        IndiceMax=abs(Lx*dh)+abs(Ly*dk)+abs(Lz*dl);
hmin=0;
if(dh<0){hmin=hmin+dh*Lx;O1=Lx;/*F1=0;*/}
if(dk<0){hmin=hmin+dk*Ly;O2=Ly;/*F2=0;*/}
if(dl<0){hmin=hmin+dl*Lz;O3=Lz;/*F3=0;*/}
hmax=0;
if(dh>=0){hmax=hmax+dh*Lx;O1=0;/*F1=Lx;*/}
if(dk>=0){hmax=hmax+dk*Ly;O2=0;/*F2=Ly;*/}
if(dl>=0){hmax=hmax+dl*Lz;O3=0;/*F3=Lz;*/}
deltah=(hmax-hmin);
hmaxint = hmax;
hminint = hmin;
deltahint = deltah;
                        fprintf(stderr,"\n O : (%d,%d,%d)",O1,O2,O3);
                        fprintf(stderr,"\n Indice minimum : %d",hminint);
                        fprintf(stderr,"\n Indice maximum : %d",hmaxint);
                        fprintf(stderr,"\n dIndice : %d",deltahint);
                        dIndx=Lx*dh+dyx*dk+dzx*dl;
                        dIndy=Ly*dk+dxy*dh+dzy*dl;
                        dIndz=Lz*dl+dxz*dh+dyz*dk;
                        dInd=pgdc3(dIndx,dIndy,dIndz);
                        fprintf(stderr,"\n Increment d'indice : %d",dInd);
if(dh!=0){tdiplx=dInd/dh;}else{tdiplx=0;}
if(dk!=0){tdiply=dInd/dk;}else{tdiply=0;}
if(dl!=0){tdiplz=dInd/dl;}else{tdiplz=0;}
/*if(tdiplx*dh==dInd){tdiply=0;tdiplz=0;}
else{if(tdiply*dk==dInd){tdiplx=0;tdiplz=0;}
else{if(tdiplz*dl==dInd){tdiplx=0;tdiply=0;}
else{fprintf(stderr,"\n PROBLEME POUR DETERMINER LE VECTEUR DE TRANSITION ENTRE PLANS");}}}*/

                        Plan=deltah/dInd;
                        fprintf(stderr,"\n Nombre de plans : %d",Plan);
                        TailDom=sqrt(Plan);
                        dPlan=sqrt((1.0*Lx*Lx)+(1.0*Ly*Ly)+(1.0*Lz*Lz))/Plan;
                        fprintf(stderr,"\n Distance entre plan : %f",dPlan);
/*fprintf(stderr,"\n t : %d %d %d",tdiplx,tdiply,tdiplz);*/

                        fprintf(stderr,"\n RECHERCHE DE LA TAILLE DU PAVAGE...%d",TailDom);

resolution:;
TDx[0]=0;TDx[1]=0;TDx[2]=0;TDx[3]=0;
TDy[0]=0;TDy[1]=0;TDy[2]=0;TDy[3]=0;
TDz[0]=0;TDz[1]=0;TDz[2]=0;TDz[3]=0;
nbT=0;
diox=Lx;
dioy=Ly;
dioz=Lz;
dmax = sqrt((diox*diox)+(dioy*dioy)+(dioz*dioz))*TailDom;
dD[0] = dmax;dD[1] = dmax;dD[2] = dmax;dD[3] = dmax;

                        for(dixmax=-TailDom;dixmax<TailDom;dixmax++){
                        for(diymax=-TailDom;diymax<TailDom;diymax++){
                        for(dizmax=-TailDom;dizmax<TailDom;dizmax++){

diox=(dixmax*Lx+diymax*dxy+dizmax*dxz);
dioy=(diymax*Ly+dixmax*dyx+dizmax*dyz);
dioz=(dizmax*Lz+dixmax*dzx+diymax*dzy);
diodist=sqrt((diox*diox)+(dioy*dioy)+(dioz*dioz));
dioprod=(dixmax*(Lx*dh+dyx*dk+dzx*dl)+
         diymax*(Ly*dk+dxy*dh+dzy*dl)+
         dizmax*(Lz*dl+dxz*dh+dyz*dk));

if(dioprod==0&&diodist!=0)
                        {
if(diodist<dD[0]&&abs(dixmax)!=abs(TDx[0])&&abs(diymax)!=abs(TDy[0])&&abs(dizmax)!=abs(TDz[0]))
{
if(nbT>=4)
{
TDx[3]=TDx[2];
TDy[3]=TDy[2];
TDz[3]=TDz[2];
dD[3]=dD[2];
}
if(nbT>=3)
{
TDx[2]=TDx[1];
TDy[2]=TDy[1];
TDz[2]=TDz[1];
dD[2]=dD[1];
}
if(nbT>=2)
{
TDx[1]=TDx[0];
TDy[1]=TDy[0];
TDz[1]=TDz[0];
dD[1]=dD[0];
}
TDx[0]=dixmax;
TDy[0]=diymax;
TDz[0]=dizmax;
dD[0]=diodist;
nbT++;
diodist=dmax;
}
if(diodist<dD[1]&&abs(dixmax)!=abs(TDx[1])&&abs(diymax)!=abs(TDy[1])&&abs(dizmax)!=abs(TDz[1]))
{
if(nbT>=4)
{
TDx[3]=TDx[2];
TDy[3]=TDy[2];
TDz[3]=TDz[2];
dD[3]=dD[2];
}
if(nbT>=3)
{
TDx[2]=TDx[1];
TDy[2]=TDy[1];
TDz[2]=TDz[1];
dD[2]=dD[1];
}
TDx[1]=dixmax;
TDy[1]=diymax;
TDz[1]=dizmax;
dD[1]=diodist;
diodist=dmax;
}
if(diodist<dD[2]&&abs(dixmax)!=abs(TDx[2])&&abs(diymax)!=abs(TDy[2])&&abs(dizmax)!=abs(TDz[2]))
{
if(nbT>=4)
{
TDx[3]=TDx[2];
TDy[3]=TDy[2];
TDz[3]=TDz[2];
dD[3]=dD[2];
}
TDx[2]=dixmax;
TDy[2]=diymax;
TDz[2]=dizmax;
dD[2]=diodist;
diodist=dmax;
}
if(diodist<dD[3]&&abs(dixmax)!=abs(TDx[3])&&abs(diymax)!=abs(TDy[3])&&abs(dizmax)!=abs(TDz[3]))
{
TDx[3]=dixmax;
TDy[3]=diymax;
TDz[3]=dizmax;
dD[3]=diodist;
diodist=dmax;
}

               }
                        }
                        }
                        }

if(TDx[3]==0&&TDy[3]==0&&TDz[3]==0)
{
TailDom=TailDom+10;
fprintf(stderr,"\n Problemes lors de la resolution des equations diophantiennes on augmente le domaine explore %d",TailDom);
goto resolution;
}
else
{
fprintf(stderr,"\n Taille (profondeur) %d %d %d",TDx[0],TDy[0],TDz[0]);
fprintf(stderr,"\n Taille (profondeur) %d %d %d",TDx[1],TDy[1],TDz[1]);
fprintf(stderr,"\n Taille (profondeur) %d %d %d",TDx[2],TDy[2],TDz[2]);
fprintf(stderr,"\n Taille (profondeur) %d %d %d",TDx[3],TDy[3],TDz[3]);
dixmix=mymax(TDx[1],TDx[0]);
dixmyx=mymax(TDx[3],TDx[2]);
diymix=mymax(TDy[1],TDy[0]);
diymyx=mymax(TDy[3],TDy[2]);
dizmix=mymax(TDz[1],TDz[0]);
dizmyx=mymax(TDz[3],TDz[2]);
dixmax=mymax(dixmix,dixmyx);
diymax=mymax(diymix,diymyx);
dizmax=mymax(dizmix,dizmyx);
diodystmin=dD[3];
/*ddixmix=dixmix-dixmyx;
ddiymix=diymix-diymyx;
ddizmix=dizmix-dizmyx;

ddixmyx=dixmix+dixmyx;
ddiymyx=diymix+diymyx;
ddizmyx=dizmix+dizmyx;*/

fprintf(stderr,"\n Taille (profondeur) %d %d %d",dixmax,diymax,dizmax);
}
fprintf(stderr,"\n TABULATION POUR LES PLANS DE LA BOITE...");

/*------------------------------------------------------------------------------------*/
   for(i=0;i<10000;i++){diplx[i]=0;diply[i]=0;diplz[i]=0;}
   posx=O1;
   posy=O2;
   posz=O3;
   for(i=0;i<=Plan;i++){
      ix_test=0;
      iy_test=0;
      iz_test=0;
      itest = 0;
      if(i!=0){
      posx=posx+tdiplx;
      if(posx<0 ){posx=posx-tdiplx;tdiplx=0;}
      if(posx>Lx){posx=posx-tdiplx;tdiplx=0;}
      if(tdiplx==0){
      posy=posy+tdiply;
      if(posy<0 ){posy=posy-tdiply;tdiply=0;}
      if(posy>Ly){posy=posy-tdiply;tdiply=0;}
      if(tdiply==0){
      posz=posz+tdiplz;
      if(posz<0 ){posz=posz-tdiplz;tdiplz=0;}
      if(posz>Lz){posz=posz-tdiplz;tdiplz=0;}
      }}}
      numDuPlan=(((posx-O1)*dh)+((posy-O2)*dk)+((posz-O3)*dl))/(1.*dInd);
      dmin = dmax;
      for(ix=-dixmax;ix<=dixmax;ix++){for(iy=-diymax;iy<=diymax;iy++){for(iz=-dizmax;iz<=dizmax;iz++)
      {
      vect_1 = posx + Lx*ix + iy*dxy +iz*dxz;
      vect_2 = posy + Ly*iy + ix*dyx +iz*dyz;
      vect_3 = posz + Lz*iz + ix*dzx +iy*dzy;
      lprodsca=(dh*vect_1)+(dk*vect_2)+(dl*vect_3);
      dtest = sqrt((vect_1*vect_1)+ (vect_2*vect_2)+(vect_3*vect_3));

/* fprintf(stderr,"\n ALORS QUOI %f %f \n",lprodsca,dtest); */

      if (lprodsca == 0 && dtest < dmin)
      {
      dmin=dtest;
      ix_test=ix;
      iy_test=iy;
      iz_test=iz;
      itest = 1;
      }
      }}}
      if(itest == 1)
      {
/*       if((i/50)*50==i){
      fprintf(stderr,"\n Pt %d %d %d",posx,posy,posz);}*/
      diplx[numDuPlan]=ix_test;
      diply[numDuPlan]=iy_test;
      diplz[numDuPlan]=iz_test;
/*if(ix_test==dixmix&&iy_test==diymix&&iz_test==dizmix){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}
if(ix_test==dixmyx&&iy_test==diymyx&&iz_test==dizmyx){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}
if(ix_test==ddixmix&&iy_test==ddiymix&&iz_test==ddizmix){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}
if(ix_test==ddixmyx&&iy_test==ddiymyx&&iz_test==ddizmyx){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}

if(ix_test==-dixmix&&iy_test==-diymix&&iz_test==-dizmix){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}
if(ix_test==-dixmyx&&iy_test==-diymyx&&iz_test==-dizmyx){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}
if(ix_test==-ddixmix&&iy_test==-ddiymix&&iz_test==-ddizmix){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}
if(ix_test==-ddixmyx&&iy_test==-ddiymyx&&iz_test==-ddizmyx){diplx[numDuPlan]=0;diply[numDuPlan]=0;diplz[numDuPlan]=0;}*/
      }
      else
      {
      fprintf(stderr,"\n %d %d %d dioph pas resolue pour le plan %d %d\n",posx,posy,posz,numDuPlan,i);
      diplx[numDuPlan]=0;
      diply[numDuPlan]=0;
      diplz[numDuPlan]=0;
      } /* pb de resolution */
      } /* sur les plans */
/*------------------------------------------------------------------------------------*/
                        }
/*		  	XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);*/

fprintf(stderr,"\n LE MODE GRAPHIQUE BASCULERA EN PLAN DEPLIE AU PROCHAIN TRACE...\n");
extinct = 0;
h=dh;
k=dk;
l=dl;
Rzoom=(dmax*0.01)/(diodystmin*1.5);
diopZ=1.0;
cristallo(&a_h,&a_v);
return(0);
break;
/* Equations diophantiennes */
/*******************************************************************************************/

/********************/
/* Zoom foret dioph */
                        case 'u' : /* appui sur la touche U */
/*
                        fprintf(stderr,"\n NOMBRE D'IMAGES (%d %d %d, 0 0 0 pour cube) ? : ",imx,imy,imz);
                        while (scanf("%d%d%d",&imx,&imy,&imz) != 3 ) getchar();
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);
                        return(0);

*/

                        fprintf(stderr,"\n Zoom pour la foret en mode diophantien (%f) ? : ",diopZ);
                        while (scanf("%f",&diopZ) != 1 ) getchar();
                        XClearWindow(dpy,win);
                        trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);
                        return (0);

                        break;
/**********/

                        case 'k' : // appui sur la touche k
                        makefilm=abs(makefilm-1);
                        fprintf(stderr,"\n CAPTURE AUTOMATIQUE : MAKEFILM NUMFILM -> %d %d \n\n",makefilm,numfilm);
                        //numfilm=0;
                        return(0);
                        break;

/*                        case 'w' : // appui sur la touche w
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        //numfilm=999;
                        fprintf(stderr,"\n UN PLAN : MAKEFILM NUMFILM -> %d %d \n\n",makefilm,numfilm);
                        traceps(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);   // trace epsf
                        return(0);
                        break;
*/

                        case 'c' : /* appui sur la touche  "c" */
                        getcool(&modecool,&fondcool,&boitecool,trcool,trR,trG,trB,tron,&cooljonc,&cooldev,&nbube);
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);*/
                        return (0);
                        break;

                        case 'n' : /* appui sur la touche  "c" */
                        prsegnum=abs(prsegnum-1);
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'v' : /* appui sur la touche  "v" */
                        fprintf(stderr,"\n VERTICAL ANGLE(%.1f) ? : ",a_v);
                        while (scanf("%f",&a_v) != 1) getchar();
/* il est demande d'entrer un nouvel angle vertical de visualisation */
/* puis on effectue un nouveau trace complet                         */
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'h' : /* appui sur la touche "h" */
/* meme sequence que pour l'angle vertical                   */
               fprintf(stderr,"\n HORIZONTAL ANGLE(%.1f) ? : ",a_h);
                        while (scanf("%f",&a_h) != 1) getchar();
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'd' : /* on a appuye sur la touche "d" */
/* on peut ainsi visualiser ou bien faire disparaitre alternativement
   (par appui sur la meme touche) un cube de 10 microns de cote. */
                        extinct = !extinct;
                        XClearWindow(dpy,win);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);

                        case 'o' : /* on a appuye sur la touche "o" */
                        segout(iseg,junct);
                        return (0);
                        break;

                        case 'g' : /* appui sur la touche "g" */
/* parametre delai */
                        fprintf(stderr,"\n DELAY (1-1000) :" );
                        while (scanf("%d",&delai) != 1) getchar();
                        return (0);
                        break;

                        case 'f' : /* appui sur la touche "f" */
                        fprintf(stderr,"\n SCALE (%.1f) ? : ",echelle);
                        while (scanf("%f",&echelle) != 1) getchar();
/* changement de la profondeur de perspective et nouveau trace */
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'z' : /* appui sur la touche "z" */
                        fprintf(stderr,"\n Zoom (%f) ? : ",Rzoom);
                        while (scanf("%f",&Rzoom) != 1) getchar();
/* changement de zoom et nouveau trace */
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'p' : /* appui sur la touche "p" */
fprintf(stderr,"\n VIEW POINT (%.0f %.0f %.0f) ? : ",pdv_x,pdv_y,pdv_z);
      while (scanf("%f%f%f",&pdv_x,&pdv_y,&pdv_z) != 3 ) getchar();
/* changement du point de vue et nouveau trace */
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'm' :/*appui sur la lettre "m" pour Miller*/
      fprintf(stderr,"\n DIRECTION INDICES (%d %d %d) ? : ",h,k,l);
                        while (scanf("%d%d%d",&h,&k,&l) != 3 ) getchar();
/* changement de representation par les indices de Miller:
   la vue de l'observateur se fait selon la direction entree */
                        cristallo(&a_h,&a_v);
            fprintf(stderr,"\n a_v = %f    ;   a_h = %f \n\n",a_v,a_h);
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 'a' : /* appui sur la touche "a" */
                                 /* appel de l'aide en ligne */
                        aide();
                        return (0);
                        break;

                        case 'l' : /* appui sur la touche "l" */
/*mise en route d'une lame mince et blocage interdiction d'effectuer
une nouvelle lame avant d'avoir restaure le cristal massif par
appui sur les touches "b" ou "c" */
                        if (!lame_matic)
                              {
                              fprintf(stderr,"\n MILLER INDICES OF THE NORMAL : ");
                              while (scanf("%d%d%d",&h,&k,&l) != 3) getchar();
                              cristallo(&a_h,&a_v);
                              lame(tr,x_in,y_in,z_in,nmax);
                              XClearWindow(dpy,win);
                              perspective(x_in,y_in,z_in,x,y);
                              segments(x,y,XO,YO,X1,Y1);
                              /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax);*/
                              lame_matic=1;
                              }
                        return (0);
                        break;

                        case 'b' : /* appui sur la touche "b" */
                                 /* restauration du cristal massif */
                        lame_matic=0;
                        for (i=0;i<*nseg;i++) *(tr+i)=0;
                        for (i=0;i<*nseg*2;i++)
                              {
                              *(x_in+i) = *(xp+i);
                              *(y_in+i) = *(yp+i);
                              *(z_in+i) = *(zp+i);
                              }
                        return (0);
                        break;

                        case 'e' : /* appui sur la touche "e" */
               /* restauration du cristal massif et nouveau trace */
                        lame_matic=0;
                        for (i=0;i<*nseg;i++) *(tr+i)=0;
                        for (i=0;i<*nseg*2;i++)
                              {
                              *(x_in+i) = *(xp+i);
                              *(y_in+i) = *(yp+i);
                              *(z_in+i) = *(zp+i);
                              }
                        XClearWindow(dpy,win);
                        perspective(x_in,y_in,z_in,x,y);
                        segments(x,y,XO,YO,X1,Y1);
                        /* trace(tr,co,iseg,XO,YO,X1,Y1,segnum,junct,nmax); */
                        return (0);
                        break;

                        case 's' : /* appui sur la touche "s" */
/* on sort de la boucle infinie apres avoir restitue la memoire */
/* ( retour au programme predraw puis au FORTRAN ) */
                        matic=0;
                        return(1);
                        break;

                        case 't' : /*appui sur la touche "t" */
/* on arrete le more revolver */
                        matic=0;
                        fprintf(stderr,"\n\nStepping mode\n\n");
                        return(0);
                        break;
                        case 'r' : /* appui sur la touche "r" */
/* on lance le mode automatique pour un nombre de pas fixe par
   l'operateur */
                        fprintf(stderr," \n AUTOMATIC MODE ");
                        fprintf(stderr," \n NUMBER OF STEPS : ");
                        while (scanf("%d",&matic) != 1) getchar();
                        if (matic<=0) {*quit = matic ;}
                        matic=*fini+matic;
                        return (1) ;
                        break;

                        case 'q' : /* appui sur la touche "q" */
                        *quit = 0 ;
                        return (1) ;
                        break;
                        }
               break;   /* ouf! */

            }
break;
}
return 0;
}
