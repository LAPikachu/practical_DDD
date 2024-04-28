/*
 ! ===================================================================================================
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
#include "extercool.h"

/****************************/
/* Equations diophantiennes */

extern	float *dsx_in,*dsy_in,*dsz_in;
static float *dbx_in,*dby_in,*dbz_in;
extern	float *ddsx_in,*ddsy_in,*ddsz_in;
extern	float *ddbx_in,*ddby_in,*ddbz_in;
extern	int *dstabSEG,*dbtabSEG;
extern	int *dstabPLUS,*dbtabPLUS;

extern	int dh,dk,dl;
extern	int x_centre,y_centre,z_centre;
extern	int dioph;
extern  int dixmax,diymax,dizmax;
extern	int dxy,dxz,dyx,dyz,dzx,dzy;
extern	int Lx,Ly,Lz;
extern	float minmax;

extern	float dmax;

extern	int SatNsegmax;
int diophyPLUS=0,diophySEG=0;

int idioph=0,idiophy=0;

int vl_1=0,vl_2=0,vl_3=0;

#ifdef ATOUM
long int tmp_x=0,tmp_y=0,tmp_z=0;
long int ttmp_x=0,ttmp_y=0,ttmp_z=0;
#else
long long int tmp_x=0,tmp_y=0,tmp_z=0;
long long int ttmp_x=0,ttmp_y=0,ttmp_z=0;
#endif

extern int dInd,Plan;
extern int O1,O2,O3;
extern int diplx[20000],diply[20000],diplz[20000];
extern float diopZ;
/* Equations diophantiennes */
/****************************/

extern  int nseg[];
extern  int tnseg[];
extern  int tnmax[];
extern  int nsmax[];
extern  int bufnmax[];
extern  int bufnseg,buf2nmax,enmic,h,k,l;
extern  unsigned int taille_fenetre;
extern  float e,d,eref,dref;

extern float *kbx_in,*kby_in,*kbz_in;
extern float *kdx_in,*kdy_in,*kdz_in;
extern float dPlan;
/**********************************************************************/

#ifdef MAC
void predraw(n,iseg,tabvois,ivec,fini,NMAX,junct,nbsegdep,nsegmax)
#else
void predraw_(n,iseg,tabvois,ivec,fini,NMAX,junct,nbsegdep,nsegmax)
#endif

/***********************************************************************
*******
******* Ce programme realise l'interface entre le programme FORTRAN de
******* Gilles Canova et le programme graphique.
******* Il a pour but de transferer les tableaux dont on a besoin pour
******* le trace depuis le FORTRAN vers le C.
******* Les parametres de passage ont meme nom en C qu'en FORTRAN.
******* Le retour du programme graphique passe par ce programme avant
******* de retourner au FORTRAN.
******* La structure gloabale du programme est la suivante:
*******
******* --> depart en C : main()
*******			  window_init() initialisation du graphisme
*******			  micromegas()  appel du programme FORTRAN
*******
******* --> boucle en FORTRAN sur ISTEP :
*******				calcul fortran
*******				call predraw --> call draw (trace)
*******	  calcul fortran    <--   predraw    <--   draw
*******
***********************************************************************/

int *n,*iseg,*tabvois,*ivec,*fini,*NMAX,*junct,*nbsegdep,*nsegmax;

{

int id,i,*j,k,*l,*dtabSEG,*dtabPLUS,m,quit=1,xnbube,ynbube,znbube,NbubeCube,base,basem,kplusun;
int inck=0,incksurdeux;
float *x_in,*y_in,*z_in;
float *dx_in,*dy_in,*dz_in;
float *ddx_in,*ddy_in,*ddz_in;
float scale,diophy_scale,un_sur_diophy_scale,t,tx,ty,tz,*y,tnmaxref,UnSurNbube;
float petx,pety,petz,decax,decay,decaz;

/****************************/
/* Equations diophantiennes */
int numDuPlan;
double numDuPlan1,FnumDuPlan;

int ix_test=0,iy_test=0,iz_test=0;
int ix=0,iy=0,iz=0,itest=0,itestold=0;
float posx=0,posy=0,posz=0;
int tposx=0,tposy=0,tposz=0;
int ttposx=0,ttposy=0,ttposz=0;
long long int scaNor;
/* Equations diophantiennes */
/****************************/

switch(nbube)
{
case 0 :
UnSurNbube=1.0;
NbubeCube=1;
break;
default:
UnSurNbube=1.0/nbube;
NbubeCube=nbube*nbube*nbube;
break;
}

/*    INFOS GEOMETRIQUES    */
/****************************/

do /* boucle de fonctionement draw predraw */
{

tnmaxref=1;          /* init value */
scale=1;             /* init value */
tnmax[0]=*NMAX;      /* Lx */
tnmax[1]=*(NMAX+1);  /* Ly */
tnmax[2]=*(NMAX+2);  /* Lz */
Lx=tnmax[0];
Ly=tnmax[1];
Lz=tnmax[2];
tnmax[3]=*(NMAX+3); /* dxy */
tnmax[4]=*(NMAX+4); /* dxz */
tnmax[5]=*(NMAX+5); /* dyx */
tnmax[6]=*(NMAX+6); /* dyz */
tnmax[7]=*(NMAX+7); /* dzx */
tnmax[8]=*(NMAX+8); /* dzy */
dxy=tnmax[3];
dxz=tnmax[4];
dyx=tnmax[5];
dyz=tnmax[6];
dzx=tnmax[7];
dzy=tnmax[8];
tnmax[9]=*(NMAX+9);   /* nbaseredmax */
tnmax[10]=*(NMAX+10); /* ntsg */
tnmax[11]=*(NMAX+11); /* nbase */
tnmax[12]=*(NMAX+12); /* 0.000001/avalue */
enmic = tnmax[12];
for (i=0;i<tnmax[11];i++){tnmax[13+i]=*(NMAX+13+i);} /* Base */
if (tnmax[0]<tnmax[1]) tnmaxref=tnmax[1];else tnmaxref=tnmax[0]; /* max de Lx Ly Lz*/
if (tnmaxref<tnmax[2]) tnmaxref=tnmax[2];
scale=((float)taille_fenetre/tnmaxref);/*on ajuste le trace */
t = (float)(taille_fenetre/2); /* 1- a la taille de la fenetre*/
tx = t*((float)tnmax[0]/tnmaxref);
ty = t*((float)tnmax[1]/tnmaxref);
tz = t*((float)tnmax[2]/tnmaxref);
petx = tx*UnSurNbube;
pety = ty*UnSurNbube;
petz = tz*UnSurNbube;
*nseg= *n+(12*NbubeCube);		       /* 2- a l'echelle du cube en angstroms */
bufnseg = (*n)+(12*NbubeCube);
*tnseg= *nbsegdep;
*nsmax= *nsegmax;
*bufnmax= tnmaxref; /* *NMAX; */
i= *nseg*8;
buf2nmax= tnmaxref; /* *NMAX; */
x_in=kbx_in;
y_in=kby_in;
z_in=kbz_in;

/************************************************************************************************************************/
/*    AFFICHAGE STANDARD    */
/****************************/

if(*nseg>SatNsegmax){fprintf(stderr,"trop de segment (mode standard)...\n");exit(EXIT_FAILURE);}

if(dioph==0) /* pas dioph */
{
for (i=0;i<*n;i++) /* sur les segments */
{
j=iseg+5*i;
k = *(j+3);
id= *(j+4);
l= ivec+(id*3-1);
vl_1=*(l-2);
vl_2=*(l-1);
vl_3=*(l);

/*
tmp_x=*j    ;
tmp_y=*(j+1);
tmp_z=*(j+2);
*/

tmp_x = *j     ;
tmp_y = *(j+1) ;
tmp_z = *(j+2) ;
tmp_x = tmp_x - x_centre;  /* pour preciser le plan a visualier */
tmp_y = tmp_y - y_centre;
tmp_z = tmp_z - z_centre;
ttmp_x=((modulo(tmp_y,Ly))-tmp_y)/Ly*dxy+((modulo(tmp_z,Lz))-tmp_z)/Lz*dxz;
ttmp_y=((modulo(tmp_x,Lx))-tmp_x)/Lx*dyx+((modulo(tmp_z,Lz))-tmp_z)/Lz*dyz;
ttmp_z=((modulo(tmp_x,Lx))-tmp_x)/Lx*dzx+((modulo(tmp_y,Ly))-tmp_y)/Ly*dzy;
ttmp_x=tmp_x+ttmp_x;
ttmp_y=tmp_y+ttmp_y;
ttmp_z=tmp_z+ttmp_z;
tmp_x=modulo(ttmp_x,Lx);
tmp_y=modulo(ttmp_y,Ly);
tmp_z=modulo(ttmp_z,Lz);

/*
tmp_x=*j    ;
tmp_y=*(j+1);
tmp_z=*(j+2);
*/

*(x_in+i)= tmp_x*scale-tx;/* tableau sauvegarde (deux fois plus long : O puis E) */
*(y_in+i)= tmp_y*scale-ty;
*(z_in+i)= tmp_z*scale-tz;
m=i+*nseg;
*(x_in+m)=(tmp_x + vl_1*k)*scale-tx;
*(y_in+m)=(tmp_y + vl_2*k)*scale-ty;
*(z_in+m)=(tmp_z + vl_3*k)*scale-tz;
}
}
/************************************************************************************************************************/
/*   AFFICHAGE DIOPHANTIEN  */
/****************************/
else
{
dx_in=dbx_in;
dy_in=dby_in;
dz_in=dbz_in;
ddx_in=ddbx_in;
ddy_in=ddby_in;
ddz_in=ddbz_in;
dtabSEG=dbtabSEG;
dtabPLUS=dbtabPLUS;
ix = dixmax;iy = diymax;iz = dizmax;

diophyPLUS=0;
diophySEG=0;

for (i=0;i<*n;i++)
{
j=iseg+5*i;
k = *(j+3);
id= *(j+4);
l= ivec+(id*3-1);
vl_1=*(l-2);
vl_2=*(l-1);
vl_3=*(l);

/* RM 2014 le point d'intersection n'est pas forcément obtenu par une translatiuon entiere du vecteur ligne de l'origine */
/* c'est le vrai problème qui explique les clignotement d'avant */
if(abs(vl_1)>abs(vl_2))
{
 diophy_scale=abs(vl_1);
}
else
{
 diophy_scale=abs(vl_2);
}
if(diophy_scale<abs(vl_3))diophy_scale=abs(vl_3);
un_sur_diophy_scale=1./diophy_scale;

scaNor=vl_1*dh+vl_2*dk+vl_3*dl; /* RM2014 produit scalaire vecteur ligne normal plan deplié pour savoir si ou pas dans le plan */

posx = *j     ;
posy = *(j+1) ;
posz = *(j+2) ;
posx =posx- x_centre;  /* pour preciser le plan a visualier */
posy =posy- y_centre;
posz =posz- z_centre;
ttposx =posx;
ttposy =posy;
ttposz =posz;

posx =posx -O1;  /* pour definir proprement les numeros de plan */
posy =posy -O2;
posz =posz -O3;
numDuPlan1=((posx*dh)+(posy*dk)+(posz*dl))/(double)dInd; /* RM2014 Numéro du plan à l'origine */

posx = ttposx + (k * vl_1) ;
posy = ttposy + (k * vl_2) ;
posz = ttposz + (k * vl_3) ;
posx =posx -O1;
posy =posy -O2;
posz =posz -O3;

posx =ttposx;
posy =ttposy;
posz =ttposz;

if(k!=0)  /* Lent mais exactement les points ou le segments perce les polygones */
{	      /* Par contre peut etre qq pbs avec sous reseau */
/* Si besoin est on peut accelerer avec qq produits scalaires pour
remplacer la boucle sur la longueur du segment*/

if(scaNor!=0) /* pas ds le plan (radical) */
{

//fprintf(stderr,"O %d %d %d\n",tmp_x,tmp_y,tmp_z);


itest = 0;
itestold = 0;

idiophy=-1;

/* RM 2014 on applique un facteur correctif d'échelle pour essayer de ne pas passer à côté d'une intersection */
/* RM 2014 (pour accélérer on pourait procéder autrement à l'aide d'un produit scalaire) */
kplusun=(k*diophy_scale)+1;
//fprintf(stderr,"SEG %d L %f Lds %f\n",i,k+1,kplusun);

inck=0;
incksurdeux=0;

for(idioph=0;idioph<kplusun;) /* RM 2013 BOUCLE SUR LA LONGUEUR DU SEGMENT*/
{

ix_test=0;iy_test=0;iz_test=0;

FnumDuPlan=(((posx-O1)*dh)+((posy-O2)*dk)+((posz-O3)*dl))/(double)dInd;
numDuPlan=FnumDuPlan;

if(numDuPlan==FnumDuPlan){ /* RM2014 Si on passe par un plan (et on suppose que c'est le cas seulement sur les points de notre réseau et a priori c'est assez raisonnable)  */

 if(numDuPlan<0||numDuPlan>Plan){

  tposx=((modulo(posy,Ly))-posy)/Ly*dxy+((modulo(posz,Lz))-posz)/Lz*dxz;/* RM2014 on calcul le Point d'impact corrigé si pb de CLP / indic plan */
  tposy=((modulo(posx,Lx))-posx)/Lx*dyx+((modulo(posz,Lz))-posz)/Lz*dyz;
  tposz=((modulo(posx,Lx))-posx)/Lx*dzx+((modulo(posy,Ly))-posy)/Ly*dzy;/* RM2014 posx,y,z st connus calculés avant entrée ds la boucle et actualisé en fin de boucle */

  tposx=posx+tposx;
  tposy=posy+tposy;
  tposz=posz+tposz;

  posx=modulo(tposx,Lx);
  posy=modulo(tposy,Lz);
  posz=modulo(tposz,Lz);

  FnumDuPlan=(((posx-O1)*dh)+((posy-O2)*dk)+((posz-O3)*dl))/(1.*dInd);
  numDuPlan=FnumDuPlan;
 }

 ix_test=diplx[numDuPlan];
 iy_test=diply[numDuPlan];
 iz_test=diplz[numDuPlan];
 if(idiophy!=idioph){itest++;} /* RM2014 idiophy est mis à zéro avant la boucle, idioph est la variable de boucle     */
 idiophy=idioph;               /*        test sert à gérer le mode itératif : pas à pas ou par tronçon entre polygone */
 tmp_x=posx;
 tmp_y=posy;
 tmp_z=posz;

// } SUPP

if(itest != itestold)
{

 ix=ix_test;
 iy=iy_test;
 iz=iz_test;
 tmp_x = tmp_x+Lx*ix+iy*dxy+iz*dxz;
 tmp_y = tmp_y+Ly*iy+ix*dyx+iz*dyz;
 tmp_z = tmp_z+Lz*iz+ix*dzx+iy*dzy;

 /* On trace le seg entre deux polygones */
 if(itest==1){ /* RM2014 On s'occupe du residu à l'origine */
  inck=idiophy;
  incksurdeux=idiophy;
 }
 else if(itest==2){ /* RM2014 puis du premier tronçon complet, on trace ce qui est sous le plan déplié c'est le plus simple */
  inck=idiophy-inck;
  incksurdeux=inck;
  if(incksurdeux==0){incksurdeux=1;inck=2;}
 }
 else if(idioph+inck<kplusun)
{ /* RM2014 puis du dernier tronçon, on ajoute le residu sortant */
 }
 else
{ /* RM2014 puis du dernier tronçon, on ajoute le residu sortant */
  incksurdeux=inck;
  inck=inck+(kplusun-idioph-1);
 }

/*fprintf(stderr,"seg=%d l=%d t=%d O=%d E=%d x!\n",i,idioph,itest,incksurdeux,inck);*/

 dtabPLUS[diophyPLUS]= i;
 tmp_x = tmp_x - (diopZ*incksurdeux*vl_1*un_sur_diophy_scale);
 tmp_y = tmp_y - (diopZ*incksurdeux*vl_2*un_sur_diophy_scale);
 tmp_z = tmp_z - (diopZ*incksurdeux*vl_3*un_sur_diophy_scale);
 *(dx_in+diophyPLUS)= tmp_x*scale-tx;
 *(dy_in+diophyPLUS)= tmp_y*scale-ty;
 *(dz_in+diophyPLUS)= tmp_z*scale-tz;
 tmp_x = tmp_x + (diopZ*inck*vl_1*un_sur_diophy_scale);
 tmp_y = tmp_y + (diopZ*inck*vl_2*un_sur_diophy_scale);
 tmp_z = tmp_z + (diopZ*inck*vl_3*un_sur_diophy_scale);
 *(ddx_in+diophyPLUS)= tmp_x*scale-tx;
 *(ddy_in+diophyPLUS)= tmp_y*scale-ty;
 *(ddz_in+diophyPLUS)= tmp_z*scale-tz;

 diophyPLUS++;
 if(diophyPLUS>(SatNsegmax*10)){fprintf(stderr,"trop de segment (mode pavage)...\n");exit(EXIT_FAILURE);}
  itestold=itest;

}

} //DEPLACE

 if(itest<2) /* RM2014 au début pas à pas */
 {
  idioph++;
  posx = posx+(float)vl_1*un_sur_diophy_scale;
  posy = posy+(float)vl_2*un_sur_diophy_scale;
  posz = posz+(float)vl_3*un_sur_diophy_scale;
 }
 else if(idioph+inck<kplusun) /* RM2014 au milieu tronçon à tronçon */
 {
  idioph=idioph+inck;
  posx = posx+(float)vl_1*inck*un_sur_diophy_scale;
  posy = posy+(float)vl_2*inck*un_sur_diophy_scale;
  posz = posz+(float)vl_3*inck*un_sur_diophy_scale;
 }
 else /* RM2014 à la fin début pas à pas*/
 {
  posx = posx+(float)vl_1*(kplusun-idioph)*un_sur_diophy_scale;
  posy = posy+(float)vl_2*(kplusun-idioph)*un_sur_diophy_scale;
  posz = posz+(float)vl_3*(kplusun-idioph)*un_sur_diophy_scale;
  idioph=idioph+(kplusun-idioph); /* Sortie */
 }


} /* iteration sur la longueur du segment secant */
}
else /* ds le plan ou dans un plan \ */
{
numDuPlan=numDuPlan1; // DESACTIVE
if(numDuPlan1 == numDuPlan)
{
ix=diplx[numDuPlan];
iy=diply[numDuPlan];
iz=diplz[numDuPlan];
tmp_x= posx+Lx*ix+iy*dxy+iz*dxz;
tmp_y= posy+Ly*iy+ix*dyx+iz*dyz;
tmp_z= posz+Lz*iz+ix*dzx+iy*dzy;
dtabSEG[diophySEG]= i;
*(x_in+i)= tmp_x*scale-tx;/* tableau sauvegarde (deux fois plus long : O puis E) */
*(y_in+i)= tmp_y*scale-ty;
*(z_in+i)= tmp_z*scale-tz;
m=i+*nseg;
*(x_in+m)=(tmp_x + vl_1*k)*scale-tx;
*(y_in+m)=(tmp_y + vl_2*k)*scale-ty;
*(z_in+m)=(tmp_z + vl_3*k)*scale-tz;
diophySEG++;
}
else // DESACTIVE
{
numDuPlan=(numDuPlan1+0.5);
ix=diplx[numDuPlan];
iy=diply[numDuPlan];
iz=diplz[numDuPlan];
tmp_x= posx+Lx*ix+iy*dxy+iz*dxz;
tmp_y= posy+Ly*iy+ix*dyx+iz*dyz;
tmp_z= posz+Lz*iz+ix*dzx+iy*dzy;
dtabSEG[diophySEG]= -i;
/* tableau sauvegarde (deux fois plus long : O puis E) */
*(x_in+i)= tmp_x*scale-tx;
*(y_in+i)= tmp_y*scale-ty;
*(z_in+i)= tmp_z*scale-tz;
m=i+*nseg;
*(x_in+m)=(tmp_x + vl_1*k)*scale-tx;
*(y_in+m)=(tmp_y + vl_2*k)*scale-ty;
*(z_in+m)=(tmp_z + vl_3*k)*scale-tz;
diophySEG++;
} /* ds le plan ? */
} /* ds un plan // */
} /* de longueur non nul */
} /* boucle sur les segments */

//fprintf(stderr,"\n dioph dans le plan %d sur %d \n",diophySEG,*n);
//fprintf(stderr,"\n dioph hors du plan %d -> %d \n",(*n-diophySEG),diophyPLUS);
//fprintf(stderr,"\n total %d \n",(diophyPLUS+diophySEG));

} /* dioph */

/************************************************************************************************************************/
/*   AFFICHAGE DE LA BOITE  */
/****************************/
switch(nbube)
{
case 0 :

y=x_in + *n; /* entree des segments formant un    */

*y = -t;     /* les plans */
*(y+1) = t;
*(y+2) = -t;
*(y+3) = t;

*(y+4) = t;
*(y+5) = -t;
*(y+6) = t;
*(y+7) = -t;

*(y+8) = 0;
*(y+9) = -t;
*(y+10) = t;
*(y+11) = -t;

y=y_in + *n;
*y = 0;
*(y+1) = t;
*(y+2) = -t;
*(y+3) = 0;

*(y+4) = 0;
*(y+5) = t;
*(y+6) = -t;
*(y+7) = 0;

*(y+8) = t;
*(y+9) = t;
*(y+10) = -t;
*(y+11) = 0;

y=z_in + *n;
*y = t;
*(y+1) = 0;
*(y+2) = 0;
*(y+3) = -t;

*(y+4) = t;
*(y+5) = 0;
*(y+6) = 0;
*(y+7) = -t;

*(y+8) = t;
*(y+9) = 0;
*(y+10) = 0;
*(y+11) = -t;

y=x_in+ *n + *nseg;
*(y+1) = -t;     /* un cube repere de taille correcte */
*(y+3) = t;
*(y) = -t;
*(y+2) = t;

*(y+5) = t;
*(y+7) = -t;
*(y+4) = t;
*(y+6) = -t;

*(y+8) = 0;
*(y+9) = -t;
*(y+10) = t;
*(y+11) = -t;

y=y_in+ *n + *nseg;
*(y+1) = 0;
*(y+3) = t;
*(y) = -t;
*(y+2) = 0;

*(y+5) = 0;
*(y+7) = t;
*(y+4) = -t;
*(y+6) = 0;

*(y+8) = -t;
*(y+9) = t;
*(y+10) = -t;
*(y+11) = 0;

y=z_in+ *n + *nseg;
*(y+1) = t;
*(y+3) = 0;
*(y) = 0;
*(y+2) = -t;

*(y+5) = t;
*(y+7) = 0;
*(y+4) = 0;
*(y+6) = -t;

*(y+8) = -t;
*(y+9) = 0;
*(y+10) = 0;
*(y+11) = -t;

break;

default:

for(xnbube=1;xnbube<(nbube+1);xnbube++){
for(ynbube=1;ynbube<(nbube+1);ynbube++){
for(znbube=1;znbube<(nbube+1);znbube++){

base = *n +(znbube-1)*12+((ynbube-1)*nbube*12)+((xnbube-1)*nbube*nbube*12);
basem = *nseg + base;

switch(nbube)
{
case 1 :
decax = 0;
decay = 0;
decaz = 0;
break;
default:
decax = -xnbube*2*petx+(nbube+1)*petx;
decay = -ynbube*2*pety+(nbube+1)*pety;
decaz = -znbube*2*petz+(nbube+1)*petz;
break;
} //switch

y=x_in + base; /* entree des segments formant     */
*y = decax+-petx;     /* les cubes */
*(y+1) = decax+petx;
*(y+2) = decax+petx;
*(y+3) = decax+-petx;
*(y+4) = decax+-petx;
*(y+5) = decax+petx;
*(y+6) = decax+petx;
*(y+7) = decax+-petx;
*(y+8) = decax+-petx;
*(y+9) = decax+petx;
*(y+10) = decax+petx;
*(y+11) = decax+-petx;
y=y_in + base;
*y = decay+-pety;
*(y+1) = decay+-pety;
*(y+2) = decay+pety;
*(y+3) = decay+pety;
*(y+4) = decay+-pety;
*(y+5) = decay+-pety;
*(y+6) = decay+pety;
*(y+7) = decay+pety;
*(y+8) = decay+-pety;
*(y+9) = decay+-pety;
*(y+10) = decay+pety;
*(y+11) = decay+pety;
y=z_in + base;
for (i=0;i<8;i++) *(y+i) = decaz+-petz;
for (i=8;i<12;i++) *(y+i) = decaz+petz;
y=x_in+ basem;
*y = decax+petx;
*(y+1) = decax+petx;
*(y+2) = decax+-petx;
*(y+3) = decax+-petx;
*(y+4) = decax+-petx;
*(y+5) = decax+petx;
*(y+6) = decax+petx;
*(y+7) = decax+-petx;
*(y+8) = decax+petx;
*(y+9) = decax+petx;
*(y+10) = decax+-petx;
*(y+11) = decax+-petx;
y=y_in+ basem;
*y = decay+-pety;
*(y+1) = decay+pety;
*(y+2) = decay+pety;
*(y+3) = decay+-pety;
*(y+4) = decay+-pety;
*(y+5) = decay+-pety;
*(y+6) = decay+pety;
*(y+7) = decay+pety;
*(y+8) = decay+-pety;
*(y+9) = decay+pety;
*(y+10) = decay+pety;
*(y+11) = decay+-pety;
y=z_in+ basem;
for (i=0;i<4;i++) *(y+i) = decaz+-petz;
for (i=4;i<12;i++) *(y+i) = decaz+petz;
}}}
break;
} //switch

} while(draw(x_in,y_in,z_in,iseg,ivec,&quit,junct,NMAX,fini));/* appel du programme de  trace */
      if(quit<=0){*fini = quit ;}
}

/**********************************************************************/
