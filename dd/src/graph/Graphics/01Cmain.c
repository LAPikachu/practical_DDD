/**
* mM (for microMegas) is an open source program of DD (Dislocation Dynamics) simulation
* originally developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
* Copyright (C) 1993, 1996, 2000, 2002, 2004   Benoit Devincre, Ladislas Kubin, Marc Condat,
* Christophe Lemarchand, Ronan Madec, Sebastien Groh, Ghiath Monnet, Julien Durinck,
* Phillipe Carrez, Christophe de Sansal, Mohamed-Gazi Tagorti.
*
* This file is part of mM.
*
* mM is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* any later version.
*
* mM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with MM; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
* - - -
* For more information see http://zig.onera.fr/mm_home_page/
* - - -
*/

/*********************************************************************//**
*
* Ceci est le corps principal du programme.Apres initialisation de
* la fenetre graphique,il donne la main au programme FORTRAN. On
* ne revient plus au programme principal sauf si il y a achevement
* du programme FORTRAN qui signifie egalement la fin de l'ensemble
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "varcool.h"
#include "external.h"

#define InuseGX 1000000
#define InuseImGX 0
#define CnuseImGX 'G'

extern int micromegas_();

/* 2 TABLEAUX : 1 tab d'affichage et 1 tab de sauvegarde */

static float *sx_in,*sy_in,*sz_in;	/* tab affichage */
static float *ksx_in,*ksy_in,*ksz_in;	/* tab sauvegarde */
static float *sx,*sy;
static char *str;
static int *sX0,*sY0,*sX1,*sY1;

float *kbx_in,*kby_in,*kbz_in;
float *bx_in,*by_in,*bz_in;
float *bx,*by;
char *btr;
int *bX0,*bY0,*bX1,*bY1;

/****************************/
/* Equations diophantiennes */
static float *dsx_in,*dsy_in,*dsz_in;
static float *dbx_in,*dby_in,*dbz_in;

static int *dsxp,*dsyp;
static int *dbxp,*dbyp;

static float *ddsx_in,*ddsy_in,*ddsz_in;
static float *ddbx_in,*ddby_in,*ddbz_in;

static int *ddsxp,*ddsyp;
static int *ddbxp,*ddbyp;

static int *dstabSEG,*dbtabSEG;
static int *dstabPLUS,*dbtabPLUS;
/* Equations diophantiennes */
/****************************/

int	argcc;
char	**argvv;
static int SatNsegmax;

int main(int argc, char **argv, char **arge)

{
static int ForNsegmax,NimagesX,NimagesY,NimagesZ;

/* RM Le 21/09/05 */
/* LIEN AVEC LA VERSION MIXTE DE MICROMEGAS */
ForNsegmax=InuseGX;
NimagesX=InuseImGX;
NimagesY=InuseImGX;
NimagesZ=InuseImGX;
/* RM Le 21/09/05 */
/* FIN DU LIEN AVEC LA VERSION MIXTE DE MICROMEGAS */

SatNsegmax=((1+2*NimagesX)*(1+2*NimagesY)*(1+2*NimagesZ))*ForNsegmax;
/* fprintf(stderr,"allocation C -> %i %i %i %i %i\n",NimagesX,NimagesY,NimagesZ,ForNsegmax,SatNsegmax);
fprintf(stderr,"[NSEGMAX]: %i segments.\n",InuseGX);
fprintf(stderr,"[Nb Images]: %i %i %i .\n", NimagesX, NimagesY, NimagesZ); */
fprintf(stderr,"memory allocation -> %i standard \n",SatNsegmax);
fprintf(stderr,"                  -> %i pavage \n",ForNsegmax*10);

sx_in=(float *)malloc(SatNsegmax*8);
if(sx_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
sy_in=(float *)malloc(SatNsegmax*8);
if(sy_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
sz_in=(float *)malloc(SatNsegmax*8);
if(sz_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}

ksx_in=(float *)malloc(SatNsegmax*8);
if(ksx_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ksy_in=(float *)malloc(SatNsegmax*8);
if(ksy_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ksz_in=(float *)malloc(SatNsegmax*8);
if(ksz_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}

/****************************/
/* Equations diophantiennes */

/* ordre de grandeur : une boite de 10 microns
1000 plans et des segments secants de 1 micron
donc 100 fois plus de segment si les segments
etaient perpendiculaires aux plans, on va mettre
10 ce qui est un peu juste mais bon...*/

dsx_in=(float *)malloc(SatNsegmax*8*10);  /* Segments secants origine 3D */
if(dsx_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dsy_in=(float *)malloc(SatNsegmax*8*10);
if(dsy_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dsz_in=(float *)malloc(SatNsegmax*8*10);
if(dsz_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dbx_in= dsx_in;
dby_in= dsy_in;
dbz_in= dsz_in;

dsxp=(int *)malloc(SatNsegmax*4*10); /* Segments secants origine 2D */
if(dsxp==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dsyp=(int *)malloc(SatNsegmax*4*10);
if(dsyp==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dbxp= dsxp;
dbyp= dsyp;

ddsx_in=(float *)malloc(SatNsegmax*8*10);  /* Segments secants extremite 3D */
if(ddsx_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ddsy_in=(float *)malloc(SatNsegmax*8*10);
if(ddsy_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ddsz_in=(float *)malloc(SatNsegmax*8*10);
if(ddsz_in==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ddbx_in= ddsx_in;
ddby_in= ddsy_in;
ddbz_in= ddsz_in;

ddsxp=(int *)malloc(SatNsegmax*4*10); /* Segments secants extremite 2D */
if(ddsxp==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ddsyp=(int *)malloc(SatNsegmax*4*10);
if(ddsyp==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
ddbxp= ddsxp;
ddbyp= ddsyp;

dstabSEG=(int *)malloc(SatNsegmax*4); /* Indice des segments dans le plan */
if(dstabSEG==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dbtabSEG=dstabSEG;

dstabPLUS=(int *)malloc(SatNsegmax*4*10); /* Indice des segments secants pour
                                          le code couleur */
if(dstabPLUS==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
dbtabPLUS=dstabPLUS;

/* Equations diophantiennes */
/****************************/

sx=(float *)malloc(SatNsegmax*8);
if(sx==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
sy=(float *)malloc(SatNsegmax*8);
if(sy==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}

str=(char *)malloc(SatNsegmax);

sX0=(int *)malloc(SatNsegmax*4);
if(sX0==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
sY0=(int *)malloc(SatNsegmax*4);
if(sY0==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
sX1=(int *)malloc(SatNsegmax*4);
if(sX1==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}
sY1=(int *)malloc(SatNsegmax*4);
if(sY1==NULL){fprintf(stderr,"l'allocation a echouee...\n");exit(EXIT_FAILURE);}

bx_in= sx_in;
by_in= sy_in;
bz_in= sz_in;
kbx_in= ksx_in;
kby_in= ksy_in;
kbz_in= ksz_in;

bx=    sx;
by=    sy;
btr=   str;
bX0=   sX0;
bY0=   sY0;
bX1=   sX1;
bY1=   sY1;

/* fprintf(stderr,"Appel X\n");*/
window_init(argcc,argvv);
/* fprintf(stderr,"Appel FORTRAN\n");*/

#ifdef MAC
micromegas();
#else
micromegas_();
#endif

free(sx_in);
free(sy_in);
free(sz_in);
free(ksx_in);
free(ksy_in);
free(ksz_in);

return 0;
}
/**********************************************************************/


