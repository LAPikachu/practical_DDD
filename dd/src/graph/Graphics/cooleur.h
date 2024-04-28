/*
!
!                                       #   #    #   #
!                                       ## ##    ## ##
!                                       # # #    # # #
!                                       #   #    #   #
!                                       #   #    #   #
!                                       #   #ICRO#   #EGAS
!__________________________________________________________________________________________
!
! MM (in the place of MicroMegas) is program of DDD (Discrete Dislocation Dynamics)
! initially developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
! Copyright (C) 1993, 1998, 2000, 2002, 2003, 2004 Benoit Devincre, Ladislas Kubin, Marc Condat,
! Ronan Madec, Ghiath Monnet
!
! This file is part of MM.
!
! MM is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! any later version.
!
! MM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MM; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
!
! For information contact by Email ronan.madec@cea.fr
!__________________________________________________________________________________________
*/

/*********************************************/
/* Lecteur du fichier de configuration de    */
/* l'affichage : couleur écran et postscript */
/* notament.                                 */
/*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "extercool.h"
// #include "varcool.h"

extern  int tnmax[];

int getcool(int *modecool,
            int *fondcool,
	           int *boitecool,
	           int *trcool,
	           double *trR,
	           double *trG,
	           double *trB,
	           int *tron,
	           int *cooljonc,
	           int *cooldev,
	           int *nbube)
{

   FILE *fichier;
   char buf[255];
   int i1,nbtss;

fichier = fopen("../in/couleur.micmeg","rt");
if (fichier == NULL)
{
 printf("\n\n\n ATTENTION !!! Pas de fichier de configuration des couleurs !\n\n\n");
 return(-1);
 exit(1);
}
else
{
 fgets(buf,255,fichier); /*Commentaires*/

 fgets(buf,255,fichier); /*Commentaires*/
 fgets(buf,255,fichier); /*Mode*/
 *modecool=atoi(buf);

 fgets(buf,255,fichier); /*Commentaires*/
 fgets(buf,255,fichier); /*Grille pour méthode des boites ou images ou encore plans pour les carto*/
 *nbube=atoi(buf);

 for(i1=0;i1<100;i1++){
  tron[i1]=0;
  trcool[i1]=0;
  trR[i1]=1.0;
  trG[i1]=1.0;
  trB[i1]=1.0;
 }; /*Initialisation*/

 fgets(buf,255,fichier); /*Commentaires*/
 fgets(buf,255,fichier); /*Fond*/
 *fondcool=atoi(buf);
 trcool[99]=*fondcool;
 fgets(buf,255,fichier);
 trR[99]=atof(buf);
 fgets(buf,255,fichier);
 trG[99]=atof(buf);
 fgets(buf,255,fichier);
 trB[99]=atof(buf);
 fgets(buf,255,fichier); /*Boite*/
 *boitecool=atoi(buf);
 trcool[98]=*boitecool;
 fgets(buf,255,fichier);
 trR[98]=atof(buf);
 fgets(buf,255,fichier);
 trG[98]=atof(buf);
 fgets(buf,255,fichier);
 trB[98]=atof(buf);

 fgets(buf,255,fichier); /*Commentaires*/
 fgets(buf,255,fichier); /*Jonction*/
 *cooljonc=atoi(buf);
 trcool[97]=*cooljonc;
 fgets(buf,255,fichier);
 trR[97]=atof(buf);
 fgets(buf,255,fichier);
 trG[97]=atof(buf);
 fgets(buf,255,fichier);
 trB[97]=atof(buf);
 fgets(buf,255,fichier); /* Glissement devie */
 *cooldev=atoi(buf);
 trcool[96]=*cooldev;
 fgets(buf,255,fichier);
 trR[96]=atof(buf);
 fgets(buf,255,fichier);
 trG[96]=atof(buf);
 fgets(buf,255,fichier);
 trB[96]=atof(buf);

 if(recool==1)
 {
 fgets(buf,255,fichier); /*Commentaires*/
 fgets(buf,255,fichier); /*Nb de type de segments / couleurs dans le fichiers de config*/
 nbtss=atoi(buf);
 if(nbtss<tnmax[10])
  fprintf(stderr,"PROBLEME DE NOMBRE DE COULEURS DANS LE FICHIER Couleur.micmeg-> %d DANS MICMEG->%d \n",nbtss,tnmax[10]);

 /* Couleur Ecran PostScript et Affichage o/N*/
 for(i1=0;i1<nbtss;i1++){
  fgets(buf,255,fichier);
  trcool[i1]=atoi(buf);

  fgets(buf,255,fichier);
  trR[i1]=atof(buf);
  fgets(buf,255,fichier);
  trG[i1]=atof(buf);
  fgets(buf,255,fichier);
  trB[i1]=atof(buf);

  fgets(buf,255,fichier);
  tron[i1]=atoi(buf);
 };
}
fclose(fichier);
fprintf(stderr,"Fichier Couleur.micmeg ok.\n");
}
return 0;
}
