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
#include <string.h>
#include "external.h"

//------------------------------------
//- GENERATION DU POSTSCRIPT COULEUR -
//------------------------------------
void EpsCouleur(FILE *fsortie,int lacouleur)
{
fprintf(fsortie, "%f %f %f setrgbcolor \n",trR[lacouleur],trG[lacouleur],trB[lacouleur]);
}

//----------------------------------
//- GENERATION DU POSTSCRIPT TRAIT -
//----------------------------------
void EpsWidth(FILE *fsortie, int width)
{
fprintf(fsortie, "%d setlinewidth\n", width*100);
}

//----------------------------------
//- GENERATION DU POSTSCRIPT LIGNE -
//----------------------------------
void EpsDrawLine(FILE *fsortie,int lex0,int ley0,int lex1,int ley1)
{
extern unsigned int taille_fenetre;
fprintf(fsortie, "%d\t%d\tmoveto\n", lex0*100, (taille_fenetre-ley0)*100);
fprintf(fsortie, "%d\t%d\tlineto\n", lex1*100, (taille_fenetre-ley1)*100);
fprintf(fsortie, "stroke \n");
}
//----------------------------------
//- GENERATION DU POSTSCRIPT CERCLE -
//----------------------------------
void EpsDrawCircle(FILE *fsortie,int lex0,int ley0,int lex1,int ley1)
{
  extern unsigned int taille_fenetre;
  int  xi,yi,xf,yf,ri;
  xi = lex0*100;
  yi = (taille_fenetre-ley0)*100;
  xf =  lex1*100;
  yf =  (taille_fenetre-ley1)*100;
  ri = sqrt((xf-xi)*(xf-xi) + (yf-yi)*(yf-yi));
// Attention, il y a un problem de format avec les 2 lignes ci dessous
// 	fprintf(fsortie, "100 setlinewidth \n", xi,yi,ri);
//  fprintf(fsortie, "0.700000 0.700000 0.700000 setrgbcolor \n", xi,yi,ri);
fprintf(fsortie, "100 setlinewidth \n");
fprintf(fsortie, "0.700000 0.700000 0.700000 setrgbcolor \n");
  fprintf(fsortie, "%d  %d  %d  0  380  arc\n", xi,yi,ri);
  //fprintf(fsortie, "stroke \n");
  fprintf(fsortie, "fill \n");
}

//--------------------------------------
//--------------------------------------
//- GENERATION DU POSTSCRIPT RECTANGLE -
//--------------------------------------
//--------------------------------------
void EpsDrawRect(FILE *fsortie,int lex0,int ley0,int lex1,int ley1)
{
extern unsigned int taille_fenetre;
fprintf(fsortie, "%d\t%d\t%d\t%d\trectfill\n",lex0*100,
(taille_fenetre-ley0)*100,lex1*100,-ley1*100);
  fprintf(fsortie, "stroke \n");
}

//-----------------------------------
//-----------------------------------
//- GENERATION DU POSTSCRIPT ENTETE -
//-----------------------------------
//-----------------------------------
void debut_postscript(FILE *fsortie,int Bx1,int By1,int Bx2,int By2)
{
extern unsigned int taille_fenetre;
fprintf(fsortie, "%%!PS-Adobe-3.0 EPSF-3.0\n");
fprintf(fsortie, "%%%%BoundingBox: %d %d %d %d\n",Bx1,By1,Bx2,By2);
fprintf(fsortie, "%%%%BeginDocument: tests.eps\n");
fprintf(fsortie, "%%%%Title: graph_micmeg.epsf\n");
fprintf(fsortie, "%%%%DocumentFonts: Helvetica\n");
fprintf(fsortie, "%%%%EndComments\n");
fprintf(fsortie, "0.01 0.01 scale\n");
}

//--------------------------------
//--------------------------------
//- GENERATION DU POSTSCRIPT FIN -
//--------------------------------
//--------------------------------
void fin_postscript(FILE *fsortie)
{
fprintf(fsortie, "showpage\n%%EOF\n%%EndDocument");
}

//-------------------------------------
//-------------------------------------
//- TRACE DANS UN FICHIERS POSTSCRIPT -
//-------------------------------------
//-------------------------------------

void traceps(	char *tr,
    int co,
    int *iseg,
    int *XO,
    int *YO,
    int *X1,
    int *Y1,
    int segnum,
    int *junct,
    int *nmax)

{
  int i,id1,Bx1=0,By1=0,Bx2,By2;
  FILE *fichier;
  char nom[22];
  extern int makefilm,numfilm;
  extern unsigned int taille_fenetre;
  int totoX,tataX,totoY,tataY;

/****************************/
/* Equations diophantiennes */
int idio;
extern  int dioph;
extern  int diophyPLUS,diophySEG;
extern  int *dbtabSEG;
extern  int *dbtabPLUS;
extern  int *dbxp,*dbyp;
extern  int *ddbxp,*ddbyp;
extern  int tnmax[];
/* Equations diophantiennes */
/****************************/

/* OUVERTURE DU FICHIER ET ENTETE */
/**********************************/
sprintf(nom, "../out/imag_%4d.epsf", numfilm);
if(nom[13]==' ')nom[13]='0';
if(nom[14]==' ')nom[14]='0';
if(nom[15]==' ')nom[15]='0';
fprintf(stderr,"\n %s \n",nom);
numfilm++;
if(numfilm==10000){numfilm=0;fprintf(stderr,"\n MAKEFILM NUMFILM -> %d %d \n\n",makefilm,numfilm);}
fichier = fopen(nom,"w");
Bx2=taille_fenetre;
By2=taille_fenetre;
debut_postscript(fichier,Bx1,By1,Bx2,By2);
EpsWidth(fichier,1);
EpsCouleur(fichier,99);
EpsDrawRect(fichier,0,0,taille_fenetre,taille_fenetre);

if(dioph==0){
/* VERSION STANDARD*/

/* TRACE DES SEGMENTS */
/**********************/
for (i=0;i<*nseg-12*nbube*nbube*nbube;i++) {

if ( *(tr+i) ) continue ;/* si tr different de 0->pas de trace */

switch(modecool)
{
  case 0 :	/* Souligne la discretisation */
          id1= *(iseg+4+5*i)-(((*(iseg+4+5*i))/tnmax[9])*tnmax[9])+1;
    if((*(junct+i))== 0)EpsCouleur(fichier,id1);
    if((*(junct+i))== 1)EpsCouleur(fichier,97);
          if((*(junct+i))== 2)EpsCouleur(fichier,96);
            if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
  break;
case 1 :	/* Souligne les systemes de glissement */
  id1= tnmax[13+*(iseg+4+5*i)-1]-1;
  if((*(junct+i)) < 3)
    {
      if((*(junct+i)) == 0) EpsCouleur(fichier,id1);
      if((*(junct+i)) == 1)
        {EpsCouleur(fichier,97);
          EpsWidth(fichier,3);}
      if(tron[id1]!=0) EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
      EpsWidth(fichier,1);
    }
  else
    {
      EpsCouleur(fichier,97);
      EpsWidth(fichier,3);
      if(tron[id1]!=0) EpsDrawCircle(fichier,XO[i],YO[i],X1[i],Y1[i]);
      EpsWidth(fichier,1);
    }

  break;
  case 2 :	/* Une couleur mais on distingue qd meme GD et JONC */
          id1= 0;
                if((*(junct+i))== 0)EpsCouleur(fichier,id1);
    if((*(junct+i))== 1)EpsCouleur(fichier,97);
          if((*(junct+i))== 2)EpsCouleur(fichier,96);
                if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
break;
case 3 :	/* Une couleur */
                id1= 0;
                EpsCouleur(fichier,id1);
                if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
  break;

  } /*switch*/
} /*for*/

}
else
{
/* VERSION PAVAGE*/
/* _________________________________________________________________________________________________________________ */
EpsWidth(fichier,2);
for (idio=0;idio<diophySEG;idio++) {

i=abs(dstabSEG[idio]);

if ( *(tr+i) ) continue ;/* si tr different de 0->pas de trace*/
switch(modecool)
{
  case 0 :	/* Souligne la discretisation */
          id1= *(iseg+4+5*i)-(((*(iseg+4+5*i))/tnmax[9])*tnmax[9])+1;
    if((*(junct+i))== *nsmax)EpsCouleur(fichier,id1);
    if((*(junct+i))!= *nsmax)EpsCouleur(fichier,97);
          if((*(junct+i))== 0)EpsCouleur(fichier,96);
            if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);

  break;
  case 1 :	/* Souligne les systemes de glissement */

          id1= tnmax[13+*(iseg+4+5*i)-1]-1;
    if((*(junct+i))== *nsmax)EpsCouleur(fichier,id1);
    if((*(junct+i))!= *nsmax)EpsCouleur(fichier,97);
          if((*(junct+i))== 0)EpsCouleur(fichier,96);
    if((*(junct+i))!= 0 && (*(junct+i))!= 1)EpsWidth(fichier,3);
            if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
    if((*(junct+i))!= 0 && (*(junct+i))!= 1)EpsWidth(fichier,2);

  break;
  case 2 :	/* Une couleur mais on distingue qd meme GD et JONC */
          id1= 0;
    if((*(junct+i))== *nsmax)EpsCouleur(fichier,id1);
    if((*(junct+i))!= *nsmax)EpsCouleur(fichier,97);
          if((*(junct+i))== 0)EpsCouleur(fichier,96);
                if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);

break;
case 3 :	/* Une couleur */
                id1= 0;
                EpsCouleur(fichier,id1);
                if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
  break;   }
} //for
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
EpsWidth(fichier,2);

for (idio=0;idio<diophyPLUS;idio++) {

i=dstabPLUS[idio];
totoX=dbxp[idio];
tataX=ddbxp[idio];
totoY=dbyp[idio];
tataY=ddbyp[idio];
switch(modecool)
{
  case 0 :	/* Souligne la discretisation */
          id1= *(iseg+4+5*i)-(((*(iseg+4+5*i))/tnmax[9])*tnmax[9])+1;
    if((*(junct+i))== *nsmax)EpsCouleur(fichier,id1);
    if((*(junct+i))!= *nsmax)EpsCouleur(fichier,97);
          if((*(junct+i))== 0)EpsCouleur(fichier,96);
            if(tron[id1]!=0)EpsDrawLine(fichier,totoX,totoY,tataX,tataY);
  break;
  case 1 :	/* Souligne les systemes de glissement */

          id1= tnmax[13+*(iseg+4+5*i)-1]-1;
    if((*(junct+i))== *nsmax)EpsCouleur(fichier,id1);
    if((*(junct+i))!= *nsmax)EpsCouleur(fichier,97);
          if((*(junct+i))== 0)EpsCouleur(fichier,96);
            if(tron[id1]!=0)EpsDrawLine(fichier,totoX,totoY,tataX,tataY);
  break;
  case 2 :	/* Une couleur mais on distingue qd meme GD et JONC */
          id1= 0;
    if((*(junct+i))== *nsmax)EpsCouleur(fichier,id1);
    if((*(junct+i))!= *nsmax)EpsCouleur(fichier,97);
          if((*(junct+i))== 0)EpsCouleur(fichier,96);
                if(tron[id1]!=0)EpsDrawLine(fichier,totoX,totoY,tataX,tataY);
break;
case 3 :	/* Une couleur */
                id1= 0;
                EpsCouleur(fichier,id1);
                if(tron[id1]!=0)EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
  break;   }
} //for


/* _________________________________________________________________________________________________________________ */
}

/* TRACE DU OU DES CUBES OU DES PLANS */
/**************************************/
  EpsWidth(fichier,1);
  EpsCouleur(fichier,98);
  i= *tnseg;
  if (i != 0){
      i=i-1;
      EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
        }
  if (extinct)
    {
      for (i= *nseg-12*nbube*nbube*nbube;i<*nseg;i++)
  {
    if ( *(tr+i) ) continue ;
    EpsDrawLine(fichier,XO[i],YO[i],X1[i],Y1[i]);
  }
    }
  EpsWidth(fichier,1);

/* COULEUR DES SYSTEMES, SEGMENTS ENTRE PLANS ET JONCTION */
/**********************************************************/
if (extinct){
for (i= 0;i<tnmax[10];i++) // Pour avoir les systemes
  {
    EpsCouleur(fichier,i);
          EpsDrawRect(fichier,(i)*10,1,10,5);
  }

  EpsCouleur(fichier,97);
        EpsDrawRect(fichier,0,7,tnmax[10]*10,2);

  EpsCouleur(fichier,96);
        EpsDrawRect(fichier,0,10,tnmax[10]*10,2);
}

/* THIS IS THE END */
/*******************/
fin_postscript(fichier);
fclose(fichier);

}