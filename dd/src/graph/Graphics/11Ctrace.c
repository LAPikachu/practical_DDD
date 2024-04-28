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

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include "extercool.h"
#include "external.h"

extern Display *dpy;  /* pour les definitions : voir window_init.c */
extern Window win;
extern GC gc;
extern unsigned long fg;
extern int nseg[1],nsmax[1],tnseg[1],extinct,tnmax[1000];

/****************************/
/* Equations diophantiennes */
int idio;
extern int dioph;
extern int *dstabSEG,*dbtabSEG;
extern int *dstabPLUS,*dbtabPLUS;
extern int diophyPLUS,diophySEG;
extern int *dbxp,*dbyp;
extern int *ddbxp,*ddbyp;
/* Equations diophantiennes */
/****************************/

/**********************************************************************/

void trace( char *tr,
  int co,
  int *iseg,
  int *XO,
  int *YO,
  int *X1,
  int *Y1,
  int segnum,
  int *junct,
  int *nmax)

/***********************************************************************
 *******
 ******* Ce programme image l'objet projete dans la fenetre de visuali-
 ******* sation.
 *******
 ***********************************************************************/

{
  int totoX,tataX,totoY,tataY;
  int i,id1,linbs,touti,tmp;
  char inbs[20];

 /* fprintf(stderr,"TRACE\n"); */

  XSetForeground(dpy,gc,fondcool);/*mise en place de la couleur de */

 if (segnum==-1) {

  if(dioph==0) /* PAS DIOPHANTIEN*/
   {

    for (i=0;i<*nseg-12*nbube*nbube*nbube;i++)
     {

      if ( *(tr+i) ) continue ;/* si tr different de 0->pas de trace*/

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      switch(modecool)
       {
       case 0 : // Souligne la discretisation MARCHE PAS SI NBASERED VARIABLE

        id1= (*(iseg+4+5*i)-(((*(iseg+4+5*i)-1)/tnmax[9])*tnmax[9]))-1;

        XSetForeground(dpy,gc,trcool[id1]);
        if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
        if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
        if(tron[id1]!=0)
         {
          XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
         }

        break;
       case 1 : // Souligne les systemes de glissement
        id1= tnmax[13+*(iseg+4+5*i)-1]-1;

        XSetForeground(dpy,gc,trcool[id1]);
        if((*(junct+i)) < 3 )
         {
          if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
          if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
          if(tron[id1]!=0)
           {
            XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
           }
         }
        else
         {
          XSetForeground(dpy,gc,cooldev);
          if(tron[id1] != 0)
           {
            // XFillArc(dpy,win,gc,XO[i],YO[i],2*(X1[i]-XO[i]),2*(Y1[i]-YO[i]),0 ,360);

            //tmp = 2*sqrt((X1[i]-XO[i])*(X1[i]-XO[i]) +(Y1[i]-YO[i])*(Y1[i]-YO[i])) ;
            tmp = *(iseg+3+5*i)/2;

                        /* fprintf(stderr,"\n i  %d ",i);
                        fprintf(stderr,"\n iseg4  %d ",*(iseg+3+5*i));
                        fprintf(stderr,"\n X1[i]  %d ",X1[i]);
                        fprintf(stderr,"\n XO[i]  %d ",XO[i]);*/

            // XFillRectangle(display,win,gc,XcoinSuperieur_gauche,
            //                 YcoinSuperieur_gauche,largeurEllipse,hauteurEllipse,
            //                 angle_start,angle_end) langle à multiplier par 64;
            XDrawArc(dpy,win,gc,XO[i]-tmp/2,YO[i]-tmp/2,tmp,tmp,0,23040);
            //XFillArc(dpy,win,gc,XO[i]-tmp/20,YO[i]-tmp/20,tmp/10,tmp/10,0,23040);
           }
         }



        break;
       case 2 :
        XSetForeground(dpy,gc,trcool[0]);
        if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
        if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
        XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
        break;
       case 3 :
        XSetForeground(dpy,gc,trcool[0]);
        XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
        break;
       case 4 : /* Souligne la discretisation MARCHE PAS SI NBASERED VARIABLE */

        id1= (*(iseg+4+5*i)-(((*(iseg+4+5*i)-1)/tnmax[9])*tnmax[9]))-1;

        XSetForeground(dpy,gc,trcool[id1]);
        if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
        if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
        if(tron[id1]!=0)
         {
          XDrawLine(dpy,win,gc,XO[i]-1,YO[i]-1,X1[i]-1,Y1[i]-1);
          XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
          XDrawLine(dpy,win,gc,XO[i]+1,YO[i]+1,X1[i]+1,Y1[i]+1);
         }

        break;
       case 5 : // Souligne les systemes de glissement

        id1= tnmax[13+*(iseg+4+5*i)-1]-1;

        XSetForeground(dpy,gc,trcool[id1]);
        if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
        if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
        if(tron[id1]!=0)
         {
          XDrawLine(dpy,win,gc,XO[i]-1,YO[i]-1,X1[i]-1,Y1[i]-1);
          XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
          XDrawLine(dpy,win,gc,XO[i]+1,YO[i]+1,X1[i]+1,Y1[i]+1);
         }

        break;
       }
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(prsegnum==1)
       {
        linbs=sprintf(inbs,"%d",i);
        XDrawString(dpy,win,gc,(XO[i]+X1[i])/2,(YO[i]+Y1[i])/2,inbs,linbs);
       }


     } //for

   }
  else /* DIOPHANTIEN-----------------------------------------------------------------------*/
   {


    for (idio=0;idio<diophySEG;idio++) {

     touti=dstabSEG[idio];
     i=abs(touti);

     if ( *(tr+i) ) continue ;/* si tr different de 0->pas de trace*/

//fprintf(stderr,"\n seg ds plan %d/%d %d modecool %d \n",idio,diophySEG,i,modecool);

     /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
     switch(modecool)
      {
      case 0 : // Souligne la discretisation MARCHE PAS SI NBASERED VARIABLE

       id1= (*(iseg+4+5*i)-(((*(iseg+4+5*i)-1)/tnmax[9])*tnmax[9]))-1;

       XSetForeground(dpy,gc,trcool[id1]);

       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
             if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
        }

       break;
      case 1 : // Souligne les systemes de glissement

       id1= tnmax[13+*(iseg+4+5*i)-1]-1;

       XSetForeground(dpy,gc,trcool[id1]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
       if(touti<0)XSetForeground(dpy,gc,boitecool);
       if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
        }


       break;
      case 2 :
       XSetForeground(dpy,gc,trcool[0]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
              XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
       break;
      case 3 :
       XSetForeground(dpy,gc,trcool[0]);
              XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
       break;
      case 4 : /* Souligne la discretisation MARCHE PAS SI NBASERED VARIABLE */

       id1= (*(iseg+4+5*i)-(((*(iseg+4+5*i)-1)/tnmax[9])*tnmax[9]))-1;

       XSetForeground(dpy,gc,trcool[id1]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
             if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,XO[i]-1,YO[i]-1,X1[i]-1,Y1[i]-1);
         XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
         XDrawLine(dpy,win,gc,XO[i]+1,YO[i]+1,X1[i]+1,Y1[i]+1);
        }

       break;
      case 5 : // Souligne les systemes de glissement

       id1= tnmax[13+*(iseg+4+5*i)-1]-1;

       XSetForeground(dpy,gc,trcool[id1]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
       if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,XO[i]-1,YO[i]-1,X1[i]-1,Y1[i]-1);
         XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
         XDrawLine(dpy,win,gc,XO[i]+1,YO[i]+1,X1[i]+1,Y1[i]+1);
        }

       break;
      }
     /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
     if(prsegnum==1)
      {
       linbs=sprintf(inbs,"%d",i);
       XDrawString(dpy,win,gc,(XO[i]+X1[i])/2,(YO[i]+Y1[i])/2,inbs,linbs);
      }

    } //for

    for (idio=0;idio<diophyPLUS;idio++) {

     i=dstabPLUS[idio];
     totoX=dbxp[idio];
     tataX=ddbxp[idio];
//     if(totoX==tataX)tataX++;// utile si GCLineWidth actif
     totoY=dbyp[idio];
     tataY=ddbyp[idio];
//     if(totoY==tataY)tataY++;// utile si GCLineWidth actif



// fprintf(stderr,"\n seg %d/%d %d modecool %d \n",idio,diophyPLUS,i,modecool);

     /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
     switch(modecool)
      {
      case 0 : // Souligne la discretisation MARCHE PAS SI NBASERED VARIABLE

       id1= (*(iseg+4+5*i)-(((*(iseg+4+5*i)-1)/tnmax[9])*tnmax[9]))-1;

       XSetForeground(dpy,gc,trcool[id1]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
             if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,totoX,totoY,tataX,tataY);
        }

       break;

      case 1 : // Souligne les systemes de glissement

       id1= tnmax[13+*(iseg+4+5*i)-1]-1;
       XSetForeground(dpy,gc,trcool[id1]);

       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
       if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,totoX,totoY,tataX,tataY);
        }


       break;
      case 2 :
       XSetForeground(dpy,gc,trcool[0]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
       XDrawLine(dpy,win,gc,totoX,totoY,tataX,tataY);
       break;
      case 3 :
       XSetForeground(dpy,gc,trcool[0]);

       XDrawLine(dpy,win,gc,totoX,totoY,tataX,tataY);
       break;
      case 4 : /* Souligne la discretisation MARCHE PAS SI NBASERED VARIABLE */

       id1= (*(iseg+4+5*i)-(((*(iseg+4+5*i)-1)/tnmax[9])*tnmax[9]))-1;

       XSetForeground(dpy,gc,trcool[id1]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
             if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,totoX-1,totoY-1,tataX-1,tataY-1);
         XDrawLine(dpy,win,gc,totoX,totoY,tataX,tataY);
         XDrawLine(dpy,win,gc,totoX+1,totoY+1,tataX+1,tataY+1);
        }

       break;
      case 5 : // Souligne les systemes de glissement

       id1= tnmax[13+*(iseg+4+5*i)-1]-1;

       XSetForeground(dpy,gc,trcool[id1]);
       if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
       if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
       if(tron[id1]!=0)
        {
         XDrawLine(dpy,win,gc,totoX-1,totoY-1,tataX-1,tataY-1);
         XDrawLine(dpy,win,gc,totoX,totoY,tataX,tataY);
         XDrawLine(dpy,win,gc,totoX+1,totoY+1,tataX+1,tataY+1);
        }

       break;
      }
     /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
     if(prsegnum==1)
      {
       linbs=sprintf(inbs,"%d",i);
       XDrawString(dpy,win,gc,(totoX+totoX)/2,(tataY+tataY)/2,inbs,linbs);
      }

    } //for

   }

 } //if


 /*                           UN SEUL SEGMENT                          */
 /**********************************************************************/
  if ((segnum!=(-1))&&(segnum-1<=*nseg-12*nbube*nbube*nbube)&&(!(*(tr+segnum-1))))
    {
      /*trace */
      i=segnum-1;
      id1= *(iseg+3+5*i);       /* !!! modif !!! */

   // **************************
   switch(modecool)
    {
    case 0 : // Souligne la discretisation

     id1= *(iseg+4+5*i)-(((*(iseg+4+5*i))/tnmax[9])*tnmax[9])+1;

     XSetForeground(dpy,gc,trcool[id1]);
     if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
     if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
     if(tron[id1]!=0) XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);

     break;
    case 1 : // Souligne les systemes de glissement

         id1= tnmax[13+*(iseg+4+5*i)-1]-1;

     XSetForeground(dpy,gc,trcool[id1]);
     if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
     if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
          if(tron[id1]!=0) XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);

     break;
    case 2 :
     XSetForeground(dpy,gc,trcool[0]);
     if((*(junct+i)) == 1) XSetForeground(dpy,gc,cooljonc);
     if ((*(junct+i))== 2) XSetForeground(dpy,gc,cooldev);
     XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
     break;
    case 3 :
     XSetForeground(dpy,gc,trcool[0]);
     XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
     break;
    case 4 : /* Plus affichage palette */
     XSetForeground(dpy,gc,trcool[0]);
     XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
     break;
    } //switch
   // **************************

      XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
    }





 /*                        trace du cube mauve                         */
 /**********************************************************************/
  XSetForeground(dpy,gc,boitecool);
  i= *tnseg;
  if (i != 0){
  i=i-1;
  XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
 }
  if (extinct)
    {
      for (i= *nseg-12*nbube*nbube*nbube;i<*nseg;i++)
    {
     if ( *(tr+i) ) continue ;
     XDrawLine(dpy,win,gc,XO[i],YO[i],X1[i],Y1[i]);
    }
    }







 /****************************************************/
 /* Diffferentes palettes :                          */
 /* 64k : i de 0 a 500 avec facteur 150 => 0 a 75000 */
 /* blanc 65535 noir 0                               */
 /* 256 : i de 0 a 255 avec facteur 1 => 0 a 255     */
 /* blanc 255 ? (1 plus probablement) noir 0         */
 /****************************************************/
 if(modecool==4){
  for (i= 0;i<500;i++) // Pour avoir la palette
   {
    XSetForeground(dpy,gc,i*1000000);
    XDrawLine(dpy,win,gc,(i-1)*5,5,(i)*5,5);
    XDrawLine(dpy,win,gc,(i-1)*5,6,(i)*5,6);
    XDrawLine(dpy,win,gc,(i-1)*5,7,(i)*5,7);
    if((i/10)*10==i)
     {
      XDrawLine(dpy,win,gc,(i-1)*5,8,(i)*5,8);
      XDrawLine(dpy,win,gc,(i-1)*5,9,(i)*5,9);
      XDrawLine(dpy,win,gc,(i-1)*5,10,(i)*5,10);
     }
   }
 }






 /* Pour avoir la legende couleur segments/systemes */
 /***************************************************/
 if (extinct){
  for (i= 0;i<tnmax[10];i++)
   {
    XSetForeground(dpy,gc,trcool[i]);
    XDrawLine(dpy,win,gc,(i)*10,15,(i+1)*10,15);
    XDrawLine(dpy,win,gc,(i)*10,16,(i+1)*10,16);
    XDrawLine(dpy,win,gc,(i)*10,17,(i+1)*10,17);
    XDrawLine(dpy,win,gc,(i)*10,18,(i+1)*10,18);
    XDrawLine(dpy,win,gc,(i)*10,19,(i+1)*10,19);
    XDrawLine(dpy,win,gc,(i)*10,20,(i+1)*10,20);
   }
  XSetForeground(dpy,gc,cooljonc);
  XDrawLine(dpy,win,gc,0,22,tnmax[10]*10,22);
  XDrawLine(dpy,win,gc,0,23,tnmax[10]*10,23);
  XDrawLine(dpy,win,gc,0,24,tnmax[10]*10,24);
  XSetForeground(dpy,gc,cooldev);
  XDrawLine(dpy,win,gc,0,26,tnmax[10]*10,26);
  XDrawLine(dpy,win,gc,0,27,tnmax[10]*10,27);
  XDrawLine(dpy,win,gc,0,28,tnmax[10]*10,28);
 }
}
