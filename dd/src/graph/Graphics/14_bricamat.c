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
#include "external.h"

int modulo(int a, int b)

{
int aa,bb,resul;

resul = -999999;

if(b>0){
  if(a>=0){
    resul = a%b;
  }
  else{
    aa=abs(a);
    resul = b-aa%b;
  }
}

if(b<0){
  bb=abs(b);
  if(a>=0){
    resul = b+a%bb;
  }
  else{
    aa=abs(a);
    resul = -aa%bb;
  }
}

return(resul);
}

/*-------------------------------------
/ -              MAX                 -
/ ------------------------------------- */

int mymax(int a, int b)

{
int aa,bb;

aa=abs(a);
bb=abs(b);

if(abs(a)>=abs(b)){
return(aa);
}
else
{
return(bb);
}
}

/*-------------------------------------
/ -              PGCD                 -
/ ------------------------------------- */

int pgdc( int aa, int bb)

{
int a,b,c,resul;

if(bb>aa){a=bb;b=aa;}else{a=aa;b=bb;}
if(a*b==0)
{
resul=mymax(a,b);
}
else
{
for(b=b;b!=0;b=c)
{
c=a-((a/b)*b);
a=b;
}
resul=abs(a);
}
return(resul);
}

/*-------------------------------------
/ -              PGCD                 -
/ ------------------------------------- */

int pgdc3(int a, int b, int c)

{
int resul;
if(a==0)
{
resul=pgdc(b,c);
}
else
{
if(b==0)
{
resul=pgdc(a,c);
}
else
{
if(c==0)
{
resul=pgdc(a,b);
}
else
{
resul=pgdc(a,pgdc(b,c));
}
}
}

return(resul);
}
