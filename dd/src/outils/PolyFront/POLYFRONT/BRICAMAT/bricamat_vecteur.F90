! MobiDiC (in the place of Mobile Dislocation Colony) is program of DD (Dislocation Dynamics)
! developed at the 'CEA, DAM, DIF' CEA, FRANCE
! Copyright (C) 2003, 2006, 2007, 2008 Ronan MADEC, Hervé TROADEC, Vincent LIARD,
! Aymeric CANTON, Anthony DYAN
!
! MobiDiC is initially based on mM (in the place of microMEGAS)
! developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
! Copyright (C) 1993, 1998, 2000 Benoit DEVINCRE, Ladislas KUBIN, Marc CONDAT,
! Ronan MADEC, Ghiath MONNET
!
! This file is part of MobiDiC
!
! MobiDiC is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! any later version.
!
! MobiDiC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MobiDiC; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
!
! For information contact by Email Ronan MADEC : ronan.madec@cea.fr 

!**********************************************************
!* valeur absolue du produit scalaire de vecteurs entiers *
!* divise par la norme des vecteurs                       *
!**********************************************************
function absprodscanorm(v1,v2)

implicit none
real(kind=DPR) :: v1(3),v2(3),absprodscanorm

absprodscanorm=abs((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/&
(sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))*&
sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))))

end function absprodscanorm

function prodscanorm(v1,v2)

implicit none
real(kind=DPR) :: v1(3),v2(3),prodscanorm

prodscanorm=(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/&
(sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))*&
sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)))

end function prodscanorm

!**********************************************************
!* valeur absolue du produit scalaire de vecteurs entiers *
!**********************************************************
function iabsprodsca(v1,v2)

implicit none
integer(kind=DPI) :: v1(3),v2(3),iabsprodsca

iabsprodsca=abs(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))

end function iabsprodsca

!###################################################
!# Calcul du carre de la norme d'un vecteur entier #
!###################################################
function inorivect(i)

implicit none

integer(kind=DPI):: i(3),inorivect

inorivect=i(1)*i(1)+i(2)*i(2)+i(3)*i(3)
return

end function inorivect

!###################################################
!# Calcul du carre de la norme d'un vecteur entier #
!###################################################
function carnorvect(i)

implicit none

real(kind=DPR):: i(3),carnorvect

carnorvect=i(1)*i(1)+i(2)*i(2)+i(3)*i(3)
return

end function carnorvect

!###################################################
!# Calcul du carre de la norme d'un vecteur entier #
!###################################################
function lnorivect(i)

implicit none

integer(kind=DPI) :: i(3)
integer(kind=DPI) :: lnorivect,r(3)

r=i
lnorivect=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)

end function lnorivect

!##########################################
!# Calcul de la norme d'un vecteur entier #
!##########################################
function norivect(i)

implicit none

integer(kind=DPI) :: i(3)
real(kind=DPR) :: r(3),norivect

r(1:3)=real(i(1:3),DPR)
norivect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
return

end function norivect

!###########################################################
!# Calcul du vecteur normalise associe a un vecteur entier #
!###########################################################
function normaivec(i)

implicit none

integer(DPI) :: i(3)
real(kind=DPR) :: normaivec(3),r(3),deno

deno=1./norivect(i)
r=dble(i)
normaivec(:)=r(:)*deno

return
end function normaivec

!########################################
!# Calcul de la norme d'un vecteur reel #
!########################################
function norvect(r)

implicit none

real(kind=DPR) :: norvect,r(3)

norvect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

return

end function norvect

!#########################################################
!# Calcul du vecteur normalise associe a un vecteur reel #
!#########################################################

function normavec(r)

implicit none

real(kind=DPR) :: normavec(3),r(3),deno

deno=1./norvect(r)            !*** division plus lourde que multiplication
normavec(:)=r(:)*deno

return
end function normavec

!##################################################
!# Produit scalaire entre deux vecteurs reels a.b #
!##################################################
function prodsca(a,b)

implicit none

real(kind=DPR) :: prodsca,a(3),b(3)

prodsca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

return
end function prodsca

!####################################################
!# Produit scalaire entre deux vecteurs entiers a.b #
!####################################################
function iprodsca(a,b)

implicit none

integer(DPI) :: iprodsca
integer(DPI)  :: a(3),b(3)

iprodsca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

return
end function iprodsca

!####################################################
!# Produit scalaire entre deux vecteurs entiers a.b #
!# Passage en DPIV                                  #
!####################################################
function viprodsca(a,b)

implicit none

integer(DPIv) :: viprodsca
integer(DPI)  :: a(3),b(3)
integer(DPIv)  :: aa(3),bb(3)

aa(:)=a(:)
bb(:)=b(:)

viprodsca=aa(1)*bb(1)+aa(2)*bb(2)+aa(3)*bb(3)

return
end function viprodsca

!#######################################################################
!# Produit scalaire entre deux vecteurs entiers a.b le résultats étant #
!# retourné dans un réel.                                              #
!#######################################################################
function prodisca(a,b)

implicit none
integer(kind=DPI) :: a(3),b(3)
real(kind=DPR) :: prodisca

prodisca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

return
end function prodisca

!######################################################################
!# Calcul du produit vectoriel a^b  (a et b sont des vecteurs entier) #
!######################################################################
function prodivec(a,b)

implicit none

integer(DPI) :: prodivec(3),a(3),b(3)

prodivec(1) = a(2)*b(3)-a(3)*b(2)
prodivec(2) = a(3)*b(1)-a(1)*b(3)
prodivec(3) = a(1)*b(2)-a(2)*b(1)

return
end function prodivec

!###################################
!# Calcul du produit vectoriel a^b #
!###################################
function prodvec(a,b)

implicit none

real(kind=DPR) :: prodvec(3),a(3),b(3)

prodvec(1) = a(2)*b(3)-a(3)*b(2)
prodvec(2) = a(3)*b(1)-a(1)*b(3)
prodvec(3) = a(1)*b(2)-a(2)*b(1)

return
end function prodvec

!*************************************
! Produit vectoriel et normalisation *
!*************************************
function nprodvec(a,b)

implicit none

real(kind=DPR) :: nprodvec(3),a(3),b(3),ndenpv

nprodvec(1) = a(2)*b(3)-a(3)*b(2)
nprodvec(2) = a(3)*b(1)-a(1)*b(3)
nprodvec(3) = a(1)*b(2)-a(2)*b(1)

ndenpv=dsqrt(nprodvec(1)**2+nprodvec(2)**2+nprodvec(3)**2)

if(ndenpv.ne.0.d0) then
nprodvec(1:3)=nprodvec(1:3)/ndenpv
else
endif

return
end function nprodvec

!##############################
!# Vecteurs perpendiculaire ? #
!##############################
function perpendi(a,b)

implicit none
logical :: perpendi

#ifdef LSD
real(kind=DPR) :: a(3),b(3)
#endif
#ifndef LSD
integer(kind=DPI) :: a(3),b(3)
#endif

#ifdef LSD
if (valeur_LSD_nulle(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))) then ! RMZa : c'est un produit scalaire
! il faudrait une traduction correct en norme i.e. projeté perpendiculairement a sur b < à precision et inversement
#endif
#ifndef LSD
if((a(1)*b(1)+a(2)*b(2)+a(3)*b(3)).eq.0)then
#endif
perpendi=.true.
else
perpendi=.false.
endif

return
end function perpendi


!########################
!# Vecteurs confondus ? #
!########################
function confondu(a,b)

implicit none
logical :: confondu
integer(DPI) :: prodivec(3),a(3),b(3)

prodivec(1) = a(2)*b(3)-a(3)*b(2)
prodivec(2) = a(3)*b(1)-a(1)*b(3)
prodivec(3) = a(1)*b(2)-a(2)*b(1)

if(prodivec(1).eq.0.and.prodivec(2).eq.0.and.prodivec(3).eq.0)then
confondu=.true.
else
confondu=.false.
endif

return
end function confondu
