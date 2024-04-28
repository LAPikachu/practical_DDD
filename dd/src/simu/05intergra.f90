
!===================================================================================================
!========================    DEBUT    MODULE  "INTERGRA"   =========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains the procedures needed to visualize the dislocation trajectory during the computations.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module INTERGRA

use VARGLOB

implicit none

contains

!##########################################################################
!# Interface Graphique                                                    #
!########################################################## 27/10/00 ######

subroutine MODEGRAF(nbsegdepin)

implicit none

integer (kind=DPI) :: nbsegdepin
integer (4) :: I,IJONC(gNSEGMAX),ISEG(5,gNSEGMAX),TABVOIS(2,gNSEGMAX),nbsegdep,tmp
integer (4) :: NsegmCP
integer (4) :: NparCP

! RM le 22/09/05
!
! ligne avant interface graphique avec dÈveloppement diophantien :
!
!integer (4) :: GIMIX,GIMIY,GIMIZ,GISEGM,GIMODUR(4),kkcp, nsegmaxcp
!
! on se sert de GIMODUR pour transferer les infos dont on a besoin,
! je fais au plus simple pour une base de maximum 200 vecteurs
!
integer (4) :: GIMIX,GIMIY,GIMIZ,GISEGM,GIMODUR(213),kkcp, nsegmaxcp

!*** Construction des variable d'echange utils

nbsegdep  = int(nbsegdepin)
BVECLINCP = int(BVECLIN)
KKCP      = int(KK)
NSEGMAXCP = int(NSEGMAX)
NsegmCP   = int(Nsegm)
NparCP    = int(Npar)

do i = 1, NsegmCP

    ISEG(1,i) = int(SEG(i)%O(1))
    if (modulo((iseg(1,i)+SEG(i)%NORME/2*bveclin(1,SEG(i)%veclin)),modur(1))                  &
    .ne.(iseg(1,i)+SEG(i)%NORME/2*bveclin(1,SEG(i)%veclin)))                                  &
    ISEG(1,i) = int(modulo((iseg(1,i)+SEG(i)%NORME/2*bveclin(1,SEG(i)%veclin)),modur(1))      &
    -SEG(i)%NORME/2*bveclin(1,SEG(i)%veclin))
    ISEG(2,i) = int(SEG(i)%O(2))
    if (modulo((iseg(2,i)+SEG(i)%NORME/2*bveclin(2,SEG(i)%veclin)),modur(2))&
    .ne.(iseg(2,i)+SEG(i)%NORME/2*bveclin(2,SEG(i)%veclin)))                                  &
    ISEG(2,i) = int(modulo((iseg(2,i)+SEG(i)%NORME/2*bveclin(2,SEG(i)%veclin)),modur(2))      &
    -SEG(i)%NORME/2*bveclin(2,SEG(i)%veclin))
    ISEG(3,i) = int(SEG(i)%O(3))
    if (modulo((iseg(3,i)+SEG(i)%NORME/2*bveclin(3,SEG(i)%veclin)),modur(3))                  &
    .ne.(iseg(3,i)+SEG(i)%NORME/2*bveclin(3,SEG(i)%veclin)))                                  &
    ISEG(3,i) = int(modulo((iseg(3,i)+SEG(i)%NORME/2*bveclin(3,SEG(i)%veclin)),modur(3))      &
    -SEG(i)%NORME/2*bveclin(3,SEG(i)%veclin))
    ISEG(4,i) = int(SEG(i)%NORME)
    ISEG(5,i) = int(SEG(i)%VECLIN)
    TABVOIS(1,i) = int(SEG(i)%Voiso)
    TABVOIS(2,i) = int(SEG(i)%Voise)

    if(seg(i)%jonc) then
       IJONC(I) = Iun
    elseif (seg(i)%gd > izero) then
       IJONC(I) = Ideux
    else
       IJONC(I) = Izero
    endif
enddo

do I = 1, NparCP
    tmp               = NsegmCP + I
    IJONC(tmp)        = itrois
    ISEG(1:3,tmp)     = int(par(i)%C(1:3))
    ISEG(4,tmp)       = nint(par(i)%R/normlin(IUN))
    ISEG(5,tmp)       = IUN
    TABVOIS(1,tmp)    = nsegmax
    TABVOIS(2,tmp)    = nsegmax
enddo

GISEGM = NsegmCP + NparCP

! GIMODUR definitions
! This table is used to transer data between the c and Fort codes

! Dimension of the simulation cell to visualise
GIMODUR(1:3) = int(MODUR(1:3)*(2*affnboimage+1))

! The shifting informations, for the moment set to zero
! POUR LE MOMENT LES DECALAGES SONT A ZERO, AU PASSAGE IL SERAIT
! PEUT ETRE BON D'UTILISER SIX PETITS SCALAIRES INITIALISES UNE BONNE
! FOIS POUR TOUTE AU LIEU DE FAIRE DES SUMS A CHAQUE OPERATION DE CLP...
GIMODUR(4:9) =0 !shiftxy

! In the following a Base of 8 vectors per slip systems is assumed
! but the graphical interface is more generale than this and works
! with more or less sofisticated simulation base of discretization
GIMODUR(10) = 8                 ! nb of vecteurs per slip systems
GIMODUR(11) = int(NBASE/8)      ! nb of slip systems taken into account
GIMODUR(12) = int(nbase)             ! total number of vectors

! The space scalling factor
GIMODUR(13) = int(0.000001/avalue)

do I = 1, int(NBASE)
 GIMODUR(13+i)=(i-1)/8+1 ! numero du systeme
enddo

!*** Le 06/06/01 : EFFET DES IMAGES ******
if (affnboimage.ne.0) then

   if(npar > izero) stop " Replica are not possible with particules"

!*** Boucles sur les images *****************
   gimix=0
   gimiy=0
   gimiz=0
   do while(gimix.ne.(2*affnboimage).or.gimiy.ne.(2*affnboimage).or.gimiz.ne.(2*affnboimage))
       if(gimix.ne.(2*affnboimage)) then
          gimix=gimix+1
       else
          if(gimiy.ne.(2*affnboimage)) then
             gimix=0
             gimiy=gimiy+1
          else
             if(gimiz.ne.(2*affnboimage)) then
                gimix=0
                gimiy=0
                gimiz=gimiz+1
             endif
          endif
       endif

       IJONC(gISEGM+1:gISEGM+NSEGM)     = int(seg(1:nsegm)%ijonc)
       ISEG(1,gISEGM+1:gISEGM+NSEGM)    = int(ISEG(1,1:NSEGM)+gimix*modur(1))
       ISEG(2,gISEGM+1:gISEGM+NSEGM)    = int(ISEG(2,1:NSEGM)+gimiy*modur(2))
       ISEG(3,gISEGM+1:gISEGM+NSEGM)    = int(ISEG(3,1:NSEGM)+gimiz*modur(3))
       ISEG(4,gISEGM+1:gISEGM+NSEGM)    = int(SEG(1:nsegm)%NORME)
       ISEG(5,gISEGM+1:gISEGM+NSEGM)    = int(SEG(1:nsegm)%VECLIN)
       TABVOIS(1,gISEGM+1:gISEGM+NSEGM) = int(SEG(1:nsegm)%Voiso)
       TABVOIS(2,gISEGM+1:gISEGM+NSEGM) = int(SEG(1:nsegm)%Voise)
       gISEGM                           = gISEGM + NsegmCP

   enddo
endif
!!*****************************************

if (gISEGM.gt.gNSEGMAX) write (*,*) "Interface graphique : debordement de tableau dans le FORTRAN"

!*** Appel a la partie graphique

call PREDRAW(gISEGM,ISEG,TABVOIS,BVECLINCP,kkcp,gImoduR,IJONC,nbsegdep,nsegmaxcp)

!*** Appel a la partie graphique

!if (kk == 50) then
!   do i = 1, gisegm
!       write(*,'(6I10)')ISEG(1:5,I),IJONC(i)
!   enddo
!endif

if (kk.eq.0) stop "kk = 0 ?????"

end subroutine MODEGRAF

end module INTERGRA
