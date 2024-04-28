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

module BRICAMAT

use constantes_precision

implicit none

real(kind=DPR) :: alfa(3,3,6)
real(kind=DPR) :: wilbeta(6,3,3)
!
!SONGER A INITIALISER LES MATRICES DE PASSAGE AVEC VECTEN

integer(kind=DPI),parameter :: dimPREM = 100000
integer(kind=DPIM),save :: tabprem(dimPREM)
!tabprem(:)=0;tabprem(1)=1;tabprem(2)=2;
!
!SONGER A INITIALISER TABPREM

REAL(KIND=DPR),save::ref_vccp(3,3),&!Changement de repere
ref_v(3),&
ref_c(3),&
ref_cp(3),&
Rinv_ref_vccp(3,3)
!
!SONGER A INITIALISER LES VECTEURS ET MATRICE DE CHANGEMENT DE REPERE

!*** Aout 2003 : VARIABLE POUR LES CL BL ou P
integer,save::  CLx, CLy, CLz ! conditions aux limites, en x,y et z
logical,save::  BLx, BLy, BLz ! conditions de type bords libre, en x,y et z POUR LA TOPOLOGIE
logical,save:: cBLx,cBLy,cBLz ! conditions de type bords libre, en x,y et z POUR LES CONTRAINTES
logical,save:: tBL,tCLP,tCLC  ! Bords libre, Condtions aux limites periodiques ou confinées
!
#ifdef LSD
REAL(KIND=DPR),save::ModuR(3)! Dimension du volume simule
REAL(KIND=DPR),save::SHIFTxy,SHIFTxz,SHIFTyx,SHIFTyz,SHIFTzx,SHIFTzy! Decalages lors des CLP
REAL(KIND=DPR),save::ModuR2(3)! Dimension du volume simule avec facteur 2
REAL(KIND=DPR),save::SHIFTxy2,SHIFTxz2,SHIFTyx2,SHIFTyz2,SHIFTzx2,SHIFTzy2! idem
REAL(KIND=DPR),dimension(3),save  :: TranslaCLP
#endif
#ifndef LSD
INTEGER(KIND=DPI),save::ModuR(3)! Dimension du volume simule
INTEGER(KIND=DPI),save::SHIFTxy,SHIFTxz,SHIFTyx,SHIFTyz,SHIFTzx,SHIFTzy! Decalages lors des CLP
INTEGER(KIND=DPI),save::ModuR2(3)! Dimension du volume simule avec facteur 2
INTEGER(KIND=DPI),save::SHIFTxy2,SHIFTxz2,SHIFTyx2,SHIFTyz2,SHIFTzx2,SHIFTzy2! idem
INTEGER(KIND=DPI),dimension(3),save  :: TranslaCLP
#endif
!
!SONGER A INITIALISER avant utilisation de CLP

!*** Generateur uniforme de nombre pseudo-aleatoire ***
INTEGER(kind=DPIS),save :: seed1, seed2, seed3
!
!SONGER A INITIALISER avec init_seeds


!*** AffMat
real(kind=DPR) :: buffy1(1000,1000)
integer*4 :: buffN,ToToNumleg
logical :: totoMin,totoMax
real*8  :: totoMinR,totoMaxR

contains

#include"bricamat_entier_fraction_premier.F90"
#include"bricamat_PartEnt.F90"
#include"bricamat_vecteur.F90"

end module bricamat
