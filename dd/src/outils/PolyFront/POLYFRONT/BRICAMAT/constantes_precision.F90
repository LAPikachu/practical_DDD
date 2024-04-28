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

module constantes_precision

implicit none

integer,parameter  :: DPR=8  !type de réels utilisés : double precision : 8
integer,parameter  :: SPR=4  !type de réels utilisés : simple precision : 4

integer,parameter  :: DPI=8  !type de entier utilisés pour les vecteurs : 4
integer,parameter  :: DPIV=8 !type de entier utilisés pour les produits scalaires : 4 ou 8 selon la base
integer,parameter  :: DPIM=8 !type d'entier utilisés pour les produits matricels : 8

integer,parameter  :: DPIG=4 !type d'entier utilisés pour la communication avec l'interface graphique en C
integer,parameter  :: DPIS=4 !type d'entier utilisés pour rand
integer,parameter  :: DPIC=4 !type d'entier utilisés pour case

real(kind=DPR),parameter ::          brica_PII = 3.1415926536D0
real(kind=DPR),parameter :: brica_PII_vers_deg = 57.295779513D0

real(kind=DPR),parameter ::          PRECISION =  0.001 ! precision servant de comparaison entre les reels
real(kind=DPR),parameter ::          MRECISION = -0.001 ! sur les normes ie /maxCoorVecBase pour coordonnées

#ifdef LSD
real(kind=DPR), parameter :: LE_ZERO = 0.0D0
real(kind=DPR), parameter :: LE_UN = 1.0D0
real(kind=DPR), parameter :: ME_UN = -1.0D0
#endif
#ifndef LSD
integer(kind=DPI),parameter :: LE_ZERO = 0
integer(kind=DPI),parameter :: LE_UN = 1
integer(kind=DPI),parameter :: ME_UN = -1
#endif

end module constantes_precision
