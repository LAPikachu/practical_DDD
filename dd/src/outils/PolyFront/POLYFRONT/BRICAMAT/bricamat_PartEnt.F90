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

!##############################################################################
!#                                                                            #
!# Retourne la partie entiere d'un reel                                       #
!#                                                                            #
!# Int et Nint retourne autre chose...                                        #
!# Int  : 0.6 -> 0 et -0.6 ->  0 Partie entiere F90-95                        #
!# Nint : 0.6 -> 1 et -0.6 -> -1 Entier le plus proche                        #
!#                                                                            #
!# PartEnt : 0.6 -> 0 et -0.6 -> -1Entier N le plus grand tel que N < R       #
!#                                                                            #
!# -2.999992 -> -3, on essaye de tenir compte de pbs de precision numerique   #
!#                                                                            #
!##############################################################################

function PartEnt(r)

implicit none

real(DPR),intent (in) :: r
integer(kind=DPI)     :: PartEnt

if (int(r,kind=DPI).ne.nint(r,kind=DPI).and.(abs(nint(r,kind=DPI)-r).le.precision)) then

if (r.ge.0.d0) then
PartEnt=nint(r,kind=DPI)
else
PartEnt=nint(r,kind=DPI)-1
endif

else

if (r.ge.0.d0) then
PartEnt=int(r,kind=DPI)
else
PartEnt=int(r,kind=DPI)-1
endif

endif

return

end function PartEnt
