
!===================================================================================================
!========================    DEBUT    MODULE  "BIGSAVE"   ==========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the simulation output writing and restart files.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module BIGSAVE

use BRICAMAT
use VARGLOB

implicit none

contains

!***********************************************************************

subroutine saveall

implicit none
# ifdef MDC
character(len=256) :: cmd
# endif
#ifdef PA
  if ((Mon_Rang_PairImpair + IUN) >  IUN) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  else
#endif

!*** CREATION/OUVERTURE-ECRASSEMENT DU FICHIER BINAIRE
# ifdef MDC

write(cmd,'("cp ../out/bigsave.bin ../out/bigsave.bin.old")')
call system(cmd)
# endif

Print *, " >>>>>>>>>>>   BIGSAVE at step KK = ",KK

open(99,FILE=file_bigsave,FORM='UNFORMATTED',STATUS='REPLACE')

Write(99) KK,nsegm,nsegmax,Npar,Nstatcontrol,Sigma,EPSO,EPSDOT,DELSIG,accutime
Write(99) SIGAPP,SIGPLUS,DELTAT,SONDMS,VITMOY,ModuR(1:3),FSIG
Write(99) Solli_Sys(1:NTSG)
Write(99) SEG(1:nsegm)
Write(99) CranSys(1:NTSG)
Write(99) RAUSYS(1:NTSG),AireSYS(1:NTSG),GAMMADOTSYS(1:NTSG),TrAppSys(1:NTSG),TrIntSys(1:NTSG)
Write(99) RAUDMO_I,RAUDIS,RAUDMO,AireVisSys(1:NTSG),AireCoinSys(1:NTSG),AireMixtSys(1:NTSG)
Write(99) EPSOLD(1:NstatControl),GAMMAOLDSYS(1:NTSG,1:Nstatcontrol)
Write(99) NTGD,NTGD1,NTGD2,OLDNTGD1,OLDNTGD2,airevis,airecoin,airemixte,Stress_fact

if (DesorientGrain) Write(99) AireSysGrain(1:nbgrains,1:NTSG)

if (npar > izero)   Write(99) par(1:npar)

if ((.not. allocation_dynamique_boites) .and. calculate_gammabox) then
  Write(99) NBoites
  Write(99) Gammabox(1:NBoites,1:NTSG)
endif

close(99)

#ifdef PA
    ! The writing process ended their task, we can now liberate everybody
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
#endif

end subroutine saveall

!***********************************************************************
!	AVEC MODIF CRADE POUR REDEFINIR LE CONTROL EN RELANCANT !!!
!***********************************************************************
subroutine loadlight

implicit none

!*** OUVERTURE DU FICHIER BINAIRE
open(99,FILE=file_bigsave,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')

!*** Generateur uniforme de nombre pseudo-aleatoire ***
read(99) seed1, seed2, seed3
!******************************************************
!*** Microstructure
read(99) NSEGM,SEG ! Nombre total de segments, Tableaux des segments

close(99)

end subroutine loadlight

!########################################################################
!#  The film.bin info are written at this step                          #
!#  Notice that the header of the film.bin file is written only at the  #
!#  the beginning of the file (first opening in Initialisation_fichiers #
!########################################################################
subroutine SAUVima

implicit none

integer,parameter   :: DPI_S=selected_int_kind(9)       ! The integer format needed for the film.bin file
integer(kind=DPI_S) :: I                                !<
integer(kind=DPI_S) :: O                                !<
integer(kind=DPI_S) :: tmp                              !<
integer(kind=DPI_S), Dimension(3) :: P                  !<
integer(kind=DPI_S), Dimension(3) :: l                  !<
integer(kind=DPI_S), Dimension(9,nsegm+npar) :: tabseg  !<

!*** Construction des variables d'echanges utiles

#ifdef PA

if ((Mon_Rang_PairImpair + IUN) >  IUN) then
  ! Only the process zero of the communicator can write in output files, the others are waiting
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
else

#endif

do I = 1,int(NSEGM)

  P(1:3) = int(SEG(i)%O(1:3))
  l(1:3) = int(SEG(i)%NORME/2*bveclin(1:3,SEG(i)%veclin))
  O = P(1) + l(1)
  if (modulo(O,modur(1)) /= O) P(1) = int(modulo(O,modur(1))- l(1))

  O = P(2) + l(2)
  if (modulo(O,modur(2)) /= O) P(2) = int(modulo(O,modur(2))- l(2))

  O = P(3) + l(3)
  if (modulo(O,modur(3)) /= O) P(3) = int(modulo(O,modur(3)) - l(3))

  tabseg(1:3,i) = P(1:3)
  tabseg(4,i)   = int(SEG(i)%NORME)
  tabseg(5,i)   = int(SEG(i)%VECLIN)

  if(seg(i)%jonc) then
    ! This segment take part to a junction, the ijonc info is saved
    tabseg(6,i) = int(seg(i)%ijonc)
  elseif (seg(i)%gd > izero) then
    ! This segment is a GD (screw segment at cross-slip intersection)
    tabseg(6,i) = -Iun
  else
    ! Nothing special to declare
    tabseg(6,i) = Izero
  endif

  tabseg(7,i) = int(SEG(i)%vnno)
  tabseg(8,i) = int(SEG(i)%vnne)
  tabseg(9,i) = int(SEG(i)%surface)

enddo

if(particules) then

  tmp = int(nsegm)

  do I = 1, int(npar)
    tmp = tmp + 1
    tabseg(1:3,tmp) = int(par(i)%C(1:3))
    tabseg(4,tmp)   = nint(par(i)%R/normlin(IUN))
    tabseg(5,tmp)   = IUN
    tabseg(6,tmp)   = Itrois
    tabseg(7,tmp)   = Izero
    tabseg(8,tmp)   = Izero
    tabseg(9,tmp)   = Izero
  enddo

endif

! The writting ...

#ifdef MDC
   write(93) int(kk,DPI_S),real(deltat,DP),int((nsegm+npar),DPI_S),1D-6,1D-6
#else
   write(93) int(kk,DPI_S),real(deltat,DP),int((nsegm+npar),DPI_S),sigma*xmu*1D-6*schmid,EPSO*100!
#endif

write(93) (int(tabseg(1:9,I),DPI_S), I=1,int(nsegm+npar))

write(93) '&'      ! This mark is used to note the end of each step info

#ifdef PA
  ! The writing process ended their task, we can now liberate everybody
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
endif
#endif

end subroutine SAUVima

!##########################################################################
!# subroutine save_config: ecrite le ficher SEG correspondant a la config #
!# a l'iteration donnee : les GD et JONC ne sont pas pris en compte       #
!###########################################################14/04/99#######
subroutine save_config

implicit none

integer(DPI)  :: i,vo,ve,vnne,vnno,ijonc
logical :: zipseg(nsegm)
character(len=40)   :: file_SEG_relax, file_SEG_save, file_SEG_surface

zipseg(1:nsegm) = .false.

61 format(I7,3I10,I7,I7,4I7,L3,I7,I3)
64 format(I7,3I10,I7,I7,4I7,L3,I7,I3,I3)

! le block suivant sert a sauver la config apres la relaxation TL dans SEG_relax

file_SEG_relax    =   '../in/SEG_relax'
file_SEG_save     =   '../in/SEG_save'
file_SEG_surface  =   '../in/SEG_surface'
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_SEG_relax    =   file_SEG_relax(1:len_trim(file_SEG_relax))//"2"
  file_SEG_save     =   file_SEG_save(1:len_trim(file_SEG_save))//"2"
  file_SEG_surface  =   file_SEG_surface(1:len_trim(file_SEG_surface))//"2"
endif

if(KK <= relax_Reac) then
  print *, "=============   SAUVGARD de configuration dans SEG_relax kk = ",KK
  open( 8,FILE=file_SEG_relax,STATUS='UNKNOWN')

  write (8,'(12F4.0)') Solli_sys(1:NTSG)
  write (8,*) NSEGM
  write (8,*) MODUR(1:3)
else
  print *, "=============   SAUVGARD de configuration dans SEG_save  KK = ",KK
  open( 8,FILE=file_SEG_save,STATUS='UNKNOWN')
  open( 9,FILE=file_SEG_surface,STATUS='UNKNOWN')

  write (8,'(12F4.0)') Solli_sys(1:NTSG)
  write (8,*) NSEGM
  write (8,*) MODUR(1:3)
  write (9,'(12F4.0)') Solli_sys(1:NTSG)
  write (9,*) NSEGM
  write (9,*) MODUR(1:3)
endif

! The loop on segments
do i = 1, nsegm

  vo = seg(i)%voiso
  vnno = seg(i)%vnno
  ve = seg(i)%voise
  vnne = seg(i)%vnne
  ijonc = seg(i)%Ijonc

  ! Sometime it is useful to restart the simulation from the SEG_save configuration and to
  ! change the value of nsegmax. Hence Nsegmax value is replaced by zero here
  if (ve == nsegmax) ve = izero
  if (vo == nsegmax) vo = izero
  if (vnno == nsegmax) vnno = izero
  if (vnne == nsegmax) vnne = izero
  if (ijonc == nsegmax) ijonc = izero

  write (8,61) i, seg(i)%O(1:3), seg(i)%norme, seg(i)%veclin, vo, vnno, ve, vnne, seg(I)%JONC, Ijonc, seg(i)%gd

  if(KK > relax_Reac) then
    write (9,64) i, seg(i)%O(1:3), seg(i)%norme, seg(i)%veclin, vo, vnno, ve, vnne, seg(I)%JONC, Ijonc, seg(i)%gd, seg(i)%surface
  endif
enddo

close(8)
close(9)

end subroutine save_config

end module BIGSAVE

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
