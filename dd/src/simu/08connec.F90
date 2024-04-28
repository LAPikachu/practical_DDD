
!===================================================================================================
!========================    DEBUT    MODULE "CONNEC"    ===========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the segments connectivity rules.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module CONNEC

use VARGLOB
use DEBUG

implicit none

contains

!############################################################################
!# subroutine init_seg(i) : initialize segment number I                     #
!# very useful when adding a new segments. All entries should be defined    #
!############################################################################
subroutine Init_seg(i)

implicit none

integer(DPI),intent(in) :: I

if (kkdebug) then
  write(379,fmt='("*** Entering Init_seg for new_segment = ",I7,":")') i
endif

seg(I)%O(:)         = Izero
seg(I)%VECLIN       = Izero
seg(I)%NORME        = Izero
seg(I)%VOISO        = NSEGMAX
seg(I)%VOISE        = NSEGMAX
seg(I)%VNNO         = NSEGMAX
seg(I)%VNNE         = NSEGMAX
seg(I)%IJONC        = NSEGMAX
seg(I)%TJONC        = Izero
seg(I)%WAIT         = Izero
seg(I)%GRAIN        = Izero
seg(I)%SURFACE      = Izero
seg(I)%DOM          = Izero
seg(I)%VarFreePlan  = Izero
seg(I)%RESDEP       = zero
seg(I)%DEPINST      = zero
seg(I)%probadev     = zero
seg(I)%JONC         = .false.
seg(I)%GD           = izero
seg(I)%DISEG        = .false.
seg(I)%BLOQUER      = .false.
seg(I)%unload       = .false.

#ifdef MDC
seg(I)%SIGFE(:)     = zero
#endif

end subroutine init_seg

!#############################################################################################
!> subroutine stopjonc : Destroy junctions when the latter have a zero length.
!! The goal of this procedure is mainly to check if the segment at O and E can be reconnected.
!! \param itemp   Index of the segment junction to eliminate
!! \param ref     The call reference
!#############################################################################################
subroutine stopjonc(itemp,ref)

implicit none

integer(kind=DPI) :: itemp    !< The segment junction index
integer(kind=DPI) :: Io       !< The first neighbor on the O side
integer(kind=DPI) :: Ie       !< The first neighbor on the E side
integer(kind=DPI) :: lo       !< The length of segment Io
integer(kind=DPI) :: le       !< The length of segment Ie
integer(kind=DPI) :: VLo      !< The veclin index of segment Io
integer(kind=DPI) :: Vle      !< The veclin index of segment Ie
integer(kind=DPI) :: icont    !< A loop counter
integer(kind=DPI) :: VTest(3) !< A vector needed for some test
integer           :: ref      !< The call reference

! Paranoid tests
if (seg(Itemp)%norme /= 0) write(379,*) "WARNING kk",kk, ", StopJonc subroutine is called with a segment of a non-zero length!"
if (.not. seg(Itemp)%jonc) write(379,*) "WARNING kk",kk, ", StopJonc subroutine is called For a not jonc segment!"

if(kkdebug) then
  write(379,*) "Subroutine StopJonc is called for segment :",itemp,"  at ref =",ref
  call conf(itemp)
!   call seginfo(itemp," At the begining of stopjonc ")
endif

! The junction definition is eliminated
SEG(itemp)%jonc   = .FALSE.
SEG(itemp)%ijonc  = NSEGMAX
SEG(itemp)%tjonc  = 0

! The first neighbors
Io = seg(itemp)%voiso
Ie = seg(itemp)%voise

! If the segment is a segment to eliminate, this is enough.
if ((Io == itemp .and. Ie == itemp) .or. out(itemp))  return

! To reconstruct the connectivity, we look for the Io=vnno and Ie=vnne
icont = iun
Do while (Io /= Nsegmax         .and. &
          seg(Io)%norme == 0    .and. &
          .not. seg(Io)%jonc    .and. &
          seg(Io)%gd < iun)

  Io    = seg(Io)%voiso
  icont = icont + iun

  if (out(Io)) then
    write(379,*) "In stopjonc ref:",ref,"a seg to eliminate is listed. Io:",Io
    call conf(Io)
    stop "In subroutine stopjonc a problem of connectivity is found!"
  endif

  if (icont > 10 ) then
    call conf(itemp)
    call seginfo(itemp," ERREUR FATALE001 ")
    write(*,*) "Step kk=",kk," stopjonc: strange connectivity on the O side of segment i= ",itemp
    stop
  endif
enddo

icont = iun
Do while (Ie /= Nsegmax             .and. &
          seg(Ie)%norme == 0        .and. &
          .not. seg(Ie)%jonc        .and. &
          seg(Ie)%gd < iun)

  Ie    = seg(Ie)%voise
  icont = icont + iun

  if (out(Ie)) then
    write(379,*) "In stopjonc ref:",ref,"a seg to eliminate is listed. Ie:",Ie
    call conf(Ie)
    stop "In subroutine stopjonc a problem of connectivity is found!"
  endif

  if (icont > 10 ) then
    call conf(itemp)
    call seginfo(itemp," ERREUR FATALE002 ")
    write(*,*) "Step kk=",kk," stopjonc: strange connectivity on the E side of segment i= ",itemp
    stop
  endif
enddo

! The vnno and vnne to reconnect
VLo = seg(Io)%veclin
VLe = seg(Ie)%veclin
Lo  = seg(Io)%norme
Le  = seg(Ie)%norme

! We look for a possible overlapping length of the segment arms attached to the junction to eliminate
! the segment arms must not be jonc or gd segments
if( Io /= nsegmax .and. .not. seg(Io)%jonc .and. (seg(Io)%gd + seg(Ie)%gd) < iun .and. &
    Ie /= nsegmax .and. .not. seg(Ie)%jonc) then

  VTest(:) = bveclin(:,VLo) + Bveclin(:,VLe)
  ! if an overlapping exist
  if ((abs(VTest(1)) + abs(VTest(2)) + abs(VTest(3))) == 0) then

    if (kkdebug) write(379,*) " Overlapping of the vnno=",Io," and vnne=",Ie," with length:", Lo, le

    if (Lo > Le) then               ! Io is longer than Ie : It is connected to the vnne of ie

      Lo            = Lo - Le
      seg(Io)%norme = Lo
      seg(Ie)%norme = Izero
      Ie            = seg(ie)%vnne

    elseif (Lo < Le) then           ! Ie is longer than Io : It is connected to the vnno of io

      Le            = Le - Lo
      seg(Ie)%norme = Le
      seg(Ie)%O(:)  = seg(Io)%O(:)
      seg(Io)%norme = Izero
      Io            = seg(io)%vnno

    else                            ! The length of Io and Ie is the same, we connect their respective vnno and vnne

      seg(Io)%norme = Izero
      seg(Ie)%norme = Izero
      Io            = seg(io)%vnno
      Ie            = seg(ie)%vnne

    endif

  endif

endif

! The connection between Io=vnno and Ie=vnne is made
call connecij(Io,Ie,1000000+ref)

! if(kkdebug) then
!   call seginfo(itemp," At the end of stopjonc ")
! endif

end subroutine stopjonc


!############################################################################
! subroutine voisinageE (i,j) : update of the definition of vnn (non zero
! neighbor) in seg(i) and its neighbors in the particular case where any new
! segments must be consider. It is applied typically when a segment length
! become zero or the opposite.
! Paranoid debug test must be avoid as much as possible in this subroutine
! since it is the most called subroutine in mM
!############################################################################
subroutine voisinage(i,ref)

implicit none

integer(kind=DPI)   :: i              !< The segment index to consider
integer(kind=DPI)   :: lon            !< The segment length
integer(kind=DPI)   :: ve             !< First neighbor at the end
integer(kind=DPI)   :: vo             !< First neighbor at the origin
integer(kind=DPI)   :: vnne           !< First non zero length neighbor at the end
integer(kind=DPI)   :: vnno           !< First non zero length neighbor at the end
integer(kind=DPI)   :: vmemo          !< An intermediate variable useful in neighbors research

integer             :: ref            !< The debug ref
integer             :: compt          !< variable to compute number of elements
integer             :: icompt         !< variable to compute number of elements

Logical             :: Ivnn           ! Segment i is not a kneecap or GD or JONC segment

! Paranoid test
if (i /= nsegmax .and. ((i==seg(i)%voiso .and. i==seg(i)%voise) .or. out(i))) then
  write(379,*)  " subroutine voisinage is called with ref=", ref,                           &
                "for a segment I to eliminate", i, seg(i)%voiso, seg(i)%voise, out(i)
  call seginfo(i,'A segment to eliminate is called in subroutine voisinage')
  stop "A segment to eliminate is called in subroutine voisinage"
endif

! Test to know if I is a segment to define as a vnn = not a pivot segment or GD or Jonc
Ivnn = ((seg(i)%norme /= izero) .or. seg(i)%gd > izero .or. seg(i)%jonc)

If (kkdebug) then
  write(379,*) "Voisinage: Update of the VNN definition for segments close to i =", i, "   REF =", ref
  if(Ivnn)  then
    write(379,*)  I, " is not a kneecap or a GD or a JONC segment "
  else
    write(379,*)  I, " is a kneecap and not a GD and not a JONC segment "
  endif
endif

!----------------------------------
! Identification of the seg(i)%vnne
ve    = seg(i)%voise
lon   = seg(ve)%norme

compt = iun
do

  ! test to identify the vnne
  if (lon /= 0 .or. seg(ve)%gd > izero .or. seg(ve)%jonc) exit

  ! the particular case of a closed loop made of segments made only of pivot segments
  ! This loop can be eliminated
  if (ve == i) then

    ve = seg(i)%voise
    do icompt = 1, compt
      vmemo             = ve
      out(vmemo)        = .true.
      seg(vmemo)%voiso  = vmemo
      ve                = seg(vmemo)%voise
      seg(vmemo)%voise  = vmemo
    enddo

    return   ! No more work to do

  endif

  if ( ve == seg(ve)%voise .or. compt > 11 ) then
    print *, " Voisinage E : fatal error 1, KK=",KK,"  REF= ",ref,"  Compt= ",compt,"  i= ",i
    call seginfo(i,'Problem in voisinage on the E side')
    stop 'see ../out/debug file !!!'
  endif

  ve    = seg(ve)%voise
  lon   = seg(ve)%norme
  compt = compt + 1

enddo

vnne        = ve
seg(i)%vnne = vnne

!----------------------------------
! Identification of the seg(i)%vnno
vo    = seg(i)%voiso
lon   = seg(vo)%norme
compt = iun

do

  ! test to identify the vnno
  if (lon /= 0 .or. seg(vo)%gd > izero .or. seg(vo)%jonc) exit

  ! the particular case of a closed loop made of segments with length zero
  ! This segments must be eliminated
  if (vo == i) then

    vo = seg(i)%voiso
    do icompt = 1, compt
      vmemo             = vo
      out(vmemo) = .true.
      seg(vmemo)%voise = vmemo
      vo = seg(vmemo)%voiso
      seg(vmemo)%voiso = vmemo
    enddo

    return   ! No more work to do

  endif

  if( vo == seg(vo)%voiso .or. compt > 11 ) then
    print *, " Voisinage O : fatal error 1, KK=",KK,"  REF= ",ref,"  Compt= ",compt,"  i= ",i
    call seginfo(i,'Problem in voisinage on the O side')
    stop 'see ../out/debug file !!!'
  endif

  vo    = seg(vo)%voiso
  lon   = seg(vo)%norme
  compt = compt + 1

enddo

vnno        = vo
seg(i)%vnno = vnno

!-----------------------------------------------------------------------------------------------
! Update of the %vnno and %vnne of the pivot segments between seg(i)%vnno and seg(i)%vnne
if(Ivnn) then   ! I is not a kneecap or a GD or a JONC segment

  if(vnne /= nsegmax)   seg(vnne)%vnno = i
  if(vnno /= nsegmax)   seg(vnno)%vnne = i

  ! the E side
  ve = seg(i)%voise
  do while (ve /= vnne)
    seg(ve)%vnne  = vnne
    seg(ve)%vnno  = i
    ve            = seg(ve)%voise
  enddo

  ! the O side
  vo = seg(i)%voiso
  do while (vo /= vnno)
    seg(vo)%vnno  = vnno
    seg(vo)%vnne  = i
    vo            = seg(vo)%voiso
  enddo

else      ! I is a kneecap and not a GD and not a JONC segment

  if(vnne /= nsegmax)   seg(vnne)%vnno = vnno
  if(vnno /= nsegmax)   seg(vnno)%vnne = vnne

  ! the E side
  ve = seg(i)%voise
  do while (ve /= vnne)
    seg(ve)%vnne  = vnne
    seg(ve)%vnno  = vnno
    ve            = seg(ve)%voise
  enddo

  ! the O side
  vo = seg(i)%voiso
  do while (vo /= vnno)
    seg(vo)%vnno  = vnno
    seg(vo)%vnne  = vnne
    vo            = seg(vo)%voiso
  enddo

endif

end subroutine voisinage

!====================================================================================================
! This procedure is used to insert a kneecap between the origin of I and a pinning point
!====================================================================================================
subroutine RotuleO(i,ref)

implicit none

integer(kind=DPI)   :: i,itemp,IVL
integer             :: ref

! Initialization
plusseg = plusseg + IUN
itemp = nsegm + plusseg

if(kkdebug)                                                                               &
write(379,'("Insert Kneecap number ",I5," between Nsegmax and ",I5," REF = ", I10)') itemp,I,ref

IVL = seg(i)%veclin

seg(itemp)%veclin       = nbrot(IVL,IVL,2)
seg(itemp)%norme        = IZERO
seg(itemp)%O(:)         = seg(I)%O(:)
seg(itemp)%voise        = I
seg(itemp)%voiso        = NSEGMAX
seg(itemp)%vnne         = I
seg(itemp)%vnno         = NSEGMAX
seg(itemp)%jonc         = .false.
seg(itemp)%Ijonc        = NSEGMAX
seg(itemp)%tjonc        = IZERO
seg(itemp)%GD           = izero
seg(itemp)%wait         = IZERO
seg(itemp)%resdep       = zero
SEG(Itemp)%bloquer      = seg(I)%bloquer
SEG(Itemp)%unload       = seg(I)%unload
seg(itemp)%grain        = seg(I)%grain
seg(itemp)%surface      = IZERO
seg(itemp)%VarFreeplan  = IZERO
seg(Itemp)%NPhase       = IZERO
seg(itemp)%anglevis     = -un

if (NbCvxDom > IUN) then
  call assign_domain(itemp,1111)
else
  seg(itemp)%dom=IUN
endif

#ifdef MDC
seg(itemp)%SIGFE(:) = zero
#endif

! mise a jours de J
seg(I)%voiso            = itemp
seg(I)%vnno             = NSEGMAX

end subroutine RotuleO

!====================================================================================================
! This procedure is used to insert a kneecap between a pinning point and the end of I
!====================================================================================================
subroutine RotuleE(i,ref)

implicit none

integer(kind=DPI)  :: i,itemp,IVL
integer            :: ref

! Initialization
plusseg = plusseg + IUN
itemp = nsegm + plusseg

if(kkdebug) &
     write(379,'("Insert Kneecap number ",I5," between ",I5," and NSEGMAX REF = ", I10)')&
     itemp,I,ref

IVL = seg(i)%veclin
seg(itemp)%veclin       = nbrot(IVL,IVL,2)
seg(itemp)%norme        = IZERO
seg(itemp)%O(:)         = seg(I)%O(:) + seg(i)%norme * bveclin(1:3,IVL)
seg(itemp)%voise        = NSEGMAX
seg(itemp)%voiso        = I
seg(itemp)%vnne         = NSEGMAX
seg(itemp)%vnno         = I
seg(itemp)%jonc         = .false.
seg(itemp)%Ijonc        = NSEGMAX
seg(itemp)%tjonc        = IZERO
seg(itemp)%GD           = izero
seg(itemp)%wait         = IZERO
seg(itemp)%resdep       = zero
SEG(Itemp)%bloquer      = seg(I)%bloquer
SEG(Itemp)%unload       = seg(I)%unload
seg(itemp)%grain        = seg(I)%grain
seg(itemp)%surface      = IZERO
seg(itemp)%VarFreeplan  = IZERO
seg(Itemp)%NPhase       = IZERO
seg(Itemp)%anglevis     = -un

if (NbCvxDom > IUN) then
  call assign_domain(itemp,1112)
else
  seg(itemp)%dom=IUN
endif

#ifdef MDC
seg(I)%SIGFE(:)         = zero
#endif

! mise a jours de J
seg(I)%voise            = itemp
seg(I)%vnne             = NSEGMAX

end subroutine RotuleE

!====================================================================================================
! This procedure is very important: it defines the way segments must be connected to each other.
! If I and J belong to the same slip system, they are directly connected with appropriate number
! of kneecap (if needed)
! If I and J are not belonging to the same slip system, kneecap CS are introduced between I and J
! unless I or J are already CS
! In the last case, if ref value is negative... then when a gd segment is introduced
! the latter is forced to be of type gd2
!====================================================================================================
subroutine connecij(i,j,ref)

implicit none

integer(kind=DPI),intent(in)    :: i        !< Segment I
integer(kind=DPI),intent(in)    :: j        !< Segment J
integer(kind=DPI)               :: itemp    !< Temporary index for segments
integer(kind=DPI)               :: rotules  !< Number of kneecap to add
integer(kind=DPI)               :: ve       !< First neighbor at End
integer(kind=DPI)               :: vo       !< First neighbor at Origin
integer(kind=DPI)               :: rvnno    !< First non null neighbor at End
integer(kind=DPI)               :: rvnne    !< First non null neighbor at Origin
integer(kind=DPI)               :: X(3)
integer(kind=DPI)               :: IVLiini  !< Backup Segment I veclin
integer(kind=DPI)               :: IVLjini  !< Backup Segment J veclin
integer(kind=DPI)               :: IVLi     !< Store Segment I veclin
integer(kind=DPI)               :: IVLj     !< Store Segment J veclin
integer(kind=DPI),dimension(3)  :: Oi       !< Segment I Origin coordinates
integer(kind=DPI),dimension(3)  :: Ei       !< Segment I End coordinates
integer(kind=DPI)               :: lgd1
integer(kind=DPI)               :: lgd2
integer(kind=DPI)               :: seggd
integer(kind=DPI)               :: lgd
integer(kind=DPI)               :: plussegini
integer(kind=DPI)               :: memoGD1,memoGD2
integer                         :: ref     !< Subroutine tracker

logical                         :: glisdev  !< Is cross-sliping
logical                         :: insert

! Paranoid tests (attention segment nsegmax can be called in connecij)
if (i /= nsegmax .and. ((i==seg(i)%voiso .and. i==seg(i)%voise) .or. out(i))) then
  write (379,*) "at step:",KK,"connecij is called with ref:",ref,"with a segment i to eliminate"
  call conf(i)
  call seginfo(i,"Problem with i at connecij begining")
  stop 'A segment to eliminate is call in connecij ?'
endif
if (j /= nsegmax .and. ((j==seg(j)%voiso .and. j==seg(j)%voise) .or. out(j))) then
  write (379,*) "at step:",KK,"connecij is called with ref:",ref,"with a segment j to eliminate"
  call conf(j)
  call seginfo(j,"Problem with i at connecij begining")
  stop 'A segment to eliminate is call in connecij ?'
endif

if (kkdebug) then
  write(379,*)  "--"
  write (379,fmt='("Enter connecij with ref:", I7,"   i=", I9,": and j=" ,I6,":")') ref,i,j
  write (379,fmt='("I | VD ",I9," Neighboring | "," VNNO =",I9,": VNNE =", I9,":")') &
                  seg(i)%veclin,seg(i)%vnno,seg(i)%vnne
  write (379,fmt='("J | VD ",I9," Neighboring | "," VNNO =",I9,": VNNE =", I9,":")') &
                  seg(j)%veclin,seg(j)%vnno,seg(j)%vnne
endif

if (GB==ITROIS) then
  if (i/=nsegmax.and.j/=nsegmax)  then
    if  (seg(i)%grain==IZERO.or.seg(j)%grain==IZERO.or.seg(i)%grain/=seg(j)%grain) then
      print*,'i is in grain',seg(i)%grain,'j is in grain',seg(j)%grain
      print*,'Is the variable %grain has been initialized ?',ref
      if (seg(i)%grain==IZERO) then
        seg(i)%grain=seg(j)%grain             ! can be put right
      else if (seg(j)%grain==IZERO) then
        seg(j)%grain=seg(i)%grain             ! can be put right
      else
        stop "A segments do not belong to its grain!!!"
      endif
    endif
  endif
endif

! Initialization
plussegini = plusseg

memoGD1 = izero    ! Var used to keep info about GD1 (no GD1 = 0)
memoGD2 = iun      ! Var used to keep info about GD2 (no GD2 = 1)

! We check if i or j is a segment GD made by cross-slip (gd = 2)
if (seg(i)%gd == ideux) memoGD2 = memoGD2 + iun
if (seg(j)%gd == ideux) memoGD2 = memoGD2 + iun
if (seg(i)%gd == iun) memoGD1 = memoGD1 + iun
if (seg(j)%gd == iun) memoGD1 = memoGD1 + iun
! Ref is negatif = a GD2 segment must be reintroduced between i and j
if (ref < izero) then
  GD2made = .false.         ! if ref < 0, we always want to build a GD2 if possible
  memoGD2 = memoGD2 + iun
endif

! if (kkdebug) write(379,*) '->1 ',i,j,memoGD1,memoGD2,oldntgd1,oldntgd2

! initialisation necessaire pour les nouvelles rotules
itemp = nsegm + plusseg
seg(itemp + 1 : Itemp + 10)%gd          = izero
seg(Itemp + 1 : Itemp + 10)%jonc        = .false.
seg(Itemp + 1 : Itemp + 10)%Ijonc       = NSEGMAX
seg(Itemp + 1 : Itemp + 10)%tjonc       = izero
seg(itemp + 1 : itemp + 10)%resdep      = zero
seg(itemp + 1 : itemp + 10)%bloquer     = seg(i)%bloquer
seg(itemp + 1 : itemp + 10)%unload      = seg(i)%unload
seg(itemp + 1 : itemp + 10)%wait        = IZERO
seg(itemp + 1 : itemp + 10)%probadev    = ZERO
seg(itemp + 1 : itemp + 10)%surface     = IZERO
seg(itemp + 1 : itemp + 10)%VarFreePlan = IZERO
seg(itemp + 1 : itemp + 10)%NPhase      = izero
seg(itemp + 1 : itemp + 10)%dom         = seg(j)%dom
seg(itemp + 1 : itemp + 10)%anglevis    = -un

#ifdef MDC
seg(itemp + 1 : itemp + 10 )%SIGFE(1) = zero
seg(itemp + 1 : itemp + 10 )%SIGFE(2) = zero
seg(itemp + 1 : itemp + 10 )%SIGFE(3) = zero
seg(itemp + 1 : itemp + 10 )%SIGFE(4) = zero
seg(itemp + 1 : itemp + 10 )%SIGFE(5) = zero
seg(itemp + 1 : itemp + 10 )%SIGFE(6) = zero
#endif

if(i == nsegmax) then
  seg(itemp + 1 : itemp + 10)%grain = seg(j)%grain
else
  seg(itemp + 1 : itemp + 10)%grain = seg(i)%grain
endif

insert = .false.

if(kkdebug)  then
  write(379,fmt='("====Connecij : i =",I9,": et j =", I9,":  / Nseg++seg =",I6,"  REF = ", I5)') i,j, itemp, ref
endif

!================================================================================
! cas ou I et J sont pts 'encrage : on ne fait rien .............................
!================================================================================
if (i == nsegmax .and. j == nsegmax) then
  write( *, fmt='(" REF = ", I5, "  KK =",I10)') ref,kk
  write( *, fmt='("connecij: Can not connect two anchor points")')
  stop
endif

ivliini = seg(i)%veclin
ivljini = seg(j)%veclin
ivli = ivliini
ivlj = ivljini

!================================================================================
! When I is a pinning point, only J must be updated
!================================================================================
if (i == nsegmax) then

  vo = seg(j)%voiso
  rotules = 0

  do while (vo /= i )

    rotules = rotules + 1
    itemp = vo
    vo = seg(vo)%voiso

    if(kkdebug) write(379,*)  "Segment ", itemp, ":  will be remove 1"

    if (seg(itemp)%norme /= 0) then
      write(379,*) " error-1 connecij KK=",kk,": REF = ", ref, " remove = ", itemp, ":  non kneecap segment"
      seg(itemp)%norme = Izero
    endif

    if (seg(itemp)%gd == ideux) then
      memoGD2 = memoGD2 + iun                 ! A seg gd is eliminated, we want to remember this info
      gd2elimin = gd2elimin + 1               ! The info about eliminated gd2 is keept for multi connecij call
      if (kkdebug) write(379,*) '>>>>>1',gd2elimin
    endif
    if (seg(itemp)%gd == iun) memoGD1 = memoGD1 + iun

!     if (kkdebug) write(379,*) '->2 ',i,j,memoGD1,memoGD2,oldntgd1,oldntgd2

    seg(itemp)%voise = itemp
    seg(itemp)%voiso = itemp
    out(itemp) = .true.
    seg(itemp)%resdep = zero
#ifdef MDC
    seg(itemp)%SIGFE(:) = zero
#endif
    IDEP(itemp) = Izero

    if (seg(itemp)%surface > IZERO .and. seg(itemp)%surface < ITROIS) then
      seg(j)%surface=seg(itemp)%surface
      seg(j)%VarFreeplan=seg(itemp)%VarFreePlan
      if (seg(itemp)%surface==IUN) seg(j)%dom=seg(itemp)%dom !A segment touching a surface with the origin has the domain of the surface
      seg(itemp)%surface=IZERO
      seg(itemp)%VarFreePlan=IZERO
    endif

    if(seg(Itemp)%jonc) then
      lgd1 = seg(Itemp)%Ijonc
      seg(lgd1)%jonc = .false.
      seg(lgd1)%Ijonc = nsegmax
      if(seg(lgd1)%norme == 0) call voisinage(lgd1,344)
!          call seginfo(lgd1, " verification de connection pour le binome ")
      seg(itemp)%jonc = .false.
      seg(itemp)%Ijonc = NSEGMAX
      write(379,*) "!> Connecij: 2: we remove the junction ",itemp,": and its pair: ", lgd1,":"
      write(379,*) "KK=",KK, "            REF = ",REF
    endif

    seg(itemp)%gd = izero
    idep(itemp)   = 0

    if (rotules > 100 ) then
      call conf(i)
      call conf(j)
      write(*,*) " REF = ", ref, "  KK =",kk
      write(*,*) "STOP kk",kk, "Infinite loop 1 in connecij"
      stop
    endif

  enddo

  ! We add a kneecap between segment J and NSEGMAX
  PlusSeg = PlusSeg +1 ; itemp = nsegm + PlusSeg
  if(KKdebug) write(379,*)  " We add a kneecap between segment J and NSEGMAX :",itemp

  seg(itemp)%veclin  = nbrot(IVLJ,IVLJ,2)
  seg(itemp)%norme   = 0
  seg(itemp)%O(:)    = seg(J)%O(:)
  seg(itemp)%voise   = j
  seg(itemp)%voiso   = i
  seg(itemp)%vnne    = j
  seg(itemp)%vnno    = i
  seg(itemp)%jonc    = .false.
  seg(itemp)%Ijonc   = NSEGMAX
  seg(itemp)%tjonc   = 0
  seg(itemp)%gd      = izero
  seg(itemp)%wait    = IZERO
  SEG(Itemp)%bloquer = seg(j)%bloquer
  SEG(Itemp)%unload  = seg(j)%unload
  seg(itemp)%grain   = seg(j)%grain
  seg(itemp)%dom     = seg(j)%dom
  seg(itemp)%anglevis= -un
  seg(itemp)%dom = seg(j)%dom

  ! mise a jours de J
  seg(j)%voiso = itemp
  seg(j)%vnno  = i

!====================================================
! When J is a pinning point, only I must be updated
!====================================================
elseif(j == nsegmax) then

  ve = seg(i)%voise
  rotules = 0

  do while (ve /= j )

    rotules = rotules + 1
    itemp = ve
    ve = seg(ve)%voise

    if(kkdebug) write(379,*)  "Segment ", itemp, ":  will be remove 2"

    if (seg(itemp)%norme /= 0 .or. seg(itemp)%jonc) then
      write(379,*) " error-2 connecij :KK=",kk," REF = ", ref, " remove = ", itemp, ":  non kneecap segment"
      seg(itemp)%norme = Izero
    endif

    if (seg(itemp)%gd == ideux) then
      memoGD2 = memoGD2 + iun         ! A seg gd is eliminated, we want to remember this info
      gd2elimin = gd2elimin + 1       ! The info about eliminated gd2 is keept for multi connecij call
      if (kkdebug) write(379,*) '>>>>>2',gd2elimin
    endif
    if (seg(itemp)%gd == iun) memoGD1 = memoGD1 + iun

!     if (kkdebug) write(379,*) '->3 ',i,j,memoGD1,memoGD2,oldntgd1,oldntgd2

    IDEP(itemp) = Izero
    seg(itemp)%voise = itemp
    seg(itemp)%voiso = itemp
    out(itemp) = .true.
    seg(itemp)%resdep = zero
#ifdef MDC
    seg(itemp)%SIGFE(:) = zero
#endif

    if (seg(itemp)%surface > IZERO .and. seg(itemp)%surface < ITROIS ) then
      seg(i)%surface=seg(itemp)%surface
      seg(i)%VarFreeplan=seg(itemp)%VarFreePlan
      if (seg(itemp)%surface==IUN) seg(i)%dom=seg(itemp)%dom !A segment touching a surface with origin has the domain of the surface
      seg(itemp)%surface=IZERO
      seg(itemp)%VarFreePlan=IZERO
    endif

    if(seg(Itemp)%jonc) then
      lgd1 = seg(Itemp)%Ijonc
      seg(lgd1)%jonc = .false. ;  seg(lgd1)%Ijonc = nsegmax
      if(seg(lgd1)%norme == 0) call voisinage(lgd1,345)
!          call seginfo(lgd1, " verification de connection pour le binome ")
      seg(itemp)%jonc = .false. ;   seg(itemp)%Ijonc = NSEGMAX
      write(379,*) "!> Connecij: 2: we remove the junction ",itemp,": and its pair: ", lgd1,":"
      write(379,*) "KK=",KK, "            REF = ",REF
    endif

    seg(itemp)%gd = izero
    idep(itemp)   = 0

    if (rotules > 100 ) then
      write(*,*) " REF = ", ref
      call conf(i)
      call conf(j)
      write(*,*) "STOP kk",kk, " Infinite loop 2 in connecij"
      stop
    endif

  enddo

  ! We add a kneecap between segment I and NSEGMAX
  PlusSeg = PlusSeg +1 ; itemp = nsegm + PlusSeg
  if(KKdebug) write(379,*)  "We add a kneecap between segment I and NSEGMAX :",itemp

  seg(itemp)%veclin  = nbrot(IVLI,IVLI,2)
  seg(itemp)%norme   = 0
  seg(itemp)%O(:)    = seg(I)%O(:)+seg(i)%norme*bveclin(:,IVLI)
  seg(itemp)%voise   = j
  seg(itemp)%voiso   = i
  seg(itemp)%vnne    = j
  seg(itemp)%vnno    = i
  seg(itemp)%jonc    = .false.
  seg(itemp)%Ijonc   = NSEGMAX
  seg(itemp)%tjonc   = 0
  seg(itemp)%gd      = izero
  seg(itemp)%wait    = IZERO
  SEG(Itemp)%bloquer = seg(i)%bloquer
  SEG(Itemp)%unload  = seg(i)%unload
  seg(itemp)%grain   = seg(i)%grain
  seg(itemp)%dom     = seg(i)%dom
  seg(itemp)%anglevis= -un

  ! mise a jours de I
  seg(I)%voise = itemp
  seg(I)%vnne = j
  call assign_domain(itemp,88)

!====================================================
! cas ou ni I ni J ne sont des point d'encrage :  ...
!====================================================
else

  ! traitement du cas general : ni i ni j est un point d'encrage
  Oi(:) = seg(i)%O(:)
  Ei(:) = Oi(:) + seg(I)%norme * bveclin(:,IVLi)
  X(:) = seg(j)%O(1:3)- EI(1:3)
  X(1) = abs(modulo(X(1),Modur(1)))
  X(2) = abs(modulo(X(2),Modur(2)))
  X(3) = abs(modulo(X(3),Modur(3)))

  if ((x(1) + x(2) + x(3)) /= izero) then
    print *, " REF = ", ref
    print * ," Oj-Ei =", seg(J)%O(:) - Ei(:)
    print *, " connecij: i and j cannot be connected",i,j
    call conf(i)
    call conf(j)
    stop
  endif

  !================================================================================
  ! supression des rotules entre i et J
  !================================================================================
  !==== Cote E de I          ======================================================

  ve = seg(i)%voise
  rotules = 0

  do while (ve /= j .and. ve /= nsegmax)

    rotules = rotules + 1
    itemp = ve
    ve = seg(ve)%voise

    if(kkdebug) write(379,*)  "Segment ", itemp, ":  will be remove 3"

    if (seg(itemp)%norme /= 0) then
      write(379,*) " error-3 connecij :KK=",kk," REF = ", ref, " remove = ", itemp, ":  non kneecap segment"
      seg(itemp)%norme = Izero
    endif

    if (seg(itemp)%gd == ideux) then
      memoGD2 = memoGD2 + iun       ! A seg gd is eliminated, we want to remember this info
      gd2elimin = gd2elimin + 1     ! The info about eliminated gd2 is keept for multi connecij call
      if (kkdebug) write(379,*) '>>>>>3',gd2elimin
    endif
    if (seg(itemp)%gd == iun) memoGD1 = memoGD1 + iun

!     if (kkdebug) write(379,*) '->4 ',i,j,memoGD1,memoGD2,oldntgd1,oldntgd2

    IDEP(itemp)         = Izero
    seg(itemp)%voise    = itemp
    seg(itemp)%voiso    = itemp
    out(itemp)          = .true.
    seg(itemp)%resdep   = Zero
#ifdef MDC
    seg(itemp)%SIGFE(:) = zero
#endif
    seg(itemp)%wait     = IZERO
    seg(Itemp)%bloquer  = seg(j)%bloquer
    seg(Itemp)%unload   = seg(j)%unload
    seg(itemp)%grain    = seg(j)%grain
    seg(itemp)%dom      = seg(j)%dom
    seg(itemp)%NPhase   = izero

    if(seg(Itemp)%jonc) then
      lgd1 = seg(Itemp)%Ijonc
      seg(lgd1)%jonc = .false.
      seg(lgd1)%Ijonc = nsegmax
      if(seg(lgd1)%norme == 0) call voisinage(lgd1,347)
!          call seginfo(lgd1, " verification de connection pour le binome ")
      seg(itemp)%jonc = .false.
      seg(itemp)%Ijonc = NSEGMAX
      if(kkdebug) then
        write(379,*) "!> Connecij: 3: we remove the junction ",itemp,": and its pair: ", lgd1,":"
        write(379,*)  "KK=",KK, "            REF = ",REF
      endif
    endif

    ! seg(itemp)%gd = izero   ! the GD info must not be initialized here since the eliminate segment can
                              ! still be called with a second connecij to reconstruct connectivity
    idep(itemp)   = 0

    if (rotules > 100 ) then
      print *, " REF = ", ref
      call conf(i)
      call conf(j)
      call seginfo(i, "Infinite loop 3 ")
      write(*,*) "STOP kk",kk, " Infinite loop 3 in connecij"
      stop
    endif

  enddo

  !         ======================================================
  ! si on sort da la boucle sur nsegmax alors il faut supprimer les rotules Cote O de J
  if (ve == nsegmax) then
    vo = seg(J)%voiso
    rotules = 0

    do while (vo /= I .and. vo /= nsegmax)

      rotules = rotules + 1
      itemp = vo
      vo = seg(vo)%voiso

      if(kkdebug) write(379,*)  "Segment ", itemp, ":  will be remove 4"

      if (seg(itemp)%gd == ideux) then
        memoGD2 = memoGD2 + iun         ! A seg gd is eliminated, we want to remember this info
        gd2elimin = gd2elimin + 1       ! The info about eliminated gd2 is keept for multi connecij call
      if (kkdebug) write(379,*) '>>>>>4',gd2elimin
      endif
      if (seg(itemp)%gd == iun) memoGD1 = memoGD1 + iun

!       if (kkdebug) write(379,*) '->5 ',i,j,memoGD1,memoGD2,oldntgd1,oldntgd2

      IDEP(itemp) = Izero
      seg(itemp)%norme  = izero
      seg(itemp)%voise  = itemp
      seg(itemp)%voiso  = itemp
      out(itemp)        = .true.
      seg(itemp)%resdep = zero
#ifdef MDC
      seg(I)%SIGFE(:)   = zero
#endif
      seg(itemp)%wait   = IZERO
      SEG(Itemp)%bloquer= seg(j)%bloquer
      SEG(Itemp)%unload = seg(j)%unload
      seg(itemp)%grain  = seg(j)%grain
      seg(itemp)%dom    = seg(j)%dom
      seg(itemp)%NPhase = izero

      if(seg(Itemp)%jonc) then
        lgd1 = seg(Itemp)%Ijonc
        seg(lgd1)%jonc  = .false.
        seg(lgd1)%Ijonc = nsegmax
        if(seg(lgd1)%norme == 0) call voisinage(lgd1,348)
!             call seginfo(lgd1, " verification de connection pour le binome ")
        seg(itemp)%jonc  = .false.
        seg(itemp)%Ijonc = NSEGMAX
        write(379,*) "!> Connecij: 4: we remove the junction ",itemp,": and its pair: ", lgd1,":"
        write(379,*)  "KK=",KK, "            REF = ",REF
      endif

      seg(itemp)%gd = izero
      idep(itemp)   = 0

      if (rotules > 100 ) then
        print *, " REF = ", ref
        call conf(i)
        call conf(j)
        write(*,*) "STOP kk",kk, " Infinite loop 33 in connecij"
        if(kkdebug) stop
      endif

    enddo

  endif

  glisdev = (sysconnec(ivli, ivlj) == 2)

  ! Attention. important:
  ! in the case where I and J are not in the same slip system, we introduce a kneecap (rotule) signed CS
  ! old procedure : if I or J are screw, it becomes CS

  ! pour utiliser le tableau nbrot, il faut mettre ivli et ivlj dans le meme system
  ! le seg GD etant immobile, on converti sont IVL dans l'autre system sans toucher
  ! au type derive seg()%veclin. C'est pour utiliser nbrot
  insert = .false.

  if (glisdev) then

    if (seg(I)%gd > izero) then

      if (seg(i)%norme /= izero .and. syseg(IVLI) == syseg(seg(seg(i)%voiso)%veclin)) then
        if (seg(i)%gd == iun) MemoGD2 = iun  ! This GD will be replaced a we must impose the same gd character
        if(KKdebug) write(379,*)  " Number of segments to add = ", rotules
        seg(i)%gd = izero
        insert = .true.
      else
        ivli = segdevIJ(ivli,ivlj)
        if (ivli == 0) then
          print *, "STOP kk",kk, " connecij : Erreur 1: segdevIJ ??"
          stop
        endif
      endif

    elseif (seg(J)%gd > izero) then

      if (seg(J)%norme /= izero .and. syseg(IVLJ) == syseg(seg(seg(J)%voise)%veclin)) then
        if (seg(j)%gd == iun) MemoGD2 = iun  ! This GD will be replaced a we must impose the same gd character
        seg(j)%gd = izero
        insert = .true.
      else
        ivlJ = segdevIJ(ivlJ,ivli)
        if (ivlJ == 0) then
          print *, "STOP kk",kk, " connecij : Erreur 2: segdevIJ ??"
          stop
        endif
      endif

    ! dans le cas ou les segments i et j ne sont pas du meme system et aucun des deux
    ! ne peut entre considere comme GD, il y a un traitement particulier asser lourd
    ! il s'agit d'introduire une rotule GD entre I et J et de faire la connectivite
    ! entre I-GD et entre GD-J (cas ou INSERT est true)
    else

      insert = .true.

    endif

  endif

  !================================================================================
  ! cas ou I et J sont du meme systems de glissement, on verifie la circulation
  !================================================================================
  if (.not. insert) then

    seg(i)%vnne = j
    seg(j)%vnno = i

    ! le nombre de rotule necessaire entre I et J
    rotules = nbrot(IVLi,IVLJ,1)
    if(KKdebug) write(379,*)  " Number of segments to add = ", rotules

    ! traitement specifique selon le cas
    select case (rotules)

        ! a ) si les segemnts sont directement connectables
      case(IZERO)
        ! mise a jours de i
        seg(i)%voise = j
        ! mise a jours des voisins de i et de j
        seg(j)%voiso = i
        ! b ) une rotule doit etre introduite

      case(IUN)
        ! d'abord la rotule
        PlusSeg = PlusSeg +1
        itemp = nsegm + PlusSeg
        if(KKdebug) write(379,*)  "CASE IUN - We add segment : ",itemp,":"
        seg(itemp)%veclin = nbrot(IVLi,IVLJ,2)
        !      if (seg(itemp)%veclin == 0)
        seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voise = j
        seg(itemp)%voiso = i
        seg(itemp)%vnne = j
        seg(itemp)%vnno = i
        out(itemp) = .false.

        !Domain attribution in case of surfaces
        if (seg(j)%surface > IUN) then
          ! J has the domain of the surface and it could be different that itemp domain
          call assign_domain(itemp,101)
        else
          seg(itemp)%dom=seg(j)%dom
        endif

        ! mise a jours de i
        seg(i)%voise = itemp
        ! mise a jours des voisins de i et de j
        seg(j)%voiso = itemp

      case(IDEUX)
        ! c ) deux rotules doivent etre introduites

        ! la 1ere rotule cote i
        PlusSeg = PlusSeg +2
        itemp = nsegm + PlusSeg -1
        if(KKdebug) write(379,*)  "CASE IDEUX - We add 1st kneecap segment : ",itemp,":"
        seg(itemp)%veclin = nbrot(IVLi,IVLJ,2)
        seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voise = itemp + 1
        seg(itemp)%voiso = i
        seg(itemp)%vnne = j
        seg(itemp)%vnno = i

        !Domain attribution in case of surfaces
        if (seg(j)%surface > IUN ) then
          ! J has the domain of the surface and it could be different that itemp domain
          call assign_domain(itemp,202)
        else
          seg(itemp)%dom=seg(j)%dom
        endif

        out(itemp) = .false.

        ! la 2eme rotule (cote j)
        itemp = nsegm + PlusSeg;
        if(KKdebug) write(379,*)  "CASE IDEUX - We add 2nd kneecap segment : ",itemp,":"
        seg(itemp)%veclin = nbrot(IVLi,IVLJ,3)
        seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = itemp -1
        seg(itemp)%voise = j
        seg(itemp)%vnne = j
        seg(itemp)%vnno = i
        seg(itemp)%dom= seg(nsegm + PlusSeg -1)%dom
        out(itemp) = .false.

        ! mise a jours de i
        seg(i)%voise = itemp-1
        ! mise a jours des voisins de j
        seg(j)%voiso = itemp
        ! si plus de deux rotule : erreur initial d'appele au programme connecIJ

      case(ITROIS)
        ! c ) trois rotules : on le tolere car inevitable
        ! la 1ere rotule cote i
        if(kkdebug) write(379,*)  " REF = ", ref," kk =",kk," Opposed segment connection - 1" &
          & ," IVLI",IVLI,"IVLJ",IVLJ,"IVLIINI",IVLIINI,"IVLJINI",IVLJINI
        PlusSeg = PlusSeg + 3
        itemp = nsegm + PlusSeg - 2
        if(KKdebug) write(379,*)  "CASE ITROIS - We add 1st kneecap segment : ",itemp,":"
        seg(itemp)%veclin = nbrot(IVLi,IVLJ,2)
        seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = i
        seg(itemp)%vnno = i
        seg(itemp)%voise = itemp + 1
        seg(itemp)%vnne = j

        !Domain attribution in case of surfaces
        if (seg(j)%surface > IUN ) then
          ! J has the domain of the surface and it could be different that itemp domain
          call assign_domain(itemp,303)
        else
          seg(itemp)%dom=seg(j)%dom
        endif

        out(itemp) = .false.


        ! 2eme rotule cote i
        itemp = nsegm + PlusSeg - 1
        if(KKdebug) write(379,*)  "CASE ITROIS - We add 2nd kneecap segment : ",itemp,":"
        seg(itemp)%veclin = nbrot(IVLi,IVLJ,3)
        seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = itemp - 1
        seg(itemp)%vnno = i
        seg(itemp)%voise = itemp + 1
        seg(itemp)%vnne = j
        seg(itemp)%dom=seg(nsegm + PlusSeg - 2)%dom
        out(itemp) = .false.


        ! la 3eme rotule (cote j)
        itemp = nsegm + PlusSeg;
        if(KKdebug) write(379,*)  "CASE ITROIS - We add 3rd kneecap segment : ",itemp,":"
        seg(itemp)%veclin = nbrot(IVLi,IVLJ,4)
        seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = itemp - 1
        seg(itemp)%vnno = i
        seg(itemp)%voise = J
        seg(itemp)%vnne = j
        seg(itemp)%dom=seg(nsegm + PlusSeg - 2)%dom
        out(itemp) = .false.

        ! mise a jours de i
        seg(i)%voise = itemp - 2

        ! mise a jours des voisins de j
        seg(j)%voiso = itemp
        ! si plus de deux rotule : erreur initial d'appele au programme connecIJ

    ENDSELECT


    ! We check the case when I and J are declare cross-sliping (GD=True) and its
    ! neighbors are in the same slip system then we want to remove the GD Flag
    ! as its not needed anymore
    if (seg(i)%gd > izero .and. (seg(i)%voiso /=i .and. seg(i)%voise /= i)) then       ! If I is  cross-sliping then

      if (kkdebug) call seginfo(i,"enter  i (1) in connecij")
      vo = seg(i)%voiso       ! recover I first neighbour at Origin
      ve = seg(i)%voise       ! recover I first neighbour at End

      if(syseg(seg(vo)%veclin) == syseg(seg(ve)%veclin)) then ! In that case, i should not cross-slip

        ! Before the gd character is removed, we remember the previous existence of two GD
        if ((memoGD1 + memoGD2) > ideux) then
          if (kkdebug) write(379,*) 'In connecij - 1, oldntgd before icrement =',  &
                                     oldntgd1,oldntgd2,memoGD1,memoGD2,' ref =',ref
          ! memoGD1 and memoGD2 are intialized because in some very particular cases the same test may
          ! be accepted in the next lines for the J segment. In such cases we don t want to increment
          ! again oldntgd1 and oldntgd2
          if (memoGD1 == ideux ) then
            oldntgd1 = oldntgd1 + 2
            memoGD1 = memoGD1 - ideux
          elseif (memoGD2 == itrois) then
            oldntgd2 = oldntgd2 + 1
            memoGD2 = memoGD2 - ideux
          endif

          if (kkdebug) write(379,*) 'In connecij - 1, oldntgd after icrement =',  &
                                     oldntgd1,oldntgd2,memoGD1,memoGD2,' ref =',ref
        endif

        seg(i)%gd     = izero
        seg(i)%veclin = ivli

        if (seg(i)%norme == izero) then

          ! As a GD is removed and seg(i)%norme=0, the vnno and vnne definition must be reconstructed around
          rvnno = seg(i)%vnno                 ! recover I first non null neighbor at Origin
          rvnne = seg(i)%vnne                 ! recover I first non null neighbor at End

          !Treatment on both side to handle pinning point problems
          !Treatment on O side
          itemp = seg(i)%voiso            ! Next segment
          rotules = izero
          do while (itemp /= rvnno)
            rotules = rotules + iun
            seg(itemp)%vnno = rvnno
            seg(itemp)%vnne = rvnne
            itemp = seg(itemp)%voiso
            if (rotules > 100) then
              write(*,*) 'A problem of infinite loop in connecij Origin-side when removing GD on I at step =',KK
              call seginfo(itemp,"A problem of infinite loop in connecij Origin-side when removing GD on I at step")
              call seginfo(i,"A problem of infinite loop in connecij Origin-side when removing GD on I at step")
              stop
            endif
          enddo

          !Treatment on E side
          itemp = seg(i)%voise            ! Next segment
          rotules = izero
          do while (itemp /= rvnne)
            rotules = rotules + iun
            seg(itemp)%vnno = rvnno
            seg(itemp)%vnne = rvnne
            itemp = seg(itemp)%voise
            if (rotules > 100) then
              write(*,*) 'A problem of infinite loop in connecij End-side when removing GD on I at step',KK
              call seginfo(itemp,"A problem of infinite loop in connecij End-side when removing GD on I at step")
              call seginfo(i,"A problem of infinite loop in connecij End-side when removing GD on I at step")
              stop
            endif
          enddo

          if (rvnno /= Nsegmax) seg(rvnno)%vnne = rvnne             ! Last segment
          if (rvnne /= Nsegmax) seg(rvnne)%vnno = rvnno             ! Last segment

          if (kkdebug) call seginfo(i,'Check after GD destruction')

        endif

        if (kkdebug) write(379,fmt='("> No more Cross-slip after connection &
          & , GD is set to .false. for i",I9,":" )') i

      endif

    endif

    if (seg(j)%gd > izero .and. (seg(j)%voiso /=j .and. seg(j)%voise /= j)) then       ! If J is  cross-sliping then

      if (kkdebug) call seginfo(j,"enter  j (2) in connecij")

      vo = seg(j)%voiso       ! recover J first neighbor at Origin
      ve = seg(j)%voise       ! recover J first neighbor at End

      if(syseg(seg(vo)%veclin) == syseg(seg(ve)%veclin)) then ! In that case, j should not cross-lip

        ! Before the gd character is removed, we remember the previous existence of two GD
        if ((memoGD1 + memoGD2) == itrois) then
          if (kkdebug) write(379,*) 'In connecij - 2, oldntgd before icrement =',  &
                                     oldntgd1,oldntgd2,memoGD1,memoGD2,' ref =',ref
          if (memoGD1 == ideux ) then
            oldntgd1 = oldntgd1 + 2
            memoGD1 = memoGD1 - ideux
          elseif (memoGD2 == itrois) then
            oldntgd2 = oldntgd2 + 2
            memoGD2 = memoGD2 - ideux
          endif

          if (kkdebug) write(379,*) 'In connecij - 2, oldntgd after icrement =',  &
                                     oldntgd1,oldntgd2,memoGD1,memoGD2,' ref =',ref
        endif

        seg(j)%gd     = izero
        seg(j)%veclin = ivlj

        if (seg(j)%norme == izero) then
          ! As a GD is removed and seg(j)%norme=0, the vnno and vnne definition must be reconstructed around
          rvnno = seg(j)%vnno                 ! recover I first non null neighbor at Origin
          rvnne = seg(j)%vnne                 ! recover I first non null neighbor at End

          !Treatment on both side to handle pinning point problems
          !Treatment on O side
          itemp = seg(j)%voiso            ! Next segment
          rotules = izero
          do while (itemp /= rvnno)
            rotules = rotules + iun
            seg(itemp)%vnno = rvnno
            seg(itemp)%vnne = rvnne
            itemp = seg(itemp)%voiso
            if (rotules > 100) then
              write(*,*) 'A problem of infinite loop in connecij Origin-side when removing GD on J at step =',KK
              call seginfo(itemp,"A problem of infinite loop in connecij Origin-side when removing GD on J at step")
              call seginfo(j,"A problem of infinite loop in connecij Origin-side when removing GD on J at step")
              stop
            endif
          enddo

          !Treatment on E side
          itemp = seg(j)%voise            ! Next segment
          rotules = izero
          do while (itemp /= rvnne)
            rotules = rotules + iun
            seg(itemp)%vnno = rvnno
            seg(itemp)%vnne = rvnne
            itemp = seg(itemp)%voise
            if (rotules > 100) then
              write(*,*) 'A problem of infinite loop in connecij End-side when removing GD on J at step',KK
              call seginfo(itemp,"A problem of infinite loop in connecij End-side when removing GD on J at step")
              call seginfo(j,"A problem of infinite loop in connecij End-side when removing GD on J at step")
              stop
            endif
          enddo

          if (rvnno /= Nsegmax) seg(rvnno)%vnne = rvnne             ! Last segment
          if (rvnne /= Nsegmax) seg(rvnne)%vnno = rvnno             ! Last segment

          if (kkdebug) call seginfo(j,'Check after GD destruction')

        endif

        if (kkdebug) write(379,fmt='("> No more Cross-slip after connection &
          & , GD is set to .false. for j",I9,":" )') j

      endif

    endif

    !================================================================================
    ! Case where I and J do not belong to the same system, we need to introduce a kneecap with cross-slip information
    !================================================================================
  else
    !================================================================================
    ! Looking for linear vector for cross-slip segment which minimize kneecaps between
    ! I-GD and GD-J
    !================================================================================
    lgd1 = assoc(ivli,1)
    lgd2 = assoc(ivli,5)

    if (nbrot(IVLi,LGD1,1) > 2 .or. nbrot(segdevIJ(LGD1,ivlj),ivlj,1) > 2) then
        lgd = lgd2
        !         if(nbrot(IVLi,LGD,1) > 2 .or. nbrot(segdevIJ(LGD,ivlj),ivlj,1) > 2) then
        !            print *, " REF = ", ref," erreur : connecIJ : connec"
        !            call conf(i) ; call conf(J)
        !            print *, "STOP kk",kk, " CONNECIJ : erreur de connection 1"
        !         endif
    elseif (nbrot(IVLi,LGD2,1) > 2 .or. nbrot(segdevIJ(LGD2,ivlj),ivlj,1) > 2) then
        lgd = lgd1
        !         if(nbrot(IVLi,LGD,1) > 2 .or. nbrot(segdevIJ(LGD,ivlj),ivlj,1) > 2) then
        !            print *, " REF = ", ref," erreur : connecIJ : connec"
        !            call conf(i) ; call conf(J)
        !            print *, "STOP kk",kk, " CONNECIJ : erreur de connection 2"
        !         endif
    else
        rotules = nbrot(IVLi,LGD1,1) + nbrot(segdevIJ(LGD1,ivlj),ivlj,1)
        itemp = nbrot(IVLi,LGD2,1) + nbrot(segdevIJ(LGD2,ivlj),ivlj,1)
        if(rotules < itemp) then
          lgd = lgd1
        else
          lgd = lgd2
        endif
    endif

    !================================================================================
    ! We add a GD kneecap between I and J
    !================================================================================
    PlusSeg = PlusSeg + 1
    seggd = nsegm + PlusSeg

    if (GD2made .and. memoGD2 > iun) memoGD2 = memoGD2 - iun  ! A GD2 was made in the previous connecij call
    if (memoGD2 > ideux) memoGD2 = ideux                      ! One more GD2 segments must be introduced

    seg(segGD)%gd           = memoGD2
    seg(segGD)%veclin       = LGD
    seg(segGD)%norme        = 0
    seg(segGD)%O(:)         = seg(J)%O(:)
    seg(segGD)%vnne         = j
    seg(segGD)%vnno         = i
    seg(segGD)%wait         = IZERO
    seg(segGD)%surface      = IZERO
    seg(segGD)%VarFreePlan  = IZERO
    seg(segGD)%bloquer      = seg(i)%bloquer
    seg(segGD)%unload       = seg(i)%unload
    seg(segGD)%grain        = seg(i)%grain
    seg(segGD)%dom          = seg(j)%dom
    seg(segGD)%anglevis     = -un

    ! Updating I and J
    seg(i)%vnne             = segGD
    seg(j)%vnno             = segGD

    !================================================================================
    ! First : Connection between I and SEGGD
    !================================================================================
    rotules = nbrot(IVLi,LGD,1)

    if(KKdebug) then
      write(379,*)  " The GD2made value before GD introduction = ",GD2made
      write(379,*)  " Cross-slip segment (GD) inserted:", seggd,": with LGD =",LGD," the memoGD2 value = ", memoGD2
      write(379,*)  " Number of kneecap to add between I and GD segment : ", rotules
    endif

    if (memoGD2 > iun) GD2made = .true.  ! A GD2 is made during this connecij call, we keep the info

    if(KKdebug) then
      write(379,*)  " The GD status GD2made, oldntgd1, oldntgd2 =",GD2made, oldntgd1, oldntgd2
    endif

    ! Select case
    select case (rotules)

    case(IZERO)         ! Segments are directly connectable

      ! Updating I
      seg(i)%voise = seggd

      ! Updating SEGGD neighbors
      seg(seggd)%voiso = i

    case(IUN)         ! We need a kneecap to connect segments

      ! First the kneecap
      PlusSeg            = PlusSeg + iun
      itemp              = nsegm + PlusSeg
      if(KKdebug) write(379,*)  "Kneecap added :",itemp,":"
      seg(itemp)%veclin  = nbrot(IVLi,LGD,2)
      seg(itemp)%norme   = 0
      seg(itemp)%O(:)    = seg(SEGGD)%O(:)
      seg(itemp)%dom     = seg(SEGGD)%dom
      seg(itemp)%voise   = seggd
      seg(itemp)%voiso   = i
      seg(itemp)%vnne    = seggd
      seg(itemp)%vnno    = i

      ! Updating I
      seg(i)%voise       = itemp

      ! Updating SEGGD neighbors
      seg(seggd)%voiso   = itemp

    case(IDEUX)         ! We need TWO kneecap to connect segments

      PlusSeg           = PlusSeg + 2

      ! First kneecap
      itemp             = nsegm + PlusSeg -1
      if(KKdebug) write(379,*)  "First kneecap added :",itemp,":"
      seg(itemp)%veclin = nbrot(IVLi,LGD,2)
      seg(itemp)%norme  = 0
      seg(itemp)%O(:)   = seg(SEGGD)%O(:)
      seg(itemp)%dom    = seg(SEGGD)%dom
      seg(itemp)%voise  = itemp + 1
      seg(itemp)%voiso  = i
      seg(itemp)%vnne   = seggd
      seg(itemp)%vnno   = i

      ! Second kneecap
      itemp             = nsegm + PlusSeg
      if(KKdebug) write(379,*)  "Second kneecap added :",itemp,":"
      seg(itemp)%veclin = nbrot(IVLi,LGD,3)
      seg(itemp)%norme  = 0
      seg(itemp)%O(:)   = seg(SEGGD)%O(:)
      seg(itemp)%dom    = seg(SEGGD)%dom
      seg(itemp)%voiso  = itemp -1
      seg(itemp)%voise  = seggd
      seg(itemp)%vnne   = seggd
      seg(itemp)%vnno   = i

      ! Updating I
      seg(i)%voise = itemp-1

      ! Updating SEGGD neighbors
      seg(seggd)%voiso = itemp

    case(ITROIS)         ! We need THREE kneecaps

      if (kkdebug) write(379,*)  " REF = ", ref," kk =",kk," Opposite segment connection"
      PlusSeg = PlusSeg + 3

      ! First kneecap ( I side)
      itemp             = nsegm + PlusSeg - 2
      if(KKdebug) write(379,*)  "First kneecap added :",itemp,":"
      seg(itemp)%veclin = nbrot(IVLi,LGD,2)
      seg(itemp)%norme  = 0
      seg(itemp)%O(:)   = seg(SEGGD)%O(:)
      seg(itemp)%dom    = seg(SEGGD)%dom
      seg(itemp)%voiso  = i
      seg(itemp)%vnno   = i
      seg(itemp)%voise  = itemp + 1
      seg(itemp)%vnne   = SEGGD

      ! Second kneecap ( I side)
      itemp             = nsegm + PlusSeg - 1
      if(KKdebug) write(379,*)  "Second kneecap added :",itemp,":"
      seg(itemp)%veclin = nbrot(IVLi,LGD,3)
      seg(itemp)%norme  = 0
      seg(itemp)%O(:)   = seg(SEGGD)%O(:)
      seg(itemp)%dom    = seg(SEGGD)%dom
      seg(itemp)%voiso  = itemp - 1
      seg(itemp)%vnno   = i
      seg(itemp)%voise  = itemp + 1
      seg(itemp)%vnne   = SEGGD


      ! Third kneecap (SEGGD side)
      itemp             = nsegm + PlusSeg;
      if(KKdebug) write(379,*)  "Third kneecap added :",itemp,":"
      seg(itemp)%veclin = nbrot(IVLi,LGD,4)
      seg(itemp)%norme  = 0
      seg(itemp)%O(:)   = seg(SEGGD)%O(:)
      seg(itemp)%dom    = seg(SEGGD)%dom
      seg(itemp)%voiso  = itemp - 1
      seg(itemp)%vnno   = i
      seg(itemp)%voise  = SEGGD
      seg(itemp)%vnne   = SEGGD

      ! Updating I
      seg(i)%voise      = itemp - 2

      ! Updationg SEGGD neighbor
      seg(SEGGD)%voiso  = itemp

    ENDSELECT

    !================================================================================
    ! Ensuite connection entre SEGGD et J
    !================================================================================
    ! le seg GD etant immobile, on converti sont IVL dans le system de J sans toucher
    ! au type derive seg()%veclin. C'est pour utiliser nbrot
    LGD = SEGDEVIJ(LGD,IVLJ)
    rotules = nbrot(LGD,IVLJ,1)
    if(KKdebug) then
        write(379,*)  "Number of kneecap segments to add between GD segment and J = ", rotules
    endif

    ! traitement specifique selon le cas
    select case (rotules)
        ! a ) si les segemnts sont directement connectables
    case(IZERO)
        ! mise a jours de segGD
        seg(segGD)%voise = J
        ! mise a jours des voisins de i et de j
        seg(J)%voiso = segGD
        ! b ) une rotule doit etre introduite
    case(IUN)
        ! d'abord la rotule
        PlusSeg = PlusSeg +1 ; itemp = nsegm + PlusSeg
        if(KKdebug) write(379,*)  "Kneecap added :",itemp,":"
        seg(itemp)%veclin = nbrot(LGD,IVLJ,2)
        seg(itemp)%norme = 0 ; seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voise = J ;    seg(itemp)%voiso = segGD
        seg(itemp)%vnne = J ;      seg(itemp)%vnno = segGD
        seg(itemp)%dom = seg(J)%dom
        ! mise a jours de segGD
        seg(segGD)%voise = itemp
        ! mise a jours des voisins de segGD et de j
        seg(J)%voiso = itemp
    case(IDEUX)
        ! c ) deux rotules doivent etre introduites

        ! la 1ere rotule cote segGD
        PlusSeg = PlusSeg + 2 ;      itemp = nsegm + PlusSeg -1
        seg(segGD)%voise = itemp
        if(KKdebug) write(379,*)  "First kneecap added :",itemp,":"
        seg(itemp)%veclin = nbrot(LGD,IVLj,2) ; seg(itemp)%norme = 0
        seg(itemp)%O(:) = seg(J)%O(:) ;      seg(itemp)%voise = itemp + 1 ;
        seg(itemp)%voiso = segGD;      seg(itemp)%vnne = J
        seg(itemp)%vnno = segGD
        seg(itemp)%dom = seg(J)%dom

        ! la 2eme rotule (cote J)
        itemp = nsegm + PlusSeg
        if(KKdebug) write(379,*)  "Second kneecap added :",itemp,":"
        seg(itemp)%veclin = nbrot(LGD,IVLJ,3)
        seg(itemp)%norme = 0;      seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = itemp -1;      seg(itemp)%voise = J
        seg(itemp)%vnne = J ;      seg(itemp)%vnno = segGD
        seg(itemp)%dom = seg(J)%dom
        ! mise a jours de segGD
        ! mise a jours des voisins de J
        seg(J)%voiso = itemp
    case(ITROIS)
        ! c ) trois rotules : on le tolere car inevitable
        ! la 1ere rotule cote i
        if (kkdebug) write(379,*)  " REF = ", ref," kk =",kk," opposite segment connection - 3"
        PlusSeg = PlusSeg + 3;
        itemp = nsegm + PlusSeg - 2
        if(KKdebug) write(379,*)  "First kneecap added :",itemp,":"
        seg(itemp)%veclin = nbrot(LGD,IVLJ,2)
        seg(itemp)%norme = 0 ;  seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = SEGGD ;  seg(itemp)%vnno = SEGGD
        seg(itemp)%voise = itemp + 1 ; seg(itemp)%vnne = j
        seg(itemp)%dom = seg(J)%dom

        ! 2eme rotule cote i
        itemp = nsegm + PlusSeg - 1
        if(KKdebug) write(379,*)  "Second kneecap added :",itemp,":"
        seg(itemp)%veclin = nbrot(LGD,IVLJ,3)
        seg(itemp)%norme = 0 ;  seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = itemp - 1 ; seg(itemp)%vnno = SEGGD
        seg(itemp)%voise = itemp + 1 ; seg(itemp)%vnne = j
        seg(itemp)%dom = seg(J)%dom

        ! la 3eme rotule (cote j)
        itemp = nsegm + PlusSeg;
        if(KKdebug) write(379,*)  "Third kneecap added :",itemp,":"
        seg(itemp)%veclin = nbrot(LGD,IVLJ,4)
        seg(itemp)%norme = 0 ;  seg(itemp)%O(:) = seg(J)%O(:)
        seg(itemp)%voiso = itemp - 1 ; seg(itemp)%vnno = SEGGD
        seg(itemp)%voise = J         ; seg(itemp)%vnne = j
        seg(itemp)%dom = seg(J)%dom

        ! mise a jours de i
        seg(SEGGD)%voise = itemp - 2
        ! mise a jours des voisins de j
        seg(j)%voiso = itemp
        ! si plus de deux rotule : erreur initial d'appele au programme connecIJ

    ENDSELECT
  endif
endif

! Remise des rotules ajoutees dans les boite correspondantes
Do rotules = plussegini + 1, plusseg

  Itemp = nsegm + rotules

  if (allocated(nsegboite)) then

    Oi(:) = seg(itemp)%O(:)
    Oi(:) = modulo(Oi(:),modur(:))
    X(1) = modulo(int(real(Oi(1))/tailleBoite(1)),NboitesX) + IUN
    X(2) = modulo(int(real(Oi(2))/tailleBoite(2)),NboitesY) + IUN
    X(3) = modulo(int(real(Oi(3))/tailleBoite(3)),NboitesZ) + IUN
    if(X(1) > NboitesX .or. X(2)> NboitesY .or. X(3) > NboitesZ) then
      print*,"error in subroutine CONNECIJ, one segment coordinates could not be defined in the multipoles doamins"
      print*,"its wrong domain coordinates", X(1:3)
      stop
    endif
    LGD = B3D_1D(X(1),X(2),X(3))

    ! For some simulation it may append that the value of ListeSegRange (defined
    ! in 01Constant.f90) is too small. A test is introduced here to detect such
    ! problems
    if (size(IndexBoite(LGD)%ListeSeg) == NsegBoite(LGD)) then
      write(379,*) 'Problem in domain :', X(1),X(2),X(3),NsegBoite(LGD)
      write(379,*) 'size(IndexBoite(LGD)%ListeSeg) =',LGD,size(IndexBoite(LGD)%ListeSeg)
      write(379,*) 'itemp =', itemp,nsegm,plussegini,plusseg,rotules
      write(379,*) 'Problem %listeseg at step =',kk, '      with Nboites = ',Nboites
      do itemp=1,Nboites
        write(379,*) itemp,NsegBoite(itemp),size(IndexBoite(itemp)%ListeSeg)
      enddo
      stop 'The ListeSegRangeIni value must be increased for this simulation! See the debug file for more info.'
    endif

    NsegBoite(LGD) =  NsegBoite(LGD) + IUN
    IndexBoite(LGD)%ListeSeg(NsegBoite(LGD)) = Itemp
    IBoite(Itemp) = LGD

  endif

enddo


end subroutine connecij

!###################################################
!# subroutine desbou(i) : destruction d'une boucle #
!###################################################
subroutine desbou(i)

implicit none

integer (kind=DPI) :: i,j
integer (kind=DPI) :: NDESTRU,DESTRU
integer (kind=DPI) :: nGDEliminate

if(kkdebug) write(379,*)  " Beginning of subroutine desbou for segment i  = ",i

! Initialization
ndestru = seg(i)%voise
destru  = ndestru
nGDEliminate = izero

! Loop on segment part of a loop to eliminate
do while (destru.ne.i)

  destru = ndestru

  if (ndestru == seg(destru)%voise) then
    write(379,*) 'exit from the desbou loop in a strange manner at step = ',KK
    exit
  endif

  ndestru = seg(destru)%voise

  if(kkdebug) write(379,*) "In desbou we eliminate seg:",destru

  ! if the loop contains gd segments, stat on the gd eliminated is made
  if (seg(destru)%gd == iun) then
    oldntgd1  = oldntgd1 + 1
    ntgd1     = ntgd1    - 1
    seg(destru)%gd  = izero
    nGDEliminate = nGDEliminate + Iun
  endif
  if (seg(destru)%gd == ideux) then
    oldntgd2  = oldntgd2 + 1
    ntgd2     = ntgd2    - 1
    seg(destru)%gd  = izero
    nGDEliminate = nGDEliminate + Ideux
    if (kkdebug) write(379,*) 'In desbou - 1, a GD elimination is made on segment =', destru,ntgd2,oldntgd2
  endif
  if (nGDEliminate == itrois) then
    ! Two GD segments with different GD values and created by connecij was eliminated,
    ! the oldntgd stat must be corrected
    oldntgd1  = oldntgd1 - 1
    oldntgd2  = oldntgd2 - 1
    nGDEliminate = izero
    if (kkdebug) write(379,*) 'In desbou - 1, a double GD elimination was made :', ntgd2,oldntgd2
  endif

  ! Part of the loop is a junction
  if (seg(destru)%jonc) then

    j = seg(destru)%ijonc

    ! We free the complementary segment of the junction
    if (seg(j)%norme /= 0) then
      seg(j)%JONC = .false.       !simple case
      seg(j)%IJONC = nsegmax
      seg(j)%tJONC = 0
    else
      if (seg(j)%vnno == seg(destru)%vnne .or. seg(j)%vnne == seg(destru)%vnno) then
        seg(j)%JONC = .false.       !simple case
        seg(j)%IJONC = nsegmax
        seg(j)%tJONC = 0
      else
        call stopjonc(j,26306)            !non trivial case the vnn must be redifined
      endif
    endif

    seg(destru)%jonc  = .false.
    seg(destru)%ijonc = nsegmax
    seg(destru)%tjonc = izero

  endif

  ! Segment destru is eliminated
  seg(destru)%voiso  = destru
  seg(destru)%voise  = destru
  out(destru)        = .true.
  seg(destru)%norme  = izero

enddo

!write (*,*) kk,i,"desbou"
!call vernsm("ici")	!*** PAS AVEC LES IBM

end subroutine desbou

!#################################################
!# subroutine vire(i) : destruction d'une source #
!#################################################
subroutine vire(i)

implicit none

INTEGER (DPI) :: i
INTEGER (DPI) :: VOIS1E,VOIS1o,NDESTRU,DESTRU

if(kkdebug) write(379,*)  " appel a  VIRE  pour   i = ",i


VOIS1e = SEG(i)%VOISe !*** Premiers voisins nuls ou pas
VOIS1o = SEG(i)%VOISo !*** Premiers voisins nuls ou pas

!*** supression des rotules
ndestru=vois1o
destru=vois1o
do while (destru /= nsegmax .and. destru /= i)
    ndestru = seg(destru)%voiso
    !write (*,*) "desbou ",destru
    seg(destru)%voiso=destru
    seg(destru)%voise=destru
    seg(destru)%norme=0
    destru=ndestru
enddo
if (destru /= i)then
   ndestru=vois1e
   destru=vois1e
   do while (destru /= nsegmax)
       ndestru=seg(destru)%voise
       !write (*,*) "desbou ",destru
       seg(destru)%voiso=destru
       seg(destru)%voise=destru
       seg(destru)%norme=0
       destru=ndestru
   enddo
endif
seg(i)%voiso=destru
seg(i)%voise=destru
seg(i)%norme=0

end subroutine vire

!#################################################

!###############################################################################
!# Procedure Corriger_config: Solve problems of connectivity on segs at O or E #
!###############################################################################
subroutine corriger_config

implicit none

integer(kind=DPI) ::  vli           !<
integer(kind=DPI) ::  vlj           !<
integer(kind=DPI) ::  compt         !<
integer(kind=DPI) ::  I             !<
integer(kind=DPI) ::  SEGI          !<
!integer(kind=DPI) ::  IVNNE         !<
integer(kind=DPI) ::  segIvnne      !<
integer(kind=DPI) ::  vnne          !<
integer(kind=DPI) ::  BO            !<
integer(kind=DPI) ::  BE            !<
integer(kind=DPI) ::  ve            !<
integer(kind=DPI) ::  vo            !<
integer(kind=DPI) ::  vnno          !<
integer(kind=DPI) ::  binome        !<
integer(kind=DPI) ::  numbvnne      !<
integer(kind=DPI) ::  long          !<
integer(kind=DPI) ::  L1            !<
integer(kind=DPI) ::  L2            !<
integer(kind=DPI) ::  tmp           !<
integer(kind=DPI) ::  CumLength     !<
integer(kind=DPI) ::  MINPARALENGTH !<
integer(kind=DPI) ::  voistest      !<
integer(kind=DPI) ::  nextvoise     !< tmp var for %voise
integer(kind=DPI) ::  itemp         !< tmp var for segi
integer(kind=DPI) ::  compteur      !< loop conter
integer(kind=DPI) ::  numbgd1       !< number of gd1 visited
integer(kind=DPI) ::  numbgd2       !< number of gd2 visited

integer(kind=DPI), DIMENSION(10) ::  memo !<
integer(kind=DPI), DIMENSION(3) ::  V     !<

!Initialization
PlusSeg = izero

if (KKdebug) then
  write(379,*)  " !!!! We start the loop Jonc in corriger_config"
!   call check_connec(" Before loop JONC .................")
endif

! loop JONC is used to eliminate asymmetrical junctions and other strange jonction configurations
JONC: Do segI = 1, nsegm

  ! A very strange test, why this update is made here ??? (BD)
  if (seg(SEGI)%probadev /= zero .and. tyseg(seg(SEGI)%veclin) /= iun) then
    print*,'strange update of %probadev in corriger_config at step =',KK,'   Segi =',segi
    seg(SEGI)%probadev = zero
  endif

  if (.not. seg(segI)%jonc) cycle JONC                                                    ! not a jonc segments
  if ((SegI == seg(SegI)%voiso .and. SegI == seg(SegI)%voise) .or. out(SegI)) Cycle JONC  ! a segment to eliminate

  ! Initialization
  GD2made = .false.
  binome = seg(SEGI)%Ijonc

  if (.not. seg(BINOME)%jonc) then

    ! A paranoid part, we should never enter this part.
    print*, "KK=",KK," corriger config: junction removal because of asymmetry between ",segI, binome

    seg(binome)%jonc = .false.
    seg(binome)%Ijonc = nsegmax

    ! in the case where binome was a kneecap, its neigbhor should be reconnected a gain
    if (seg(binome)%norme == 0) call connecij(seg(binome)%vnno,seg(binome)%vnne,11345)

    seg(segI)%jonc = .false.
    seg(segI)%Ijonc = NSEGMAX

    ! in the case where binome was a kneecap, its neigbhor should be reconnected a gain
    if (seg(segI)%norme == 0) call connecij(seg(segI)%vnno,seg(segI)%vnne,12345)

  else

    !elimination of two segments in collinear 'junction' separated by a half loop with a GD
    if (((seg(seg(segi)%vnno)%GD * seg(seg(binome)%vnne)%GD) > izero .and.        &
         seg(segi)%vnno == seg(binome)%vnne) .or.                                 &
        ((seg(seg(segi)%vnne)%GD * seg(seg(binome)%vnno)%GD) > izero .and.        &
         seg(segi)%vnne == seg(binome)%vnno)) then

      ! The particular case of two segments in collinear 'junction' closed on themself with gd segments at both ends
      if (seg(segi)%vnno == seg(binome)%vnne .and. seg(segi)%vnne == seg(binome)%vnno) then

        ! Determination of hands of old junction for the cutting of segments : iO side
        nextvoise = seg(segi)%voise
        compteur  = izero

        do while (nextvoise /= segi .and. compteur < 100)
          itemp             = nextvoise
          nextvoise         = seg(itemp)%voise
          out(itemp)        = .true.
          seg(itemp)%voiso  = itemp
          seg(itemp)%voise  = itemp
          seg(itemp)%jonc   = .false.
          seg(itemp)%Ijonc  = NSEGMAX
          seg(itemp)%tjonc  = IZERO
          compteur  = compteur + iun

          ! Info about the eliminated GD segments must be saved
          if (seg(itemp)%gd == iun)   oldntGD1 = oldntGD1 + 1
          if (seg(itemp)%gd == ideux) oldntGD2 = oldntGD2 + 1

          if (kkdebug) write(379,*) " The segment ", itemp,                          &
             ": is eliminated in loop Jonc in corriger_config. oldntGD2 =", oldntGD2
        enddo

        if (compteur > 99) then
          call seginfo(segi,'Problem in corrier_config and the Jonc loop')
          call seginfo(binome,'Problem in corrier_config and the Jonc loop')
          stop 'Problem in corrier_config and the Jonc loop'
        endif

        out(segi)       = .true.
        seg(segi)%voiso = segi
        seg(segi)%voise = segi
        seg(segi)%jonc  = .false.
        seg(segi)%Ijonc = NSEGMAX
        seg(segi)%tjonc = IZERO

        ! Info about the eliminated GD segments must be saved
        if (seg(segi)%gd == iun)   oldntGD1 = oldntGD1 + 1
        if (seg(segi)%gd == ideux) oldntGD2 = oldntGD2 + 1

        if (kkdebug) write(379,*) " The segment ", segi,                          &
           ": is eliminated in loop Jonc in corriger_config. oldntGD2 =", oldntGD2

      else

        seg(segI)%jonc = .false.
        seg(segI)%Ijonc = NSEGMAX
        seg(segI)%tjonc = IZERO

        if (seg(segI)%norme == 0)   call connecij(seg(segI)%vnno,seg(segI)%vnne,14345)

        seg(binome)%jonc = .false.
        seg(binome)%Ijonc = NSEGMAX
        seg(binome)%tjonc = IZERO

        if (seg(binome)%norme == 0) call connecij(seg(binome)%vnno,seg(binome)%vnne,15345)

      endif

    endif

  endif

enddo JONC

nsegm = nsegm + PLUSSEG

if (nsegm > nsegmax) then
  write(*,*) 'Problem in (5) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif

PlusSeg = izero
tmp = nsegmax

if (KKdebug) then
  write(379,*)  " !!!! We start the loop NETGD in corriger_config"
  !   call check_connec(" Before loop NETGD .................")
endif

! the loop NETGD is used to eliminate unnecessary GD characters
NETGD: do SegI = 1 , Nsegm

  if (seg(SegI)%gd < iun) cycle NETGD                                                         ! SegI is not a GD
  if ((SegI == seg(SegI)%voiso .and. SegI == seg(SegI)%voise) .or. out(SegI)) Cycle  NETGD    ! SegI is to eliminate

  ! Initialization
  GD2made = .false.
  Bo  = seg(SegI)%vnno
  Be  = seg(SegI)%vnne

  ! if SegI is in between two pining points, no need to be a GD
  if (Bo == nsegmax .and. Be == nsegmax) then

    seg(SegI)%gd = izero
    if(seg(SegI)%norme /= izero) then
      call connecIJ(Bo,SegI,144)
      GD2made = .false.
      call connecIJ(SegI,Be,145)
    else
      ! this case is crazy, SegI is to eliminate
      seg(SegI)%voiso = SEGI
      seg(SegI)%voise = SEGI
    endif

  ! si I est voisin de point d'encrage ou d'un GD: annuler GD et connecter directement son vnn
  elseif (Bo == nsegmax .or. seg(Bo)%gd > izero) then

    if(seg(SegI)%norme /= izero) then
      call connecIJ(SegI,Be,147)
      GD2made = .false.
      call connecIJ(Bo,SegI,146)
    else
      call connecIJ(Bo,Be,1467)
    endif

    seg(SegI)%gd = izero

  ! si I est voisin de point d'encrage: annuler GD et connecter directement son vnn
  elseif (Be == nsegmax .or. seg(Be)%gd > izero) then

    if(seg(SegI)%norme /= izero) then
      call connecIJ(Bo,SegI,148)
      GD2made = .false.
      call connecIJ(SegI,Be,149)
    else
      call connecIJ(Bo,Be,1489)
    endif

    seg(SegI)%gd = izero

  ! cette condition sert a liberer le seg I GD des qu'il n'est plus rotule
  ! la subtilite est dans le choix de la rotule GD et son emplacement:
  ! 1) on attribue a I le veclin du systeme le plus charge
  ! 2) on introduit la rotule cote system le moin actif
  else

    VLI =  seg(SegI)%veclin

    if (Bo /= nsegmax .and. Be /= nsegmax) then
      L1 = syseg(seg(Bo)%veclin)
      L2 = syseg(seg(Be)%veclin)
    else
      ! The exceptional case of GD segments ending at a surface with a kneecap as only neighbor
      L1 = syseg(seg(seg(SegI)%voiso)%veclin)
      L2 = syseg(seg(seg(SegI)%voise)%veclin)
    endif

    if (L1 /= L2) then

      if (seg(SegI)%norme /= izero) then
        ! ici on sait que I est gd non nul et que ses voisins ne sont pas points d'encrage
        ! definition of the loading factor for Bo and vnne

        if (DesorientGrain) then
          vo = int(solli_sys(L1) *  SchmidSysGrain(seg(SegI)%grain,L1),DPI)
          ve = int(solli_sys(L2) *  SchmidSysGrain(seg(SegI)%grain,L2),DPI)
        else
          vo = int(solli_sys(L1) *  SchmidSys(L1),DPI)
          ve = int(solli_sys(L2) *  SchmidSys(L2),DPI)
        endif

        ! cas ou le system cote O est plus sollicite que le system cote E: on introduit une rotule GD, cote Se si il est vis
        if (vo > ve) then
          ! on ajoute une rotule entre I et son Vnne
          !  Bo -- I -- i1 -- Be
          if(KKdebug) write(379,*)  " Transfer du caractere GD de: ",SEGI,"  a la rotule entre I et Vnne"
          VLI = segdevIJ(VLI,seg(Bo)%veclin)
          !      print *, " Creation de rotule GD cote E. Le VLI de ancien GD est devevu ", seg(I)%veclin
        elseif (vo < ve) then
          ! cas ou le Be est plus sollicite que Bo, on passe le caractere GD a Bo, si il est vis
          ! on ajoute une rotule entre Bo et I
          !  Bo -- i1 -- I -- Be
          if(KKdebug) write(379,*)  " Transfer du caractere GD de: ",SEGI,"  a la rotule entre Vnno et I"
          VLI = segdevIJ(VLI,seg(Be)%veclin)
        !else : vo==ve :  Loading on slip systems is the same and thus we do not change GD Veclin
          !      print *, " Creation de rotule GD cote O. Le VLI de ancien GD est devevu ", seg(I)%veclin
        endif

        ! determination si il faut ajouter une rotule gd de chaque cote de I
        SEG(SEGI)%VECLIN = VLI
        L1 = seg(seg(SegI)%voiso)%veclin
        L2 = seg(seg(SegI)%voise)%veclin

        if (syseg(VLI) /= syseg(L1))  call connecIJ(Bo,SegI,76232)
        if (syseg(VLI) /= syseg(L2))  call connecIJ(SegI,Be,7635)

        seg(SegI)%gd = izero    ! SegI is not any more a gd and was replaced

        if (KKdebug) call seginfo(SegI," Transfer du caractere GD cote vnno ou vnne")

      endif

    ! if bnvvo and Be belong to the same system: no need for GD
    else

      seg(SegI)%gd = izero

      if (seg(SegI)%norme /= izero) then
        SEG(SegI)%VECLIN = segdevIJ(Vli,seg(Be)%veclin)
      else
        call connecIJ(BO,BE,151)
      endif

    endif

  endif

enddo  NETGD

NsegM = Nsegm + PlusSeg
if(nsegm > nsegmax) then
  write(*,*) 'Problem in (6) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif
PlusSeg = izero

if (KKdebug) then
  write(379,*)  " !!!! We start the loop side_O in corriger_config"
  call check_connec(" Before loop SIDE_O .................")
endif

!First loop function is to eliminate accumulated kneecaps at O when the side segment is a fix point.
SIDE_O: Do segI = 1, nsegm

  if ((SegI == seg(SegI)%voiso .and. SegI == seg(SegI)%voise) .or. out(SegI)) Cycle SIDE_O  ! a segment to eliminate

  if (seg(segI)%GD > izero .and. SEG(SEGI)%JONC) then               !Only non zero segments are of interest
    PRINT *, "CORRIGER_CONFIG: fatal error: I=",SEGI," is junc and GD, KK=",kk
    stop
  endif

  ! Initialization
  GD2made = .false.

  if(seg(segI)%norme /= izero) then               !Only non zero segments are of interest

    vo = seg(segI)%voiso
    ve = seg(segI)%voise

    if (seg(vo)%norme /= izero) cycle SIDE_O   !Only zero length segments are now of interest
    if (segI == ve .and. segI == vo) cycle SIDE_O    !To be eliminate anyway

    !Who is next
    do compt = 1, 100

      !the following segs are not interesting, we can stop the research loop
      if (seg(vo)%jonc .or. (seg(vo)%norme /= izero .and. vo /= nsegmax)) CYCLE side_O

      if (vo == nsegmax) then                !This seg is a pinning point
        if (compt > IDEUX) then
            EXIT
        else
            Cycle side_O         !There is less than three kneecaps, then
        endif                   !everything is OK
      endif
      vo = seg(vo)%voiso

    enddo

    if (compt > 99) then           !Something stupid is found, better stop
      print *, "corrig_conf: infinite loop 1: I= ",SEGI,"   kk = ",KK
      if (kkdebug)  stop
      cycle  side_O
      !stop
    endif

    !This VNNO is a pinning point, but there is to many kneecaps between I and VNNO
    if (kkdebug) write(379,*)  " KK ",kk," Connection correction between i =",segi," and vnno=NSEGMAX"
    call connecIJ(vo,SegI,2564)     !A correct connection is rebuild

  endif

enddo SIDE_O


if (KKdebug) then
  write(379,*)  " !!!! We start the loop RECOUVREMENT in corriger_config"
endif

!Second loop function is to eliminate possible overlapping segs at the end of a step
RECOUVREMENT: Do segi = 1, nsegm

  ! Initialization
  GD2made = .false.

  ! One treat only segments of finite length
  if (seg(segI)%norme /= izero) then

    if ((SegI == seg(SegI)%voiso .and. SegI == seg(SegI)%voise) .or. out(SegI)) Cycle RECOUVREMENT  ! a segment to eliminate

    ! Initialization
    numbgd1 = izero
    numbgd2 = izero

    long  = seg(segI)%norme
    Vli   = seg(segI)%veclin
    vo    = seg(segI)%voiso
    ve    = seg(segI)%voise

    numbvnne = izero

    !Research of the next true VNNE, without considering intermediate VNNE of GD and JONC type if of zero length
    do compt = 1, 100

      if (ve == nsegmax) cycle RECOUVREMENT          !A pining point, see after

      if (seg(ve)%norme /= izero) then               !ve is the VNNE

        V(:)  = bveclin(:,vli) + bveclin(:,seg(ve)%veclin)

        !We check if i and its VNNE can geometrically overlap, V(:)=zero
        tmp = abs(V(1)) + abs(V(2)) + abs(V(3))
        if(tmp /= izero) then
          cycle RECOUVREMENT   !The segs are not averlaping, it s OK
        else
          exit                 !The segs are overlaping
        endif

      else

        ! During the vnne search if we pass too many GD or JONC segments of zero length
        ! better to wait for simplification before doing the modification
        if (seg(ve)%gd > izero .or. seg(ve)%jonc) then
          numbvnne = numbvnne + iun
          if (numbvnne > iun) cycle RECOUVREMENT
        endif

      endif

      !We check that the neighbour declaration is OK
      if(seg(ve)%voise == ve .or. seg(ve)%voise == SEGI) then
        print *, " fatal erreur : bad neighbour side E for i=",SEGI," kk = ",KK
        stop
      endif

      ve = seg(ve)%voise

    enddo

    if (compt > 99) then           !Something stupid is found, better stop
      print *, "corrig_conf: infinite loop 2: I= ",SEGI,"   kk = ",KK
      if (kkdebug) then
          call seginfo(segI, " corriger_config: erreur ..............")
          stop
      endif
      cycle RECOUVREMENT
    endif

    !Let see who are the overlapping segs
    vnne  = ve  ! last ve found is either gd, junc or non-zero length
    VlJ   = seg(vnne)%veclin

    if (kkdebug) write(379,*)  " KK",kk," Overlapping between i =",segi," and its vnne",vnne

    !We check at the origin of i, if several segs with the same VLD
    !are connected to form a longer straight line section
    BO    = SEGI
    L1    = Long
    V(:)  = seg(BO)%O(:)

    do compt = 1, 100

      BO =  seg(BO)%voiso
      ! We do not reconstruct lines including jonc or GD or touching a free surface
      !if (seg(BO)%gd > izero .or. seg(BO)%jonc .or. seg(BO)%surface>0) cycle RECOUVREMENT
      if (seg(BO)%jonc .or. seg(BO)%surface>0) cycle RECOUVREMENT

      if (BO == vnne) then   !This is a closed flat loop to be deleted in subroutine net
        if (kkdebug) call seginfo(BO,"A small closed loop eliminated on the BO side")
        call DESBOU(BO)
        cycle RECOUVREMENT
      endif

      if (seg(BO)%norme == izero) cycle      !This seg is a kneecap

      if (BO /= nsegmax .and. seg(BO)%veclin == vli) then
        L1 = L1 + seg(BO)%norme            !The total segs length
        V(:) = seg(BO)%O(:)                !The new Origin coordinates
      else
        EXIT   !This seg has a diferent VLD, we can stop the do loop
      endif

    enddo

    if (compt > 99) then          !Something stupid is found, better stop
      if (GB==IZERO) then  ! problem with barriers. to see later
        if (.not.seg(segI)%bloquer) then
          print *, "corrig_conf: boucle infinie 3 kk = ",KK
          if (kkdebug) then
              call seginfo(segI, " corriger_config: erreur ..............")
              stop
          endif
        endif
      endif
      cycle RECOUVREMENT
    endif

    !We check at the extremity of vnne, if several segs with the same VLD
    !are connected to form a longer straight line section
    BE = vnne
    L2 = seg(BE)%norme

    do compt = 1, 100

      BE = seg(BE)%voise

      ! We do not reconstruct lines including jonc or GD or touching a free surface
      if (seg(BE)%gd > izero .or. seg(BE)%jonc .or. seg(BE)%surface>0) cycle RECOUVREMENT

      if (BE == BO) then  !This is a closed flat loop to be deleted in subroutine net
        if (kkdebug) call seginfo(BE,"A small closed loop eliminated on the BE side")
        call DESBOU(BE)
        cycle RECOUVREMENT
      endif

      if (seg(BE)%norme == izero) cycle      !This seg is a kneecap

      if (BE /= nsegmax .and. seg(BE)%veclin == vlJ) then
        L2 = L2 + seg(BE)%norme            !The total segs length
      else
        EXIT      !This seg has a diferent VLD, we can stop the do loop
      endif

    enddo

    if (compt > 99) then          !Something stupid is found, better stop
      if (GB==IZERO) then  ! problem with barriers. to see later
        if (.not.seg(vnne)%bloquer) then
          print *, "corrig_conf: boucle infinie 4 kk = ",KK
          if (kkdebug) then
              call seginfo(segI, " corriger_config: erreur ..............")
              stop
          endif
        endif
      endif

      if (kkdebug)  stop
      cycle RECOUVREMENT  ! la il s'agit d'une microboucle a supprimer dans net
    endif

    !All the segs between BO and BE are now eliminated
    !reconnection is made latter
    if (BO /= nsegmax) then

      tmp = seg(BO)%voise

      do compt = 1,200

        if (tmp == BE) exit
        if(kkdebug) write(379,*)  "  !> This segment will be removed  (1) : ",tmp,": with gd =",seg(tmp)%gd
        ve = seg(tmp)%voise

        !The seg tmp is eliminated
        seg(tmp)%norme = izero
        seg(tmp)%voiso = tmp
        seg(tmp)%voise = tmp
        out(tmp)       = .true.
        if (seg(tmp)%GD == iun)   numbgd1 = numbgd1 + iun
        if (seg(tmp)%GD == ideux) numbgd2 = numbgd2 + iun
        seg(tmp)%GD    = izero

        !If tmp is jonc, we must liberate its binome
        if(seg(tmp)%jonc)then
          if(kkdebug) write(379,*)  "corriger_config: juntion removal 10,KK",kk
          binome = seg(tmp)%Ijonc
          seg(binome)%Ijonc = nsegmax
          seg(binome)%jonc  =.false.

          if(seg(binome)%norme == izero)  call connecIJ(seg(binome)%vnno,seg(binome)%vnne,2812)

          seg(tmp)%Ijonc = nsegmax
          seg(tmp)%jonc =.false.
        endif

        tmp = ve           !Next one

      enddo

      if (compt > 199) then       !Something stupid is found, better stop
        print *, "corrig_conf: boucle infinie 5 kk = ",KK
        if (kkdebug) then
          call seginfo(segI, " corriger_config: erreur ..............")
          stop
        endif
        if (kkdebug)  stop
        cycle RECOUVREMENT
      endif

    else

      tmp = seg(BE)%voiso

      do compt = 1,200

        if (tmp == BO) exit
        if(kkdebug) write(379,*)  "  !> This segment will be removed  (2) : ",tmp,": with gd =",seg(tmp)%gd
        vo = seg(tmp)%voiso

        !The seg tmp is eliminated
        seg(tmp)%norme = izero
        seg(tmp)%voiso = tmp
        seg(tmp)%voise = tmp
        out(tmp)       = .true.
        if (seg(tmp)%GD == iun)   numbgd1 = numbgd1 + iun
        if (seg(tmp)%GD == ideux) numbgd2 = numbgd2 + iun
        seg(tmp)%gd    = izero

        !If tmp is jonc, we must liberate its binome
        if(seg(tmp)%jonc)then
          if(kkdebug) write(379,*)  " corriger_config: juntion removal 11,kk",kk
          binome = seg(tmp)%Ijonc
          seg(binome)%Ijonc = nsegmax
          seg(binome)%jonc  =.false.

          if(seg(binome)%norme == izero)  call connecIJ(seg(binome)%vnno,seg(binome)%vnne,2812)

          seg(tmp)%Ijonc = nsegmax
          seg(tmp)%jonc =.false.
        endif

        !Transfer of surface informations
        if (seg(tmp)%surface > IZERO)then
          seg(SEGI)%surface      = seg(tmp)%surface
          seg(SEGI)%VarFreePlan  = seg(tmp)%VarFreePlan
        endif

        tmp = vo            !Next one

      enddo

      if (compt > 199) then      !Something stupid is found, better stop
        print *, "corrig_conf: infinite loop 6 kk = ",KK
        if (kkdebug) then
          call seginfo(segI, " corriger_config: error ..............")
          stop
        endif
        if (kkdebug)  stop
        cycle RECOUVREMENT
      endif

    endif

    ! During the segs elimination between BO and BE we may have eliminated two GD seg
    ! These GD cannot be reconstructructed
    if (numbgd1 == ideux) oldntGD1 = oldntGD1 + 2
    if (numbgd2 == ideux) oldntGD2 = oldntGD2 + 2
    if (numbgd1 == iun .and. numbgd2 == iun) then
      oldntGD1 = oldntGD1 + 1
      oldntGD2 = oldntGD2 + 1
    endif

    if (kkdebug) then
      write(379,*) " Avant reconstructionnumbgd1, numbgd2, oldntGD1, oldntGD2 = ", numbgd1, numbgd2, oldntGD1, oldntGD2
      write(379,*) "Overlapping between L1= ", L1," et L2=",L2, " KK =",kk
      write(379,*) "BO =",BO, ": BE =",BE,":"
    endif

    !Now the reconstruction, we start with the case when L1<L2
    if (L1 < L2) then

      !Segs are added to connect BO and BE
      plusseg = plusseg + 1
      tmp = nsegm + plusseg
      if (kkdebug) write(379,fmt='("L1 < L2, we add the segment",I7)') tmp

      seg(tmp)%O(:)       = V(:)
      seg(tmp)%veclin     = VLJ
      seg(tmp)%norme      = L2 - L1
      seg(tmp)%voiso      = BO
      seg(tmp)%voise      = BE
      seg(tmp)%grain      = seg(SEGI)%grain
      seg(tmp)%wait       = IZERO
      seg(tmp)%bloquer    = seg(SEGI)%bloquer
      seg(tmp)%unload     = seg(SEGI)%unload
      seg(tmp)%anglevis   = seg(SEGI)%anglevis
      seg(tmp)%dom        = seg(SEGI)%dom
      seg(tmp)%surface    = seg(SEGI)%surface
      seg(tmp)%VarFreePlan  = seg(SEGI)%VarFreePlan

      if(BO /= nsegmax) seg(BO)%voise = tmp
      if(BE /= nsegmax) seg(BE)%voiso = tmp

      if(seg(vnne)%surface /= IZERO) then
        seg(tmp)%surface      = seg(vnne)%surface
        seg(tmp)%VarFreePlan  = seg(vnne)%VarFreePlan
        if (seg(vnne)%surface==IUN) seg(tmp)%dom = seg(vnne)%dom
        seg(vnne)%surface     = IZERO
        seg(vnne)%VarFreePlan = IZERO
      endif

      ! if a gd2 segment was eliminated, we need to recover the information after the connecij between vnno and vnne
      if (numbgd2 == ideux) then
        GD2made = .false.
        call connecIJ(BO,tmp, -66966)
        if (GD2made) oldntgd2 = oldntgd2 - 1
        GD2made = .false.
        call connecIJ(tmp,BE, -66967)
        if (GD2made) oldntgd2 = oldntgd2 - 1
      elseif (numbgd2 == iun) then
        if (BO == NSEGMAX .or. BE == NSEGMAX) then
          GD2made = .false.
          call connecIJ(BO,tmp, 34961)
          GD2made = .false.
          call connecIJ(tmp,BE, 34962)
        elseif (syseg(seg(BO)%veclin) /= syseg(seg(tmp)%veclin)) then
          GD2made = .false.
          call connecIJ(BO,tmp, -55966)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
          GD2made = .false.
          call connecIJ(tmp,BE, 55967)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
        elseif (syseg(seg(tmp)%veclin) /= syseg(seg(BE)%veclin)) then
          GD2made = .false.
          call connecIJ(BO,tmp, 44966)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
          GD2made = .false.
          call connecIJ(tmp,BE, -44967)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
        endif
      else
        GD2made = .false.
        call connecIJ(BO,tmp, 33966)
        GD2made = .false.
        call connecIJ(tmp,BE, 33967)
      endif

    elseif (L1 > L2) then

      !Segs are added to connect BO and BE
      plusseg = plusseg + 1
      tmp = nsegm + plusseg
      if (kkdebug) write(379,fmt='("L1 > L2, we add the segment",I7)') tmp

      seg(tmp)%O(:)     = V(:)
      seg(tmp)%veclin   = VLI
      seg(tmp)%norme    = L1 - L2
      seg(tmp)%voiso    = BO
      seg(tmp)%voise    = BE

      !Domain:
      if (L1 > Long) then !meaning that itemp origin may has been displaced due to segment merging on BO side
        call assign_domain(tmp,1112)
      else
        seg(tmp)%dom      = seg(SEGI)%dom
      endif

      if(BO /= nsegmax) seg(BO)%voise = tmp
      if(BE /= nsegmax) seg(BE)%voiso = tmp
      if(seg(Segi)%surface /= IZERO) then
        seg(tmp)%surface      = seg(Segi)%surface
        seg(tmp)%VarFreePlan  = seg(Segi)%VarFreePlan
        seg(Segi)%surface     = IZERO
        seg(Segi)%VarFreePlan = IZERO
      endif

      seg(tmp)%grain          = seg(SEGI)%grain
      if (kkdebug) call seginfo(tmp,'L1 > L2 : After surface decision')
      seg(tmp)%wait           = IZERO
      seg(tmp)%bloquer        = seg(SEGI)%bloquer
      seg(tmp)%unload         = seg(SEGI)%unload
      seg(tmp)%anglevis       = seg(SEGI)%anglevis

      ! if a gd2 segment was eliminated, we need to recover the information after the connecij between vnno and vnne
      if (numbgd2 == ideux) then
        GD2made = .false.
        call connecIJ(BO,tmp, -67968)
        if (GD2made) oldntgd2 = oldntgd2 - 1
        GD2made = .false.
        call connecIJ(tmp,BE, -67969)
        if (GD2made) oldntgd2 = oldntgd2 - 1
      elseif (numbgd2 == iun) then
        if (BO == NSEGMAX .or. BE == NSEGMAX) then
          GD2made = .false.
          call connecIJ(BO,tmp, 34963)
          GD2made = .false.
          call connecIJ(tmp,BE, 34964)
        elseif (syseg(seg(BO)%veclin) /= syseg(seg(tmp)%veclin)) then
          GD2made = .false.
          call connecIJ(BO,tmp, -57968)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
          GD2made = .false.
          call connecIJ(tmp,BE, 57969)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
        elseif (syseg(seg(tmp)%veclin) /= syseg(seg(BE)%veclin)) then
          GD2made = .false.
          call connecIJ(BO,tmp, 47968)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
          GD2made = .false.
          call connecIJ(tmp,BE, -47969)
          if (GD2made .and. numbgd1 > izero) oldntgd2 = oldntgd2 - 1
        endif
      else
        GD2made = .false.
        call connecIJ(BO,tmp, 37968)
        GD2made = .false.
        call connecIJ(tmp,BE, 37969)
      endif

    else

      !The simplest case L1 = L2
      !BO and BE can be directly reconnected
      if(BO /= nsegmax) seg(BO)%voise = BE
      if(BE /= nsegmax) seg(BE)%voiso = BO

      if (numbgd2 /= izero) then
        GD2made = .false.
        call connecIJ(BO,BE, -39699)
      else
        GD2made = .false.
        call connecIJ(BO,BE, 39699)
      endif

    endif

    if (kkdebug) write(379,*) " Apres reconstruction numbgd1, numbgd2, oldntGD1, oldntGD2 = ", &
                              numbgd1, numbgd2, oldntGD1, oldntGD2

  endif

enddo RECOUVREMENT

if (KKdebug) then
  write(379,*)  " !!!! We start the loop side_E in corriger_config"
!   call check_connec(" Before loop SIDE_E .................")
endif

!Third loop function is to eliminate accumulated kneecaps at E side
SIDE_E: Do segi = 1, nsegm

  if ((SegI == seg(SegI)%voiso .and. SegI == seg(SegI)%voise) .or. out(SegI)) Cycle SIDE_E  ! a segment to eliminate

  ! Initialization
  GD2made = .false.
  long = seg(segI)%norme
  tmp = nsegmax

  if(long /= izero .or. seg(segi)%jonc) then        ! Only segs of finite length are of interest
    ve = seg(segI)%voise

    do compt = 1, 100
      if(ve == nsegmax) then            !ve is a pinning point
        if(compt < 3) then
            cycle SIDE_E                !The connectivity is a priori OK
        else
            call connecIJ (SEGI,ve,9874)   !Clean connection is made
            cycle SIDE_E
        endif
      elseif(seg(ve)%norme /= izero .or. seg(ve)%jonc) then
        EXIT
      endif
      ve = seg(ve)%voise

    enddo

    if (compt > 99) then         !Something stupid is found, better stop
      print *, "corrig_conf: boucle infinie 22: I= ",SEGI,"   kk = ",KK
      call seginfo(segi, " corriger_config: erreur ..............")
      if (kkdebug)  stop
      CYCLE side_E
    endif

    vnne = ve
    Vli = seg(segI)%veclin
    VlJ = seg(ve)%veclin

    !case of connectivity between segs of differents slip systemes
    if(syseg(VLI) /= syseg(VlJ)) then
      L1 = segdevIJ(VLJ,VLI)
      ! calculation of the optimal nomber of kneecaps between I and vnne including a CS kneecap.
      if(nbrot(VLi,assoc(VLI,1),1) + nbrot(assoc(VLI,1),L1,1) <         &
         nbrot(VLi,assoc(VLI,5),1) + nbrot(assoc(VLI,5),L1,1)) then
        L2 = nbrot(VLi,assoc(VLI,1),1) + nbrot(assoc(VLI,1),L1,1)
      else
        L2 = nbrot(VLi,assoc(VLI,5),1) + nbrot(assoc(VLI,5),L1,1)
      endif

      if((seg(segI)%gd + seg(vnne)%gd) > izero) then

        if(compt > IUN + nbrot(VLI,L1,1)) then
          if (kkdebug) write(379,*)  "corrig_conf: correction type1 : I= ",SEGI,"   kk = ",KK
          !call seginfo(i, " corriger_config: before correction type 1 ......")
          call connecIJ (SEGI,vnne,9274)
          !call seginfo(segI, " corriger_config: aftere correction type 1 ......")
        endif

      elseif(compt > L2 + IDEUX ) then

        if (kkdebug) write(379,*)  "corrig_conf: correction type2 : I= ",SEGI,"   kk = ",KK
        !call seginfo(segI, " corriger_config: before correction type 2 ......")
        call connecIJ (SEGI,vnne,9275)
        !call seginfo(i, " corriger_config: aftere correction type 2 ......")

      endif

    !case of connectivity between segs of the same slip systeme
    else

      if (compt >  IUN + nbrot(VLI,VLj,1)) then
        if (kkdebug) write(379,*)  "corrig_conf: correction type3 : I= ",SEGI,"   kk = ",KK
        !call seginfo(segI, " corriger_config: before correction type 3 ......")
        call connecIJ (SEGI,vnne,9276)
        !call seginfo(i, " corriger_config: aftere correction type 3 ......")
      endif

    endif

  endif

enddo SIDE_E

Nsegm = Nsegm + PlusSeg
if(nsegm > nsegmax) then
  write(*,*) 'Problem in (7) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif
PlusSeg = izero

if (KKdebug) then
  write(379,*)  " !!!! We start the loop MINIPARA in corriger_config"
!   call check_connec(" Before loop MINIPARA .................")
endif

MINPARALENGTH = IDEUX * MODEP_MAX
memo(:) = izero

! The spacial case of tiny segments in parallel direction and connected with a kneecap
! When possible those segments are reconnected to reduce the total number of segments and
! to help segments displacement on the simulation lattice.
MINPARA: do I = 1,Nsegm

  if (seg(I)%norme == IZERO)                                                    cycle MINPARA ! I must not be a kneecap
  if ((I == seg(I)%voiso .and. I == seg(I)%voise) .or. out(I))                  cycle MINPARA ! a segment to eliminate
  if (seg(I)%norme > MINPARALENGTH)                                             cycle MINPARA ! Are the segments small enough
  if (seg(seg(i)%vnne)%norme > MINPARALENGTH)                                   cycle MINPARA
  if (seg(I)%veclin /= seg(seg(I)%vnne)%veclin)                                 cycle MINPARA ! Are the two segments parallel ?
  if (seg(I)%jonc .or. seg(I)%gd > izero)                                       cycle MINPARA ! There is no jonc or gd in the vicinity
  if (seg(seg(I)%voise)%jonc .or. seg(seg(I)%voise)%gd > izero)                 cycle MINPARA
  if (seg(seg(I)%vnne)%jonc .or. seg(seg(I)%vnne)%gd > izero)                   cycle MINPARA
  if (seg(seg(I)%vnne)%voise == NSEGMAX  .or. seg(seg(I)%vnne)%vnne == NSEGMAX) cycle MINPARA ! No pinning points
  if (seg(I)%surface /= 0 .or. seg(seg(I)%vnne)%surface /= 0)                   cycle MINPARA ! The segments are not touching a surface

  if(kkdebug) then
    write(379,*) "Segment", I, "is recombined with ",seg(i)%vnne, "in the MINPARALENGTH loop!"
    call seginfo(I,"recombination in MINPARALENGTH")
  endif

  ! Initialization
  GD2made = .false.
  compt = izero

  ! The recombination
  segIvnne = seg(I)%vnne
  seg(I)%norme = seg(I)%norme + seg(segIvnne)%norme
  ! Kneecaps between I and segIvnne  + segIvnne are eliminated
  voistest = seg(I)%voise

  do

    compt = compt + iun
    memo(compt) = voistest            ! We keep info about segments to be eliminated
    if (voistest /= segIvnne) then
      !print*,'voistest1=',voistest,segIvnne,out(voistest),compt
      voistest = seg(voistest)%voise
    else
      !print*,'voistest2=',voistest,segIvnne,out(voistest),compt
      exit
    endif
    if (compt > 9) then
      write(*,*) "problem in the MinPara loop of NET !",compt
      call seginfo(segIvnne,"problem in the MinPara loop of NET")
      stop
    endif

  enddo

  ! Connectivity of the seg(I)%vnne segment is rebuild
  seg(seg(segIvnne)%vnne)%vnno = I

  if (seg(segIvnne)%vnne == seg(segIvnne)%voise) then
    ! good for us there is no kneecap
    seg(seg(segIvnne)%vnne)%voiso = I
  else
    ! The kneecap connectivity must be reconstructed
    voistest = seg(seg(segIvnne)%vnne)%voiso

    do

      seg(voistest)%vnno = I
      if (seg(voistest)%voiso /= segIvnne) then
        voistest = seg(voistest)%voiso
      else
        seg(voistest)%voiso = I
        exit
      endif

    enddo

  endif

  ! Connectivity of I is reconstructed
  seg(I)%voise = seg(seg(I)%vnne)%voise
  seg(I)%vnne  = seg(seg(I)%vnne)%vnne

  ! Useless segments are now eliminated
  do segi = 1,compt
    seg(memo(segi))%voiso = memo(segi)
    seg(memo(segi))%voise = memo(segi)
    out(memo(segi))       = .true.
  enddo

enddo MINPARA

if (KKdebug) then
  write(379,*)  " !!!! We start the loop SUPERPOSITION in corriger_config"
!   call check_connec(" Before loop SUPERPOSITION .................")
endif

! Case of segment superposition close to gd
SUPERPOSITION: do i = 1, nsegm

  ! Initialization
  numbgd2 = izero
  plusseg = IZERO

  if (seg(i)%norme == IZERO)                                      cycle
  if ((I == seg(I)%voiso .and. I == seg(I)%voise) .or. OUT(I))    cycle
  vnno = seg(i)%vnno
  vnne = seg(i)%vnne
  if (vnno == Nsegmax .or. vnne == Nsegmax)         cycle   ! This test must not be made if we are close to a pinning point

  ! Initialization
  GD2made = .false.

  ! The superposition test
  if ((seg(vnno)%gd + seg(seg(vnne)%vnne)%gd) > izero                   .and. &
      seg(i)%veclin /= seg(vnne)%veclin .and. seg(vnne)%norme /= IZERO  .and. &
      tyseg(seg(i)%veclin) == tyseg(seg(vnne)%veclin))                               then

    if (kkdebug) then
      write(379,*) 'i',i
      write(379,*) 'vnne', vnne
      call seginfo(i,"superposition")
    endif

    ! If it is a small loop, we do nothing since the segments are going to be eliminated in NET latter
    tmp       = vnne
    CumLength = seg(i)%norme

    do compt = 1, 100
      ! The MicroBoucle test
      if (tmp == i .and. CumLength < MicroBoucle) then
        if (kkdebug) write(379,*) 'A small loop to be eliminate in Net is found in superposition. i = ',i
        cycle SUPERPOSITION
      endif

      tmp = seg(tmp)%vnne
      CumLength = CumLength + seg(tmp)%norme
    enddo

    ! The two superposing segments have not the same length
    if (seg(i)%norme /= seg(seg(i)%vnne)%norme) then

      tmp = seg(i)%voise

      do compt = 1,10

        if(tmp == seg(i)%vnne) then
            exit
        endif
        if(kkdebug) write(379,*)  "KK=",kk,"coriger_config GD: elimination 1 du seg =",tmp
        ve = seg(tmp)%voise
        !The seg tmp is eliminated
        seg(tmp)%norme = izero
        seg(tmp)%voiso = tmp
        seg(tmp)%voise = tmp
        if (seg(tmp)%GD > iun) numbgd2 = numbgd2 + 1
        seg(tmp)%GD    = izero
        if (seg(tmp)%jonc) then
          seg(seg(tmp)%Ijonc)%Ijonc=nsegmax
          seg(seg(tmp)%Ijonc)%jonc=.false.
          if (.not. OUT(seg(seg(tmp)%Ijonc)%voiso) .and. &
              .not. OUT(seg(seg(tmp)%Ijonc)%voise)) call voisinage(seg(tmp)%Ijonc,81)
          seg(tmp)%Ijonc = nsegmax
          seg(tmp)%jonc =.false.
        endif
        OUT(tmp) = .true.
        tmp = ve            !Next one

      enddo

      if (seg(i)%norme > seg(seg(i)%vnne)%norme) then

        tmp=seg(seg(i)%vnne)%voise

        do compt = 1,10

          if(seg(tmp)%norme /= Izero) exit

          if(kkdebug) write(379,*)  "KK=",kk,"coriger_config GD: elimination 2 du seg =",tmp

          ve = seg(tmp)%voise

          !The seg tmp is eliminated
          seg(tmp)%norme = izero
          seg(tmp)%voiso = tmp
          seg(tmp)%voise = tmp
          if (seg(tmp)%GD > iun) numbgd2 = numbgd2 + 1
          seg(tmp)%GD    = izero

          if (seg(tmp)%jonc) then
            seg(seg(tmp)%Ijonc)%Ijonc=nsegmax
            seg(seg(tmp)%Ijonc)%jonc=.false.
            if (.not. OUT(seg(seg(tmp)%Ijonc)%voiso).and. &
                .not. OUT(seg(seg(tmp)%Ijonc)%voise)) call voisinage(seg(tmp)%Ijonc,82)
            seg(tmp)%Ijonc = nsegmax
            seg(tmp)%jonc =.false.
          endif

          OUT(tmp) = .true.
          tmp = ve            !Next one

        enddo

        vnne=tmp

        seg(i)%norme           = seg(i)%norme - seg(seg(i)%vnne)%norme
        seg(seg(i)%vnne)%norme = izero
        seg(seg(i)%vnne)%voiso = seg(i)%vnne
        seg(seg(i)%vnne)%voise = seg(i)%vnne
        seg(seg(i)%vnne)%GD    = izero

        if (seg(seg(i)%vnne)%jonc) then
          seg(seg(seg(i)%vnne)%Ijonc)%Ijonc=nsegmax
          seg(seg(seg(i)%vnne)%Ijonc)%jonc=.false.
          if (.not. OUT(seg(seg(seg(i)%vnne)%Ijonc)%voiso).and. &
              .not. OUT(seg(seg(seg(i)%vnne)%Ijonc)%voise)) call voisinage(seg(seg(i)%vnne)%Ijonc,821)
          seg(seg(i)%vnne)%Ijonc = nsegmax
          seg(seg(i)%vnne)%jonc =.false.
        endif

        OUT(seg(i)%vnne) = .true.
        seg(i)%voise=vnne
        seg(i)%vnne=vnne
        seg(vnne)%voiso=i
        seg(vnne)%vnno=i

        if (seg(i)%Ijonc == vnne) then
          seg(i)%Ijonc = nsegmax
          seg(i)%jonc = .false.
          seg(vnne)%Ijonc = nsegmax
          seg(vnne)%jonc = .false.
        elseif(seg(vnne)%jonc .and. seg(seg(vnne)%ijonc)%ijonc/=vnne) then
          seg(vnne)%Ijonc = nsegmax
          seg(vnne)%jonc = .false.
        elseif(seg(i)%jonc .and. seg(seg(i)%ijonc)%ijonc/= i) then
          seg(i)%Ijonc = nsegmax
          seg(i)%jonc = .false.
        endif

        ! If a gd2 segment was eliminated, we need to recover the information after the connecij between vnno and vnne
        if (numbgd2 /= izero) then
          GD2made = .false.
          call connecIJ (i,vnne,-9273)
        else
          GD2made = .false.
          call connecIJ (i,vnne,9273)
        endif



      else

        tmp=seg(i)%voiso

        do compt = 1,10
          if(seg(tmp)%norme /= Izero ) then
              !seg(tmp)%Ijonc = nsegmax
              !seg(tmp)%jonc =.false.
              exit
          endif
          if(kkdebug) write(379,*)  "KK=",kk,"coriger_config GD: elimination 3 du seg =",tmp
          vo = seg(tmp)%voiso
          !The seg tmp is eliminated
          seg(tmp)%norme = izero
          seg(tmp)%voiso = tmp
          seg(tmp)%voise = tmp
          if (seg(tmp)%GD > iun) numbgd2 = numbgd2 + 1
          seg(tmp)%GD    = izero
          if (seg(tmp)%jonc) then
            seg(seg(tmp)%Ijonc)%Ijonc=nsegmax
            seg(seg(tmp)%Ijonc)%jonc=.false.
            if (.not. OUT(seg(seg(tmp)%Ijonc)%voiso) .and. &
                .not. OUT(seg(seg(tmp)%Ijonc)%voise))call voisinage(seg(tmp)%Ijonc,83)
            seg(tmp)%Ijonc = nsegmax
            seg(tmp)%jonc =.false.
          endif
          OUT(tmp) = .true.
          tmp = vo            !Next one
        enddo
        vnno = tmp


        vnne=seg(i)%vnne
        seg(vnne)%norme= seg(vnne)%norme-seg(i)%norme

        seg(vnne)%O     = seg(i)%O
        seg(vnne)%dom   = seg(i)%dom

        seg(i)%norme = izero
        seg(i)%voiso = i
        seg(i)%voise = i
        seg(i)%GD    = izero

        if (seg(i)%jonc) then
          seg(seg(i)%Ijonc)%Ijonc=nsegmax
          seg(seg(i)%Ijonc)%jonc=.false.
          if (.not. OUT(seg(seg(i)%Ijonc)%voiso) .and. &
              .not. OUT(seg(seg(i)%Ijonc)%voise))call voisinage(seg(i)%Ijonc,831)
            seg(i)%Ijonc = nsegmax
            seg(i)%jonc =.false.
        endif

        OUT(i) = .true.

        seg(vnno)%voise = vnne
        seg(vnno)%vnne  = vnne
        seg(vnne)%voiso = vnno
        seg(vnne)%vnno  = vnno

        if (seg(vnne)%Ijonc == vnno) then
          seg(vnne)%Ijonc = nsegmax
          seg(vnne)%jonc = .false.
          seg(vnno)%Ijonc = nsegmax
          seg(vnno)%jonc = .false.
        elseif(seg(vnne)%jonc .and. seg(seg(vnne)%ijonc)%ijonc/=vnne) then
          seg(vnne)%Ijonc = nsegmax
          seg(vnne)%jonc = .false.
        elseif(seg(vnno)%jonc .and. seg(seg(vnno)%ijonc)%ijonc/= vnno) then
          seg(vnno)%Ijonc = nsegmax
          seg(vnno)%jonc = .false.
        endif

        ! If a gd2 segment was eliminated, we need to recover the information after the connecij between vnno and vnne
        if (numbgd2 /= izero) then
          GD2made = .false.
          call connecIJ (vnno,vnne,-8274)
        else
          GD2made = .false.
          call connecIJ (vnno,vnne,8274)
        endif

      endif

    else

      tmp = seg(i)%voiso

      !eliminating rotule on vnno and i side
      do compt = 1,10

        if(seg(tmp)%norme /= Izero ) exit
        if(kkdebug) write(379,*)  "KK=",kk,"coriger_config GD: elimination 4 du seg =",tmp

        vo = seg(tmp)%voiso

        !The seg tmp is eliminated
        seg(tmp)%norme = izero
        seg(tmp)%voiso = tmp
        seg(tmp)%voise = tmp
        if (seg(tmp)%GD > iun) numbgd2 = numbgd2 + 1
        seg(tmp)%GD    = izero

        if (seg(tmp)%jonc) then
          seg(seg(tmp)%Ijonc)%Ijonc=nsegmax
          seg(seg(tmp)%Ijonc)%jonc=.false.
          if (.not. OUT(seg(seg(tmp)%Ijonc)%voiso) .and. &
              .not. OUT(seg(seg(tmp)%Ijonc)%voise))         call voisinage(seg(tmp)%Ijonc,91)
          seg(tmp)%Ijonc = nsegmax
          seg(tmp)%jonc =.false.
        endif

        OUT(tmp) = .true.
        tmp = vo            !Next one

      enddo

      vnno = tmp

      tmp = seg(i)%voise

      do compt = 1,10

        if(tmp == seg(i)%vnne) exit
        if(kkdebug) write(379,*)  "KK=",kk,"coriger_config GD: elimination 5 du seg =",tmp
        ve = seg(tmp)%voise

        !The seg tmp is eliminated
        seg(tmp)%norme = izero
        seg(tmp)%voiso = tmp
        seg(tmp)%voise = tmp
        if (seg(tmp)%GD > iun) numbgd2 = numbgd2 + 1
        seg(tmp)%GD    = izero

        if (seg(tmp)%jonc) then
          seg(seg(tmp)%Ijonc)%Ijonc=nsegmax
          seg(seg(tmp)%Ijonc)%jonc=.false.
          if (.not. OUT(seg(seg(tmp)%Ijonc)%voiso) .and. &
              .not. OUT(seg(seg(tmp)%Ijonc)%voise))call voisinage(seg(tmp)%Ijonc,92)
          seg(tmp)%Ijonc = nsegmax
          seg(tmp)%jonc =.false.
        endif

        OUT(tmp) = .true.
        tmp = ve            !Next one

      enddo

      tmp=seg(seg(i)%vnne)%voise

      do compt = 1,10

        if(seg(tmp)%norme /= Izero) exit
        if(kkdebug) write(379,*)  "KK=",kk,"coriger_config GD: elimination 6 du seg =",tmp

        ve = seg(tmp)%voise

        !The seg tmp is eliminated
        seg(tmp)%norme = izero
        seg(tmp)%voiso = tmp
        seg(tmp)%voise = tmp
        if (seg(tmp)%GD > iun) numbgd2 = numbgd2 + 1
        seg(tmp)%GD    = izero

        if (seg(tmp)%jonc) then
          seg(seg(tmp)%Ijonc)%Ijonc=nsegmax
          seg(seg(tmp)%Ijonc)%jonc=.false.
          if (.not. OUT(seg(seg(tmp)%Ijonc)%voiso) .and. &
              .not. OUT(seg(seg(tmp)%Ijonc)%voise))call voisinage(seg(tmp)%Ijonc,93)
          seg(tmp)%Ijonc = nsegmax
          seg(tmp)%jonc =.false.
        endif

        OUT(tmp) = .true.
        tmp = ve            !Next one

      enddo

      vnne = tmp

      seg(seg(i)%vnne)%norme = izero
      seg(seg(i)%vnne)%voiso = seg(i)%vnne
      seg(seg(i)%vnne)%voise = seg(i)%vnne
      seg(seg(i)%vnne)%GD    = izero

      if (seg(seg(i)%vnne)%jonc) then
        seg(seg(seg(i)%vnne)%Ijonc)%Ijonc=nsegmax
        seg(seg(seg(i)%vnne)%Ijonc)%jonc=.false.
        if (.not. OUT(seg(seg(seg(i)%vnne)%Ijonc)%voiso) .and.      &
            .not. OUT(seg(seg(seg(i)%vnne)%Ijonc)%voise))             call voisinage(seg(seg(i)%vnne)%Ijonc,94)
        seg(seg(i)%vnne)%Ijonc = nsegmax
        seg(seg(i)%vnne)%jonc =.false.
      endif

      OUT(seg(i)%vnne)= .TRUE.

      seg(i)%norme = izero
      seg(i)%voiso = i
      seg(i)%voise = i
      seg(i)%GD    = izero

      if (seg(i)%jonc) then
        seg(seg(i)%Ijonc)%Ijonc=nsegmax
        seg(seg(i)%Ijonc)%jonc=.false.
        if (.not. OUT(seg(seg(i)%Ijonc)%voiso) .and.    &
            .not. OUT(seg(seg(i)%Ijonc)%voise))           call voisinage(seg(i)%Ijonc,95)
        seg(i)%Ijonc = nsegmax
        seg(i)%jonc =.false.
      endif

      OUT(i)=.TRUE.

      seg(vnno)%vnne  = vnne
      seg(vnno)%voise = vnne
      seg(vnne)%vnno  = vnno
      seg(vnne)%voiso = vnno

      if (seg(vnne)%Ijonc == vnno) then
        seg(vnne)%Ijonc = nsegmax
        seg(vnne)%jonc = .false.
        seg(vnno)%Ijonc = nsegmax
        seg(vnno)%jonc = .false.
      elseif(seg(vnne)%jonc .and. seg(seg(vnne)%ijonc)%ijonc/=vnne) then
        seg(vnne)%Ijonc = nsegmax
        seg(vnne)%jonc = .false.
      elseif(seg(vnno)%jonc .and. seg(seg(vnno)%ijonc)%ijonc/= vnno) then
        seg(vnno)%Ijonc = nsegmax
        seg(vnno)%jonc = .false.
      endif

      ! If a gd2 segment was eliminated, we need to recover the information after the connecij between vnno and vnne
      if (numbgd2 /= izero) then
        GD2made = .false.
        call connecIJ (vnno,vnne,-9275)
      else
        GD2made = .false.
        call connecIJ (vnno,vnne,9275)
      endif

    endif

    nsegm = nsegm + plusseg

    if (kkdebug) call seginfo(vnno,"superposition end")

  endif

enddo SUPERPOSITION

end subroutine corriger_config


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assign_domain(i,ref)

implicit none

integer (kind=DPI) :: OldDom    !< Domain when segment arrive in loop
integer (kind=DPI) :: NewDom    !< Domain finally attributed to the sgment
integer (kind=DPI) :: i   !< Segment number
integer (kind=DPI) :: jj  !< Counter
integer (kind=DPI) :: ii  !< Counter
!integer (kind=DPI) :: surfcount !< Surface counter

integer (kind=DPI) :: plane   !< Store the number of transparent surface and grain boundary surface checked through loops
integer (kind=DPI) :: planeFP !< Store the number of free surface checked through loops

real (kind=DP)     :: Oi(3)       !< Coordinates of the segment Origin
real (kind=DP)     :: DotProdO    !<  Scalar Product between the vector [Oi-Intersect0] and the surface normal

real (kind=DP)     :: normint(3)  !< Surface normal (Miller indices , reals)

logical,dimension(NbCvxDom) :: bloquerD     !< 'True' if the segment origin is NOT in the domain, 'False' otherwise
!logical                    :: IsInDom(2)   !< 'True' if the segment origin is in the domain, 'False' otherwise
integer                     :: ref          !< Case tracker

Oi(:) = real(modulo(seg(i)%O(:),modur(:)),DP)

plane     = IZERO
planeFP   = IZERO
OldDom    = seg(i)%dom
NewDom    = IZERO

bloquerD(:) = .false.

if (kkdebug) then
  write(379,*)  ""
  write(379,fmt='("-- Assign domain ",I6," (Old segment domain : ", I3,") for segment ",I7,":")') ref,OldDom, I
endif

!Segment or kneecap touching a surface

if     ( seg(i)%surface==IUN                                            ) then
  seg(i)%dom = Plane_dom(seg(i)%VarFreePlan)
  if (kkdebug) then
    write(379,fmt='("Seg(i)%surface ",I4,"")') seg(i)%surface
    write(379,fmt='("-- Assign domain n1 ",I6," (New segment domain : ", I3,")")') ref,seg(i)%dom
  endif
  return

elseif ( seg(i)%norme == IZERO .and. seg(seg(I)%voise)%surface == IUN   ) then
  seg(i)%dom = Plane_dom(seg(seg(I)%voise)%VarFreePlan)
  if (kkdebug) then
    write(379,fmt='("seg(I)%voise ",I7,": Surface ",I4)') seg(I)%voise, seg(seg(I)%voise)%surface
    write(379,fmt='("-- Assign domain n2 ",I6," (New segment domain : ", I3,")")') ref,seg(i)%dom
  endif
  return

elseif ( seg(i)%norme == IZERO .and. seg(seg(I)%voiso)%surface == IDEUX ) then
  seg(i)%dom = Plane_dom(seg(seg(I)%voiso)%VarFreePlan)
    if (kkdebug) then
      write(379,fmt='("seg(I)%voiso ",I7,": Surface ",I4)') seg(I)%voiso, seg(seg(I)%voiso)%surface
      write(379,fmt='("-- Assign domain n3 ",I6," (New segment domain : ", I3,")")') ref,seg(i)%dom
    endif
  return

endif

do ii=1,NbCvxDom !Loop over the domains

  if (kkdebug) write(379,fmt='("Looking if segment", I9,  ": belongs to domain ",I4," (Origin : ",3I15," )")') i,ii,int(Oi,DPI)

  !The segment position is first tested with respect to transparent and grain boundary surfaces of domain ii

  if (NbPlanCvxDom(1,ii)==IZERO .and. nbcvxDom > 1) then
    write(*,fmt='("Stop - The Domain",I3, "has no transparent and grain boundary surfaces")') ii
    !stop
  endif

  if (kkdebug) write(379,fmt='(3x,"Looking at closed boundary surface ")')

  do jj= plane+1,plane+NbPlanCvxDom(1,ii) ! Loop over transparent and grain boundary surfaces


    normint(:) = Plane_MillerR(1:3,jj)
    DotProdO=Oi(1)*normint(1)+Oi(2)*normint(2)+Oi(3)*normint(3)-Plane_Pos(jj)

    if (kkdebug) write(379,fmt='(x,"VFP",I3,x,"Miller",3I5,x,"Dotprod", F20.3)')jj,Plane_MillerI(1:3,jj),DotProdO

    if (DotProdO > zero) then

      bloquerD(ii) = .true.
      exit
    endif

  enddo ! Loop End over transparent and grain boundary surfaces

  if (bloquerD(ii)) then

      plane=plane+NbPlanCvxDom(1,ii)
      planeFP=planeFP+NbPlanCvxDom(2,ii)

      cycle

  else

    !If the considered segment end is inside the domain when just looking at grain boundaries surfaces,
    !we check if is always true when looking at the domain free surfaces
    if (NbPlanCvxDom(2,ii)/=IZERO) then ! Check if there is free surfaces for domain ii
      if (kkdebug) write(379,fmt='(3x,"Looking at free surfaces ")')
      do jj= NbPlan+planeFP+1,NbPlan+planeFP+NbPlanCvxDom(2,ii) ! Loop over free surfaces

        ! Normale projection of OiP on the surface
         normint(:) = Plane_MillerR(1:3,jj)
         DotProdO=Oi(1)*normint(1)+Oi(2)*normint(2)+Oi(3)*normint(3)-Plane_Pos(jj)

         if (kkdebug) write(379,fmt='(x,"VFP",I3,x,"Miller",3I5,x,"Dotprod", F20.3)')jj,Plane_MillerI(1:3,jj),DotProdO

         if (DotProdO > zero) then

            bloquerD(ii) = .true.
            exit

         endif

      enddo ! Loop over free surfaces

    else
      if (kkdebug) write(379,fmt='("No Free surface to check")')
    endif
  endif

  if (bloquerD(ii)) then

    if (ii /= NbCvxDom) then

      plane=plane+NbPlanCvxDom(1,ii)
      planeFP=planeFP+NbPlanCvxDom(2,ii)

      cycle
    endif

  else
    seg(i)%dom = ii !Assign New Domain

    if (kkdebug) write(379,fmt='(">> Segment Domain is at the end of this procedure : ",I3)') seg(i)%dom
    if (NbCvxDom > IZERO .and. seg(i)%dom == IZERO) stop 'Assign Domain - Domain == 0 at the end of the subroutine'
   return
  endif

  plane=plane+NbPlanCvxDom(1,ii)
  planeFP=planeFP+NbPlanCvxDom(2,ii)

enddo !End of the loop over domains

end subroutine assign_domain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_domain(i,ref)

implicit none

integer (kind=DPI) :: Domi    !< Domain when segment arrive in loop
integer (kind=DPI) :: i   !< Segment number
integer (kind=DPI) :: jj  !< Counter
integer (kind=DPI) :: ii  !< Counter

integer (kind=DPI) :: plane   !< Store the number of transparent surface and grain boundary surface checked through loops
integer (kind=DPI) :: planeFP !< Store the number of free surface checked through loops

real (kind=DP)     :: DotProdO !<  Scalar Product between the vector [Oi-Intersect0] and the surface normal
real (kind=DP)     :: Oi(3)       !< Coordinates of the segment Origin

real (kind=DP)     :: normint(3)    !< Surface normal (Miller indices , reals)
real (kind=DP)     :: interS(3)

integer            :: ref !< Case tracker

Oi(:) = real(modulo(seg(i)%O(:),modur(:)))

if (seg(i)%surface == IUN) then
  normint(:)    = Plane_MillerR(1:3,seg(i)%VarfreePlan)
  InterS(:)  = InterPlanSeg(normint(:),Plane_pos(seg(i)%VarfreePlan),Oi,Bveclin(:,seg(i)%veclin)) !intersection with free surface after displacement!!
  Oi(:) =InterS
elseif(seg(seg(i)%voiso)%surface == IDEUX) then
  normint(:)    = Plane_MillerR(1:3,seg(seg(i)%voiso)%VarfreePlan)
  InterS(:)  = InterPlanSeg(normint(:),Plane_pos(seg(seg(i)%voiso)%VarfreePlan),Oi,Bveclin(:,seg(seg(i)%voiso)%veclin)) !intersection with free surface after displacement!!
  Oi(:) =InterS
elseif( seg(seg(i)%voise)%surface  == IUN) then
  normint(:)    = Plane_MillerR(1:3,seg(seg(i)%voise)%VarfreePlan)
  InterS(:)  = InterPlanSeg(normint(:),Plane_pos(seg(seg(i)%voise)%VarfreePlan),Oi,Bveclin(:,seg(seg(i)%voise)%veclin)) !intersection with free surface after displacement!!
  Oi(:) =InterS
endif



plane=IZERO
planeFP=IZERO
Domi=seg(i)%dom

if (kkdebug) then
  write(379,*)  ""
  write(379,fmt='("-- Check domain ",I6,A2,3I12,A2,I6)') i, " (",int(Oi,DPI),") ",ref
endif

if (NbCvxDom > IZERO .and. Domi<=IZERO) then
  call seginfo(i,"check domain stop")
  stop 'Check Domain - No Domain for this segment'
endif
do ii=1,NbCvxDom !Loop over the domains

  if (ii/=domi) then

    plane=plane+NbPlanCvxDom(1,ii)     !Update plane for domain ii+1
    planeFP=planeFP+NbPlanCvxDom(2,ii) !Update planeFP for domain ii+1
    if (kkdebug) write(379,fmt='(x,"Skipping domain ",I4)') ii

  else
    if (kkdebug) write(379,fmt='(x,"Testing Domain ",I4)') ii

    if (kkdebug) write(379,fmt='(3x,"Looking at closed boundary surface ")')
    !The segment position is first tested with respect to transparent and grain boundary surfaces of domain ii
    do jj= plane+1,plane+NbPlanCvxDom(1,domi) ! Loop over transparent and grain boundary surfaces

      normint(:) = Plane_MillerR(1:3,jj)
      DotProdO=Oi(1)*normint(1)+Oi(2)*normint(2)+Oi(3)*normint(3)-Plane_Pos(jj)

      if (kkdebug) write(379,fmt='(x,"VFP",I3,x,"Miller",3I5,x,"Dotprod", F20.3)')jj,Plane_MillerI(1:3,jj),DotProdO

      if (DotProdO > numtol_dotp) then

        print *, "Origin ", Oi
        print *, "Plan ", Plane_MillerI(1:3,jj),Plane_pos(jj), seg(i)%varfreeplan
        call seginfo(I,"Segment origin outside FIXED boundary")
        stop 'Segment origin outside FIXED boundary'

      endif

    enddo ! Loop End over transparent and grain boundary surfaces

    !If the considered segment end is inside the domain when just looking at grain boundaries surfaces,
    !we check if is always true when looking at the domain free surfaces
    if (kkdebug) write(379,fmt='(3x,"Looking at free surface ")')

    if (NbPlanCvxDom(2,ii)/=IZERO) then ! Check if there is free surfaces for domain ii
      do jj= NbPlan+planeFP+1,NbPlan+planeFP+NbPlanCvxDom(2,ii) ! Loop over free surfaces

        ! Normale projection of OiP on the surface
        normint(:) = Plane_MillerR(1:3,jj)
        DotProdO=Oi(1)*normint(1)+Oi(2)*normint(2)+Oi(3)*normint(3)-Plane_pos(jj)

        if (kkdebug) write(379,fmt='(x,"VFP",I3,x,"Miller",3I5,x,"Dotprod", F20.3)') jj,Plane_MillerI(1:3,jj),DotProdO

        if (DotProdO > numtol_dotp) then

          print *, "Origin ", Oi
          print *, "Plan ", Plane_MillerI(1:3,jj),Plane_pos(jj), seg(i)%varfreeplan
          call seginfo(I,"Segment origin outside Free boundary")
          stop "Segment origin outside Free boundary"

        endif

      enddo ! Loop over free surfaces

    else
      if (kkdebug) write(379,fmt='("No Free surface to check")')
    endif

  endif
enddo !End of the loop over domains

!If we went until here then domain is correct
if (kkdebug) then
  write(379,fmt='("== > Domain ", I4," is consistent for seg ",I9,":")') seg(i)%dom, i
  write(379,fmt='("")')
endif

end subroutine check_domain

end module CONNEC

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
