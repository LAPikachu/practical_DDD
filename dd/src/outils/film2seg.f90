program main

use CONSTANTES
use VARGLOB
use BRICAMAT
use VARBASE
use DEBUG
use INIT
use CONNEC

implicit none

real(kind=DP)       :: xstress

! Initialization
CFC = .true.

call GenSeg

! Variables defined globally in VARBASE or VARGLOB in DPI=13:
! BVECLIN; SEG; MODUR;
! In Film.bin, integers are coded in DPI=8.
! Intermediate local variables are therefore required

!#######################
        contains
!#######################

subroutine GenSeg

  implicit none

  integer,parameter     :: DPI_S=selected_int_kind(9)   ! The integer format needed for the film.bin file
  character             :: fichier*50,fichbase*50,caractest
  integer(kind=DPI_S)   :: icristallo,nsegmaxf, modurloc(3)
  integer(kind=DPI_S)   :: ijk,debutaffich
  integer(kind=DPI_S)   :: KODERR, kk, nsegmloc !nsegm
  integer(kind=DPI)     :: i,j
  integer(kind=DPI_S)   :: minseg
  integer(kind=DPI_S)   :: tabseg(9,nsegmax)
  integer(kind=DPI)     :: vo,ve,vnne,vnno
  integer(kind=DPI)     :: ijonc,itemp
  character(len=40)     :: file_seg_read, file_seg_film
  logical               :: kkdebug
  logical               :: corr(nsegmax), Fire

  integer(kind=DPI_S), allocatable,save ::  bvecnorloc(:,:)
  integer(kind=DPI_S), allocatable,save ::  bveclinloc(:,:)

  Print *, " ==========================================================================="
  Print *, "                                CopyRight"
  Print *, " "
  Print *, " mM (for microMegas) is an open source program of DD (Dislocation Dynamics)"
  Print *, " simulation originally developed at the 'Laboratoire d'Etude"
  Print *, " des Microstructure), CNRS-ONERA, FRANCE."
  Print *, " Copyright (C) 1993, 1996, 2000, 2002, 2004   B. Devincre, L. Kubin, M. Condat,"
  Print *, " C. Lemarchand, R. Madec, S. Groh, G. Monnet, J. Durinck, P. Carrez, C. de Sansal,"
  Print *, " M-G. Tagorti, S. Queyreau, S. Naamane, A. Roos, A. Vattre."
  Print *, "  "
  Print *, " This is a free software, and you are welcome to redistribute it under the"
  Print *, " conditions of the GNU General Public License (see the 'Copyright' file in"
  Print *, " the program distribution)."
  Print *, " ==========================================================================="
  Print *, "              "


  ! We first select the configuration of the film.bin file we want to read
  write(*,*) "Step number (in film.bin) you want to extract a SEG configuration ?"
  read (*,*) debutaffich

!   write(*,*) 'PLEASE NOTE THE FOLLOWING LIMITATIONS OF THE FILM.BIN FILE:'
!   write(*,*) 'one phase only, no surface, 12 slip systems'

  !3) LOOP TO FIND THE SELECTED STEP
  kk=1
  plusseg = 0
  Fire=.false.

  ! The file we want to analyze
  fichier = "film.bin"
  open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR,action='read')

  ! Lecture of informations contained on the file, not used for the calculation
  read(50,IOSTAT=KODERR) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
  write (*,*) icristallo
  read(50,IOSTAT=KODERR) nsegmaxf     ! The nsegmax value used in this simulation
  write (*,*) nsegmaxf
  read(50,IOSTAT=KODERR) avalue       ! The simulation lattice parameter (the simulation scaling)
  write (*,*) avalue
  read(50,IOSTAT=KODERR) nbase        ! Nb of vectors direction used for the material crystal symmetry
  write (*,*) nbase

  NBsysdev        = 1
  NTSG            = 12
  NstatControl    = 1
  corr(:)   = .false.
  call Allocation
  solli_sys(1:12) = 1
  allocate ( bveclinloc(1:3,1:nbase))
  allocate ( bvecnorloc(1:3,1:nbase))

  read(50,IOSTAT=KODERR) (bveclinloc(:,I),bvecnorloc(:,I), I=1,nbase)    ! table of vector direction and displacement for this sim

  write (*,*) "The BVD file we use for this extraction!"
  do I=1,nbase
      write (*,*) I, bveclinloc(1:3,I), bvecnorloc(1:3,I)
  enddo

  bveclin(:,:)=bveclinloc(:,:)
  bvecnor(:,:)=bvecnorloc(:,:)

  write(*,*) " Dimension of the 3D simulated volume in a value !"
  write(*,*) " "
  read(50,IOSTAT=KODERR) modurloc(:)     ! Dimension of the 3D simulated volume in a value
  write (*,*) modurloc(1:3)
  modur(:) = modurloc(:)

  kkdebug=.true.

  ! Test on the nsegmax definition
  if (nsegmax /= nsegmaxf) then
      write(*,*) 'The standard value for nsegmax was modified in this film !', nsegmax, nsegmaxf
      write(*,*) 'You must compile again histo.f90 and toolsbox.f90 with the correct value : STOP !'
      stop
  endif

  ! Selection of the BVD list of vector
  if (icristallo == 1) fichbase = "CS"
  if (icristallo == 2) fichbase = "BCC"
  if (icristallo == 3) fichbase = "CFC"
  if (icristallo == 4) fichbase = "HCP"
  if (icristallo == 5) fichbase = "ORT"
  if (icristallo == 6) fichbase = "MGO"
  if (icristallo == 7) fichbase = "DC"

  if (icristallo > 7 .or. icristallo < 1) then
      print*, "cristallo : Unknown value"
      stop
  endif
  print *, "  The crystal symmetry is of type ",fichbase

  ! LOOP ON THE SIMULATION STEPS
  ijk = -1
  Steps_loop_17: do while (kk.gt.0)

    ijk = ijk + 1

    ! READING OF THE SEGMENT DATA FROM FILM.BIN
    read(50,IOSTAT=KODERR) kk,deltat,NsegMloc,XSTRESS,EPSO
    write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I5)') KK,EPSO,XSTRESS, nsegmloc
    nsegm=nsegmloc
    if(KODERR < 0) then
      print*,"  End of the film.bin file"
      exit Steps_loop_17
    elseif(KODERR > 0) then
      print*,"Reading problem 1",KODERR,ijk
      kk=1
      cycle Steps_loop_17
    endif

    ! The segments micro structure is loaded
    read(50,IOSTAT=KODERR) (tabseg(1:9,J), J=1,nsegm)

    if(KODERR < 0) then
      print*,"  End of the film.bin file"
      exit Steps_loop_17
    elseif(KODERR > 0) then
      print*,"Reading Error 2",KODERR,nsegm,j
      cycle Steps_loop_17
    endif

    ! The steps we want to use for statistic
    if(kk >= debutaffich .and. .not. Fire) then ! Begin of the "Big If"

      ! SEG DATA ARE OVERIDEN BY FILM.BIN DATA
      do j=1,nsegm

        seg(j)%O(1:3)   = tabseg(1:3,j)
        seg(j)%norme    = tabseg(4,j)
        seg(j)%veclin   = tabseg(5,j)

        seg(j)%voiso    = tabseg(7,j)
        seg(j)%vnno     = tabseg(7,j)
        seg(j)%voise    = tabseg(8,j)
        seg(j)%vnne     = tabseg(8,j)

        seg(j)%jonc     = .false.
        seg(j)%Ijonc    = nsegmax
        seg(j)%gd       = izero
        SEG(j)%surface  = tabseg(9,j)

        if (tabseg(6,j) < 0) then

          seg(j)%gd     = iun

        elseif (tabseg(6,j) > 0) then

          seg(j)%jonc   = .true.
          seg(j)%Ijonc  = tabseg(6,j)

        endif

      enddo

61        format(I7,3I10,I7,I7,4I7,L3,I7,I3)

      ! UNMODIFIED SEG DATA FROM FILM.BIN
      write (*,*) " "
      write (*,*) "Read the segment configuration a step  KK = ",KK
      file_seg_read    =   '../out/SEG_Read'
      open  (9,FILE=file_seg_read,STATUS='UNKNOWN')
      write (9,'(12F4.0)') solli_sys(1:12)
      write (9,*) NsegM
      write (9,*) modur(1:3)
      do i=1,nsegm
        write (9,61) i, seg(i)%O(1:3), seg(i)%norme, seg(i)%veclin, seg(i)%vnno,seg(i)%vnno, &
                        seg(i)%vnne, seg(i)%vnne, seg(I)%JONC, seg(i)%Ijonc, seg(i)%gd
      enddo
      close (9)

      ! DEFINITION AND TABULATION OF ARRAYS
      call Tabulation
!         write (*,*) 'We pass subroutine Tabulation'

      ! SINCE FILM.BIN DOES NOT CONTAIN NON ZERO LENGTH NEIGHBOOR
      ! INFORMATION, WE NEED TO REDO THE CONNECTION USING CONNECIJ
      do i = 1, nsegm

!         print*,"I=", i, seg(i)%O(1:3), seg(i)%norme, seg(i)%veclin, seg(i)%vnno,seg(i)%vnno, &
!                         seg(i)%vnne, seg(i)%vnne, seg(I)%JONC, seg(i)%Ijonc, seg(i)%gd

        vo = seg(i)%voiso
        ve = seg(i)%voise

        if (vo == i .and. ve==i) cycle !Segment that will be removed

        vnno = seg(i)%vnno
        vnne = seg(i)%vnne

        ! search for ordinary pivot segments for which we lost the voiso and voise info
        ! Those segments are eliminated
        if ( seg(i)%norme==0 .and. seg(i)%gd < iun) then

          if ((seg(vnno)%voise /= i .and. seg(vnno)%vnne /= i) .and.    &
              (seg(vnne)%voiso /= i .and. seg(vnne)%vnno /= i))           then

            seg(i)%voiso=i ! we dont know how to connect i with its vnn
            seg(i)%voise=i ! i will be removed

          else

            call connecij(int(i,8), int(vnne,8),1)

          endif

        else

          call connecij(int(i,8), int(vnne,8),2)

        endif

      enddo

      write (*,*) "The segment connectivity was build !"

      ! NEXT WE HANDLE SEG DATA FOR SEGMENTS THAT WILL BE DELETED (PROCEDURE COPIED FROM CONNEC)
      minseg=0
      do i = 1, nsegm + plusseg

        vo = seg(i)%voiso
        ve = seg(i)%voise

        if (vo == i .and. ve ==i) then

!           write (*,*) 'suppress seg', i
          itemp=i
          IDEP(itemp) = Izero
          seg(itemp)%norme = izero
          seg(itemp)%voise = itemp
          seg(itemp)%voiso = itemp
          seg(itemp)%vnno = itemp
          seg(itemp)%vnne = itemp
          out(itemp) = .true.
          seg(itemp)%resdep = zero
          seg(itemp)%wait = IZERO
          SEG(Itemp)%bloquer=seg(j)%bloquer
          seg(itemp)%grain=seg(j)%grain
          seg(itemp)%dom=seg(j)%dom
          seg(itemp)%NPhase=izero

          if(seg(Itemp)%jonc) then
!             write (*,*) 'Suppressed seg is junct', itemp
            ijonc = seg(Itemp)%Ijonc
            seg(ijonc)%jonc = .false. ;  seg(ijonc)%Ijonc = nsegmax
            if(seg(ijonc)%norme == 0) call voisinage(int(ijonc,8),348)
            !             call seginfo(ijonc, " verification de connection pour le binome ")
            seg(itemp)%jonc = .false. ;   seg(itemp)%Ijonc = NSEGMAX
            print*," connecij: 4:on suprime une jonction et son binom ",itemp, ijonc
            print*, "KK=",KK, "            REF = ",2
            !        stop
          endif

          seg(itemp)%GD = izero
          idep(itemp) = 0

        endif

      enddo

      write (*,*) "The segment configuration was cleanup ! nsegm =", nsegm, "+", plusseg

      ! WRITING OF THE CLEANED SEG FILE
      file_seg_film    =   '../in/SEG_film.bin'
      open( 18,FILE=file_seg_film,STATUS='UNKNOWN')
      write (18,'(12F4.0)') solli_sys(1:12)
      write (18,*) NsegM + plusseg
      write (18,*) modur(1:3)

      do i = 1, nsegm + plusseg
        vo = seg(i)%voiso
        vnno = seg(i)%vnno
        ve = seg(i)%voise
        vnne = seg(i)%vnne
        ijonc = seg(i)%Ijonc
        ! Sometime it is useful to restart the simulation from the seg_save configuration and to
        ! change the value of nsegmax. Hence Nsegmax value is replaced by zero here
        if (ve == nsegmax) ve = 0
        if (vo == nsegmax) vo = 0
        if (vnno == nsegmax) vnno = 0
        if (vnne == nsegmax) vnne = 0
        if (ijonc == nsegmax) ijonc = 0
        write (18,61) i, seg(i)%O(1:3), seg(i)%norme, seg(i)%veclin, vo, vnno, ve, vnne, seg(I)%JONC,&
            Ijonc, seg(i)%gd
!         write (*,61)  i, seg(i)%O(1:3), seg(i)%norme, seg(i)%veclin, vo, vnno, ve, vnne, seg(I)%JONC,&
!             Ijonc, seg(i)%gd
      enddo

      close(18)

      write (*,*) "The segment configuration was saved in file: ", file_seg_film

      exit Steps_loop_17

      Fire = .true.
      !END OF STEPS OF INTEREST
    endif

    ! The mark for the end of step is read
    read(50,IOSTAT=KODERR) caractest

    if(KODERR < 0) then
      print*,"  End of the film.bin file"
      stop
    elseif(KODERR > 0) then
      print*,"Reading error 3",caractest
      cycle Steps_loop_17
    endif

    !END OF LOOP ON ITERATIONS

  enddo Steps_loop_17

  close(50)

end subroutine GenSeg

end program main
