module film2para_module

implicit none

integer, parameter  :: DPI=selected_int_kind(9)           ! Definition of the : double precision integer used in the simulation
! integer, parameter  :: DPI=selected_int_kind(13)          ! Definition of the : double precision integer used in the simulation
integer, parameter  :: DPI_S=selected_int_kind(9)         ! The integer format needed for the film.bin file
integer, parameter  :: DP=selected_real_kind(p=14)        ! Definition of the : double precision real used in the simulation
integer, parameter  :: TabDim=1000000                     ! Max number of segments accounted for in paraview
integer, parameter  :: nsegmax=100000
integer(kind=DPI_S) :: iiseg, ijunctt, tabvois, surface
integer(kind=DPI)   :: iseg, junctt,slipsys
real(kind=DP)       :: Vcourb

real(kind=DP)       :: curv_vect,seg_center

logical             :: foil_cut

dimension iiseg(5,nsegmax), ijunctt(nsegmax), tabvois(2,nsegmax), surface(nsegmax)
dimension iseg(5,TabDim), junctt(TabDim)
dimension Vcourb(3,nsegmax)
dimension curv_vect(3,nsegmax),seg_center(3,nsegmax)


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine reading the film.bin file (cf camera.f90) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine read_bin_file_camera(UnitEnsightCase, UnitEnsightGeo)
subroutine read_bin_file_camera()

  implicit none

!  integer, intent(in) :: UnitEnsightCase, UnitEnsightGeo

  real(kind=DP)       :: xstress,epso,avalue,deltat
  real(kind=DP)       :: nveclame(3)
  real(kind=DP)       :: HalfModur(3),thickness,DistCenter
  real(kind=DP)       :: UperSurf, LowerSurf, Altitude, IsegO(3)
  integer             :: cristallo,nsegmaxf,l,icompt
  integer(kind=DPI)   :: i,ijk,j,startstep,endstep
  integer(kind=DPI)   :: modur_ref(3)
  integer(kind=DPI_S) :: nbase, bveclin(3,96),bvecnor(3,96),modur(3),kk,nsegm, freqpause
  integer(kind=DPI)   :: GIMIX,GIMIY,GIMIZ,GISEGM,GIMODUR(213)
  integer(kind=DPI)   :: affnboimagex,affnboimagey,affnboimagez
  integer(kind=DPI)   :: KODERR
  integer(kind=DPI)   :: veclame(3)

  character(len=50)   :: fichbase,filename,filename1
  character(len=1)    :: caractest,switch
  character(len=6)    :: file_index                             !< Index of the debug_loop file in parallel computation
  character(len=6)    :: jacquesadit                            !< string of 6 zeros

  integer             :: nb_time_steps, previous_kk,deltakk

  ! Variables used in the vector curvature calculation
  integer(kind=DPI)   :: P1(3),P2(3),P3(3),Oi(3),Ei(3),VNN,compteur
  integer(kind=DPI)   :: CompteurG, CompteurD   !< for radius of curvature statistics
  integer(kind=DPI)   :: ivli, vli, long,npoints
  integer(kind=DPI)   :: Rlocal, Rlocal2
  integer(kind=DPI)   :: BUFTLRO(1:3,3), I1O, I1E, NBDIM(3), DECAL(3)
  real(kind=DP)       :: Ldis_act, alpha
  real(kind=DP)       :: norme        !< for radius of curvature statistics


  logical,dimension(NSEGMAX)          :: tabGD    !< The GD marker
  logical,dimension(NSEGMAX)          :: tabjonc  !< The jonction marker

  integer (kind=DPI)                  :: XLOMAX           !< The discretization length
  integer(kind=DPI),  dimension(3)    :: ROO, trAB, Tr
  real(kind=DP),      dimension(3)    :: RA, RB

  real, parameter     :: fac_rlocal = 1.5   !< The factor applied to the rlocal cste (see vector curvature calculation)

  integer :: status_input !< Existence of the input file
  integer :: line !< Line counter
  integer :: NR=-1   !< Number of line in input file
  integer :: ios  !< IOSTAT return

  character(LEN=1) :: junk !< junk variable

  integer(kind=DPI)   :: div, incr, incr2
  div = 1
  incr = 0
  incr2 = 0
  jacquesadit ='00000'

  nveclame=(/0.,0.,0./)
  NR=0
  endstep=0
  foil_cut = .false.      ! Default solution is no thin foil cut

  switch='n' !Initialisation : we do not want to reuse values of f2p_def file

  status_input=access('../in/mM_Data_Files/f2p_dat','r')
  IF (status_input == 0) THEN
    OPEN(89,file='../in/mM_Data_Files/f2p_dat',STATUS='OLD')
    NR=0
    DO line=1,50
      READ(89,*,IOSTAT=ios) junk
      IF(ios/=0)EXIT
      IF(line==50)THEN
        write(*,*) "Error: Maximum number of records exceeded..."
        write(*,*) "Exiting program now..."
        STOP
      ENDIF
      NR=NR+1
    ENDDO
    CLOSE(89)
  ENDIF

  IF ((status_input == 0) .and. (NR>=9)) THEN
    OPEN(89,file='../in/mM_Data_Files/f2p_dat',STATUS='OLD')
    read (89,*) affnboimagex,affnboimagey,affnboimagez
    read (89,*) veclame(1),veclame(2),veclame(3)

    nveclame(1) = real(veclame(1))/sqrt(real(veclame(1))**2+real(veclame(2))**2+real(veclame(3))**2)
    nveclame(2) = real(veclame(2))/sqrt(real(veclame(1))**2+real(veclame(2))**2+real(veclame(3))**2)
    nveclame(3) = real(veclame(3))/sqrt(real(veclame(1))**2+real(veclame(2))**2+real(veclame(3))**2)

    read (89,*) thickness

    read (89,*) DistCenter

    read (89,*) freqpause

    read (89,*) startstep

    read (89,*) endstep

    read (89,*) Ldis_act

    read (89,*) slipsys

    close(89)

    write(*,*)
    write(*,'(A)') "--------------------------------------------------------------------"
    write(*,'(A)') " Parameters in file ../in/mM_Data_Files/f2p_dat (Last used settings)"
    write(*,'(A)') "--------------------------------------------------------------------"
    write(*,'(A62,3I3)')  "> Number of periodic replica used in X,  Y and Z direction:", affnboimagex,affnboimagey,affnboimagez
    write(*,'(A62,3I3)')  "> Thin foil direction vector:                              ", veclame(1),veclame(2),veclame(3)
    write(*,'(A62,F9.3)') "> Thin foil thickness (micron):                            ",thickness
    write(*,'(A62,F9.3)') "> Distance from the reference cell center (micron):        ", DistCenter
    write(*,'(A62,I9)')   "> Pause frequency:                                         ",freqpause
    write(*,'(A62,I9)')   "> Starting step:                                           ",startstep
    write(*,'(A62,I9)')   "> Ending step:                                             ",endstep
    write(*,'(A62,F9.3)') "> Discretization length:                                   ",Ldis_act
    write(*,'(A62,I9)')   "> Curvature radius are computed for slip system n#         ",slipsys

    write(*,*)
    write(*,*) " < Do you want to change these settings?? (y/n) >"

    read(*,*) switch

  ENDIF

  IF ((status_input > 0) .or. (switch=='y') .or. (NR<9)) THEN

    ! The number of replica
    write (*,*) "Number of periodic replica used in direction: X  Y  Z  ?"
    read (*,*,iostat=ios) affnboimagex,affnboimagey,affnboimagez

    DO WHILE ((affnboimagex<0) .or. (affnboimagez<0) .or. (affnboimagez<0) .or. (ios/=0))
      write (*,*) '!> Numbers of periodic replica must be equal or greater than zero <!'
      write (*,*) "Number of periodic replica used in direction: X  Y  Z  ?"
      read (*,*,iostat=ios) affnboimagex,affnboimagey,affnboimagez
    END DO

    ! The thin foil normal direction
    write(*,*) "Thin foil direction vector: X  Y  Z  ?"
    read (*,*,iostat=ios) veclame(1),veclame(2),veclame(3)

    nveclame(1) = real(veclame(1))/sqrt(real(veclame(1))**2+real(veclame(2))**2+real(veclame(3))**2)
    nveclame(2) = real(veclame(2))/sqrt(real(veclame(1))**2+real(veclame(2))**2+real(veclame(3))**2)
    nveclame(3) = real(veclame(3))/sqrt(real(veclame(1))**2+real(veclame(2))**2+real(veclame(3))**2)

    ! The thin foil thickness
    write(*,*) "Thin foil thickness (micron): X ?"
    read (*,*,iostat=ios) thickness

    DO WHILE ((thickness<=0) .or. (ios/=0))
      write (*,*) '!> Thin foil thickness must be greater than zero <!'
      write(*,*) "Thin foil thickness (micron): X ?"
      read (*,*,iostat=ios) thickness
    END DO

    ! The DistCenter
    write(*,*) "Distance from the reference cell center (micron): X ?"
    read (*,*,iostat=ios) DistCenter

    DO WHILE ((DistCenter<0) .or. (ios/=0))
      write (*,*) '!> Distance from the reference cell center must be equal or greater than zero <!'
      write(*,*) "Distance from the reference cell center (micron): X ?"
      read (*,*,iostat=ios) DistCenter
    ENDDO

    ! The freqpause
    write (*,*) "Pause frequency: ? "
    read (*,*,iostat=ios) freqpause

    DO WHILE ((freqpause<=0) .or. (ios/=0))
      write (*,*) '!> Pause frequency must be greater than zero <!'
      write (*,*) "Pause frequency:  ?"
      read (*,*,iostat=ios) freqpause
    END DO

    ! The startstep
    write (*,*) "Starting step:  ?"
    read (*,*,iostat=ios) startstep

    DO WHILE ((startstep<0) .or. (ios/=0))
      write (*,*) '!> Starting step must be equal or greater than zero <!'
      write (*,*) "Starting step:  ?"
      read (*,*,iostat=ios) startstep
    END DO

    ! The endstep
    write (*,*) "Ending step:  ?"
    read (*,*,iostat=ios) endstep

    DO WHILE ((endstep - startstep <=0) .or. (ios/=0))
      write (*,*) '!> Ending step must be greater than the starting step <!'
      write (*,*) "Ending step:  ?"
      read (*,*,iostat=ios) endstep
    END DO

    ! The Ldis_act
    write (*,*) "The discretization length used during the simulation (micron):  ? (if 0 : No curvature vector calculation)"
    !write (*,'(A,A)') "This parameter is needed for the curvature vector calculation",&
    !"and must be set to 0 if the latter calculation is not needed!"
    read (*,*,iostat=ios) Ldis_act

    DO WHILE ((Ldis_act<0) .or. (ios/=0))
      write (*,*) "The discretization length used during the simulation (micron):  ? (if 0 : No curvature vector calculation)"
      read (*,*,iostat=ios) Ldis_act
    END DO

    ! The slipsys
    write (*,*) "Slip system you want to compute curvature radius vector (>12 all slip systems)"
    read (*,*,iostat=ios) slipsys

    DO WHILE ((slipsys <=0) .or. (ios/=0))
      write (*,*) '!> Slip system number must be greater than zero <!'
      write (*,*) "Slip system you want to compute curvature radius vector (>12 all slip systems)"
      read (*,*,iostat=ios) slipsys
    END DO

    open(89,file='../in/mM_Data_Files/f2p_dat',STATUS='UNKNOWN')

    write(89,*) affnboimagex,affnboimagey,affnboimagez
    write(89,*) veclame(1),veclame(2),veclame(3)
    write(89,*) thickness
    write(89,*) DistCenter
    write(89,*) freqpause
    write(89,*) startstep
    write(89,*) endstep
    write(89,*) Ldis_act
    write(89,*) slipsys

    close(89)

  ENDIF

  !** Ouverture de l unite contenant le film
  !      print *, " Nom du fichier film dans le repertoire dd/out ?"
  !      read(*,*) fichier

  open(49,file='../out/film.bin',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)

  read(49,IOSTAT=KODERR) cristallo
  read(49,IOSTAT=KODERR) nsegmaxf
  read(49,IOSTAT=KODERR) avalue       ! unite de simulation: relation espace reel-valeurs entieres (real)
  read(49,IOSTAT=KODERR) nbase        ! nombre de vecteurs de discretisation (entier)
  read(49,IOSTAT=KODERR) (bveclin(:,I),bvecnor(:,I), I=1,nbase) ! les 3 composantes du vecteur ligne et plan
  read(49,IOSTAT=KODERR) modur(:)     ! dimensions de la boite (3 entier)

  if(KODERR < 0) then
     print*,"  End of the film.bin file"
     stop
  elseif(KODERR > 0) then
     print*,"Reading problem 1",KODERR
     stop
  endif

  if (nsegmax /= nsegmaxf) then
    write(*,*) 'The standard value of nsegmax was modified in this film : STOP !'
    read(*,*)
  endif

  ! en cas de probleme avec les redemarages par exemple
  if (cristallo > 4 .or. cristallo < 1) then
    write(*,*) 'le fichier de la base de vecteur n est pas connu'
    write(*,*) ' Cristall lu est : ',cristallo
    write(*,*) 'sans doute un probleme de redemarage'
    write(*,*) 'CS=1, CC=2, CFC=3, HC=4'
    write(*,*) 'quelle base de vecteur doit on prendre ?'
    read (*,*) cristallo
  endif

  if (cristallo == 1) fichbase = "CS"
  if (cristallo == 2) fichbase = "BCC"
  if (cristallo == 3) fichbase = "CFC"
  if (cristallo == 4) fichbase = "HCP"

  if (cristallo > 4 .or. cristallo < 1) then
    print*, "cristallo : valeur inconnue"          ! No way
    stop
  endif

  print *, "  Base de vecteur utilisee est de type ",fichbase

  ! Defintion of the regularization length used in the curvature definition for
  ! the local line tension. This length is scaled by the discretization length
  ! All the segments participating to the local line tension will be excludded
  ! from the interactions (self-interaction)
  Ldis_act  = Ldis_act * 1D-6               ! the discretization length in meter
  Xlomax    = NINT(Ldis_act/avalue,DPI)     ! the discretization length in the simulation units
  rlocal    = INT(xlomax * fac_rlocal,DPI)  ! the reference length used for the curvature (line tension) calculation
  rlocal2   = rlocal * rlocal

  ijk = -1

  ! Quantities needed for the thin foil process
  ! Simulation Units
  thickness = thickness * 1.e-6 / avalue
  DistCenter = DistCenter * 1.e-6 / avalue
  ! Reference cell center coordinates
  HalfModur(1:3) = real(Modur(1:3))/2.
  ! Limits of the thin foil
  UperSurf  = (dot_product(HalfModur,nveclame) + DistCenter) + (0.5 * thickness)
  LowerSurf = (dot_product(HalfModur,nveclame) + DistCenter) - (0.5 * thickness)

  !**************************************************
  !***  Start of the time loop
  nb_time_steps = 0
  previous_kk   = 0
  kk = 0
  l=0

  ! The header analysis is finished, we know close film.bin to restart again
  close(49)


  open(50,file='../out/film.bin',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)

  read(50,IOSTAT=KODERR) cristallo
  read(50,IOSTAT=KODERR) nsegmaxf
  read(50,IOSTAT=KODERR) avalue       ! unite de simulation: relation espace reel-valeurs entieres (real)
  read(50,IOSTAT=KODERR) nbase        ! nombre de vecteurs de discretisation (entier)
  read(50,IOSTAT=KODERR) (bveclin(:,I),bvecnor(:,I), I=1,nbase) ! les 3 composantes du vecteur ligne et plan
  read(50,IOSTAT=KODERR) modur(:)     ! dimensions de la boite (3 entier)

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  Steps_loop: do while (kk >= 0)

    ijk = ijk + 1

    !*** The time step header
    read(50,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO

    if (kk > endstep) exit Steps_loop  ! No need to go further we have finished

    ! This step must be analysed since
    if ( (kk >= startstep) .and. (modulo(kk, freqpause)==0) ) then

      ! Progress of the extraction is given by printing the header at each time step
      write(*,'(I7,", time:",ES10.2E2,",Eps:",ES10.2E2,"%, Sigma =",ES10.2E2," MPa, Nseg = ",I8)')     &
            KK, deltat * kk, EPSO, XSTRESS, nsegm

      if(KODERR < 0) then
          print*,"  Fin de fichier FILM"
          goto 1000
          ! stop
      elseif(KODERR > 0) then
          print*,"probleme lecture 1",KODERR,ijk
          cycle Steps_loop
      endif

      ! GIMODUR definitions
      ! This table is used to transfer data between the c and Fort codes

      ! Dimension of the simulation cell to visualize
      GIMODUR(1)= MODUR(1)*(affnboimagex+1)
      GIMODUR(2)= MODUR(2)*(affnboimagey+1)
      GIMODUR(3)= MODUR(3)*(affnboimagez+1)

      ! The shifting informations, for the moment set to zero
      ! POUR LE MOMENT LES DECALAGES SONT A ZERO, AU PASSAGE IL SERAIT
      ! PEUT ETRE BON D'UTILISER SIX PETITS SCALAIRES INITIALISES UNE BONNE
      ! FOIS POUR TOUTE AU LIEU DE FAIRE DES SUMS A CHAQUE OPERATION DE CLP...
      GIMODUR(4)=0 !shiftxy
      GIMODUR(5)=0 !shiftxz
      GIMODUR(6)=0 !shiftyx
      GIMODUR(7)=0 !shiftyz
      GIMODUR(8)=0 !shiftzx
      GIMODUR(9)=0 !shiftzy

      ! In the following a Base of 8 vectors per slip systems is assumed
      ! but the graphical interface is more generale than this and works
      ! with more or less sofisticated simulation base of discretization
      GIMODUR(10)=8                 ! nb of vecteurs per slip systems
      GIMODUR(11)=NBASE/8           ! nb of slip systems taken into account
      GIMODUR(12)=nbase             ! total number of vectors

      ! The space scalling factor
      GIMODUR(13)=int(0.000001/avalue,DPI)

      do I = 1,NBASE
        GIMODUR(13+i) = (i-1)/8 + 1     ! Index of slip systems (needed to the
      enddo                             ! graphical information function

      ! The segments information is loaded
      READ(50,IOSTAT=KODERR) (iiseg(1:5,J),ijunctt(j),tabvois(1:2,j),surface(j), J=1,nsegm)

      if(KODERR < 0) then
        print*,"  Fin de fichier FILM"
        goto 1000
        ! stop
      elseif(KODERR > 0) then
        print*,"probleme lecture 2"
        kk=1
        cycle Steps_loop
      endif

      ! As not special procedure exist today for particles, the later are subtracted from the paraview analysis
      do j=1,nsegm
        if (tabvois(1,j) == 0 .and. tabvois(2,j) == 0) then
          nsegm = j-1   ! The real number of segments
          exit
        endif
      enddo

      ! We save the exact shape (3) of the simulation box
      moduR_ref(3) = moduR(3)

      ! Test to build segments replica in periodic boxes
      if (affnboimagex.ne.0 .or. affnboimagey.ne.0 .or. affnboimagez.ne.0) then

        ! Initialization
        GISEGM=0
        gimix=0
        gimiy=0
        gimiz=0

        ! The replicas loops
        do gimix = -affnboimagex, affnboimagex
        do gimiy = -affnboimagey, affnboimagey
        do gimiz = -affnboimagez, affnboimagez

          icompt = 0

          do i = 1, nsegm

            IsegO(1) = iISEG(1,i) + gimix * modur(1)
            IsegO(2) = iISEG(2,i) + gimiy * modur(2)
            IsegO(3) = iISEG(3,i) + gimiz * modur(3)

            ! The thin foil test
            altitude = dot_product(IsegO,nveclame)
            if (Altitude > Upersurf .or. Altitude < LowerSurf) then
              foil_cut = .true.
              cycle
            endif


            ! The segment is in the thin foil - > we keep the information
            icompt = icompt + 1
            ISEG(1,gISEGM+icompt) = int(IsegO(1),DPI)
            ISEG(2,gISEGM+icompt) = int(IsegO(2),DPI)
            ISEG(3,gISEGM+icompt) = int(IsegO(3),DPI)
            ISEG(4,gISEGM+icompt) = iiSEG(4,i)
            ISEG(5,gISEGM+icompt) = iiSEG(5,i)
            JUNCTT(gISEGM+icompt) = iJUNCTT(i)

          enddo

          gISEGM = gISEGM + icompt    ! The new number of segments

        enddo
        enddo
        enddo

      else

        icompt = 0

        do i = 1,nsegm

          IsegO(1) = iISEG(1,i)
          IsegO(2) = iISEG(2,i)
          IsegO(3) = iISEG(3,i)

          ! The thin foil test
          altitude = dot_product(IsegO,nveclame)
          if (Altitude > Upersurf .or. Altitude < LowerSurf) then
              foil_cut = .true.
              cycle
          endif
          icompt = icompt + 1
          ISEG(1:5,icompt) = iISEG(1:5,i)
          JUNCTT(icompt)   = iJUNCTT(i)

        enddo

        gISEGM = icompt

      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculation of the vector radius of curvature for the segments
      ! This calculation is essentially a copy and paste of the calculation made in mM
      ! to calculate the local line tension.
      ! The calculation is made only in the reference box without replicas
      VCourb(1:3,1:nsegm)   = 0.0            ! The curvature vector is initialized to zero

      if (Ldis_act /= 0.0) then

        if (affnboimagex == 0 .and. affnboimagey == 0 .and. affnboimagez == 0) then

          ! The junction and GD information is recovered
          tabGD(1:nsegm)   = .false.
          tabjonc(1:nsegm) = .false.

          do i = 1, nsegm
            if (JUNCTT(i) < 0) tabGD(i)   = .true.   ! This segment is a GD segment
            if (JUNCTT(i) > 0) tabjonc(i) = .true.   ! this segment is part of a junction
          enddo

          do i = 1, nsegm

            long          = ISEG(4,I)
            VLI           = ISEG(5,I)
            i1o           = tabvois(1,I)
            i1e           = tabvois(2,I)

            if (long == 0) cycle     ! The curvature vector cannot be calculated on pivot segments

            ! DETERMINATION of the points needed for the calculation of the local radius of curvature:
            ! P2 is always the center of I
            P2(1:3) = iseg(1:3,i) + (long * bveclin(1:3,VLI)) / 2

            ! P1 is :
            ! - the origin of I, if I1o is a pinning point
            ! - the center of I1o, if I1o is a real segment (not a part of a long segment discretized in small pieces)
            ! - the center of the effective neighbor before discretization
            if (i1o  == nsegmax .or. tabjonc(i1o)) then

              P1(1:3) = iseg(1:3,i)

            else
              ! Oi : is the beginning of the effective I1o
              ! Ei : is the end of the effective I1o
              Ei(1:3) = iseg(1:3,i)

              ! to avoid PBC problem, all coordinates are calculated relatively to the origin of I
              Tr(1:3) = iseg(4,I1o) * bveclin(1:3,ISEG(5,i1o)) / 2
              Oi(1:3) = Ei(1:3) - Tr(1:3) * 2

              compteur = 1

              if (RLOCAL2 < 0.0) then
                ! P1 is the center of the vnn segment
                P1(1:3) = Oi(1:3) + Tr(1:3)
              else

                ! The vector connecting P2 and the center of i1o
                ROO(1:3) = (Oi(1:3) + Tr(1:3)) - P2(1:3)

                ! Do we need to look for more distant segment i10 to regularize the line tension curvature ?
                if (SUM(real(ROO,DP)*real(ROO,DP)) < Rlocal2) then

                  !The half length of i1o segment
                  TrAB(1:3) = Tr(1:3)

                  ! Ra and Rb are vectors between P2 and the extremities of i1o
                  RA(1:3)=real(ROO(1:3) + TrAB(1:3),DP)
                  RB(1:3)=real(ROO(1:3) - TrAB(1:3),DP)

                  vnn = i1o

                  ! exact calculation of the line tension impose that the local radius of curvature is calculated
                  ! with points P1 and P3 not very close from P2. A minimum distance of Rlocal is imposed!

                  do while (tabvois(1,vnn) /= nsegmax  .and.          &
                            .not. tabjonc(vnn)       .and.          &
                            .not. tabgd(vnn)         .and.          &
                            SUM(RA*RA) < Rlocal2     .and.          &
                            SUM(RB*RB) < Rlocal2)

                    ! The next i1o segment along the dislocation line is tested
                    vnn = tabvois(1,vnn)
                    compteur = compteur + 1

                    ! The vector connecting P2 and the center of previous vnn
                    Tr(1:3) = TrAB(1:3)

                    !The half length of the new vnn segment
                    IVLI = iseg(5,vnn)
                    TrAB(1:3) = (iseg(4,vnn)*bveclin(1:3,IVLI)) / 2

                    !The center of the new vnn segment
                    ROO(1:3) = ROO(1:3) - (Tr(1:3) + TrAB(1:3))

                    ! Ra and Rb are vectors between P2 and the extremities of new i1o
                    RA(:)=real(ROO(:)+TrAB(:),DP)
                    RB(:)=real(ROO(:)-TrAB(:),DP)

                    ! We exit the loop if more than 10 segments were tested
                    if(compteur > 10) exit

                  enddo

                  ! P1 is the center of the vnn segment
                  P1(1:3) = P2(1:3) + ROO(1:3)
                endif

              endif

              ! Now we look for the problem of long straight lines cut in many segments,
              ! straight lines section must be treated as one long segment, whatever rlocal
              if (compteur == 1) then

                ! The initial vnn is redefined
                vnn = tabvois(1,i1o)

                do while (vnn /= nsegmax .and. iseg(5,vnn) == iseg(5,i1o))

                  Oi(1:3) =  Oi(1:3) - iseg(4,vnn) * bveclin(1:3,iseg(5,vnn))
                  compteur = compteur + 1
                  vnn = tabvois(1,vnn)

                  ! no more than 10 aligned segments are considered
                  if(compteur > 10) exit

                enddo

                P1(1:3) = (Oi + Ei) / 2     ! The new P1

              endif

              compteurG = compteur  ! the compteur value is saved for statistics

            endif  ! End of the tests on the P1 side

            ! P3 is :
            ! - the end     of I if I1e is pinning point
            ! - the center  of I1e if I1e is a real segment (not a part of a long segment beeing descretised)
            ! - the center of the effctive neighbor before descretisation
            if (i1e  == nsegmax .or. tabjonc(i1e)) then

              P3(1:3) = iseg(1:3,i) + long * bveclin(1:3,VLI)  ! end of I

            else

              ! begining of effective I1e = end of I
              Oi(1:3) = iseg(1:3,i) + long * bveclin(1:3,VLI)

              ! starting end of I1e
              Tr(:) =   iseg(4,i1e) * bveclin(1:3,iseg(5,i1e)) / 2
              Ei(1:3) = Oi(1:3) + Tr(1:3) * 2

              compteur = 1

              if (Rlocal2 < 0.0) then

                ! P3 is the center of the vnn segment
                P3(1:3) = Oi(1:3) + Tr(:)

              else

                ! The vector connecting P2 and the center of i1e
                ROO(1:3) = (Oi(1:3) + Tr(1:3)) - P2(1:3)

                ! Do we need to look for more distant segment i1e to regularize the line tension curvature ?
                if (SUM(real(ROO,DP)*real(ROO,DP)) < Rlocal2) then

                  !The half length of i1e segment
                  TrAB(:) = Tr(1:3)

                  ! Ra and Rb are vectors between P2 and the extremities of i1o
                  RA(:) = real(ROO(:) + TrAB(:),DP)
                  RB(:) = real(ROO(:) - TrAB(:),DP)

                  vnn = i1e

                  ! exact calculation of the line tension impose that the local radius of curvature is calculated
                  ! with points P1 and P3 not very close from P2. A minimum distance of Rlocal is imposed!
                  do while (tabvois(2,vnn) /= nsegmax       .and.     &
                            .not.tabjonc(vnn)               .and.     &
                            .not.tabgd(vnn)                 .and.     &
                            SUM(RA*RA) < Rlocal2            .and.     &
                            SUM(RB*RB) < Rlocal2)

                    ! The next i1e segment along the dislocation line is tested
                    vnn = tabvois(2,vnn)
                    compteur = compteur + 1

                    ! The vector connecting P2 and the center of previous vnn
                    Tr(1:3) = TrAB(1:3)
                    !The half length of the new vnn segment
                    IVLI = iseg(5,vnn)
                    TrAB(1:3) = (iseg(4,vnn)*bveclin(1:3,IVLI)) / 2

                    ! The vector connecting P2 and the center of vnn
                    ROO(1:3) = ROO(1:3) + (Tr(1:3) + TrAB(1:3))

                    ! Ra et Rb distances entre le centre de i et les extremites de vnn
                    RA(:)=real(ROO(:) - TrAB(:),DP)
                    RB(:)=real(ROO(:) + TrAB(:),DP)

                    if(compteur > 10) exit

                  enddo

                  ! P3 is the center of the vnn segment
                  P3(1:3) = P2(1:3) + ROO(1:3)

                endif

              endif

              ! Now we look for the problem of long straight lines cut in many segments
              ! straight lines section must be treated as one long segment, whatever rlocal
              if (compteur == 1) then

                ! The initial vnn is redefined
                vnn = tabvois(2,i1e)

                do while(vnn /= nsegmax .and. iseg(5,vnn) == iseg(5,i1e))

                  EI(1:3) = Ei(1:3) + iseg(4,vnn) * bveclin(1:3,iseg(5,vnn)) ! end of effective I1e
                  compteur = compteur + 1

                  vnn = tabvois(2,vnn)
                  if(compteur > 10) exit

                enddo

                P3(1:3) = (Oi + Ei) / 2

              endif

              CompteurD = Compteur ! We save compteur value for statistics

            endif     ! End of the P3 tests

            if (surface(i) == 0) then    ! The curvature vector has no meaning for the segments touching a surface or interface

              ! test on the existence of a local curvature (VCourb /= infinity)
              Decal(1:3) = P1(1:3) - P2(1:3)
              NBdim(1:3) = P3(1:3) - P2(1:3)

              ! the following test is true only when P1,P2,P3 are alligned
              if ((Decal(2)*NBdim(3) /= Decal(3)*NBdim(2)) .or.          &
                  (Decal(3)*NBdim(1) /= Decal(1)*NBdim(3)) .or.          &
                  (Decal(1)*NBdim(2) /= Decal(2)*NBdim(1)))    then

                BUFTLRO(1:3,1) = P1(1:3)
                BUFTLRO(1:3,2) = P2(1:3)
                BUFTLRO(1:3,3) = P3(1:3)

                ! Calcul du vecteur rayon de courbure
                VCOURB(1:3,i) = VECCOURBCRAMER(real(Decal),real(NBdim),alpha)

                norme = sqrt(VCOURB(1,i)*VCOURB(1,i) + &
                             VCOURB(2,i)*VCOURB(2,i) + &
                             VCOURB(3,i)*VCOURB(3,i) ) ! avec la fonction intrinsèque dot_product

                !the following lines are commented, they are useful to visualize curvature
                !vector in paraview: in case of large curvature vectors you can impose a threshold (default=10 times the discretization length)
                !to reduce the maximum size of the vectors in the graphical interface, to improve the visualization.

                !if (norme > Xlomax*20) then
                !   VCOURB(1:3,i)=VCOURB(1:3,i)*Xlomax*20/norme
                !endif


              endif

            endif

!             213  format(2x,3(I4,2x),2x, E13.5,2x, L3, 9(I4,2x), 2(E13.5,2x))
!             write(28,213) i,              &           ! numéro du noeud
!                           compteurG,      &           ! nb voisins gauche
!                           compteurD,      &           ! nb voisins droite
!                           norme                       ! norme of the curvature_radius vector

          enddo     ! End of the segment loop

        endif     ! End of the replica test for the curvature vector calculation

      endif     ! End of the test on the curvature vector calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! The integrity of the film file is tested
      read(50,IOSTAT=KODERR) caractest

      if(KODERR < 0) then
          print*,"  Fin de fichier FILM"
          goto 1000
          ! stop
      elseif(KODERR > 0) then
          print*,"probleme lecture 3",caractest
          kk=1
          cycle Steps_loop
      endif

      print*,'exporting data for kk=',kk
      deltakk = kk - previous_kk

      !call write_geo_ascii_file(UnitEnsightGeo, gisegm, bveclin,bvecnor,cristallo)
      nb_time_steps = nb_time_steps + 1

      open(10,STATUS='SCRATCH')
        write  (10,'(I0)') nb_time_steps
        rewind (10)
        read   (10,*) file_index
      close(10)


      incr = int(nb_time_steps/div)
      if (incr == 10) then
        div = div*10
        incr2 = incr2+1
      endif

      filename= "../out/vtk/film_"//trim(jacquesadit(1:5-incr2))//trim(file_index)//".vtp"

      ! vtk construction is possible form the moment only in the absence of thin foil cut
      if (.not. foil_cut) then

        open (unit=3, file=filename, status='unknown')

        call reorder_microstructure(gisegm, bveclin, npoints)

        close (3)

        if (Ldis_act /= 0.0) then


          filename1= "../out/vtk/curvature_"//trim(jacquesadit(1:5-incr2))//trim(file_index)//".vtp"

          open (unit=4, file=filename1, status='unknown')

          call curv_radius(npoints)

          close(4)

        endif

      endif

      previous_kk = kk

    else

      read(50,IOSTAT=KODERR)

      read(50,IOSTAT=KODERR)

      if(KODERR < 0) then
         print*,"  End of the film.bin file"
         stop
      elseif(KODERR > 0) then
         print*,"Reading problem 666",KODERR
         stop
      endif

      previous_kk = kk

    endif

  enddo Steps_loop

  close (50)

1000  continue

  write (*,*)
  write (*,*) "************************************************************************************"
  write (*,*) "***** Details for dislocation view settings in paraview in out/f2pcolorscales ******"
  write (*,*) "************************************************************************************"
  write (*,*)

  ! Apres lecture des donnes relatives a la construction du film
  ! on appelle write_case_file pour ecrire le fichier.geo correspondant a la
  ! geometrie
  call write_pvd_file(nb_time_steps,deltakk,deltat,startstep)
  !call write_case_file(UnitEnsightCase, "film.geo", nb_time_steps,deltakk,deltat)

end subroutine read_bin_file_camera


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine permettant d ecrire le film.case !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_pvd_file(nb_time_steps,deltakk,deltat,startstep)

  implicit none

  integer, intent(in)       ::nb_time_steps,deltakk,startstep
  real(kind=DP),intent(in)  :: deltat

  integer                   :: itime
  real(kind=DP)             :: dtime
  character(len=40)         :: cdtime,citime,filename
  character(len=6)          :: jacquesadit

  integer(kind=DPI)   :: div, incr, incr2
  div = 1
  incr = 0
  incr2 = 0
  jacquesadit ='00000'

  open (unit=4, file='../out/film.pvd', status='unknown')
  write (4,'(a)') '<?xml version="1.0"?>'
  write (4,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
  write (4,'(a)') '<Collection>'

  do itime = 1 , nb_time_steps
    dtime = ((itime-1) * deltat * deltakk) + (startstep * deltat)
    ! write (cdtime, '(E23.15)') dtime
    open(10,STATUS='SCRATCH')
      write  (10,'(E23.15)') dtime
      rewind (10)
      read   (10,*) cdtime
    close(10)

    open(10,STATUS='SCRATCH')
      write  (10,'(I0)') itime
      rewind (10)
      read   (10,*) citime
    close(10)

    incr = int(itime/div)
    if (incr == 10) then
      div = div*10
      incr2 = incr2+1
    endif

    filename= "vtk/film_"//trim(jacquesadit(1:5-incr2))//trim(citime)//".vtp"

    write (4,'(a)',ADVANCE='NO') '<DataSet timestep="'//trim(adjustl(cdtime))//''
    write (4,'(a)') '" part="0" file="'//trim(adjustl(filename))//'"/>'
  end do

  write (4,'(a)') '</Collection>'
  write (4,'(a)') '</VTKFile>'
  close(4)

end subroutine write_pvd_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reorder_microstructure(gisegm, bveclin ,nptsctr)

  implicit none

  integer(kind=DPI), intent(in) :: gisegm   ! , kk
  integer(kind=DPI_S), intent(in) :: bveclin(3,96)!,bvecnor(3,96)
  integer(DPI)  :: i,Ei(3),dislo,j,ntotpoints,ntotlines,next
  integer(DPI)  :: linepoints,totelem,limloop,limlinepoints!points(3,nsegmax),
  integer(DPI)  :: slipsystem(gisegm),junctype(gisegm),juncbinom(gisegm),slipsystemjunc(gisegm),linesign(gisegm)
  integer(DPI)  :: is,js,nptsctr
  integer(DPI),dimension(12,12) ::tabjonc
  integer(DPI)  :: juncttold
  integer(DPI)  :: testsign,disloold,dislonext
  real(DP) :: points(3,nsegmax),seg_center(3,nsegmax),curv_vect(3,nsegmax)
  logical :: zipseg(nsegmax), firstloop
  integer(DPI), allocatable :: loop(:,:),segfirst(:)
  ! real(DP) :: start, finish

  !call cpu_time(start)
  ! reduction of the maximum number of loops and segments to reduce the size memory allocation
  limloop  = int(gisegm*0.35) + 1
  !limloop  = int(gisegm)
  limlinepoints = int(nsegmax*0.2)

  allocate(loop(limloop,limlinepoints))
  allocate(segfirst(limloop))

  segfirst(1:limloop) = 0

  zipseg(1:nsegmax) = .false.
  slipsystem(1:gisegm) = 0
  junctype(1:gisegm) = 0
  juncbinom(1:gisegm) = 0
  slipsystemjunc(1:gisegm) = 0
  linesign(1:gisegm) = 0

  !loop(1:limloop,1:limlinepoints) = 0

  tabjonc(:,:) = 0

  tabjonc(1,:)=(/0,4,2,3,2,0,3,2,0,2,1,1/)
  tabjonc(2,:)=(/4,0,0,2,3,2,2,0,2,3,1,1/)
  tabjonc(3,:)=(/2,0,0,4,2,3,2,0,1,1,2,3/)
  tabjonc(4,:)=(/3,2,4,0,0,2,3,2,1,1,0,2/)
  tabjonc(5,:)=(/2,3,2,0,0,4,1,1,2,3,0,2/)
  tabjonc(6,:)=(/0,2,3,2,4,0,1,1,0,2,2,3/)
  tabjonc(7,:)=(/3,2,2,3,1,1,0,4,2,0,2,0/)
  tabjonc(8,:)=(/2,0,0,2,1,1,4,0,3,2,3,2/)
  tabjonc(9,:)=(/0,2,1,1,2,0,2,3,0,4,3,2/)
 tabjonc(10,:)=(/2,3,1,1,3,2,0,2,4,0,2,0/)
 tabjonc(11,:)=(/1,1,2,0,0,2,2,3,3,2,0,4/)
 tabjonc(12,:)=(/1,1,3,2,2,3,0,2,2,0,4,0/)

  Do i=1,12
    Do j=1,12
      if (tabjonc(i,j) /= tabjonc(j,i)) then
        print *, i, j
        stop "Problem here "
      endif
    enddo
  enddo

  is =0
  js =0

  i = 1

  !to write the .vtp file the connectivity is rebuilt
  if (gisegm < nsegmax) then

    ! 1. We discard segments with zero lenght
    do dislo = 1, gisegm
      if(iseg(4,dislo) == 0 .and. junctt(dislo)==0)  zipseg(dislo) = .true. !kneecaps are not considered, they have zero length so we don't display it
      if(iseg(4,dislo) == 0 .and. junctt(dislo)/=0) then  !GD segments are treated as kneecaps, but in ths case we need to rebuild the connectivity
        zipseg(dislo)=.true.
        tabvois(1,tabvois(2,dislo)) = tabvois(1,dislo)
        tabvois(2,tabvois(1,dislo)) = tabvois(2,dislo)
      endif
      if (tabvois(1,dislo) > gisegm ) then
        i = i+1
        segfirst(i)= dislo
      endif
    enddo

    segfirst(1) = i
    firstloop   = .false.
    if (segfirst(1) > 1) firstloop = .true.

    ntotpoints = 0
    ntotlines = 0
    nptsctr = 0

    disloold = 0
    dislonext = 0
    dislo = 1

    j = 1

    ! Begin main loop
    Do while(disloold < gisegm )!.not. ALL(zipseg(1:gisegm))) !.and. (dislo < gisegm)

      if (firstloop) then !Pinning point first

        dislo = segfirst(j+1)
        j = j+1
        if (j==segfirst(1)) firstloop=.false.

      else ! Other segments - The algorithm follows PBC segments before going on other segments

        if (dislonext /= 0) then
          dislo = dislonext
          dislonext = 0
        else
          dislo = disloold + 1
          disloold = dislo
        endif

      endif

      if (zipseg(dislo) .or. dislo > gisegm) cycle ! Segment has already taken into account

      linepoints = 0 ! initialize line points
      i = dislo

      ntotlines=ntotlines+1 ! A new line is created
      juncttold = junctt(i) ! Marker for junctions

      !Origin of the line
      ntotpoints=ntotpoints+1 ! Total points is increased
      linepoints=linepoints+1 ! Line point is increased
      points(1,ntotpoints)=iseg(1,i)
      points(2,ntotpoints)=iseg(2,i)
      points(3,ntotpoints)=iseg(3,i)

      loop(ntotlines,linepoints+1)=ntotpoints-1

      !--------------------------------------
      !Definition of the scalars for the line
      !--------------------------------------
      is = int((iseg(5,i)-1)/8)+1

      slipsystem(ntotlines) = is! int((iseg(5,i)-1)/8)+1  !iseg
      slipsystemjunc(ntotlines) =  is

      ! For junctype
      ! Hirth = 13
      ! Glissiles = 14 (No distinctions between both glissiles here)
      ! Lomer = 15

      if (junctt(i) /= 0) then
        js = int((iseg(5,JUNCTT(i))-1)/8)+1
        junctype(ntotlines)       =  tabjonc(is,js)               !junction type
        slipsystemjunc(ntotlines) =  13                           !Change to junction (13)
        juncbinom(ntotlines)  =  min(is,js)*100 + max(is,js)  !junction binom. e.g call 911 for junc between sys 9 and 11
        linesign(ntotlines) =  13

      else
        linesign(ntotlines) = int((iseg(5,i)-1)/8)+1
        testsign= iseg(5,i)-int((iseg(5,i)-1)/8)*8
        if (testsign <= 3 .or. testsign == 8) then
            linesign(ntotlines) =(int((iseg(5,i)-1)/8)+1)*2-1
        else
            linesign(ntotlines) = (int((iseg(5,i)-1)/8)+1)*2
        endif

      endif

      if (slipsys == int((iseg(5,i)-1)/8)+1 .or. slipsys >12) then
        nptsctr=nptsctr+1
        seg_center(1:3,nptsctr)=real((2.d0*iseg(1:3,i)+iseg(4,i)*bveclin(1:3,iseg(5,i)))*0.5d0,DPI)
        curv_vect(1:3,nptsctr) =VCOURB(1:3,i)
      endif

      !--------------------------------------
      !Look for the line end
      !--------------------------------------
      next = 0
      Do while (i <= gisegm .and. .not. zipseg(i))!.and. juncttold==junctt(i)) ! loop over segments

        zipseg(i) = .true. !Segment is taken into account

        next = tabvois(2,i)
        Do while (zipseg(next) .and. next <= gisegm .and. iseg(4,next) == 0 ) !Looking for next segment
          next = tabvois(2,next)
        enddo

        if ((int((iseg(5,i)-1)/8)+1 /= int((iseg(5,next)-1)/8)+1) &                ! System change
          & .or. (junctt(i)==0 .and. junctt(next) > 0) &                            ! Move from segment to junction
          & .or. ((junctt(i) > 0 .and. junctt(next) == 0) )) then   ! Move from junction to segment

          Ei(1:3) = iseg(1:3,i) + iseg(4,i)*bveclin(1:3,iseg(5,i))! End point is calculated

          ntotpoints=ntotpoints+1
          linepoints=linepoints+1
          points(1:3,ntotpoints)=Ei(1:3)


          loop(ntotlines,linepoints+1)=ntotpoints-1

          dislonext = tabvois(2,i)
          i = gisegm+1 ! loop is ended

        else

          if (next > gisegm) then !nsegmax point

            Ei(1:3) = iseg(1:3,i) + iseg(4,i)*bveclin(1:3,iseg(5,i))

            ntotpoints=ntotpoints+1
            linepoints=linepoints+1

            points(1,ntotpoints)=Ei(1)
            points(2,ntotpoints)=Ei(2)
            points(3,ntotpoints)=Ei(3)

            i = gisegm+1 ! loop is ended

            loop(ntotlines,linepoints+1)=ntotpoints-1

          else

            Ei(1:3) = iseg(1:3,i) + iseg(4,i)*bveclin(1:3,iseg(5,i))

            if (Ei(1) /= iseg(1,next) .or. Ei(2) /= iseg(2,next) .or. Ei(3) /= iseg(3,next)) then ! Periodic boundary conditions

              ntotpoints=ntotpoints+1
              linepoints=linepoints+1

              points(1,ntotpoints)=Ei(1)
              points(2,ntotpoints)=Ei(2)
              points(3,ntotpoints)=Ei(3)

              i = gisegm+1 ! loop is ended
              dislonext = next

              loop(ntotlines,linepoints+1)=ntotpoints-1

            else

              if (iseg(5,i) /= iseg(5,next) .or. zipseg(next) ) then !Reduction of point number for straight lines
                ntotpoints=ntotpoints+1
                linepoints=linepoints+1
                points(1,ntotpoints)=Ei(1)
                points(2,ntotpoints)=Ei(2)
                points(3,ntotpoints)=Ei(3)

                loop(ntotlines,linepoints+1)=ntotpoints-1
              endif

              i = next

            endif

          endif
        endif

      enddo

      loop(ntotlines,1)=linepoints !Total points of the line

      if (linepoints==1) then
        print *, "Problem here"
        stop
      endif

  enddo

    !print *, "     > Points | Lines produced : ",ntotpoints,ntotlines

    do i=1,ntotlines
      totelem=totelem + loop(i,1)
    enddo


    write (3,'(a)') '<?xml version="1.0"?>'
    write (3,'(a)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
    write (3,'(x,a)') '<PolyData>'
    write (3,'(2x,a,I8,a,I8,a)') '<Piece NumberOfPoints="',ntotpoints,'" NumberOfVerts="0" NumberOfLines="',ntotlines,'"'
    write (3,'(9x,a)') 'NumberOfStrips="0" NumberOfPolys="0">'

    write (3,'(3x,a)') '<Points>'
    write (3,'(4x,a)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
    do i = 1, ntotpoints
      write (3, '(5x,g14.6,x,g14.6,x,g14.6)' ) points(1:3,i)
    end do
    write (3,'(4x,a)') '</DataArray>'
    write (3,'(3x,a)') '</Points>'

    write (3,'(3x,a)') '<Lines>'
    write (3,'(4x,a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
    do i = 1, ntotlines
      do j= 2, loop(i,1) + 1
          if (j==loop(i,1) + 1) then
            write ( 3, '(I8)') loop(i,j)
          else
            write ( 3, '(I8)', ADVANCE='NO') loop(i,j)
          endif
        enddo
    enddo
    !write (3,'(5x,a)')'0 1 2 3 4 5 6 7 8'
    write (3,'(4x,a)') '</DataArray>'

    write (3,'(4x,a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
    !write (3, '(5x,a)')'3 6 9'
    totelem=0
    do i = 1, ntotlines
      totelem =totelem + loop(i,1)
      if (i==ntotlines) then
      write ( 3, '(I8)') totelem
      else
      write ( 3, '(I8)', ADVANCE='NO') totelem
      endif
    enddo
    write (3,'(4x,a)') '</DataArray>'
     write (3,'(3x,a)') '</Lines>'

    write (3,'(3x,a)') '<CellData Scalars="Systems and Junctions, Systems only ,&
& Junctions type,Junction binoms,Dislocation line sign">' !,Junction binom

    write (3,'(4x,a)') '<DataArray type="Int32" Name="Systems and Junctions" format="ascii">'
    do i = 1, ntotlines
      write (3, '(5x,I4)') slipsystemjunc(i)
    enddo
    write (3,'(4x,a)') '</DataArray>'

    write (3,'(4x,a)') '<DataArray type="Int32" Name="Systems only" format="ascii">'
    do i = 1, ntotlines
      write (3, '(5x,I3)') slipsystem(i)
    enddo
    write (3,'(4x,a)') '</DataArray>'

    write (3,'(4x,a)') '<DataArray type="Int32" Name="Junctions type" format="ascii">'
    do i = 1, ntotlines
      write (3, '(5x,I3)') junctype(i)
    enddo
    write (3,'(4x,a)') '</DataArray>'

    write (3,'(4x,a)') '<DataArray type="Int32" Name="Junction binoms" format="ascii">'
    do i = 1, ntotlines
      write (3, '(5x,I4)') juncbinom(i)
    enddo
    write (3,'(4x,a)') '</DataArray>'

    write (3,'(4x,a)') '<DataArray type="Int32" Name="Dislocation line sign" format="ascii">'
    do i = 1, ntotlines
      write (3, '(5x,I4)') linesign(i)
    enddo
    write (3,'(4x,a)') '</DataArray>'


    write (3,'(3x,a)') '</CellData>'


    write (3,'(2x,a)') '</Piece>'
    write (3,'(x,a)') '</PolyData>'
    write (3,'(a)') '</VTKFile>'
  else
    print *, " Segments information are not listed since Nsegm > Nsegmax"
  endif

!call cpu_time(finish)
!print '("Time = ",f6.3," seconds.")',finish-start

end subroutine reorder_microstructure



subroutine curv_radius(ntotpoints)
implicit none

  integer(DPI)  :: i,ntotpoints

!   do i=1,ntotpoints
!   write(32,*) seg_center(:,i),curv_vect(:,i)
!   enddo

  write (4,'(a)') '<?xml version="1.0"?>'
  write (4,'(a)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
  write (4,'(x,a)') '<PolyData>'
  write (4,'(2x,a,I8,a,I8,a)') '<Piece NumberOfPoints="',ntotpoints,'" NumberOfVerts="0" NumberOfLines="0"'
  write (4,'(9x,a)') 'NumberOfStrips="0" NumberOfPolys="0">'

  write (4,'(3x,a)') '<Points>'
  write (4,'(4x,a)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  do i = 1, ntotpoints
    write (4, '(5x,g14.6,x,g14.6,x,g14.6)' ) seg_center(1:3,i)
  end do
  write (4,'(4x,a)') '</DataArray>'
  write (4,'(3x,a)') '</Points>'

  write (4,'(3x,a)') '<PointData Vectors="curvature_radius">'
  write (4,'(4x,a)') '<DataArray type="Float32" Name="curvature_radius" NumberOfComponents="3" format="ascii">'
  do i = 1, ntotpoints
  write (4, '(5x,g14.6,x,g14.6,x,g14.6)') curv_vect(1:3,i)
  enddo
  write (4,'(4x,a)') '</DataArray>'
  write (4,'(3x,a)') '</PointData>'


  write (4,'(2x,a)') '</Piece>'
  write (4,'(x,a)') '</PolyData>'
  write (4,'(a)') '</VTKFile>'

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function cose2v(a,b)

implicit none

real(DP) :: cose2v,a(3),b(3),na,nb

na=norvect(a)
nb=norvect(b)
if (na.ne.(0.).and.nb.ne.(0.)) then
   cose2v = dot_product(a,b)/(na*nb)
   if (abs(cose2v) > 1.0) cose2v = 1.0 * sign(1.0D0,cose2v) !avoid numerical errors
else
   cose2v = 0.
endif

end function cose2v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function norvect(r)

implicit none

real(kind=DP) :: norvect,r(3)

norvect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

end function norvect

!###########################################################################
!> \brief
!!  Calcul du rayon de courbure pour un cercle passant par trois points
!!  M,P et Q
!!  Pour un rayon de courbure infinie le vecteur retourne par la
!!  fonction est nul.
!!  Alpha : angle entre M et Q sur le cercle
!!  T : vecteur tangant ie normal a R dans MPQ
!############################################################### 30/07/99 ##
!> \ingroup Vect_calc
!> \todo    Has to be translated in english / description of parameters is needed
!! \date    30/07/99
!! \param   OM        Distance between the circle center and one of the points of left vnn for the calculation of curvature radius
!! \param   OP        Distance between the circle center and the center of the considered segment for the calculation of curvature radius
!! \param   OQ        Distance between the circle center and one of the points of right vnn for the calculation of curvature radius
!! \param   alpha     Angle between OM and OQ

FUNCTION VECCOURBCRAMER(PMint,PQint,alpha) RESULT(VECCOURB)

implicit none

real(kind=DP), dimension(3)            :: VECCOURB
real(kind=DPI),dimension(3),intent(in) :: PMint,PQint
real(kind=DP)                          :: ALPHA1, ALPHA2,normVC2,invnormVC2
real(kind=DP)                          :: DET,DETx,DETy,DETz,arg1,arg2
real(kind=DP),intent(out)              :: ALPHA
real(kind=DP), dimension(4)            :: PM,PQ,PMvecPQ
real(kind=DP), dimension(3)            :: A,B,OM,OQ,OR
real(kind=DP), parameter               :: half = 0.5D0
real(kind=DP), parameter               :: PII  = 4. * DATAN(1.D0)

PM(1:3) = PMint
A(1:3) = half * PM(1:3)
PM(4) = PM(1)*A(1) + PM(2)*A(2) + PM(3)*A(3)


PQ(1:3) = PQint
B(1:3) = half * PQ(1:3)
PQ(4) = PQ(1)*B(1) + PQ(2)*B(2) + PQ(3)*B(3)


PMvecPQ(1) = (PM(2)*PQ(3) - PM(3)*PQ(2))
PMvecPQ(2) = (PM(3)*PQ(1) - PM(1)*PQ(3))
PMvecPQ(3) = (PM(1)*PQ(2) - PM(2)*PQ(1))
PMvecPQ(4)= 0.d0

!Cramer's rule to determinate the center of the circumference defined by the three
!points P,M,Q (P1,P2,P3 in the subroutine FORCE)

DET = PM(1)*PQ(2)*PMvecPQ(3) - PM(1)*PQ(3)*PMvecPQ(2) &
    - PM(2)*PQ(1)*PMvecPQ(3) + PM(2)*PQ(3)*PMvecPQ(1) &
    + PM(3)*PQ(1)*PMvecPQ(2) - PM(3)*PQ(2)*PMvecPQ(1)



if (DET < 1.d-15) then
  !The points P,M,Q are aligned, no curvature radius
  if (DET<0.d0) print *, DET
  VECCOURB(1) = 0.d0
  VECCOURB(2) = 0.d0
  VECCOURB(3) = 0.d0
  alpha=0.d0
else

  DETx = PM(4)*PQ(2)*PMvecPQ(3) - PM(4)*PQ(3)*PMvecPQ(2) &
       - PM(2)*PQ(4)*PMvecPQ(3) + PM(2)*PQ(3)*PMvecPQ(4) &
       + PM(3)*PQ(4)*PMvecPQ(2) - PM(3)*PQ(2)*PMvecPQ(4)

  DETy = PM(1)*PQ(4)*PMvecPQ(3) - PM(1)*PQ(3)*PMvecPQ(4) &
       - PM(4)*PQ(1)*PMvecPQ(3) + PM(4)*PQ(3)*PMvecPQ(1) &
       + PM(3)*PQ(1)*PMvecPQ(4) - PM(3)*PQ(4)*PMvecPQ(1)

  DETz = PM(1)*PQ(2)*PMvecPQ(4) - PM(1)*PQ(4)*PMvecPQ(2) &
       - PM(2)*PQ(1)*PMvecPQ(4) + PM(2)*PQ(4)*PMvecPQ(1) &
       + PM(4)*PQ(1)*PMvecPQ(2) - PM(4)*PQ(2)*PMvecPQ(1)

  OR(1) = DETx/DET     !x coordinate of the center of the circumference
  OR(2) = DETy/DET     !y coordinate of the center of the circumference
  OR(3) = DETz/DET     !z coordinate of the center of the circumference
  VECCOURB(1:3)=OR(1:3)  !curvature radius corresponds to the center of the circumference,
                         ! indeed the point P is located at (0,0,0) and the curvature radius is OR-P

  normVC2 = VECCOURB(1)*VECCOURB(1)+VECCOURB(2)*VECCOURB(2)+VECCOURB(3)*VECCOURB(3)
  invnormVC2 = 1/normVC2

  OM(1:3)=OR(1:3)-PM(1:3)
  OQ(1:3)=OR(1:3)-PQ(1:3)
!   if (kkdebug) then
!   write(379,*)OM(1)*OM(1)+OM(2)*OM(2)+OM(3)*OM(3)
!   write(379,*)OQ(1)*OQ(1)+OQ(2)*OQ(2)+OQ(3)*OQ(3)
!   write(379,*)OM
!   write(379,*)OQ
!   write(379,*) VECCOURB
!   write(379,*) normVC2
!   write(379,*) (OM(1)*VECCOURB(1)+OM(2)*VECCOURB(2)+OM(3)*VECCOURB(3))
!   write(379,*) (OQ(1)*VECCOURB(1)+OQ(2)*VECCOURB(2)+OQ(3)*VECCOURB(3))
!   write(379,*) (OM(1)*VECCOURB(1)+OM(2)*VECCOURB(2)+OM(3)*VECCOURB(3))/normVC2
!   write(379,*) (OQ(1)*VECCOURB(1)+OQ(2)*VECCOURB(2)+OQ(3)*VECCOURB(3))/normVC2
!   call FLUSH()
!   endif
  arg1=(OM(1)*VECCOURB(1)+OM(2)*VECCOURB(2)+OM(3)*VECCOURB(3))
  arg2=(OQ(1)*VECCOURB(1)+OQ(2)*VECCOURB(2)+OQ(3)*VECCOURB(3))
  !stop
  if (abs(arg1+normVC2) < 1.d-6) then
     ALPHA1=Pii*dsqrt(normVC2)
  else
     ALPHA1=Acos(arg1*invnormVC2)
  endif

  if (abs(arg2+normVC2) < 1.d-6) then
     ALPHA2=Pii*dsqrt(normVC2)
  else
     ALPHA2=Acos(arg2*invnormVC2)
  endif

  ALPHA=ALPHA1+ALPHA2
endif

END FUNCTION VECCOURBCRAMER


end module film2para_module


!!!!!!!!!!!!!!!!
! main program !
!!!!!!!!!!!!!!!!
program film2para


use film2para_module

implicit none

character(len=256) :: cmd



cmd = "rm -f ../out/vtk/film*.vtp"
call system(cmd)
cmd = "rm -f ../out/vtk/curvature*.vtp"
call system(cmd)
cmd = "rm -f ../out/film.pvd"
call system(cmd)

Print *, " "
Print *, " "
Print *, " "
Print *, "                                    MM    M  EEEE   GGGG    AAA    SSSS"
Print *, "                                    M M M M  EE    G       A   A  SS   "
Print *, "         MM MM  I  CCC  RRRR  OOOO  M  M  M  EEEE  G  GGG  AAAAA   SSS "
Print *, "         M M M  I  C    RRR   O  0  M     M  EE    GG  GG  A   A     SS"
Print *, "         M   M  I  CCC  R  R  0000  M     M  EEEE   GGGG   A   A  SSSS "
Print *, " "
Print *, "                                             --> VTK image creation    "
Print *, " "
Print *, " ========================================================================================"


call read_bin_file_camera()

end program film2para
