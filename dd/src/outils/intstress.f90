include '../../src/simu/00copyright.f90'

PROGRAM intstress

!$ use OMP_LIB

IMPLICIT NONE

INTEGER,PARAMETER  :: DP=selected_real_kind(p=14)   !< Definition of real number kind
integer,parameter  :: DPI=selected_int_kind(9)      !< Definition of integer number kind

integer(kind=DPI),parameter   :: nsegmax=21000 !Maximum number of segments
integer(kind=DPI),parameter   :: nsegmax2=200000 !Maximum number of segments

integer(kind=DPI),parameter  :: nbasered = 8  !< Dimension de la base reduite

real(kind=DP), parameter  :: PI     = 3.14159265358979D0

type segment
    real(kind=dp)                   :: resdep   !< fraction non entiere du deplacement
    integer(kind=dpi),dimension(3)  :: o        !< coordonnees de l'origine des segments

    integer(kind=dpi)    :: veclin, &   !< indice du vecteur ligne et deplacement
    norme,                        &   !< longueur du segment (positif)
    voiso,                        &   !< indice du voisin de l'origine du segment
    voise,                        &   !< indice du voisin de l'extremite du segment
    vnno,                         &   !< indice du voisin non nul de l'origine du segment
    vnne,                         &   !< indice du voisin non nul de l'extremite du segment
    ijonc,                        &   !< indice du segment binome lors des jonctions
    tjonc,                        &   !< nb d'iterations d'existence de la jonction
    wait,                         &   !< nombre de pas immobiles
    grain                             !< localisation du segment dans un polycristal

    logical ::  jonc,             &   !< clef d'immobilisation pour les jonctions
                gd,               &   !< clef d'immobilisation pour le gd
                diseg,            &   !< clef de discretization du segment
                bloquer,          &   !< clef de bloquage des segments sur une barriere
                dedant                !< clef pour les seg. entier dans la boite
end type segment

type(SEGMENT),dimension(nsegmax) :: seg   !< The segments tab

integer(kind=dpi)   :: icristallo   !<

integer(kind=dpi)   :: nb_step_stat !<
integer(kind=dpi)   :: Nplane       !<
integer(kind=dpi)   :: I            !<
integer(kind=dpi)   :: J            !<
integer(kind=dpi)   :: indexB       !<
real(kind=dpi)      :: useless      !<

integer(kind=DPI) :: nbase        !<
integer(kind=DPI) :: k            !<
integer(kind=DPI) :: kk           !<
integer(kind=DPI) :: kkk          !<
integer(kind=DPI) :: debutaffich  !<
integer(kind=DPI) :: freqaffich   !<
integer(kind=DPI) :: Nnodes       !<
integer(kind=DPI) :: slipsys      !< Considered slip system
integer(kind=DPI) :: KODERR       !<
integer(kind=DPI) :: nsegm        !< Number of segments
integer(kind=DPI) :: nsegmaxf     !<
integer(kind=DPI) :: finaffich    !<
integer(kind=DPI) :: NR           !< Number of line in input file
integer(kind=DPI) :: nline        !< Line counter

integer :: status_input !< Existence of the input file

integer(kind=DPI),Dimension(3) :: modur     !<
integer(kind=DPI),Dimension(3) :: BoxCenter !<
integer(kind=DPI),Dimension(3) :: NRep(1:3) !<

integer(kind=DPI), dimension(5,nsegmax)   :: iseg               !< The segments information
integer(kind=DPI), dimension(2,nsegmax)   :: TABVOIS            !< The segment first neighbors
integer(kind=DPI), dimension(nsegmax)     :: junctt             !< The junction pair-segment
integer(kind=DPI), dimension(nsegmax)     :: ID_surface         !< the segments surface ID

integer(kind=DPI),dimension(3,96) :: bveclin       !<
integer(kind=DPI),dimension(3,96) :: bvecnor

real(kind=dp)  :: GridDist       !<
real(kind=dp)  :: GridNormalNorm !<
real(kind=dp)  :: widthp         !<
real(kind=dp)  :: coeftaille     !<

real(kind=dp),dimension(96)    :: normlin !<

real(kind=DP) :: PlaneSp      !<
real(kind=DP) :: bdiva        !<
real(kind=DP) :: bdivpa       !<
real(kind=DP) :: VecBurgers   !<
real(kind=DP) :: dpoiss       !<
real(kind=DP) :: avalue       !<
real(kind=DP) :: epso         !<
real(kind=DP) :: deltat       !<
real(kind=DP) :: mu           !<
real(kind=DP) :: rr1          !<
real(kind=DP) :: rr2          !<
real(kind=DP) :: rr3          !<
real(kind=DP) :: norm         !<
real(kind=DP) :: xstress      !<
real(kind=DP) :: core_radius  !< Core radius spreading factor used in short range stress calculation

real(kind=DP),Dimension(3)  :: X            !<
real(kind=DP),Dimension(3)  :: Y            !<
real(kind=DP),Dimension(3)  :: GridNormal   !<
real(kind=dp),dimension(3)  :: r            !< Coordinate of the stress calculation
real(kind=dp),dimension(6)  :: dombounds    !< Boundaries if restriction on the exploration domain
real(kind=DP),Dimension(3)  :: MapCenter    !<
real(kind=DP),Dimension(3)  :: maxbounds     !< Maximum value in X, Y and Z direction
real(kind=DP),Dimension(3)  :: minbounds     !< Minimum value in X, Y and Z direction

real(kind=DP),Dimension(6):: sigtmp         !<
real(kind=DP),Dimension(6):: sigtmprep     !<

real(kind=DP),Dimension(3,nsegmax)  :: T      !< vecteur ligne normalise
real(kind=DP),Dimension(3,nsegmax)  :: B      !< vecteur de Burgers normalise
real(kind=DP),Dimension(3,nsegmax)  :: BPVT   !< vecteur normale a B et T

real(kind=DP),dimension(6,100,100)   :: sigint      !<
real(kind=DP),dimension(6,100,100)   :: sigintrep   !<
real(kind=DP),dimension(1000,3)   :: mycoords         !<  Store the X,Y,Z coordinates from a file ( mydata('Number of lines',3) )

logical                     :: Film         !< True if Film.bin is used
logical                     :: IsSingular   !< True if we want the microMegas singular field | False allows to choose core radius spreading factor
logical                     :: Islimits     !< True if we want limit the exploration to a subdomain
logical                     :: StationFile     !< True if we want to use a file containing X,Y,Z Coordinates

logical,dimension(nsegmax) :: hide !<
character           :: caractest  !<
character(len=1)    :: carac      !<
character(len=3)    :: crystal_structure !<
character(len=70)   :: filename   !< The results file name
character(len=50)   :: fichier    !<
character(len=60)   :: fichbase   !<
character(len=60)   :: materiau   !<
character(len=60)   :: control    !<
character(len=60)   :: segments   !<
character(len=60)   :: cristallo  !<
character(len=256)  :: command    !< Unix Command used in subroutine system
character(len=256)  :: filepath   !< File path
character(len=8)    :: date       !<
character(len=10)   :: time       !<
character(len=15)   :: mydate     !<
character(len=15)   :: core_radius_str !< String to store core radius for write statements

Print *, " "
Print *, " "
Print *, " "
Print *, "                                    MM    M  EEEE   GGGG    AAA    SSSS"
Print *, "                                    M M M M  EE    G       A   A  SS   "
Print *, "         MM MM  I  CCC  RRRR  OOOO  M  M  M  EEEE  G  GGG  AAAAA   SSS "
Print *, "         M M M  I  C    RRR   O  0  M     M  EE    GG  GG  A   A     SS"
Print *, "         M   M  I  CCC  R  R  0000  M     M  EEEE   GGGG   A   A  SSSS "
Print *, " "
Print *, "                                             --> histo    "
Print *, " "
Print *, " ========================================================================================"
Print *, " "


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (5)                                                   !
  ! The dislocation microstructure (segment output file) is used      !
  ! to make a map of internal stress at the nodes of a plane grid.    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Whatever append before, deallocate variables
  !deallocate(bveclin)
  !deallocate(bvecnor)
  !deallocate(seg)
  !deallocate(hide)
  !deallocate(normlin)

  ! The simulation parameters

  open(1,file='../../in/input.dd',STATUS='OLD')

  carac = " "
  do while (carac /= "-")
    read (1,*) carac
  enddo

  read (1,*) materiau
  read (1,*) control
  read (1,*) segments

  close (1)

  ! Identification of the crystal symetry
  open(1,file="../../in/"//materiau,STATUS='OLD')

  ! some useless data are by passed first
  read(1,*) useless
  read(1,*) useless
  read(1,*) mu        ! The shear modulus
  mu = mu * 1.D9      ! now mu is in Pa
  read(1,*) useless
  read(1,*) dpoiss    ! Poisson ratio
  read(1,*) useless
  read(1,*) useless
  read(1,*) useless
  read(1,*) crystal_structure(1:3) ! the crystal symetry
  read(1,*) useless
  read(1,*) VecBurgers    ! The burgers vector in Angstrom
  VecBurgers = VecBurgers * 1.d-10
  close (1)

  ! Choice between Seg_save and film.bin
  write(*,*) "Use Film.bin (F = SEGSAVE) ? "
  read (*,*) Film

  if (Film) then

    print *, " > We use film.bin informations "
    ! The file we want to analyze
    fichier = "film.bin"
    open(50,file='../../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
    ! Important information at the begining of the film.bin file
    read(50) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
    read(50) nsegmaxf     ! The nsegmax value used in this simulation
    read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
    read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry

    !allocate (bvecnor(3,nbase))
    !allocate (bveclin(3,nbase))

    read(50) (bveclin(:,I),bvecnor(:,I), I=1,nbase)    ! table of vector direction and displacement for this simulation
    read(50) modur(:)     ! Dimension of the 3D simulated volume in a value

    ! Test on the nsegmax definition
    if (nsegmax /= nsegmaxf) then
      write(*,*) "Number of segments in film.bin = ",nsegmaxf,"but in this programm = ",nsegmax
      write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
      stop
    endif

  else

    print *, " > We use SEGSAVE informations "

    ! The simulation elementary lattice parameter
    cristallo = '../../out/BVD.'//crystal_structure(1:3)

    open (1,FILE=cristallo,STATUS='OLD')

    read (1,*) nbase , avalue
    !allocate (bveclin(3,nbase))

    do j = 1,nbase
        read (1,*) i, bveclin(1:3,j)
        if(j /= i) print *, "stop, a problem appears when reading BVD"
    enddo

    close(1)

    open(50,file='../../in/SEG_save',STATUS='OLD')

    read(50,*) useless                     ! Parameter of effective loading on the systems
    read(50,*) Nsegm                       ! Nb of segments

    ! Test on the nsegmax definition
    if (nsegm > nsegmax) then
      write(*,*) "Number of segments in SEGSAVE = ",nsegmaxf,"but in this program = ",nsegmax
        write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
        stop
    endif

    read(50,*) modur(:)  ! Dimensions of the simulation cell

    !allocate (seg(nsegm))
    !allocate (hide(nsegm))

    do I = 1,nsegm

      read (50,*) j,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
                 seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%Ijonc,seg(i)%gd

      ! Identification of pinning point is made in conformity with the mM simulation rule
      if(seg(i)%vnno  == 0) seg(i)%vnno  = nsegmax
      if(seg(i)%voiso == 0) seg(i)%voiso = nsegmax
      if(seg(i)%vnne  == 0) seg(i)%vnne  = nsegmax
      if(seg(i)%voise == 0) seg(i)%voise = nsegmax
      if(seg(i)%Ijonc == 0) seg(i)%ijonc = nsegmax

      if(j /= i) then
        print *,"A problem appears when loading SEG_save"
        stop
      endif

    enddo

  endif

  ! The center of the simulation cell
  BoxCenter(1:3) = modur(1:3)/2

  !allocate (normlin(nbase))

  do j = 1,nbase
      normlin(j)=sqrt(real(bveclin(1,j))**2+  &
                      real(bveclin(2,j))**2+  &
                      real(bveclin(3,j))**2)
  enddo

  nb_step_stat = 0

  If (Film) then
    ! We first select the part of the film.bin file we want to use for the stats
    write(*,*) "Statistics start at simulation step:  ? "
    read (*,*) debutaffich
    write(*,*) "Statistics end at simulation step:  ? "
    read (*,*)  finaffich
    write(*,*) "Snapshot frequency for the film analysis: (Warning : must be coherent with KSTATS) ? "
    read (*,*) freqaffich
    if (freqaffich  <= 0)freqaffich = 1
  else !set 1 for all parameters in case of SEG_save
    debutaffich =1
    finaffich   =1
    freqaffich  =1
  endif

  ! Choice between Seg_save and film.bin
  write(*,*) "Use a file containing point coordinates (F = Node mapping) ? "
  read (*,*) StationFile

  if (.not. StationFile) then
    ! We load information about the stress Grid
    write(*,*) "How nany nodes per side shall we take for the squared stress grid "
    read(*,*) Nnodes

    DO WHILE (Nnodes < 1 .or. Nnodes > 500 )
      write(*,*) "<! Please consider a number of nodes between 5 and 500 !> "
      write(*,*) "How nany nodes per side shall we take for the squared stress grid "
      read(*,*) Nnodes
    ENDDO

    write(*,*) "Aspect ratio between stress mapping and the simulation cell "
    read(*,*) coeftaille

    minbounds(:) = 0.

    maxbounds(1)=modur(1)
    maxbounds(2)=modur(2)
    maxbounds(3)=modur(3)

  else
    DO WHILE (status_input/=0)
      NR=-1
      status_input=1

      filepath='../../out/location.dat'
      call get_file_line_number(filepath,status_input,NR) !Read the station file
    ENDDO

    OPEN(90,file=filepath,STATUS='OLD')
    !Now we can allocate data variables
    !ALLOCATE(mycoords(NR,3))

    !Now read data into mydata
    DO J=1,NR
      READ(90,*) mycoords(J,1) ,  mycoords(J,2), mycoords(J,3)
    ENDDO

    CLOSE(90)
    Nnodes = int(sqrt(real(NR)))
    maxbounds(:)= MAXVAL (mycoords, DIM=1)
    minbounds(:)= MINVAL (mycoords, DIM=1)
  endif

  write(*,*) "Number of replicates per side of the reference volume Nx, Ny, Nz ? "
  read(*,*) NRep(1:3)

  write(*,*) "Please enter the grid normal direction, x, y, z "
  read(*,*) GridNormal(1:3)

  GridNormalNorm = sqrt(GridNormal(1)**2+GridNormal(2)**2+GridNormal(3)**2)
  GridNormal(1:3) = GridNormal(1:3)/GridNormalNorm

  Nplane=0
  DO WHILE (Nplane < 1)
    write(*,*) "Enter the number of planes you want to analyze. "
    read(*,*) Nplane
  ENDDO

  if (Nplane > 1) then
    write(*,*) "Enter the spacing between the planes in microns. "
    read(*,*) PlaneSp
    PlaneSp = PlaneSp * 1.D-6 / avalue
  else
    PlaneSp = 0
  endif

  write(*,*) "Enter the normal distance from the top most plane to the simulation cell center (microns) "
  read(*,*) GridDist
  GridDist = GridDist * 1.D-6 / avalue


  !slipsys
  write (*,*) "Which slip system do you want to consider (>12 all slip systems) ? "
  read (*,*) slipsys

  DO WHILE ((slipsys <=0))
    write (*,*) '!> Slip system number must be greater than zero <!'
    write (*,*) "Slip system you want to compute curvature radius vector (>12 all slip systems)"
    read  (*,*) slipsys
  END DO

  write(*,*) "Do you want the singular stress field (T/F) ? [T = Singular stress field | F = Regularized stress field] "
  read(*,*) IsSingular

  If  (IsSingular) then
    core_radius = 0.5 !This Value should be the same as in microMegas
  else
    write(*,*) "Which core radius spreading factor for short range stress calculation ('halfthickness' in DCM) ? "
    read(*,*) core_radius
  endif

  write(*,fmt='("Do you want to explore only a subdomain (T/F) ?  [Default : ",6I9," ] ")') &
  & int(minbounds(1)),int(maxbounds(1)),int(minbounds(2)),int(maxbounds(2)),int(minbounds(3)),int(maxbounds(3))
  read(*,*) Islimits

  If (Islimits) then

  dombounds(:)=-100.

    DO WHILE (dombounds(1) < 0 .or. dombounds(1) > modur(1))
      write (*,*) ' xmin  >'
      read  (*,*) dombounds(1)
    END DO

    DO WHILE (dombounds(2) < 0 .or. dombounds(2) > modur(1) .or. dombounds(2) < dombounds(1))
      write (*,*) ' xmax >'
      read  (*,*) dombounds(2)
    END DO

    DO WHILE (dombounds(3) < 0 .or. dombounds(3) > modur(2))
      write (*,*) ' ymin >'
      read  (*,*) dombounds(3)
    END DO

    DO WHILE (dombounds(4) < 0 .or. dombounds(4) > modur(2) .or. dombounds(4) < dombounds(3))
      write (*,*) ' ymax >'
      read  (*,*) dombounds(4)
    END DO

    DO WHILE (dombounds(5) < 0 .or. dombounds(5) > modur(3))
      write (*,*) ' zmin >'
      read  (*,*) dombounds(5)
    END DO

    DO WHILE (dombounds(6) < 0 .or. dombounds(6) > modur(3) .or. dombounds(6) < dombounds(5))
      write (*,*) ' zmax >'
      read  (*,*) dombounds(6)
    END DO

  else

  dombounds(1)= real(minbounds(1),DP)
  dombounds(2)= real(maxbounds(1),DP)
  dombounds(3)= real(minbounds(2),DP)
  dombounds(4)= real(maxbounds(2),DP)
  dombounds(5)= real(minbounds(3),DP)
  dombounds(6)= real(maxbounds(3),DP)

  endif

  if (Film) then
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

  else
    fichbase = crystal_structure
  endif

  print *, "--"
  print *, "  The crystal symmetry is of type ",fichbase
  print *, " "

  ! sigint: only segments on the box are used to caluclate stress

  !allocate (sigint(3,3,Nnodes,Nnodes))
  !allocate (sigintrep(3,3,Nnodes,Nnodes))

  ! Setup the result directory
  call date_and_time(date,time)!,zone,values)
  mydate=date//'_'//time(1:6)

  command = 'mkdir ../../out/SigInt_dir'//'_'//mydate
  call system(command)

  ! WE START THE LOOP ON SIMULATION STEPS SAVED IN FILM.BIN

  print *, "-------------------------"
  print *, " Beginning of step loop "
  print *, "-------------------------"
  print *, ""

  ! INITIALIZATION ON THE LOOP ON THE SELECTED STEP
  kk=0

  Steps_loop_4: do while (kk.ge.0 .and. kk .le. finaffich)

    sigint(:,:,:)     = 0. !Initialization
    sigintrep(:,:,:)  =  0.

    if (Film) then
      ! A new simulation step
      read(50,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO

      if(kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0) then
        print*," > Selected step "
        write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I5)') KK,EPSO,XSTRESS, nsegm
      endif

      if(KODERR < 0) then
        print*,"  End of the film.bin file "
        exit Steps_loop_4
      elseif(KODERR > 0) then
        print*,"Reading problem 1",KODERR,kk
        stop
        !cycle Steps_loop_4
      endif

      ! The segments micro structure is loaded
      READ(50,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),ID_surface(j), J=1,nsegm)

      if(KODERR < 0) then
        print*,"  End of the film.bin file"
        exit Steps_loop_4
      elseif(KODERR > 0) then
        print*,"Reading Error 2",KODERR,nsegm,j
        stop
        !cycle Steps_loop_4
      endif
    else
      kk=1 !Segsave
    endif

    ! The steps we want to use for statistic
    if(kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0) then ! Begin of the "Big If"

      if (Film) then

        !nsegm (number of segments) might be different for each step
        !allocate (seg(nsegm))
        !allocate (hide(nsegm))

        ! The film.bin information is used
        do j=1,nsegm

          seg(j)%O(1:3) = Iseg(1:3,j)
          seg(j)%norme  = Iseg(4,j)
          seg(j)%veclin = Iseg(5,j)

          seg(j)%vnno   = tabvois(1,j)
          seg(j)%vnne   = tabvois(2,j)

          seg(j)%jonc   = .false.
          seg(j)%Ijonc  = nsegmax
          seg(j)%gd     = izero

          if (Junctt(j) < 0) then
            seg(j)%gd     = iun
          elseif (Junctt(j) > 0) then
            seg(j)%jonc   = .true.
            seg(j)%Ijonc  = Junctt(j)
          endif
       enddo

        print*," > film.bin informations loaded !"

      endif

      ! Start of the calculations

      !Initialization
      hide(1:nsegm) = .false.

      !Important segments information is predefined
      do K = 1,nsegm
        ! rotule segments can be eliminated from the calculations
        if (seg(k)%norme == 0) hide(k)=.true.

        ! The segment vector direction
        T(1:3,K) = bveclin(1:3,seg(k)%veclin)/normlin(seg(k)%veclin)
        ! The Burgers vector of i
        indexB = (int((seg(k)%veclin-1)/8))*8 + 1
        B(1:3,K) = bveclin(1:3,indexB)/normlin(indexB)

        !B produit vectoriel T
        bpvt(1,k) = b(2,k)*t(3,k)-b(3,k)*t(2,k)
        bpvt(2,k) = b(3,k)*t(1,k)-b(1,k)*t(3,k)
        bpvt(3,k) = b(1,k)*t(2,k)-b(2,k)*t(1,k)
      enddo

      !1) CARTOGRAPHY AXIS DETERMINATION
      ! A direct frame is built starting from the thinfoil normal
      if(GridNormal(1)+GridNormal(2) == 0.) then
        ! first axis: in plan (0y0)
        X(1) = Gridnormal(3)
        X(2) = 0
        X(3) = -Gridnormal(1)
        ! second axis normal to X and N
        Y(1) = Gridnormal(2)*X(3)-Gridnormal(3)*X(2)
        Y(2) = Gridnormal(3)*X(1)-Gridnormal(1)*X(3)
        Y(3) = Gridnormal(1)*X(2)-Gridnormal(2)*X(1)
        ! axis are normalized
        norm = sqrt(X(1)**2+ X(2)**2+ X(3)**2)
        X(:) = X(:) / norm
        norm = sqrt(Y(1)**2+ Y(2)**2+ Y(3)**2)
        Y(:) = Y(:) / norm
      else
        ! first axis: in plan (00z)
        X(1) = Gridnormal(2)
        X(2) = - Gridnormal(1)
        X(3) = 0
        ! second axis normal to X and N
        Y(1) = Gridnormal(2)*X(3)-Gridnormal(3)*X(2)
        Y(2) = Gridnormal(3)*X(1)-Gridnormal(1)*X(3)
        Y(3) = Gridnormal(1)*X(2)-Gridnormal(2)*X(1)
        ! axis are normalized
        norm = sqrt(X(1)**2+ X(2)**2+ X(3)**2)
        X(:) = X(:) / norm
        norm = sqrt(Y(1)**2+ Y(2)**2+ Y(3)**2)
        Y(:) = Y(:) / norm
      endif

      ! Mapping increment
      IF (.not. StationFile) widthp= coeftaille * MAX(Modur(1),Modur(2),Modur(3)) / real(Nnodes)

      ! Calculations can be done on multiple (Nplane) planes automatically.
      do kkk=1,Nplane

        ! The Center of the stress mapping
        MapCenter(1) = BoxCenter(1)+GridNormal(1)*(GridDist-(PlaneSp*(kkk-1)))
        MapCenter(2) = BoxCenter(2)+GridNormal(2)*(GridDist-(PlaneSp*(kkk-1)))
        MapCenter(3) = BoxCenter(3)+GridNormal(3)*(GridDist-(PlaneSp*(kkk-1)))

        If (StationFile) nline=1

        do i=1,Nnodes
            print *,"   > Computing line : ", i
            If (.not. StationFile) then
              rr1 = mapcenter(1) + (i-real(Nnodes,DP)*0.5 - 0.5)*X(1)*widthp
              rr2 = mapcenter(2) + (i-real(Nnodes,DP)*0.5 - 0.5)*X(2)*widthp
              rr3 = mapcenter(3) + (i-real(Nnodes,DP)*0.5 - 0.5)*X(3)*widthp
            endif

          do j=1,Nnodes
            If (.not. StationFile) then
              r(1) = rr1 + (j-real(Nnodes,DP)*0.5 - 0.5)*Y(1)*widthp
              r(2) = rr2 + (j-real(Nnodes,DP)*0.5 - 0.5)*Y(2)*widthp
              r(3) = rr3 + (j-real(Nnodes,DP)*0.5 - 0.5)*Y(3)*widthp
            else
              r(1) = mycoords(nline,1)
              r(2) = mycoords(nline,2)
              r(3) = mycoords(nline,3)
              nline=nline+1 ! increment the number of lines
            endif
            if ((r(1) >= dombounds(1)) .and. (r(2) >= dombounds(3)) .and. (r(3) >= dombounds(5)) .and.  &
              & (r(1) <= dombounds(2)) .and. (r(2) <= dombounds(4)) .and. (r(3) <= dombounds(6))) then

              call sigma_int(nsegm,r,sigtmp,sigtmprep,NRep) ! Calculation of the stress tensor at coordinate r

              sigint(:,i,j) = sigint(:,i,j) + sigtmp(:)
              sigintrep(:,i,j) = sigintrep(:,i,j) + sigtmprep(:)

            else

              sigint(:,i,j) = 0.0
              sigintrep(:,i,j) = 0.0

            endif

          end do

        end do

      end do   ! Nplane


      ! Management of the calculation outputs
      nb_step_stat = nb_step_stat +1                    ! The filename number

    OPEN (7,STATUS = 'SCRATCH') ! Get core radius as a string
    WRITE(7,*) core_radius
    REWIND(7)
    READ(7,*) core_radius_str
    CLOSE(7)

  114 format ('../../out/SigInt_dir','_',A15,'/sigint_h',A6,'_incr',I0, '.txt')         ! The file name construction

      write (filename, 114) mydate, core_radius_str, nb_step_stat                ! definition of a different filename for each loop
      open (unit=601, file=filename, status='unknown')  ! The file name is open
      write(601,'(2x,A2,2x,A2,6(2x,A5))') 'X','Y','Sig11','Sig22','Sig33',&
      'Sig23','Sig13','Sig12'

      do i = 1, Nnodes
        do j = 1, Nnodes
          write(601,*) i, j, sigint(1,i,j) * mu * 1D-6 / Nplane,    &
                             sigint(2,i,j) * mu * 1D-6 / Nplane,    &
                             sigint(3,i,j) * mu * 1D-6 / Nplane,    &
                             sigint(4,i,j) * mu * 1D-6 / Nplane,    &
                             sigint(5,i,j) * mu * 1D-6 / Nplane,    &
                             sigint(6,i,j) * mu * 1D-6 / Nplane
        end do
      end do
      close(601)

  115 format ('../../out/SigInt_dir','_',A15,'/sigintrep_h',A6,'_incr',I0, '.txt')         ! The file name construction

      write (filename, 115) mydate, core_radius_str, nb_step_stat                ! definition of a different filename for each loop
      open (unit=602, file=filename, status='unknown')  ! The file name is open
      write(602,'(2x,A2,2x,A2,6(2x,A5))') 'X','Y','Sig11','Sig22','Sig33',&
      'Sig23','Sig13','Sig12'

      do i = 1, Nnodes
        do j = 1, Nnodes
          write(602,*) i, j, sigintrep(1,i,j) * mu * 1D-6 / Nplane,    &
                             sigintrep(2,i,j) * mu * 1D-6 / Nplane,    &
                             sigintrep(3,i,j) * mu * 1D-6 / Nplane,    &
                             sigintrep(4,i,j) * mu * 1D-6 / Nplane,    &
                             sigintrep(5,i,j) * mu * 1D-6 / Nplane,    &
                             sigintrep(6,i,j) * mu * 1D-6 / Nplane
        end do
      end do
      close(602)

      ! END OF THE CALCULATIONS MADE ON THIS STEP
      !  ======>>>>>

      !deallocate (seg)
      !deallocate (hide)
    endif ! End of the "big if" selecting the step

    If (Film) then
      ! The mark for the end of step is read
      read(50,IOSTAT=KODERR) caractest

      if(KODERR < 0) then
          print*,"  End of the film.bin file 3"
          exit Steps_loop_4
      elseif(KODERR > 0) then
          print*,"Reading error - Could not find the mark for the end of step ",caractest
          stop
          !cycle Steps_loop_4
      endif

    else
      kk=kk+1
    endif

  enddo Steps_loop_4

  print *, ""
  print *, "-------------------------"
  print *, " End of step loop "
  print *, "-------------------------"
  print *, ""

  print *,' > Results are in /out/SigInt_dir'

  close(50)

  print*," > Internal stress field calculation achieved with success"

   !If (StationFile) deallocate(mycoords)

CONTAINS

!##########################################################################
!# Calculation of the internal stress associate to a dislocation          #
!# microstructure with nsegm segments at point r                          #
!##########################################################################
subroutine sigma_int(nsegm,r,sigtmp,sigtmprep,NRep)

implicit none

integer (kind=dpi),intent(in)            :: nsegm             !< Number of segments
real(kind=dp),dimension(3),intent(in)    :: r                 !< Coordinate of the current calculation
real(kind=dp),dimension(6),intent(out) :: sigtmp            !< The stress tensor at r
real(kind=dp),dimension(6),intent(out) :: sigtmprep         !< The stress tensor at r
real(kind=dp)                            :: core2             !< (Core radius / avalue)^2
!real(kind=dp),dimension(6)         :: sigint
!real(kind=dp),dimension(6)         :: sigintrep
!real(kind=dp),dimension(6,nsegm)         :: sigint2
!real(kind=dp),dimension(6,nsegm)         :: sigintrep2
real(kind=dp),dimension(3,nsegm)         :: ra,rb, ra0, rb0
integer(kind=dpi)                        :: i,carac,NX,NY,NZ, Nrep(1:3)
integer(kind=dpi),dimension(3)           :: Tr,ab,centrei

real(kind=DP) :: d2,s1,s2,a2_p_d2,inv_a2_p_d2
real(kind=DP) :: s1Ra3inv,Ra1_3inv,Ra1_inv
real(kind=DP) :: s2Ra3inv,Ra2_3inv,Ra2_inv
real(kind=DP) :: s1_03,s1_13,s1_05,s1_15,s1_25
real(kind=DP) :: s2_03,s2_13,s2_05,s2_15,s2_25
real(kind=DP) :: s_03,s_13,s_05,s_15,s_25
real(kind=DP) :: m4pn,mn4pn,a2m4pn,a2m8p
real(kind=DP) :: d_pv_b_ps_t,commun,R_ps_t,Ra2
real(kind=DP),dimension(3) :: cai_d,d_pv_b,t_pv_b
real(kind=DP),dimension(6) :: t_pt_t
real(kind=DP),dimension(6) :: d_pt_d,t_pt_d,t_pt_t_pv_b,d_pt_t_pv_b,t_pt_d_pv_b
real(kind=DP),dimension(6) :: I_03,I_13,I_05,I_15,I_25

integer(kind=dpi)   :: IndexSeg
integer :: threadnum !< Thread number if openMP multithreading

!sigint2(:,:)   =  0.0d0
sigtmp(:) =  0.0d0
!sigintrep2(:,:)   =  0.0d0
sigtmprep(:) =  0.0d0

BdivA = vecburgers / Avalue

core2 = (core_radius * BdivA)*(core_radius * BdivA)

bdivpa = bdiva/pI

m4pn=0.25d0/(1.0d0-dpoiss)
mn4pn=m4pn*dpoiss
a2m4pn=core2*m4pn
a2m8p=core2*0.125d0

! The translation needed for periodic boundaries
! With this translation r is shifted at the center of the simulation box
Tr(1:3)   = BoxCenter(1:3)- nint(r(1:3))

!print *, '#initialisation'
do i = 1,nsegm
  ! Those segments don t need to be considered
  if (hide(i)) cycle
  ab(:) = seg(i)%norme*bveclin(:,seg(i)%veclin)
  centrei(:) = seg(i)%o(:) + ab(:)/2
  centrei(:) = modulo(centrei(:)+Tr(:), ModuR(:))
  ra0(:,i) = real(BoxCenter(:) - (centrei(:) - ab(:)/2),DP)
  rb0(:,i) = real(BoxCenter(:) - (centrei(:) + ab(:)/2),DP)
enddo

!print *,"threadnum :", threadnum

  do NX = -Nrep(1),Nrep(1)
  do NY = -Nrep(2),Nrep(2)
  do NZ = -Nrep(3),Nrep(3)

!---------------------------------------------------------------
! Multi threading work only with nsegmax arround 10000 segments
!---------------------------------------------------------------
!$OMP PARALLEL NUM_THREADS(4) DEFAULT(NONE)  , &
! !$OMP PARALLEL DEFAULT(NONE)  , &
!$OMP& PRIVATE(d2,s1,s2,a2_p_d2,inv_a2_p_d2,s1Ra3inv,Ra1_3inv,Ra1_inv,s2Ra3inv,Ra2_3inv,Ra2_inv) , &
!$OMP& PRIVATE(s1_03,s1_13,s1_05,s1_15,s1_25,s2_03,s2_13,s2_05,s2_15,s2_25,s_03,s_13,s_05,s_15,s_25) , &
!$OMP& PRIVATE(d_pv_b_ps_t,commun,R_ps_t,Ra2,cai_d,d_pv_b,t_pv_b) , &
!$OMP& PRIVATE(t_pt_t,d_pt_d,t_pt_d,t_pt_t_pv_b,d_pt_t_pv_b,t_pt_d_pv_b,I_03,I_13,I_05,I_15,I_25) , &
!$OMP& PRIVATE(NX,NY,NZ,IndexSeg,carac,threadnum) , &
!$OMP& SHARED(ra,rb,nsegm,NREP,slipsys,hide,ra0,rb0,modur) , &
!$OMP& SHARED(T,B,core2,bpvt,bdivpa,seg,m4pn,mn4pn,a2m4pn,a2m8p) , &
!$OMP& REDUCTION(+:sigtmp) , &
!$OMP& REDUCTION(+:sigtmprep)
!$ threadnum = OMP_GET_THREAD_NUM()

!$OMP DO
do i = 1,nsegm

    ! a test to keep only the contribution of one slip system
    if (slipsys < 13) then
       IndexSeg = (int((seg(i)%veclin-1)/8)) + 1
       if (seg(i)%jonc) IndexSeg = 0
       if (IndexSeg /= slipsys) cycle
    endif

    ! Those segments don t need to be considered
    if (hide(i)) cycle

    ! information relative to the segments
    ! The segment type (the tyseg type 1 <-> 4!)
    carac = seg(i)%veclin - (((seg(i)%veclin-1)/(nbasered/2))*(nbasered/2))
    ! The vetors connecting r to O and E of segment i

    ra(1,i) = ra0(1,i) + NX*MODUR(1)
    ra(2,i) = ra0(2,i) + NY*MODUR(2)
    ra(3,i) = ra0(3,i) + NZ*MODUR(3)
    rb(1,i) = rb0(1,i) + NX*MODUR(1)
    rb(2,i) = rb0(2,i) + NY*MODUR(2)
    rb(3,i) = rb0(3,i) + NZ*MODUR(3)

    R_ps_t=ra(1,i)*T(1,i)+ra(2,i)*T(2,i)+ra(3,i)*T(3,i)
    s1=-R_ps_t
    cai_d=RA(:,i)-R_ps_t*T(:,i)
    d2=cai_d(1)*cai_d(1)+cai_d(2)*cai_d(2)+cai_d(3)*cai_d(3) ! distance au segment
    R_ps_t=rb(1,i)*T(1,i)+rb(2,i)*T(2,i)+rb(3,i)*T(3,i)
    s2=-R_ps_t

    a2_p_d2=core2+d2
    inv_a2_p_d2=1.d0/a2_p_d2

    Ra2=a2_p_d2+s1*s1
    Ra1_inv=1.d0/sqrt(Ra2)
    Ra1_3inv=Ra1_inv/Ra2 ! 1/ + cher que 2* ? OUI
    s1Ra3inv=s1*Ra1_3inv

    Ra2=a2_p_d2+s2*s2
    Ra2_inv=1.d0/sqrt(Ra2)
    Ra2_3inv=Ra2_inv/Ra2 ! 1/ + cher que 2* ? ENFIN PAREIL (dans la barre d'erreur)
    s2Ra3inv=s2*Ra2_3inv

    s1_03=s1*Ra1_inv*inv_a2_p_d2
    s2_03=s2*Ra2_inv*inv_a2_p_d2
    s_03=s2_03-s1_03

    s1_13=-Ra1_inv
    s2_13=-Ra2_inv
    s_13=s2_13-s1_13

    s1_05=(2.d0*s1_03+s1Ra3inv)*inv_a2_p_d2
    s2_05=(2.d0*s2_03+s2Ra3inv)*inv_a2_p_d2
    s_05=s2_05-s1_05

    s1_15=-Ra1_3inv
    s2_15=-Ra2_3inv
    s_15=s2_15-s1_15

    s1_25=s1_03-s1Ra3inv
    s2_25=s2_03-s2Ra3inv
    s_25=s2_25-s1_25

    t_pv_b=-BPVT(:,i)

    d_pv_b(1)=cai_d(2)*B(3,i)-cai_d(3)*B(2,i)
    d_pv_b(2)=cai_d(3)*B(1,i)-cai_d(1)*B(3,i)
    d_pv_b(3)=cai_d(1)*B(2,i)-cai_d(2)*B(1,i)

    d_pv_b_ps_t=d_pv_b(1)*T(1,i)+d_pv_b(2)*T(2,i)+d_pv_b(3)*T(3,i)

    d_pt_d(1) = cai_d(1)*cai_d(1)
    d_pt_d(2) = cai_d(2)*cai_d(2)
    d_pt_d(3) = cai_d(3)*cai_d(3)
    d_pt_d(4) = cai_d(2)*cai_d(3)
    d_pt_d(5) = cai_d(1)*cai_d(3)
    d_pt_d(6) = cai_d(1)*cai_d(2)

    t_pt_t(1) = T(1,i)*T(1,i)
    t_pt_t(2) = T(2,i)*T(2,i)
    t_pt_t(3) = T(3,i)*T(3,i)
    t_pt_t(4) = T(2,i)*T(3,i)
    t_pt_t(5) = T(1,i)*T(3,i)
    t_pt_t(6) = T(1,i)*T(2,i)

    t_pt_d(1) = (T(1,i)*cai_d(1))*2.d0
    t_pt_d(2) = (T(2,i)*cai_d(2))*2.d0
    t_pt_d(3) = (T(3,i)*cai_d(3))*2.d0
    t_pt_d(4) = T(2,i)*cai_d(3)+cai_d(2)*T(3,i)
    t_pt_d(5) = T(1,i)*cai_d(3)+cai_d(1)*T(3,i)
    t_pt_d(6) = T(1,i)*cai_d(2)+cai_d(1)*T(2,i)

    t_pt_t_pv_b(1) = (T(1,i)*t_pv_b(1))*2.d0
    t_pt_t_pv_b(2) = (T(2,i)*t_pv_b(2))*2.d0
    t_pt_t_pv_b(3) = (T(3,i)*t_pv_b(3))*2.d0
    t_pt_t_pv_b(4) = T(2,i)*t_pv_b(3)+t_pv_b(2)*T(3,i)
    t_pt_t_pv_b(5) = T(1,i)*t_pv_b(3)+t_pv_b(1)*T(3,i)
    t_pt_t_pv_b(6) = T(1,i)*t_pv_b(2)+t_pv_b(1)*T(2,i)

    d_pt_t_pv_b(1) = (cai_d(1)*t_pv_b(1))*2.d0
    d_pt_t_pv_b(2) = (cai_d(2)*t_pv_b(2))*2.d0
    d_pt_t_pv_b(3) = (cai_d(3)*t_pv_b(3))*2.d0
    d_pt_t_pv_b(4) = cai_d(2)*t_pv_b(3)+t_pv_b(2)*cai_d(3)
    d_pt_t_pv_b(5) = cai_d(1)*t_pv_b(3)+t_pv_b(1)*cai_d(3)
    d_pt_t_pv_b(6) = cai_d(1)*t_pv_b(2)+t_pv_b(1)*cai_d(2)

    t_pt_d_pv_b(1) = (T(1,i)*d_pv_b(1))*2.d0
    t_pt_d_pv_b(2) = (T(2,i)*d_pv_b(2))*2.d0
    t_pt_d_pv_b(3) = (T(3,i)*d_pv_b(3))*2.d0
    t_pt_d_pv_b(4) = T(2,i)*d_pv_b(3)+d_pv_b(2)*T(3,i)
    t_pt_d_pv_b(5) = T(1,i)*d_pv_b(3)+d_pv_b(1)*T(3,i)
    t_pt_d_pv_b(6) = T(1,i)*d_pv_b(2)+d_pv_b(1)*T(2,i)

    commun=m4pn*d_pv_b_ps_t

    I_03(1)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(1))-0.25d0*t_pt_d_pv_b(1)
    I_03(2)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(2))-0.25d0*t_pt_d_pv_b(2)
    I_03(3)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(3))-0.25d0*t_pt_d_pv_b(3)
    I_03(4)=m4pn*d_pt_t_pv_b(4)-0.25d0*t_pt_d_pv_b(4)
    I_03(5)=m4pn*d_pt_t_pv_b(5)-0.25d0*t_pt_d_pv_b(5)
    I_03(6)=m4pn*d_pt_t_pv_b(6)-0.25d0*t_pt_d_pv_b(6)

    if (carac==1) then !test if screw segment | quicker calculation
      if (nx/=0.or.ny/=0.or.nz/=0)then
        sigtmprep(:)=sigtmprep(:)+ I_03*s_03*bdivpa
        !sigintrep2(:,I) = I_03*s_03*bdivpa
      else
        sigtmp(:)=sigtmp(:)+ I_03*s_03*bdivpa
        sigtmprep(:)=sigtmprep(:)+ I_03*s_03*bdivpa
        !sigint2(:,I) = I_03*s_03*bdivpa
        !sigintrep2(:,I) = I_03*s_03*bdivpa
      endif
      cycle
    endif

    I_13=-mn4pn*t_pt_t_pv_b(:)

    I_15(1)=a2m8p*t_pt_t_pv_b(1)-commun*t_pt_d(1)
    I_15(2)=a2m8p*t_pt_t_pv_b(2)-commun*t_pt_d(2)
    I_15(3)=a2m8p*t_pt_t_pv_b(3)-commun*t_pt_d(3)
    I_15(4)=a2m8p*t_pt_t_pv_b(4)-commun*t_pt_d(4)
    I_15(5)=a2m8p*t_pt_t_pv_b(5)-commun*t_pt_d(5)
    I_15(6)=a2m8p*t_pt_t_pv_b(6)-commun*t_pt_d(6)

    I_05(1)=commun*(core2+d_pt_d(1))-a2m8p*t_pt_d_pv_b(1)
    I_05(2)=commun*(core2+d_pt_d(2))-a2m8p*t_pt_d_pv_b(2)
    I_05(3)=commun*(core2+d_pt_d(3))-a2m8p*t_pt_d_pv_b(3)
    I_05(4)=commun*d_pt_d(4)-a2m8p*t_pt_d_pv_b(4)
    I_05(5)=commun*d_pt_d(5)-a2m8p*t_pt_d_pv_b(5)
    I_05(6)=commun*d_pt_d(6)-a2m8p*t_pt_d_pv_b(6)

    I_25=commun*t_pt_t

    if (nx/=0.or.ny/=0.or.nz/=0)then !If any Replica
      sigtmprep(:)=sigtmprep(:)+ (I_03*s_03+I_13*s_13+I_05*s_05+I_15*s_15+I_25*s_25)*bdivpa
      !sigintrep2(:,I) = (I_03*s_03+I_13*s_13+I_05*s_05+I_15*s_15+I_25*s_25)*bdivpa
    else
      sigtmp(:)=sigtmp(:)+ (I_03*s_03+I_13*s_13+I_05*s_05+I_15*s_15+I_25*s_25)*bdivpa
      sigtmprep(:)=sigtmprep(:)+ (I_03*s_03+I_13*s_13+I_05*s_05+I_15*s_15+I_25*s_25)*bdivpa
      !sigint2(:,I) = (I_03*s_03+I_13*s_13+I_05*s_05+I_15*s_15+I_25*s_25)*bdivpa
      !sigintrep2(:,I) =  (I_03*s_03+I_13*s_13+I_05*s_05+I_15*s_15+I_25*s_25)*bdivpa
    endif


enddo !i = 1,nsegm

!$OMP END DO
!$OMP END PARALLEL


  enddo !NVOL
  enddo !NVOL
  enddo !NVOL

!sigtmp(6) = SUM (sigint2(6,:))
!sigtmp(4) = SUM (sigint2(4,:))
!sigtmp(5) = SUM (sigint2(5,:))
!sigtmp(1) = SUM (sigint2(1,:))
!sigtmp(2) = SUM (sigint2(2,:))
!sigtmp(3) = SUM (sigint2(3,:))

!sigtmprep(6) = SUM (sigintrep2(6,:))
!sigtmprep(4) = SUM (sigintrep2(4,:))
!sigtmprep(5) = SUM (sigintrep2(5,:))
!sigtmprep(1) = SUM (sigintrep2(1,:))
!sigtmprep(2) = SUM (sigintrep2(2,:))
!sigtmprep(3) = SUM (sigintrep2(3,:))

end subroutine sigma_int

subroutine get_file_line_number(filepath,status_input,NR)

implicit none

integer             :: status_input !< Existence of the input file
integer(kind=DPI)   :: ios  !< IOSTAT return
integer(kind=DPI)   :: line         !< Line counter
integer(kind=DPI)   :: NR           !< Number of line in input file
character(len=256)  :: filepath     !< File path
character(LEN=1) :: junk !< junk variable

  status_input=access(filepath,'r')

  write(*,*) " > File path ", filepath
  write(*,*) " > Input statut " ,status_input
  IF (status_input == 0) THEN
    OPEN(89,file=filepath,STATUS='OLD')
    NR=0
    DO line=1,1000000
      READ(89,*,IOSTAT=ios) junk
      IF(ios/=0)EXIT
      IF(line==1000000)THEN
        write(*,*) "Error: Maximum number of records exceeded, File has more than 1000000 lines"
        write(*,*) "Exiting program now..."
        STOP
      ENDIF
      NR=NR+1
    ENDDO
    CLOSE(89)
  write(*,*) " > Number of lines ", NR
  ENDIF

end subroutine get_file_line_number

END PROGRAM intstress


