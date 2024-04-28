include '../../src/simu/00copyright.f90'

!*************************************************************************
! This program is part of the mM project. It isnsegmaxf
! devoted to the analysis of dislocation microstructures
! made with the mM simulation code. The default output file analyzed is
! the ASCI file "SEG_save".
! The different statistical analysis programmed are the folowing:
! 1) A global junction length analysis. Length of junctions are measured
!    and a file named "histojonc" of the junction s length is saved. This
!    data file can be used to construct histograms of junctions length.
! 2) A statistical measure of dislocations line roughness (or fractal
!    dimension). This measure is based on a calculation of dislocations
!    length included in a sphere of increasing radius R.
! 3) The dislocation microstructure (segment output file) is sorted
!    from small to long dislocations debris and saved in a file named LOOP_save.
! 4) The dislocation microstructure (segment output file) is analysed
!    to make statistics on jogs. Jogs are classified according to
!    a mobility criteria (balance between the resolved PK on each side
! 5) The dislocation microstructure (segment output file) is used
!    to make a map of internal stress at the nodes of a plane grid.
! 6) Analysis of the orientation of the junctions with respect to
!    the deplacement sens of the primary segments
! 7) Analysis of the distribution of the dislocation lines character
!    in the dislocation microstructure. Outpyt of this calculation is
!    a list of of segment character and length in file 'statcarac'.
! 8) The Nye tensor calculation and approximate GND density
! 9) The total dislocation density calculation in a thin foil
! 10) Resolved shear strain mapping from gammabox.bin
! 11) 3D -> 2D extraction of segments coordinates in a reference plane
!*************************************************************************
Program histo

!!!!!!!!!!!!!
implicit none
!!!!!!!!!!!!!

integer,parameter   :: DP=selected_real_kind(p=14)   !< Definition of real number kind
integer,parameter   :: DPI=selected_int_kind(13)     !< Definition of integer number kind
integer,parameter   :: DPI_S=selected_int_kind(9)    !< The integer format needed for the film.bin file

integer(kind=DPI),parameter   :: nsegmax=100000 !Maximum number of segments

integer(kind=DPI),parameter  :: nbasered = 8  !< Dimension de la base reduite
integer(kind=DPI),parameter  :: izero   = 0
integer(kind=DPI),parameter  :: iun     = 1

real(kind=DP), parameter  :: PI     = 3.14159265358979D0

type segment
    real(kind=dp)                  :: resdep   !< fraction non entiere du deplacement
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
    gd,                           &   !< The type of GD segments 0 = not GD, 1 = GD made by cross-slip, 2 = GD made by colli annihilation
    grain                             !< localisation du segment dans un polycristal

    logical ::  jonc,             &   !< clef d'immobilisation pour les jonctions
                diseg,            &   !< clef de discretization du segment
                bloquer,          &   !< clef de bloquage des segments sur une barriere
                dedant                !< clef pour les seg. entier dans la boite
end type segment

type(SEGMENT),allocatable :: seg(:),   &   !< The segments tab
                             segl(:),  &   !< The segment debris tab
                             segll(:)      !< The segment debris tab in order



integer(kind=dpi)   :: nb_step_stat !<
integer(kind=dpi)   :: I            !<
integer(kind=dpi)   :: J            !<
integer(kind=dpi)   :: indexB       !<
real(kind=dpi)      :: useless      !<


integer(kind=DPI) :: ijk          !<
integer(kind=DPI) :: k            !<
integer(kind=DPI) :: debutaffich  !<
integer(kind=DPI) :: freqaffich   !<
integer(kind=DPI) :: finaffich    !<

! Variables that must be of dim DPI_S
integer(kind=DPI_S)                   :: nbase          !< size of the mM segment list
integer(kind=DPI_S)                   :: nsegmaxf       !< maximum number of segments
integer(kind=DPI_S),allocatable       :: bveclin(:,:)   !< list of segments line direction
integer(kind=DPI_S),allocatable       :: bvecdep(:,:)   !< list of segment displacement direction
integer(kind=DPI_S),allocatable,save  :: bvecnor(:,:)   !< list of segment glide plane normal direction
integer(kind=DPI_S)                   :: icristallo     !< type of crystallo
integer(kind=DPI_S),Dimension(3)      :: modur          !< dimension of the simulation volume (a units)
integer(kind=DPI_S)                   :: kk             !< Simulation steps
integer(kind=DPI_S)                   :: nsegm          !< Number of segments
integer(kind=DPI_S)                   :: KODERR         !< ISTAT error code
integer(kind=DPI_S),allocatable       :: iseg(:,:)      !< The segments list
integer(kind=DPI_S),allocatable       :: TABVOIS(:,:)   !< The segment first neighbors
integer(kind=DPI_S),allocatable       :: junctt(:)      !< The junction pair-segment
integer(kind=DPI_S),allocatable       :: ID_surface(:)  !< the segments surface ID

integer(kind=DPI),Dimension(3) :: BoxCenter !<

real(kind=dp),allocatable    :: normlin(:) !<

real(kind=DP) :: VecBurgers   !<
real(kind=DP) :: dpoiss       !<
real(kind=DP) :: avalue       !<
real(kind=DP) :: epso         !<
real(kind=DP) :: deltat       !<
real(kind=DP) :: mu           !<
real(kind=DP) :: norm         !<
real(kind=DP) :: xstress      !<

integer(kind=dpi),dimension(3)  :: Lattice_O              !< Coordinate of the lattice origin for density analysis
integer(kind=dpi),dimension(3)  :: v1,v2,v3               !< 3 vectors to define the lattice elementary cell
integer(kind=dpi)               :: n1,n2,n3               !< 3 integers to define the lattice dimension in v1, v2, v3
integer(kind=dpi)               :: c1,c2,c3               !< 3 integers to define the ra vector coordinate
real(kind=dp), dimension(3)     :: A, B, C ,D             !< Vectors for the Cramer's rule
real(kind=dp)                   :: sysdet,xdet,ydet,zdet  !< Determinants in the Cramer's rule
real(kind=dp), dimension(3)     :: coord

logical                         :: Film         !< True if Film.bin is used
logical                         :: TheGoodStep  !< The Good step logical test

logical,allocatable             :: hide(:)  !<

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
character(len=8)    :: date       !<
character(len=10)   :: time       !<
character(len=15)   :: mydate     !<

real(kind=8),dimension(3,3) :: sigapp         !< The applied stress tensor
real(kind=8),dimension(3)   :: z              !< The applied stress axis
real(kind=8),dimension(3)   :: btmp,vectmp    !< Burgers Vector, tmp vector
real(kind=8),dimension(3)   :: PK             !< Peach-Koeler vector
real(kind=8)                :: TauBO,TauBE    !< Resolved shear stress at O and E

integer(kind=8),dimension(12,12) ::tabjonc

real(kind=8)                :: density,ampli
real(kind=8)                :: totjunclength,junclength
real(kind=8)                :: AvJuncLength
real(kind=8)                :: lengthO,lengthE
real(kind=8)                :: radius2,incradius,length(1:10000,1:3)
real(kind=8)                :: normvecO2,normvecE2
real(kind=DP),allocatable   :: normdep(:),normnor(:),Solli_sys(:),sysnormlin(:,:),sysnormnor(:,:),tens(:,:,:)
real(kind=8),allocatable    :: looplength(:),lengthl(:)
integer(kind=8)             :: LessSeg,LessSeg2,ratio,jsys
integer(kind=DPI_S)         :: ntsg,Nsegloop
integer(kind=8)             :: Oi(1:3),Ei(1:3),Oj(1:3),Ej(1:3),TOi(1:3),decal(1:3)
integer(kind=8)             :: O(1:3),E(1:3),vecO(1:3),vecE(1:3)
integer(kind=8)             :: IBV,PtJonc(1:4,4),CNN,buf1jonc,buf2jonc
integer(kind=8)             :: bufkjonc,bufojonc,bufejonc,bufbufjonc,indexvois,indexi,indexvoiso
integer(kind=8)             :: lengthinc,idO,idE,ibvO,ibvE,normeO,normeE,normeOinc,normeEinc
integer(kind=8)             :: question
logical                     :: planar
logical,allocatable         :: done(:)
real(kind=dp)               :: realtemp
logical                     :: startdebris

integer(kind=dpi)           :: numbdebris           !< Number of loop debris in case(3)
real(kind=dp),allocatable   :: listlength(:), ListLengthPrim(:), ListLengthNormal(:)    !< list of loop information in case(3)
real(kind=dp),allocatable   :: LengthPrim(:), LengthNormal(:)                           !< list of loop information in case(3)
real(kind=dp),dimension(3)  :: VecPrim, VecNormal   !< The two reference directions in case(3)
real(kind=dp)               :: VecPrimNorm, VecNormalNorm, MaxPrim, MinPrim, MaxNormal, MinNormal

integer(kind=DPI)           :: Nb_phase, Index_phase  !< number of phases used, index of the phase
                                                      !< Present even in non parallel mode on the material input file
integer(kind=DPI)           :: slipsys                !< Considered slip system

! More variables for case 1&10
real(kind=DP)                               :: Lccusp         !< The critical (maxi) length of cusp
real(kind=DP)                               :: Lcuspcumul     !< The cumulated length of cusp
integer(kind=DPI)                           :: nbjuncnonnul   !< Number of non zero junctions
integer(kind=DPI)                           :: totjuncsegs    !< Total number of junctions segments
integer(kind=DPI)                           :: nbjonc         !< Total number of junctions
integer(kind=DPI)                           :: nbjuncnul      !< Total number of zero junctions
integer(kind=DPI)                           :: juncsegs       !< Number of segments involve in a junction
integer(kind=DPI)                           :: idebjonc       !< Indice of the segment starting a junction
integer(kind=DPI)                           :: is             !< Slip system of the segment starting a junction
integer(kind=DPI)                           :: js             !< Junction binome slip system of the segment starting a junction
integer(kind=DPI)                           :: juncnum        !< Junction indice 1: Hirth, 2: Glissile, 3: Lomer

real(kind=DP), dimension(3)                 :: juncstat       !< Lenght per junction type 1: Hirth, 2: Glissile, 3: Lomer

integer(kind=DPI)                           :: icusp          !<
logical, dimension(:), allocatable          :: zip            !<
logical                                     :: vononj         !<
logical                                     :: cusp           !<
logical                                     :: venonj         !<

! variables for case 7
real(kind=DP)               :: RAYON2,ANGLE, normvli, signtest, TWOPII
real(kind=DP)               :: normalvec(1:3), normallinvec(1:3)
real(kind=DP),allocatable   :: anglevis(:)
real(kind=DP), parameter    :: HALF    = 0.5D0 ,              &
                               PII     = 3.14159265358979D0 , &
                               UN      = 1.
real(kind=DP),allocatable,save      :: VECNORLIN(:,:)    !< Base de vecteurs ligne normalises
integer(kind=DPI),allocatable,save  :: assoc(:,:)

! variables for case 8, 9 & 12
real(kind=DP),allocatable   :: nye_tensor(:,:,:,:,:)  !< Tab for the nye tensor on the lattice
real(kind=DP),allocatable   :: dens_GND(:,:,:)        !< Tab for the GND density on the lattice
integer(kind=DPI)           :: normi                  !< Norm of segment i
real(kind=DP), dimension(3) :: center                 !< Center of the considered elementary segment

real(kind=DP)               :: voltot, volb, tmp, tmp2, vltmp(3)
real(kind=DP)               :: ra(3), lensecmin, lengthcumul
real(kind=DP)               :: b1(3), b2(3),bj1(3), bj2(3)
real(kind=DP)               :: OiOjscaN,OiejscaN, vcourb(3)
integer(kind=DPI)           :: l, m, vli(3), sys, shortseg
integer(kind=DPI)           :: f(3), G(3), bi(3)
integer(kind=DPI)           :: p1(3),p3(3),buftlro(3,3),p2(3)
integer(kind=DPI)           :: ii, vdi(3),nvec(3), jj, sens(96)
integer(kind=DPI)           :: modurx,modury,modurz, lj, vlj(3), li
integer(kind=DPI)           :: ij, Oa(3)
integer(kind=DPI)           :: vnn, sensTlj(4), vnno,vnne, conf(4), compteur, sensTLjMemo
logical                     :: testpmpq
real(DP)                    :: alpha

integer(kind=DPI)   :: nbsegdep

real(kind=DP),allocatable :: dislocation_density(:,:,:,:)       !< Tab for the density on the lattice
real(kind=DP),allocatable :: junction_density(:,:,:,:)          !< For the junction density on the lattice

! Variables for case 12
integer(kind=DPI)          :: JK, ni, max_length,IX,IY,IZ
integer(kind=DPI_S)        :: NBoites, NBoitesX, NBoitesY, NBoitesZ
real(kind=DP)              :: tens_sum(3,3)
real(kind=DP)              :: Px,Py,Pz            !   <  Gammabox center location
real(kind=DP),allocatable  :: Gammabox_table(:,:),Gamma_res(:,:),    &
                              somme(:),tens_box(:,:,:)

! Variables for case 13
Real(kind=DP)                                 :: Pos        !< Plane position (d)
Real(kind=DP),dimension(3)                    :: Point      !< One point in the 3D space
Integer(kind=DPI),dimension(3)                :: Miller     !< Plane normal vector (a,b,c)
Real(kind=DP),dimension(3)                    :: Mil        !< Plane normal vector (a,b,c)
Integer(kind=DPI),dimension(:), allocatable   :: Sys2D      !< System number
Integer(kind=DPI),dimension(:), allocatable   :: Jonc2D     !< The junction binome
Integer(kind=DPI),dimension(:), allocatable   :: Sign2D     !< Sign of the line (Burgers vector direction in 2D)
Real(kind=DP),dimension(:,:), allocatable     :: rxy        !< One point along the tested segment
Real(kind=DP),dimension(3)                    :: V2Dx,V2Dy  !< The two principale directions in the extraction plane
Real(kind=DP),dimension(3)                    :: VecOi,VecEi!< Two vectors between point and O or E
Real(kind=DP)                                 :: ScaO,ScaE  !< Two scalar product
Real(kind=DP)                                 :: normivect  !<
Real(kind=DP),dimension(3)                    :: Oi2D,Ei2D  !<
Integer(kind=DPI),dimension(3)                :: btmp2D     !<

Integer(kind=DPI)           :: Veclin       !<

!!!!!!!!!!!!!!!!!!!
! Initialisations !
!!!!!!!!!!!!!!!!!!!
density=0.
totjunclength=0.
nbjonc=0
nbjuncnul=0
totjuncsegs=0

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The simulation parameters  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1,file='../in/input.dd',STATUS='OLD')

carac = " "
do while (carac /= "-")
  read (1,*) carac
enddo

read (1,*) materiau
! read (1,*) control
! read (1,*) segments

close (1)

! Identification of the crystal symetry
open(1,file="../in/"//materiau,STATUS='OLD')

! some useless data are by passed first
read(1,*) Nb_Phase
read(1,*) Index_Phase
read(1,*) mu        ! The shear modulus
mu = mu * 1.D9      ! now mu is in Pa
read(1,*) useless
read(1,*) dpoiss    ! Poisson ratio
read(1,*) useless
read(1,*) useless
read(1,*) useless

! the crystal symetry
read(1,*) crystal_structure(1:3)

read(1,*) useless
read(1,*) VecBurgers    ! The burgers vector in Angstrom
VecBurgers = VecBurgers * 1.d-10
close (1)

! The simulation elementary lattice parameter
 cristallo = '../out/BVD.'//crystal_structure(1:3)

open(1,FILE=cristallo,STATUS='OLD')

read (1,*) nbase , avalue

allocate (bveclin(3,nbase))
allocate (normlin(nbase))
allocate (bvecdep(3,nbase))
allocate (normdep(nbase))
allocate (bvecnor(3,nbase))
allocate (normnor(nbase))
allocate (sysnormlin(3,nbase/8))
allocate (sysnormnor(3,nbase/8))
allocate (tens(3,3,nbase/8))


jsys = 0
do j = 1,nbase

  read (1,*) i, bveclin(1:3,j),bvecdep(1:3,j),bvecnor(1:3,j)
  if(j /= i) print *, "stop, a problem appears when reading BVD"
  normlin(j)=sqrt(real(bveclin(1,j))**2+  &
                  real(bveclin(2,j))**2+  &
                  real(bveclin(3,j))**2)
  normdep(j)=sqrt(real(bvecdep(1,j))**2+  &
                  real(bvecdep(2,j))**2+  &
                  real(bvecdep(3,j))**2)
  normnor(j)=sqrt(real(bvecnor(1,j))**2+  &
                  real(bvecnor(2,j))**2+  &
                  real(bvecnor(3,j))**2)
  if (modulo(int(j-1),8)==0) then
    jsys = jsys +1
    sysnormlin(1:3,jsys) = bveclin(1:3,j)/normlin(j)
    sysnormnor(1:3,jsys) = bvecnor(1:3,j)/normlin(j)
  endif

enddo

tens(:,:,:) = 0.

do jsys=1,nbase/8
  do i=1,3
    do j=1,3
      tens(i,j,jsys) = (1./2.)*(sysnormlin(i,jsys)*sysnormnor(j,jsys)+sysnormlin(j,jsys)*sysnormnor(i,jsys))
    enddo
  enddo
enddo

close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loading of the segments microstructure !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1,file='../in/SEG_save',STATUS='OLD')

ntsg=nbase/8
allocate (solli_sys(ntsg))

read(1,*) Solli_Sys(1:ntsg)           ! Parameter of effective loading on the systems
read(1,*) Nsegm                       ! Nb of segments
read(1,*) modur(1),modur(2),modur(3)  ! Dimensions of the simulation cell

! The center of the simulation cell
BoxCenter(1:3) = modur(1:3)/2

allocate (seg(nsegm))
allocate (hide(nsegm)) ! used for case 5 and 13

do I = 1,nsegm

  read (1,*) j,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
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

close(1)

! The different calculation to proceed are now selected
question = 1

do while (question > 0 .and. question < 14)

  write(*,*) ""
  write(*,*) ""
  write(*,*) "What type of calculation do you want to proceed ?"
  write(*,*) "(the default segment file analysed is SEG_save!)"
  write(*,*) "1  =  Junctions length analysis from film.bin"
  write(*,*) "2  =  Dislocation lines roughness (fractal dimension)"
  write(*,*) "3  =  Dislocation loops distribution"
  write(*,*) "4  =  Collinear superjogs analysis"
  write(*,*) "5  =  Internal stress mapping"
  write(*,*) "6  =  Junctions symmetry analysis"
  write(*,*) "7  =  Segments character analysis"
  write(*,*) "8  =  Calculation of Nye tensor and approximate GND density"
  write(*,*) "9  =  Dislocation density analysis"
  write(*,*) "10 =  Internal resolved shear strain  mapping from gammabox.bin"
  write(*,*) "11 =  3D -> 2D extraction from film.bin of segments coordinates"
  write(*,*) "0  =  Exit"
  write(*,*) ""
  write(*,*) "Your decision ?"

  read(*,*) question
  write(*,*) ""
  write(*,*) ""

  select case (question)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (1)                                                !
  ! Computation of the true length of junctions, i.e. section where !
  ! the two segments making the junction are overlapping.           !
  ! Same as 1) but stats are made with the film.bin file            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(1)

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

    ! selected slipsys
    write (*,*) "Which slip system do you want to analyse (1-12), Number > 13 means all systems ? "
    read (*,*) slipsys
    do while (slipsys < 1)
      write (*,*) '!> Slip system number must be > 0 !'
      write (*,*) "Which slip system do you want to analyse (1-12), Number > 13 means all systems ? "
      read  (*,*) slipsys
    enddo

    ! We first select the part of the film.bin file we want to use for the stats
    write (*,*) "statistics start at step:  ?"
    read (*,*) debutaffich
    if (debutaffich  <= 0) debutaffich = 1

    write (*,*) "Snapshot frequency:  ?"
    read (*,*) freqaffich
    if (freqaffich  <= 0) freqaffich = 1

    ! Initializations
    nbsegdep  = 0
    kk        = 1

    ! Warning: The critical length for cusps can be changed
    Lccusp = 100 * normlin(1) * avalue
    print *, 'Maximum length of cups along junctions =',Lccusp,'normlin(1)=',normlin(1)

    ! Results are saved in files histojonc and statjonc
    open(1,file='../out/histojonc',STATUS='UNKNOWN')
    open(2,file='../out/statjonc',STATUS='UNKNOWN')

    write(2,*) '--> kk    nbjonc    totjunclength   avjonclength      &
              &  nbjoncnul       Hirth length       Glissile length       Lomer length'

    ! The file we want to analyze
    fichier = "film.bin"
    open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
!
    ! Important information at the begining of the film.bin file
    read(50) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
    read(50) nsegmaxf     ! The nsegmax value used in this simulation
    read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
    read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry
    read(50) (bveclin(:,I),bvecnor(:,I), I=1,nbase) ! table of vector direction and displacement for this simulation
    read(50) modur(:)     ! Dimension of the 3D simulated volume in a value

    allocate  (iseg(5,nsegmax))
    allocate  (TABVOIS(2,nsegmax))
    allocate  (junctt(nsegmax))
    allocate  (ID_surface(nsegmax))

    ! Test on the nsegmax definition
    if (nsegmax /= nsegmaxf) then
      write(*,*) 'The default value of nsegmax=100000 was modified in this film !'
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

    ! We start the loop on simulation steps saved in film.bin
    ijk = 0

    Steps_loop: do while (kk.gt.0)

      juncstat(:) = 0
      ijk = ijk + 1

      ! A new simulation step
      read(50,IOSTAT=KODERR) kk, deltat, NSEGM, XSTRESS, EPSO
      write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I7)') KK,EPSO,XSTRESS,nsegm

      if(KODERR < 0) then
        print*,"  End of the film.bin file"
        stop
      elseif(KODERR > 0) then
        print*,"Reading problem 2",KODERR,ijk
        kk=1
        cycle Steps_loop
      endif

      if ( (kk >= debutaffich) .and. (modulo(kk, freqaffich)== 0) ) then

        ! The segments micro structure is loaded
        READ(50,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),ID_surface(j), J=1,nsegm)

        if(KODERR < 0) then
          print*,"  End of the film.bin file"
          stop
        elseif(KODERR > 0) then
          print*,"Reading Error 2",KODERR,nsegm,j
          cycle Steps_loop
        endif

        ! We transfer the information from film.bin
        !  ======>>>>>
        print*,'kk=',kk

        deallocate (seg)
        allocate (seg(1:nsegm))

        ! The film.bin information is used to over write the SEG_save information
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

        ! print*,'out do loop',KK

        ! Start of the stats made on this step
        !  ======>>>>>
        allocate (zip(1:nsegm))
        zip(1:nsegm) = .false.

        ! The file header
        write(1,*) kk, EPSO
        write(1,*) '--> Index     Nb_segs     Sys     Junclength'

        seg_loop: do k = 1, nsegm

          ! Only segment with %jonc=true are tested
          if (.not. seg(k)%jonc) cycle seg_loop

          ! This segment has already been treated
          if (zip(k)) cycle seg_loop

          ! we select the segments to be analysed if needed
          if (slipsys < 13) then
            if (int((seg(seg(k)%ijonc)%veclin-1)/8)+1 /= slipsys .and.  &
                int((seg(k)%veclin-1)/8)+1 /= slipsys)  cycle seg_loop
          endif

          ! Initialisation
          i           = k         ! The segment we test
          vononj      = .true.
          Lcuspcumul  = 0.

          ! We look for the segment origin i of the junction
          do while (vononj)

            ! The vnno of i is tested
            if (seg(seg(i)%vnno)%jonc) then

              ! This is not the origine of the junction and we must look for the next vnno
              if (seg(seg(i)%vnno)%veclin == seg(i)%veclin) then
                i = seg(i)%vnno
              else
                ! This segment is part of a junction, but not in the same direction
                vononj = .false.
              endif

            else                 ! if vnno not junction

              icusp = seg(i)%vnno
              cusp  = .false.

              ! We look for junction continuity in length smaller than Lccusp ?
              do while (Lcuspcumul < Lccusp)

                Lcuspcumul = Lcuspcumul + seg(icusp)%norme * normlin(seg(icusp)%veclin) * avalue

                !If we find a seg in junction configuration, we continue the junction
                if (seg(icusp)%vnno > nsegm) exit

                if ( seg(seg(icusp)%vnno)%jonc .and. seg(seg(icusp)%vnno)%veclin == seg(k)%veclin ) then
                  cusp = .true.
                  i    = seg(icusp)%vnno
                  exit
                endif

                icusp = seg(icusp)%vnno

              enddo

              vononj = .false.

              ! A cusp was found along the junction
              if (cusp) vononj=.true.

            endif

          end do

          ! The segment where the real junction is starting
          idebjonc  = i

          ! A paranoid test
          if (seg(i)%ijonc == 0) then
            print *,'problem with i'
            STOP
          endif

          nbjonc      = nbjonc + 1
          juncsegs    = 0
          junclength  = 0.
          venonj      = .true.
          Lcuspcumul  = 0

          ! the loop to find the end of the junction
          do while (venonj)

            if (seg(i)%jonc) then

              if (seg(i)%veclin == seg(k)%veclin) then

                j         = seg(i)%ijonc
                juncsegs  = juncsegs + 1
                zip(i)    = .true.

                ! Segments are shifted at the center of the simulation cell to avoid problems
                ! with the periodic boundary conditions.
                IBV=seg(i)%veclin
                TOi(1:3)=(seg(i)%O(1:3)+seg(i)%norme/2*Bveclin(1:3,IBV)-modur(1:3)/2)
                Oi(1:3) = modulo(int(seg(i)%o(1:3)-toi(1:3),DPI),modur(1:3))
                Oj(1:3) = modulo(int(seg(J)%o(1:3)-toi(1:3),DPI),modur(1:3))
                Ei(1:3) = modulo(int(seg(i)%o(1:3)+seg(i)%norme*bveclin(1:3,IBV)-toi(1:3),DPI),modur(1:3))
                Ej(1:3) = modulo(int(seg(J)%o(1:3)+seg(J)%norme*bveclin(1:3,seg(J)%veclin)-toi(1:3),DPI),modur(1:3))

                ! We check that the two segments are realy along the same line
                ! cnn is the index of one bveclin dimension indices different from zero
                decal(1:3) = Oj(1:3) - Oi(1:3)

                if (bveclin(1,IBV) /= 0) then
                  cnn=1
                  if (decal(1) /= 0 .and. decal(2) /= 0 .and.          &
                  (decal(1)*bveclin(2,seg(i)%veclin))/(decal(2)*bveclin(1,seg(i)%veclin)) /= 1) then
                    print*,'Segments',i,'and',j,'are disconnected 1-1'
                    print*,seg(i)%o,seg(i)%norme,seg(i)%veclin
                    print*,seg(j)%o,seg(j)%norme,seg(j)%veclin
                    print*,oi,oj
                    print*,ei,ej
                    stop
                  elseif (decal(1) /= 0 .and. decal(3) /= 0 .and.          &
                  (decal(1)*bveclin(3,seg(i)%veclin))/(decal(3)*bveclin(1,seg(i)%veclin)) /= 1) then
                    print*,'Segments',i,'and',j,'are disconnected 1-2'
                    print*,seg(i)%o,seg(i)%norme,seg(i)%veclin
                    print*,seg(j)%o,seg(j)%norme,seg(j)%veclin
                    print*,oi,oj
                    print*,ei,ej
                    stop
                  endif
                elseif (bveclin(2,IBV) /= 0) then
                  cnn=2
                  if (decal(1) /= 0 .and. decal(2) /= 0 .and.          &
                  (decal(1)*bveclin(2,seg(i)%veclin))/(decal(2)*bveclin(1,seg(i)%veclin)) /= 1) then
                    print*,'Segments',i,'and',j,'are disconnected 2-1'
                    print*,seg(i)%o,seg(i)%norme,seg(i)%veclin
                    print*,seg(j)%o,seg(j)%norme,seg(j)%veclin
                    print*,oi,oj
                    print*,ei,ej
                    stop
                  elseif (decal(1) /= 0 .and. decal(3) /= 0 .and.          &
                  (decal(1)*bveclin(3,seg(i)%veclin))/(decal(3)*bveclin(1,seg(i)%veclin)) /= 1) then
                    print*,'Segments',i,'and',j,'are disconnected 2-2'
                    print*,seg(i)%o,seg(i)%norme,seg(i)%veclin
                    print*,seg(j)%o,seg(j)%norme,seg(j)%veclin
                    print*,oi,oj
                    print*,ei,ej
                    stop
                  endif
                elseif (bveclin(3,IBV) /= 0) then
                  cnn=3
                  if (decal(1) /= 0 .and. decal(2) /= 0 .and.          &
                  (decal(1)*bveclin(2,seg(i)%veclin))/(decal(2)*bveclin(1,seg(i)%veclin)) /= 1) then
                    print*,'Segments',i,'and',j,'are disconnected 3-1'
                    print*,seg(i)%o,seg(i)%norme,seg(i)%veclin
                    print*,seg(j)%o,seg(j)%norme,seg(j)%veclin
                    print*,oi,oj
                    print*,ei,ej
                    stop
                  elseif (decal(1) /= 0 .and. decal(3) /= 0 .and.          &
                  (decal(1)*bveclin(3,seg(i)%veclin))/(decal(3)*bveclin(1,seg(i)%veclin)) /= 1) then
                    print*,'Segments',i,'and',j,'are disconnected 3-2'
                    print*,seg(i)%o,seg(i)%norme,seg(i)%veclin
                    print*,seg(j)%o,seg(j)%norme,seg(j)%veclin
                    print*,oi,oj
                    print*,ei,ej
                    stop
                  endif
                else
                  print *, "I end up here in histo, maybe a case was not taken into account - what about cnn value ?"
                  stop
                endif

                ! Here the algorithm used in the simulation (11topolo.f90) is used as it is to
                ! calcul the jonction length. In the folowing:
                ! ptjonc(n,m)
                !
                !   n  ,  m -->
                !         1           2           3            4
                !   1    Oix         Eix         Ojx          Ejx
                !   2    Oiy         Eiy         Ojy          Ejy
                !   3    Oiz         Eiz         Ojz          Ejz
                !   4     0        Eix-Oix      Ojx-Oix     Ejx-Oix  (where x is the ccn dimension)
                !
                ! The tab ptjonc defines the 4 points appearing in the junction construction
                ! respectively Oi Ei Oj Ej
                ! buf1junc : the point with the smallest (negative) coordinate
                ! bufojunc : the point at the begining of the junction
                ! bufejunc : the point at the end of the junction
                ! buf2junc : the point with the largest (positive) coordinate

                PtJonc(1:3,1) = Oi(1:3)
                PtJonc(1:3,2) = Ei(1:3)
                PtJonc(1:3,3) = Oj(1:3)
                PtJonc(1:3,4) = Ej(1:3)
                PtJonc(4,1)   = 0
                PtJonc(4,2)   = PtJonc(CNN,2)-PtJonc(CNN,1)
                PtJonc(4,3)   = PtJonc(CNN,3)-PtJonc(CNN,1)
                PtJonc(4,4)   = PtJonc(CNN,4)-PtJonc(CNN,1)

                ! Identification of the two extremities
                buf1jonc = 1
                buf2jonc = 1
                do bufkjonc = 2, 4, 1
                  if (PtJonc(4,bufkjonc) <= PtJonc(4,buf1jonc)) buf1jonc = bufkjonc
                  if (PtJonc(4,bufkjonc) >  PtJonc(4,buf2jonc)) buf2jonc = bufkjonc
                enddo

                ! Identification of the points starting and ending the junction
                bufojonc = 0
                bufejonc = 0
                do bufkjonc = 1, 4, 1
                  if (bufkjonc /= buf1jonc .and. bufkjonc /= buf2jonc .and. bufojonc == 0)  &
                    bufojonc = bufkjonc
                  if (bufkjonc /= buf1jonc .and. bufkjonc /= buf2jonc .and. bufkjonc /= bufojonc .and. bufejonc == 0) &
                    bufejonc = bufkjonc
                enddo

                if (PtJonc(4,bufojonc) > PtJonc(4,bufejonc)) then
                  bufbufjonc = bufojonc
                  bufojonc   = bufejonc
                  bufejonc   = bufbufjonc
                endif

                if (juncsegs == 1) then
                  ! The particular case of two binome junction segments with no overlap
                  ! This particular case may exist temporarily since we test junction
                  ! zipping-unzipping during a few steps
                  if ((buf1jonc == 1 .and. bufojonc == 2) .or.    &
                      (buf1jonc == 2 .and. bufojonc == 1) .or.    &
                      (buf1jonc == 3 .and. bufojonc == 4) .or.    &
                      (buf1jonc == 4 .and. bufojonc == 3)) then
                    ratio = 0  ! The junction length is set to zero
                  else
                    if (seg(i)%vnne < nsegm) then
                      if (seg(seg(i)%vnne)%jonc) then
                        if (bveclin(CNN,IBV) > 0) then
                          ! The (Ei - bufojonc) length
                          ratio =abs((ptjonc(CNN,2) - (ptjonc(CNN,1)+ptjonc(4,bufojonc))) / bveclin(CNN,IBV))
                          print*, "ration 1 ", ratio
                        else
                          ! The (Ei - bufojonc) length
                          ratio =abs((ptjonc(CNN,2) - (ptjonc(CNN,1)+ptjonc(4,bufejonc))) / bveclin(CNN,IBV))
                          print*, "ration 2 ", ratio
                        endif
                      else
                        ! The overlaping length
                        ratio =abs((ptjonc(4,bufejonc) - ptjonc(4,bufojonc)) / bveclin(CNN,IBV))
                        print*, "ration 3 ", ratio
                        print*, buf1jonc, bufojonc, bufejonc, buf2jonc
                      endif
                    endif
                  endif
                else
                  if (seg(i)%vnne < nsegm) then
                    if ( seg(seg(i)%vnne)%jonc ) then
                      ! The total i segment length
                      ratio = seg(i)%norme
                    else
                      if (bveclin(CNN,IBV) > 0) then
                        ! The (Bufejonc - Oi) length
                        ratio = abs((ptjonc(CNN,1) + ptjonc(4,bufejonc) - ptjonc(CNN,1)) / bveclin(CNN,IBV))
                      else
                        ! The (Bufejonc - Oi) length
                        ratio = abs((ptjonc(CNN,1) + ptjonc(4,bufojonc) - ptjonc(CNN,1)) / bveclin(CNN,IBV))
                      endif
                    endif
                  endif
                endif

                ! The junction length
                junclength = junclength +  ratio * normlin(IBV) * avalue

                !go to the next e neighbor
                i=seg(i)%vnne
                if (i > nsegm)  venonj = .false.

              else

                venonj = .false.

              endif

            else

              icusp = i
              cusp  = .false.
              ! Do we have cusp along the Lccusp length ?
              do while ( Lcuspcumul < Lccusp)
                Lcuspcumul = Lcuspcumul + seg(icusp)%norme*normlin(seg(icusp)%veclin)*avalue
                !If we have a seg in junction configuration, we continue the junction
                if (seg(icusp)%vnne > nsegm) exit
                if ( seg(seg(icusp)%vnne)%jonc .and. seg(seg(icusp)%vnne)%veclin == seg(k)%veclin ) then
                  cusp=.true.
                  i=seg(icusp)%vnne
                  exit
                endif
                icusp=seg(icusp)%vnne
              enddo
              venonj=.false.
              if (cusp) venonj=.true.

            endif

          end do ! Endloop for the segs owned by the junction

          ! Junction type selection
          is = int((seg(idebjonc)%veclin-1)/8)+1
          js = int((seg(seg(idebjonc)%ijonc)%veclin-1)/8)+1
          juncnum = tabjonc(is,js)

          ! Paranoid tests
          ! check if colli is zero length
          if (juncnum > 3 .and. junclength /= 0.) then
            print*, " A collinear reaction with non zero overlap length is found"
            print*, " Overlap length = ",junclength
            print*, " juncsegs = ", juncsegs
            print*, "Segment i = ", idebjonc, is, seg(idebjonc)%norme
            print*, bveclin(1:3, seg(idebjonc)%veclin)
            print*, "Segment j = ", seg(idebjonc)%ijonc, js, seg(seg(idebjonc)%ijonc)%norme
            print*, bveclin(1:3, seg(seg(idebjonc)%ijonc)%veclin)
            print*, Oi
            print*, Ei
            print*, Ej
            print*, Oj
            Stop
          endif
           ! Non defined junctions are found
          if (juncnum == 0) then
            print *, idebjonc, seg(idebjonc)%ijonc, is, js, juncnum, junclength, juncsegs, &
                    seg(idebjonc)%norme,seg(seg(idebjonc)%ijonc)%norme
            stop "Problem junction type"
          endif

          ! The solution we found
          if (juncnum < 4 .and. juncnum > 0 ) then

            juncstat(juncnum) = juncstat(juncnum) + junclength

            ! Progress of the calculation are saved
            totjunclength = totjunclength + junclength
            totjuncsegs   = totjuncsegs   + juncsegs
            if ( junclength == 0) nbjuncnul = nbjuncnul + 1

            ! We save the information regarding this junction
            write(1,63) idebjonc,juncsegs,int( ((seg(idebjonc)%veclin-1)/8)+1),                   &
                        int(((seg(seg(idebjonc)%ijonc)%veclin-1)/8)+1),                           &
                        seg(idebjonc)%veclin,seg(seg(idebjonc)%ijonc)%veclin,junclength,juncnum
            63 format(1X,6I6,E13.5,I6)
          endif

        enddo seg_loop

        deallocate (zip)

        nbjuncnonnul = nbjonc - nbjuncnul
        ! The Average length of a junction
        if (nbjuncnonnul /= 0) then
          AvJuncLength = totjunclength/nbjuncnonnul
        else
          AvJuncLength = 0.
        endif

        write(1,*) 'Number of junctions =',nbjonc
        write(1,*) 'Number of junctions with zero length =',nbjuncnul

        write(1,*) 'Average length of a junction (without zero length junctions in the calculation =', AvJuncLength
        write(1,*) 'Total segs treated =',totjuncsegs
        write(1,*) 'Hirth total length =',juncstat(1)
        write(1,*) 'Glissile total length =',juncstat(2)
        write(1,*) 'Lomer total length =',juncstat(3)

        ! End of the stats made on this step
        !  ======>>>>>
        write(2,44) kk,nbjonc,totjunclength,AvJuncLength,nbjuncnul, juncstat(:)
        44 format(I10,I6,2E13.5,I6,3E13.5)

        ! On remet les variables a 0
        nbjonc        = 0
        nbjuncnul     = 0
        nbjuncnonnul  = 0
        totjuncsegs   = 0
        totjunclength = 0

        ! The mark for the end of step is read
        read(50,IOSTAT=KODERR) caractest

        if(KODERR < 0) then
          print*,"  End of the film.bin file"
          stop
        elseif(KODERR > 0) then
          print*,"Reading error 3",caractest
          cycle Steps_loop
        endif

      else

        READ(50,IOSTAT=KODERR)

        caractest = '@'
        do while (caractest /= '&')
          read(50,IOSTAT=KODERR) caractest

          if(KODERR < 0) then
            print*,"  End of the film.bin file"
            stop
          elseif(KODERR > 0) then
            print*,"Reading problem 666",KODERR
            stop
          endif
        enddo

      endif

    enddo Steps_loop

    print *,'results are in /out/histojonc & /out/statjonc'
    read *

    deallocate  (iseg)
    deallocate  (TABVOIS)
    deallocate  (junctt)
    deallocate  (ID_surface)

    close(50)
    close(1)
    close(2)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (2)                                                 !
  ! Computation of the dislocations line roughness (or fractal      !
  ! dimension). This measure is based on the calculation of         !
  ! dislocation length included in a sphere of increasing radius R. !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(2)

    length(:,:) = 0.

    ! If we must restrict the calculation to the glide plane, planar = true
    planar = .true.
    !planar = .false.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Length is measured form the origin of each segment in the microstructure
    loopseg: do i=1,nsegm

      if (seg(i)%norme == 0) cycle    ! Only non-zero segments are tested

      print*,"i=",i

      ! the point where the measure of length start
      Oi(1:3) = seg(i)%o(1:3)
      !print*,"segment Oi= ",i,Oi
      O(1:3) = Oi(1:3)         ! The length mark along the O side
      E(1:3) = Oi(1:3)         ! The length mark along the E side
      ! initializations
      lengthO = 0.    ! Cumulated length on the O side
      lengthE = 0.    ! Cumulated length on the E side
      lengthinc = 0   ! Number of raduis increment tested

      ! 1st seg to explore on the "O" side
      idO = seg(i)%vnno            ! Identity of the segment considered to measure the length
      if (idO == 0) cycle loopseg  ! A end point, we move to the next Oi
      ibvO = seg(idO)%veclin       ! the character od idO
      normeO = seg(idO)%norme      ! the norme of segment idO
      normeOinc = 0                ! Increment of norme tested

      ! 2d seg to explore on the "E" side
      idE = i
      ibvE = seg(idE)%veclin
      normeE = seg(idE)%norme
      normeEinc = 0

      !  print*,ibvO,normeO,ibvE,normeE

      if (planar) then
        !The two segments are not in the same glide plane
        if((ibvO-1)/8 .ne. (ibvE-1)/8) cycle loopseg  ! Segments are not in the same glide plane
      endif

      ! A test to check that every thing is OK
      if (idO == idE) then
        write(*,*) "A micro-loop of non-zero length is found in the microstructure!"
        stop
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! the loop on the increasing sphere of radius = sqrt(radius2)
      radius2 = 0
      incradius = 10                     ! The increment size for increasing radius
      surradius : do while (radius2 >= 0)
        radius2 = (sqrt(radius2) + incradius)**2     ! The radius tested in the loop
        lengthinc = lengthinc + 1                    ! One more larger radius is tested
        !print*,"radius2 =", lengthinc, radius2

          ! The "O" side
          normeOinc = normeOinc + 1        ! Since a new radius is tested an increment along the line is added
          if (normeOinc > normeO) then
            ! The total length of the segment is reached, we move to the next one on the O side
            if (planar) then
              ! the following line is needed to consider only the roughness in the initial glide plane
              if((seg(idO)%veclin-1)/8 .ne. (seg(seg(idO)%vnno)%veclin-1)/8) cycle loopseg  ! Segments are not in the same glide plane
            endif
            idO = seg(idO)%vnno          ! Identity of the segment considered to measure the length
            if (idO == 0) cycle loopseg  ! A end point, we move to the next Oi
            if (idO == idE) then
              ! The all procedure is stoped since the dislocation loop is closed
              ! we move to the next Oi
              exit surradius
            endif
            ibvO = seg(idO)%veclin   ! the character od idO
            normeO = seg(idO)%norme  ! the norme of segment idO
            normeOinc = 1            ! Increment of norme tested
          endif
          O(:) = O(:) + bveclin(:,ibvO)                       ! The mark on the O side is moved
          vecO(:) = (O(:) - Oi(:))                            ! A new vector Oi-O is defined
          normvecO2 = vecO(1)**2 + vecO(2)**2 + vecO(3)**2    ! The norme of Oi-O is calculated

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! The loop on segments inside the radius : The O side
          surnormvecO : do while (normvecO2 < radius2)

            lengthO = lengthO + normlin(ibvO)*avalue    ! The length of line on the O side is increased

            !print*,"lengthO",(lengthO/avalue)**2,normvecO2,radius2
            !read(*,*)

            normeOinc = normeOinc + 1        ! A new increment of length is tested
            if (normeOinc > normeO) then
              ! The total length of the segment is reached, we move to the next one on the O side
              if (planar) then
                ! the following line is needed to consider only the roughness in the initial glide plane
                if((seg(idO)%veclin-1)/8 .ne. (seg(seg(idO)%vnno)%veclin-1)/8) cycle loopseg  ! Segments are not in the same glide plane
              endif
              idO = seg(idO)%vnno          ! Identity of the segment considered to measure the length
              if (idO == 0) cycle loopseg  ! A end point, we move to the next Oi
              if (idO == idE) then
                ! The all procedure is stoped since the dislocation loop is closed
                ! we move to the next Oi
                exit surradius
              endif
              ibvO = seg(idO)%veclin   ! the character od idO
              normeO = seg(idO)%norme  ! the norme of segment idO
              normeOinc = 0            ! Increment of norme tested
            endif
            O(:) = O(:) + bveclin(:,ibvO)                       ! The mark on the O side is moved
            vecO(:) = (O(:) - Oi(:))                            ! A new vector Oi-O is defined
            normvecO2 = vecO(1)**2 + vecO(2)**2 + vecO(3)**2    ! The norme of Oi-O is calculated

          enddo surnormvecO   ! the normvecO2 loop
          O(:) = O(:) - bveclin(:,ibvO)   ! The last move in the do loop was not OK
          ! The end of line inside the radius is reached on the O side

          ! The "E" side
          normeEinc = normeEinc + 1        ! Since a new radius is tested an increment along the line is added
          if (normeEinc > normeE) then
            ! The total length of the segment is reached, we move to the next one on the E side
            if (planar) then
              ! the following line is needed to consider only the roughness in the initial glide plane
              if((seg(idE)%veclin-1)/8 .ne. (seg(seg(idE)%vnne)%veclin-1)/8) cycle loopseg  ! Segments are not in the same glide plane
            endif
            idE = seg(idE)%vnne          ! Identity of the segment considered to measure the length
            if (idE == 0) cycle loopseg  ! A end point, we move to the next Oi
            if (idO == idE) then
              ! The all procedure is stoped since the dislocation loop is closed
              ! we move to the next Oi
              exit surradius
            endif
            ibvE = seg(idE)%veclin   ! the character od idO
            normeE = seg(idE)%norme  ! the norme of segment idO
            normeEinc = 1            ! Increment of norme tested
          endif
          E(:) = E(:) + bveclin(:,ibvE)                       ! The mark on the E side is moved
          vecE(:) = (E(:) - Oi(:))                            ! A new vector Oi-E is defined
          normvecE2 = vecE(1)**2 + vecE(2)**2 + vecE(3)**2    ! The norme of Oi-E is calculated

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! The loop on segments inside the radius : The E side
          surnormvecE : do while (normvecE2 < radius2)

            lengthE = lengthE + normlin(ibvE)*avalue    ! The length of line on the E side is increased

            normeEinc = normeEinc + 1        ! A new increment of length is tested
            if (normeEinc > normeE) then
              ! The total length of the segment is reached, we move to the next one on the E side
              if (planar) then
                ! the following line is needed to consider only the roughness in the initial glide plane
                if((seg(idE)%veclin-1)/8 .ne. (seg(seg(idE)%vnne)%veclin-1)/8) cycle loopseg  ! Segments are not in the same glide plane
              endif
              idE = seg(idE)%vnne          ! Identity of the segment considered to measure the length
              if (idE == 0) cycle loopseg  ! A end point, we move to the next Oi
              if (idO == idE) then
                ! The all procedure is stoped since the dislocation loop is closed
                ! we move to the next Oi
                exit surradius
              endif
              ibvE = seg(idE)%veclin   ! the character od idO
              normeE = seg(idE)%norme  ! the norme of segment idO
              normeEinc = 0            ! Increment of norme tested
            endif
            E(:) = E(:) + bveclin(:,ibvE)                       ! The mark on the E side is moved
            vecE(:) = (E(:) - Oi(:))                            ! A new vector Oi-E is defined
            normvecE2 = vecE(1)**2 + vecE(2)**2 + vecE(3)**2    ! The norme of Oi-E is calculated

          enddo surnormvecE   ! the normvecE2 loop
          E(:) = E(:) - bveclin(:,ibvE)   ! The last move in the do loop was not OK
          ! The end of line inside the radius is reached on the E side

          ! Before moving to the next radius, we save results at this step
          length(lengthinc,1) = sqrt(radius2)*avalue                    ! The sphere radius
          length(lengthinc,2) = length(lengthinc,2) + lengthO + lengthE ! Length of dislocation starting from Oi in this radius
          length(lengthinc,3) = length(lengthinc,3) + 1                 ! The number of measures done at this radius

      enddo surradius   ! radius loop

    enddo loopseg ! i seg loop

    print*,"The computation proceeded whithout problem!"

    ! Results are saved in the file roughness
    open(1,file='../out/roughness',STATUS='UNKNOWN')
    ! The file header
    write(1,*) '-->  Index     radius     length     stat'
    i=1
    do while (i /= 0)
      write(1,34) i,length(i,1),length(i,2)/length(i,3),length(i,3)
      34 format(1X,I6,3E13.5)
      if (length(i+1,1) /= 0) then
        i = i + 1
      else
        exit
      endif
    enddo

    close(1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (3)                                                 !
  ! The dislocation microstructure (segment output file) is sorted  !
  ! from small to long dislocations debris and saved in a file      !
  ! named LOOP_save.                                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(3)

    allocate (done(nsegm))           ! To remember if a segment was tested
    allocate (looplength(nsegm))     ! To cmulate the debris length
    allocate (lengthl(nsegm))
    allocate (segl(nsegm))           ! A seg tab to save the loop debris result
    allocate (segll(nsegm))          ! A seg tab to save the loop debris sorting result from smaller to bigger
    allocate (listlength(nsegm))        ! The list of loop debris length
    allocate (listlengthPrim(nsegm))    ! The list of loop debris length in the Prim direction
    allocate (listlengthNormal(nsegm))  ! The list of loop debris length in the Normal direction
    allocate (lengthPrim(nsegm))        ! The list of loop debris length in the Prim direction
    allocate (lengthNormal(nsegm))      ! The list of loop debris length in the Normal direction

    ! Initialization
    done(:) = .false.
    looplength(:) = 0
    Nsegloop = 0
    numbdebris = 0

    write(*,*) "What is the loop length primary direction ?"
    read(*,*) VecPrim(1:3)
    VecPrimNorm = sqrt(VecPrim(1)**2+VecPrim(2)**2+VecPrim(3)**2)
    VecPrim(1:3) = VecPrim(1:3)/VecPrimNorm
    write(*,*) "What is the loop length normal direction ?"
    read(*,*) VecNormal(1:3)
    VecNormalNorm = sqrt(VecNormal(1)**2+VecNormal(2)**2+VecNormal(3)**2)
    VecNormal(1:3) = VecNormal(1:3)/VecNormalNorm

    ! The segment loop
    loopi: do i=1,nsegm

      ! This segment is now tested, no reason to come back later
      done(i) = .true.

      ! Dislocation line including junction are not considered as loop debris
      if (seg(i)%jonc) cycle loopi

      ! Calculation of the loop length
      indexvois = seg(i)%voise
      looplength(i) = seg(i)%norme*normlin(seg(i)%veclin)*avalue     ! The first increment of length is saved
      ! Calculation of the loop dimensions in two reference directions
      O(1:3) = seg(i)%O(1:3)
      E(1:3) = O(1:3) + seg(i)%norme*bveclin(1:3,seg(i)%veclin)
      MaxPrim = max(dot_product(O(:), VecPrim(1:3)), dot_product(E(:), VecPrim(1:3)))
      MinPrim = min(dot_product(O(:), VecPrim(1:3)), dot_product(E(:), VecPrim(1:3)))
      MaxNormal = max(dot_product(O(:), VecNormal(1:3)), dot_product(E(:), VecNormal(1:3)))
      MinNormal = min(dot_product(O(:), VecNormal(1:3)), dot_product(E(:), VecNormal(1:3)))

      j = 2
      ! The hand checking
      do while (j > 0)

        indexi = indexvois

        ! This dislocation has already been tested
        if (done(indexi)) cycle loopi
        done(indexi) = .true.

        ! Dislocation line including junction are not considered as loop debris
        if (seg(indexi)%jonc) cycle loopi
        ! Dislocation ending with a pinning point are not considered as loop debris
        if (seg(indexi)%voise == 0) cycle loopi

        ! Calculation of the loop length
        indexvois = seg(indexi)%voise
        looplength(i) = looplength(i) + seg(indexi)%norme*normlin(seg(indexi)%veclin)*avalue
        ! Calculation of the loop dimensions in two reference directions
        O(1:3) = E(1:3)
        E(1:3) = O(1:3) + seg(indexi)%norme*bveclin(1:3,seg(indexi)%veclin)
        MaxPrim = max(dot_product(O(:), VecPrim(1:3)), dot_product(E(:), VecPrim(1:3)), MaxPrim)
        MinPrim = min(dot_product(O(:), VecPrim(1:3)), dot_product(E(:), VecPrim(1:3)), MinPrim)
        MaxNormal = max(dot_product(O(:), VecNormal(1:3)), dot_product(E(:), VecNormal(1:3)), MaxNormal)
        MinNormal = min(dot_product(O(:), VecNormal(1:3)), dot_product(E(:), VecNormal(1:3)), MinNormal)
        LengthPrim(i) = (MaxPrim - MinPrim) * avalue
        LengthNormal(i) = (MaxNormal - MinNormal) * avalue

        !the dislocation loop is closed (indexvois = i)
        if (indexvois == i) exit

        ! One more segment must be tested
        j = j + 1

      enddo

      ! We must do one stupid loop to determine the number of pivot segments
      indexi = i
      LessSeg = 0
      do k = nsegloop + 1, nsegloop + j
        if (seg(indexi)%norme == 0 .and. seg(indexi)%gd < iun) LessSeg = LessSeg + 1
        indexi = seg(indexi)%voise
      enddo

      ! Segments part of debris are saved in a new tab named segl
      indexi = i
      LessSeg2 = 0

      ! The debris loop informations are save in a tab Listxxx
      numbdebris = numbdebris + 1
      listlength(numbdebris) = looplength(i)
      ListLengthPrim(numbdebris) = LengthPrim(i)
      ListLengthNormal(numbdebris) = LengthNormal(i)

      do k = nsegloop+1, nsegloop + j

        if (seg(indexi)%norme == 0 .and. seg(indexi)%gd < iun) then
          LessSeg2 = LessSeg2 + 1
        else
          segl(k-LessSeg2) = seg(indexi)
          lengthl(k-LessSeg2) = looplength(i)
          ! Attention we must intialize the neighbourg info
          if (k - LessSeg2 == nsegloop+1) then
              segl(k- LessSeg2)%voiso = -(j - 1 - LessSeg)
              segl(k- LessSeg2)%vnno  = -(j - 1 - LessSeg)
              segl(k- LessSeg2)%voise = -1
              segl(k- LessSeg2)%vnne  = -1
          elseif (k -LessSeg2 == nsegloop +j - LessSeg) then
              segl(k- LessSeg2)%voiso = -1
              segl(k- LessSeg2)%vnno  = -1
              segl(k- LessSeg2)%voise = -(j - 1 - LessSeg)
              segl(k- LessSeg2)%vnne  = -(j - 1 - LessSeg)
          else
              segl(k- LessSeg2)%voiso = -1
              segl(k- LessSeg2)%vnno  = -1
              segl(k- LessSeg2)%voise = -1
              segl(k- LessSeg2)%vnne  = -1
          endif
        endif

        indexi = seg(indexi)%voise

      enddo

      ! The total number of segments saved is segl is saved
      nsegloop = nsegloop + j - LessSeg

    enddo loopi ! i seg loop

    ! The listlength array is sorted in ascending order of magnitude
L1: do i = 1, numbdebris-1
L2:   do j = i+1, numbdebris
        if (listlength(i) > listlength(j)) then
          realtemp = listlength(j)
          listlength(j) = listlength(i)
          listlength(i) = realtemp

          realtemp = listlengthPrim(j)
          listlengthPrim(j) = listlengthPrim(i)
          listlengthPrim(i) = realtemp

          realtemp = listlengthNormal(j)
          listlengthNormal(j) = listlengthNormal(i)
          listlengthNormal(i) = realtemp
        endif
      enddo L2
    enddo L1

    ! The segl list is reconstructed in ascending order
    done(:) = .false.
    k = 0
L3: do i = 1, numbdebris
      startdebris = .false.
L4:   do j = 1,nsegloop
        if (.not. done(j) .and. lengthl(j) == listlength(i)) then
          k = k + 1
          startdebris = .true.
          done(j) = .true.
          segll(k) = segl(j)
        endif
        if (startdebris .and. lengthl(j) /= listlength(i)) exit L4
      enddo L4
    enddo L3

    ! Paranoid test
    if (k /= nsegloop) then
      print*,'Problem in the loop list reconstruction', k, nsegloop
      stop
    endif

    ! Results are saved in the file roughness
    open(1,file='../in/LOOP_save',STATUS='UNKNOWN')

    ! The file header
    write(1,'(12F4.0)') Solli_Sys(1:ntsg)           ! Parameter of effective loading on the systems
    write(1,* ) Nsegloop                    ! Nb of segments
    write(1,*) modur(1),modur(2),modur(3)  ! Dimensions of the simulation cell

    do i = 1,Nsegloop
      if (segll(i)%Ijonc == nsegmax) segll(i)%Ijonc = 0
      write(1,61) i,segll(i)%O(1:3),segll(i)%norme,segll(i)%veclin,segll(i)%voiso,segll(i)%vnno, &
                segll(i)%voise,segll(i)%vnne,segll(i)%JONC,segll(i)%Ijonc,segll(i)%gd,lengthl(i)
    enddo
    61 format(I5,3I10,I5,I7,4I7,L3,I7,I3,3X,E13.6)

    close(1)

    ! Results are saved in the file roughness
    open(1,file='../in/LOOP_list',STATUS='UNKNOWN')
    do i = 1, numbdebris
      write(1,62) i,listlength(i),ListLengthPrim(i),ListLengthNormal(i)
    enddo
    62 format(I5,3E13.6)
    close(1)

    print*,"The computation proceeded whithout problem!"

    deallocate (done)           ! To remember if a segment was tested
    deallocate (looplength)     ! To cmulate the debris length
    deallocate (lengthl)
    deallocate (segl)           ! A seg tab to save the loop debris result
    deallocate (segll)          ! A seg tab to save the loop debris sorting result from smaller to bigger
    deallocate (listlength)       ! The list of loop debris length
    deallocate (listlengthPrim)   ! The list of loop debris length in Prim direction
    deallocate (listlengthNormal) ! The list of loop debris length in Normal direction
    deallocate (lengthPrim)       ! The list of loop debris length in Prim direction
    deallocate (lengthNormal)     ! The list of loop debris length in Normal direction


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (4)                                                   !
  ! The dislocation microstructure (segment output file) is analysed  !
  ! to make statistics on jogs. Jogs are classified according to      !
  ! a mobility criteria (balance between the resolved PK on each side !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(4)

    ! Initialization
    Nsegloop = 0
    write(*,*) "What is the uniaxial applied stress axis 'X Y Z'?"
    read(*,*) Z(1:3)
    ampli = 1.
    ! THe reference applied stress tensor is defined
    call  ORIENTDIEG(AMPLI,Z,SIGAPP,.TRUE.)

    ! Results are saved in the file Jogs_stat
    open(1,file='../out/jogs_stat',STATUS='UNKNOWN')
    ! The file header
    write(1,*) '-->  Index     TauBO     TauBE     sum'

    ! The segment loop
    do i=1,nsegm

      ! We are looking fro screew segements ending jogs
      if (seg(i)%gd < iun) cycle

      ! The number of screw GD found
      Nsegloop = Nsegloop + 1

      ! First non zero length segment at O and E of i
      indexvois  = seg(i)%vnne
      indexvoiso = seg(i)%vnno

      ! The Burgers vector of i, vnno and vnne
      indexB = (int((seg(i)%veclin-1)/8))*8 + 1
      btmp(1:3) = bveclin(1:3,indexB)/normlin(indexB)

      ! The resolved applied shear stress at the O of seg GD i
      ! The PK force vector
      vectmp(1:3) = 0.
      do J = 1,3
        vectmp(1) = vectmp(1) + SIGapp(1,J)*Btmp(J)
        vectmp(2) = vectmp(2) + SIGapp(2,J)*Btmp(J)
        vectmp(3) = vectmp(3) + SIGapp(3,J)*Btmp(J)
      enddo
      PK(1) = (vectmp(2)*bveclin(3,seg(indexvoiso)%veclin) &
            - vectmp(3)*bveclin(2,seg(indexvoiso)%veclin))&
            / normlin(seg(indexvoiso)%veclin)
      PK(2) = (vectmp(3)*bveclin(1,seg(indexvoiso)%veclin) &
            - vectmp(1)*bveclin(3,seg(indexvoiso)%veclin))&
            / normlin(seg(indexvoiso)%veclin)
      PK(3) = (vectmp(1)*bveclin(2,seg(indexvoiso)%veclin) &
            - vectmp(2)*bveclin(1,seg(indexvoiso)%veclin))&
            / normlin(seg(indexvoiso)%veclin)
      ! The resolved force in the screw (GD) direction
      TauBO =((PK(1)*bveclin(1,indexB))+        &
              (PK(2)*bveclin(2,indexB))+        &
              (PK(3)*bveclin(3,indexB)))/normlin(indexB)
      if(TauBO < 1.d-10) TauBO= 0.

      ! The resolved applied shear stress at the E of seg GD i
      ! The PK force vector
      ! vectmp is unchanged
      PK(1) = (vectmp(2)*bveclin(3,seg(indexvois)%veclin) &
            - vectmp(3)*bveclin(2,seg(indexvois)%veclin))&
            / normlin(seg(indexvois)%veclin)
      PK(2) = (vectmp(3)*bveclin(1,seg(indexvois)%veclin) &
            - vectmp(1)*bveclin(3,seg(indexvois)%veclin))&
            / normlin(seg(indexvois)%veclin)
      PK(3) = (vectmp(1)*bveclin(2,seg(indexvois)%veclin) &
            - vectmp(2)*bveclin(1,seg(indexvois)%veclin))&
            / normlin(seg(indexvois)%veclin)
      ! The resolved force in the screw (GD) direction
      TauBE =((PK(1)*bveclin(1,indexB))+        &
              (PK(2)*bveclin(2,indexB))+        &
              (PK(3)*bveclin(3,indexB)))/normlin(indexB)
      if(TauBE < 1.d-10) TauBE= 0.

      ! Forces on this GD are saved
      ! Negative value means segments are both working but in opposite directions
      ! zero means only one segment is working on the GD segment, to move it
      ! Positive value means segments are both working in the same direction
      if (max(abs(TauBO),abs(TauBE)) /= 0.) then
        write(1,39) nsegloop,TauBO,TauBE,(abs(TauBO+TauBE)-max(abs(TauBO),abs(TauBE)))
      else
        write(1,39) nsegloop,TauBO,TauBE,0.
      endif

    enddo ! i seg loop

    print*,"The computation proceeded whithout problem!"

    39 format(1X,I6,3E13.5)
    close(1)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (5)                                                   !
  ! The dislocation microstructure (segment output file) is used      !
  ! to make a map of internal stress at the nodes of a plane grid.    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(5)

      call sigma_int ! Calculation of the stress tensor at coordinates

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (6)                                                   !
  ! analysis of the orientation of the junctions with respect to      !
  ! the displacement sens in a tested slip system                     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(6)

    ! The quantities calculated in this section
    !sensTlj(:): Statistic on segments touching a junction
    !  sensTlj(1) : number of segments having a junction as neighbor
    !  sensTlj(2) : number of segments, which zip the junction as an effect of the applied stress
    !  sensTlj(3) : number of segments, which unzip the junction as an effect of the applied stress
    !  sensTlj(4) : number of segments, which are neutral at the junction triple node as an effect of the applied stress
    !conf(:): configuration of the line when attached to 2 junctions
    !  conf(1) the two junctions are zipped by the mobile segment in between (it s harder and harder to glide)
    !  conf(2) the two junctions are unzipped by the mobile segment in between ( it s easier and easier to glide)
    !  conf(3) one junction is zipped,one unzipped ( the neutral case, the mobile dislocation glide at constant stress)
    !  conf(4) the 2 junctions are align (the mobile segment in between move for free in both directions)


    ! INITIALIZATION

    sensTlj(1:4) = 0
    conf(1:4)    = 0
    shortseg     = 0

    ! Selection of the system considered for this analysis
    write(*,*) "System of interest ?"
    read(*,*) sys

    sigapp(:,:) = 0.
    !stress tensor to define the sens of displacement of each segment
    write(*,*) "Please define the  applied stress tensor used to build the microstructure."
    write(*,*) "sig11, sig22, sig33, sig32, sig31, sig21   (Unit vector)."
    read(*,*) sigapp(1,1),sigapp(2,2),sigapp(3,3),sigapp(3,2),sigapp(3,1),sigapp(2,1)
    ! This tensor is symmetric
    sigapp(1,2) = sigapp(2,1)
    sigapp(1,3) = sigapp(3,1)
    sigapp(2,3) = sigapp(3,2)
    write(*,*) ' '
    write(*,*) 'sigapp', sigapp(1,:)
    write(*,*) 'sigapp', sigapp(2,:)
    write(*,*) 'sigapp', sigapp(3,:)

    ! The minimum number of connected segments to consider in the dislocation section between two junctions
    write(*,*) " "
    write(*,*) "Minimum Length of dislocation line section between two junctions considered in the stats (units micron) !"
    read(*,*) lensecmin
    if (lensecmin < 0) then
      write(*,*) "The lensecmin value is stupid !"
      stop
    endif
    lensecmin = lensecmin * 1.e-6 / avalue  ! lensecmin is redefined in simulation units

    ! Calculation of the Burgers vector associated to the system 'sys'
    ij=(sys-1)*8+1
    bi(1:3)=bveclin(1:3,ij)
    write(*,*) ' '
    write(*,*) 'tested slip system =', ij, 'Burgers vector =',bi(:)
    write(*,*) ' '

    ! The systems displacement sens
    sens(1:96) = 0
    do ij = 1,12
      do ii= 1,8
        li = (ij-1) * 8 + ii      ! Segment index
        vli(:) = bveclin(:,li)    ! Segment direction vector
        VDi(:) = Bvecdep(:,Li)    ! Segment displacement vector
        !write(*,*) ' '
        !write(*,*) '##ii', ii, 'li', li
        !write(*,*) 'vli', vli(:)
        !write(*,*) 'vdi', vdi(:)
        ! The P-K resolved stress calculation
        do l = 1,3
          G(l) = int(bi(1)*sigapp(1,l) + bi(2)*sigapp(2,l) + bi(3)*sigapp(3,l))
        enddo
        F(1) = G(2)*vli(3) - G(3)*vli(2)
        F(2) = G(3)*vli(1) - G(1)*vli(3)
        F(3) = G(1)*vli(2) - G(2)*vli(1)
        OiOjscaN = F(1)*vdi(1) + F(2)*vdi(2) + F(3)*vdi(3)
        sens(li) = 1
        if (OiOjscaN < 0) sens(li)=-1
        !write(*,*) 'f', f(:)
        !write(*,*) 'g', g(:)
        !write(*,*) 'sca',OiOjscaN
        !write(*,*) 'sens',sens(li)
      enddo
    enddo

    ! modulo are applied to segment origin to center the microstructure if needed
    modurx = modur(1)
    modury = modur(2)
    modurz = modur(3)
    seg(1:nsegm)%o(1) = modulo(seg(1:nsegm)%o(1),modurx)
    seg(1:nsegm)%o(2) = modulo(seg(1:nsegm)%o(2),modury)
    seg(1:nsegm)%o(3) = modulo(seg(1:nsegm)%o(3),modurz)

    ! Beginning of the calculation

    ! A loop on all slip system
    do jj = 1, 12

      ! A loop on all segments
      do ii = 1, nsegm

        ! We first eliminate all the useless segments
        if (seg(ii)%norme ==0) cycle                  ! The segment length is zero
        if (seg(ii)%jonc) cycle                       ! The segment is a junction segment
        indexB = (int((seg(ii)%veclin-1)/8))*8 + 1
        tmp = indexB/8
        if ((floor(tmp) + 1) /= sys) cycle            ! The segment is not part of the slip system we want to test

        ! The vnno and vnne of the tested segment
        vnno = seg(ii)%vnno
        vnne = seg(ii)%vnne

        ! We keep only segment ii having a junction as neighbor
        ! The analysis start with segments having neighbor on the O side
        if (seg(vnno)%jonc) then

          ! The tested segment properties
          Li      = SEG(ii)%VECLIN
          VDi(:)  = Bvecdep(:,Li) * sens(li)
          vli(:)  = Bveclin(:,li)
          Nvec(:) = Bvecnor(:,Li)
          ! The O and E coordinates of the junction segment touching segment ii
          Ej(:)   = seg(ii)%O
          OJ(:)   = Ej(:) - seg(vnno)%norme*bveclin(:,seg(vnno)%veclin)

          !!!!!!!!! Part 1 !!!!!!!!!!

          ! A new segment touching a junction is found
          sensTlj(1) = sensTlj(1) + 1
          sensTLjMemo = 0

          ! Calculation of P1 the junction center
          ! The initial vnn is redefined
          compteur = 1
          vnn = seg(vnno)%vnno
          !search for aligned neighboors
          do while(vnn /= nsegmax .and. seg(vnn)%veclin == seg(vnno)%veclin)
            compteur = compteur + 1
            OJ(:) = OJ(:) - seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)
            vnn = seg(vnn)%vnno
            if(compteur > 10) then
              exit    ! This straight segment is made with more than 10 segments (stupid)
            endif
          enddo
          P1(1:3) = (Oj + Ej) / 2

          ! Calculation of P2 the center of segment ii
          P2(:) = seg(ii)%o(:) + int(0.5*float(seg(ii)%norme*vli(:)))

          ! We explore the other side of segment to calculate the curvature at ii
          ! The properties of the first segment on the E side of ii
          Lj      = SEG(vnne)%VECLIN
          VLJ(:)  = Bveclin(:,Lj)
          OJ(:)   = seg(ii)%O + seg(ii)%norme*vli(:)
          EJ(:)   = Oj(:) + seg(vnne)%norme*bveclin(:,seg(vnne)%veclin)

          ! calculation of P3 the center of the straight line section on the E side of ii
          ! The initial vnn is redefined
          compteur = 1
          vnn = seg(vnne)%vnne
          !search for aligned neighboors
          do while(vnn /= nsegmax .and. seg(vnn)%veclin == seg(vnne)%veclin)
            compteur = compteur + 1
            EJ(:) = EJ(:) + seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)
            vnn = seg(vnn)%vnne
            if(compteur > 10) then
              exit    ! This straight segment is made with more than 10 segments (stupid)
            endif
          enddo
          P3(1:3) = (Oj + Ej) / 2

          ! Calculation of the curvature vector VCOURB from P1,P2,P3
          ! The aim of this calculation is to calculate the real vector normal to the dislocation line section
          ! discretized with segment ii. The discrete segment character is different from the dislocation character
          BUFTLRO(1:3,1) = P1(1:3)
          BUFTLRO(1:3,2) = P2(1:3)
          BUFTLRO(1:3,3) = P3(1:3)
          VCOURB(1:3) = VECCOURB(BUFTLRO(1:3,1),BUFTLRO(1:3,2),BUFTLRO(1:3,3),TESTPMPQ,alpha)

          ! The spacial case where VCOURB(1:3) cannot correctly be defined (see mM simulation code for more details)
          ! Better to stop the analysis here for segment ii
          if (TESTPMPQ) then
            sensTlj(1) = sensTlj(1) - 1
            cycle
          endif

          ! The VCOURB vector is forced to be oriented in the applied stress direction of displacement
          if ((VCOURB(1)*VDI(1) + VCOURB(2)*VDI(2) + VCOURB(3)*VDI(3)) < 0) then
            !if sca vcourb with vdi <0 then vcourb is not correctly
            !oriented with respect to the applied stress tensor
            ! print *, 'PB sens VCOURB', vcourb(:),'vdi',vdi(:)
            VCOURB(:) = -VCOURB(:)
          else if (VCOURB(1)==0 .and. VCOURB(2)==0 .and. VCOURB(3)==0) then
            ! if vcourb=0 then P1,P2,P3 are aligned
            ! print *, 'seg aligne vcourbe = vdi'
            VCOURB(:) = vdi(:)
          endif

          ! The vector // to the junction direction.
          ! This vector sens is signed to decrease the junction length
          Lj        = SEG(vnno)%VECLIN
          VLJ(:)    = -Bveclin(:,Lj)

          ! Is the action of the applied stress on ii wants to zip or unzip the junction
          OiOjscaN  = VCOURB(1)*VLJ(1) + VCOURB(2)*VLJ(2) + VCOURB(3)*VLJ(3)
          if (OiOjscaN < 0) then
            sensTlj(2) = sensTlj(2) + 1   ! applied stress on ii works to zip the junction
            sensTLjMemo = 2
          elseif (OiOjscaN > 0) then
            sensTlj(3) = sensTlj(3) + 1   ! applied stress on ii works to unzip the junction
            sensTLjMemo = 3
          else
            sensTlj(4) = sensTlj(4) + 1   ! applied stress on ii is neutral
            sensTLjMemo = 4
          endif
!           print *, '>vnno jonc', vnno,Lj,VLJ
!           print *, 'P1',P1(:),'P2',P2(:),'p3',p3(:)
!           print *, 'vcourb',VCOURB, 'vlj',vlj
!           print *, 'scaN',OiOjscaN

          !!!!!!!!! Part 2 !!!!!!!!!!

          !Now we follow the line in the E direction of ii to find the next junction
          ij = 0
          lengthcumul = seg(Li)%norme*NormLin(seg(Li)%veclin)   ! The length calculation is initiated to the first segment length
          !we are only interested by junction segments
          do while (.not.seg(vnne)%jonc .and. vnne/=ii .and. vnne/=nsegmax)
            ! The dislocation length between two junctions is summed
            lengthcumul = lengthcumul + seg(vnne)%norme*NormLin(seg(vnne)%veclin)
            vnne=seg(vnne)%vnne
            ij = ij + 1
            if (ij > nsegmax) then
              write(*,*) 'a problem is found in the segment connectivity'
              stop
            endif
          enddo

          ! Spacial cases we do not want to consider
          if (vnne == ii) then
            sensTlj(1) = sensTlj(1) - 1 ! the previous segment on the O side is finaly not considered
            cycle                       ! A closed dislocation loop
          endif
          if (vnne == nsegmax) then
            sensTlj(1) = sensTlj(1) - 1 ! the previous segment on the O side is finaly not considered
            cycle                  ! pinning point
          endif

          ! A new segment touching a junction is found
          sensTlj(1) = sensTlj(1) + 1

          ! Properties of the new tested segment touching a junction on the E side
          Lj = SEG(vnne)%VECLIN
          VLJ(:) = Bveclin(:,Lj)
          OJ(:) = seg(vnne)%O
          EJ(:) = Oj(:) + seg(vnne)%norme*bveclin(:,seg(vnne)%veclin)

          !definition of P3 the center of the junction on the E side
          ! The initial vnn is redefined
          compteur = 1
          vnn = seg(vnne)%vnne
          !search for aligned neighbors
          do while(vnn /= nsegmax .and. seg(vnn)%veclin == seg(vnne)%veclin)
            compteur = compteur + 1
            EJ(:) = EJ(:) + seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)
            vnn = seg(vnn)%vnne
            if(compteur > 10) then
              exit    ! This straight segment is made with more than 10 segments (stupid)
            endif
          enddo
          P3(1:3) = (Oj + Ej) / 2

          ! P2 is the center of seg(vnne)%vnno
          vnn = seg(vnne)%vnno
          Lj = SEG(vnn)%VECLIN
          VLJ(:) = Bveclin(:,Lj)
          P2(:) = seg(vnne)%o(:) - int(0.5*float(seg(vnn)%norme*VLJ(:)))

          ! We explore the other side of segment touching a junction on the E side
          Ej(:) = P2(:) - int(0.5*float(seg(vnn)%norme*vlj(:)))
          vnn = seg(seg(vnne)%vnno)%vnno
          Lj = SEG(vnn)%VECLIN
          VLJ(:) = Bveclin(:,Lj)
          OJ(:) = Ej(:) - seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)

          !definition of P1 the center of the vnno
          ! The initial vnn is redefined
          compteur = 1
          vnno = seg(seg(vnne)%vnno)%vnno
          vnn = seg(vnno)%vnno
          !search for aligned neighboors
          do while(vnn /= nsegmax .and. seg(vnn)%veclin == seg(vnno)%veclin)
            compteur = compteur + 1
            OJ(:) = OJ(:) - seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)
            vnn = seg(vnn)%vnno
            if(compteur > 10) then
              exit    ! This straight segment is made with more than 10 segments (stupid)
            endif
          enddo
          P1(1:3) = (Oj + Ej) / 2

          Li= SEG(seg(vnne)%vnno)%VECLIN
          VDi(:) = Bvecdep(:,Li) * sens(li)

          ! Calculation of the curvature vector VCOURB from P1,P2,P3
          ! The aim of this calculation is to calculate the real vector normal to the dislocation line section
          BUFTLRO(1:3,1) = P1(1:3)
          BUFTLRO(1:3,2) = P2(1:3)
          BUFTLRO(1:3,3) = P3(1:3)
          ! calculation of the curvature vector from 3 pt P1,P2,P3
          VCOURB(1:3) = VECCOURB(BUFTLRO(1:3,1),BUFTLRO(1:3,2),BUFTLRO(1:3,3),TESTPMPQ,alpha)

          ! The spacial case where VCOURB(1:3) cannot correctly be defined (see mM simulation code for more details)
          ! Better to stop the analysis here
          if (TESTPMPQ) then
            sensTlj(1) = sensTlj(1) - 2    ! The last two segment are not considered
            sensTLj(sensTLjMemo) = sensTLj(sensTLjMemo) -1  ! We remove the information on the O side
            cycle
          endif

          ! The VCOURB vector is forced to be oriented in the applied stress direction of displacement
          if ((VCOURB(1)*VDI(1) + VCOURB(2)*VDI(2) + VCOURB(3)*VDI(3)) < 0) then
            !if sca vcourb with vdi <0 then vcourb is not correctely
            !oriented with respect to the applied stress tensor
            ! print *, 'PB sens VCOURB', vcourb(:),'vdi',vdi(:)
            VCOURB(:) = -VCOURB(:)
          else if (VCOURB(1)==0.and.VCOURB(2)==0.and.VCOURB(3)==0) then
            ! if vcourb=0 then P1,P2,P3 are aligned
            ! print *, 'seg aligne vcourbe=vdi'
            VCOURB(:) = vdi(:)
          endif

          ! The vector // to the junction direction.
          ! This vector sens is signed to decrease the junction length
          Lj = SEG(vnne)%VECLIN
          VLJ(:) = Bveclin(:,Lj)

          OiEjscaN = VCOURB(1)*VLJ(1)+VCOURB(2)*VLJ(2)+VCOURB(3)*VLJ(3)
          if (OiEjscaN < 0) then
            sensTlj(2) = sensTlj(2) + 1     ! applied stress on ii works to zip the junction
          elseif (OiEjscaN > 0) then
            sensTlj(3) = sensTlj(3) + 1     ! applied stress on ii works to unzip the junction
          else
            sensTlj(4) = sensTlj(4) + 1     ! applied stress on ii works to unzip the junction
          endif
!           print *, '>vnne jonc', vnne,Lj,VLJ
!           print *, 'P1',P1(:),'P2',P2(:),'p3',p3(:)
!           print *, 'vcourb',VCOURB, 'vlj',vlj
!           print *, 'scaN',OiEjscaN

          ! The finale analysis is not made on very short line section ( not realy representative )
          if (lengthcumul > lensecmin) then
            if (OiEjscaN==0 .and. OiOjscaN==0) then
                 conf(4) = conf(4) + 1                                ! two // junctions in a raw
            else if (OiojscaN==0 .and. OiEjscaN/=0) then
              if (OiEjscaN < 0) conf(1) = conf(1) + 1                   ! A hardening geometry
              if (OiEjscaN > 0) conf(2) = conf(2) + 1                   ! A softening geometry
            else if (OiEjscaN==0 .and. OiOjscaN/=0) then
              if (OiOjscaN < 0) conf(1) = conf(1) + 1                   ! A hardening geometry
              if (OiOjscaN > 0) conf(2) = conf(2) + 1                   ! A softening geometry
            else if (OiEjscaN/=0 .and. OiOjscaN/=0) then
              if (OiOjscaN < 0 .and. OiEjscaN < 0) conf(1) = conf(1) + 1  ! A hardening geometry
              if (OiOjscaN > 0 .and. OiEjscaN > 0) conf(2) = conf(2) + 1  ! A softening geometry
              if (OiOjscaN*OiEjscaN<0) conf(3) = conf(3) + 1          ! The indeterminate case
            endif
          else
            shortseg = shortseg + 1
          endif

        endif   ! test on segment connected to a junction
      enddo   ! loop on segments
    enddo   ! loop on slip systems

    tmp = real((conf(1)+conf(2)+conf(3)+conf(4)), DP)
    if ((conf(1)+conf(2)+conf(3)+conf(4)+shortseg) /= sensTlj(1)/2) then
      write(*,*) " A problem in found in the case(6) calculation)"
      stop
    endif

    write(*,*) "##############################################################"
    write(*,*) "total number of tested line sections pinned by two junctions =", sensTlj(1)/2
    write(*,*) " "
    write(*,*) "Number of triple node configurations tested      =", sensTlj(1)
    write(*,*) "Number of segments working to zip a junction     =", sensTlj(2)
    write(*,*) "Number of segments working to unzip a junction   =", sensTlj(3)
    write(*,*) "Number of segments neutral on a junction         =", sensTlj(4)
    write(*,*) " "
    write(*,*) "Hard forest configurations                       =", conf(1)
    write(*,*) "Soft forest configurations                       =", conf(2)
    write(*,*) "Neutral or undetermined forest configurations    =", conf(3:4)
    write(*,*) "Number of short segment not analysed             =", shortseg
    write(*,*) " "
    write(*,*) "Forest configuration statistic :    ", "Hard",conf(1)/tmp,"Soft",conf(2)/tmp,"Others",(conf(3)+conf(4))/tmp
    write(*,*) "##############################################################"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (7)                                                   !
  ! Analysis of the distribution of the dislocation lines character   !
  ! in the dislocation microstructure. Output of this calculation is  !
  ! a list of of segment character and length in file statcarac       !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(7)

    !modulo are applied to segment origin to reconstruct
    !properly the microstructure
    modurx = modur(1)
    modury = modur(2)
    modurz = modur(3)
    seg(1:nsegm)%o(1) = modulo(seg(1:nsegm)%o(1),modurx)
    seg(1:nsegm)%o(2) = modulo(seg(1:nsegm)%o(2),modury)
    seg(1:nsegm)%o(3) = modulo(seg(1:nsegm)%o(3),modurz)

    open(111,FILE='../out/statcarac',STATUS='UNKNOWN')

    !initialization
    allocate (anglevis(nsegm))
    allocate (vecnorlin(3,nbase))
    allocate (assoc(nbase,8))
    TWOPII = 2.*PII
    do II=1,NBASE
      do IJ=1,8
        ASSOC(II,IJ) = IJ + ((II-1)/8)*8 ! vecteur IJ associe
      enddo
    enddo

    ! The normalized vector in the segment direction
    do ij= 1, nbase
      VecNorLin(1:3,ij) = normaivec(int(BVECLIN(1:3,ij),DPI))
    enddo

    !LOOP ON SEGMT II
    do ii = 1,nsegm
      ! if ii is null or jonc or do not belong to the system of interest jj => cycle
      if (seg(ii)%norme ==0) cycle
      indexB = (int((seg(ii)%veclin-1)/8))*8 + 1
      tmp = indexB/8
      sys = floor(tmp) + 1

      vnno = seg(ii)%vnno
      vnne = seg(ii)%vnne

      Li= SEG(ii)%VECLIN
      VDi(:) = Bvecdep(:,Li)
      vli(:) = bveclin(:,li)
      normvli= dsqrt(real(vli(1)**2+vli(2)**2+vli(3)**2,DP))
      normalvec(:) = normaivec(int(Bvecnor(:,Li),DPI))
      !we look for the neighboor at the origin of ii
      compteur = 1
      Ej(:) = seg(ii)%O
      OJ(:) = Ej(:) - seg(vnno)%norme*bveclin(:,seg(vnno)%veclin)

      IF (seg(vnno)%JONC) THEN
        ! The special case where the vnno is junction
        P1(1:3) =seg(ii)%o(:)
      ELSE
        !definition of P1 the center of the vnno
        ! The initial vnn is redefined
        vnn = seg(vnno)%vnno
        !search for aligned neighboors
        do while(vnn /= nsegmax .and. seg(vnn)%veclin == seg(vnno)%veclin)
          compteur = compteur + 1
          OJ(:) = OJ(:) - seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)
          vnn = seg(vnn)%vnno
          if(compteur > 10) then
            exit
          endif
        enddo
        P1(1:3) = (Oj + Ej) / 2
      ENDIF

      ! P2 is the center of ii
      P2(:) = seg(ii)%o(:)+ int(0.5*float(seg(ii)%norme*vli(:)))

      ! we now look for the neighboor at the extremity of ii
      compteur = 1
      Lj = SEG(vnne)%VECLIN
      VLJ(:) = Bveclin(:,Lj)
      OJ(:) = seg(ii)%O(:)+seg(ii)%norme*vli(:)
      EJ(:) = Oj(:) + seg(vnne)%norme*bveclin(:,seg(vnne)%veclin)

      IF (seg(vnne)%JONC) THEN
        ! The special case where the vnno is junction
        P3(1:3) =seg(ii)%o(:)+seg(ii)%norme*vli(:)
      ELSE
        !definition of P3 the center of the vnne
        ! The initial vnn is redefined
        vnn = seg(vnne)%vnne
        !search for aligned neighboors
        do while(vnn /= nsegmax .and. seg(vnn)%veclin == seg(vnne)%veclin)
          compteur = compteur + 1
          EJ(:) = EJ(:) + seg(vnn)%norme*bveclin(:,seg(vnn)%veclin)
          vnn = seg(vnn)%vnne
          if(compteur > 10) then
            exit
          endif
        enddo
        P3(1:3) = (Oj + Ej) / 2
      ENDIF


      ! calculation of the curvature vector from 3 pt P1,P2,P3
      BUFTLRO(1:3,1) = P1(1:3)
      BUFTLRO(1:3,2) = P2(1:3)
      BUFTLRO(1:3,3) = P3(1:3)
      VCOURB(1:3)=VECCOURB(BUFTLRO(1:3,1),BUFTLRO(1:3,2),BUFTLRO(1:3,3),TESTPMPQ,alpha)

      if (TESTPMPQ) then
        ! Le rayon de courbure est suppose etre le demi diametre du cercle
        ! passant par P1-P3 (ou M-Q)
        VCOURB(1:3)= VDI(1:3)
      endif

      ! La norme de ce vecteur est le rayon de courbure
      RAYON2=VCOURB(1)*VCOURB(1)+VCOURB(2)*VCOURB(2)+VCOURB(3)*VCOURB(3)

      ! Si les trois points sont alignes, vcourb = vdi
      if (rayon2==0) then
        vcourb(:) = vdi(:)
        rayon2 = VCOURB(1)*VCOURB(1)+VCOURB(2)*VCOURB(2)+VCOURB(3)*VCOURB(3)
      endif

      ! Angle de la tangente locale a la direction vis
      ! (calcul utile seulement dans le cas anisotrope)
      ANGLE=VCOURB(1)*vecnorlin(1,assoc(Li,1))+   &
            VCOURB(2)*vecnorlin(2,assoc(Li,1))+   &
            VCOURB(3)*vecnorlin(3,assoc(Li,1))

      ANGLE=ANGLE/dsqrt(rayon2)

      ! To define the angle between -Pi and Pi an additional test with respect to the
      ! slip plane normale is needed
      normallinvec(1) = VCOURB(2)/dsqrt(rayon2)*vecnorlin(3,assoc(Li,1)) -     &
                    VCOURB(3)/dsqrt(rayon2)*vecnorlin(2,assoc(Li,1))
      normallinvec(2) = VCOURB(3)/dsqrt(rayon2)*vecnorlin(1,assoc(Li,1)) -     &
                    VCOURB(1)/dsqrt(rayon2)*vecnorlin(3,assoc(Li,1))
      normallinvec(3) = VCOURB(1)/dsqrt(rayon2)*vecnorlin(2,assoc(Li,1)) -     &
                    VCOURB(2)/dsqrt(rayon2)*vecnorlin(1,assoc(Li,1))

      signtest = normallinvec(1)*normalvec(1)+normallinvec(2)*normalvec(2)+normallinvec(3)*normalvec(3)

      ! angle between the normal to I and b
      if (dabs(angle) > UN) then
        print*,'strange',angle
        angle = dsign(UN,angle)
      endif

      anglevis(ii) = modulo(dble(PII) * HALF + sign(UN,signtest)*dacos(ANGLE),TWOPII)
!       ANGLE       = dble(PII)      * HALF - dacos(ANGLE)


      31 format(1X,1I6,2X, E13.6,2X, E13.6)
      write (111,31) ii, anglevis(ii),seg(ii)%norme*avalue*normvli
      print*, ii, anglevis(ii),seg(ii)%norme*avalue*normvli

    enddo

    deallocate (anglevis)
    close(111)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (11)                                                  !
  ! NYE TENSOR AND GND approximate density on a lattice               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(8)

    ! The nye tensor and approximate gnd density pattern is calculated on a lattice
    ! The lattice is used to defines boxes where the dislocation density is locally calculated
    ! The lattice is defined with :
    ! 1) Origine coordinate - Lattice_O
    ! 2) 3 vectors defining the elementary lattice volume shape and dimension - v1, v2, v3
    ! 3) 3 integer numbers defining the dimension of the lattice - n1, n2, n3


    !Whatever append before, deallocate variables
    deallocate(bveclin)
    deallocate(bvecnor)
    deallocate(seg)
    deallocate(normlin)

    !---------------------------------------------------------

    ! The simulation input files are defined
    open(1,file='../in/input.dd',STATUS='OLD')

    carac = " "
    do while (carac /= "-")
      read (1,*) carac
    enddo

    read (1,*) materiau
    read (1,*) control
    read (1,*) segments

    close (1)

    !---------------------------------------------------------

    ! Some material parameter must be loaded
    open(1,file="../in/"//materiau,STATUS='OLD')

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
    read(1,*) crystal_structure(1:3) ! the crystal symmetry
    read(1,*) useless
    read(1,*) VecBurgers    ! The burgers vector in Angstrom
    VecBurgers = VecBurgers * 1.d-10

    close (1)

    !---------------------------------------------------------

    ! Shall we work with 'Seg_save' or 'film.bin'
    write(*,*) "Do you want to analyse the -Seg_save- or -film.bin- file"
    write(*,*) "F = Seg_save, T = film.bin"
    read (*,*) Film

    if (Film) then

      print *, " > We are going to analyse dislocation configuration in -film.bin-"
      ! The file we want to analyze
      fichier = "film.bin"
      open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
      ! Important information at the begining of the film.bin file
      read(50) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
      read(50) nsegmaxf     ! The nsegmax value used in this simulation
      read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
      read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry

      allocate  (bvecnor(3,nbase))
      allocate  (bveclin(3,nbase))

      read(50) (bveclin(:,I),bvecnor(:,I), I=1,nbase)   ! table of vector direction and displacement for this simulation
      read(50) modur(:)                                 ! Dimension of the 3D simulated volume in a value

      ! Test on the nsegmax definition
      if (nsegmax /= nsegmaxf) then
        write(*,*) "Number of segments in film.bin = ",nsegmaxf,"but in this programm = ",nsegmax
        write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
      !  stop
      endif

    else

      print *, " > We are going to analyse dislocation configuration in -Seg_save-"

      ! The simulation elementary lattice parameter
      cristallo = '../out/BVD.'//crystal_structure(1:3)

      open (1,FILE=cristallo,STATUS='OLD')

      read (1,*) nbase , avalue

      allocate (bvecnor(3,nbase))
      allocate (bveclin(3,nbase))

      do j = 1,nbase
          read (1,*) i, bveclin(1:3,j)
          if(j /= i) print *, "stop, a problem appears when reading BVD"
      enddo

      close(1)

      open(50,file='../in/SEG_save',STATUS='OLD')

      read(50,*) useless
      read(50,*) Nsegm                       ! Nb of segments

      ! Test on the nsegmax definition
      if (nsegm > nsegmax) then
        write(*,*) "Number of segments in -Seg_save- = ",nsegmaxf,"but in this program = ",nsegmax
        write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
       ! stop
      endif

      read(50,*) modur(:)  ! Dimensions of the simulation cell

      allocate (seg(nsegm))
      allocate  (iseg(5,nsegm))
      allocate  (tabvois(2,nsegm))
      allocate  (junctt(nsegm))
      allocate  (ID_surface(nsegm))

      do I = 1,nsegm

        read (50,*) j,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
                   seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%Ijonc,seg(i)%gd

        print*,j,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
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

    !---------------------------------------------------------

    ! The normal of slip planes is defined
    allocate (normlin(nbase))
    normlin(:) = 0
    do j = 1,nbase
        normlin(j)=sqrt(real(bveclin(1,j))**2 + real(bveclin(2,j))**2 + real(bveclin(3,j))**2)
    enddo

    !---------------------------------------------------------

    ! The configuration(s) we want to analyse is defined

    If (Film) then
      ! We first select the part of the film.bin file we want to use for the stats
      write(*,*) "Analyse start at simulation step:  ? "
      read (*,*) debutaffich
      write(*,*) "Analyse end at simulation step:  ? "
      read (*,*) finaffich
      write(*,*) "Snapshot frequency for the film analysis: (Warning : must be coherent with KSTATS) ? "
      read (*,*) freqaffich
      if (freqaffich  <= 0)freqaffich = 1
    else !set 1 for all parameters in case of SEG_save
      debutaffich = 1
      finaffich   = 1
      freqaffich  = 1
    endif

    nb_step_stat = 0    ! The variable counting the number of configuration analysed is initialized

    ! The slip system symmetry is defined
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

    !---------------------------------------------------------

    ! Definition of the calculation lattice
    write(*,*) "We now defined the parameters used to defined tha calculation lattice!"
    write(*,*) " "
    write(*,*) " Please give the coordinates of the Lattice origine in a unit : Ox  Oy  Oz"
    read (*,*) Lattice_O(1:3)
    write(*,*) " >>>> ", Lattice_O(1:3)
    write(*,*) " "
    write(*,*) " Please give the 3 vectors defining the Lattice elementary cell (must be in a unit)!"
    write(*,*) " "
    write(*,*) " V1x  V1y  V1y ?"
    read (*,*) V1(1:3)
    write(*,*) " V2x  V2y  V2y ?"
    read (*,*) V2(1:3)
    write(*,*) " V3x  V3y  V3y ?"
    read (*,*) V3(1:3)
    write(*,*) " >>>> V1 = ", V1(1:3)
    write(*,*) " >>>> V2 = ", V2(1:3)
    write(*,*) " >>>> V3 = ", V3(1:3)
    write(*,*) " "
    write(*,*) " Please give the 3 integer number defining the lattice dimension : n1  n2  n3"
    read (*,*) n1, n2, n3
    write(*,*) " "
    write(*,*) "Please enter the system of interest"
    write(*,*) "0 = all slip systems"
    read (*,*) sys
    write(*,*) " "
    !---------------------------------------------------------

    allocate (nye_tensor(1:n1,1:n2,1:n3,1:3,1:3))
    allocate (dens_GND(1:n1,1:n2,1:n3))

    ! The simulation box volume and lattice elementary cell
    A(1:3) = real(v1(1:3),dp)
    B(1:3) = real(v2(1:3),dp)
    C(1:3) = real(v3(1:3),dp)

    volb    = Determinant(A(:), B(:), C(:)) * avalue**3
    voltot  = volb * n1 * n2 * n3

    !---------------------------------------------------------

    ! Setup the result directory
    call date_and_time(date,time)!,zone,values)
    mydate=date//'_'//time(1:6)

    command = 'mkdir ../out/Nye_GND_dir'//'_'//mydate
    call system(command)

    ! WE START THE LOOP ON SIMULATION STEPS SAVED IN FILM.BIN

    print *, "-------------------------"
    print *, " Beginning of step loop "
    print *, "-------------------------"
    print *, ""

    ! INITIALIZATION ON THE LOOP ON THE SELECTED STEP
    kk=0

    Steps_loop_33: do while (kk >= 0 .and. kk <= finaffich)

      TheGoodStep = .false.

      if (Film) then

        ! A new simulation step
        read(50,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO

        TheGoodStep = (kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0)

        if (TheGoodStep) then

          print*," > Selected step "
          write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I5)') KK,EPSO,XSTRESS,nsegm

          allocate (seg(nsegm))
          allocate (iseg(5,nsegm))
          allocate (junctt(nsegm))
          allocate (tabvois(2,nsegm))
          allocate (ID_surface(nsegm))

          ! The segments micro structure is loaded
          read(50,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),ID_surface(j), J=1,nsegm)

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

        else

          ! We pass this step
          read(50,IOSTAT=KODERR)

        endif

        if(KODERR < 0) then
          print*,"  End of the film.bin file "
          exit Steps_loop_33
        elseif(KODERR > 0) then
          print*,"Reading problem 3",KODERR,kk
          stop
        endif

        ! The mark for the end of step is read
        read(50,IOSTAT=KODERR) caractest

        if(KODERR < 0) then
            print*,"  End of the film.bin file 3"
            exit Steps_loop_33
        elseif(KODERR > 0) then
            print*,"Reading error - Could not find the mark for the end of step ",caractest
            stop
        endif

      else

        kk = 1 ! Seg_save

      endif

    !---------------------------------------------------------

      ! The steps we want to use for statistic
      if (TheGoodStep) then ! Begin of the "Big If"

        ! Start of the calculations

        Nye_tensor(:,:,:,:,:) = 0
        dens_GND(:,:,:) = 0

          !Loop on the segments
          do i = 1, nsegm

            normi = seg(i)%norme

            ! Only segments with a non-zero length are of interest
            if (normi /= 0) then

              ! Identification of the segment slip system
              veclin = seg(i)%veclin
              indexB = (int((veclin-1)/8))*8 + 1
              tmp = indexB/8
              sys = floor(tmp)+1
              if (sys/=0 .and. floor(tmp)+1/=sys) cycle !Select slip system

              Oa(1:3) = seg(i)%O(:)

              ! precise computation of the density imposes to decompose the length of segments
              do J = 1,normi
                !J = normi
                vli(1:3) = bveclin(1:3,veclin)

                center(:) = real(Oa(:)) + (real(J)-0.5) * real(vli(:))
                center(:) = modulo(center(:), real(modur(:)))

                ! Vector connecting the density latice origin and the points along the segment length
                ra(:)   = center(:) - real(Lattice_O(:))

                D(1:3)  = ra(1:3)

                ! Calculation of the coordinate of ra() in the density grid
                ! Using the Cramer's rule we solve a linear system :
                ! v1 * fac1 + v2 * fac2 + v3 * fac3 = ra
                sysdet = Determinant(A, B, C)

                if (Abs(sysdet) < 10e-15) Then
                  coord = (/-300000,-300000,-300000/)
                else
                  xdet = Determinant(D, B, C)
                  ydet = Determinant(A, D, C)
                  zdet = Determinant(A, B, D)

                  coord(1)=xdet/sysdet
                  coord(2)=ydet/sysdet
                  coord(3)=zdet/sysdet
                endif

                ! coordinate of the segment section in the density tab
                c1 = nint(coord(1)) + 1
                c2 = nint(coord(2)) + 1
                c3 = nint(coord(3)) + 1

                ! We must test if the segment section is outside the grid
                if (c1 < 1 .or. c2 < 1 .or. c3 < 1 .or. c1 > n1 .or. c2 > n2 .or. c3 > n3) cycle

                !Determination of vectors L (line vector) and B (Burgers vector)
                !tmp = (seg(i)%veclin-1)/8
                !indexB = (int(tmp))*8 + 1 ! Burgers vector identification
                ! Junction special case
                if (seg(i)%jonc) then
                  b1(:)= bveclin(:,indexB)
                  !j=seg(i)%ijonc
                  tmp2 = (seg(seg(i)%ijonc)%veclin-1)/8
                  !indexB = (int(tmp))*8 + 1
                  b2(:)=bveclin(:,(int(tmp2))*8 + 1)
                  bj1(:)=b1(:)+b2(:)
                  bj2(:)=b1(:)-b2(:)
                  if (bj1(1)**2+bj1(2)**2+bj1(3)**2 < bj2(1)**2+bj2(2)**2+bj2(3)**2) then
                    btmp(:)=bj1(:)
                  else
                    btmp(:)=bj2(:)
                  endif
                  norm=sqrt(real(btmp(1))**2+ real(btmp(2))**2+ real(btmp(3))**2)
                  if (norm == 0) cycle
                  btmp(1:3) = btmp(:)/norm
                  !indexb = 0

                  ! Unjoined segments
                else
                  norm=sqrt(real(bveclin(1,indexb))**2+  &
                  real(bveclin(2,indexb))**2+  &
                  real(bveclin(3,indexb))**2)
                  btmp(1:3) = bveclin(1:3,indexB)/norm
                endif


                ! In the Nye tensor calculation, the dislocation line vector is unitary.
                ! The burgers vector is not unitary
                normvli=sqrt(real(vli(1))**2 + real(vli(2))**2 + real(vli(3))**2)
                vltmp(:)=vli(:)/normvli

                !Calculation of the Nye tensor
                do l = 1,3
                  do m = 1,3
                    if (seg(i)%jonc) then
                        !Junction: as we have 2 segments
                        !he contribution is divided by 2
                        Nye_tensor(c1,c2,c3,l,m) = Nye_tensor(c1,c2,c3,l,m) &
                                                 & + 0.5*btmp(l)*VecBurgers*vltmp(m)*(NormLin(veclin)*avalue)/volb
                    else

                        Nye_tensor(c1,c2,c3,l,m) = Nye_tensor(c1,c2,c3,l,m) &
                                                 &  +  btmp(l)*VecBurgers*vltmp(m)*(NormLin(veclin)*avalue)/volb
                    endif
                  enddo ! m
                enddo ! l


              enddo       ! loop on the segment length (j = 1,normi)

            endif         ! Non zero test

          enddo           ! loop on the segment number (i = 1, nsegm)

          ! Dislocation GND density evaluation based on the Nye tensor L1 norm (see Ultramicroscopy 133 (2013) 8–15)
          do c1 = 1,n1
            do c2 = 1,n2
              do c3 = 1,n3
dens_GND(c1,c2,c3) =  ( abs(Nye_tensor(c1,c2,c3,1,1))+abs(Nye_tensor(c1,c2,c3,1,2)) + abs(Nye_tensor(c1,c2,c3,1,3)) + &
                        abs(Nye_tensor(c1,c2,c3,2,1))+abs(Nye_tensor(c1,c2,c3,2,2)) + abs(Nye_tensor(c1,c2,c3,2,3)) +  &
                        abs(Nye_tensor(c1,c2,c3,3,1))+abs(Nye_tensor(c1,c2,c3,3,2)) + abs(Nye_tensor(c1,c2,c3,3,3))) &
                      & / VecBurgers
              enddo
            enddo
          enddo

        ! Management of the calculation outputs
        nb_step_stat = nb_step_stat +1                    ! The filename number

    1153 format ('../out/Nye_GND_dir','_',A15,'/dislocations_nye_incr',I0, '.txt') ! The file name construction

        write (filename, 1153) mydate, nb_step_stat       ! definition of a different filename for each loop
        open (unit=212, file=filename, status='unknown')  ! The file name is open

        ! the total density results

        !Writing of results for this step
        ! column 1: indice of box along the v1 direction
        ! column 2: indice of box along the v2 direction
        ! column 3: indice of box along the v3 direction
        ! column 4->ntsg+4: dislocation/junction density on the box

        print*,"Writing of density in the output file out/Nye_GND_dir_xxx/"
        write(212,'(2x,A2,2x,A2,2x,A2,2x, 2x,14A9)') 'X','Y','Z','Nye11','Nye12','Nye13',&
& 'Nye21','Nye22','Nye23','Nye31','Nye32','Nye33','rho GND'

        do c1 = 1, n1
          do c2 = 1,n2
            do c3 = 1,n3
            write(212,*) c1, c2, c3, Nye_tensor(c1,c2,c3,1,:),Nye_tensor(c1,c2,c3,2,:),Nye_tensor(c1,c2,c3,3,:),dens_GND(c1,c2,c3)
            enddo
          enddo
        enddo

        close(212)

        ! END OF THE CALCULATIONS MADE ON THIS STEP
        !  ======>>>>>

        deallocate (seg)
        deallocate (iseg)
        deallocate (junctt)
        deallocate (tabvois)
        deallocate (ID_surface)

      endif ! End of the "big if" selecting the step

      ! In the Seg_Save case only one pass is made
      if (.not. Film) kk = kk + 1

    enddo Steps_loop_33

    print *, ""
    print *, "-------------------------"
    print *, " End of step loop "
    print *, "-------------------------"
    print *, ""

    print *,' > Results are in /out/Nye_GND_dir'

    close(50)

    print*," > Nye tensor and APPROXIMATED GND density calculation achieved with success"

    deallocate (nye_tensor)
    deallocate (dens_GND)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (9)                                                   !
  ! TOTAL DISLOCATION DENSITY IN A SIMULATION THINFOIL                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(9)

    ! The dislocation density pattern is calculated on a lattice
    ! The lattice is used to defines boxes where the dislocation density is locally calculated
    ! The lattice is defined with :
    ! 1) Origine coordinate - Lattice_O
    ! 2) 3 vectors defining the elementary lattice volume shape and dimension - v1, v2, v3
    ! 3) 3 integer numbers defining the dimension of the lattice - n1, n2, n3


    !Whatever append before, deallocate variables
    deallocate(bveclin)
    deallocate(bvecnor)
    deallocate(seg)
    deallocate(normlin)

    !---------------------------------------------------------

    ! The simulation input files are defined
    open(1,file='../in/input.dd',STATUS='OLD')

    carac = " "
    do while (carac /= "-")
      read (1,*) carac
    enddo

    read (1,*) materiau
    read (1,*) control
    read (1,*) segments

    close (1)

    !---------------------------------------------------------

    ! Some material parameter must be loaded
    open(1,file="../in/"//materiau,STATUS='OLD')

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
    read(1,*) crystal_structure(1:3) ! the crystal symmetry
    read(1,*) useless
    read(1,*) VecBurgers    ! The burgers vector in Angstrom
    VecBurgers = VecBurgers * 1.d-10

    close (1)

    !---------------------------------------------------------

    ! Shall we work with 'Seg_save' or 'film.bin'
    write(*,*) "Do you want to analyse the -Seg_save- or -film.bin- file"
    write(*,*) "F = Seg_save, T = film.bin"
    read (*,*) Film

    if (Film) then

      print *, " > We are going to analyse dislocation configuration in -film.bin-"
      ! The file we want to analyze
      fichier = "film.bin"
      open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
      ! Important information at the begining of the film.bin file
      read(50) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
      read(50) nsegmaxf     ! The nsegmax value used in this simulation
      read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
      read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry

      allocate  (bvecnor(3,nbase))
      allocate  (bveclin(3,nbase))

      read(50) (bveclin(:,I),bvecnor(:,I), I=1,nbase)   ! table of vector direction and displacement for this simulation
      read(50) modur(:)                                 ! Dimension of the 3D simulated volume in a value

      ! Test on the nsegmax definition
      if (nsegmax /= nsegmaxf) then
        write(*,*) "Number of segments in film.bin = ",nsegmaxf,"but in this programm = ",nsegmax
        write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
        stop
      endif

    else

      print *, " > We are going to analyse dislocation configuration in -Seg_save-"

      ! The simulation elementary lattice parameter
      cristallo = '../out/BVD.'//crystal_structure(1:3)

      open(1,FILE=cristallo,STATUS='OLD')

      read (1,*) nbase , avalue

      allocate (bvecnor(3,nbase))
      allocate (bveclin(3,nbase))

      do j = 1,nbase
          read (1,*) i, bveclin(1:3,j)
          if(j /= i) print *, "stop, a problem appears when reading BVD"
      enddo

      close(1)

      open(50,file='../in/SEG_save',STATUS='OLD')

      read(50,*) useless
      read(50,*) Nsegm                       ! Nb of segments

      ! Test on the nsegmax definition
      if (nsegm > nsegmax) then
        write(*,*) "Number of segments in -Seg_save- = ",nsegmaxf,"but in this program = ",nsegmax
        write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
        stop
      endif

      read(50,*) modur(:)  ! Dimensions of the simulation cell

      allocate (seg(nsegm))
      allocate  (iseg(5,nsegm))
      allocate  (tabvois(2,nsegm))
      allocate  (junctt(nsegm))
      allocate  (ID_surface(nsegm))

      do I = 1,nsegm

        read (50,*) j,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
                   seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%Ijonc,seg(i)%gd

        print*,j,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
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

    !---------------------------------------------------------

    ! The normal of slip planes is defined
    allocate (normlin(nbase))
    normlin(:) = 0
    do j = 1,nbase
        normlin(j)=sqrt(real(bveclin(1,j))**2 + real(bveclin(2,j))**2 + real(bveclin(3,j))**2)
    enddo

    !---------------------------------------------------------

    ! The configuration(s) we want to analyse is defined

    If (Film) then
      ! We first select the part of the film.bin file we want to use for the stats
      write(*,*) "Analyse start at simulation step:  ? "
      read (*,*) debutaffich
      write(*,*) "Analyse end at simulation step:  ? "
      read (*,*) finaffich
      write(*,*) "Snapshot frequency for the film analysis: (Warning : must be coherent with KSTATS) ? "
      read (*,*) freqaffich
      if (freqaffich  <= 0)freqaffich = 1
    else !set 1 for all parameters in case of SEG_save
      debutaffich = 1
      finaffich   = 1
      freqaffich  = 1
    endif

    nb_step_stat = 0    ! The variable counting the number of configuration analysed is initialized

    ! The slip system symmetry is defined
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

    !---------------------------------------------------------

    ! Definition of the calculation lattice
    write(*,*) "We now defined the parameters used to defined tha calculation lattice!"
    write(*,*) " "
    write(*,*) " Please give the coordinates of the Lattice origin in a unit : Ox  Oy  Oz"
    read (*,*) Lattice_O(1:3)
    write(*,*) " >>>> ", Lattice_O(1:3)
    write(*,*) " "
    write(*,*) " Please give the 3 vectors defining the Lattice elementary cell (must be in a unit)!"
    write(*,*) " "
    write(*,*) " V1x  V1y  V1y ?"
    read (*,*) V1(1:3)
    write(*,*) " V2x  V2y  V2y ?"
    read (*,*) V2(1:3)
    write(*,*) " V3x  V3y  V3y ?"
    read (*,*) V3(1:3)
    write(*,*) " >>>> V1 = ", V1(1:3)
    write(*,*) " >>>> V2 = ", V2(1:3)
    write(*,*) " >>>> V3 = ", V3(1:3)
    write(*,*) " "
    write(*,*) " Please give the 3 integer number defining the lattice dimension : n1  n2  n3"
    read (*,*) n1, n2, n3
    write(*,*) " "

    !---------------------------------------------------------

    allocate (dislocation_density(1:n1,1:n2,1:n3,1:ntsg+1))
    allocate (junction_density   (1:n1,1:n2,1:n3,1:ntsg+1))

    ! The simulation box volume and lattice elementary cell
    A(1:3) = real(v1(1:3),dp)
    B(1:3) = real(v2(1:3),dp)
    C(1:3) = real(v3(1:3),dp)

    volb    = Determinant(A(:), B(:), C(:)) * avalue**3
    voltot  = volb * n1 * n2 * n3

    !---------------------------------------------------------

    ! Setup the result directory
    call date_and_time(date,time)!,zone,values)
    mydate=date//'_'//time(1:6)

    command = 'mkdir ../out/Density_dir'//'_'//mydate
    call system(command)

    ! WE START THE LOOP ON SIMULATION STEPS SAVED IN FILM.BIN

    print *, "-------------------------"
    print *, " Beginning of step loop "
    print *, "-------------------------"
    print *, ""

    ! INITIALIZATION ON THE LOOP ON THE SELECTED STEP
    kk=0

    Steps_loop_22: do while (kk >= 0 .and. kk <= finaffich)

      TheGoodStep = .false.

      if (Film) then

        ! A new simulation step
        read(50,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO

        TheGoodStep = (kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0)

        if (TheGoodStep) then

          print*," > Selected step "
          write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I5)') KK,EPSO,XSTRESS,nsegm

          allocate (seg(nsegm))
          allocate (iseg(5,nsegm))
          allocate (junctt(nsegm))
          allocate (tabvois(2,nsegm))
          allocate (ID_surface(nsegm))

          ! The segments micro structure is loaded
          read(50,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),ID_surface(j), J=1,nsegm)

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

        else

          ! We pass this step
          read(50,IOSTAT=KODERR)

        endif

        if(KODERR < 0) then
          print*,"  End of the film.bin file "
          exit Steps_loop_22
        elseif(KODERR > 0) then
          print*,"Reading problem 1",KODERR,kk
          stop
        endif

        ! The mark for the end of step is read
        read(50,IOSTAT=KODERR) caractest

        if(KODERR < 0) then
            print*,"  End of the film.bin file 3"
            exit Steps_loop_22
        elseif(KODERR > 0) then
            print*,"Reading error - Could not find the mark for the end of step ",caractest
            stop
        endif

      else

        kk = 1 ! Seg_save

      endif

    !---------------------------------------------------------

      ! The steps we want to use for statistic
      if (TheGoodStep) then ! Begin of the "Big If"

        ! Start of the calculations

!        ib(:) = 0
        dislocation_density(:,:,:,:) = 0
        junction_density(:,:,:,:) = 0

        !Loop on the segments
        do i = 1, nsegm

          ! Only segments with a non-zero length are of interest
          if (seg(i)%norme /= 0) then

            ! Identification of the segment slip system
            veclin = seg(i)%veclin
            indexB = (int((veclin-1)/8))*8 + 1
            tmp2 = indexB/8
            sys = floor(tmp2)+1

            Oa(1:3) = seg(i)%O(:)

            ! precise computation of the density imposes to decompose the length of segments
            do J = 1,seg(i)%norme

              ! Vector connecting the density latice origin and the points along the segment length
              ra(:)   = (real(Oa(:)) + (real(J)-0.5) * real(bveclin(:,veclin))) - real(Lattice_O(:))
              D(1:3)  = ra(1:3)

              ! Calculation of the coordinate of ra() in the density grid
              ! Using the Cramer's rule we solve a linear system :
              ! v1 * fac1 + v2 * fac2 + v3 * fac3 = ra
              sysdet = Determinant(A, B, C)

              if (Abs(sysdet) < 10e-15) Then
                coord = (/-300000,-300000,-300000/)
              else
                xdet = Determinant(D, B, C)
                ydet = Determinant(A, D, C)
                zdet = Determinant(A, B, D)

                coord(1)=xdet/sysdet
                coord(2)=ydet/sysdet
                coord(3)=zdet/sysdet
              endif

              ! coordinate of the segment section in the density tab
              c1 = nint(coord(1)) + 1
              c2 = nint(coord(2)) + 1
              c3 = nint(coord(3)) + 1

              ! We must test if the segment section is outside the grid
              if (c1 < 1 .or. c2 < 1 .or. c3 < 1 .or. c1 > n1 .or. c2 > n2 .or. c3 > n3) cycle

              ! density calculation
              ! junction case: divided by 2 because junctions are made of 2 segments

              if (seg(i)%jonc) then

                dislocation_density(c1,c2,c3,sys) =  dislocation_density(c1,c2,c3,sys) + 0.5 * NormLin(veclin)
                junction_density(c1,c2,c3,sys)    =  junction_density(c1,c2,c3,sys) + 0.5 * NormLin(veclin)
                !Total
                dislocation_density(c1,c2,c3,ntsg+1) =  dislocation_density(c1,c2,c3,ntsg+1) + 0.5 * NormLin(veclin)
                junction_density(c1,c2,c3,ntsg+1)    =  junction_density(c1,c2,c3,ntsg+1) + 0.5 * NormLin(veclin)

              else

                dislocation_density(c1,c2,c3,sys) =  dislocation_density(c1,c2,c3,sys) + NormLin(veclin)
                !Total
                dislocation_density(c1,c2,c3,ntsg+1) =  dislocation_density(c1,c2,c3,ntsg+1) + NormLin(veclin)

              endif

            enddo       ! loop on the segment length (j = 1,seg(i)%norme)

          endif         ! Non zero test

        enddo           ! loop on the segment number (i = 1, nsegm)

        dislocation_density(:,:,:,:)  = dislocation_density(:,:,:,:) * avalue / volb
        junction_density(:,:,:,:)     = junction_density(:,:,:,:) * avalue / volb

        ! Management of the calculation outputs
        nb_step_stat = nb_step_stat +1                    ! The filename number

    1151 format ('../out/Density_dir','_',A15,'/dislocations_dens_incr',I0, '.txt') ! The file name construction

        write (filename, 1151) mydate, nb_step_stat       ! definition of a different filename for each loop
        open(unit=212, file=filename, status='unknown')  ! The file name is open

        ! the total density results

        !Writing of results for this step
        ! column 1: indice of box along the v1 direction
        ! column 2: indice of box along the v2 direction
        ! column 3: indice of box along the v3 direction
        ! column 4->ntsg+4: dislocation/junction density on the box

        print*,"Writing of density in the output file out/Density_dir_xxx/"
        write(212,'(2x,A2,2x,A2,2x,A9,2x,A9,2x,A20)') 'X','Y','Z','dislocation densities'

        do c1 = 1, n1
          do c2 = 1,n2
            do c3 = 1,n3
            write(212,*) c1, c2, c3, dislocation_density(c1,c2,c3,1:ntsg+1)
            enddo
          enddo
        enddo

        close(212)

        ! the junction density results
    1152 format ('../out/Density_dir','_',A15,'/junctions_dens_incr',I0, '.txt')         ! The file name construction

        write (filename, 1152) mydate, nb_step_stat                ! definition of a different filename for each loop
        open(unit=212, file=filename, status='unknown')  ! The file name is open

        print*,"Writing of an approximate junction density in the output file out/Density_dir_xxx"
        write(212,'(2x,A2,2x,A2,2x,A9,2x,A9,2x,A32)') 'X','Y','Z','approximate junction densities'

        do c1 = 1, n1
          do c2 = 1,n2
            do c3 = 1,n3
            write(212,*) c1, c2, c3, junction_density(c1,c2,c3,1:ntsg+1)
            enddo
          enddo
        enddo

        close(212)


        ! END OF THE CALCULATIONS MADE ON THIS STEP
        !  ======>>>>>

        deallocate (seg)
        deallocate (iseg)
        deallocate (junctt)
        deallocate (tabvois)
        deallocate (ID_surface)

      endif ! End of the "big if" selecting the step

      ! In the Seg_Save case only one pass is made
      if (.not. Film) kk = kk + 1

    enddo Steps_loop_22

    print *, ""
    print *, "-------------------------"
    print *, " End of step loop "
    print *, "-------------------------"
    print *, ""

    print *,' > Results are in /out/density_dir'

    close(50)

    print*," > Density calculation achieved with success"

    deallocate (dislocation_density)
    deallocate (junction_density)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (10)                                                  !
  ! Mean shear strain evolution on a given axis direction             !
  ! Uses gammabox.bin                                                 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(10)

    ! 1) INPUTS
    ! We first select the part of the film.bin file we want to use for the stats
    write(*,*) "statistics start at step:  ?"
    read (*,*) debutaffich
    write(*,*) "statistics end at step:  ?"
    read (*,*) finaffich
    write(*,*) "Snapshot frequency:  ?"
    read (*,*) freqaffich
    if (freqaffich  <= 0)freqaffich = 1

    ! The file we want to analyze
    fichier = "gammabox.bin"
    open(61,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
    read(61) NBoites
    read(61) NboitesX
    read(61) NboitesY
    read(61) NboitesZ
    read(61) NTSG
    print *, 'NboitesX', NboitesX
    print *, 'NboitesY', NboitesY
    print *, 'NboitesZ', NboitesZ
    print *, 'NTSG', NTSG

    max_length = max(NboitesX,NBoitesY,NBoitesZ)

    allocate(somme(NTSG))
    allocate(Gammabox_Table(NBoites,NTSG))
    allocate(Gamma_res(2*max_length,NTSG))
    allocate(tens_box(3,3,NBoites))

    ! Initializations
    somme(:)                 = 0.
    Gammabox_Table(:,:)      = 0.
    tens_box(:,:,:)          = 0.

    ! 2) LECTURE OF FILM.BIN
    ! The results files
    open(62,FILE='../out/gammabox_result.txt',STATUS='UNKNOWN')
    open(63,FILE='../out/straintensor_result.txt',STATUS='UNKNOWN')

    command = 'mkdir ../out/GammaboxVTK'
    call system(command)

    ! READING OF GAMMABOX.BIN
    ! Initialization of parameters for the loop
    kk=1

    ! start of the loop
    ijk = -1

    Steps_loop_5: do while (kk.gt.0)

      ijk=ijk+1

      ! A new simulation step
      read(61,IOSTAT=KODERR) kk
      read(61,IOSTAT=KODERR) (somme(JK), JK=1,NTSG)

      ! The segments micro structure is loaded
      do ni=1,NBoites
        read(61,IOSTAT=KODERR) Px,Py,Pz,(Gammabox_Table(ni,JK),JK=1,NTSG)
      enddo

      ! Calculation of the strain tensor for all boxes combined
      do JK=1,NTSG
        do i=1,3
          do j=1,3
            tens_sum(i,j) = tens_sum(i,j) + (somme(JK)*tens(i,j,JK))
          enddo
        enddo
      enddo

      ! Calculation of the strain tensor for each box
      do ni =1,NBoites
        do i=1,3
          do j=1,3
            do JK=1,NTSG
              tens_box(i,j,ni) = tens_box(i,j,ni) + (Gammabox_Table(ni,JK) * tens(i,j,JK))
            enddo
          enddo
        enddo
      enddo

      !3) ANALYSIS ON THE DESIRED STEPS
      !The steps we want to use for statistic

      if(kk >= debutaffich .and. kk <=finaffich .and. modulo(ijk, freqaffich) == 0) then ! Begin of the "Big If"


        !VTK FILE FOR GAMMABOX SUM over steps
115     format ('../out/GammaboxVTK/GammaboxVTK_kk',I0, '.vtk')  ! The file name construction
        write (filename, 115) kk                                 ! definition of a different filename for each loop

        open(unit=10, status='unknown',file=filename)

        write(10,'(A)') '# vtk DataFile Version 2.0'
        write(10,'(A)') 'Unstructured Grid Field'
        write(10,'(A)') 'ASCII'
        write(10,'(A)') 'DATASET UNSTRUCTURED_GRID'

        ! The field coordinates
        write(10,'(A,I10,A)') 'POINTS', NBoitesX*NBoitesY*NBoitesZ, '  float'

        do IZ = 1,NBoitesZ
          do IY = 1,NBoitesY
            do IX = 1,NBoitesX

              !Formule which gives the box number function of IX,IY,IZ (similar to B1D_3D)
              !I = IX + (IY-1) * NBoitesX + (IZ-1)* (NBoitesX*NBoitesY)
              write(10,*) IX,IY,IZ
              !Gammabox_Table(I,JK)

            enddo
          enddo
        enddo

        write(10,'(A,I10)') 'POINT_DATA', NBoitesX*NBoitesY*NBoitesZ

        do JK = 1,NTSG

          write (filename, *) JK
          ! The field amplitude in real

          write(10,'(A)') 'SCALARS sys'//ADJUSTL(TRIM(filename))//' float'
          write(10,'(A)') 'LOOKUP_TABLE default'

          do IZ = 1,NBoitesZ
            do IY = 1,NBoitesY
              do IX = 1,NBoitesX

                !Formule which gives the box number function of IX,IY,IZ (similar to B1D_3D)
                I = IX + (IY-1) * NBoitesX + (IZ-1)* (NBoitesX*NBoitesY)
                write(10,*) Gammabox_Table(I,JK)

              enddo
            enddo
          enddo
        enddo

        close(10)

        ! The output file is written
        write(*,*) "Writting of results on file gammabox_result.txt | kk =",kk
        write(62,'(a,I6)') 'kk=',kk
        write(62,'(12(1x,a,I2))') ('Gam on sys:',JK,JK=1,NTSG)
        write(62,'(12(2x,e11.4))') (somme(JK), JK=1,NTSG)
        write(62,*)

        write(62,'(12(4x,a,I2))') ('Gam box:',JK, JK=1,NTSG)
        do I=1,NBoites
          write(62,'(12(6x,e11.4))') (Gammabox_Table(I,JK), JK=1,NTSG)
        enddo

        ! The tensor output file is written
        write(*,*) "Writting of results on file straintensor_result | kk =",kk
        write(63,'(a,I6)') 'kk=',kk
        write(63,'(a,I6)') ('Strain Tensor on simulation box (e11, e12, e13, e21, e22, etc.:')
        do i=1,3
            write(63,'(6(6x,e11.4))') (tens_sum(i,j), j=1,3)
        enddo

        write(63,'(a,I6)') ('Strain Tensor on sub boxes:')

        do ni=1,NBoites
            write(63,*) &
          ni,tens_box(1,1,ni),tens_box(2,2,ni),tens_box(3,3,ni),tens_box(1,2,ni),tens_box(2,3,ni),tens_box(3,1,ni)
        enddo

      endif ! End of the "big if" selecting the step

      ! 4) END
      read(61,IOSTAT=KODERR) caractest
      if(KODERR /= 0) then
        print*,"  End of the gammabox.bin file"
        stop
      endif

    enddo Steps_loop_5

    deallocate(somme)
    deallocate(Gammabox_Table)
    deallocate(Gamma_res)

    close(61)
    close(62)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCULATION (11)                                                  !
  ! 3D -> 2D extraction from film                                     !
  ! 2D coordinate of segments cutting a reference plane are extracted !
  ! from a 3D microstructure for 2D analysis.                         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(11)

    ! We first setup the result directory
    command = 'rm -fr ../out/2D_extract'
    call system(command)
    command = 'mkdir ../out/2D_extract'
    call system(command)
    nb_step_stat = 0

    ! We first select the part of the film.bin file we want to use for the stats
    write(*,*) "statistics start at simulation step:  ?"
    read (*,*) debutaffich
    write(*,*) "statistics end at simulation step:  ?"
    read (*,*) finaffich
    write(*,*) "Snapshot frequency for the film analysis:  ?"
    read (*,*) freqaffich
    if (freqaffich  <= 0) freqaffich = 1

    ! Definition of the reference plane for the 2D extraction
    write(*,*) "Enter the 3 Miller indices m(1:3) of the reference plane for the 2D extraction"
    read (*,*) Miller(1:3)
    write(*,*) "Enter a reel value defining the plane location p=(m(1)x+m(2)y+m(3)z)."
    read (*,*) Pos
    write(*,*) "Enter the 3 Miller indices of one reference direction in the extraction plane"
    read (*,*) V2Dx(1:3)

    ! Calculation of the two reference directions in the extraction plane V2Dx and V2Dy
    normivect=dsqrt(real(Miller(1)*Miller(1)+Miller(2)*Miller(2)+Miller(3)*Miller(3),DP))
    Mil(1:3) = Miller(1:3) / normivect
    normivect=dsqrt(V2Dx(1)*V2Dx(1)+V2Dx(2)*V2Dx(2)+V2Dx(3)*V2Dx(3))
    V2Dx(1:3) = V2Dx(1:3) / normivect
    V2Dy(1) = Mil(2)*V2Dx(3)-Mil(3)*V2Dx(2)
    V2Dy(2) = Mil(3)*V2Dx(1)-Mil(1)*V2Dx(3)
    V2Dy(3) = Mil(1)*V2Dx(2)-Mil(2)*V2Dx(1)

    ! Loop to find the selected step
    kk=1

    ! The file we want to analyze
    fichier = "film.bin"
    open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)

    ! Lecture of informations contained on the file, not used for the calculation
    read(50) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
    read(50) nsegmaxf     ! The nsegmax value used in this simulation
    read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
    read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry
    read(50) (bveclin(:,I),bvecnor(:,I), I=1,nbase)    ! table of vector direction and displacement for this simulation
    read(50) modur(:)     ! Dimension of the 3D simulated volume in a value

    allocate(iseg(5,nsegmax))
    allocate(tabvois(2,nsegmax))
    allocate(junctt(nsegmax))
    allocate(ID_surface(nsegmax))

    ! We start the loop at the simulation steps saved in film.bin
    ijk = -1

    Steps_loop_6: do while (kk.gt.0)

      ijk = ijk + 1

      ! A new simulation step
      read(50,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO
      write(*,'("kk = ",I7)') KK

      if(KODERR < 0) then
        print*,"  End of the film.bin file"
        stop
      elseif(KODERR > 0) then
        print*,"Reading problem 4",KODERR,ijk
        kk=1
        cycle Steps_loop_6
      endif

      read(50,IOSTAT=KODERR) (Iseg(1:5,k),Junctt(k),tabvois(1:2,k),ID_surface(k), k=1,nsegm)

        if(KODERR < 0) then
          print*,"  End of the film.bin file"
          stop
        elseif(KODERR > 0) then
          print*,"Reading Error 2",KODERR,nsegm,j
          cycle Steps_loop_6
        endif

      ! Extraction at the selected steps
      if(kk >= debutaffich .and. kk <=finaffich .and. modulo(ijk, freqaffich) == 0) then ! Begin of the "Big If"

        !print*,'kk=',kk

        deallocate (seg)
        allocate (seg(nsegm))

        ! The film.bin information is used to overwrite the SEG_save information
        Do j=1,nsegm

          seg(j)%O(1:3) = Iseg(1:3,j)
          seg(j)%norme  = Iseg(4,j)
          seg(j)%veclin = Iseg(5,j)

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

        ! Allocation of the useful list
        allocate (Sys2D(nsegm))
        allocate (Sign2D(nsegm))
        allocate (rxy(2,nsegm))
        allocate (Jonc2D(nsegm))

        Jonc2D(:) = 0    ! The default solution = segment is not forming a junction
        j = 0            ! Number of segments cutting the plane

        do i=1, nsegm

          if (seg(i)%norme /= 0) then

            Oi2D(:) = real(seg(i)%O(:))
            Ei2D(:) = real(seg(i)%O(:) + seg(i)%norme*bveclin(:,seg(i)%veclin))

            VecOi(:) = Oi2D(:) - Pos*Mil(:)
            VecEi(:) = Ei2D(:) - Pos*Mil(:)
            ScaO = (Mil(1)*VecOi(1) + Mil(2)*VecOi(2) + Mil(3)*VecOi(3))
            ScaE = (Mil(1)*VecEi(1) + Mil(2)*VecEi(2) + Mil(3)*VecEi(3))

            ! Test to keep only segments with O and E on each side of the extraction plane
            if (ScaO*ScaE > 0.) cycle

            j = j + 1     ! one more segment

            !Determination of the segment Sys2D and sign
            tmp       = (seg(i)%veclin-1)/8 +1
            Sys2D(j)  = floor(tmp)
            btmp2D(:) = bveclin(:,seg(i)%veclin)
            Sign2D(j) = sign(1,int(Miller(1)*btmp2D(1)+Miller(2)*btmp2D(2)+Miller(3)*btmp2D(3)))

            ! Junction case
            if (seg(i)%jonc) Jonc2D(j) = seg(i)%ijonc

            ! The intersection coordinates is calculated with the function InterPLanSeg
            Point(:) = InterPlanSeg(Miller,Pos,Oi2D,btmp2D)

            ! The point coordinates are resolved in the 2D plane
            rxy(1,j) = (Point(1)*V2Dx(1)+Point(2)*V2Dx(2)+Point(3)*V2Dx(3))
            rxy(2,j) = (Point(1)*V2Dy(1)+Point(2)*V2Dy(2)+Point(3)*V2Dy(3))

          endif

        enddo

        print*,'out loop',KK,'with 2D_nsegm =',j

        ! management of the calculation outputs
        nb_step_stat = nb_step_stat +1                    ! The filename number
        write (filename, 112) nb_step_stat                ! definition of a different filename for each loop
112     format ('../out/2D_extract/2DSeg', I0, '.txt')    ! The file name construction
        open(unit=212, file=filename, status='unknown')  ! The file name is open

        !Writing of results for this step
        print*,"Writing of 2D microstructure in the output file out/2D_extract/2DSegXXXX.txt"
        write(212,*) 'kk=',kk

        do i = 1,j
          write(212,*) i,rxy(1,i),rxy(2,i),Sys2D(i),Sign2D(i),Jonc2D(i)
        enddo

        close(212)

        ! END OF THE CALCULATIONS MADE ON THIS STEP
        !  ======>>>>>

        deallocate (Sys2D)
        deallocate (Sign2D)
        deallocate (rxy)
        deallocate (Jonc2D)

      endif ! End of the "big if" selecting the step

      ! The mark for the end of step is read
      read(50,IOSTAT=KODERR) caractest

      if(KODERR < 0) then
          print*,"  End of the film.bin file"
          stop
      elseif(KODERR > 0) then
          print*,"Reading error 3",caractest
          cycle Steps_loop_6
      endif

    enddo Steps_loop_6

    print *,'results are in /out/2D_extract'
    read *

    close(50)


  case default
  print*,"bye bye ..."

  end select

enddo  ! End of the main do while loop

stop


contains

!###########################################################################
!#    Calcul du rayon de courbure pour un cercle passant par trois points  #
!# M,P et Q                                                                #
!#    Pour un rayon de courbure infinie le vecteur retourne par la         #
!# fonction est nul.                                                       #
!#    Alpha : angle entre M et Q sur le cercle                             #
!#    T : vecteur tangant ie normal a R dans MPQ                           #
!#                                                                         #
!############################################################### 30/07/99 ##
FUNCTION VECCOURB(OM,OP,OQ,TESTPMPQ,alpha)

implicit none

real(kind=DP), dimension(3)            :: VECCOURB,resul
Integer(kind=DPI), dimension(3),intent(in) :: OM,OP,OQ
logical, intent(out)                  :: TESTPMPQ
real(kind=DP)                          :: ALPHA
real(kind=DP), dimension(3)            :: PM,PQ
real(kind=DP)                          :: PMPQ,PMPQ2,PM2,PQ2,PM2PQ2
real(kind=DP)                          :: FAC1,FAC2,FAC3
real(kind=DP), parameter               :: HALF=0.5D0, un= 1.0D0

PM(1:3)=real(OM(1:3)-OP(1:3),DP)
PQ(1:3)=real(OQ(1:3)-OP(1:3),DP)
PMPQ=PM(1)*PQ(1)+PM(2)*PQ(2)+PM(3)*PQ(3)
PMPQ2=PMPQ*PMPQ
PM2=PM(1)*PM(1)+PM(2)*PM(2)+PM(3)*PM(3)
PQ2=PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3)
PM2PQ2=PM2*PQ2
FAC1=(PM2PQ2-PMPQ2)


IF(FAC1.NE.0.0D0) THEN
   FAC1=HALF/FAC1
   FAC2=PQ2*(PM2-PMPQ)
   FAC3=PM2*(PQ2-PMPQ)
   RESUL(1:3)=FAC1*(FAC2*PM(1:3)+FAC3*PQ(1:3))
   alpha = cose2v(-resul(1:3)+PM(1:3),-resul(1:3)+PQ(1:3))
   if (dabs(alpha).gt. real(1,DP)) alpha = sign(real(1,DP),alpha)
   alpha=acos(alpha)
   VECCOURB(1:3)=RESUL(1:3)
ELSE
   VECCOURB(1:3)=0.0D0
   ALPHA=0
ENDIF

if (PMPQ.gt.0) then
   TESTPMPQ=.true.
   ALPHA=4.0_DP*atan(un)   ! Pi
   VECCOURB(1:3)=PM(1:3)+PQ(1:3) ! ici seul la direction de veccourb importe
else
   TESTPMPQ=.false.
endif

END FUNCTION VECCOURB

!###########################################################################
!# Calcul de la norme d'un vecteur entier                                  #
!############################################################### 06/11/98 ##
function norivect(i)

implicit none

integer(DPI)   :: i(3)
real(kind=DP)    :: r(3),norivect

r(1:3)=real((i(1:3)),DP)
norivect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

end function norivect

!###########################################################################
!# Cosinus entre deux vecteurs                                             #
!############################################################### 06/11/98 ##
function cose2v(a,b)   ! cos2v(a,b) serait sans doute mieux

implicit none

real(DP) :: cose2v,a(3),b(3),na,nb

na=(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))**0.5
nb=(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))**0.5
if (na.ne.(0.).and.nb.ne.(0.)) then
   cose2v = dot_product(a,b)/(na*nb)
else
   cose2v = 0.
endif

end function cose2v

!###########################################################################
!# Calcul du vecteur normalise associe a un vecteur entier                 #
!############################################################### 06/11/98 ##
function normaivec(i)

implicit none

integer(DPI) :: i(3)
real(kind=DP) :: normaivec(3),r(3),deno

deno = UN/norivect(i)
r=real(i,DP)
normaivec(:)=r(:)*deno

end function normaivec

!###########################################################################
!# Pour extraire d'un tenseur, l'amplitude d'un champ dans une direction   #
!# donnee. Et vis versa.                                                   #
!# SENS=.True.   /  Amplitude uniaxiale + axe -> tenseur.                  #
!# SENS=.False.  /  Tenseur + axe             -> amplitude uniaxiale       #
!###########################################################################
subroutine orientdieg(ampli,z,xtens,sens)

implicit none

logical        :: sens
integer        :: i
real(kind=8)   :: ampli,xtens(3,3)
real(kind=8)   :: z(3)

if(sens) then

   xtens(1:3,1:3) = 0.
   do i=1,3
       xtens(i,1) = z(i)*z(1)*ampli
       xtens(i,2) = z(i)*z(2)*ampli
       xtens(i,3) = z(i)*z(3)*ampli
   enddo

else

   ampli = 0.
   do i=1,3
       ampli=ampli+xtens(i,1)*z(i)*z(1)
       ampli=ampli+xtens(i,2)*z(i)*z(2)
       ampli=ampli+xtens(i,3)*z(i)*z(3)
   enddo

endif

end subroutine orientdieg

!##########################################################################
!# Calculation of the internal stress associate to a dislocation          #
!# microstructure with nsegm segments at point r                          #
!##########################################################################
subroutine sigma_int

!$ use OMP_LIB

IMPLICIT NONE

INTEGER,PARAMETER  :: DP=selected_real_kind(p=14)   !< Definition of real number kind
integer,parameter  :: DPI=selected_int_kind(9)      !< Definition of integer number kind

integer(kind=dpi),allocatable  :: origin(:,:)        !< coordonnees de l'origine des segments
integer(kind=dpi),allocatable    :: veclin(:)      !< indice du vecteur ligne et deplacement
integer(kind=dpi),allocatable    :: norme(:)       !< longueur du segment (positif)
logical,allocatable ::  jonc(:)                   !< clef d'immobilisation pour les jonctions

integer(kind=dpi)   :: icristallo   !<
integer, dimension(3) :: timer
integer, dimension(3) :: montemps


integer(kind=dpi)   :: nb_step_stat !<
integer(kind=dpi)   :: I            !<
integer(kind=dpi)   :: J            !<
integer(kind=dpi)   :: indexB       !<
real(kind=dp)       :: useless      !< Useless real variable
integer(kind=dpi)   :: uselessi     !< Useless integer variable

integer(kind=DPI) :: nbase        !<
integer(kind=DPI) :: k            !<
integer(kind=DPI) :: kk           !<
integer(kind=DPI) :: kkk          !<
integer(kind=DPI) :: ns          !<
integer(kind=DPI) :: debutaffich  !<
integer(kind=DPI) :: freqaffich   !<
integer(kind=DPI) :: Nnodes(2)       !<
integer(kind=DPI) :: KODERR       !<
integer(kind=DPI) :: nsegm        !< Number of segments
integer(kind=DPI) :: nsegmaxf     !<
integer(kind=DPI) :: finaffich    !<
integer(kind=DPI) :: NR           !< Number of line in input file
integer(kind=DPI) :: nline        !< Line counter
integer(kind=dpi)   :: IndexSeg

integer(kind=DPI) :: modur1     !<
integer(kind=DPI) :: modur2     !<
integer(kind=DPI) :: modur3     !<

integer(kind=DPI) :: BoxCenter1 !<
integer(kind=DPI) :: BoxCenter2 !<
integer(kind=DPI) :: BoxCenter3 !<

integer(kind=DPI) :: NRep1 !<
integer(kind=DPI) :: NRep2 !<
integer(kind=DPI) :: NRep3 !<

integer(kind=DPI),allocatable     :: junctt(:)             !< The junction pair-segment
integer(kind=DPI), allocatable     :: Veclin_T(:)            !< the segments veclin_number for selected segments
integer(kind=DPI), allocatable     :: Norme_T(:)             !< the segments norme for selected segments
integer(kind=DPI), allocatable     :: Origin_T(:,:)           !< the segments origin for selected segments

integer(kind=DPI), allocatable :: bveclin(:,:)       !<

real(kind=dp)  :: GridNormalNorm !<
real(kind=dp)  :: GridDist!<
real(kind=dp)  :: widthp(3)      !< The 3 grid step
real(kind=dp)  :: coeftaille     !<

real(kind=dp), allocatable    :: normlin(:) !<

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
real(kind=dp)  :: r1            !< Coordinate of the stress calculation
real(kind=dp)  :: r2            !< Coordinate of the stress calculation
real(kind=dp)  :: r3            !< Coordinate of the stress calculation
real(kind=dp),dimension(6)  :: dombounds    !< Boundaries if restriction on the exploration domain
real(kind=DP),Dimension(3)  :: MapCenter    !<
real(kind=DP),Dimension(3)  :: maxbounds     !< Maximum value in X, Y and Z direction
real(kind=DP),Dimension(3)  :: minbounds     !< Minimum value in X, Y and Z direction

real(kind=DP) :: sigtmp1,sigtmp2,sigtmp3,sigtmp4,sigtmp5,sigtmp6                        !<
real(kind=DP) :: sigtmprep1,sigtmprep2,sigtmprep3,sigtmprep4,sigtmprep5,sigtmprep6     !<

real(kind=DP), allocatable   :: sigint1(:),sigint2(:),sigint3(:),sigint4(:),sigint5(:),sigint6(:)   !<
real(kind=DP), allocatable   :: sigintrep1(:),sigintrep2(:),sigintrep3(:),sigintrep4(:),sigintrep5(:),sigintrep6(:)   !<
real(kind=DP), allocatable   :: mycoords(:,:)      !<  Store the X,Y,Z coordinates from a file ( mydata('Number of lines',3) )

integer                     :: Film         !< True if Film.bin is used
logical                     :: sqgrid       !< True if we want to draw a square grid (same grid step in X and Y direction)
logical                     :: IsSingular   !< True if we want the microMegas singular field | False allows to choose core radius spreading factor
logical                     :: Islimits     !< True if we want limit the exploration to a subdomain
logical                     :: StationFile  !< True if we want to use a file containing X,Y,Z Coordinates
logical,dimension(2)         :: conditions   !< Condition for segments filtering

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
character(len=15)   :: argc       !< Number of threads in OMP
character(len=15)   :: core_radius_str !< String to store core radius for write statements

integer(kind=dpi)           :: ii
integer(kind=dpi)           :: caraci
integer(kind=dpi)           :: veclini
integer(kind=dpi)           :: normei
integer(kind=dpi)           :: NX,NY,NZ
integer(kind=dpi)           :: Tr1,Tr2,Tr3
integer(kind=dpi)           :: ab1,ab2,ab3
integer(kind=dpi)           :: centrei1,centrei2,centrei3

real(kind=dp) :: T1i,T2i,T3i
real(kind=dp) :: B1i,B2i,B3i
real(kind=dp) :: BPVT1i,BPVT2i,BPVT3i
real(kind=dp) :: ra01,ra02,ra03
real(kind=dp) :: rb01,rb02,rb03
real(kind=dp) :: ra_1,ra_2,ra_3,rb_1,rb_2,rb_3
real(kind=DP) :: d2,s1,s2,a2_p_d2,inv_a2_p_d2
real(kind=DP) :: s1Ra3inv,Ra1_3inv,Ra1_inv
real(kind=DP) :: s2Ra3inv,Ra2_3inv,Ra2_inv
real(kind=DP) :: s1_03,s1_13,s1_05,s1_15,s1_25
real(kind=DP) :: s2_03,s2_13,s2_05,s2_15,s2_25
real(kind=DP) :: s_03,s_13,s_05,s_15,s_25
real(kind=DP) :: m4pn,mn4pn,a2m4pn,a2m8p
real(kind=DP) :: d_pv_b_ps_t,commun,R_ps_t,Ra2
real(kind=DP) :: cai_d1,cai_d2,cai_d3
real(kind=DP) :: d_pv_b1,d_pv_b2,d_pv_b3
real(kind=DP) :: t_pv_b1,t_pv_b2,t_pv_b3
real(kind=DP) :: I_25_1,I_25_2,I_25_3,I_25_4,I_25_5,I_25_6
real(kind=DP) :: I_15_1,I_15_2,I_15_3,I_15_4,I_15_5,I_15_6
real(kind=DP) :: I_05_1,I_05_2,I_05_3,I_05_4,I_05_5,I_05_6
real(kind=DP) :: I_13_1,I_13_2,I_13_3,I_13_4,I_13_5,I_13_6
real(kind=DP) :: I_03_1,I_03_2,I_03_3,I_03_4,I_03_5,I_03_6
real(kind=DP) :: t_pt_d_pv_b1,t_pt_d_pv_b2,t_pt_d_pv_b3,t_pt_d_pv_b4,t_pt_d_pv_b5,t_pt_d_pv_b6
real(kind=DP) :: d_pt_t_pv_b1,d_pt_t_pv_b2,d_pt_t_pv_b3,d_pt_t_pv_b4,d_pt_t_pv_b5,d_pt_t_pv_b6
real(kind=DP) :: t_pt_t_pv_b1,t_pt_t_pv_b2,t_pt_t_pv_b3,t_pt_t_pv_b4,t_pt_t_pv_b5,t_pt_t_pv_b6
real(kind=DP) :: t_pt_d1,t_pt_d2,t_pt_d3,t_pt_d4,t_pt_d5,t_pt_d6
real(kind=DP) :: t_pt_t1,t_pt_t2,t_pt_t3,t_pt_t4,t_pt_t5,t_pt_t6
real(kind=DP) :: d_pt_d1,d_pt_d2,d_pt_d3,d_pt_d4,d_pt_d5,d_pt_d6

real(kind=DP) :: bdivpa       !<
real(kind=DP) :: bdiva        !<
real(kind=dp) :: core2        !< (Core radius / avalue)^2

! cutoff variables
real(kind=DP), dimension(11)        :: buf
real(kind=DP)                       :: realnorme,dt2,dn2,dt
real(kind=DP)                       :: dist2,ExtR2
logical                             :: ied

integer :: threadnum !< Thread number if openMP multithreading


CALL getarg(1, argc)

If (argc=='-smp') then
  CALL getarg(2, argc)
  read( argc, '(i10)' ) threadnum
  IF (threadnum > 0) then
    write(*,*) "-------------------------------------------- "
    write(*,*) "Number of Threads asked = ", argc
    write(*,*) "-------------------------------------------- "
write(*,*) "                                             "
  ELSE
    print *, "Bad number of threads = ", argc
  ENDIF
else
  threadnum=4
endif

if (threadnum==4) then
  write(*,*) "-------------------------------------------- "
  write(*,*) " Threads number is set by default to 4       "
  write(*,*) " To change the number of threads please use :"
  write(*,*) "                                             "
  write(*,*) "       ./histo -smp  nbthreads               "
  write(*,*) "                                             "
  write(*,*) "---------------------------------------------"
  write(*,*) "                                             "
endif

    ! The simulation parameters

  open(1,file='../in/input.dd',STATUS='OLD')

  carac = " "
  do while (carac /= "-")
    read (1,*) carac
  enddo

  read (1,*) materiau
  read (1,*) control
  read (1,*) segments

  close (1)

  ! Identification of the crystal symetry
  open(1,file="../in/"//materiau,STATUS='OLD')

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
  film = -1
  do while (film < 1 .or. film > 2)
    write(*,*) "Use Film.bin (1) or use SEGSAVE (2) ? "
    read (*,*) film
  enddo

  if (Film== 1) then

    print *, " > We use film.bin informations "
    ! The file we want to analyze
    fichier = "film.bin"
    open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
    ! Important information at the begining of the film.bin file
    read(50) icristallo   ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
    read(50) nsegmaxf     ! The nsegmax value used in this simulation
    read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
    read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry

    allocate (bveclin(3,nbase))
    allocate (origin(3,nsegmax))
    allocate (Jonc(nsegmax))
    allocate (Junctt(nsegmax))
    allocate (Norme(nsegmax))
    allocate (veclin(nsegmax))

    read(50) (bveclin(:,I),uselessi,uselessi,uselessi, I=1,nbase)    ! table of vector direction and displacement for this simulation
    read(50) modur1,modur2,modur3     ! Dimension of the 3D simulated volume in a value

    ! Test on the nsegmax definition
    if (nsegmax /= nsegmaxf) then
      write(*,*) "Number of segments in film.bin = ",nsegmaxf,"but in this programm = ",nsegmax
      write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
      stop
    endif

  else

    print *, " > We use SEGSAVE informations "

    ! The simulation elementary lattice parameter
    cristallo = '../out/BVD.'//crystal_structure(1:3)

    open(1,FILE=cristallo,STATUS='OLD')

    read (1,*) nbase , avalue
    allocate (bveclin(3,nbase))

    do j = 1,nbase
        read (1,*) i, bveclin(1:3,j)
        if(j /= i) print *, "stop, a problem appears when reading BVD"
    enddo

    close(1)

    open(50,file='../in/SEG_save',STATUS='OLD')

    read(50,*) useless                     ! Parameter of effective loading on the systems
    read(50,*) Nsegm                       ! Nb of segments

    allocate (origin(3,Nsegm))
    allocate (Jonc(Nsegm))
    allocate (Norme(Nsegm))
    allocate (veclin(Nsegm))

    ! Test on the nsegmax definition
    if (nsegm > nsegmax) then
      write(*,*) "Number of segments in SEGSAVE = ",nsegmaxf,"but in this program = ",nsegmax
      write(*,*) 'Please built program with the correct value of nsegmax : STOP !'
      stop
    endif

    read(50,*) modur1,modur2,modur3  ! Dimensions of the simulation cell

    do I = 1,nsegm

      read (50,*) j,origin(:,i),norme(i),veclin(i),uselessi,uselessi, &
                 uselessi,uselessi,JONC(i),uselessi,uselessi

    enddo

  endif

  ! Bdiva
  BdivA = vecburgers / Avalue

  ! The center of the simulation cell
  BoxCenter1 = modur1/2
  BoxCenter2 = modur2/2
  BoxCenter3 = modur3/2

  allocate (normlin(nbase))

  do j = 1,nbase
      normlin(j)=sqrt(real(bveclin(1,j))**2+  &
                      real(bveclin(2,j))**2+  &
                      real(bveclin(3,j))**2)
  enddo

  nb_step_stat = 0

  If (Film == 1) then
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
    write(*,*) "How nany nodes per side shall we take for the grid "
    read(*,*) Nnodes(1)

    DO WHILE (Nnodes(1) < 10 .or. Nnodes(2) > 500 )
      write(*,*) "<! Please consider a number of nodes between 10 and 500 !> "
      write(*,*) "How nany nodes per side shall we take for the grid "
      read(*,*) Nnodes(1)
    ENDDO

    NNodes(2)=Nnodes(1)

    write(*,*) "Aspect ratio between stress mapping and the simulation cell "
    read(*,*) coeftaille

    write(*,*) "Do you want a squared grid (T) or exact grid (F)"
    read(*,*) sqgrid

    minbounds(:) = 0.

    maxbounds(1)=modur1
    maxbounds(2)=modur2
    maxbounds(3)=modur3

  else

    filepath='../in/mM_Data_Files/locations_dat'
    KODERR=0
    NR=-1
    OPEN(90,file=filepath,STATUS='OLD')
      read(90,*,IOSTAT=KODERR)
      DO WHILE (KODERR==0)
        NR=NR+1
        read(90,*,IOSTAT=KODERR)
      ENDDO
    CLOSE(90)

    print *, "Number of line in the file locations_dat : ", NR

    OPEN(90,file=filepath,STATUS='OLD')
    !Now we can allocate data variables
    ALLOCATE(mycoords(NR,3))

    !Now read data into mydata
    READ(90,*)
    DO J=1,NR
      READ(90,*,IOSTAT=KODERR) mycoords(J,1), mycoords(J,2), mycoords(J,3)
    ENDDO
    CLOSE(90)
    Nnodes(1) = 1
    Nnodes(2) = NR
    maxbounds(:)= MAXVAL (mycoords, DIM=1)
    minbounds(:)= MINVAL (mycoords, DIM=1)

  endif

  write(*,*) "Number of replicates per side of the reference volume Nx, Ny, Nz ? "
  read(*,*) NRep1,Nrep2,Nrep3

  if (.not. StationFile) then

    write(*,*) "Please enter the 2D grid normal direction, x, y, z "
    read(*,*) GridNormal(1:3)

    GridNormalNorm = sqrt(GridNormal(1)**2+GridNormal(2)**2+GridNormal(3)**2)
    GridNormal(1:3) = GridNormal(1:3)/GridNormalNorm

    write(*,*) "Enter the normal distance of the 2D grid from the simulation cell center (microns) "
    read(*,*) GridDist
    GridDist = GridDist * 1.D-6 / avalue

  endif !if not StationFile

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
    core_radius = 0.5*BdivA !This Value should be the same as in microMegas - BdivA for simulation value
  else
    write(*,*) "Which core radius spreading factor for short range stress calculation ('halfthickness' in simulation unit) ? "
    read(*,*) core_radius !Halfthickness already in simulation unit
  endif

  write(*,fmt='("Do you want to explore only a subdomain (T/F) ?  [Default : ",6I9," ] ")') &
  & int(minbounds(1)),int(maxbounds(1)),int(minbounds(2)),int(maxbounds(2)),int(minbounds(3)),int(maxbounds(3))
  read(*,*) Islimits

  If (Islimits) then

    dombounds(:)=-100.

    DO WHILE (dombounds(1) < 0 .or. dombounds(1) > modur1)
      write (*,*) ' xmin  >'
      read  (*,*) dombounds(1)
    END DO

    DO WHILE (dombounds(2) < 0 .or. dombounds(2) > modur1 .or. dombounds(2) < dombounds(1))
      write (*,*) ' xmax >'
      read  (*,*) dombounds(2)
    END DO

    DO WHILE (dombounds(3) < 0 .or. dombounds(3) > modur2)
      write (*,*) ' ymin >'
      read  (*,*) dombounds(3)
    END DO

    DO WHILE (dombounds(4) < 0 .or. dombounds(4) > modur2 .or. dombounds(4) < dombounds(3))
      write (*,*) ' ymax >'
      read  (*,*) dombounds(4)
    END DO

    DO WHILE (dombounds(5) < 0 .or. dombounds(5) > modur3)
      write (*,*) ' zmin >'
      read  (*,*) dombounds(5)
    END DO

    DO WHILE (dombounds(6) < 0 .or. dombounds(6) > modur3 .or. dombounds(6) < dombounds(5))
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

  if (Film == 1) then
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

  if (StationFile) then
    nline = nr
  else
    nline = NNodes(1)*Nnodes(2)
    allocate(mycoords(nline,3))
  endif

  allocate (sigint1(nline))
  allocate (sigint2(nline))
  allocate (sigint3(nline))
  allocate (sigint4(nline))
  allocate (sigint5(nline))
  allocate (sigint6(nline))

  allocate (sigintrep1(nline))
  allocate (sigintrep2(nline))
  allocate (sigintrep3(nline))
  allocate (sigintrep4(nline))
  allocate (sigintrep5(nline))
  allocate (sigintrep6(nline))

  ! Setup the result directory
  call date_and_time(date,time)!,zone,values)
  mydate=date//'_'//time(1:6)

  command = 'mkdir ../out/SigInt_dir' !//'_'//mydate
  call system(command)

  ! WE START THE LOOP ON SIMULATION STEPS SAVED IN FILM.BIN

  print *, "-------------------------"
  print *, " Beginning of step loop "
  print *, "-------------------------"
  print *, ""

  ! INITIALIZATION ON THE LOOP ON THE SELECTED STEP
  kk=0
  call itime(timer)
  montemps = timer

  if (.not. StationFile) then
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
    if (sqgrid) then  ! the square grid solution
      widthp(:) = coeftaille * max(Modur1,Modur2,Modur3) / real(Nnodes(1))
    else              ! the exact shape grid
      widthp(1) = coeftaille * Modur1 / real(Nnodes(1))
      widthp(2) = coeftaille * Modur2 / real(Nnodes(1))
      widthp(3) = coeftaille * Modur3 / real(Nnodes(1))
    endif

    kkk = 0

    MapCenter(1) = BoxCenter1+GridNormal(1)*GridDist
    MapCenter(2) = BoxCenter2+GridNormal(2)*GridDist
    MapCenter(3) = BoxCenter3+GridNormal(3)*GridDist


    do i=1,Nnodes(1)

      rr1 = mapcenter(1) + (i-real(Nnodes(1),DP)*0.5 - 0.5)*X(1)*widthp(1)
      rr2 = mapcenter(2) + (i-real(Nnodes(1),DP)*0.5 - 0.5)*X(2)*widthp(2)
      rr3 = mapcenter(3) + (i-real(Nnodes(1),DP)*0.5 - 0.5)*X(3)*widthp(3)

      do j=1,Nnodes(2)
        kkk = kkk + 1
        mycoords(kkk,1) = rr1 + (j-real(Nnodes(2),DP)*0.5 - 0.5)*Y(1)*widthp(1)
        mycoords(kkk,2) = rr2 + (j-real(Nnodes(2),DP)*0.5 - 0.5)*Y(2)*widthp(2)
        mycoords(kkk,3) = rr3 + (j-real(Nnodes(2),DP)*0.5 - 0.5)*Y(3)*widthp(3)

      enddo

    enddo

  endif

  !Parameters for the calculation of sigma
  core2 = (core_radius)*(core_radius)

  bdivpa = bdiva/pI

  m4pn=0.25d0/(1.0d0-dpoiss)
  mn4pn=m4pn*dpoiss
  a2m4pn=core2*m4pn
  a2m8p=core2*0.125d0
  !####################################



  Steps_loop_4: do while (kk.ge.0 .and. kk .le. finaffich)

    sigint1(:)     = 0. !Initialization
    sigint2(:)     = 0.
    sigint3(:)     = 0.
    sigint4(:)     = 0.
    sigint5(:)     = 0.
    sigint6(:)     = 0.

    sigintrep1(:)  =  0.
    sigintrep2(:)  =  0.
    sigintrep3(:)  =  0.
    sigintrep4(:)  =  0.
    sigintrep5(:)  =  0.
    sigintrep6(:)  =  0.

    if (Film == 1) then

      ! A new simulation step
      read(50,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO
      if (modulo(kk, 1000) == 0) print *, kk
      if(kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0) then
        print*," > Selected step "
        write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I5)') KK,EPSO,XSTRESS, nsegm
      endif

      if(KODERR < 0) then
        print*,"  End of the film.bin file "
        exit Steps_loop_4
      elseif(KODERR > 0) then
        print*,"Reading problem 5",KODERR,kk
        stop
        !cycle Steps_loop_4
      endif


      if(kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0) then

        read (50,IOSTAT=KODERR) (origin(1:3,j),norme(j),veclin(j),Junctt(j), uselessi,uselessi,uselessi, J=1,nsegm)
      !READ(50,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),ID_surface(j), J=1,nsegm)

        if(KODERR < 0) then
          print*,"  End of the film.bin file"
          exit Steps_loop_4
        elseif(KODERR > 0) then
          print*,"Reading Error 2",KODERR,nsegm,j
          stop
          !cycle Steps_loop_4
        endif
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

        cycle Steps_loop_4

      endif

    else

      kk=1 !Segsave

    endif


    ! The steps we want to use for statistic
    if(kk >= debutaffich .and. kk <=finaffich .and. modulo(kk, freqaffich) == 0) then ! Begin of the "Big If"

      if (Film == 1) then

      print*," > film.bin informations loaded !"
      ! Start of the calculations
      endif
      !Initialisation
      ns = 0

      !First pass to know segment of interrest number
      do K = 1,nsegm

        conditions(:)=.False.

        ! Not a kneecap
        conditions(1) = (norme(k) /= 0)

        !In the slipsystem(s) requested
        conditions(2)= .TRUE.
        if (slipsys < 13) then
           IndexSeg = (int((veclin(k)-1)/8)) + 1
           if (Junctt(k) > 0) IndexSeg = 0
           if (IndexSeg /= slipsys) conditions(2)= .FALSE.
        endif

        if (All(conditions)) then
          ns=ns+1
        endif

      enddo

      !---------------------
      !Allocation
      allocate(veclin_T(ns))
      allocate(Norme_T(ns))
      allocate(Origin_T(3,ns))

      !Initialsation
      veclin_T(:) = 0
      Norme_T(:) = 0
      Origin_T(:,:) = 0

      !Reinitialisation of ns
      ns=0
      !---------------------

      !Store only segments of interrest - Reduction of the arrays
      do K = 1,nsegm

        ! Not a kneecap
        conditions(1) = (norme(k) /= 0)

        !In the slipsystem(s) requested
        conditions(2)= .TRUE.
        if (slipsys < 13) then
           IndexSeg = (int((veclin(k)-1)/8)) + 1
           if (junctt(k)>0) IndexSeg = 0
           if (IndexSeg /= slipsys) conditions(2)= .FALSE.
        endif

        if (All(conditions)) then
          ns=ns+1
          !attributions
          Veclin_T(ns)  = veclin(k)
          Norme_T(ns)   = norme(k)
          Origin_T(:,ns)  = origin(:,k)

        endif

      enddo


      !---------------------------------------------------------------
      ! Multi threading
      !---------------------------------------------------------------
      !$OMP PARALLEL NUM_THREADS(threadnum) DEFAULT(PRIVATE)  , &
      !$OMP& SHARED(normlin)                        , &
      !$OMP& SHARED(BoxCenter1,BoxCenter2,BoxCenter3)   , &
      !$OMP& SHARED(core2,core_radius)                              , &
      !$OMP& SHARED(m4pn,mn4pn,a2m4pn,a2m8p)            , &
      !$OMP& SHARED(bdivpa)                             , &
      !$OMP& SHARED(ns)                                 , &
      !$OMP& SHARED(nline)                                 , &
      !$OMP& SHARED(origin_t)      , &
      !$OMP& SHARED(modur1,modur2,modur3)               , &
      !$OMP& SHARED(dombounds)                          , &
      !$OMP& SHARED(veclin_t,norme_t)                   , &
      !$OMP& SHARED(bveclin)                            , &
      !$OMP& SHARED(Nrep1,Nrep2,Nrep3)                  , &
      !$OMP& SHARED(mycoords)                  , &
      !$OMP& SHARED(sigint1) , &
      !$OMP& SHARED(sigint2)   , &
      !$OMP& SHARED(sigint3)   , &
      !$OMP& SHARED(sigint4)   , &
      !$OMP& SHARED(sigint5)   , &
      !$OMP& SHARED(sigint6) ,&
      !$OMP& SHARED(sigintrep1)        , &
      !$OMP& SHARED(sigintrep2)    , &
      !$OMP& SHARED(sigintrep3)   , &
      !$OMP& SHARED(sigintrep4)   , &
      !$OMP& SHARED(sigintrep5)   , &
      !$OMP& SHARED(sigintrep6)
       !
      !$OMP DO
      do kkk=1,nline

        r1 = mycoords(kkk,1)
        r2 = mycoords(kkk,2)
        r3 = mycoords(kkk,3)

        if ((r1 >= dombounds(1)) .and. (r2 >= dombounds(3)) .and. (r3 >= dombounds(5)) .and.  &
          & (r1 <= dombounds(2)) .and. (r2 <= dombounds(4)) .and. (r3 <= dombounds(6))) then

          sigtmp1 =  0.0d0
          sigtmp2 =  0.0d0
          sigtmp3 =  0.0d0
          sigtmp4 =  0.0d0
          sigtmp5 =  0.0d0
          sigtmp6 =  0.0d0

          sigtmprep1 =  0.0d0
          sigtmprep2 =  0.0d0
          sigtmprep3 =  0.0d0
          sigtmprep4 =  0.0d0
          sigtmprep5 =  0.0d0
          sigtmprep6 =  0.0d0


          ! The translation needed for periodic boundaries
          ! With this translation r is shifted at the center of the simulation box
          Tr1   = BoxCenter1- nint(r1)
          Tr2   = BoxCenter2- nint(r2)
          Tr3   = BoxCenter3- nint(r3)

          do ii = 1,ns

            do NX = -Nrep1,Nrep1
            do NY = -Nrep2,Nrep2
            do NZ = -Nrep3,Nrep3

              veclini = Veclin_T(ii)
              normei  = Norme_T(ii)

              ! The segment vector direction
              T1i = bveclin(1,veclini)/normlin(veclini)
              T2i = bveclin(2,veclini)/normlin(veclini)
              T3i = bveclin(3,veclini)/normlin(veclini)

              ! The Burgers vector of i
              indexB = (int((veclini-1)/8))*8 + 1
              B1i = bveclin(1,indexB)/normlin(indexB)
              B2i = bveclin(2,indexB)/normlin(indexB)
              B3i = bveclin(3,indexB)/normlin(indexB)

              !B produit vectoriel T
              BPVT1i = B2i*T3i-B3i*T2i
              BPVT2i = B3i*T1i-B1i*T3i
              BPVT3i = B1i*T2i-B2i*T1i

              ! information relative to the segments
              ! The segment type (the tyseg type 1 <-> 4!)
              caraci = int(veclini - (((veclini-1)/(nbasered/2))*(nbasered/2)))
              ! The vetors connecting r to O and E of segment i

              ab1 = normei*bveclin(1,veclini)
              ab2 = normei*bveclin(2,veclini)
              ab3 = normei*bveclin(3,veclini)

              centrei1 = Origin_T(1,ii) + ab1/2
              centrei2 = Origin_T(2,ii) + ab2/2
              centrei3 = Origin_T(3,ii) + ab3/2

              centrei1 = modulo(centrei1+Tr1, ModuR1)
              centrei2 = modulo(centrei2+Tr2, ModuR2)
              centrei3 = modulo(centrei3+Tr3, ModuR3)

              ra01 = real(BoxCenter1 - (centrei1 - ab1/2),DP)
              ra02 = real(BoxCenter2 - (centrei2 - ab2/2),DP)
              ra03 = real(BoxCenter3 - (centrei3 - ab3/2),DP)

              rb01 = real(BoxCenter1 - (centrei1 + ab1/2),DP)
              rb02 = real(BoxCenter2 - (centrei2 + ab2/2),DP)
              rb03 = real(BoxCenter3 - (centrei3 + ab3/2),DP)

              centrei1 = centrei1 + NX*MODUR1
              centrei2 = centrei2 + NY*MODUR2
              centrei3 = centrei3 + NZ*MODUR3

              buf(1)  =real(BoxCenter1 - centrei1,DP) ! v
              buf(2)  =real(BoxCenter2 - centrei2,DP)
              buf(3)  =real(BoxCenter3 - centrei3,DP)

              realnorme = real(normei,DP)
              dist2=buf(1)*buf(1)+buf(2)*buf(2)+buf(3)*buf(3)

              ! if the point is at a distance smaller then cutoff radius from the center of ns segment
              ! no calculation of the stress field

              if (dist2 < 4*core2) cycle

              ExtR2=core_radius + realnorme*0.5d0
              ExtR2=ExtR2*ExtR2
              !if the distance between segment centers <= than the cut-off radius we already know we must compute the short-range correction.
              !On the other hand, if the distance > we need a test to verify if the center of the segment (R) in which we want to compute the force is
              !inside the homogeneization region of the segment i
              if (dist2 < ExtR2) then
                buf(4)    =buf(1)*t1i+buf(2)*t2i+buf(3)*t3i! v.t
                buf(5)    =buf(4)/(t1i*t1i+t2i*t2i+t3i*t3i)! (v.t)/(t.t)
                buf(6)  =buf(5)*t1i ! vt=(v.t)/(t.t)*t
                buf(7)  =buf(5)*t2i
                buf(8)  =buf(5)*t3i
                buf(9:11) =buf(1:3)-buf(6:8) ! vn=v-vt
                dn2=buf(9)*buf(9)+buf(10)*buf(10)+buf(11)*buf(11)
                dt2=buf(6)*buf(6)+buf(7)*buf(7)+buf(8)*buf(8)
                dt=dsqrt(dt2)
                ied=((dt-realnorme*0.5d0)**2+dn2 < 4*core2) ! test if the point is inside a sphere at the extremity
                if ((ied .or. (dn2 < 4*core2 .and. dt < realnorme*0.5d0))) cycle! R is inside the homogeneizetion area of segment I !
              endif

              ra_1 = ra01 + NX*MODUR1
              ra_2 = ra02 + NY*MODUR2
              ra_3 = ra03 + NZ*MODUR3
              rb_1 = rb01 + NX*MODUR1
              rb_2 = rb02 + NY*MODUR2
              rb_3 = rb03 + NZ*MODUR3

              R_ps_t=ra_1*T1i+ra_2*T2i+ra_3*T3i
              s1=-R_ps_t

              cai_d1=ra_1-R_ps_t*T1i
              cai_d2=ra_2-R_ps_t*T2i
              cai_d3=ra_3-R_ps_t*T3i

              d2=cai_d1*cai_d1+cai_d2*cai_d2+cai_d3*cai_d3 ! distance au segment
              R_ps_t=rb_1*T1i+rb_2*T2i+rb_3*T3i
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

              t_pv_b1=-BPVT1i
              t_pv_b2=-BPVT2i
              t_pv_b3=-BPVT3i


              d_pv_b1=cai_d2*B3i-cai_d3*B2i
              d_pv_b2=cai_d3*B1i-cai_d1*B3i
              d_pv_b3=cai_d1*B2i-cai_d2*B1i

              d_pv_b_ps_t=d_pv_b1*T1i+d_pv_b2*T2i+d_pv_b3*T3i

              d_pt_d1 = cai_d1*cai_d1
              d_pt_d2 = cai_d2*cai_d2
              d_pt_d3 = cai_d3*cai_d3
              d_pt_d4 = cai_d2*cai_d3
              d_pt_d5 = cai_d1*cai_d3
              d_pt_d6 = cai_d1*cai_d2

              t_pt_t1 = T1i*T1i
              t_pt_t2 = T2i*T2i
              t_pt_t3 = T3i*T3i
              t_pt_t4 = T2i*T3i
              t_pt_t5 = T1i*T3i
              t_pt_t6 = T1i*T2i

              t_pt_d1 = (T1i*cai_d1)*2.d0
              t_pt_d2 = (T2i*cai_d2)*2.d0
              t_pt_d3 = (T3i*cai_d3)*2.d0
              t_pt_d4 = T2i*cai_d3+cai_d2*T3i
              t_pt_d5 = T1i*cai_d3+cai_d1*T3i
              t_pt_d6 = T1i*cai_d2+cai_d1*T2i

              t_pt_t_pv_b1 = (T1i*t_pv_b1)*2.d0
              t_pt_t_pv_b2 = (T2i*t_pv_b2)*2.d0
              t_pt_t_pv_b3 = (T3i*t_pv_b3)*2.d0
              t_pt_t_pv_b4 = T2i*t_pv_b3+t_pv_b2*T3i
              t_pt_t_pv_b5 = T1i*t_pv_b3+t_pv_b1*T3i
              t_pt_t_pv_b6 = T1i*t_pv_b2+t_pv_b1*T2i

              d_pt_t_pv_b1 = (cai_d1*t_pv_b1)*2.d0
              d_pt_t_pv_b2 = (cai_d2*t_pv_b2)*2.d0
              d_pt_t_pv_b3 = (cai_d3*t_pv_b3)*2.d0
              d_pt_t_pv_b4 = cai_d2*t_pv_b3+t_pv_b2*cai_d3
              d_pt_t_pv_b5 = cai_d1*t_pv_b3+t_pv_b1*cai_d3
              d_pt_t_pv_b6 = cai_d1*t_pv_b2+t_pv_b1*cai_d2

              t_pt_d_pv_b1 = (T1i*d_pv_b1)*2.d0
              t_pt_d_pv_b2 = (T2i*d_pv_b2)*2.d0
              t_pt_d_pv_b3 = (T3i*d_pv_b3)*2.d0
              t_pt_d_pv_b4 = T2i*d_pv_b3+d_pv_b2*T3i
              t_pt_d_pv_b5 = T1i*d_pv_b3+d_pv_b1*T3i
              t_pt_d_pv_b6 = T1i*d_pv_b2+d_pv_b1*T2i

              commun=m4pn*d_pv_b_ps_t

              I_03_1=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b1)-0.25d0*t_pt_d_pv_b1
              I_03_2=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b2)-0.25d0*t_pt_d_pv_b2
              I_03_3=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b3)-0.25d0*t_pt_d_pv_b3
              I_03_4=m4pn*d_pt_t_pv_b4-0.25d0*t_pt_d_pv_b4
              I_03_5=m4pn*d_pt_t_pv_b5-0.25d0*t_pt_d_pv_b5
              I_03_6=m4pn*d_pt_t_pv_b6-0.25d0*t_pt_d_pv_b6


              if (caraci==1) then

                if (nx/=0.or.ny/=0.or.nz/=0)then !If any Replica

                  sigtmprep1=sigtmprep1+ I_03_1*s_03*bdivpa
                  sigtmprep2=sigtmprep2+ I_03_2*s_03*bdivpa
                  sigtmprep3=sigtmprep3+ I_03_3*s_03*bdivpa
                  sigtmprep4=sigtmprep4+ I_03_4*s_03*bdivpa
                  sigtmprep5=sigtmprep5+ I_03_5*s_03*bdivpa
                  sigtmprep6=sigtmprep6+ I_03_6*s_03*bdivpa

                else

                  sigtmp1=sigtmp1+ I_03_1*s_03*bdivpa
                  sigtmp2=sigtmp2+ I_03_2*s_03*bdivpa
                  sigtmp3=sigtmp3+ I_03_3*s_03*bdivpa
                  sigtmp4=sigtmp4+ I_03_4*s_03*bdivpa
                  sigtmp5=sigtmp5+ I_03_5*s_03*bdivpa
                  sigtmp6=sigtmp6+ I_03_6*s_03*bdivpa

                  sigtmprep1=sigtmprep1+ I_03_1*s_03*bdivpa
                  sigtmprep2=sigtmprep2+ I_03_2*s_03*bdivpa
                  sigtmprep3=sigtmprep3+ I_03_3*s_03*bdivpa
                  sigtmprep4=sigtmprep4+ I_03_4*s_03*bdivpa
                  sigtmprep5=sigtmprep5+ I_03_5*s_03*bdivpa
                  sigtmprep6=sigtmprep6+ I_03_6*s_03*bdivpa

                endif

              else

                if (nx/=0.or.ny/=0.or.nz/=0)then !If any Replica

                  I_13_1=-mn4pn*t_pt_t_pv_b1
                  I_13_2=-mn4pn*t_pt_t_pv_b2
                  I_13_3=-mn4pn*t_pt_t_pv_b3
                  I_13_4=-mn4pn*t_pt_t_pv_b4
                  I_13_5=-mn4pn*t_pt_t_pv_b5
                  I_13_6=-mn4pn*t_pt_t_pv_b6

                  I_15_1=a2m8p*t_pt_t_pv_b1-commun*t_pt_d1
                  I_15_2=a2m8p*t_pt_t_pv_b2-commun*t_pt_d2
                  I_15_3=a2m8p*t_pt_t_pv_b3-commun*t_pt_d3
                  I_15_4=a2m8p*t_pt_t_pv_b4-commun*t_pt_d4
                  I_15_5=a2m8p*t_pt_t_pv_b5-commun*t_pt_d5
                  I_15_6=a2m8p*t_pt_t_pv_b6-commun*t_pt_d6

                  I_05_1=commun*(core2+d_pt_d1)-a2m8p*t_pt_d_pv_b1
                  I_05_2=commun*(core2+d_pt_d2)-a2m8p*t_pt_d_pv_b2
                  I_05_3=commun*(core2+d_pt_d3)-a2m8p*t_pt_d_pv_b3
                  I_05_4=commun*d_pt_d4-a2m8p*t_pt_d_pv_b4
                  I_05_5=commun*d_pt_d5-a2m8p*t_pt_d_pv_b5
                  I_05_6=commun*d_pt_d6-a2m8p*t_pt_d_pv_b6

                  I_25_1=commun*t_pt_t1
                  I_25_2=commun*t_pt_t2
                  I_25_3=commun*t_pt_t3
                  I_25_4=commun*t_pt_t4
                  I_25_5=commun*t_pt_t5
                  I_25_6=commun*t_pt_t6

                  sigtmprep1=sigtmprep1+ (I_03_1*s_03 +I_13_1*s_13+I_05_1*s_05+I_15_1*s_15+I_25_1*s_25)*bdivpa
                  sigtmprep2=sigtmprep2+ (I_03_2*s_03 +I_13_2*s_13+I_05_2*s_05+I_15_2*s_15+I_25_2*s_25)*bdivpa
                  sigtmprep3=sigtmprep3+ (I_03_3*s_03 +I_13_3*s_13+I_05_3*s_05+I_15_3*s_15+I_25_3*s_25)*bdivpa
                  sigtmprep4=sigtmprep4+ (I_03_4*s_03 +I_13_4*s_13+I_05_4*s_05+I_15_4*s_15+I_25_4*s_25)*bdivpa
                  sigtmprep5=sigtmprep5+ (I_03_5*s_03 +I_13_5*s_13+I_05_5*s_05+I_15_5*s_15+I_25_5*s_25)*bdivpa
                  sigtmprep6=sigtmprep6+ (I_03_6*s_03 +I_13_6*s_13+I_05_6*s_05+I_15_6*s_15+I_25_6*s_25)*bdivpa

                else

                  I_13_1=-mn4pn*t_pt_t_pv_b1
                  I_13_2=-mn4pn*t_pt_t_pv_b2
                  I_13_3=-mn4pn*t_pt_t_pv_b3
                  I_13_4=-mn4pn*t_pt_t_pv_b4
                  I_13_5=-mn4pn*t_pt_t_pv_b5
                  I_13_6=-mn4pn*t_pt_t_pv_b6

                  I_15_1=a2m8p*t_pt_t_pv_b1-commun*t_pt_d1
                  I_15_2=a2m8p*t_pt_t_pv_b2-commun*t_pt_d2
                  I_15_3=a2m8p*t_pt_t_pv_b3-commun*t_pt_d3
                  I_15_4=a2m8p*t_pt_t_pv_b4-commun*t_pt_d4
                  I_15_5=a2m8p*t_pt_t_pv_b5-commun*t_pt_d5
                  I_15_6=a2m8p*t_pt_t_pv_b6-commun*t_pt_d6

                  I_05_1=commun*(core2+d_pt_d1)-a2m8p*t_pt_d_pv_b1
                  I_05_2=commun*(core2+d_pt_d2)-a2m8p*t_pt_d_pv_b2
                  I_05_3=commun*(core2+d_pt_d3)-a2m8p*t_pt_d_pv_b3
                  I_05_4=commun*d_pt_d4-a2m8p*t_pt_d_pv_b4
                  I_05_5=commun*d_pt_d5-a2m8p*t_pt_d_pv_b5
                  I_05_6=commun*d_pt_d6-a2m8p*t_pt_d_pv_b6

                  I_25_1=commun*t_pt_t1
                  I_25_2=commun*t_pt_t2
                  I_25_3=commun*t_pt_t3
                  I_25_4=commun*t_pt_t4
                  I_25_5=commun*t_pt_t5
                  I_25_6=commun*t_pt_t6

                  sigtmp1=sigtmp1+ (I_03_1*s_03+I_13_1*s_13+I_05_1*s_05+I_15_1*s_15+I_25_1*s_25)*bdivpa
                  sigtmp2=sigtmp2+ (I_03_2*s_03+I_13_2*s_13+I_05_2*s_05+I_15_2*s_15+I_25_2*s_25)*bdivpa
                  sigtmp3=sigtmp3+ (I_03_3*s_03+I_13_3*s_13+I_05_3*s_05+I_15_3*s_15+I_25_3*s_25)*bdivpa
                  sigtmp4=sigtmp4+ (I_03_4*s_03+I_13_4*s_13+I_05_4*s_05+I_15_4*s_15+I_25_4*s_25)*bdivpa
                  sigtmp5=sigtmp5+ (I_03_5*s_03+I_13_5*s_13+I_05_5*s_05+I_15_5*s_15+I_25_5*s_25)*bdivpa
                  sigtmp6=sigtmp6+ (I_03_6*s_03+I_13_6*s_13+I_05_6*s_05+I_15_6*s_15+I_25_6*s_25)*bdivpa

                  sigtmprep1=sigtmprep1+ (I_03_1*s_03 +I_13_1*s_13+I_05_1*s_05+I_15_1*s_15+I_25_1*s_25)*bdivpa
                  sigtmprep2=sigtmprep2+ (I_03_2*s_03 +I_13_2*s_13+I_05_2*s_05+I_15_2*s_15+I_25_2*s_25)*bdivpa
                  sigtmprep3=sigtmprep3+ (I_03_3*s_03 +I_13_3*s_13+I_05_3*s_05+I_15_3*s_15+I_25_3*s_25)*bdivpa
                  sigtmprep4=sigtmprep4+ (I_03_4*s_03 +I_13_4*s_13+I_05_4*s_05+I_15_4*s_15+I_25_4*s_25)*bdivpa
                  sigtmprep5=sigtmprep5+ (I_03_5*s_03 +I_13_5*s_13+I_05_5*s_05+I_15_5*s_15+I_25_5*s_25)*bdivpa
                  sigtmprep6=sigtmprep6+ (I_03_6*s_03 +I_13_6*s_13+I_05_6*s_05+I_15_6*s_15+I_25_6*s_25)*bdivpa

                endif

              endif

            enddo !NVOL
            enddo !NVOL
            enddo !NVOL

          enddo !i = 1,nsegm

          sigint1(kkk) = sigtmp1
          sigint2(kkk) = sigtmp2
          sigint3(kkk) = sigtmp3
          sigint4(kkk) = sigtmp4
          sigint5(kkk) = sigtmp5
          sigint6(kkk) = sigtmp6

          sigintrep1(kkk) =  sigtmprep1
          sigintrep2(kkk) =  sigtmprep2
          sigintrep3(kkk) =  sigtmprep3
          sigintrep4(kkk) =  sigtmprep4
          sigintrep5(kkk) =  sigtmprep5
          sigintrep6(kkk) =  sigtmprep6


        else

          sigint1(kkk)     = 0. !Nothing
          sigint2(kkk)     = 0.
          sigint3(kkk)     = 0.
          sigint4(kkk)     = 0.
          sigint5(kkk)     = 0.
          sigint6(kkk)     = 0.

          sigintrep1(kkk)  =  0.
          sigintrep2(kkk)  =  0.
          sigintrep3(kkk)  =  0.
          sigintrep4(kkk)  =  0.
          sigintrep5(kkk)  =  0.
          sigintrep6(kkk)  =  0.

        endif

      enddo   ! Nline

      !$OMP END DO
      !$OMP END PARALLEL

      !#########

      ! Management of the calculation outputs
      nb_step_stat = nb_step_stat +1                    ! The filename number

    OPEN(7,STATUS = 'SCRATCH') ! Get core radius as a string
    WRITE(7,*) core_radius
    REWIND(7)
    READ(7,*) core_radius_str
    CLOSE(7)

  114 format ('../out/SigInt_dir/sigint_h',A6,'_incr',I0, '.txt')         ! The file name construction

      write (filename, 114) core_radius_str, nb_step_stat                ! definition of a different filename for each loop
      open(unit=601, file=filename, status='unknown')  ! The file name is open

      if (StationFile) then
        write(601,'(2x,A2,2x,A2,2x,A2,6(2x,A5))') 'X','Y','Z','Sig11','Sig22','Sig33',&
           'Sig23','Sig13','Sig12'

        do j=1,nline
           write(601,*) mycoords(j,1:3), sigint1(j) * mu * 1D-6 ,    &
                                         sigint2(j) * mu * 1D-6 ,    &
                                         sigint3(j) * mu * 1D-6 ,    &
                                         sigint4(j) * mu * 1D-6 ,    &
                                         sigint5(j) * mu * 1D-6 ,    &
                                         sigint6(j) * mu * 1D-6
        enddo
      else
        write(601,'(2x,A2,2x,A2,6(2x,A5))') 'X','Y','Sig11','Sig22','Sig33',&
          'Sig23','Sig13','Sig12'
          kkk=0
          do i = 1, Nnodes(1)
            do j = 1, Nnodes(2)
              kkk= kkk+1
              write(601,*) i, j, sigint1(kkk) * mu * 1D-6 ,    &
                                 sigint2(kkk) * mu * 1D-6 ,    &
                                 sigint3(kkk) * mu * 1D-6 ,    &
                                 sigint4(kkk) * mu * 1D-6 ,    &
                                 sigint5(kkk) * mu * 1D-6 ,    &
                                 sigint6(kkk) * mu * 1D-6
            enddo
          enddo

     endif


      close(601)
      if (Nrep1 /= 0 .or. NRep2 /= 0 .or. Nrep3 /=0) then  !Stress files including replicas are wiritten only in the case that replicas exists!!!

  115 format ('../out/SigInt_dir/sigintrep_h',A6,'_incr',I0, '.txt')         ! The file name construction

        write (filename, 115) core_radius_str, nb_step_stat                ! definition of a different filename for each loop
        open(unit=602, file=filename, status='unknown')  ! The file name is open


        if(StationFile) then

          write(602,'(2x,A2,2x,A2,2x,A2,6(2x,A5))') 'X','Y','Z','Sig11','Sig22','Sig33',&
         'Sig23','Sig13','Sig12'

          do j=1,nline
            write(602,*) mycoords(j,1:3), sigintrep1(j) * mu * 1D-6 ,    &
                                          sigintrep2(j) * mu * 1D-6 ,    &
                                          sigintrep3(j) * mu * 1D-6 ,    &
                                          sigintrep4(j) * mu * 1D-6 ,    &
                                          sigintrep5(j) * mu * 1D-6 ,    &
                                          sigintrep6(j) * mu * 1D-6
          enddo

        else
          write(601,'(2x,A2,2x,A2,6(2x,A5))') 'X','Y','Sig11','Sig22','Sig33',&
            'Sig23','Sig13','Sig12'
            kkk=0
            do i = 1, Nnodes(1)
              do j = 1, Nnodes(2)
                kkk= kkk+1
                write(601,*) i, j, sigintrep1(kkk) * mu * 1D-6 ,    &
                                   sigintrep2(kkk) * mu * 1D-6 ,    &
                                   sigintrep3(kkk) * mu * 1D-6 ,    &
                                   sigintrep4(kkk) * mu * 1D-6 ,    &
                                   sigintrep5(kkk) * mu * 1D-6 ,    &
                                   sigintrep6(kkk) * mu * 1D-6
              enddo
            enddo

        endif

        close(602)

      endif !replicas
      ! END OF THE CALCULATIONS MADE ON THIS STEP
      !  ======>>>>>

      !de-allocation of temporary arrays to store info about selected segments
      deallocate(veclin_t)
      deallocate(Norme_T)
      deallocate(Origin_T)

    endif ! End of the "big if" selecting the step

    If (Film == 1) then
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

  deallocate(mycoords)

  deallocate (bveclin)
  deallocate (origin)
  deallocate (Jonc)
  deallocate (Norme)
  deallocate (veclin)
  deallocate (normlin)

  deallocate (sigint1)
  deallocate (sigint2)
  deallocate (sigint3)
  deallocate (sigint4)
  deallocate (sigint5)
  deallocate (sigint6)

  deallocate (sigintrep1)
  deallocate (sigintrep2)
  deallocate (sigintrep3)
  deallocate (sigintrep4)
  deallocate (sigintrep5)
  deallocate (sigintrep6)

  call itime(timer)
  write(*,*) "Execution time : ", (timer(1)*3600 + timer(2)*60 + timer(3)) &
             &-(montemps(1)*3600 + montemps(2)*60 + montemps(3)), "seconds"

end subroutine sigma_int

!################################################################################
!> \brief Calculation of the intersection between a plane and a line direction. #
!>     Miller  = Miller indices of the plane                                    #
!>     Pos     = Index of the plane                                             #
!>     OrigSeg = One point of the line                                          #
!>     VecSeg  = The line vector                                                #
!################################################################################
function InterPlanSeg (Miller,Pos,OrigSeg,VecSeg)

implicit none

real(kind=DP),intent(in)                  :: pos    !< Plane position (d)
real(kind=DP),dimension(3),intent(in)     :: OrigSeg
Integer(kind=DPI),dimension(3),intent(in) :: VecSeg
Integer(kind=DPI),dimension(3),intent(in) :: Miller !< Plane normal vector (a,b,c)

real(kind=DP),dimension(3)                :: InterPlanSeg
real(kind=DP),dimension(3)                :: VecSegreal
real(kind=DP)                             :: Interm
real(kind=DP),dimension(3)                :: CoefPlan
real(kind=DP)                             :: Coeff
real(kind=DP)                             :: scalaire


CoefPlan(1:3) = real(Miller(1:3),DP)
VecSegReal(1:3) = real(VecSeg(1:3),DP)

scalaire =   CoefPlan(1)*VecSegReal(1) &
           + CoefPlan(2)*VecSegReal(2) &
           + CoefPlan(3)*VecSegReal(3)


if (abs(scalaire) < 1.d-12) then
  ! The line is parallel to the plane
  InterplanSeg(:)=(/0.,0.,0./)                  ! Il n y a pas d intersection, donc zero
else

  Interm = + OrigSeg(1)* Coefplan(1)   &
           + OrigSeg(2)* Coefplan(2)   &
           + OrigSeg(3)* Coefplan(3)   &
           - Pos

  Coeff = -Interm/scalaire

  ! Coordinates of the line and plane intersection point
  InterPlanSeg(:) = OrigSeg(1:3) + Coeff*VecSegReal(:)
endif

end function InterPlanSeg

Function Determinant(A, B, C)

  Real(kind=DP), Dimension (3), Intent (in) :: A, B, C
  Real(kind=DP) :: Determinant

  Determinant = - A(3)*B(2)*C(1) + A(2)*B(3)*C(1) + A(3)*B(1)*C(2) &
                - A(1)*B(3)*C(2) - A(2)*B(1)*C(3) + A(1)*B(2)*C(3)

End Function Determinant

end program
