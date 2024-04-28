
!===================================================================================================
!========================    DEBUT   MODULE  "INIT"   ==============================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains many procedures needed to initialize (prepare) different steps of calculations
!> of the simulation. This include simulation data loading.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module INIT

use VARGLOB
use BRICAMAT
use DEBUG
use cartographie

#ifdef MDC
use mdctools
#endif

implicit none


real(DP)             :: echelle,Ldis_act,Ldis_nact,dmusurdT,moduleG0
real(DP)             :: facteurcoin(10),facteurvis(10),facteurmixt(10),rapport
real(DP)             :: det       ! determinant of the rotation matrix defined in matrotation_def
integer(DPI)         :: Nloi(NTSG_max,nbasered/2), facteur_boite_rot,Xlomax_nact
character(len=40)    :: fichier,fichbase,fichtab,materiau,control
character(len=34)    :: segments
character(len=40)    :: file_mat_rot
character(len=40)    :: fichseg
Integer (Kind=DPI)   :: shift_rotation


!#######
contains
!#######

!###########################################################################
! This subroutine load all the input parameters of the simulation
!
! Simulation Data are shared in between 4 different files:
!  1)  File 'materiaux' : The material parameters
!  2)  File 'control'   : Parameters used to control the simulaton
!  3)  File 'seg3D'     : The inital segment microstrucure

!###########################################################################
subroutine lire_donnees

implicit none

!# Declarations

real(kind=DP)       :: valeur,Burgerstemp
Integer (Kind=DPI)  :: i,itemp,ktemp,syst
Integer (Kind=DPI)  :: iphase
character(len=1)    :: carac
logical             :: test_PA

!Crack problem variables
character(len=256)  :: cmd !> String to run shell commands
character(len=256)  :: fact_type_str,Plate_a_str,Plate_H_str,Plate_W_str
Integer (Kind=DPI)  :: fact_type  ! Geometrical factor method
Integer (Kind=DPI)  :: Plate_H    ! Plate HALF height
Integer (Kind=DPI)  :: Plate_W    ! Plate weigth

#ifdef PA
Integer             :: Clef_PairImpair
character(len=40)   :: file_debug         !< Name of the debug file
character(len=2)    :: file_index         !< index of the file_debug in parallel computation
#endif

#ifdef MDC

! Variable declaration for MPI
integer                   :: ZebCheckTid=0
character(len=1024)       :: CurrentWorkingDirectory

#endif

! 99 format(I6,I6,I6)  !n de format mis arbitriarement sylvain

! Initialisation des loi de vitesse brute
Loi(1:NLV_MAX)%Arrhenius  = izero
Loi(1:NLV_MAX)%h          = zero
Loi(1:NLV_MAX)%V0         = zero
Loi(1:NLV_MAX)% deltaG0   = zero
Loi(1:NLV_MAX)% tau0      = zero
Loi(1:NLV_MAX)% coef_p    = zero
Loi(1:NLV_MAX)% coef_q    = zero
Loi(1:NLV_MAX)% BF        = zero
Loi(1:NLV_MAX)% friction  = zero

Matrice_ROTBVD (1:3,1)  = (/iun,izero,izero/)
Matrice_ROTBVD (1:3,2)  = (/izero,iun,izero/)
Matrice_ROTBVD (1:3,3)  = (/izero,izero,iun/)
SHIFT(1:3,1:3)          = izero
SHIFTBOITE(1:3,1:3)     = izero
DECALCLP(1:3,1:3)       = izero
Facteur_boite_rot       = iun
Homothetie              = iun
nbgrains                = izero

#if defined(PA) || defined(MDC)
 call MPI_INIT(IERR)
 !
 if (IERR /= 0) stop "Erreur d'initialisation de MPI"
#endif

#ifdef PA

! The process rank in MPI_COMM_WORLD is defined
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MON_RANG, IERR)
IF (IERR /= 0) stop "Erreur de la requete du rang"

! The size of the MPI_COMM_WORD communicator is defined
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, TAILLE, IERR)
IF (IERR /= 0) stop "Erreur de la requete de la taille"

!Initialisation des variables globales
VITSIGINT=0
VITSIGLP=0

! The default solution (Nb_Phase = 1)
TAILLE_PAIRIMPAIR = TAILLE
Ma_couleur = IZERO
Mon_Rang_PairImpair = Mon_Rang

#endif


! Definition and opening of the debug output files unit
#ifdef PA
open(10,STATUS='SCRATCH')
  write  (10,'(I2)') Mon_Rang
  rewind (10)
  read   (10,*) file_index
close(10)
file_debug= "../out/debug/proc_"//trim(file_index)//".txt"
open (unit=379, file=file_debug, status='unknown')
#else
open (unit=379, file='../out/debug/proc_0.txt', status='unknown')
#endif

#ifdef MDC

#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
#endif

! Initialisation of communication between Zeb and mM
print *, ''
print *, '(Z<->DD 0) STARTED initialisation of message passing'

call zebmpiinit

call zebmpipackrecv

call zebmpiunpack(zebchecktid)

if (zebchecktid == zebulontid) then
  call zebmpisend(magic_dislo)
else
  stop  'FATAL ERROR: mM tid is not equal to Zebulon tid'
endif

! mM is repositionned at Zebulon currentworkingdirectory
call zebmpiunpack(stringlength)
!- call zebmpiunpack(currentworkingdirectory, stringlength)
call zebmpiunpack_chararray(currentworkingdirectory, stringlength)
call  chdir(currentworkingdirectory(1:stringlength))

print *, '(Z<->DD 0) FINISHED Initialisation of message passing'
print *, ''


#if defined(PA) && defined(MDC)
    ! Le proc zero a fini d ecrire, on libere tout le monde
  endif
#endif


#endif




!#########################################################################
!# Chargement des donnees propres au materiau : Fichier 'materiau'       #
!#########################################################################

#ifdef MDC
open(1,file='../in/DCM/input.dcm',STATUS='OLD')
#else
open(1,file='../in/input.dd',STATUS='OLD')
#endif

carac = " "
do while (carac /= "-")
    read (1,*) carac
enddo

read (1,*) materiau
read (1,*) control
read (1,*) segments


close (1)

fichseg = "../in/"//segments
#ifdef MDC
print *, " Name of the MATERIAL file : ../in/DCM/"//materiau
print *, " Name of the CONTROL file : ../in/DCM/"// control
print *, " Name of the SEGMENTS file : "//fichseg

open(1,file="../in/DCM/"//materiau,STATUS='OLD')

#else

print *, " Name of the MATERIAL file : ../in/"//materiau
print *, " Name of the CONTROL file : ../in/"// control
print *, " Name of the SEGMENTS file : "//fichseg

open(1,file="../in/"//materiau,STATUS='OLD')
#endif

!**   How many phases must be considered in this material
read(1,*) Nb_phase
#ifdef MDC

if (Nb_phase > 1) print *, "To this day The DCM works only with a single phase"

#endif
! Simulation with more than one phase can only be run with mmp and pair number of nodes
test_PA = .FALSE.
#ifdef PA
test_PA = .TRUE.
if ((taille/2)*2 /= taille) test_PA = .FALSE.   ! taille must be pair
if (Nb_phase > 2) test_PA = .FALSE.             ! Only two phases can be considered for the moment
#endif
if (.not. test_PA .and. Nb_phase /= 1) then
  write(*,*) '  '
  write(*,*) 'Simulation in multi-phase mode can only be run with mmp and with a pair number of procs !!'
  stop
endif

#ifdef PA
! when Nb_phase = 1 :
! A new communicator COMM_PAIRIMPAIR is defined which is identical to MPI_COMM_WORD (Ma_Couleur = 0)
! when Nb_phase = 2 :
! A new communicator COMM_PAIRIMPAIR is defined to split MPI_COMM_WORD into pair and impair process (Ma_Couleur = 0 or 1)
if (Nb_phase == 2) then
  ! The rank in COMM_PAIRIMPAIR is in the same order than in MPI_COM_WORLD
  Clef_PairImpair = Mon_rang

  ! The process color is defined to split MPI_COM_WORLD in two equal parts
  Ma_Couleur = mod(Mon_rang,2)
!endif

! The creation of COMM_PAIRIMPAIR
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,Ma_Couleur,Clef_PairImpair,COMM_PAIRIMPAIR,IERR)

if (IERR /= 0) stop "COMM_PAIRIMPAIR initialization error"

! The process rank in COMM_PAIRIMPAIR is redefined
CALL MPI_COMM_RANK(COMM_PAIRIMPAIR, MON_RANG_PAIRIMPAIR, IERR)
! The size of COMM_PAIRIMPAIR is redefined
CALL MPI_COMM_SIZE(COMM_PAIRIMPAIR, TAILLE_PAIRIMPAIR, IERR)
endif
print*,"MPI_COMM_WORLD", Taille, taille_PAIRIMPAIR, Ma_couleur, Mon_rang, Mon_rang_PAIRIMPAIR
#endif

loop_phase : do iphase = 1,Nb_phase

  ! The pair phase of the material is treated by the pair procs and the impair phase are treated by the impair procs
  if (Nb_phase /= 1 .and. iphase > 1) then
    if ( (Ma_couleur + 1) < iphase) cycle loop_phase
  endif

  !**   The phase index (=1 if only one phase)
  read(1,*) Index_phase
#ifndef MDC
! (!!!DCM!!!) Elastic constant are defined in the FE code

  !**   Module de cisaillement MU a T = 0 K
  read(1,*) ModuleG0

  !**   pente de la courbe Module de cisaillement = fonction (Temperature)
  read(1,*) dmusurdT

  !** Module de POISSON
  read (1,*) DPOISS
#endif

  !**  tauIII = contrainte de debut du stade III (restauration)
  read(1,*) TAUIII

  !** parametre du Glissement devie
  read(1,*) BETA
  read(1,*) XALPHA

  !**  lecture de la structure cristallographique
  read(1,*) crystal_structure(1:3)

  read (1,*) Nb_slip_types

  syst=0
  do itemp = 1, Nb_slip_types
      read (1,*) Burgerstemp        ! norme vecteur Buergers en Angstrom de la famille
      write (*, '( " The Burgers Vector  = ", F6.3, " A ")') Burgerstemp
      read (1,*) Slip(itemp)%planes(1:3)
      read (1,*) Slip(itemp)%directions(1:3)

      if (etat( Slip(itemp)%planes(1:3),Slip(itemp)%directions(1:3)) /= 0) &
      stop "Problem with the elementary vectors direction!!!"

      read (1,*) Slip(itemp)%Nsystemes   ! nb de systemes de glissement

  ! Chaque systeme de glissement a sa norme vecteur de vecteur de Burgers
      do ktemp=(itemp-1)*syst+1,itemp*Slip(itemp)%Nsystemes
          Slip(ktemp)%VecBurgers=Burgerstemp
      enddo
      syst=Slip(itemp)%Nsystemes
      read (1,*)   &
      Nloi(itemp,1),  & ! numero de la loi de vitesse associe au vis
      Nloi(itemp,2),  & ! numero de la loi de vitesse associe au mixte 1
      Nloi(itemp,3),  & ! numero de la loi de vitesse associe au coin
      Nloi(itemp,4)   ! numero de la loi de vitesse associe au mixte 2
      if (Nloi(itemp,3)-Nloi(itemp,2) /= 0 .or. Nloi(itemp,4)-Nloi(itemp,2) /= 0) then
        print *, "Currently, mM does not give the possibility to use different thermally activated"
        print *, "velocity law for non screw segments !!! Additional work is needed for that!"
        print *, "Nevertheless, if you want to run a simulation, your responsibility ..."
        read(*,*)
      endif
  enddo

  ! lecture du nombre de lois de vitesse a lire
  ! Attention par defaut la loi n 1 est celle de la relaxation
  read (1,*) NLV
  do i = 1, NLV
      read(1,*) itemp  ! numero de la loi
      if (itemp /= i ) stop " Erreur fichier materiau : n loi /= ordre de la loi, FATALE !"

      read(1,*) Loi(itemp)%Arrhenius  ! clef d'activation thermique
  !   0 :  loi visqueseue ;
  !   1 : activation thermique sans effet de lognueur ;
  !   2 : double decrochement

  ! decrire la friction du reseau
      if (Loi(itemp)%arrhenius == IDEUX) then
        read(1,*) valeur
        Loi(itemp)% h = valeur     ! attac frenquency for thermal activated mobility (s-1)
        read(1,*) valeur
        Loi(itemp)% deltag0 = valeur ! energie d'activation totale ( pour T= 0 K)
        read(1,*) valeur
        Loi(itemp)%friction = valeur        ! contyrainte au plateau athermique
        read(1,*) valeur
        Loi(itemp)%tau0 = valeur
        read(1,*) valeur
        Loi(itemp)% coef_p = valeur ! exposant p de la loie de vitesse (pour les vis)
        read(1,*) valeur
        Loi(itemp)% coef_q = valeur ! exposant q de la loie de vitesse (pour les vis)
      elseif (Loi(itemp)%arrhenius == 1) then
        read(1,*) rapport
      elseif (Loi(itemp)%arrhenius == 0) then
        read(1,*) valeur
        Loi(itemp)%BF = valeur
        read(1,*) valeur
        Loi(itemp)%friction = valeur !  donnees pour les coin
      else
        stop " Velocity law unknown!!!"
      endif

  enddo

enddo loop_phase

close (1)

!########################################################
!# Loading of the simulation input parameters
!########################################################

#ifdef PA
! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  control = control(1:len_trim(control))//"2"
endif
#endif

#ifdef MDC
open(2,file="../in/DCM/"//control,STATUS='OLD')
#else
open(2,file="../in/"//control,STATUS='OLD')
#endif

#ifndef MDC
read(2,*) SIDEJA            ! Key controlling the simulation starting status (new start, restart, ..)
#endif

#ifdef MDC
read(2,*) Mode_deformation_key ! key to add deformation not seen by FE (Zebulon)
#endif

read(2,*) Mode_deformation  ! Deformation mode imposed in the simulation

#ifdef MDC
 if (mode_deformation > 9  .or. mode_deformation  < 0) &
  stop "only mode_deformation equal to 0,1,2,3,4,5,6,7,8,9 implemented in the DCM!!!"
#else
if (mode_deformation > 9 .or. mode_deformation < 0) stop "Unknown loading model!!!"
#endif


read (2,*) echelle          ! Reference scale, i.e. size of elementary screw vectors in BVD.XXX (Burgers vector unit)

#ifdef MDC
shear=.false.
#else
read(2,*) shear             !  Stress controlling key. If .true. the control is made by considering the resolved shear stress
#endif
                            ! on the slip system with the highest Schmid factor
read(2,*) Sigma0            ! The starting stress in MPa

#ifndef MDC
read(2,*) DELTAT0           ! Elementary time step of the simulation
#endif

read(2,*) SigmaPoint        ! Elementary stress increments (Pa s-1)       (mode_deformation=3,6,7)

read(2,*) EpsilonPoint      ! Imposed strain rate (s-1)                   (mode_deformation=0,5,6,7)
Epspointsign = sign(un,Epsilonpoint)

read(2,*) RAID              ! Apparent Young modulus                      (mode_deformation=0,5,6,7)
!#ifndef MDC
! In fatigue mode (5) we need also
read(2,*) EPSMAX            ! Maximum plastic strain per cycles           (mode_deformation=5)

read(2,*) FSIG   ! Signe initial (fatigue)
!#endif

!*** Options on the dislocations elementary properties
read(2,*) LINTEN                ! Key defining the type of line tension used

select case(linten)            ! control for the key validity
  case(:-1)
     Stop 'The line tension key: linten, is not correct !'
case(5:)
     Stop 'The line tension key: linten, is not correct !'
case (3)
     ! Discreet values of the anisitropic line tension are defined in the file 'disdi'
     ! Such file can be defined thanks to the program 'DisDi' of Joel Douin
     open(222,file="../in/mM_Data_Files/linetension_dat",STATUS='OLD')
       Do I=1,91
         read(222,*) Disd(1,I),Disd(2,I),Disd(3,I)
       enddo
     close(222)
endselect

read(2,*) GLDEV             !< key to switch on/off the cross-slip

read(2,*) key_nucleation    !< key to switch on/off the nucleation process

read(2,*) key_infline       !< key to switch on/off the applied stress on PBC induced infinite lines

read(2,*) key_crack         !< key to add the stress field of a sharp crack in the simulated volume

read (2,*) Z(1:3)           ! Tensile or compression uniaxial test direction (Miller indices)

! For the uniaxiale test control, the vector Z must be a unit vector
if (norvect(Z) > 0.01) Z(1:3) = Z(1:3)/DSQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))   ! Uniaxiale case only


read(2,*) TEMPERATURE       ! Temperature (Kelvin)

read(2,*) Facteur_Depmax    ! Maximum segment displacement in a simulation step (Burgers vector unit)


read(2,*) NstatControl      ! Number of steps accounted for in the simulation control procedure (mode_deformation = 0,5,6,7)


read(2,*) relax_TL          ! Number of steps ascribed to initial relaxation when involving line tension only
read(2,*) relax_INT         ! Number of steps ascribed to initial relaxation under no applied load and with no contact reactions allowed
read(2,*) relax_reac        ! Number of steps ascribed to initial relaxation under no applied loading

#ifdef MDC
! Initial relaxation steps can not be applied with MDC
if ((relax_TL + relax_INT + relax_reac) /= 0 ) then
  write(*,*) "ERROR RELAXATION in the control file"
  write(*,*) "  relax_TL   = ", relax_TL
  write(*,*) "  relax_INT  = ", relax_INT
  write(*,*) "  relax_reac = ", relax_reac
endif
#endif

#ifndef MDC
read(2,*) NSTEP             ! Number of steps of simulation
#endif

read(2,*) Ldis_act          ! Segment length (micron) at which line discretization is systematically tested  (micron)
read(2,*) Ldis_nact         ! Same as above for segments belonging to an inactive slip systems (micron)

read(2,*) Period            !  Number of steps before force on waiting (quasi-immobile) segments is recalculated
if (NLV > 2 .and. Period > 2) then
  print*, ' '
  print*, 'The Wait algorithm cannot be used with a simulation using more than one'
  print*, 'mobility law since the segments character must be identify at every steps !'
  print*, 'The Wait period must be set to 1 !'
  stop
endif

#ifndef MDC
read(2,*) KRC               ! Number of waiting steps before long-range contribution to internal stress is calculated (if Greengard method used)
#endif

read (2,*) L_Boite          ! Linear mean size of domains defined used in Greengard's method (micron)
! if negative --> Dynamic mod, an optimum number of domains is calculated preiodically during the simulation
! if larger than 1/4 of the simulated volume in the smallest linear direction = Greengard method is not used
! else, a constant number of domains is used with dimensions as cloased as possible of L_boite

#ifdef PA
if (Nb_Phase > IUN .and. L_Boite < IZERO) &
  stop "Simulation with Nb_Phase = 2 must be run with a constant number of multipoles domains"
#endif

#ifndef MDC
read(2,*) PBCIxDIM,PBCIyDIM,PBCIzDIM  ! Number of replicas used in Greengard method (in x,y and z directions),
                                      ! symmetric long rang solution is imposed if negative

read(2,*) AB_LENGTH_SEUIL   ! Maximum segment length neglected in the long range contribution (if Greengard method used) (Echelle unit see line 3)
#endif

read(2,*) DCFJ              ! Maximum Distance at which stress is calculated on the segments connected to a junction (unit: Burgers vector)

read(2,*) GB                ! interface and surface definitions: 0-> inactivates, 5-> planes, 2->spherical, 3->regular 3D network (polycristal)

read(2,*) TauINT_LIMITE     ! Critical stress at which segments are considered as in a singular field (MPa)

read(2,*) KISAUVE           ! Write periodicity of segment configurations and information needed to restart computation

read(2,*) KSTAT             ! Write periodicity of results

read(2,*) KKIM              ! Write periodicity of the trajectory film
if (mod(kisauve,kkim) /= izero) &
stop ' Stop !!!! the parameters kisauve and kkim must verify mod(kisauve,kkim) /= 0'

read(2,*) KPREDRAW          ! Periodicity of refreshment of the graphical interface in gmm mode

read(2,*) Shift_rotation    ! Key for translation and rotation of the simulated volume (see the shift_rotation file)

! Debug parameters
read(2,*) iterinfo          ! Step of debugging (no debugging if negative)
read(2,*) sysinfo           ! Slip system of interest in debugging procedure (needed in simulations with many segments)

#ifdef MDC

 read(2,*) halfthickness         !half-thickness of the inclusion in which the eignestrain is defined

#endif


close(2)

#ifdef MDC

#if defined(MDC) && defined(PA)
 ! dpoiss_sum = zero
 ! xmu_sum = zero
  if (Mon_Rang == IZERO) then
#endif

  open(66,FILE='../out/debug/MDCinfo.txt',STATUS='unknown')
  !** Gets the elastic tensor from Z
  call elasticity_matrix ! Contains Z ->DD No. 3

  !** Calculates the shear modulus (Reuss average for anisotropic elasticity)
  call FEshear_modulus

  !** Gets a restart key from Z
  call restart ! Contains Z ->DD No. 4

  close(66)
#if defined(MDC) && defined(PA)
  endif

CALL MPI_BCAST(sideja,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(dpoiss,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)


#endif

#endif

!########################################################
!# Loading of the simulation input parameters - crack
!########################################################
if(key_crack) then

  open(666,file="../in/crack_def",STATUS='OLD')

  do
    read(666,'(a1)') carac
    if (carac == "!") exit
  enddo

  read(666,*) CrackTip(1:3)     ! Crack tip coordinates
  print *, "###################################### IMPORTANT ######################################"
  print *, "Crack Key has been activated - A crack field is added to the homogeneous applied stress!"
  print *, "###################################### IMPORTANT ######################################"
  if (Z(1)/=IZERO .or. Z(2)/=IZERO) &
  stop "STOP - So far crack field can only be activate with a load along z only"
  if (mode_deformation < 3 .or. mode_deformation > 4) &
  stop "STOP - So far crack field can only be activate with mode deformation 3 or 4"

  read(666,*) fact_type   ! Geometrical factor method
  read(666,*) acrack      ! Crack length
  read(666,*) plate_H     ! Plate HALF height
  read(666,*) plate_W     ! Plate weigth

!  write (fact_type_str,'(I15)') fact_type
  open(10,STATUS='SCRATCH')
    write  (10,'(I15)') fact_type
    rewind (10)
    read   (10,*) fact_type_str
  close(10)
!  write (plate_a_str,'(I15)') acrack
  open(10,STATUS='SCRATCH')
    write  (10,'(I15)') acrack
    rewind (10)
    read   (10,*) plate_a_str
  close(10)
!  write (plate_H_str,'(I15)') plate_H
  open(10,STATUS='SCRATCH')
    write  (10,'(I15)') plate_H
    rewind (10)
    read   (10,*) plate_H_str
  close(10)
!  write (plate_W_str,'(I15)') plate_W
  open(10,STATUS='SCRATCH')
    write  (10,'(I15)') plate_W
    rewind (10)
    read   (10,*) plate_W_str
  close(10)

  !Calculate geometric factor
  cmd = "python ../src/outils/crack_problem/SIF_fact_calculation.py"//' '// &
      & trim(fact_type_str)//' '//trim(plate_a_str)//' '//trim(plate_H_str)//' '//trim(plate_W_str)

  write(*,*) '---> ',cmd

  call system(cmd)

  open(667,file="../out/siffactresult.dat",STATUS='OLD')
    read(667,*) SIF_fact
  close(667)

  print *, "######################################    END    ######################################"
endif

close(666)
!########################################################

#ifdef PA
! Steps of relaxations with only the line tension are not possible in parallel computations
if (relax_TL /= izero) stop "relax_TL /= zero is not possible with mmp or mdcp!"

! Simulation with the APB procedure can not be applied in the mmp mode for the moment
if (GB == IQUATRE) stop "Simulation with the APB procedure can not be applied in the mmp mode for the moment!"
#endif


if (Shift_ROTATION > itrois .or. Shift_ROTATION < izero) then
   stop " Unknown boundary conditions!!!"
else
  select case(Shift_Rotation)
    ! Shift_Rotation = 0 : CLP sans  shift et sans rotation
    case (IZERO)
      decalageCPL    = .false.
      RotationBVD = .false.

    ! Shift_Rotation = 1 : CLP avec  decalageCPL et sans rotation
    case (IUN)
      open(2,file="../in/shift",STATUS='OLD')

!*** Shifts pour CLP shiftees
      read(2,*) DECALCLP(1:3,1) ! Shift en face a
      read(2,*) DECALCLP(1:3,2) ! Shift en face b
      read(2,*) DECALCLP(1:3,3) ! Shift en face c
      close(2)
      decalageCPL    = .true.
      RotationBVD = .false.

    ! Shift_Rotation = 2 :  CLP sans  shift et avec rotation
    case (IDEUX)
      decalageCPL    = .false.
      RotationBVD = .true.

      file_mat_rot = "../in/matrotation_def"

      ! A second set of input files is needed when Nb_phase = 2
#ifdef PA
      if (Nb_phase == 2 .and. Ma_Couleur == 1) then
        file_mat_rot = file_mat_rot(1:len_trim(file_mat_rot))//"2"
      endif
#endif

      open (33,file=file_mat_rot,status='old')     !n de fichier mis arbitriarement
      read (33,*) Matrice_rotBVD(1,1:3)
      read (33,*) Matrice_rotBVD(2,1:3)
      read (33,*) Matrice_rotBVD(3,1:3)
      read (33,*) facteur_boite_rot   ! en tournant la BVD , on agrandit les vecteur elementaire
                                      ! il faut modifier facteur_boite
      read (33,*) homothetie          ! en tournant la BVD , on agrandit les vecteur elementaire: facteur d'ohomothetie

      ! Matrix determinant calculation: if the determinant is equal to -1, the matrix is not correct, because gamma will
      ! thus be inverted
      det =(  Matrice_rotBVD(1,1) * (Matrice_rotBVD(2,2) * Matrice_rotBVD(3,3) - Matrice_rotBVD(3,2) * Matrice_rotBVD(2,3))  &
            - Matrice_rotBVD(1,2) * (Matrice_rotBVD(2,1) * Matrice_rotBVD(3,3) - Matrice_rotBVD(3,1) * Matrice_rotBVD(2,3))  &
            + Matrice_rotBVD(1,3) * (Matrice_rotBVD(2,1) * Matrice_rotBVD(3,2) - Matrice_rotBVD(2,2) * Matrice_rotBVD(3,1))) &
            /(homothetie*homothetie*homothetie)

      if (det < 0.99 .or. det > 1.01) then
          print *, "the determinant of the rotation matrix must be equal to 1"
          print *, "the determinant calculated =",det
          stop
      endif

      print *, " ============================================================================"
      print *, " ===============        Rotation des vecteur de base de discretisation"
      print *, " Facteur de boite :"     ,facteur_boite_rot
      print *, " Facteur d homothetie :" ,homothetie
      print *, " The rotation matrix :"
      write (*, '( " |  ",3I7,"  |   ")') matrice_rotBVD(1,1:3)
      write (*, '( " |  ",3I7,"  |   ")') matrice_rotBVD(2,1:3)
      write (*, '( " |  ",3I7,"  |   ")') matrice_rotBVD(3,1:3)
      print *, " ======================================================================="

    ! Shift_Rotation = 3 :  CLP avec  shift et avec rotation
    case (ITROIS)
      decalageCPL    = .true.
      RotationBVD = .true.

  endselect

endif
!==========================================================================================
!==========================================================================================

end subroutine lire_donnees

!################################################################################
!# This subroutine is used to read a reference stress tensor defined in the
!# input file "../in/tensapp". This solution is considered when UNIAXIALE loading
!# conditions are not used.
!################################################################################
subroutine lire_tenseur

implicit none

real(kind=DP) rval,Ztrans(3)
character(len=40)    :: file_tensapp

file_tensapp = "../in/tensapp_def"

#ifdef PA
! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_tensapp = file_tensapp(1:len_trim(file_tensapp))//"2"
endif
#endif

open(2,file=file_tensapp,STATUS='OLD')

! The loading tensor
read(2,*) tensapp(1,1:3)
read(2,*) tensapp(2,1:3)
read(2,*) tensapp(3,1:3)

! A direction used to calculate deformation in a reference axis direction
! Such direction is used as a macro strain gauge to control plastic strain
! in a particular direction
read(2,*) Z(1:3)
Z(1:3) = Z(1:3)/DSQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))

! The loading tensor must be length conservative in the strain gauge direction
Ztrans=matmul(tensapp,Z)
rval = DSQRT(Ztrans(1)*Ztrans(1)+Ztrans(2)*Ztrans(2)+Ztrans(3)*Ztrans(3))

if(rval == zero) then
  write(*,*) "The strain gauge direction given in the file tensapp_def"
  write(*,*) "is incompatible with the applied stress tensor!!!"
  stop
else
  ! The applied stress tensor is corrected to be length conservative in the Z direction
  tensapp = tensapp/rval
endif

close(2)

end subroutine lire_tenseur

!################################################################################
!# /brief  The file containing the imposed stress field is loaded
!################################################################################
subroutine lire_FE_SigappS

implicit none

Integer(kind=DPI)     :: loopx,loopy,loopz
integer(kind=DPI)     :: Numb_points
integer(kind=DPI)     :: npi
real(kind=DP)         :: realdimFE_Sig(3)
real(kind=DP)         :: Stress_fact_tmp
character(len=40)     :: file_FE_Sig

!PrintFieldVTK tabs
real(kind=DPI),allocatable      :: XCooR(:),YCoor(:),ZCoor(:)   !< Box Center Coordinates
real(kind=DPI),allocatable      :: FieldAmp(:,:)                !< The field amplitude we want to trace
character(len=24)               :: file2print                   !< The name and place of the file to print

file_FE_Sig = "../in/FE_Sig_def"

#ifdef PA
! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_FE_Sig = file_FE_Sig(1:len_trim(file_FE_Sig))//"2"
endif
#endif

open(2,file=file_FE_Sig,STATUS='OLD')

! The stress fact to applied to the stress field
read(2,*) Stress_fact_tmp

! If stress_fact is set to zero, this mean that we want to control the loading with the FE_Sig file definition
! but at constant strain rate. Stress_fact must then be initialized to one
if (Stress_fact_tmp == zero) then
  StrainRateMode8 = .true.
else
  StrainRateMode8 = .false.
endif

! The initial value of Stress_fact must not be load if we restart the simulation
if (sideja /= IUN)  Stress_fact = Stress_fact_tmp

! The parameters for StrainRateMode8 = .true. are always loaded
read(2,*) ref_deform        ! The reference deformation amplitude need for normalization in strain rate control
read(2,*) stress_fmin       ! The min factor applied to the stress_fact value in strain rate control mode
read(2,*) stress_fmax       ! The mac factor applied to the stress_fact value in strain rate control mode

!region of the simulation box in which the FE stress field is defined
read(2,*) FEgrid_or(1:3)  !coordinate of the corner closer to simulation box origin
read(2,*) FEgrid_size(1:3) !size of the grid in X,Y and Z direction in simulation box unit

! The Tab dimensions
read(2,*) realdimFE_Sig(1:3)
dimFE_Sig (1)= int(realdimFE_Sig(1))
dimFE_Sig (2)= int(realdimFE_Sig(2))
dimFE_Sig (3)= int(realdimFE_Sig(3))

allocate (FE_SigappS(dimFE_Sig(1),dimFE_Sig(2),dimFE_Sig(3),6))

print *,'reading external stress field from the file FE_Sig'
print *,'dimension of FE_Sig',dimFE_Sig(1),dimFE_Sig(2),dimFE_Sig(3)
! The FE_SigappS tab is loaded
do loopz=1,dimFE_Sig(3)
  do loopy=1,dimFE_Sig(2)
    do loopx=1,dimFE_Sig(1)
      read(2,*) FE_SigappS(loopx,loopy,loopz,1:6)
      ! Stress field is defined in reduce unit
      FE_SigappS(loopx,loopy,loopz,1:6) = FE_SigappS(loopx,loopy,loopz,1:6) * 1.D6 / XMU
    enddo
  enddo
enddo

close(2)

! Lets print the stress field components in VTK files to do some visualization if needed
Numb_points = dimFE_Sig(1) * dimFE_Sig(2) * dimFE_Sig(3)

allocate(XCooR(Numb_points))
allocate(YCooR(Numb_points))
allocate(ZCooR(Numb_points))
allocate(FieldAmp(Numb_points,7))

npi = 0
do loopz=1,dimFE_Sig(3)
  do loopy=1,dimFE_Sig(2)
    do loopx=1,dimFE_Sig(1)

      npi = npi + 1
      XCooR(npi) = real(FEgrid_or(1)) + real(FEgrid_size(1)) * ((real(loopx)-0.5)/real(dimFE_Sig(1)))
      YCooR(npi) = real(FEgrid_or(2)) + real(FEgrid_size(2)) * ((real(loopy)-0.5)/real(dimFE_Sig(2)))
      ZCooR(npi) = real(FEgrid_or(3)) + real(FEgrid_size(3)) * ((real(loopz)-0.5)/real(dimFE_Sig(3)))
      FieldAmp(npi,1:6) = FE_SigappS(loopx,loopy,loopz,1:6) * XMU / 1.D6
      FieldAmp(npi,7) = 1./sqrt(6.)* sqrt(                                                          &
                        (FE_SigappS(loopx,loopy,loopz,1)-FE_SigappS(loopx,loopy,loopz,2))**2 +      &
                        (FE_SigappS(loopx,loopy,loopz,2)-FE_SigappS(loopx,loopy,loopz,3))**2 +      &
                        (FE_SigappS(loopx,loopy,loopz,3)-FE_SigappS(loopx,loopy,loopz,1))**2 +      &
                        (FE_SigappS(loopx,loopy,loopz,4)**2+                                        &
                         FE_SigappS(loopx,loopy,loopz,5)**2+                                        &
                         FE_SigappS(loopx,loopy,loopz,6)**2) * 6.)
    enddo
  enddo
enddo

file2print = '../out/FE_SigappS_11.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,1),file2print)
file2print = '../out/FE_SigappS_22.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,2),file2print)
file2print = '../out/FE_SigappS_33.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,3),file2print)
file2print = '../out/FE_SigappS_32.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,4),file2print)
file2print = '../out/FE_SigappS_31.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,5),file2print)
file2print = '../out/FE_SigappS_21.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,6),file2print)
file2print = '../out/FE_SigappS_VM.vtk'
call PrintFieldVTK(Numb_points,XCooR,YCoor,ZCoor,FieldAmp(:,7),file2print)

deallocate(XCooR)
deallocate(YCooR)
deallocate(ZCooR)
deallocate(FieldAmp)


end subroutine lire_FE_SigappS

!################################################################################
!> \ingroup creep
!> In this subroutine, the parameters controlling the dislocation climb velocity are defined.
!! The input data file name is "../in/climb_def".
!################################################################################
subroutine lire_climb

implicit none

real(kind = DP)     :: Umig         !< var vacancy migration energy (eV)
real(kind = DP)     :: Dv0          !< diffusivity at T = 0K
real(kind = DP)     :: Dv           !< the vacancy bulk diffusion coefficient
real(kind = DP)     :: Uform        !< vacancy formation energy (eV)
real(kind = DP)     :: c0           !< equilibrium vacancy concentration in a defect-free crystal

character(len=40)   :: file_climb   !< the name of the file containing the creep data

! initialization
file_climb    = "../in/climb_def"
TwoPi         = deux * Pii
InvVecBurgers = un / VecBurgers
UnKT_Joule    = UnKT / unev

open(2,file=file_climb,STATUS='OLD')

! parameters for equilibrium vacancy concentration:
!   Uform : vacancy formation energy
read(2,*) Uform
  c0 = exp(-UnkT * Uform)       ! equilibrium vacancy concentration in a defect-free crystal


! parameters for vacancy bulk diffusion:
!   Umig  : vacancy migration energy
!   Dv0   : diffusivity at T = 0K
read(2,*) Umig
read(2,*) Dv0
  Dv = Dv0 * exp(-UnkT * Umig)  ! the vacancy bulk diffusion coefficient

! satfac : defined through cinfinity = satfac * c0
!          approximate model assumes constant remote concentration of vacancies
read(2,*) satfac

! the atomic volume
read(2,*) omega

! lnrrc  : logarithm term appearing in the velocity equation
read(2,*) lnrrc

! Calculation of the prefactor uses in the vclim equation
prefac2 = TwoPi * InvVecBurgers * Dv * c0 / lnrrc

close(2)

end subroutine lire_climb

!#############################################################################
!# determine les parametre principaux de la simulation
!#############################################################################
subroutine Parametrage

Implicit NONE

INTEGER(DPI) :: i
real(DP)     :: maille_mat
#ifndef MDC
real(DP)     :: resolution
#endif

print *, " "
print*, "Essential parameters fixing the stack size"
print*, "  REAL Double precision are coded with ", DP, " Octets"
print*, "  INTEGER long are coded with ", DPI," Octets "
print*, "  Maximum number of segments NSEGMAX equal ", nsegmax

! Initialization of the time steps
if (sideja /= IUN) then
  Relax_INT = Relax_INT + Relax_TL
  Relax_Reac = relax_Reac + Relax_INT
  NStep = NStep + Relax_Reac
  StepControl = Relax_Reac + NstatControl
  StepNextTimeBreak = StepControl

  print *, " "
  print *, "Relaxation steps"
  print *, " "
  print *, "  First relaxation run (LT only)                 : ", Relax_TL
  print *, "  Second relaxation run (+ elastic interactions) : ", Relax_INT
  print *, "  Third relaxation run (+ contact reactions)     : ", relax_Reac
  print *, "  Total number of simulation steps               : ", Nstep
  if (relax_Reac /= 0) then
    print*, " "
    print *, "  During relaxation, the velocity law is /1/ for all the segments"
  endif

  if(sideja == ideux) then
    print *, " "
    print *, "=============== Starting simulation with particles ================"
    print *, " "
    particules = .true.
  else
    npar = Izero
  endif

  if(sideja == itrois) then
    print *, " "
    print *, "=============== Starting simulation with small prismatic loops ================"
    print *, " "
    Loops =  .true.
endif

endif

! there are only 3 cases for cartograph
if (mode_deformation == 1) then
! the code is used to generate the carto files seg???? in in/carto/
   cartograph = 1
elseif(mode_deformation == 2) then
! the code is used through a script, carto or gcarto, to scan the carto files in in/carto/
! to make the appropriate statistics
   cartograph = 2
else
! the code doesn t do anything related to cartographies
   cartograph = 0
endif

Deltat = Deltat0

write(*,*) " "
write(*,'(" Simulation temperature = ", I5," K")') nint(temperature)
VecBurgers = Slip(1)%Vecburgers

#ifndef MDC
! =======================================================================
!  Calcule du module de cisaillement a la temperature de simulation TEMP
if (temperature /= zero) then
UnkT = UN / Boltzmann / temperature
else
   UnkT = 1000000.0
endif

XMU = ModuleG0 + dmusurdT * TEMPERATURE      ! module en GPa
XMU = XMU * 1.D9                             ! module en Pa

!The isotropic Stiffness tensor is defined
call MakeStiffTens
#endif

!==========================================================================================
!====== Attention : ici on definit toutes variable cles du programme =====================
!==========================================================================================
! =======================================================================
! Definition de tous les parametres relies a la cristallographie
! FACTEUR_CRIS le rapport entre le modul du vecteur vis (ex. [440] et le module du
! vecteur de Burgers dans le repere cristallo(ex. 1/2 [110])

! Initilisation des clefs choix de la cristallo
DC = .FALSE.; CS = .false.; BCC = .false.; CFC = .false.; HCP = .false.; ORT = .false.; MGO = .false.
If (crystal_structure(1:2) == 'CS') then
   CS = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_CS'
   facteur_cris = FACTEUR_CS
   fichbase = '../out/base.cs'
   fichtab  = '../out/jonction.cs'
   facteur_boite = facteur_boite_CS
   angle_vis = angle_vis_HC
   NbSysDev =  NbSysDev_CS
elseif (crystal_structure(1:3) == 'CFC')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_CFC'
   CFC = .true.
   facteur_cris = FACTEUR_CFC
   fichbase = '../out/base.cfc'
   fichtab  = '../out/jonction.cfc'
   facteur_boite = facteur_boite_CFC
   angle_vis = angle_vis_CFC
   NbSysDev =  NbSysDev_CFC
elseif (crystal_structure(1:3) == 'DC')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_DC'
   CFC = .true.
   facteur_cris = FACTEUR_DC
   fichbase = '../out/base.dc'
   fichtab  = '../out/jonction.dc'
   facteur_boite = facteur_boite_DC
   angle_vis = angle_vis_DC
   NbSysDev =  NbSysDev_DC
elseIf (crystal_structure(1:3) == 'BCC')   then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_BCC'
   BCC = .true.
   facteur_cris =  FACTEUR_BCC
   fichbase = '../out/base.bcc'
   fichtab  = '../out/jonction.bcc'
   facteur_boite = facteur_boite_BCC
   angle_vis = angle_vis_BCC
   NbSysDev =  NbSysDev_BCC
elseIf (crystal_structure(1:3) == 'HCP')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_HC'
   HCP = .true.
   facteur_cris =  FACTEUR_HC
   fichbase = '../out/base.hc'
   fichtab  = '../out/jonction.hc'
   facteur_boite = facteur_boite_HC
   angle_vis = angle_vis_HC
   NbSysDev =  NbSysDev_HCP
elseif (crystal_structure(1:3) == 'ORT') then
   ORT = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_ORT'
   facteur_cris = FACTEUR_ORT
   fichbase = '../out/base.ort'
   fichtab  = '../out/jonction.ort'
   facteur_boite = facteur_boite_ORT
   angle_vis = angle_vis_ORT
   NbSysDev =  NbSysDev_ORT
elseIf (crystal_structure(1:3) == 'MGO')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_MGO'
   MGO = .true.
   facteur_cris =  FACTEUR_MGO
   fichbase = '../out/base.MGO'
   fichtab  = '../out/jonction.MGO'
   facteur_boite = facteur_boite_MGO
   angle_vis = angle_vis_MGO
   NbSysDev =  NbSysDev_MGO
else
   stop " Crystal symetry not yet defined,  BY"
endif

if(RotationBVD) then
   Facteur_boite = INT(facteur_boite_rot)
   Facteur_cris = INT(Facteur_cris * homothetie)
endif

write(*,*) "Facteur_Boite = ", Facteur_boite

angle_vis = angle_vis * PII / 180.0D0
if(.not. CS .and. .not. BCC .and. .not. HCP .and. .not. CFC .and. .not. ORT .and. .not. MGO .and. .not. DC) &
   stop "The parameter angle_vis is not defined !"

! Definition of the critical angle used to defined the screw, edge and mixte statistics
Critic_angle = quatre * angle_vis

!==========================================================================================
write(*, '("  The crystal structure factor, FACTEUR_CRIS = ", I5)') facteur_cris
write(*,*)
write(*,*) "Important simulation lengths"
!==========================================================================================

! =======================================================================
! Dimensionnement des tableaux de base de vecteurs
! =======================================================================

! Calcul du nombre total de systemes de glissement
NTSG = 0
do i = 1, Nb_slip_types
    NTSG = NTSG + Slip(i)%Nsystemes
enddo

! ATTENTION ICI definition de la dimension de la base de vecteurs de discreti
! calcul de la dimension des tableaux des vecteurs unitaires de discretisation
NBASE = NBASERED * NTSG

slip(1:ntsg)%VecBurgers = slip(1:ntsg)%VecBurgers * 1.D-10    ! mise a l'echelle reelle de vecteur de burgers
if (CS) VecBurgers=slip(2)%VecBurgers
if (ORT .or. CFC .or. HCP .or. BCC .or. MGO .or. DC) VecBurgers=slip(1)%VecBurgers

if (ORT .or. CS) then
   maille_mat=VecBurgers
elseif (HCP .or. CFC .or. MGO .or. DC) then
   maille_mat = VecBurgers * 1.4142135 ! calcul du paramatre de maille physique
elseif (BCC) then
   maille_mat = VecBurgers * 1.1547D0 ! calcul du paramatre de maille physique
else
   stop ' relation Burgers-parametre de maille non connue'
endif

! =======================================================================
! =======================================================================
! ====================  Attention : tres important      =================
! =======================================================================
#ifndef MDC
resolution = echelle * VecBurgers   ! resolution effective de la simulation
vitesse_unitaire = resolution/deltat
#endif
!==========================================================================================
!==========================================================================================
!==========================================================================================
! Calcul de l'unite de la simulation pour avoir la taille des vecteurs
! vis de la simulation (e.g.|[022]|*avalue) egale a la resolution physique
! souhaitee, c-a-d egale a echelle*Vecteur Burgers reel du materiau
! Ceci permet de conserver tous les parametre de la simulation inchanges
avalue = echelle * maille_mat / facteur_cris    ! donne Burgers(simu)= echelle*Burgers(vrai)
write (*,'("  The elementary reference length, AVALUE = ", F8.5," nm" )') avalue * 1D9

BdivA = vecburgers / Avalue    ! The definition of VecBurgers norme in avalue unit
Bspread2 = (spread_core_radius * BdivA)*(spread_core_radius * BdivA)
write (*,'("  The important length, BdivA = ", F8.5, " (avalue unit)")') BdivA
write (*,'("  The spreading radius used in shrt range interactions, Bspread = ", F8.5, " (avalue unit)")') sqrt(Bspread2)

!==========================================================================================
! initialisation des deplacements maximaux
!==========================================================================================

! =========================================================
! Definition de la longuerur de discretisation homogeneisee
write(*, '("  Loaded systems discretization length (micron)     = ", F12.6)') Ldis_act
write(*, '("  Unloaded systems discretization length (micron) = ", F12.6)') Ldis_nact
Ldis_act = Ldis_act * 1D-6   ! longueur de discretisation en metre (sys actif)
Ldis_nact = Ldis_nact * 1D-6   ! longueur de discretisation en metre (sys non actif)
Xlomax    = NINT(Ldis_act/avalue,DPI) ! longueur de discretisation en (a)
Xlomax_nact  = NINT(Ldis_nact/avalue,DPI) ! longueur de discretisation en (a) (sys non actif)

! Ecriture a  l ecran pour verification des comformite
write(*, '("  Loaded   systems discretization length (in AVALUE) = ", I10)') Xlomax
write(*, '("  Unloaded systems discretization length (in AVALUE) = ", I10)') Xlomax_nact
!write(*, '("  Jogs dragging effectiveness  = ", F5.0," %" )') Efficacite_crans

kkdebug = (iterinfo == IUN)

!The anisotric line tension need to be expressed in reduce units
if (linten==ITROIS) then
Do I=1,91
    Disd(3,I) = Disd(3,I)* 1.e-9 /(Xmu*avalue*avalue)
enddo
  DisdInc = 1./Disd(1,2) ! The reverse of the angle increment used in Disd
endif

end subroutine parametrage

!#########################################################################
!# Sous-programme necessaire pour l'initialisation des fichier en fonction
!# de la clef (SIDEJA) : demarrage/redemarrage de la simulation
!#########################################################################
Subroutine initialisation_fichiers

implicit none

integer,parameter   :: DPI_S=selected_int_kind(9)   ! The integer format needed for the film.bin file
integer(DPI)             :: cristallo,JK

#ifdef PA
if ((Mon_Rang_PairImpair + IUN) >  IUN) then
  ! Only the process zero of the communicator can write in output files, the others are waiting
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
else
#endif

!#######################################################################
!# Lecture des fichiers contenants les segments et de la table des
!# voisins lorsque l'on repart d'une precedante simulation
!#######################################################################

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! opening of the output files

      file_graph    =   "../out/graph"
      file_stat     =   "../out/stat"
      file_rau      =   "../out/rau"
      file_raujonc  =   "../out/raujonc"
      file_gamma    =   "../out/gamma"
      file_gammap   =   "../out/gammap"
      file_eps      =   "../out/eps"
      file_travapp  =   "../out/travapp"
      file_travint  =   "../out/travint"
      file_bigsave  =   "../out/bigsave"
      file_film     =   "../out/film"
      file_timebreak=   "../out/timebreak"
      ! gammabox will be written only if allocation dynamic for Greengard method is false
      ! (not defined at this step)
      if (L_boite > ZERO) file_gammabox = "../out/gammabox"


#ifdef PA
      ! A second set of input files is needed when Nb_phase = 2
      if (Nb_phase == 2 .and. Ma_Couleur == 1) then
        file_graph    =   file_graph(1:len_trim(file_graph))//"2"
        file_stat     =   file_stat(1:len_trim(file_stat))//"2"
        file_rau      =   file_rau(1:len_trim(file_rau))//"2"
        file_raujonc  =   file_raujonc(1:len_trim(file_raujonc))//"2"
        file_gamma    =   file_gamma(1:len_trim(file_gamma))//"2"
        file_gammap   =   file_gammap(1:len_trim(file_gammap))//"2"
        file_travapp  =   file_travapp(1:len_trim(file_travapp))//"2"
        file_travint  =   file_travint(1:len_trim(file_travint))//"2"
        file_bigsave  =   file_bigsave(1:len_trim(file_bigsave))//"2"
        file_film     =   file_film(1:len_trim(file_film))//"2"
        if (L_boite > ZERO) file_gammabox = file_gammabox(1:len_trim(file_gammabox))//"2"
      endif
#endif

      file_graph    =   file_graph(1:len_trim(file_graph))//".txt"
      file_stat     =   file_stat(1:len_trim(file_stat))//".txt"
      file_rau      =   file_rau(1:len_trim(file_rau))//".txt"
      file_eps      =   file_eps(1:len_trim(file_eps))//".txt"
      file_raujonc  =   file_raujonc(1:len_trim(file_raujonc))//".txt"
      file_gamma    =   file_gamma(1:len_trim(file_gamma))//".txt"
      file_gammap   =   file_gammap(1:len_trim(file_gammap))//".txt"
      file_travapp  =   file_travapp(1:len_trim(file_travapp))//".txt"
      file_travint  =   file_travint(1:len_trim(file_travint))//".txt"
      file_bigsave  =   file_bigsave(1:len_trim(file_bigsave))//".bin"
      file_film     =   file_film(1:len_trim(file_film))//".bin"
      file_timebreak=   file_timebreak(1:len_trim(file_timebreak))//".txt"
      if (L_boite > ZERO) file_gammabox = file_gammabox(1:len_trim(file_gammabox))//".bin"


select case (SIDEJA)

  !############################
  case (IZERO,ideux)  ! Demarrage standart

    !#### The new graph output file
    open( 8,FILE=file_graph,STATUS='REPLACE')
    ! The file header
    write(8,'(13(4x,a))') 'KK','EPS_TOT','EPS_PLA','SIGMA','GAMMA','TAU','EPS_DOT', &
                        & 'RHO_TOT','ALPHA','ZP_X','ZP_Y','ZP_Z','TIME'
    close( 8)

    !#### The new stat output file
    open(99,FILE=file_stat,STATUS='REPLACE')
    ! The file header
    write(99,'(22(4x,a))') 'EPS_PLA','JONC','GD','LL','SWEEP_S','SWEEP_E','SWEEP_M','RHO_1',    &
                         & 'RHO_2','LENGTH_1','LENGTH_2','TAU_1','TAU_2','DensLoop','DensInf', &
                         & 'DensFreeInf','Nloop','NInfline','GD1','GD2', 'OldGD1', 'OldGD2'
    close(99)

    !#### The new rau output file
    open(28,FILE=file_rau,STATUS='REPLACE')
    ! The file header
    write(28,'(4x,a,4x,a,21(4x,i2,a))') 'EPS_PLA','RHO_GROSEG',(JK,'_RHO', JK=1,NTSG)
    close(28)

    !#### The new eps output file
    open(16,FILE=file_eps,STATUS='REPLACE')
    ! The file header
    write(16,'(13A14)') 'kk','ep11','ep22','ep33','ep23','ep13','ep12', &
                      & 'sig11','sig22','sig33','sig23','sig13','sig12'
    close(16)

    !#### The new raujonc output file
    open(29,FILE=file_raujonc,STATUS='REPLACE')
    ! The file header
    write(29,'(4x,a,40(4x,i2,a))') 'EPS_PLA',(JK,'_RHO_J', JK=1,NTSG),(JK,'_N_J', JK=1,NTSG)
    close(29)

    !#### The new gamma output file
    open(38,FILE=file_gamma,STATUS='REPLACE')
    ! The file header
    write(38,'(4x,a,20(4x,i2,a))') 'EPS_PLA',(JK,'_GAM', JK=1,NTSG)
    close(38)

    !#### The new travapp output file
    open(48,FILE=file_travapp,STATUS='REPLACE')
    ! The file header
    write(48,'(4x,a,20(4x,i2,a))') 'EPS_PLA',(JK,'_TRA_APP', JK=1,NTSG)
    close(48)

    !#### The new travint output file
    open(49,FILE=file_travint,STATUS='REPLACE')
    ! The file header
    write(49,'(4x,a,20(4x,i2,a))') 'EPS_PLA',(JK,'_TRA_INT', JK=1,NTSG)
    close(49)

    !#### The new gammap output file
    open(25,FILE=file_gammap,STATUS='REPLACE')
    ! The file header
    write(25,'(4x,a,20(4x,i2,a))') 'EPS_PLA',(JK,'_GAM_DOT', JK=1,NTSG)
    close(25)

    !#### The new timebreak output file
    if (mode_deformation == 9) then
      open(79,FILE=file_timebreak,STATUS='REPLACE')
      ! The file header
      write(79,'(2x,3(4x,a))') 'step','timebreak','vitesse'
      close(79)
    endif

    !#### The new gammabox output file
    if (L_boite > ZERO) then
      open(61,FILE=file_gammabox,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=KODERR)
      close(61)
    endif

    !&&&& multicristals outfiles
    if (GB /= 0) then
#ifdef PA
      if (Nb_Phase == 2) stop  "GB = 3 is not possible with Nb_Phase =2"
#endif

      open(47,FILE='../out/polydens.txt',STATUS='REPLACE')
      close(47)

      open(46,FILE='../out/polyrotation.txt',STATUS='REPLACE')
      close(46)

      open(45,FILE='../out/polygamma.txt',STATUS='REPLACE')
      close(45)

    endif

    !&&&& multicristals outfiles
#ifdef MDC
    open(67,FILE='../out/debug/MDCstepstat.txt',STATUS='REPLACE')
    write(67,'(2x,3(4x,a))') 'kk','Nb_segm_dep','Nb_segm_tot'
    close(67)
#endif

    !#### The new bigsave output file
    open(99,FILE=file_bigsave,FORM='UNFORMATTED',STATUS='REPLACE')
    close(99)

    !#### The new film output file
    open(93,FILE=file_film,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=KODERR)

    if (CS)  cristallo = 1
    if (BCC) cristallo = 2
    if (CFC) cristallo = 3
    if (HCP) cristallo = 4
    if (ORT) cristallo = 5
    if (MGO) cristallo = 6
    if (DC)  cristallo = 7

    write (93) int(cristallo,DPI_S)
    write (93) int(nsegmax,DPI_S)

  !############################
  case (IUN) ! The restart case

    ! All the output files are suppose to exist, so nothing to do except for the film output file
    ! The film.bin file must be re-open and positioned at the restart position
    ! This operation is made latter in the subroutine generation_microstructure

    ! Ouverture des unites d'ecriture
    open(93,FILE=file_film,FORM='UNFORMATTED',STATUS='old',ACTION="READWRITE",IOSTAT=KODERR)

  end select

#ifdef PA
  ! Proc zero ended is task, we can now liberate everybody
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
endif
#endif

end subroutine Initialisation_fichiers


!#############################################################################
!# Dimensionnement de la base des segments unitaire-VECTEUIRS ALLOUABLES
!# SOUS-PROGRAMME DIMENSIONNER_BASERED
!################################################################# 05/11/98 ##
subroutine Allocation

implicit none

allocate (BVECLIN(3,NBASE))
allocate (BVECLINCP(3,NBASE))
allocate (BVECDEP(3,NBASE))
allocate (BVECNOR(3,NBASE))
allocate (courbuO(NBASE,2))
allocate (courbuE(NBASE,2))
allocate (ASSOC(NBASE,NBASERED))
allocate (INVSI(NBASERED))
allocate (DEASSOC(NBASE))
allocate (TYSEG(NBASE))
allocate (SISEG(NBASE))
allocate (SEGdev(NBASE,NBsysdev))
allocate (SEGdevIJ(NBASE,NBase))
allocate (CONEC(NBASE,2))
allocate (DECONEC(NBASE,NBASE))
allocate (SC(nbase,nbase))
allocate (CoefSC(nbase,nbase))
allocate (AXEJONCi(NBASE,NBASE))
allocate (AXEJONCj(NBASE,NBASE))
allocate (SYSOBST(NBASE,NBASE))
allocate (SYSCONNEC(NBASE,NBASE))
allocate (EtatSYS(NBASE,NBASE))    !*** type d'interaction
allocate (NBROT(NBASE,NBASE,4))
allocate (GDROT(NBASE,NBASE))
allocate (INVCONEC(NBASE,NBASE))
allocate (PartVis(NBASE))
allocate (ProjVis(NBASE))
allocate (Projcoin(NBASE))
allocate (SYSEG(NBASE))
allocate (LOMA(NBASE))
allocate (JSLOMA(NBASE))
allocate (FAC1DEP(NBASE,NBASERED))
allocate (FAC2DEP(NBASE,NBASERED,NBASERED))
allocate (MODDEPLA(NBASE,NBASERED,NBASERED))
allocate (VECNORLIN(3,NBASE))
allocate (VECNORDEP(3,NBASE))
allocate (VECNORNOR(3,NBASE))
allocate (NORMLIN(NBASE))
allocate (NORMDEP(NBASE))
allocate (tau0(nbase))
allocate (Bf(nbase))
allocate (RAUSYS(NTSG))
allocate (RAUSYS_precipite(NTSG))
allocate (NJONCSYS(NTSG))
allocate (RAUSYSjonc(NTSG))
allocate (CranSys(NTSG))
allocate (AIRESYS(NTSG))
allocate (AIRESYSInst(NTSG))
allocate (GAMSYS(NTSG))
allocate (AireVisSys(NTSG))
allocate (AireCoinSys(NTSG))
allocate (AireMixtSys(NTSG))
allocate (TrAppSys(NTSG))
allocate (TrAppSysInst(NTSG))
allocate (TrIntSys(NTSG))
allocate (TrIntSysInst(NTSG))
allocate (SchmidSys(NTSG))
allocate (Signe_TAUAPP(NTSG))
allocate (TensRotSys(NTSG,3,3))
allocate (TensEpsSys(NTSG,3,3))
allocate (TensDistSys(NTSG,3,3))
allocate (GAMMADOTSYS(NTSG))
allocate (Numero_loi(NBASE))
allocate (Inn(NBASE))
allocate (EPSOLD(NstatControl))
allocate (GAMMAOLDSYS(NTSG,NstatControl))
allocate(SysDev(NTSG,NbSysDev))

end subroutine Allocation

!************************************************************************
! Ce sous-programme ecrit dans le fichier (BVD.cfc, cc,cs, etc.)
! les vecteurs de la base de simulation
!  Ghiath Monnet -  07/02/01
!************************************************************************
subroutine READ_basedevecteurs

implicit none

integer (kind=DPI) :: i,j
character(len=40)  :: cristallo

cristallo = '../in/CRISTALLOGRAPHIE/BVD.'//crystal_structure(1:3)

open (15,FILE=cristallo,STATUS='OLD')  ! fichier de sortie
! d'abord on ecrit la liste brut des vecteur de discretisation
! pour l'utilisation d'autre programmes
read (15,*) nbase       ! unique de programme : lu par camera ....

do j = 1,nbase
  read (15,*) i, bvecnor(1:3,j), bveclin(1:3,j),bvecdep(1:3,j)
  if (j /= i) print *, "stop probleme de lecture de BVD"
enddo

close (15)

end subroutine READ_basedevecteurs

!############################################################################
!#  Ce sous progamme effectue une rotation+homothetie de tout les vecteurs  #
!#  n, l et d pour les materiaux CC                                         #
!#  Sylvain Queyreau -24janv05-                                             #
!############################################################################
subroutine rotation_BVD

!declaration des variables
integer(DPI) :: i
!integer, dimension (3,nbase) :: bvecnortourn, bveclintourn, bvecdeptourn

!rotation sur tous les bvecnor
do i = 1 , nbase
    Bvecnor(1:3,I) = MatMul(Matrice_rotBVD,Bvecnor(1:3,I))
    Bveclin(1:3,I) = MatMul(Matrice_rotBVD,Bveclin(1:3,I))
    Bvecdep(1:3,I) = MatMul(Matrice_rotBVD,Bvecdep(1:3,I))
enddo
! rotation of the tensile axis
Z =  MatMul(Matrice_rotBVD,Z)
Z(1:3) = Z(1:3)/DSQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))

end subroutine rotation_BVD

!#############################################################################
!############           initialisation des variables globales       ##########
!#############################################################################
subroutine initialiser

implicit none

!**** PARANOIAC SETTING ****
RAUDIS = ZERO
RAUDIS_matrice = ZERO
RAUDIS_precipite = ZERO
RAUDMO_I = ZERO
VITMOY = ZERO
EPSDOT = ZERO
kk = 0
npar = izero
accutime = zero
airevis = zero
airecoin = zero
airemixte = zero
OLDNTGD1 = izero
OLDNTGD2 = izero
AireSys(1:NTSG) = zero
AireSysInst(1:NTSG) = zero
AireVisSys(1:NTSG)  = zero
AireCoinSys(1:NTSG) = zero
AireMixtSys(1:NTSG) = zero
RAUSYSjonc(1:NTSG) = zero
TrAppSys(1:NTSG) = zero
TrAppSysInst(1:NTSG) = zero
TrIntSys(1:NTSG) = zero
TrIntSysInst(1:NTSG) = zero
TrApp = zero
TrAppInst = zero
TrInt = zero
TrIntInst = zero
IDEP(1:NSEGMAX)=0
VL_CS(1:NSEGMAX) = 0
IBoite(1:NSEGMAX) = IZERO
QUIDEP(1:NSEGMAX)=0
DAIRE(1:NSEGMAX)=zero
OUT(1:NSEGMAX)=.false.
TensRot(1:3,1:3)=zero
TensDist(1:3,1:3)=zero
EPSDOT = ZERO           !*** Initialisations ...
Solli_sys(1:NTSG_max) = UN
RAUSYS(:)               = zero
RAUSYS_precipite(:)     = zero
CRANSYS(:)              = zero
GAMSYS(:)               =zero
GAMMADOTSYS(:)          =0.d0
kk_debut_jonc = 0  ! iteration a la quelle commence la premiere jonction
kkjonc = 0         ! nb d'iteration ou la jonction existe
Sigma = zero
MICROBOUCLE = izero
Nboites = IUN
Nbplan =IZERO
NbplanMax =IZERO
NbFreePlan = IZERO
NbplanDom = IZERO
NbCvxDom = IZERO
InclExcl=.false.            ! A default value is given, but the real value is set in module microstructure
DesorientGrain=.false.      ! The default solution is no disorientation
Nb_GroSeg = izero
GroSeg(1:Nsegmax/idix) = izero
seg(:)%zerotl = .False.

end subroutine initialiser

!#########################################################################
!# Procedure de tabulation des deplacements et autres tableau            #
!############################################################## 02/2001 ##
subroutine Tabulation

implicit none

integer(kind=DPI)  :: v1(3),v2(3),v3(3),d1(3),d2(3),Jindice,alpha,beta,gamma
integer(kind=DPI)  :: i,j,jj,ii,k,qq,bi(3), bj(3),li(3),lj(3),bjonc(3),facboite
integer(kind=DPI)  :: buf,buf1,indj,indk,systeme(3,3)
integer(kind=DPI)  :: k2,br,ni(3),nj(3),Ldim(3),Ldimtemp(3)
integer(kind=DPI)  :: calcsys, edgedir
integer(kind=DPI)  :: totoBUGv1(3),totoBUGv2(3),coupleIJ(2,2),coupleIJc(2,2)
real(kind=DP)      :: x,rv(3),rv1(3),dimensions(3),dimtemp(3),inconnues(3)
real(kind=DP)      :: nrv,nrv1
logical            :: CtiUneVis
character          :: calculboite

real(kind=DP), dimension(3,3) ::   MatIdent = RESHAPE((/1.,0.,0.,0.,1.,0.,0.,0.,1./),SHAPE=(/3,3/))

br = nbasered      ! pour faciliter l'ecriture

!*** TABULATION ANTI SPIRALE EN PT D'ANCRAGE :
!     RENVOIE LE CARACTERE REMPLACANT LES DEUX INITIAUX
!      (Juste la traduction en indice de la base reduite d'une somme vectorielle)
SC(:,:)=0
CoefSC(:,:)=0
!allocate (SC(nbase,nbase))
!allocate (CoefSC(nbase,nbase))
qq = 0
if (HCP) then
   facteurvis (:)  = (/TROIS,TROIS,TROIS,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/Quatre,un,un,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/un,half,half,half,zero,zero,zero,zero,zero,zero/)
elseif (BCC) then
   facteurvis (:)  = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/deux,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/0.334D0,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
elseif (ORT) then
   facteurvis (:)  = (/6.0D0,cinq,6.0D0,cinq,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/DIX,dix,cinq,6.0D0,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/Un,Un,Un,Un,Un,Un,zero,zero,zero,zero/)
elseif(CFC) then
   facteurvis (:)  = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/half,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
elseif(CS) then
   facteurvis (:)  = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
elseif(MGO) then
   facteurvis (:)  = (/un,un,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/un,un,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/un,half,zero,zero,zero,zero,zero,zero,zero,zero/)
elseif(DC) then
   facteurvis (:)  = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/un,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/half,zero,zero,zero,zero,zero,zero,zero,zero,zero/)
endif
do i = 1,nb_slip_types
    do j = 1, slip(i)%nsystemes
        k = qq*nbasered
!        print *, i,j,k
        if(facteurmixt(i) == 1) then
           SC(k + 1,k + 3)= k + 2 !GENESME
           CoefSC(k + 1,k + 3)= 1
           SC(k + 1,k + 4)= k + 3
           CoefSC(k + 1,k + 4)= 1
           SC(k + 1,k + 6)= k + 7
           CoefSC(k + 1,k + 6)= 1
           SC(k + 1,k + 7)= k + 8 !GENESME
           CoefSC(k + 1,k + 7)= 1

           SC(k + 3,k + 1)= k + 2 !GENESME
           CoefSC(k + 3,k + 1)= 1
           SC(k + 3,k + 5)= k + 4 !GENESME
           CoefSC(k + 3,k + 5)= 1
           SC(k + 3,k + 6)= k + 5
           CoefSC(k + 3,k + 6)= 1
           SC(k + 3,k + 8)= k + 1
           CoefSC(k + 3,k + 8)= 1

           SC(k + 5,k + 2)= k + 3
           CoefSC(k + 5,k + 2)= 1
           SC(k + 5,k + 3)= k + 4 !GENESME
           CoefSC(k + 5,k + 3)= 1
           SC(k + 5,k + 7)= k + 6 !GENESME
           CoefSC(k + 5,k + 7)= 1
           SC(k + 5,k + 8)= k + 7
           CoefSC(k + 5,k + 8)= 1

           SC(k + 7,k + 1)= k + 8 !GENESME
           CoefSC(k + 7,k + 1)= 1
           SC(k + 7,k + 2)= k + 1
           CoefSC(k + 7,k + 2)= 1
           SC(k + 7,k + 4)= k + 5 ! ATTENTION : il faut une base de type PM VCM1M2
           CoefSC(k + 7,k + 4)= 1
           SC(k + 7,k + 5)= k + 6 !GENESME
           CoefSC(k + 7,k + 5)= 1

        elseif (facteurmixt(i) == 0.5) then

           SC(k + 1,k + 3)= k + 2 !GENESME
           CoefSC(k + 1,k + 3)= 2 !GENESME
           SC(k + 1,k + 4)= k + 2
           CoefSC(k + 1,k + 4)= 1
           SC(k + 1,k + 6)= k + 8
           CoefSC(k + 1,k + 6)= 1
           SC(k + 1,k + 7)=  k + 8 !GENESME
           CoefSC(k + 1,k + 7)= 2 !GENESME

           SC(k + 3,k + 1)= k + 2 !GENESME
           CoefSC(k + 3,k + 1)= 2 !GENESME
           SC(k + 3,k + 5)=  k + 4 !GENESME
           CoefSC(k + 3,k + 5)= 2
           SC(k + 3,k + 6)=  k + 4
           CoefSC(k + 3,k + 6)= 1
           SC(k + 3,k + 8)= k +  2
           CoefSC(k + 3,k + 8) = 1

           SC(k + 5,k + 2)= k +  4
           CoefSC(k + 5,k + 2)= 1
           SC(k + 5,k + 3)= k +  4 !GENESME
           CoefSC(k + 5,k + 3)= 2
           SC(k + 5,k + 7)= k +  6 !GENESME
           CoefSC(k + 5,k + 7)= 2
           SC(k + 5,k + 8)= k +  6
           CoefSC(k + 5,k + 8)= 1

           SC(k + 7,k + 1)= k +  8 !GENESME
           CoefSC(k + 7,k + 1)= 2
           SC(k + 7,k + 2)= k +  8
           CoefSC(k + 7,k + 2)= 1
           SC(k + 7,k + 4)= k +  6
           CoefSC(k + 7,k + 4)= 1
           SC(k + 7,k + 5)= k +  6 !GENESME
           CoefSC(k + 7,k + 5)= 2
        elseif(facteurmixt(i) == 0.334D0) then
           SC(k + 1,k + 4)= k + 2
           CoefSC(k + 1,k + 4)= 1
           SC(k + 1,k + 6)= k + 8
           CoefSC(k + 1,k + 6)= 1
           SC(k + 1,k + 7)= k + 8 !GENESME
           CoefSC(k + 1,k + 7)= 3

           SC(k + 2,k + 5)= k + 4 !GENESME
           CoefSC(k + 2,k + 5)= 1
           SC(k + 2,k + 7)= k + 8 !GENESME
           CoefSC(k + 2,k + 7)= 2
           SC(k + 2,k + 8)= k + 1
           CoefSC(k + 2,k + 8)= 1

           SC(k + 3,k + 5)= k + 4 !GENESME
           CoefSC(k + 3,k + 5)= 3
           SC(k + 3,k + 6)= k + 4 !GENESME
           CoefSC(k + 3,k + 6)= 2

           SC(k + 4,k + 1)= k + 2
           CoefSC(k + 4,k + 1)= 1
           SC(k + 4,k + 6)= k + 5 !GENESME
           CoefSC(k + 4,k + 6)= 1

        else
           stop  " facteur mixt non pris en compte pour la tabulation de SC"
        endif
        qq = qq + 1
    enddo
enddo


! verification des tableau Sc et coefSC
do i =1,nbase
    do j =1,nbase
        if (SC(i,j) /= 0) then
           v1(:) = bveclin(:,I) + bveclin(:,j)
           v2(:) = coefSC(i,j) * bveclin(:,SC(i,j))
           if (etat(v1,v2) /= 3) then
              print *, "i,j ", i ,j
              print *, "Sc, Coef ", sc(i,j),coefsc(i,j)
              print *, v1
              print *, v2
              stop "Stop !!!!!"
           endif
        endif
    enddo
enddo

!*** Construction du tableau des vecteurs normalises
do I=1,NBASE
    if(modulo(bveclin(1,i),ideux) /= izero .or. &
         modulo(bveclin(2,i),ideux)  /= izero .or. &
         modulo(bveclin(3,i),ideux)  /= izero)  then
       print *, "force : vecteur BVD non pair ",I," Interdit pour calcul de force"
    endif

    if (bveclin(1,i) /= 0) then
       Inn(i) = 1
    elseif (bveclin(2,i) /= 0) then
       Inn(i) = 2
    elseif (bveclin(3,i) /= 0) then
       Inn(i) = 3
    else
       stop  " Base de vecteur erronee"
    endif

    NormLin(I) = norivect(BVECLIN(1:3,I))
    VecNorLin(1:3,I) = normaivec(BVECLIN(1:3,I))
    NormDep(I) = norivect(BVECDEP(1:3,I))
    VecNorDep(1:3,I) = normaivec(BVECDEP(1:3,I))
    VecNorNor(1:3,I) = normaivec(BVECNOR(1:3,I))
enddo

!*** Construction des tableaux utiles a la manipulation de la base

do I=1,NBASERED

    INVSI(I)=(I+(NBASERED/2))

    if (INVSI(I).gt.NBASERED) INVSI(I)=INVSI(I)-NBASERED

enddo

do I=1,NBASE
    do J=1,NBASERED
        ASSOC(I,J) = J + ((I-IUN)/NBASERED)*NBASERED ! vecteur J associe
    enddo

    DEASSOC(I) = I-ASSOC(I,IUN)+1    ! type "etendu" (1 a 8)
    TYSEG(I)   = I-( ((I-IUN)/(NBASERED/IDEUX))*(NBASERED/IDEUX) )   ! type "light"  (1 a 4)
    totoBUGv1(1:3)=bveclin(1:3,I)
    totoBUGv2(1:3)=bveclin(1:3,assoc(I,1))
    SISEG(I)   = sicose2v(totoBUGv2,totoBUGv2)
    if (siseg(I).eq.0) SISEG(I)=siseg(I-1)
!    BISYSEG(I) = ((I-IUN)/(IDEUX*NBASERED))+IUN    ! bisysteme de glissement
    SYSEG(I)   = ((I-IUN)/(NBASERED))+IUN      !   systeme de glissement
enddo

! Calculation of the applied stress tensor
if (uniaxiale) then
  call  ORIENTDIEG(un,Z,SIGAPP,.TRUE.)         !The applied stress tensor
else
  call  TensIncr(un,tensapp,SIGAPP,.TRUE.)     !The applied stress tensor
endif

! Calculation of the initial Schmid factor for the slip systems. Note that
! SchmidSys may be smoothly modified during the simulation due to plastic rotation
call GlidCompPK(sigapp, nbase, ntsg, VecNorLin, VecNorNor, SchmidSys)

! Calculation of the elementary rotation tensor for each slip system, such tensor is needed
! to compute the total rigid body rotation due to each slip system contriution
do ii=1,NTSG
    I=(ii-IUN)*nbasered+IUN    ! Index of the first vector of each slip system
    rv(:) =real(bveclin(:,I),DP)
    rv1(:)=real(bvecnor(:,I),DP)
    nrv =norvect(rv)
    nrv1=norvect(rv1)
    do j=1,3
        do k=1,3
            ! The system plastic rotation tensor
            TensRotSys(ii,j,k)=half*(real(bveclin(j,I),DP)*real(bvecnor(k,I),DP)-    &
                                     real(bveclin(k,I),DP)*real(bvecnor(j,I),DP))    &
                              /(nrv*nrv1)
            ! The systems plastic deformation tensor
            TensEpsSys(ii,j,k)=half*(real(bveclin(j,I),DP)*real(bvecnor(k,I),DP)+    &
                                     real(bveclin(k,I),DP)*real(bvecnor(j,I),DP))    &
                              /(nrv*nrv1)
            ! The systems plastic distortion tensor
            Z=normavec(Z)     ! In the above equation Z must be a unit vector
            TensDistSys(ii,j,k)=(real(bveclin(j,I),DP)*real(bvecnor(k,I),DP))/(nrv*nrv1)- ( MatIdent(j,k)*                       &
                                (real(bveclin(1,I),DP)*Z(1)+real(bveclin(2,I),DP)*Z(2)+real(bveclin(3,I),DP)*Z(3))*              &
                                (real(bvecnor(1,I),DP)*Z(1)+real(bvecnor(2,I),DP)*Z(2)+real(bvecnor(3,I),DP)*Z(3))/(nrv*nrv1))
        enddo
    enddo
enddo

! Definition of slip system particular projection directions
do I=1,NBASE

    totoBUGv1(1:3) = bveclin(1:3,i)

    ! The screw direction common to all collinear slip systems
    ! Partvis is the normalized positive screw length of the Nbase vectors
    ! Projvis is the not normalized screw length of the Nbase vectors
    totoBUGv2(1:3) = bveclin(1:3,assoc(I,1))
    PartVis(I) = Real(abs(dot_product(totoBUGv1(1:3),totoBUGv2(1:3))),DP)
    PartVis(I) = PartVis(I)/normlin(i)/normlin(assoc(I,1))
    ProjVis(I) = Real(dot_product(totoBUGv1(1:3),totoBUGv2(1:3)),DP)

    ! The same edge direction must be used for all collinear slip systems (systems with the same b)
    ! Projcoin is the not normalized egde length of the Nbase vectors
    calcsys = syseg(I) - modulo((syseg(I)-iun),(NBSysDev+iun))
    edgedir = assoc(((calcsys-iun)*nbasered+iun),3)
    totoBUGv2(1:3) = bveclin(1:3,edgedir)
    Projcoin(I) = Real(dot_product(totoBUGv1(1:3),totoBUGv2(1:3)),DP)

    do J = 1, 2, 1
      indj = i+((j-1)*2-1)
      if (indj.gt.assoc(i,NBASERED)) indj=indj-NBASERED
      if (indj.lt.assoc(i,1)) indj=indj+NBASERED
      CONEC(I,J) = indj
      DECONEC(I,indj)=J
    enddo

    INVCONEC(I,CONEC(I,1)) = CONEC(I,2)
    INVCONEC(I,CONEC(I,2)) = CONEC(I,1)

enddo

do I=1,NBASE
    do J=1,NBASE
        nbrot(I,J,:) = nsegmax
        if(syseg(i) /=syseg(j) ) CYCLE
        if (J.eq.conec(I,1).or.J.eq.conec(I,2)) then
           NBROT(I,J,1)=0 !*** Segment directement connectable
           NBROT(I,J,2)=0
           NBROT(I,J,3)=0
           NBROT(I,J,4)=0
        else
           if (I.eq.J) then
              NBROT(I,J,1)=1 !*** 1 rotule
              NBROT(I,J,2)=conec(I,1)
              NBROT(I,J,3)=conec(I,2)
              NBROT(I,J,4)=0
           else
              if (conec(I,2).eq.conec(J,1)) then
                 NBROT(I,J,1)=1 !*** 1 rotule
                 NBROT(I,J,2)=conec(I,2)
                 NBROT(I,J,3)=conec(I,2)
                 NBROT(I,J,4)=0
              else
                 if (conec(I,1).eq.conec(J,2)) then
                    NBROT(I,J,1)=1 !*** 1 rotule
                    NBROT(I,J,2)=conec(I,1)
                    NBROT(I,J,3)=conec(I,1)
                    NBROT(I,J,4)=0
                 else
                    if (conec(conec(I,1),1).eq.conec(j,2).and.(nbasered.ge.8)) then
                       NBROT(I,J,1)=2 !*** 2 rotule
                       NBROT(I,J,2)=conec(I,1)
                       NBROT(I,J,3)=conec(J,2)
                       NBROT(I,J,4)=0
                    else
                       if (conec(conec(I,2),2).eq.conec(j,1).and.(nbasered.ge.8)) then
                          NBROT(I,J,1)=2 !*** 2 rotule
                          NBROT(I,J,2)=conec(I,2)
                          NBROT(I,J,3)=conec(J,1)
                          NBROT(I,J,4)=0
                       else
                          NBROT(I,J,1)=(nbasered/2)-1 !*** 3 rotules : recouvrement !!!
                          if(I < J) then
                             NBROT(I,J,2) = I + 1
                             NBROT(I,J,3) = I + 2
                             NBROT(I,J,4) = I + 3
                          else
                             if(I/nbasered == (i-1)/nbasered) then
                                NBROT(I,J,2) = I + 1
                             else
                                NBROT(I,J,2) = I + 1 - nbasered
                             endif
                             if((I+1)/nbasered == (i-1)/nbasered) then
                                NBROT(I,J,3) = I + 2
                             else
                                NBROT(I,J,3) = I + 2 - nbasered
                             endif
                             if((I+2)/nbasered == (i-1)/nbasered) then
                                NBROT(I,J,4) = I + 3
                             else
                                NBROT(I,J,4) = I + 3 - nbasered
                             endif
                          endif
                       endif
                    endif
                 endif
              endif
           endif
        endif
    enddo
enddo


!**************** TABULATION POUR LE GLISSEMENT DEVIE ************
EtatSYS(:,:) = IDIX
do I=1,NBASE
    do J=1,NBASE
        bi =  bveclin(:,assoc(i,1))       ! vecteur de burgers pour systeme i
        bj =  bveclin(:,assoc(j,1))       ! vecteur de burgers pour systeme j
        ni = bvecnor(:,i)       ! normal au plan i
        nj = bvecnor(:,j)
        if(abs(etat(ni,nj)) > 1 .and. abs(etat(bi,bj)) > 1) then
           EtatSys(I,J) = izero  ! meme systeme
        elseif(abs(etat(ni,nj)) > 1 .and. abs(etat(bi,bj)) < 2) then
           EtatSys(I,J) = IUN ! systemes coplanairs
        elseif(abs(etat(ni,nj)) < 2 .and. abs(etat(bi,bj)) > 1) then
           EtatSys(I,J) = IDEUX ! systemes coliniaires
        else
           EtatSys(I,J) = ITROIS ! systemes coplanairs
        endif

        GDROT(I,J)=0

!*** appel du type gdrot(vect de sysi,vect de sysj'=dev(sysj)=sysi)
        if (ASSOC(I,1).eq.ASSOC(J,1)) then
           if (tyseg(i).ne.1.and.tyseg(j).ne.1) then
              CtiUneVis=.false.
              if ((NBROT(I,J,1).ne.(nbasered/2-1)).and.(NBROT(I,J,1).ne.0)) then

                 do Jindice=1,NBROT(I,J,1),1
!*** recherche d'une vis ds les rotules standards
                     if (tyseg(NBROT(I,J,Jindice+1)).eq.1) CtiUneVis=.true.
                 enddo
                 if (CtiUneVis) then
!*** On a deja une vis => pas de pb
                    GDROT(I,J)=0
                 else
!*** Autrement on choisit la vis la plus adaptee
                    totoBUGv1(1:3)=bveclin(1:3,i)+bveclin(1:3,j)
                    totoBUGv2(1:3)=bveclin(1:3,assoc(i,1))
                    if(sicose2v(totoBUGv1(1:3),&
                    totoBUGv2(1:3)).gt.0) then
                       GDROT(I,J)=ASSOC(I,1)
                    else
                       GDROT(I,J)=ASSOC(I,(nbasered/2)+1)
                    endif
                 endif

              else
                 GDROT(I,J)=ASSOC(I,1)
              endif !*** PAS VRAIMENT DU RECOUVREMENT SI PAS VIS PUISQUE PLAN et PLAN DEV
           endif !if (tyseg(i).ne.1.and.tyseg(j).ne.1) then
!*** On a deja une vis => pas de pb
        endif !if (ASSOC(I,1).ne.ASSOC(J,1)) then
!*** Sys differents
! write (*,*) 'gdrot',i,j,GDROT(I,J)
    enddo
enddo

!*** Ronan MADEC Le 06/02/01 *************************************

!*** Tabulation des facteurs de la formule generale traitant le deplacement d'un
!*  segment connaissant son vecteur ligne et ceux de ses premiers voisins ainsi
!*  que son vecteur deplacement.
!* indice 1 : voisin en o
!* indice 2 : le segemnt i
!* indice 3 : voisin en e
!* fac1dep (i,j) renvoie la variation de la norme du voisin j quand i avance de 1 pas
!* fac2dep (i,j) renvoie la variation de la norme du voisin I quand i avance de 1 pas
!* dL1 = (v1.d2) / abs(v1.d2) n   ou   n = v dt / d2
!*     = fac1dep(v2,v1) n
!*fac1dep(v2,v1) = (v1.d2) / abs(v1.d2)
!* dO2 = dL1 v1 = fac1dep n v1
!*
!* dL3 = -(v3.d2) / abs(v3.d2) n
!*     = - fac1dep(v2,v3) n
!*
!* dO3 = - dL3 v3 = fac1dep n v3
!*
!* dv2 = dO3-dO2
!*     = - (- fac1dep(v2,v3) n v3 + fac1dep(v2,v1) n v1)
!*     = - n (- fac1dep(v2,v3) v3 + fac1dep(v2,v1) v1)
!*     = n fac2dep(v2,v1,v3) v2

!*** boucle sur le segment a deplacer

fac1dep(1:NBASE,1:NBASERED) = real(nsegmax,DP)
fac2dep(1:NBASE,1:NBASERED,1:NBASERED) = real(nsegmax,DP)
MODDEPLA(1:NBASE,1:NBASERED,1:NBASERED) = izero
modep_max = izero
do i=1,NBASE,1
    v2(1:3) = bveclin(1:3,i)
    d2(1:3) = bvecdep(1:3,i)

!*** fac1dep
!*** boucle sur son "VOISO" ou son "VOISE"
    do j=1,NBASERED,1
!*** indice j
        indj = ASSOC(I,J)
        v1(1:3) = bveclin(1:3,indj)
! procedure de Ronan
        d1(1:3) = bvecdep(1:3,indj)
        buf = dot_product(v1,d2)
        buf1 = int(abs(buf),DPI)

        if (buf1.ne.0) then
! procedure Ghiath
           fac1dep(I,J) = real(inorivect(d2(:)),DP)/real(buf,DP)  ! = depsurprojec
!*** tabulation de la courbure
           if (indj.eq.conec(i,1).or.indj.eq.conec(i,2)) then
              if (buf.gt.zero) then
                 courbuO(i,2) = indj    ! pour le deplacement negatif de i
                 courbuE(i,1) = indj
              else
                 courbuO(i,1) = indj ! pour le deplacement positif de i
                 courbuE(i,2) = indj
              endif
           endif
           x = abs(1.0/fac1dep(I,J))
           if(x/= UN .and. x/= DEUX .and. x/= TROIS .and. x/= QUATRE) then
              print*,  " i =",i," v2  =",v2(:)
              print*,  " d2 =", d2(:)
              print*,  "  v1=", v1(:)
              print*,  "  fac1dep=",fac1dep (i,j)
              stop 'Stop !!!!'
           endif
        endif
    enddo


!*** boucle sur son "VOISO"
    do j=1,NBASERED
!*** indice j
! Ghiath : je multiplie tous les vecteur par 8 pour s'assurer
! qu'il n'y pas de 0.5000 qui train
        indj = ASSOC(I,J)
        v1(1:3) = bveclin(1:3,indj)
!  d1(1:3) = bvecdep(1:3,indj)
!*** boucle sur son "VOISE"
        do k=1,NBASERED
!*** indice k
            indk = ASSOC(I,K)
            v3(1:3) = bveclin(1:3,indk)
!   d3(1:3) = bvecdep(1:3,indk)

!*** fac2dep
            if(deassoc(i).ne.j.and.deassoc(i).ne. k.and.deassoc(i).ne.invsi(j).and. &
            deassoc(i).ne.invsi(k)) then
               rv(1:3) = fac1dep(i,k) * Real(v3(1:3),DP) - fac1dep(i,j) * real(v1(1:3),DP)
! les variations de longueur de vo et ve pour dep = 1 de I, doit garder I dans la
! meme direction
               x = rv(Inn(I)) / real(v2(Inn(I)),DP)
               rv1(:) = real(d2(:),DP)
               if(.not. egalite(dot_product(rv,rv1),ZERO)) then
                  print *, " j,i,k ", j,i,k
                  print *, "v1 =", v1,fac1dep(i,j)
                  print *, "v2 =", v2
                  print *, "v3 =", v3,fac1dep(i,k)
                  print *, " I apres dep=",rv
                  print *, "prod =", dot_product(rv,rv1)
                  print *, " fac2dep =", x
                  print*, " fac1dep de vo et ve incompatible"
                  stop
               endif
! verification
               fac2dep(i,j,k) = x

! ghiath : afin de generaliser la procedure mixt/mixt, on fait le test sur les segmnt i
! et ces voisins. L'idee est de garantir un longueur entier des segments en modulant
! le deplacement : absdep doit etre multiple de MODEP, qui est une variable integer
! calculee en fonction de caractersitiques geometriques des segmnts, ex cas des mixt/mixt
! dans les cfc MODEP=2

! initialisation
               jj = 1
               qq = 1
               k2 = 1
               x = dabs(fac1dep(i,j))  ! compatibilite avec le VO
               if (.not. entier(x)) then
                  if (x < 1) then
                     jj = nint(ppem(x)/x)
                  else
                     stop "INIT:ERREUR 1: il faut changer la procedure de MODDEPLA"
                  endif
               endif

               x = dabs(fac1dep(i,k))    ! compatibilite avec le VE
               if (.not. entier(x)) then
                  if (x < 1) then
                     qq = nint(ppem(x)/x)
                  else
                     stop "INIT:ERREUR 2: il faut changer la procedure de MODDEPLA"
                  endif
               endif

               x = dabs(fac2dep(i,j,k))    ! compatibilite avec I
               if (.not. entier(x)) then
                  k2 =  nint(ppem(x)/x)
               endif
               moddepla(i,j,k) = ppcm3(jj,int(qq,DPI),k2)
               if(moddepla(i,j,k) > modep_max) modep_max = moddepla(i,j,k)
            endif

!            write (*,*) j,i,k," moddepla ",moddepla(i,j,k)
!           write (*,*)dabs(fac1dep(i,j)) ,dabs(fac1dep(i,k)),dabs(fac2dep(i,j,k))
        enddo
    enddo
enddo

write(*,*) " "
write(*,*) "Modep_Max of this Discretisation Vector Base : ", modep_max

!###########################################
!##### tabulation des axes de jonction #####
!###########################################
axejoncj(:,:) = IZERO         !   INITIALISATION
axejonci(:,:) = IZERO         !   INITIALISATION
SYSOBST(:,:) = 999        !   INITIALISATION
SYSCONNEC(:,:) = 999   !POUR CONNECTER LES SEGMENTS ENTRE EUX

! allocation of tables necessary for cross slip
SegDev(1:NBASE,1:NBsysdev) = Izero ! to get corresponding cros slip segments in all cross slip systems
SegDevIJ(1:NBASE,1:NBase) = Izero  ! to get corresponding cros slip segment in the slip system of J
SysDev(1:NTSG,1:NbSysDev)  = Izero ! to get the all cros slip systems for the given system

do ii=1,NTSG
    gamma = izero
    do jj = 1 , NTSG
        if(ii == jj ) CYCLE
        I = (II-IUN) * nbasered + IUN
        J = (JJ-IUN) * nbasered + IUN    ! Index of the first vector in the slip system JJ

        ni = bvecnor(:,i)       ! normal au plan i
        nj = bvecnor(:,j)

        bi =  bveclin(:,i)       ! vecteur de burgers pour systeme i
        bj =  bveclin(:,j)       ! vecteur de burgers pour systeme j
!*** TABULATION DES TYPES D'OBSTACLES ET DE CONNECTIVITES
        if (int(abs(etat(ni,nj)),DPI) < 2 .and. int(abs(etat(bi,bj)),DPI) >= 2) then   ! system colineaire
           if (etat(bi,bj) /=  3) then
              print *, " Fatal error:INIT: systems collinear and first vectors are different"
              print *, " Not supported option"
              stop
           endif
           ! system colineaire
           if(gamma == nbsysdev) stop " Fatal error:INIT: detected NbSysDev > given Nbsysdev"
           gamma = gamma + IUN
           SysDev(II,gamma)  = JJ
!            print *, "SysDev(II,gamma) ",II,gamma,SysDev(II,gamma)
           do k =  0 , nbasered - iun
               SegDev(I+K,gamma) = J + K
!            print *, "SegDev(I+K,gamma) ",I+K, gamma,SegDev(I+K,gamma)
           enddo
        endif
    enddo
enddo

write(*,*) "Detected number of cross slip systems : ", gamma

do i=1,nbase,1
    do j=1,nbase,1
        coupleIJ(:,:) = nsegmax
        coupleIJc(:,:) = nsegmax
        li = bveclin(:,i)       ! vecteur de ligne i
        lj = bveclin(:,j)       ! vecteur de ligne j
        ni = bvecnor(:,i)       ! normal au plan i
        nj = bvecnor(:,j)
        bi =  bveclin(:,assoc(i,1))       ! vecteur de burgers pour systeme i
        bj =  bveclin(:,assoc(j,1))       ! vecteur de burgers pour systeme j
        bjonc = bi(:) + bj(:)  ! vecteur de jonction si les vecteurs de binom de jonction tont PARALLELs

!*** TABULATION DES TYPES D'OBSTACLES ET DE CONNECTIVITES

        if (int(abs(etat(ni,nj)),DPI) >= 2) then   ! meme plan de gliss
           if (int(abs(etat(bi,bj)),DPI) >= 2) then     !meme systeme de gliss
              sysconnec(i,j) = 1
              sysOBST(i,j) = 1
              segdevIJ(i,J) = assoc(J,1) + deassoc(i) -1
           else
              sysconnec(i,j) = 0     ! systems differents mais coplanaires
              sysOBST(i,j) = 0
           endif
        else             ! plan differents
           if (int(abs(etat(bi,bj)),DPI) >= 2) then       ! meme b
              segdevIJ(i,J) = assoc(J,1) + deassoc(i) -1
              sysconnec(i,j) = 2     ! vecteur de BUrgers : configuration de GD
! nouveau : 10/10 : devseg(i,J) renvoie l'indice du VL correspondant a I dans
! on maintient la forme initial de DEVSEG pour pouvoir continuer a fonctioner pour les CFC
              if ((int(abs(etat(li,bi)),DPI) >= 2) .and. (int(abs(etat(lj,bj)),DPI) >= 2)) then
                 sysOBST(i,j) = 1   ! les segments i ET J correspondent a des vis
              else
                 sysOBST(i,j) = 6
              endif
           else
              sysconnec(i,j) = 3     ! Burgers differents et non coplanaires
              sysOBST(i,j) = 6
           endif
        endif

! 1) d'abord : on ne travaille pas avec meme systemes de glissement
        if (Etatsys(i,j) == izero) cycle      ! ni,nj : meme system de glissement
! modification : dans les bcc, les jonction entre coplanaires sont cessile
! donc il faut kles prendre en compte
        if (Etatsys(i,j) == IUN .and. .not. BCC) cycle      ! meme plan : sauf BCC: jonction coplanaire cissile
! des vecteur de Burgers dans la base reduite
! Attention : on verifie plus tard le signe des vecteur lignes de la jonction pour i et pour j
! Cette modification est introduite pour faire un choix physique des vecteurs de jonction qui doit
! conduire a une jonction mixte de Hirth : (bi + bj) n'est pas perpendiculaire a l'axe de jonction
!
! determination des couples de vecteur parallels en meme temps dans les deux systems(4vecteurs)
        surII : do ii = assoc(i,1),assoc(i,br)
            surJJ : do jj = assoc(j,1), assoc(j,br)

! On ne garde que les couples parallels, forcement parallels aussi a l"axe d'intersection si ni et nj sont differents
! ce choix suppose implicitement que bi.bj < 0.
                if (etat(bveclin(:,ii),bveclin(:,jj)) > 1) then
! Reste les 4 couples possibles
!(coupleIJ(1,1) == indice du premier vecteur dans sys I parallel a un vecteur de sys J
!(coupleIJ(1,2) == indice de vecteur dans sys J parallel au vecteur coupleIJ(1,1)
!(coupleIJ(2,1) == indice de 2eme vecteur dans sys I parallel a un vecteur de sys J
!(coupleIJ(2,2) == indice de vecteur dans sys J parallel au vecteur coupleIJ(2,1)
! dans le cas des sys copla dans les BCC , il y a deux autres couples a considerer apres ceux de coupleIJ
                   if(coupleIJ(1,1) == nsegmax) then
                      coupleIJ(1,1) = ii
                   elseif(coupleIJ(2,1) == nsegmax) then
                      coupleIJ(2,1) = ii
                   elseif(coupleIJc(1,1) == nsegmax) then
                      coupleIJc(1,1) = ii
                   else
                      coupleIJc(2,1) = ii
                   endif
                   if(coupleIJ(1,2) == nsegmax) then
                      coupleIJ(1,2) = jj
                   elseif(coupleIJ(2,2) == nsegmax) then
                      coupleIJ(2,2) = jj
                   elseif(coupleIJc(1,2) == nsegmax) then
                      coupleIJc(1,2) = jj
                   else
                      coupleIJc(2,2) = jj
                   endif
                endif

            enddo surJJ
        enddo surII
        if(coupleIJ(1,1) == nsegmax .or. coupleIJ(1,2) == nsegmax .or. &
        coupleIJ(2,1) == nsegmax .or. coupleIJ(2,2) == nsegmax ) then
           if (.not. HCP) then
              stop  "INIT : erreur de determineation des axes d'intersection ..."
           else
              cycle
           endif
        endif
        if(etatsys(i,j) == IUN) then
           if(coupleIJc(1,1) == nsegmax .or. coupleIJc(1,2) == nsegmax .or.   &
           coupleIJc(2,1) == nsegmax .or. coupleIJc(2,2) == nsegmax )         &
              stop "INIT : erreur de determineation des axes d'intersection ..."
        endif
!********************************************************************************************
!*********** on selectionne les deux couples de jonction parmis les 4 de COUPLEIJ ***********
!*********** selon un critere energitique : exe Frank                             ***********
!********************************************************************************************
! Frank : choix des couple de jonction
! regle de Frank  (bj.bj) * (li.lj) < 0
! rappel : etat = -3,-2,-1,0,1,2,3 == aniparallels egaux, antiparallels, optu,,perpendiculaire, aigui,parralels,egaux
        alpha = etat(bi,bj)
        beta = etat(bjonc,bveclin(:,coupleIJ(1,1)))   ! relation entre (bj+bi) et l'axe d'intersection
        if(alpha < 0) then
! les couples deja selectionnes sont parallels : c-a-d correspondent au cas ou bj.bi < 0

        elseif(alpha > 0 .or. (alpha == izero .and. beta == Izero)) then
!  si l'angle entre bi bj est aigu, on choisit les couples antiparallels
! si bi.bj = 0 : il y a deux cas bjonc est perpendiculaire ou pas a l'axe de jonction
! A la place de Frank, on applique la regle energitique entre la Hirth coin ou mixte
! Resultat : il faut choisir les couples qui fournissent la configuration mixte
! Donc si bjonc.axe = 0 il faut prendre le couple antiparallels pour imposer bjonc = bi-bj
           ii = coupleIJ(1,2)
           if (deassoc(ii) < ICinq) then
              coupleIJ(1,2) = II + Iquatre
           else
              coupleIJ(1,2) = II - Iquatre
           endif
           ii = coupleIJ(2,2)
           if (deassoc(ii) < ICinq) then
              coupleIJ(2,2) = II + Iquatre
           else
              coupleIJ(2,2) = II - Iquatre
           endif
           if(etatsys(i,j) == IUN) then
              ii = coupleIJc(1,2)
              if (deassoc(ii) < ICinq) then
                 coupleIJc(1,2) = II + Iquatre
              else
                 coupleIJc(1,2) = II - Iquatre
              endif
              ii = coupleIJc(2,2)
              if (deassoc(ii) < ICinq) then
                 coupleIJc(2,2) = II + Iquatre
              else
                 coupleIJc(2,2) = II - Iquatre
              endif
           endif
        endif

!  verifixcation du critere de Frank bi.bj * li.lj < 0
!           print *,i,j,coupleIJ(1,1),coupleIJ(1,2)
!           print *,i,j,coupleIJ(2,1),coupleIJ(2,2)
        ii = coupleIJ(1,1) ; jj = coupleIJ(1,2)
        if(etat(bi,bj) * etat(bveclin(:,ii),bveclin(:,jj)) > izero ) &
           stop " erreur 1 INIT : mauvais chois coupleIJ "

        ii = coupleIJ(2,1) ; jj = coupleIJ(2,2)
        if(etat(bi,bj) * etat(bveclin(:,ii),bveclin(:,jj)) > izero ) &
           stop " erreur 2 INIT : mauvais chois coupleIJ "

        if(etatsys(i,j) == IUN) then
           ii = coupleIJc(1,1) ; jj = coupleIJc(1,2)
           if(etat(bi,bj) * etat(bveclin(:,ii),bveclin(:,jj)) > izero ) &
              stop " erreur 3 INIT : mauvais chois coupleIJ "
           ii = coupleIJc(2,1) ; jj = coupleIJc(2,2)
           if(etat(bi,bj) * etat(bveclin(:,ii),bveclin(:,jj)) > izero ) &
              stop " erreur 3 INIT : mauvais chois coupleIJ "
        endif
!  write(78,'("I =",I3," J=",I3," II1 =", I7," JJ1= ",I7," II2 = ", I7," JJ2= ",I7," II1 =", I7," JJ1= ",I7," II2 = ", I7," JJ2= ",I7)') &
!  i,j,coupleIJ(1,1),coupleIJ(1,2),coupleIJ(2,1),coupleIJ(2,2),coupleIJc(1,1),coupleIJc(1,2), coupleIJc(2,1),coupleIJc(2,2)
!***************************************************************************************************
!*********** maintenant qu'on est sur que le binome de jonction est stable energitiquement
!*********** on selectionne un des 2 couples qui minimise les les segments I et J
!***************************************************************************************************
        ii = CoupleIJ(1,1) ;       jj = CoupleIJ(1,2)
        x = icose2v(li,bveclin(:,ii)) + icose2v(lj,bveclin(:,jj))
        ii = CoupleIJ(2,1) ;       jj = CoupleIJ(2,2)
        rv(1) = icose2v(li,bveclin(:,ii)) + icose2v(lj,bveclin(:,jj))
! cas des systemes secants
        if(etatsys(i,j) > IUN ) then
           if(rv(1) > x) then
              ii = coupleIJ(2,1) ; jj = coupleIJ(2,2)
           else
              ii = coupleIJ(1,1) ; jj = coupleIJ(1,2)
           endif

! Derniere condition: il ne faut pas que l'un des binomes de jonction soit antiparallel au segments initial
           if(nbrot(i,ii,1) < itrois .and. nbrot(j,jj,1) < itrois) then
              axejonci(i,j) = ii ; axejoncj(i,j) = jj
           endif
! si coplanaire il y aura 2 autres couple a prendre en compte :coupleIJc
        else
           ii = CoupleIJc(1,1) ;       jj = CoupleIJc(1,2)
           rv(2) = icose2v(li,bveclin(:,ii)) + icose2v(lj,bveclin(:,jj))
           ii = CoupleIJc(2,1) ;       jj = CoupleIJc(2,2)
           rv(3) = icose2v(li,bveclin(:,ii)) + icose2v(lj,bveclin(:,jj))
           if(rv(3) > rv(2) .and. rv(3) > rv(1) .and. rv(3) > x) then
              ii = coupleIJc(2,1) ; jj = coupleIJc(2,2)
           elseif(rv(2) > rv(1) .and. rv(2) > x) then
              ii = coupleIJc(1,1) ; jj = coupleIJc(1,2)
           elseif(rv(1) > x) then
              ii = coupleIJ(2,1) ; jj = coupleIJ(2,2)
           else
              ii = coupleIJ(1,1) ; jj = coupleIJ(1,2)
           endif
           if(nbrot(i,ii,1) < itrois .and. nbrot(j,jj,1) < itrois) then
              axejonci(i,j) = ii ; axejoncj(i,j) = jj
           endif
        endif
    enddo
enddo

! extension du tableau nbrot(...) au systemes devies
do i=1,nbase,1
    do j=1,nbase,1
!        if(SegDevIJ(I,J) /= izero) then
!          print *, "I,J,SegDevIJ",I,J,segdev(I,IUN),segdev(I,Ideux),SegDevIJ(I,J)
!        endif
        if(segdevIJ(i,j) /= 0) then
           nbrot(i,j,1) = nbrot(segdevIJ(i,j),j,1)
        endif
    enddo
enddo


! ICI on dit definit le numero de la loi de vitesse associee a chaque vecteur de la base
ii = izero
do i = 1,Nb_slip_types          ! boucle sur les type de systemes
    do j=1,Slip(i)%Nsystemes ! boucle sur le nombre de systeme dans chaque type
        do k=1,br            ! boucle sur Nbasered
            ii = ii + 1        ! donne l'indice dans nbase

! Nloi(i,j) a ete lu dans le fichier materiau. Il donne pour le type de systeme
! de glissement (i) la loi associee au caractere vis, m1, coin, m2 (j = 1,2,3,4)
            Numero_loi(ii) = Nloi(i,modulo(ii-1,nbasered/IDEUX) + 1 )
!            print *, ii, modulo(ii-1,nbasered/IDEUX)+1,Numero_loi(ii)
        enddo
    enddo
enddo


! =======================================================================
! =========  Pour les problemes de sous-reseau  ========================
! =======================================================================
!  Definition des vecteur de translation unitaire pour paver l'espace
!  Ceci est utilise seulement pour les HCP
! =======================================================================
!  elles servent a detecter les probleme des sous reseau. le reper est constitue de :
! ici on defini les trois translations elementaires pour la construction de
! l'espace, c-a-d le reseau principal.
! ces dierctions sont :
!         a1 : la direction a1 : 1/2 [ 0  1 -1 ] = [0 66 -66]
!         a2 : la direction a2 : 1/2 [-1  1  0 ] = [-66 66 0]
!         c : la direction normal : 2/3 [1 1 1]  = [88 88 88]
! les points de l'espace doivent etre accessibles apartir du point [000] et de ces
! translations uniquement.

! pour les HCP, les trois direction entieres doivent etre :
! a1 = [0 66 -66]
! a2 = [-66 0 66]
! c = [-88 -88 -88]

if (HCP) then
   x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
   y_reseau(:) = bveclin (:,nbasered+1)   !========  vecteur vis 2eme systeme
   z_reseau(:) = bveclin (:,3)            !========  vecteur coin 1er systeme
elseif (ORT) then
   x_reseau(:) = bveclin (:,nbasered+1)   !========  vecteur vis 2eme systeme
   y_reseau(:) = bveclin (:,3)            !========  vecteur coin 1er systeme
   z_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
elseif (BCC) then
!  Attention: si les vecteurs x,y,z ne sont pas egaux
! aux bveclin 1,4,25 la procedure de decomposition en vecteurs elementaires
! est fausse
   x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
   y_reseau(:) = bveclin (:,4)            !========  vecteur mixte 1er systeme
   z_reseau(:) = bveclin (:,25)           !========  vecteur coin 1er systeme
elseif (CFC) then
! Pour le CFC on se sert des x_reseau, y_reseau et z_reseau uniquement pour le calcul du facteur de boite
   x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
   y_reseau(:) = bveclin (:,3)            !========  vecteur coin 1er systeme
   z_reseau(:) = bveclin (:,2*nbasered+1) !========  vecteur vis 3eme systeme
elseif (CS) then
! Pour le CS on se sert des x_reseau, y_reseau et z_reseau uniquement pour le calcul du facteur de boite
   x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
   y_reseau(:) = bveclin (:,3)            !========  vecteur coin 1er systeme
   z_reseau(:) = bveclin (:,nbasered+7) !========  vecteur vis 3eme systeme
elseif (MGO) then
   x_reseau(:) = bveclin (:,2*nbasered+3)   !6 0 0
   y_reseau(:) = bveclin (:,4*nbasered+7)   !0 6 0
   z_reseau(:) = bveclin (:,3)    !0 0 6
elseif (DC) then
! Pour le CFC on se sert des x_reseau, y_reseau et z_reseau uniquement pour le calcul du facteur de boite
   x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
   y_reseau(:) = bveclin (:,3)            !========  vecteur coin 1er systeme
   z_reseau(:) = bveclin (:,2*nbasered+1) !========  vecteur vis 3eme systeme

else
   stop  " vecteurs de translation elementaires non definis pour cette cristallo, by"
endif

! =======================================================================
!  rappel des unites de translation pour la definition du reseau spatial
! ====    elles sont definies dans le modules varglob    ==============
! =======================================================================
write(*,*)
write(*,'("The elementary translation vectors ")')
write(*, '("  X  = [ ",3(I6,X)," ]")') x_reseau(1:3)
write(*, '("  Y  = [ ",3(I6,X)," ]")') y_reseau(1:3)
write(*, '("  Z  = [ ",3(I6,X)," ]")') z_reseau(1:3)

write(*,*)
write(*,'("The elementary segments length")')
write(*,'("  Screw     = ",F8.4," nm")')  normlin(1)*avalue*1D9
write(*,'("  Mixted 1  = ",F8.4," nm")')  normlin(2)*avalue*1D9
write(*,'("  Edge      = ",F8.4," nm")')  normlin(3)*avalue*1D9
write(*,'("  Mixted 2  = ",F8.4," nm")')  normlin(4)*avalue*1D9

!#######################################################################################
!Some time it is necessary to test calculate the "facteur_boite" directly in mm
calculboite="N"   ! The default mode is No

if (calculboite == "Y") then
   print*,"Calcul du facteur_boite !"
!   read(*,*) calculboite

! Trois longueurs de reference vis coin et mixte pour la discretisation
!** ici on va calculer les dimensions du plus petit volume simule et le facteur_boite
! Pour trouver les dimensions du plus petit volume simule on resout les 3 systemes d equations:
! alpha1*(xreseau)+beta1*(yreseau)+gamma1*(zreseau)= (L1,0,0)
! alpha2*(xreseau)+beta2*(yreseau)+gamma2*(zreseau)= (0,L2,0)
! alpha3*(xreseau)+beta3*(yreseau)+gamma3*(zreseau)= (0,0,L3)
! L1,L2,L3 dimensions du plus petit volume simule
! ce calcul iteratif n est pas optimise et est assez long

      systeme(1:3,1)=x_reseau
      systeme(1:3,2)=y_reseau
      systeme(1:3,3)=z_reseau

      if (HCP) then
         systeme(1:3,1)=2*x_reseau
         systeme(1:3,2)=2*y_reseau
         systeme(1:3,3)=2*z_reseau
      endif

      Ldim(:)=0
      Boucle:do i=1,3
          Ldimtemp(:)=0
          dimensions(1:3)=0
          do alpha=-200,200
              do beta=-200,200
                  do gamma=-200,200
                      inconnues(1)=alpha
                      inconnues(2)=beta
                      inconnues(3)=gamma
                      dimensions=matmul(systeme,inconnues)
                      dimtemp(:)=0 ; dimtemp(i)=dimensions(i)
                      if (dimensions(1)==dimtemp(1) .and. dimensions(2)==dimtemp(2) .and. dimensions(3)==dimtemp(3) &
                      .and. entier(dimensions(i))) then
                         if (Ldimtemp(1)==0 .and. Ldimtemp(2)==0 .and. Ldimtemp(3)==0) Ldimtemp(:)=int(dimensions(:),DPI)
                         if (norvect(dimensions(:)) < norvect(real(Ldimtemp(:),DP))) Ldimtemp(:)=int(dimensions(:),DPI)
                      endif
                  enddo
              enddo
          enddo
          Ldim(:)=abs(Ldimtemp(:))+Ldim(:)
      enddo Boucle
      print*,Ldim(:)

! calcul du facteur_boite qui est le ppcm des dimensions du plus petit volume
! methode differente de la fonction ppcm3: on cherche les entiers multiplicateurs des dimensions
! qui permettent d obtenir le ppcm (evite de faire une boucle sur un nombre gigantesque d iterations
! lorsque les dimensions sont importantes).
      bppcm:do i=1,200
          do j=1,200
              do k=1,200
                  if (Ldim(1)*i==Ldim(2)*j .and. Ldim(1)*i==Ldim(3)*k) then
                     facboite=Ldim(1)*i
                     exit bppcm
                  endif
              enddo
          enddo
      enddo bppcm
      print*,"le facteur_boite vaut: ",facboite
      read(*,*)
endif

if (decalageCPL) then
! CLP shiftees: tableau de vecteurs elementaire de shifts
! on fait en sorte que les decalages soient des vecteurs du reseau
! i=1 ---> sortie face a
! i=2 ---> sortie face b
! i=3 ---> sortie face c
   do i=1,3,1
       SHIFT(1:3,i)=nint(real(DECALCLP(1:3,i),DP)/Facteur_boite,DPI)*Facteur_boite
   enddo
endif
end subroutine Tabulation

!************************************************************************
! Ce sous-programme ecrit dans le fichier (BVD.cfc, cc,cs, etc.)
! les vecteurs de la base de simulation
!  Ghiath Monnet -  07/02/01
!************************************************************************
subroutine write_basedevecteurs

implicit none

integer (kind=DPI) :: i,j,itemp
character(len=40)  :: file_BVD

61 format(I3,': (',3(I5,1X),')[',3(I6,1X),']. Dep: [',3(I6,1X),']')
63 format(I5,"  ",3(1X,I6),"  ",3(I6,1X),"      ",3(I3,1X))

file_BVD = '../out/BVD'

#ifdef PA
! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_BVD = file_BVD(1:len_trim(file_BVD))//"2"
endif
#endif

file_BVD = file_BVD(1:len_trim(file_BVD))//"."//crystal_structure(1:3)

open (14,FILE=file_BVD,STATUS='UNKNOWN')  ! fichier de sortie

write (14,*) nbase , avalue         ! unique de programme : lu par camera ....
write (14,*) "   "
do i = 1,nbase
    write (14,63) i, bveclin(1:3,i),bvecdep(1:3,i),bvecnor(1:3,i)
    if(modulo(i,ihuit) == izero) write (14,*) "   "
enddo
write(14,*) "  "
write(14,*) "  "

write(14,*)  NTSG, NbSysDev,"       !   NTSG ;    nb of combinations for cros slip "
write(14,*) " "
do i = 1,NTSG
    write(14,*)  i,sysdev(i,1:NbSysDev)
enddo



write(14,*) "  "
write(14,*) "  "
if (HCP) then
   write(14,*) "Indice        Plan            ligne (/22)                   deplacement "
elseif (BCC) then
   write(14,*) "Indice        Plan            ligne (/6)                   deplacement "
else
   write(14,*) "Indice        plan       ligne                 deplacement "
endif
! ensuite il y a l'ecriture formatee de la BVD pour une lecture agreable
do itemp = 1,NTSG
    write (14,*) " Systeme numero : ",itemp
    do j = 1,8
        i = (itemp - 1) * nbasered + j
        if (HCP) then
           write (14,61) i, bvecnor(1:3,i), bveclin(1:3,i)/22,bvecdep(1:3,i)
        elseif (BCC) then
           write (14,61) i, bvecnor(1:3,i), bveclin(1:3,i)/6,bvecdep(1:3,i)
        else
           write (14,61) i, bvecnor(1:3,i), bveclin(1:3,i),bvecdep(1:3,i)
        endif
    enddo
enddo

close(14)

end subroutine write_basedevecteurs

!**************************************************************************************
!> In this subroutine we load the initial segment configuration and we setup the
!! simulation loading parameters in agreement with the key (mode_deformation)
!**************************************************************************************
subroutine Configuration

implicit none

integer(kind=DPI) :: I,compt
real(kind=DP)     :: temp

! =======================================================================
#ifndef MDC
if(sideja /= ideux) then
  write(*,*)
  write(*,*) "Elementary displacement velocity per time step"
   temp = normdep(1)*avalue/deltat
   write (*,'("  Screw direction  = ",F15.6," m/s")') temp
   temp = normdep(2)*avalue/deltat
   write (*,'("  Mixte direction  = ",F15.6," m/s")') temp
   temp = normdep(3)*avalue/deltat
   write (*,'("  Edge direction   = ",F15.6," m/s")') temp
endif
#endif


DCFJ = DCFJ * BdivA  ! DCFJ = Distance for Calculation of Force on Junction arms (in avalue unit)

! The compatibility between DCFJ, xlomax and echelle is tested
write(*,*)
write (*,'(" The junction force calculation DCFJ length  = ",I5," BVD elementary vector")') INT(DCFJ / maxval(Normlin(:)))
if ((INT(DCFJ / maxval(Normlin(:))) < 1) .or. ((XLOMAX*half/DCFJ) < 1)) then
  print*,INT(DCFJ / maxval(Normlin(:)))
  print*,XLOMAX*half/DCFJ
  write (*,*) "Stop: the DCFJ value does not seems to be consistant with the ECHELLE or LDIS_ACT values"
  stop
endif

Tauint_limite = Tauint_limite *1.0D6 / Xmu
write (*,'("Stress amplitude fixing the velocity screening in stress shaft = ",F15.2," MPa")') &
TauINT_LIMITE * 1.0D-6 * Xmu
print*, " "

#ifndef MDC
if(krc == 0)stop 'KRC == 0. Interdit , Bye'
#endif

! Some variables must be initialized as a function of the loading conditions
STRATE       = .false.       ! Clef de control en vitesse de def = cste
METFUT       = .false.       ! Clef de control en metallofute
METFUT2      = .false.       ! Clef de control en metallofute2
FATIG        = .false.       ! Clef de control en fatigue
CREEP        = .false.       ! Clef de control en fluage \ingroup creep

#ifdef MDC
if (mode_deformation_key) then
#endif
!> if select case \a mode_deformation equals
select case(mode_deformation)

!> - case (0) the strain rate deformation mode (uniaxial loading)
case (0)
   STRATE       = .True.
   FSIG         = 0
   EPSMAX       = 0
   write(*, '("Loading type: Uniaxial constant strain rate",/)')
   relative_rate = .false.

!> - case (1:2) the carto mM loading mode
case (1:2)

!> - case (3) the stress rate deformation mode (uniaxial loading)
case (3)
   FSIG         = 0
   EPSMAX       = 0
   write(*,'("Loading type: Uniaxial constant stress rate",/)')

!> - case (4) the constant stress deformation mode (uniaxial loading)
case (4)
   SigmaPoint   = 0.0
   DeltaEpsilon = 0.0
   FSIG         = 0
   EPSMAX       = 0
   write(*, '("Loading type: Uniaxial constant stress (sigma0)",/)')

!> - case (5) the fatigue deformation mode (uniaxial loading)
case (5)
   FATIG        = .true.
   DeltaEpsilon = 0
   write(*, '("Loading type: Fatigue ",/)')

!> - case (6) the MetalloFute deformation mode (uniaxial loading)
case (6)
   METFUT       = .true.
   EPSMAX       = 0
   write(*, '("Loading type: MetalloFute",/)')

!> - case (7) the MetalloFute2 deformation mode (uniaxial loading)
case (7)
   METFUT2      = .true.        ! Clef de control en metallofute
   EPSMAX       = 0
   write(*, '("Loading type: MetalloFute2",/)')

!> - case (8) the Finite Element tabulated stress loading mode (heterogeneous loading)
case (8)
   FSIG         = 0
   EPSMAX       = 0
   shear        = .false.
   write(*, '(" Loading type: imposed constant tabulated stress field",/)')

!> - case (9) the creep mode loading defintion (uniaxial loading) \ingroup creep
case (9)
   CREEP        = .true.
   SigmaPoint   = 0.0
   DeltaEpsilon = 0.0
   FSIG         = 0
   EPSMAX       = 0
   write(*, '("Loading type: Uniaxial creep mode at constant stress (sigma0)",/)')

   if (SiDeja /= 2) then
    write(*, '(" ",/)')
    write(*, '("The mode_deformation = 9, i.e. Creep mode can be applied only with SiDeja = 2 !!!",/)')
    stop
   endif

end select
#ifdef MDC
endif
#endif

!=============================================
!> Loading of the initial segments configuration
!=============================================

file_fichseg = fichseg
#ifdef PA
! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_fichseg = file_fichseg(1:len_trim(file_fichseg))//"2"
endif
#endif

! The segment initial configuration is loaded from file_fichseg
open(3,file=file_fichseg,STATUS='OLD')

! solli_sys(i) = 0 means a forest segment
! solli_sys(i) = 1 means a standard segment
read(3,*) Solli_Sys(1:NTSG)

! Nombre de segments dans la configuration initiale
read(3,*) Nsegm


read(3,*)  modur(1),modur(2),modur(3)

close(3)

if(Nsegm == izero) then

#ifdef PA
   if (Nb_Phase == 2) stop "The junction input mode cannot be used with Nb_phase = 2"
#endif
    !> When (Nsegm == izero), an alternative segment file is loaded to test junction reactions
   Solli_Sys(1:NTSG)= zero
   open (3,FILE='../in/junction_def',STATUS='OLD')
   read(3,*) Nsegm  ! nombre de segemnts initiaux
   read(3,*) temp
   read(3,*) temp
   read(3,*) temp
   read(3,*) I
   read(3,*) temp
   read(3,*) temp
   read(3,*) temp
   read(3,*) temp
   solli_sys(syseg(I)) = temp
   print *,"1", I, syseg(I),solli_sys(syseg(I))

   if(nsegm == 2) then
      read(3,*) I
      read(3,*) temp
      read(3,*) temp
      read(3,*) temp
      read(3,*) temp
      solli_sys(syseg(I)) = temp
      print *,"2", I, syseg(I),solli_sys(syseg(I))
   endif
endif
close(3)

effet_longueur = .false.

Do I = 1, NLV
    if(Loi(i)%arrhenius == 2) effet_longueur = .true.
enddo

write(*,*) "Length effect in mobility laws = ", effet_longueur

! The maximum displacement is first defined with respect to the discretization
! length with respect to all segment directions
compt = Izero
Dlimite_max = izero

do I=1,NBASE
    LOMA(I)     = INT(XLOMAX/NormLin(I),DPI)       !*** PLUS RAPIDE
    if (effet_longueur .and. tyseg(I) == 1 ) LOMA(I)  = LOMA(I) * Icinq
    if(Dlimite_max < loma(i)) Dlimite_max = Loma(i)
    JSLOMA(I)   = INT(XLOMAX/(NormLin(I)*4.0),DPI) !*** PLUS RAPIDE
    if (solli_sys(syseg(i)) < numtols) then
       LOMA(I) = INT(XLOMAX_nact/NormLin(I),DPI)     !*** PLUS RAPIDE
       if (effet_longueur .and. tyseg(I) == 1 ) LOMA(I)  = LOMA(I) * Icinq
       compt = compt + IUN
    endif

  if(loma(i) < 10) then
    print *, " Stop : The discretization length of the VecLin :",i
    print *, " is too small! The simulation ECHELLE shoulb be decreased  or "
    print *, " the discretization length XLOMAX should be increased."
    print *, " LOMAI = ", loma(i)
    print *, " XLOMAX = ", XLOMAX
    print *, " NormLin = ", NormLin(I)
    stop
  endif
enddo

! The maximum displacement length is in a second time compared to the input variable Facteur_depmax
! Facteur_depmax defined a second upper limit if we consider a mean displacement of one at each steps
if(Dlimite_max > facteur_depmax) Dlimite_max = INT(Facteur_depmax * un,DPI)

! The Dlimite value is defined for all the mobility laws
do i = 1, NLV
    Dlimite(I) = Dlimite_max
enddo

if (compt /= Izero) then
   write(*,'(/," Discretisation de systemes non actifs pour :", I8, " systemes.")') compt/nbasered
endif

write(*,*)
write(*, '(" Important simulation parameters (in AVALUE)")')
write(*, '("  Maximum displacement length per step = ",I8)') Dlimite_max
write(*, '("  Screw discretization length = ",I8)') loma(1)
write(*, '("  Mixte discretization length = ",I8)') loma(2)
write(*, '("  Edge discretization length  = ",I8)') loma(3)

!=====================
! Some important variables initialization
!=====================
! Longueur de discretisation de la courbure en unite a ???
ICROSMAX = INT(XLOMAX/IQUATRE,DPI)

! Size of the tiny loops that must be eliminated during the dynamic (hence modeling diffusion effect)
! The size of such loops is a priori defined with the parameter Microboucle_ratio (in 01Constant),
! but this length is modified in the subroutine generer_particule if particules are considered.
Microboucle = int((float(XLOMAX) * Microboucle_ratio),DPI)
write(*,*)
write (*,'("  Dimension of the collapsing loops (in a) = ", I6," BVD  = ",F5.2," nm")') &
Microboucle, Microboucle * avalue * 1D9

! Defintion of the regularization length used in the curvature definition for
! the local line tension. This length is scaled by the discretization length
! All the segments participating to the local line tension will be excludded
! from the interactions (self-interaction)
rlocal = INT(float(xlomax) * fac_rlocal,DPI)
write (*,'("  Dimension of the line tension regularisation length (in microns) = ",F5.2)') &
      &  real((deux * rlocal * avalue * 1.D6),DP)
rlocal2 = rlocal * rlocal

! =============================================
!>  Definition of the slip systems mobility law
! =============================================

write(*,*)
write(*,*) "Number of velocity law taken into account = ", NLV
If(NLV < ideux) stop "NLV must be > 1, The velocity law /1/ is dedicated to the relaxation steps!!!"

! Using IS unites for the mobility laws
Do I = 1, NLV
    if(i== 1 .and. Loi(i)%arrhenius /= 0) then
       print *, " For the sake of simplicity, the first velocity law Schid be Viscous"
       print *, " Please use a law with arrhenius type = 0                           "
       stop
    endif
    if(Loi(i)%arrhenius == 2) then
       Loi(i)%friction = Loi(i)%friction  * 1E+6        ! contyrainte au plateau athermique
       Loi(i)%tau0 = Loi(i)%tau0 * 1E+6
    elseif(Loi(i)%arrhenius == 1) then
       if(Loi(i-1)%arrhenius /= 2) then
          print *, " Please insert the mobility law of type 2 (length effect) before law number: ",I
          stop
       endif
! the following procedure is define mobility law for non screw dislocations:
! the idea is to use the same Boltzmann factor but with a prefactor length-independant
! Therefore : the prefactor v0 = ratio (introduced in Materiau) * h(screw) * 1m
! the division by avalvue is to put put the velocity in simulation unite, ie. a/s
       Loi(i)%V0 = rapport * Loi(i-1)%h * 1E-6 / avalue    ! attac frenquency for thermal activated mobility (s-1)
       Loi(i)% deltag0 = Loi(i-1)% deltag0         ! energie d'activation totale ( pour T= 0 K)
       Loi(i)%friction =  Loi(i-1)%friction         ! contyrainte au plateau athermique
       Loi(i)%tau0 = Loi(i-1)%tau0
       Loi(i)% coef_p = Loi(i-1)% coef_p  ! exposant p de la loie de vitesse (pour les vis)
       Loi(i)% coef_q = Loi(i-1)% coef_q ! exposant q de la loie de vitesse (pour les vis)

       print *, " Velocities ratio Edge/Vis(1m) :", Loi(i)%V0*avalue / Loi(i-1)%h/1E-6

       elseif(Loi(i)%arrhenius == izero) then
          Loi(i)%BF = Loi(i)%BF * avalue
          Loi(I)%friction  = Loi(i)%friction * 1E+6
    endif
enddo

If ( NLV > NLV_max) stop " nombre de lois de vitesse trop grand"

end subroutine configuration


!#############################################################################
!# Desallocation des tableaux de la simulation (allacated by dimentionner_base)
!# SOUS-PROGRAMME ANNULER_DIMENSIONNER_BASERED
!################################################################# 05/11/98 ##
subroutine Desallocation

implicit none

deallocate (BVECLIN)
deallocate (BVECDEP)
deallocate (BVECNOR)
deallocate (courbuO)
deallocate (courbuE)
deallocate (ASSOC)
deallocate (INVSI)
deallocate (DEASSOC)
deallocate (TYSEG)
deallocate (SISEG)
deallocate (SEGDev)
deallocate (CONEC)
deallocate (DECONEC)
deallocate (SC)
deallocate (CoefSC)
deallocate (AXEJONCi)
deallocate (AXEJONCj)
deallocate (SYSOBST)
deallocate (SYSCONNEC)
!deallocate (LIEN)
!deallocate (LIENS)
deallocate (NBROT)
deallocate (GDROT)
deallocate (INVCONEC)
deallocate (PartVis)
deallocate (ProjVis)
deallocate (Projcoin)
deallocate (SYSEG)
deallocate (LOMA)
deallocate (JSLOMA)
deallocate (FAC1DEP)
deallocate (FAC2DEP)
deallocate (MODDEPLA)
deallocate (VECNORLIN)
deallocate (VECNORDEP)
deallocate (NORMLIN)
deallocate (NORMDEP)
deallocate (RAUSYS)
deallocate (RAUSYS_precipite)
deallocate (NJONCSYS)
deallocate (CRANSYS)
deallocate (GAMSYS)
deallocate (SchmidSys)
deallocate (Signe_TAUAPP)
deallocate (TensRotSys)
deallocate (TensEpsSys)
deallocate (TrAppSys)
deallocate (TrIntSys)
deallocate (GAMMADOTSYS)

if ((.not. allocation_dynamique_boites) .and. calculate_gammabox ) deallocate (Gammabox)

if (desorientgrain) then
  if (GB == 3) then
    deallocate (PlanMiller)
    deallocate (PlanPos)
  endif
  deallocate (MatRotGrain)
  deallocate (rauGRAINsys)
  deallocate (rausysjoncGRAIN)
  deallocate (RAUGRAINjonc)
  deallocate (RAUGRAIN)
  deallocate (AIRESYSGRAIN)
  deallocate (TensRotGrain)
  deallocate (TensRotSysGrain)
  deallocate (Gamsysgrain)
  deallocate (SchmidSysGrain)
  deallocate (SigAppOriente)
endif

deallocate (Bf)
deallocate (Tau0)

#ifdef MDC
#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
if (Mon_Rang ==  IZERO) then
#endif

 deallocate(sweptsurfdata)

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
endif
#endif

#endif
end subroutine Desallocation

!#############################################################################
!> This subroutine is used to load and define all the simulation parameters
!#############################################################################
subroutine LOAD

implicit none

!#########################################################################
!# Parametres de controle : chargement et reconstitution des donnees     #
!# physiques ou topologiques utile a la simulation                       #
!#########################################################################
call lire_donnees

#ifndef MDC
! if the loading axis is set to (0,0,0), a reference stress tensor is used
! in the loading procedure. We are not in uniaxiale condition !
if (norvect(Z) < 0.01) then

  UNIAXIALE = .FALSE.

  if (shear) then
    write(*,*) 'The shear mode control can be used only with uniaxiale test'
    write(*,*) 'The shear control parameter must be .false.'
    stop
  endif

  call Lire_tenseur

else

  UNIAXIALE = .TRUE.

endif
#else
 if (mode_deformation_key) UNIAXIALE = .TRUE.
#endif

call parametrage ! Contains (indirectly) DD<- Z No. 2

! With this mode we need to initialize (load) the tab FE_SigappS
if (mode_deformation == 8) then
  call Lire_FE_SigappS
  if (shear) stop 'Simulation with mode_deformation = 8 cannot be used with the (shear) key option'
endif

if (mode_deformation == 9) then
  call Lire_climb
endif

if (key_crack) then
  allocate (Sigcrack(6))
endif

!#############################################################################
! Initialisation des fichiers en fonction du mode de demarage de la simulation
!#############################################################################
!call initialisation_fichiers

!==========================================================================================
!======== Connaissant desormais le nombre total des systemes de glissement, on allouer=====
!======== les tableaux relatives a la base de discretisation ============================
! =============================================================
! Allocation de la dimention des tableaux utils a la simulation
call allocation

!==========================================================================================
!========= ensuite on genere la base de vecterurs de discretisation ======================
!==========================================================================================
! Generation et definition des caracteristiques des segments de dislocations
!call generer_basedevecteurs

! Ecriture de la base de vecteurs de discretisation format brut et clean

call read_basedevecteurs

if(rotationBVD) call rotation_BVD

call initialiser

!#######################################################################
! tabulations annexes decrivant les correlations de la base des vecteurs
!#######################################################################
call tabulation

call write_basedevecteurs

if(cartograph == 1) THEN
   call carto
   write(*,'(///)')
   Print *, " GENERATION DES FICHIERS DE CARTO DANS ../IN/CARTO/ TERMINEE"
   Print *, " SORTIE NORMALE DU PROGRAMME"
   PRINT*,  " POUR EXECUTER DD, VEUILLEZ METTRE /cartograph/ a UN ou DEUX"
   write(*,'(///)')
   STOP
ENDIF


end subroutine LOAD

!##########################################################################
!# It is sometime usefull to define the material isotropic stiffness tensor
!# In this subroutine the Lame constante definition is used.
!##########################################################################
Subroutine MakeStiffTens

implicit none

! Initialization
Lambda = (2.* DPoiss * XMu) / (1 - (2 * DPoiss))
InvXMu = 1. / XMu

! A siffness matrix is defined with the Lame constantes
Stiffness(:,:) = zero
Stiffness(1,1) = (Lambda + 2. * XMu)
Stiffness(2,2) = (Lambda + 2. * XMu)
Stiffness(3,3) = (Lambda + 2. * XMu)
Stiffness(1,2) = Lambda
Stiffness(1,3) = Lambda
Stiffness(2,3) = Lambda
Stiffness(2,1) = Lambda
Stiffness(3,1) = Lambda
Stiffness(3,2) = Lambda
Stiffness(4,4) = XMu
Stiffness(5,5) = XMu
Stiffness(6,6) = XMu

end Subroutine MakeStiffTens

end module INIT

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
