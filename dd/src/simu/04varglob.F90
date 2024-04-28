
!===================================================================================================
!========================    DEBUT    MODULE  "VARGLOB"   ==========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains the definition of the main constants used in the microMegas modules.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module VARGLOB

use constantes  ! les constantes de la simulation
use VARBASE     ! tableau de la base

implicit none

integer(kind=DPI) :: domboundsX(2),domboundsY(2),domboundsZ(2)

logical       :: key_infline      = .false.   !< PBC induced infinite lines can be artificialy pinned at one point
logical       :: key_nucleation   = .false.   !< Nucleation activation key
logical       :: key_Rot_Dom      = .false.   !< key to impose crystal rotation in domains

real(kind=DP),allocatable,save :: TAU0(:),  & !< Contrainte resolue de Peierls
                                  BF(:)       !< coefficient de frottement par VBD

!**** POUR LES STAT...
Integer(kind=DPI),allocatable,save ::       &
njoncsys(:)          !< Number of junction for each system

real(kind=DP),allocatable,save ::           &
     RAUSYS(:),          &!< Densite par systeme de glissement
     RAUSYS_precipite(:),&!< Densite par systeme de glissement
     RAUSYSjonc(:),      &!< Densite par systeme de glissement des jonctions
     CRANSYS(:),         &!< crans formes  par unite de longueur de dislo pour chaque systeme
     GamSYS(:),          &!< Gamma par systeme de glissement
     Gammabox(:,:),      &!< Gamma for each slip system on each multipole box
     TensRotSys(:,:,:),  &!< Elementary rotation tensor for each slip system (Wij)
     TensEpsSys(:,:,:),  &!< Elementary deformation tensor for each slip system (Eij)
     TensDistSys(:,:,:), &!< Elementary distortion tensor for each slip system (Bij)
     AireSYS(:),         &!< aire ballaye en a2 par systeme de glissement
     AireSYSInst(:),     &!< aire ballaye en a2 par systeme de glissement par iteration
     AireVisSYS(:),      &!< aire ballaye en a2 par systeme de glissement
     AireCoinSYS(:),     &!< aire ballaye en a2 par systeme de glissement
     AireMixtSYS(:),     &!< aire ballaye en a2 par systeme de glissement
     SchmidSys(:),       &!< facteur de schmid signe par systeme
     TrAppSys(:),        &!< Travail applique cumule par systeme
     TrAppSysInst(:),    &!< Travail applique instantane par systeme
     TrIntSys(:),        &!< Travail interne cumule par systeme
     TrIntSysInst(:),    &!< Travail interne instantanne par systeme
     GAMMADOTSYS(:),     &!< Instantaneous strain rate on each slip systems
     Signe_TAUAPP(:)      !< Sign of the applied stress on each slip systems (use in APB calculations)

real(kind=DP),save ::   L_boite    ! Length of the Greengard Domains

Logical :: Normalise_SR ! clef de normalisation de strain rate par rau pour 1012 m-2
Logical :: RCMAKE  ! clef de reallocation dynamic
Logical :: Methode_boites  ! clef d'utilisation de la methode des boites (NBOITES = 1 ou non)
Logical :: Particules  ! clef de presence de particules dans la siomulation
Logical :: Loops          ! logical variable indicating the presence of small loops in the simulation box
Logical :: relative_rate  ! if T, the "fixed" strain rate is reevaluated to keep the mean dislocation
                          ! velocity constant
real(kind=DP) :: strainrate_ini ! vitesse de deformation initiale
real(kind=DP) :: rau_ini ! initial dislocation density
real(kind=DP) :: L_average_screw ! average length of screw segments
real(kind=DP) :: tau_average_screw ! the average effective stress on screw segments
real(kind=DP) :: L_average_edge ! average length of edge segments
real(kind=DP) :: tau_average_edge ! the average effective stress on edge segments
real(kind=DP) :: rau_screw ! total density of segments considered as screw
real(kind=DP) :: rau_edge !total density of segments considered as non-screw
real(kind=DP) :: volboite, BdivAboite ! For Gammabox calculation
real(kind=DP) :: Critic_angle   !< Critical angle to define screw, edge and mixte line curvature

!ATTENTION pour debuggage
integer (kind=DPI) :: iterinfo  ! : l'iteraction de debuggage du code
integer (kind=DPI) :: sysinfo  ! : si > 0, les info de debugage seront donnees seuelement sur ce systeme
integer (kind=DPI) :: MicroBoucle  ! taille limite (en a) pour la suppression des boucles
integer (kind=DPI) :: NBASE     ! dimension de la base des vecteur utilises pour
integer (kind=DPI) :: Modep_max ! le modep max relatif a la morphologie et la crystallographie
integer (kind=DPI) :: VL_CS(Nsegmax) ! table giving the new line vector for a cros slip segment
integer*4          :: koderr    ! pour eviter un arret du calcul quand la lecture d un fichier depasse sa taille


!==================================================================
!===== Variables used in the different control modes ==============
!==================================================================
real(kind=DP)::         &
       TrApp           ,&! Travail applique cumule
       TrAppInst       ,&! Travail applique instantane
       TrInt           ,&! Travail interne cumule
       TrIntInst       ,&! Travail interne instantanne
       Z(3)            ,&! Axe de sollicitation
       ZPrime(3)       ,&! Axe de sollicitation modified by the rigid body rotation
       TensRot(3,3)    ,&! Plastic rotation tensor
       TensEps(3,3)    ,&! Plastic deformation tensor
       TensDist(3,3)   ,&! Plastic distortion tensor
       EPSMAX          ,&! Amplitude de la deformation par cycle de fat.
       FSIG            ,&! Signe du premier cycle de fatigue
       RAID            ,&!< The effective stiffness used in constant strain rate control
       DeltaEpsilon    ,&! Amplitude de Deformation pour le calcul de la vitesse de deformation
       SigmaPoint      ,&! Vitesse d increment de contrainte = cste
       EpsilonPoint    ,&! Vitesse d increment de deformation = cste
       Epspointsign    ,&! Sign of the strain rate
       Angle_vis       ,&! angle entre la ligne moyenne et la direction de b
       stress_fact     ,&!< The factor applied to the stress value in control mode 8 (FE_sig)
       stress_fmin     ,&!< The min factor applied to the stress value in control mode 8 (FE_sig)
       stress_fmax     ,&!< The max factor applied to the stress value in control mode 8 (FE_sig)
       ref_deform        !< The reference deformation amplitude need for normalization in control mode 8 (FE_sig)

logical::                    &
      STRATE                ,&  !< Constant strain rate control key
      METFUT                ,&  !< Clef de control en metallofute
      METFUT2               ,&  !< Clef de control en metallofute2 (materiaux a forte friction de reseau)
      FATIG                 ,&  !< Fatigue control key
      CRACK                 ,&  !< Crack control key
      UNIAXIALE             ,&  !< false quand il faut lire le tenseur dans un fichier separe 9par defaut true)
      decalageCPL           ,&  !< CPL shift control key
      RotationBVD           ,&  !< clef de rotation de CPL
      StrainRateMode8           !< Constant strain rate control key in mode 8

real(kind=DP),dimension(6,6) :: Stiffness       !< The material elastic stiffness tensor
real(kind=DP),dimension(3,3) :: TestSig         !< a stress tensor to do test
real(kind=DP),dimension(3,3) :: SigAppS         !< The applied stress on a segment
real(kind=DP)                :: Lambda,InvXMu   !< The Lame constants

!*** Variables pour le controle : MISE EN OEUVRE

integer (kind=DPI)::         &
        Period              ,& ! Nombre de pas immobiles
        NstatControl        ,& ! Number of steps considered for moving averages in simulation monitoring
        StepControl         ,& ! Time step when loading control start to have some meaning
        Matrice_ROTBVD(3,3) ,& ! Matrice de rotation pour CLP
        SHIFT(3,3)          ,& ! vecteur decalage du reseau
        SHIFTBOITE(3,3)     ,& ! vecteur decalage pour methode des boites
        DECALCLP(3,3)          ! vecteur decalage lu dans cont

real(kind=DP),allocatable::  &
    EPSOLD(:)               ,&! NstatControl precedantes valeures de epsilon
    GAMMAOLDSYS(:,:)          ! par systeme cette fois ci

real(kind=DP)::      &
    Sigma           ,&!< Contrainte de tracation dans la direction de l'axe de traction
    Sigma0          ,&!< Contrainte de raction INITIATE dans la direction de l'axe de traction
    EPSO            ,&!< Deformation totale dans la dir. de l axe de traction
    EPS_pl          ,&!< Deformation plastique dans la dir. de l axe de traction
    EPS_pl_inst     ,&!< Deformation plastique dans la dir. de l axe de traction par iteration
    EPS_el          ,&!< Deformation elastic dans la dir. de l axe de traction
    EPSDOT          ,&!< instantaneous strain rate calculated in the tensile axis direction
    DELSIG            !< Increment de contrainte par pas

!*** Crack variables ***************
logical           ::  key_crack            !< Crack field on/off
integer(kind=DPI) ::  CrackTip(3)         !< Crack tip coordinates (burger vector)
real(kind=DP)     ::  SIF_fact            !< Plate geometric factor
integer(kind=DPI) ::  acrack              !< Crack length

real(kind=DP)     ::  slope               !< Percentage of Applied Stress
real(kind=DP), allocatable, save ::     Sigcrack(:)     !< The stress tensor for 2D crack stress field
!***********************************

logical::  REVERS ! Clef de changement de signe en fatigue
logical::  kkdebug ! Clef de debugage de la simation

!*** Variables Conditions Periodique Aux Limites ***************
integer(kind=DPI) :: ModuR(3)           !< The integer dimension of the simulated periodic volume
integer(kind=DPI) :: ModurMax           !< The largest dimension of the simulated periodic volume
real(kind=DP)     :: InvmoduR(3)        !< The inverse of the dimension of the simulated periodic volume
real(kind=DP)     :: HalfModuR(3)       !< Half dimension of the simulated periodic volume
!**************************************************************

logical            :: OUT(NSEGMAX),OUTO(NSEGMAX),OUTE(NSEGMAX) !*** NOS INDEX DES SEGMENT ALONGES
integer(kind=DPI)  :: TABLONG(NSEGMAX) !*** CORRECTION SUR LES LONGUEURS DES SEGMENTS EN SURFACE
integer(kind=DPI)  :: buftablong,bufptablong(3),bufvtablong(3),buftlad,buftlsi !*** BUFFERS

integer(kind=DPI)  :: CSOED(NSEGMAX,2) ! On coupe le segment en O ou E de D

! The variables connected to the Greengard algo
integer(kind=DPI), allocatable, save  :: IBOITEP(:)          !< Domain where a particle is localized
integer(kind=DPI), allocatable, save  :: NsegBOITE(:)        !< Number of segments in each domains
integer(kind=DPI), allocatable, save  :: ListeparBOITE(:,:)  !< liste des particule par boite
integer(kind=DPI), allocatable, save  :: NparBOITE(:)        !< nombre de particules par boite
integer(kind=DPI), allocatable, save  :: B3D_1D(:,:,:)       !< renvoie le numero de la boite en fonction de ces coordonnees
integer(kind=DPI), allocatable, save  :: B1D_3D(:,:)         !< renvoie les coordonnees de la boite a partir du numero
integer(kind=DPI), allocatable, save  :: IndBoiteVois(:,:)   !< tableau des 27 boites voisines a une boite B
integer(kind=DPI), allocatable, save  :: NBox2Shift(:)       !< total number of boxes to shift for a given box BI
                                                             !< at the border of the simulation box
integer(kind=DPI), allocatable, save  :: Box2Shift(:, :)     !< index of boxes to shift for a given box BI at the
                                                             !< border of the simulation box
real(kind=DP), allocatable, save      :: SIGBOX(:,:,:)       !< contrainte aucentre des boites
real(kind=DP)                         :: TailleBoite(3)      !< dimensions of the domains
logical, allocatable, save            :: Key_Box_2Border(:)  !< key is true when a box is on a border or at the
                                                             ! second position near a border

real(kind=DP)           :: dom_effective_deb(3)         ,& !< boundaries of the effective simulation volume
                           dom_effective_end(3)            ! used for the multipole domains

logical                     ::  allocation_dynamique_boites    !< Key to activate automatic allocation for the multopoles algo
integer(kind=DPI)           ::  NBoites                     ,& !< Number of domains
                                NboitesX,NboitesY,NboitesZ  ,& !< Number of domains in the X, Y, Z directions
                                NSegmax_Boite               ,& !< maximum number of segments in a domain
                                Nparmax_Boite               ,& !< maximum number of particles in a domain
                                IBOITE(NSEGMAX)             ,& !< index of the domain attached to segment I
                                LengthLimit                    !< Maximum length of a segment in "BOITE"

! Building of the special array IndexBoite(:)%ListeSeg(:) with two allocatable dimensions
type intarray2D
  integer(kind=DPI), dimension(:), allocatable  :: ListeSeg
end type intarray2D
type (intarray2D) , dimension(:), allocatable   :: IndexBoite


integer(kind=DPI) ::  IDEP(NSEGMAX)     ,&! List of displacement for each segments I during a step
                      QUIDEP(NSEGMAX)   ,&! List of moving segments during a step
                      NBDEP             ,&! Number of moving segments
                      IJONC(NSEGMAX)      ! Index of the segment binome of a junction segment

integer           ::  LoiSeg(NSEGMAX)     ! Index of the velocity law attached to each segment
                                          ! to minimize parallel computations, loiseg must be a simple integer

real(kind=DP)      :: RBOX(3)     ,&! Coordinate of the center of one multipole box
                      normsurf(3)   ! unitary plan normal (derived from Plane_Miller)

integer (kind=DPI) :: KK            ,&    !< The time step increments
                      KKjonc        ,&    !<
                      nsegm         ,&    !< The number of segments
                      nsegINI       ,&    !< The number of segments considered for the last Greengard optimization
                      plusseg       ,&    !< Number of segment to add
                      kk_debut_jonc ,&    !<
                      axe1_carto    ,&    !<
                      axe2_carto    ,&    !<
                      KKGD          ,&    !<
                      kk_debut_GD         !<

Real(DP)      :: phi1_carto, phi2_carto,dist_carto_ini
Character     :: fichier_carto*10

! The input output file name definition
character(len=40)   :: file_graph, file_stat, file_eps, file_rau, file_raujonc, file_gamma, file_travapp, file_travint,  &
                       file_gammap, file_bigsave, file_film, file_b_plan, file_particules, file_defparticules, &
                       file_fichseg, file_gammabox, file_timebreak

INTEGER  Ma_Couleur       ! Colors of process for simulations with Nd_Phase = 1 or 2

real(kind=DP) :: GB2_size        ! rayon de la barriere spherique
integer (kind=DPI) ::    &
        KRC             ,&! Periodicite du calcul des boites
        ! IMAXI         ,&! Longueur ref. = Deplacement maximal par step
        ICROSMAX        ,&! Longueur minimum pour devier
        relax_TL        ,&! nb de pas de relaxation avec tension de ligne seulement
        relax_INT       ,&! nb de pas de relaxation avec tension de ligne, sigma_INT mais sans reaction
        relax_reac      ,&! nb de pas de simulation sans chargement
        NSTEP           ,&! Nombre d'iteration total (relaxation et chargement)
        Nombre_FR       ,&! Nombre total dessource de franK_Read
        KISAUVE         ,&! Periodicite de sauvegarde du bigsave
        KKIM            ,&! Periodicite de sauvegarde du bigsave
        KPREDRAW        ,&! Periodicite des appels a l'interface graph.
        KSTAT           ,&! Periodicite de sauvegarde des mesures
        IFORET          ,&! Nombre de dislocations immobiles
        CLPNBOIMAGE     ,&! nombre d'image de l'echantillon pour les Conditions aux Limites Periodiques
        AFFNBOIMAGE     ,&! Nombre d'images affiche
        segments_fous     ! nombre de segments depassant la vitesse du sont


real(kind=DP)::       &
     accutime        ,&!TEMPS DE LA SIMULATION
     SIGAPP(3,3)     ,&! Sigma appliquee
     TensApp(3,3)    ,&! tenseur applique lu dans le ficher (../in/tensapp) quand l'axe de traction est nul
     SIGPLUS(3,3)    ,&! Sigma increment
     TENS(3,3)       ,&! ? Tens
     TauTOT(NSEGMAX) ,&! Total resolved shear stress
     TauTL (NSEGMAX) ,&! Line tension
     TauAPP(NSEGMAX) ,&! Applied resolved shear stress
     TauINT(NSEGMAX) ,&! Internal resolved shear stress
     TauMon(NSEGMAX) ,&! climb resolved shear stress
     DAIRE(NSEGMAX)  ,&! Aire Balayee par un segment par step
     !DTABV(10000,ntsg_max)  ,&! Table des solutions de la TL vis
     !DTABC(10000,ntsg_max)	,&! Table des solutions de la TL coin
     dep_moyen_R    ! deplacement moyen reel par iteration

integer (kind=DPI) ::                   dimFE_Sig(3),FEgrid_or(3),FEgrid_size(3) !< Dimensions of the FE mesh used to imposed a stress gradient
real(kind=DP), allocatable, save ::     FE_SigappS(:,:,:,:)     !< The stress tensor as defined in the FE mesh

integer (kind=DPI) :: &
     Rlocal          ,&! parameter length use to regularize the local line tension
     Rlocal2           ! RLOCAL**2

integer(kind=DPI)  :: NSEGMO,NSEGMJONC , &
     Dlimite_max  ! deplacement maximal des seg/iteration remplace sound

real(kind=DP)::      &
Dlimite(NLV_max)    ,&!< deplacement maximal cale sur les la vitesse moyenne
Temperature         ,&!< Temperature
UnkT                ,&!< 1/kT   in eV units
UnkT_Joule          ,&!< 1/kT   in Joule units
RAUDMO_I            ,&!< Densite de dislocations mobiles
DENSITEMOB          ,&!< Densite de dislocations mobiles
RAUDIS              ,&!< Densite totale de dislocations
RAUDIS_matrice      ,&!< Densite de dislocations dans la matrice
RAUDIS_precipite    ,&!< Densite de dislocations dans le precipite
RAUDISjonc          ,&!< Densite totale de dislocations
RAUDISjonc_precipite,&!< Densite totale de dislocations
RAUDISvis           ,&!< Densite totale de dislocations vis
RAUDISnonvis        ,&!< Densite totale de dislocations coin et mixte
RAUDMO              ,&!< Densite totale de dislocations mobiles
XMU                 ,&!< Module de cisaillement
VITA                ,&!< Vitesse des dislocation
VITMOY              ,&!< Vitesse moyenne des dislocations
VITMAX              ,&!< Vitesse maxi (utile dans metallofute2)
TAUIII              ,&!< Contrainte (exp.) de debut de stade III
XALPHA              ,&!< Ratio sigma devie et primaire pour le CS
ARGEXP              ,&!< Constante hors de l exp de la proba de CS
AVALUE              ,&!< Valeur du parametre de longueur de la sim.
BdivA               ,&!< Norme de Burger exprime en unite avalue
Bspread2            ,&!< Spreading radius used in short range interaction formulas
Bdivpa              ,&!< Norme de Burger exprime en unite avalue divided by pi
VOLUME              ,&!< Volume de la boite de simulation
DELTAT              ,&!< Pas de temps de la simulation
Facteur_Depmax      ,&!< rapport deplacement maximal/deplacement moyen des segments
SONDMS              ,&!< Idem SOUND, mais en metre par seconde
XLSPLIT             ,&!< Taille des boites pour le calcul des sig
XLSMOIT             ,&!< Cste => XLSPLIT / 2
XLINV               ,&!< Cste => XLINV=UN/XLSPLIT
FA                  ,&!< Cste => FA=(UN-DPOISS)
FACT                ,&!< Cste => FACT=-1*HALF/(UN-DPOISS)
Solli_Sys(NTSG_MAX) !< A prefactor used to artificially increase or decrease the Schmid factor

integer (kind=DPI) :: SIDEJA ! Clef de redemarage d'un calcul precedant
! 0 -> New simulation
! 1 -> Restart of a previous simulation
! 2 -> Simulation with precipitates

integer (kind=DPI) :: XLOMAX    ! Longueur de discretisation

integer (kind=DPI) :: GB     ! Clef d'activation des barrieres
! 0 -> inexistantes
! 1 -> planes
! 2 -> spheriques
! 3 -> pavage

integer (kind=DPI) :: LINTEN     ! Line tension control key
! 0 -> Cste energy approximation
! 1 -> Isotropic solution
! 2 -> Anisotrpic solution

real(kind=DP), dimension(3,91) :: Disd     ! Table contaning anisotripic value of line tension taken from Disdi (c.f. J. Douin)
real(kind=DP)                  :: DisdInc  ! The Disd elementary Increment

logical ::       &
GLDEV           ,&! Clef d'activation du glissement devie
Singulier(NSEGMAX)! le segemnt est au voisinage d'un pui de force

character(len=10) :: DAT ,&! Date
TIM   ! Temps

real(kind=DP):: BETA ! Contrainte resolue de Peierls

integer(kind=DPI) :: IMAXCRIS  !Longeur des segments touchant les bords de la boite

!=================================================================================
! Variables ajoutees pour la definition des donnees cristallographiques
! des systemes de glissement ------------------------Ghiath , 10/01/01

character (3)      :: crystal_structure
Integer (kind=DPI) :: Nb_slip_types       ! nombre de famille de systeme de glissement
integer  :: NbSysDev !*** nb de systemes de deviation possible effectif (varie en fonction de la cristallo
Integer (kind=DPI) :: NTSG                ! nombre total de systeme de glissement = somme de nb_slip_systems
Integer (kind=DPI) :: NLV                 ! le nombre des loie de vitesse utilise
Real    (kind=DP)  :: vitesse_unitaire    ! vitesse unitaire en m/s des segments coin
type SLIP_DATA
    integer (kind=DPI)               :: Nsystemes           ! nombre total de systemes de glissements
    integer (kind=DPI),dimension(3)  :: planes
    integer (kind=DPI),dimension(3)  :: directions
    real (kind=DPI)                  :: Vecburgers
end type SLIP_DATA
type(SLIP_DATA),dimension(NTSG_MAX) :: Slip  ! type derive des systeme de glissement

type LOI_Vitesse_brute
   integer(kind=DPI)              :: Arrhenius
   real(kind=DP)                  :: h       ! Attac frequency in an Arrhenuis law type of mobility law
   real(kind=DP)                  :: deltaG0 ! Energie d'activation totale en eV (pour T = 0 K pour les vis)
   real(kind=DP)                  :: tau0    ! tau effectif a zero K et sans oxygene (pour les vis)
   real(kind=DP)                  :: coef_p  ! exposant p de la loie de vitesse (pour les vis)
   real(kind=DP)                  :: coef_q  ! exposant q de la loie de vitesse (pour les vis)
   real(kind=DP)                  :: V0      ! ratio V(edg) / Vscrew)
   real(kind=DP)                  :: BF
   real(kind=DP)                  :: friction
end type LOI_vitesse_brute

type(LOI_vitesse_brute),dimension(NLV_max) :: Loi! type derive des systeme de glissement

integer (kind=DPI) :: Mode_deformation !< 1:sigma imposee; 2: epsilon imp.; 3 fluage;DPI(fatigue) ; 5metallofute;5(fatigue); 10: crack
#ifdef MDC
logical              :: Mode_deformation_key
#endif
integer :: cartograph     !< 0:aucun traitement ; 1: stat sur les caerto; 2 generer les ficiers
integer :: FACTEUR_CRIS   !< indique le facteur crystallographic lu
integer :: FACTEUR_BOITE  !< facteur multiplicatif des dimensions de la boite
integer :: NTGD           !< Nombre total de segments GD
integer :: NTGD1          !< Nombre total de segments GD = 1 (made by collinear annihilation)
integer :: NTGD2          !< Nombre total de segments GD = 2 (made by cross-slip)
integer :: OLDNTGD1       !< Nombre total de segments GD = 1 (made by collinear annihilation) eliminated in NET
integer :: OLDNTGD2       !< Nombre total de segments GD = 2 (made by cross-slip) eliminated in NET
integer :: gd2elimin      !< A var to remember the number of gd2 eliminated
integer :: NTFS           !< Total number of segment touching a freesurface

logical :: GD2made        !< Did we introduced a GD2 between two call to connecij

integer (kind=DPI) :: homothetie, nbjonction
integer (kind=DPI) :: x_reseau(3), y_reseau(3), z_reseau(3) ! vecteur de translations unitaire pour paver l'espace

REAL(kind=DP) ::                 &
VecBurgers                      ,& !< The Burgers vector size
deltat0                         ,& !< The simulation time step
TauINT_Limite                      !< The critical stress defining the critical limit where the dynamics cannot be solved correctly
                                   !! For such segments in a singular field we impose arbitrary IDEP = 2

REAL(DP) :: D_brute(NLV_max), &   !< The predicted mean displacement per unit length for each type of velocity law
            D_vraie(NLV_max), &   !< The effective mean displacement per unit length for each type of velocity law
            Lseg(NLV_max),    &   !< The predicted total length of mobile segments for each type of velocity law
            Lseg_V(NLV_max)       !< The effective total length of mobile segments for each type of velocity law

logical  :: MGO,ORT,CFC,DC,CS,BCC,HCP,effet_longueur,shear

real(DP) :: schmid,airevis,airecoin,airemixte
real(DP) :: FrankLoop, FrankLine, FrankFreeLine
integer  :: FreeLoopNumb, LoopNumb

!===========================================================================
!================SEGMENT VARIABLES =========================================
!===========================================================================
type SEGMENT
   real(kind=DP) ::   RESDEP,     &   !< Non Integer part of segment displacement
                      DEPINST,    &   !< Instantaneous displacement saved for wait procedure
                      probadev,   &   !< Probability of deviation used for waiting segment
                      taudev,     &   !< Resolved shear stress in the deviation plane
                      anglevis        !< The segment character, angle to the screw direction
#ifdef MDC
   real(kind=DP), dimension(6) :: SigFE ! FE stress field
#endif

    integer (kind=DPI) ,dimension(3)  :: O  !< Coordinates at segments Origin
    integer (kind=DPI) :: VECLIN, &         !< Indice of Line vector and displacement
    NORME,      & !< Segment length (Must be positive)
    VOISO,      & !< First neighboring segment from segment origin (with velocity)
    VOISE,      & !< First neighboring segment from segment extremity (with velocity)
    VNNO,       & !< First non-zero neighboring segment from segment origin (with velocity)
    VNNE,       & !< First non-zero neighboring segment from segment extremity (with velocity)
    IJONC,      & !< Index of binomial segment if junction
    TJONC,      & !< Junction life time (Number of Iteration)
    WAIT,       & !< Number of step with a displacement smaller than unity
    GRAIN,      & !< Localization of the segments in the polycrystal microstructure
    DOM,        & !< Localization of the segment when facing the problem with multi-domains (GB=5)
    SURFACE,    & !< Key to identify segments touching free surfaces (0=no,1=origin,2=end,3=both sides)
    VARFREEPLAN,& !< Index of free surface touched by the segment
    Nsysdev,    & !< The corresponding deviation system
    GD,         & !< The type of GD segments 0 = not GD, 1 = GD made by cross-slip, 2 = GD made by colli annihilation
    NPhase        !< In multi-phase materials the phase number

logical ::  JONC,     &    !< Hold key for junctions
            DISEG,    &    !< Discretization key for segment
            BLOQUER,  &    !< Hold key for a barrier
            unload,   &    !< True if a segment is artificially not affected by the applied load
            ZEROTL         !< True for segments on surface, changing from domain
end type SEGMENT

!##### Variables #######

type(SEGMENT),dimension(NSEGMAX) :: SEG ! Tableaux des segments
! ...................................................................................
! variables for long segments, incompatible with the multipole domains!
integer (DPI) :: Nb_GroSeg  ! nb of very long segment
integer (DPI) :: GroSeg(Nsegmax/idix)  ! table of long segment (NbGroSeg)

logical :: LongSeg(NSEGMAX)    ! true for segments > domaine size

!=====================================================================================
!================  The particles mode variables  =====================================
!=====================================================================================

type precipite
  real(kind=DP)                   :: tau      ! effective Stress needed to penetrate the precipitate
  integer(kind=DPI) ,dimension(3) :: C        ! Coordinates of the particles center
  integer(kind=DPI)               :: R        ! radius of the particles
end type precipite

type(precipite)  ,allocatable,save    :: par(:)                 !< The list of particles \ingroup particles

integer(kind=DPI)                     :: npar                   !< Nb of particles \ingroup particles
integer(kind=DPI), allocatable, save  :: NparTester(:)          !< Number of particles to be detected for each box \ingroup particles
integer(kind=DPI), allocatable, save  :: listeParTester(:,:)    !< list of particles to be detected for each box \ingroup particles

!=====================================================================================
!================  The creep control mode variables  =================================
!=====================================================================================

logical::             creep             !< The switch variable to define the creep loading mode \ingroup creep
logical::             timebreak         !< The switch variable to activate a climb timebreak event \ingroup creep

integer (kind=DPI)                    :: StepNextTimeBreak      !< Time step where the next climb timebreak event can take place \ingroup creep
integer(kind=DPI)                     :: npar_cut               !< Nb of cutting planes \ingroup creep
integer(kind=DPI)                     :: npar_touch             !< Nb of dislocation sections touching a particles before the cutting process \ingroup creep
integer(kind=DPI), allocatable        :: par_cut(:,:)           !< A flag tab to keep information on part cutting \ingroup creep
                                                                !! - 1dim   = the cutting plane label
                                                                !! - 2dim_1 = the particle label
                                                                !! - 2dim_2 = the cutting slip system
                                                                !! - 2dim_3 = the cutting plane coordinate
                                                                !! - 2dim_4 = number of segments found in this cutting plane
                                                                !! - 2dim_5 = number of steps with zero segments found in the cutting plane
integer(kind=DPI), allocatable        :: par_touch(:,:)         !< A flag tab to keep information on segments stopped at particles contour \ingroup creep
                                                                !! - 1dim   = the touching plane label
                                                                !! - 2dim_1 = the particle label
                                                                !! - 2dim_2 = the touching slip system
                                                                !! - 2dim_3 = the touching plane coordinate
                                                                !! - 2dim_4 = number of segments found in this touching plane
                                                                !! - 2dim_5 = number of steps with zero segments found in the touching plane

real(kind=DP)                         :: TwoPi                  !< \f$ 2. \times \pi \f$ \ingroup creep
real(kind=DP)                         :: Omega                  !< atomic volume  (m**3) \ingroup creep
real(kind=DP)                         :: prefac2                !< a variable needed for the climb velocity equation \ingroup creep
real(kind=DP)                         :: satfac                 !< saturate concentration of vacancies \ingroup creep
real(kind=DP)                         :: InvVecBurgers          !< \f$ 1 / b \f$ \ingroup creep
real(kind=DP)                         :: lnrrc                  !< logarithm term appearing in the climb velocity equation \ingroup creep

real (kind=DP),    allocatable        :: h_par_touch(:)         !< The remaining climb distance for segments listed in par_touch(:,:) to bypass their blocking particle \ingroup creep
real (kind=DP),    allocatable        :: angle_par_touch(:)     !< The angle to edge for segments in a par_touch(:,:) \ingroup creep
real (kind=DP),    allocatable        :: tau_par_touch(:)       !< The resolved climb stress on the climbing segment in a par_touch(:,:) \ingroup creep

!========================================================================================

real(kind=DP) :: DPOISS   ! Coeficient de Poisson

real(kind=DP) :: AB_LENGTH_SEUIL             !longueur significative
! Sert  a zapper des segment ds la boucle la plus interne du programme ***

Real(kind=DP) :: DCFJ ! Distance pour le calcul de force sur les bras de jonction

#if defined(PA) || defined(MDC)
!====================================================
!=========    Variables pour le parallelisme  =======
!====================================================
INCLUDE 'mpif.h'
INTEGER  IERR                 ! retour des fonctions du MPI; 1 en cas d'erreur, 0 RAS
#endif

#ifdef PA
!integer                     :: MPI_COMM_WORLD
INTEGER  TAILLE               ! Processes number in MPI_COM_WORLD
INTEGER  TAILLE_PAIRIMPAIR    ! Processes number in COMM_PAIRIMPAIR
INTEGER  Mon_Rang             ! Rank of Process in MPI_COMM_WORD
INTEGER  Mon_Rang_PairImpair  ! Rank of Process in MPI_COMM_WORD
INTEGER  COMM_PAIRIMPAIR      ! The Pair-Impair Process Communicator used when Nb_Phase =2
INTEGER  STATUS(MPI_STATUS_SIZE)

!Calcul des contraintes
REAL (KIND=DP)  VITSIGINT   !la vitesse du process au pas n+1
REAL (KIND=DP)  VITSIGLP    !la vitesse du process au pas n+1 pour le calcul LP

!Repartition des segments
LOGICAL MONSEG(NSEGMAX)

#endif

!================================================================
!==============variables for the internal planes=================
!================================================================

integer(DPI) :: NbPlan,NbplanDom,NbPlanMax    ! Numbers of different type of planes
logical :: InclExcl

type plans
    Integer(kind=DPI) :: Miller(3)
    Integer(kind=DPI) :: Dom
    real(kind=DP)     :: Pos
    Integer(kind=DPI) :: Knd  !boundary type
end type plans

!new definition of planes!! derived type was slow better to define arrays...
integer(kind=DPI), allocatable :: Plane_MillerI(:,:)
real(kind=DP), allocatable     :: Plane_MillerR(:,:)
real(kind=DP), allocatable     :: Plane_pos(:)
integer(kind=DPI), allocatable :: Plane_knd(:)
integer(kind=DPI), allocatable :: Plane_dom(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


real(kind=DP)    :: TAU_APB,SIG_APB

type(plans),dimension(100) :: Varplan

!================================================================
!==============variables for the free surface planes=============
!================================================================

integer (DPI) :: NbFreePlan            !Number of free planes to consider in the simulation

!================================================================
!==============variables for Concave domains=============
!================================================================

integer(DPI)              :: NbCvxDom   !Number of convex domain used to compose a cancave geometry or obstacles
integer(DPI),allocatable  :: NbPlanCvxDom(:,:),ListCvxDom(:),ListCvxDomSegi(:)  !Number of plains belonging to each convex object
Logical,allocatable  :: DomConnec(:,:) !

!==================================================================
!=============== Variables needed for polycrystal computation =====
!=============== based on a rotation of the applied stress in =====
!=============== each grain.                                  =====
!==================================================================

Integer(kind=DPI)             ::  NbGrains,                 & !< Number of grains
                                  DimSyst,                  & !<
                                  longminseg                  !< The minimum length of segments touching grain boundaries
Integer(kind=DPI),allocatable ::  PlanMiller(:,:,:)
logical                       ::  DesorientGrain              !< The logical key used to activate the grain disorientation in multi-grains computations
Real(kind=DP),allocatable     ::  SigApporiente(:,:,:),     & !< The stress tensor applied to each grains (after rotation)
                                  MatRotGrain(:,:,:),       & !< The rotation matrix associated to each grains
                                  RauGrain(:),              & !< Dislocation density per grains
                                  RauGrainJonc(:),          & !< Junction density per grains
                                  RauGrainSys(:,:),         & !< Dislocation density per slip system
                                  RauSysJoncGrain(:,:),     & !< Junction density per slip systems
                                  TensRotGrain(:,:,:),      & !< Elementary rotation tensor for each grain
                                  TensRotSysGrain(:,:,:,:), & !< Elementary rotation tensor for each slip system (Wij)for each grain
                                  GamSysgrain(:,:),         & !< Gamma per grain and per slip systems
                                  AireSysGrain(:,:),        & !< slipped area per grain and per slip systems
                                  SchmidSysGrain(:,:),      & !< Schmid factor per grain and per slip systems
                                  PlanPos(:,:)                !<

!=============================================================================================
!===============variables for the long-rang replica  with the Greengard method ===============
!=============================================================================================

integer (kind=DPI)             :: PBCIxDIM,PBCIyDIM,PBCIzDIM,MAXMODUR
logical                        :: PBCIstatus

!========================================================================
!===============variables for the multi-phase simulations ===============
!========================================================================

integer (kind=DPI)             :: Nb_phase, Index_phase  ! number of phases used, index of the phase
                                                         ! Present even in non parallel mode on the material input file

!===============================================================================
!===============variables for passing arguments to mm executable ===============
!===============================================================================


logical::domvtx=.FALSE.
logical::initconf=.FALSE.




#ifdef PA
integer(kind=DPI),allocatable  :: liste_RO(:,:),      & ! list of the boxes and segment centers close to the boundary
                                  liste_RO_Recv(:,:)    ! equivalent to liste_RO, but received from the other phase

integer(kind=DPI),allocatable  :: list_boxes(:)         ! list of boxes close to the boundary

integer(kind=DPI)  :: NBoites_boundary,       & ! number of boxes close to the boundary
                      Nb_seg                    ! number of segments in the total box list close to the boundary
#endif

#ifdef MDC
!==================================================================
!=============== Variables for the coupling DD-EF =================
!==================================================================

!*** MPI include and variable

logical                     :: zebmpiinitialized=.FALSE.
integer                     :: zebmpierr
integer                     :: zebmpistatus(MPI_STATUS_SIZE)
integer                     :: zebmpicomm
character, pointer          :: zebmpisendbuf(:)
integer                     :: zebmpisendbufpos
integer                     :: zebmpisendbufsize=0
character, pointer          :: zebmpirecvbuf(:)
integer                     :: zebmpirecvbufpos
integer                     :: zebmpirecvbufsize=0

INTEGER              :: ZebulonTid=0    !AR- TaskID of Zebulon
INTEGER              :: MAGIC_DISLO=666 !AR- Magic number used for sending acknowledgements
INTEGER              :: StringLength    !AR- Length of a string to be transferred

!****

integer(kind=DPI)                                  :: halfthickness
real(kind=DP)                                      :: halfthickness2,SR_cutoff,SR_cutoff2
integer                                            :: mdc_timestep_ratio
logical                                            :: zcall
integer                                            :: nbswept
integer(kind=DPI)                                  :: nsegm_INI
integer,dimension(:),allocatable                   :: sweptsurfdata

integer,dimension(:,:,:),allocatable                 :: fsurfext !array to compute surface extesion
real,dimension(:,:),allocatable                        :: fsurfext_norm !norm


#endif


end module VARGLOB
!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
