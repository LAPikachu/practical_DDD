!===================================================================================================
!========================    DEBUT    MODULE  "CONSTANTES"  ========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains the definition of all the constants used in the microMegas.
!> \todo Some comments in this module have to be more specific and/or translated in English
module CONSTANTES

implicit none

integer,parameter  :: DP=selected_real_kind(p=14) !< Size for the variables with real type
integer,parameter  :: DPI=selected_int_kind(13)   !< Size for the variables with integer type
integer,parameter  :: DPI_S=selected_int_kind(9)  !< Size for the variables with simple integer type
integer,parameter  :: NLV_max=10                  !< nombre total de lois de vitesse utilisees
integer,parameter  :: NTSG_MAX=25                 !< nombre maximal de systemes de glissement
integer(kind=DPI),parameter  :: NSEGMAX=100000    !< dimension des tableaux
integer(kind=DPI),parameter  :: GNSEGMAX=100000   !< dimension des tableaux avec images
integer,parameter  :: FACTEUR_CS=2                !< facteur d'echelle realle  ---> Base de vecteur
integer,parameter  :: FACTEUR_CFC=4               !< facteur d'echelle realle  ---> Base de vecteur
integer,parameter  :: FACTEUR_DC=4                !< facteur d'echelle realle  ---> Base de vecteur
integer,parameter  :: FACTEUR_BCC=12              !< facteur d'echelle realle  ---> Base de vecteur
integer,parameter  :: FACTEUR_HC=132              !< facteur d'echelle realle  ---> Base de vecteur
integer,parameter  :: FACTEUR_ORT=2               !< facteur d'echelle realle  ---> Base de vecteur 6222
integer,parameter  :: FACTEUR_MGO=12              !< facteur d'echelle realle  ---> Base de vecteur
integer,parameter  :: Facteur_boite_CFC = 8       !< facteur multiplicatif de la boite de simu
integer,parameter  :: Facteur_boite_DC = 8        !< facteur multiplicatif de la boite de simu
integer,parameter  :: Facteur_boite_CS = 8        !< facteur multiplicatif de la boite de simu
integer,parameter  :: Facteur_boite_BCC = 24      !< facteur multiplicatif de la boite de simu
integer,parameter  :: Facteur_boite_HC = 1584     !< facteur multiplicatif de la boite de simu
integer,parameter  :: Facteur_boite_ORT = 8       !< facteur multiplicatif de la boite de simu 31110
integer,parameter  :: Facteur_boite_MGO = 6       !< facteur multiplicatif de la boite de simu
integer,parameter  :: Nb_max_sys_dev = 2          !< maximum Number of collinear interactions all crystallographies

integer,parameter  :: ListeSegRangeIni = 200      !< The initial segment range used in the ListeSeg tab dim optimization
real,parameter     :: ListeSegRangeFac = 0.2      !< An increase factor used in the ListeSeg tab dim optimization algo
                                                  !< the default value is 20%

real(kind=dp) , parameter :: unpeu = 0.3d-6    !< Limit for stress in interaction calculations

real          , parameter :: numtols = 1.e-6   !< Machine simple precision
real(kind=dp) , parameter :: numtold = 1.d-12  !< Machine double precision
real(kind=dp) , parameter :: numtol_dotp = 1.d-5  !< Machine double precision
!< precision servant de comparaison entre les reelles.
!! Attention : ce parametre est tres important il sert a comparer les nombres reels entre eux.
!! Le cas le plus frequent est un codage en 32 bits ce qui fournit, en principe, une
!! resolution superieure a 1E-9. De toute facon la precision doit etre superieur a 2puissance15,
!! qui represente la resolution entiere (codage en 16 bits)

real, parameter :: Boltzmann = 0.0000868      !< Boltzmann constant \f$ k_B \f$ (eV/K)
real, parameter :: unev = 1.602176565E-19     !< 1 electron Volt \f$ e_V \f$ in Joule units (J)
real, parameter :: Debye = 8.8E13             !< Debye Frequency in s-1
real, parameter :: Angle_vis_CFC = 2.0        !< Error in degree on the screw dislocation character definition
real, parameter :: Angle_vis_DC  = 2.0
real, parameter :: Angle_vis_HC  = 2.0
real, parameter :: Angle_vis_BCC = 0.5
real, parameter :: Angle_vis_ORT = 1.0
real, parameter :: Angle_vis_MGO = 2.0
real, parameter :: Ldedou = 1.0E-6            !< distance en m de dedoublement de sourcesquand cartograph est mis = 5
real, parameter :: Effica_crans = 0.0         !< Ratio d efficacite des cran (compris entre 0 et 1)
real, parameter :: fac_rlocal = 0.5           !< The prefector applied to the rlocal cste (see local line tension)
real, parameter :: spread_core_radius = 0.5d0 !< core radius spreading factor used in short range stress calculation

integer,dimension(4),parameter  :: IPERM=(/2,3,1,2/)

!#####################################################
!# Parameters related to the discretization process  #
!#####################################################

integer,parameter  :: NBASERED = 8      !< Number of vectors (directions) defining dislocation character per slip system

integer,parameter  :: NBSysDev_CFC = 1  !<*** nb de systemes de deviation possible
integer,parameter  :: NBSysDev_DC  = 1  !<*** nb de systemes de deviation possible
integer,parameter  :: NBSysDev_BCC = 2  !<*** nb de systemes de deviation possible
integer,parameter  :: NBSysDev_HCP = 3  !<*** nb de systemes de deviation possible
integer,parameter  :: NBSysDev_CS  = 1  !<*** nb de systemes de deviation possible
integer,parameter  :: NBSysDev_MGO = 1  !<*** nb de systemes de deviation possible
integer,parameter  :: NBSysDev_ORT = 1  !<*** nb de systemes de deviation possible

real   ,parameter  :: Microboucle_ratio = 1.0 !< Parameter defining the length of loops to be eliminate during the dynamics
                                              !! This length is a priori proportional to the discretization length (xlomax)

integer,parameter  :: Npar_max   = 40000  !<*** dimension de la base reduite
integer,parameter  :: NbItJseuil = 10
integer,parameter  :: tjoncmax   = 3      !<*** The germination time for junctions

integer (kind=DPI),parameter :: IZERO   = 0   ,&
                                IUN     = 1   ,&
                                IDEUX   = 2   ,&
                                ITROIS  = 3   ,&
                                IQUATRE = 4   ,&
                                ICINQ   = 5   ,&
                                ISIX    = 6   ,&
                                ISEPT   = 7   ,&
                                IHUIT   = 8   ,&
                                INEUF   = 9   ,&
                                IDIX    = 10  ,&
                                IONZE   = 11  ,&
                                IDOUZE   = 12  ,&
                                ITREIZE  = 13  ,&
                                ICENT   = 100

real(kind=DP), parameter ::     ZERO    = 0.0D0                 ,&
                                SQRT2   = 1.41421356D0          ,&
                                SQRT3   = 1.73205081D0          ,&
                                UNQSQ3  = 0.14433757D0          ,&
                                INVDIX  = 0.1D0                 ,&
                                QUART   = 0.25D0                ,&
                                TIERS   = 0.33333333D0          ,&
                                HALF    = 0.5D0                 ,&
                                PERCEN  = 0.6D0                 ,&
                                UN      = 1.0D0                 ,&
                                DEUX    = 2.0D0                 ,&
                                TROIS   = 3.0D0                 ,&
                                QUATRE  = 4.0D0                 ,&
                                CINQ    = 5.0D0                 ,&
                                DIX     = 10.0D0                ,&
                                CENT    = 100.0D0               ,&
                                MILLION = 1.0D6                 ,&
                                PII     = QUATRE * DATAN(UN)    ,&
                                HALFPII = (HALF * PII)

real(kind=DP), dimension(3,3), parameter  :: MatIdent = RESHAPE((/1.,0.,0.,0.,1.,0.,0.,0.,1./),SHAPE=(/3,3/)) !< Identity 3*3 Matrix

integer (kind=DPI), parameter :: RIGIDBODY = 2  !< total volume rotation key for the loading Calculs
                                                !! 1-> no body rotation
                                                !! 2 -> body rotation (the default solution)

logical                     ::  calculate_gammabox = .False.    !< Key to write the gammabox file

end module constantes

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
