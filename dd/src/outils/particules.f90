
!===================================================================================================
!========================    DEBUT    MODULE     ===================================================
!===================================================================================================



!#########################################################################
!#                                                                       #
!# Module INIT : flux de donnees entrant et initialisation des           #
!#                                                                       #
!#  structures de donnees                                                #
!#                                                                       #
!# Function noeud(Oi) : fournit le noeud du reseau le plus proche de Oi  #
!#                                                                       #
!# subroutine lecture_donnees : donnes materiaux et de control           #
!#                                                                       #
!# subroutine parametrage : determine tous les parametres de la simulation
!#                          sauf les parametres relies a la boit de      #
!#                          simulation et au segments initiaux           #
!#                                                                       #
!# subroutine initialiser_fichiers : initialisation des fichiers sortis  #
!#                                                                       #
!# subroutine Allocation : Dimensionnement des tableaux de la base       #
!#                                                                       #
!# subroutine generer_basedevecteurs : c'est la procedure qui generer    #
!#            la base de vecteur de discretisation, verifie la compatib. #
!#            en fonction de la cristallographie en question             #
!#                                                                       #
!# subroutine initialiser : initialisation des tableaux et variables     #
!#                                                                       #
!# subroutine Tabulation : Construction des tableaux utils pour la       #
!#                         topologie de la simulation                    #
!#                                                                       #
!#                                                                       #
!# subroutine lire_segments : lit le fichier des segemnts initial, et si #
!#            le nombre de segemnts lu vaut 0 alors il lit les infos     #
!#            dans le fichier jonction pour la generation de configura.  #
!#            typique pour l'etude de jonction                           #
!#                                                                       #
!# subroutine configuration : determine tous les parametre associe aux   #
!#            sources de dislocation initiales, boite, dicretisation,etc #
!#                                                                       #
!# subroutine write_basedevecteurs : ecrit de fichier BVD.? contenant les#
!#            vecteurs generer de discretisation                         #
!#                                                                       #
!# subroutine write_tableaux : ecrite quelque tableux pour information   #
!#                                                                       #
!# subroutine LOAD : le sous-programme appele dans le main               #
!#                                                                       #
!#                                                                       #
!# subroutine TABTLIN : Tabulation du calcul de la tension de ligne      #
!#                                                                       #
!# subroutine Desallocation                                              #
!#                                                                       #
!############################################################# 28/1/02 ###

module INIT

use VARGLOB
use BRICAMAT
use DEBUG
use cartographie

implicit none

type LOI_Vitesse_brute
    integer                        :: Arrhenius
    real(kind=DP)                  :: V0  ! constante de la loi de vitesse d'Arrhenius en m/s
    real(kind=DP)                  :: L0  ! longueur de reference des vis (pour la loi de vitesse)
    real(kind=DP)                  :: deltaG0 ! Energie d'activation totale en eV (pour T = 0 K pour les vis)
    real(kind=DP)                  :: tau0eff !tau effectif a zero K et sans oxygene (pour les vis)
    real(kind=DP)                  :: tauath!conrainte au plateau athermique
    real(kind=DP)                  :: coef_p  ! exposant p de la loie de vitesse (pour les vis)
    real(kind=DP)                  :: coef_q  ! exposant q de la loie de vitesse (pour les vis)
    real(kind=DP)                  :: coef_visqueux
    real(kind=DP)                  :: max_friction
end type LOI_vitesse_brute
type(LOI_vitesse_brute),dimension(NLV_max)	:: Loi_brute	! type derive des systeme de glissement


real(DP):: echelle,Ldis_act,Ldis_nact,Xlomax_nact,dmusurdT,moduleG0
integer (DPI) :: Nloi(NTSG_max,nbasered/2)
character (len=60) :: fichier,fichbase,fichtab,fichseg,materiau, segments, control
real  (DP)   :: facteurcoin(10),facteurvis(10),facteurmixt(10),DisReacloc



!#######
contains
!#######



!##########################################################################
!# subroutine noeud(o) : renvoie le point le plus proche de O situe sur  #
!# le reseau principal de la simulation                                   #
!###########################################################14/04/99#######
Function noeud(Oi,homo)

integer (kind=DPI) :: Oi(3),I,a1(3),a2(3),c(3),noeud(3),O(3)
integer (kind=DPI) :: coefa1,coefa2,coefc,homo
real(DP) :: c2


if (CS .or. CFC) then
   noeud(:) = (Oi(:) / Facteur_boite) * Facteur_boite
   return
endif

if (HC) then
! attention cde calcul n'est valable que pour les vecteurs explicites a1, a2, c
   a1(1:3) = homo * x_reseau(1:3) ! unite de translation dans la direction coin
   a2(1:3) = homo * y_reseau(1:3) ! unite de translation dans la direction vis
   c(1:3)  = homo * z_reseau(1:3) ! unite de translation dans la direction
! de la normal au plan de glissement ! Ici on respecte egalement le meme facteur
! d'homothetie enre l'espace reel et de de simulation.
   c2 = inorivect(c)
   coefc = int(real(iprodsca(Oi,c),DP) / c2,DPI)

! d'abord : conformite des dimensions de la boite au reseau principal
   O(1:3) = Oi(1:3) - coefc * c(1:3) ! la projection de Oi sur le plan (111)
! donc O() = coefa1 * a1 + coefa2 * a2 =
! = (coefa2*a2(1) ;  coefa1*a1(2) ; coefa1*a1(3))
   coefa1 = int(real(O(2),DP)/real(a1(2),DP),DPI)
   coefa2 = int(real(O(1),DP)/real(a2(1),DP),DPI)
   noeud(1) = coefa1*a1(1) + coefa2*a2(1) + coefc*c(1)
   noeud(2) = coefa1*a1(2) + coefa2*a2(2) + coefc*c(2)
   noeud(3) = coefa1*a1(3) + coefa2*a2(3) + coefc*c(3)
   return
endif
if (BCC) then
   a1(1:3) = homo * x_reseau(1:3) ! unite de translation dans la direction vis1 du sys 1
   a2(1:3) = homo * y_reseau(1:3) ! unite de translation dans la direction mixt2 du sys1
   c(1:3)  = homo * z_reseau(1:3) ! unite de translation dans la direction vis1 du sys7

! C est normal au plan de glissement du sys1 ! Ici on respecte egalement le meme facteur
! d'homothetie enre l'espace reel et de de simulation.
   coefa1 = nint(HALF*real(Oi(3)-Oi(2)) / real(a1(1),DP))
   coefa2 = nint(HALF*real(Oi(3)-Oi(1)) / real(a1(1),DP))
   coefc = nint(HALF*real(Oi(1)+Oi(2)) / real(a1(1),DP))
! fournir noeud(:) : le point le plus proche de Oi
   noeud(1) = coefa1*a1(1) + coefa2*a2(1) + coefc*c(1)
   noeud(2) = coefa1*a1(2) + coefa2*a2(2) + coefc*c(2)
   noeud(3) = coefa1*a1(3) + coefa2*a2(3) + coefc*c(3)
   return
endif
RETURN
end function noeud


!###########################################################################
! Ce sous_ programme lit toutes les donnees du programme principal
!
! Ces donnees sont separees en 3 groupes correspondants a 3 fichiers
!
!  1)  Fcihier 'materiaux' : Donnees physiques propres au materiau
!  2)  Fichier 'control'   : Donnees et parametres de la simulation
!  3)  Fichier 'seg3D'     : regroupe les caracteristiques des segments
!                             presents au demarage de la simulation
!
!   Ghiath Mohamed 10/01/01
!###########################################################################
subroutine lire_donnees

implicit none

logical :: loi_relax
!# Declarations

real(kind=DP)       :: valeur
Integer (Kind=DPI)  :: i,itemp,jtemp,type_index,ios
character(len=3)    :: Xr*3
character(len=1)    :: carac
character(len=60)   :: ligne

! Initialisation des loi de vitesse brute
Loi_brute(1:NLV_MAX)%Arrhenius = 0
Loi_brute(1:NLV_MAX)%V0        = zero
Loi_brute(1:NLV_MAX)% L0    = zero
Loi_brute(1:NLV_MAX)% deltaG0    = zero
Loi_brute(1:NLV_MAX)% tau0eff    = zero
Loi_brute(1:NLV_MAX)% tauath   = zero
Loi_brute(1:NLV_MAX)% coef_p    = zero
Loi_brute(1:NLV_MAX)% coef_q    = zero
Loi_brute(1:NLV_MAX)% coef_visqueux   = zero
Loi_brute(1:NLV_MAX)% max_friction   = zero


!*** Pour eviter ce type de message etrange sur les DEC :
!    Unaligned access pid=8803 <Pm_Dec_GX>
!     va=0x14272e19c pc=0x3ff81a091b0 ra=0x3ff81a091a8 inst=0x9c090000

!#########################################################################
!# Chargement des donnees propres au materiau : Fichier 'materiau'       #
!#########################################################################
open(1,file='../in/input.dd',STATUS='OLD')
carac = " "
do while (carac /= "-")
    read (1,*) carac
enddo
read (1,*) materiau
read (1,*) control
read (1,*) segments
close (1)

open(1,file="../in/"//materiau,STATUS='OLD')

read (1,*) VecBurgers        ! vecteur Buergers en Angstrom

!**   Module de cisaillement MU a T = 0 K
read(1,*) ModuleG0

!**   pente de la courbe Module de cisaillement = fonction (Temperature)
read(1,*) dmusurdT

!*** Module de POISSON
read (1,*) DPOISS

!**  tauIII = contrainte de debut du stade III (restauration)
read(1,*) TAUIII

!**  lecture de la structure cristallographique
read(1,*) crystal_structure(1:3)

!***  Coefficients de friction visqueux BF=5.E-5(Pa.s)/(mu=42GPa)
!**   Friction seche (Peierls) TAU0*b=b*1.0D0(MPa)/(mu=42GPA)
read (1,*) Nb_slip_types

do itemp = 1, Nb_slip_types
    read (1,*) Slip(itemp)%planes(1:3)
    read (1,*) Slip(itemp)%directions(1:3)

    if (etat( Slip(itemp)%planes(1:3),Slip(itemp)%directions(1:3)) /= 0) &
    stop " Erreur de donne crystallo : vis n'est pas perpendiculaire au plan"
    read (1,*) Slip(itemp)%Nsystemes   ! nb de systemes de glissement

    read (1,*)   &
    Nloi(itemp,1),  & ! numero de la loi de vitesse associe au vis
    Nloi(itemp,2),  & ! numero de la loi de vitesse associe au mixte 1
    Nloi(itemp,3),  & ! numero de la loi de vitesse associe au coin
    Nloi(itemp,4)   ! numero de la loi de vitesse associe au mixte 2
    if (Nloi(itemp,3)-Nloi(itemp,2) /= 0 .or. Nloi(itemp,4)-Nloi(itemp,2) /= 0) then
       print *, " Le code ne permet pas de gerer des lois de vitesses differentes"
       print *, " pour les coins, mixte1 et mixte2 .... "
       print *, " Veuillez mettre ajour le calcul du tableau anglevis dans force...., By"
       stop
    endif
ENDDO

! lecture du nombre de lois de vitesse a lire
read (1,*) NLV
loi_relax = .false.
do i = 1, NLV
    read(1,*) itemp  ! numero de la loi
    if (itemp == 0 ) then
       loi_relax = .true.
       itemp = NLV_MAX
    endif
    if (i > NLV) stop "erreur de lecture : n de la loi de vitesse faux"
    read(1,*) Loi_brute(itemp)%Arrhenius  ! clef d'activation thermique
!    print *, "  )%Arrhenius =",Loi_brute(itemp)%Arrhenius

! decrire la friction du reseau
    if (Loi_brute(itemp)%arrhenius == 2) then
       read(1,*) valeur
       Loi_brute(itemp)%tau0eff = valeur! tau effectif a T= 0 K sans oxygene
       read(1,*) valeur
       Loi_brute(itemp)%tauath = valeur! contyrainte au plateau athermique
       read(1,*) valeur
       Loi_brute(itemp)% V0 = valeur      ! constante de la loi de vitesse
       read(1,*) valeur
       Loi_brute(itemp)% L0 = valeur      ! longueur de reference des vis
       read(1,*) valeur
       Loi_brute(itemp)% deltag0 = valeur ! energie d'activation totale ( pour T= 0 K)
       read(1,*) valeur
       Loi_brute(itemp)% coef_p = valeur ! exposant p de la loie de vitesse (pour les vis)
       read(1,*) valeur
       Loi_brute(itemp)% coef_q = valeur ! exposant q de la loie de vitesse (pour les vis)
    elseif (Loi_brute(itemp)%arrhenius == 1) then
       read(1,*) valeur
       Loi_brute(itemp)%tau0eff = valeur! tau effectif a T= 0 K sans oxygene
       read(1,*) valeur
       Loi_brute(itemp)%tauath = valeur! contyrainte au plateau athermique
       read(1,*) valeur
       Loi_brute(itemp)% V0 = valeur      ! constante de la loi de vitesse
       read(1,*) valeur
       Loi_brute(itemp)% deltag0 = valeur ! energie d'activation totale ( pour T= 0 K)
       read(1,*) valeur
       Loi_brute(itemp)% coef_p = valeur ! exposant p de la loie de vitesse (pour les vis)
       read(1,*) valeur
       Loi_brute(itemp)% coef_q = valeur ! exposant q de la loie de vitesse (pour les vis)
    elseif (Loi_brute(itemp)%arrhenius == 0) then
       read(1,*) valeur
       Loi_brute(itemp)%coef_visqueux = valeur
       read(1,*) valeur
       Loi_brute(itemp)%max_friction = valeur ! donnees pour les coin
    else
       stop " type de loi de vitesse inconnu "
    endif

    if (.not. loi_relax) then
       Loi_brute(NLV_MAX) =  Loi_brute(1)
       print *, " Pas de loi de relaxation detectee."
       print *, " Si une relaxation est prevue, c'est la loi numero 1"
       print *, " qui sera utilisee."
     endif
End do
close (1)

!########################################################INIT
!# Chargement des parametres de control de la simulation
!########################################################
open(2,file="../in/"//control,STATUS='OLD')

read(2,*) Mode_deformation  ! mode de deformation dans la simulation
read(2,*) cartograph  ! clef pour la generation de carto de reaction locales
if (mode_deformation > 5 .or. mode_deformation < 1) stop " mode de deformation inconnu"
read (2,*) echelle     ! longueur de reference de la simulation exprime en b

! Valeurs Input relatives au mode de deformation selectionne
read(2,*) Sigma0    ! contrainte initial imposee directement sans iteration en MPa
read(2,*) DELTAT0 !*** Pas de temps
read(2,*) pasvariable !*** clef pour ajuster ou non le pas de temps par iteration
read(2,*) depmoyen  ! deplacement moyen a imposer pour le calage du pas de temps
read(2,*) SigmaPoint        ! increment de contrainte par iteration en KPa
read(2,*) EpsilonPoint      ! vitesse de deformation
read(2,*) DeltaEpsilon      ! interval ou on calcul la moyenne de deformation

!*** Options sur les mecanismes a prendre en compte
read(2,*) LINTEN 		!	tension de ligne seulement (!! Plus actif pour le moment)
read(2,*) GLDEV			!	GLDEV	: glissement devie prit en compte ?
read(2,*) LTSIMPLE		!	LTSIMPLE	: quel type de tension de ligne  (isotrope ou pas)

! Attention LINTEN et GLDEV ne peuvent etre actif en meme temps
if (LINTEN .and. GLDEV) then
  Write(*,*) "attention les clefs LINTEN et GLDEV ne peuvent etre active en meme temps"
  stop
endif


!*** Axe de traction dans le repere crystallographique
read (2,*) Z(1:3)

!*** parametre physiques de la simulation Module d'elasticite apparent
read(2,*) RAID   !Module d'elasticite apparent
read(2,*) TEMPERATURE   ! temperature en Kelvin

!*** rapport entre deplacement maximal autorise/deplacement moyen par iteration
read(2,*) Facteur_Depmax

!*** parametres propres aux simulations
read(2,*) relax_TL  !    nb de pas de relaxation avec tension de ligne seulement
read(2,*) relax_INT  !  nb de pas de relaxation avec tension de ligne, sigma_INT mais sans reaction
read(2,*) relax_reac !  nb de pas de simulation sans chargement
read(2,*) NSTEP      !  nb maximal de pas de simulation avec chargement apres diverses relaxations
read(2,*) Ldis_act !*** longuer de discretisation en microns
read(2,*) Ldis_nact !*** longuer de discretisation en microns
!*** Periodicite du renouvellement des interactions a longues portees
read(2,*) KRC
!*** Nombre de tranches pour la methode des boites
read (2,*) ISPLIT
!*** Nombre d'image de l'echantillon pour les Conditions aux Limites Periodiques
read (2,*) CLPNBOIMAGE
!*** Nombre d'image de l'echantillon pour les Conditions aux Limites Periodiques
read (2,*) AFFNBOIMAGE

!*** Longueur de segment ajoute pour les problemes de conditions aux limites
read (2,*) IMAXCRIS
!*** Longueur significative des segments pour la methode des boites
read(2,*) AB_LENGTH_SEUIL
read(2,*) DCFJ        ! Distance pour le calcul de force sur les bras de jonction

! pour le mode FATIGUE : 5
read(2,*) EPSMAX ! Amplitude maximale de deformation(fatigue)
read(2,*) FSIG   ! Signe initial (fatigue)

! parametre du Glissement devie
read(2,*) XALPHA
read(2,*) BETA

! clef obsoletes
read(2,*) SIDEJA			!	reprise d'une ancienne simulation ?
read(2,*) GB                ! Activer des murs infranchissable dans la simulation


read(2,*) Disreacloc  ! Distacne servant pour le calcul de TauINT_limite

!*** Periodicite des sauvegardes d'une configuration
read(2,*) KISAUVE   ! commentaire
!*** Periodictie des calculs statistique
read(2,*) KSTAT
!*** Periodictie du film
read(2,*) KKIM
!*** Periodicite des appels au programme gaphique de Marc
read(2,*) KPREDRAW
print *, " 6",KPREDRAW

! variables de debuggage
read(2,*) iterinfo
close(2)


end subroutine lire_donnees



!###########################################################################
! Ce sous_ programme lit le tenseur de sollicitation lorse que l'axe de traction
! est nul. Il s'agit du cas ou la sillicitation n'est pas uniaxiale
! Fichiera lire  'tensapp' ecrire les trois lignes du tenseur dans
!                          trois ligne du fichier
!
!   Ghiath Mohamed 10/01/01
!###########################################################################
subroutine lire_tenseur

implicit none

open(2,file="../in/tensapp_def",STATUS='OLD')

read(2,*) tensapp(1,1:3)
read(2,*) tensapp(2,1:3)
read(2,*) tensapp(3,1:3)

close(2)

end subroutine lire_tenseur


!#############################################################################
!##############determine les parametre principaux de la simulation   #########
!#############################################################################
subroutine Parametrage

Implicit NONE

INTEGER(DPI) :: i
real(DP) :: resolution, maille_mat
! initialisation de la vrai variable de pas de temps elementaire

print *, "Precision REAL : ", DP
print *, "Precision INTEGER : ", DPI
if(sideja == 0) then
   print *, "Nb iterations de relaxation (tension de ligne)            : ", Relax_TL
   print *, "Nb iterations de relaxation (Tension de ligne,sigma_Int)  : ", Relax_INT
   print *, "Nb iterations de relaxation sans chargement               : ", relax_Reac
   print *, "Nombre d'iterations avec chargement externe               : ", Nstep
   Relax_INT = Relax_INT + Relax_TL
   Relax_Reac = relax_Reac + Relax_INT
   NStep = NStep + Relax_Reac
   print *, "Nombre TOTAL d'iterations de la simulation                : ", Nstep
endif
Deltat = Deltat0

! =======================================================================
!  Calcule du module de cisaillement a la temperature de simulation TEMP

XMU = ModuleG0 + dmusurdT * TEMPERATURE      ! module en GPa
write (*,'(/," Module de cisaillement = ",F5.2," GPa")') XMU

XMU = XMU * 1.D9                       ! module en Pa

!==========================================================================================
!====== Attention : ici on definit toutes variable cles du programme ====================
!==========================================================================================
! =======================================================================
! Definition de tous les parametres relies a la cristallographie
! FACTEUR_CRIS le rapport entre le modul du vecteur vis (ex. [440] et le module du
! vecteur de Burgers dans le repere cristallo(ex. 1/2 [110])

! Initilisation des clefs choix de la cristallo
CS = .false. ; BCC = .false. ; CFC = .false. ; HC = .false. ; ORT = .false.
If (crystal_structure(1:2) == 'CS') then
   CS = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_CS'
   facteur_cris = FACTEUR_CS
   fichbase = '../out/base.cs'
   fichtab  = '../out/jonction.cs'
   facteur_boite = facteur_boite_CS
   angle_vis = angle_vis_HC
elseif (crystal_structure(1:3) == 'CFC')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_CFC'
   CFC = .true.
   facteur_cris = FACTEUR_CFC
   fichbase = '../out/base.cfc'
   fichtab  = '../out/jonction.cfc'
   facteur_boite = facteur_boite_CFC
   angle_vis = angle_vis_CFC
elseIf (crystal_structure(1:3) == 'BCC')   then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_BCC'
   BCC = .true.
   facteur_cris =  FACTEUR_BCC
   fichbase = '../out/base.bcc'
   fichtab  = '../out/jonction.bcc'
   facteur_boite = facteur_boite_BCC
   angle_vis = angle_vis_BCC
!   stop " facteur boite indetermine pour les BCC"
elseIf (crystal_structure(1:3) == 'HCP')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_HC'
   HC = .true.
   facteur_cris =  FACTEUR_HC
   fichbase = '../out/base.hc'
   fichtab  = '../out/jonction.hc'
   facteur_boite = facteur_boite_HC
   angle_vis = angle_vis_HC
elseif (crystal_structure(1:3) == 'ORT') then
   ORT = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_ORT'
   facteur_cris = FACTEUR_ORT
   fichbase = '../out/base.ort'
   fichtab  = '../out/jonction.ort'
   facteur_boite = facteur_boite_ORT
   angle_vis = angle_vis_HC
else
   stop " structure crystallographique non identifiee ,  BY"
endif

if(GLDEV .and. .not. CFC) then
   print *, " GD active pour ",crystal_structure(1:3)
   print *, " le calcul de force de deviation est indetermine car, en general, "
   print *, " le plan de deviation n'est pas unique et la selection n'est pas encore "
   print *, " effectuee dans elasti. Desole, BY."
   stop
endif

if (iterinfo <= 0) iterinfo = 100000000
if(.not. BCC .and. .not. HC .and. .not. CFC) then
   print*, " la constante angle_vis n'est pas definie ?, By !"
   stop
endif
!==========================================================================================
!==========================================================================================
!==========================================================================================
write(*, '(/,"Le Facteur homogeneisation spatiale = ", I5)') facteur_cris

VecBurgers = VecBurgers * 1.D-10    ! mise a l'echelle reelle de vecteur de burgers
if (ORT .or. CS) then
   maille_mat=VecBurgers
elseif (HC .or. CFC) then
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

resolution = echelle * VecBurgers   ! resolution effective de la simulation
vitesse_unitaire = resolution/deltat


! =======================================================================
! ====================  calcul de TauINT_limite       =================
! ===C'est la contrainte interne au dela de laquelle on considere qu'on =
! === est proche d'un reaction local et que la prediction de la vitesse =
! === ne doit plus suivre ===============================================
! =======================================================================

TauINT_LIMITE = XMU / (2.0*PII*Echelle*DisReacloc)
write (*,'(/," Contrainte interne limite pour les reactions locales = ",F6.2," MPa")') &
TauINT_LIMITE * 1.0D-6


!==========================================================================================
!==========================================================================================
!==========================================================================================
! Calcul de l'unite de la simulation pour avoir la taille des vecteurs
! vis de la simulation (e.g.|[022]|*avalue) egale a la resolution physique
! souhaitee, c-a-d egale a echelle*Vecteur Burgers reel du materiau
! Ceci permet de conserver tous les parametre de la simulation inchanges
avalue = echelle * maille_mat / facteur_cris    ! donne Burgers(simu)= echelle*Burgers(vrai)
write (*,'(/," Unite de la simulation (avalue) = ", F8.5," A" )') avalue * 1D10
!==========================================================================================
! initialisation des vitesses maximales

SOUND(1:NLV_max) = QUATRE * vitesse_unitaire/avalue
UNSOUND(1:NLV_max) = UN/SOUND(1:NLV_max)
!==========================================================================================
!==========================================================================================


! Un petit message pour me prevenir du cas limite ou avalue est plus petit que 10b
if (avalue < 10 * maille_mat) then
   write (*,'(/,"Attention : UNITE DE SILULATION = ", F6.3," du seuil elastique")') &
   avalue/(10 * maille_mat)
endif

Bdiva = VecBurgers / avalue
write(*, '(/,"le parametre BdivA = ", F10.5)') bdiva


! =========================================================
! Definition de la longuerur de discretisation homogeneisee
write(*, '(/," Discretisation des systemes actifs (micron)     = ", F6.3)') Ldis_act
write(*, '(/," Discretisation des systemes non actifs (micron) = ", F6.3)') Ldis_nact
Ldis_act = Ldis_act * 1D-6   ! longueur de discretisation en metre (sys actif)
Ldis_nact = Ldis_nact * 1D-6   ! longueur de discretisation en metre (sys non actif)
Xlomax    = DNINT(Ldis_act/avalue) ! longueur de discretisation en (a)
Xlomax_nact  = DNINT(Ldis_nact/avalue) ! longueur de discretisation en (a) (sys non actif)

! Ecriture a  l ecran pour verification des comformite
write(*, '(/,"Longueur de discretisation (en a) = ", F8.1)') Xlomax

write(*, '(/,"Efficacite des crans  = ", F5.0," %" )') Efficacite_crans
Efficacite_crans = Efficacite_crans / 100.0

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

RETURN

end subroutine parametrage

!#########################################################################
!   Sous-programme necessaire pour l'initialisation des fichier en fonction
!   de la clef (SIDEJA) : demarrage/redemarrage de la simulation
!#########################################################################
Subroutine initialisation_fichiers

implicit none

integer :: cristallo

!#######################################################################
!# Lecture des fichiers contenants les segments et de la table des
!# voisins lorsque l'on repart d'une precedante simulation
!#######################################################################
select case (SIDEJA)

!############################
case (0)  ! Demarage standart

! initialisation des fichiers d'ecriture
   open( 8,FILE='../out/stat',STATUS='REPLACE')
   close( 8)
   open(18,FILE='../out/sigeps',STATUS='REPLACE')
   close(18)
   open(28,FILE='../out/rau',STATUS='REPLACE')
   close(28)
   open(38,FILE='../out/gamma',STATUS='REPLACE')
   close(38)
   open(48,FILE='../out/travapp',STATUS='REPLACE')
   close(48)
   open(49,FILE='../out/travint',STATUS='REPLACE')
   close(49)
   open(25,FILE='../out/gammap',STATUS='REPLACE')
   close(25)

   open(99,FILE='../out/bigsave.bin',FORM='UNFORMATTED',STATUS='REPLACE')

   open(95,FILE='../out/histo.txt',STATUS='UNKNOWN')

   open(92,FILE='../out/film',FORM='UNFORMATTED',STATUS='REPLACE')
!   open(94,FILE='../out/SEGO3D',FORM='FORMATTED',STATUS='REPLACE')
   if (CS)  cristallo = 1
   if (BCC) cristallo = 2
   if (CFC) cristallo = 3
   if (HC)  cristallo = 4
   if (ORT) cristallo = 5

   write (92) cristallo
!############################################################
case (1) ! Cas ou on repart d'une configuration deja utilisee


! Ouverture des unites d'ecriture
   open(92,FILE='../out/film',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')

!##################################################
case (2) ! redemarage avec changement de parametres

! Ouverture des unites d'ecriture
   open( 8,FILE='../out/stat',STATUS='REPLACE')
   open(18,FILE='../out/sigeps',STATUS='REPLACE')
   open(28,FILE='../out/rau',STATUS='REPLACE')
   open(38,FILE='../out/gamma',STATUS='REPLACE')
   open(48,FILE='../out/trav',STATUS='REPLACE')
   open(49,FILE='../out/lj',STATUS='REPLACE')
   open(25,FILE='../out/gammap',STATUS='REPLACE')
   open(92,FILE='../out/film',FORM='UNFORMATTED',STATUS='REPLACE')
!   open(94,FILE='../out/SEGO3D',FORM='FORMATTED',STATUS='REPLACE')
   close( 8)
   close(18)
   close(28)
   close(38)
   close(48)
   close(49)
   close(25)
   close(99)
end select

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
allocate (DEVSEG(NBASE,NBASE))
allocate (DEVSEGtemp(NBASE))
allocate (CONEC(NBASE,2))
allocate (DECONEC(NBASE,NBASE))
allocate (SC(nbase,nbase))
allocate (CoefSC(nbase,nbase))
allocate (AXEJONCi(NBASE,NBASE))
allocate (AXEJONCj(NBASE,NBASE))
allocate (SYSOBST(NBASE,NBASE))
allocate (SYSCONNEC(NBASE,NBASE))
allocate (SYSSTAT(NBASE,NBASE))    !*** type d'interaction
allocate (LIEN(NBASE,NBASE))
allocate (LIENS(NBASE,NBASE,2))
allocate (NBROT(NBASE,NBASE,3))
allocate (GDROT(NBASE,NBASE))
allocate (INVCONEC(NBASE,NBASE))
allocate (PartVis(NBASE))
allocate (SYSEG(NBASE))
allocate (LOMA(NBASE))
allocate (JSLOMA(NBASE))
allocate (LOCLMA(NBASE))
allocate (FAC1DEP(NBASE,NBASERED))
allocate (FAC2DEP(NBASE,NBASERED,NBASERED))
allocate (MODDEPLA(NBASE,NBASERED,NBASERED))
allocate (VECNORLIN(3,NBASE))
allocate (VECNORDEP(3,NBASE))
allocate (NORMLIN(NBASE))
allocate (NORMDEP(NBASE))
allocate (tau0(nbase))
allocate (Bf(nbase))
allocate (RAUSYS(NTSG))
allocate (CranSys(NTSG))
allocate (AIRESYS(NTSG))
allocate (AIRESYSInst(NTSG))
allocate (GAMSYS(NTSG))
allocate (AireVisSys(NTSG))
allocate (AireCoinSys(NTSG))
allocate (TrAppSys(NTSG))
allocate (TrAppSysInst(NTSG))
allocate (TrIntSys(NTSG))
allocate (TrIntSysInst(NTSG))
allocate (SchmidSys(NTSG))
allocate (GAMMADOTSYS(NTSG))
allocate (Numero_loi(NBASE))
allocate (Inn(NBASE))


end subroutine Allocation


!************************************************************************
! Ce sous-programme sert a generer la base des vecteurs utilises pour
! la discretisation des dislocations: vecteur ligne, plan et deplacement
!  Ghiath Mohanmed -  15/01/01
!************************************************************************
subroutine generer_basedevecteurs

implicit none

! Declarations
real(kind=DP)  :: fac(nbasered/2),taille,facr,projection,depsurprojecvis,tailledep,preci
real(kind=DP)  :: projecsurdep,depsurprojeccoin,depsurprojec

integer(DPI)  :: matrix (3,3), vecteur(3), plan (3),direction(3),i,j,inverse,correcvis
integer(DPI)  :: itemp,vis_dir(3),coin_dir (3),plan_ref(3), direction_ref (3),ivo,ive,mixte(3)
integer(DPI)  :: normeseg(nbasered/2-1),facteur,facteurdep(nbasered),essai(3)
integer(DPI)  :: correcbase1,transi(3),correcdep1,correcbase2,br
integer(DPI)  :: indicesys, uniformisation,limitej, correcdep2
integer(DPI)  :: facteur1, facteur2,depmixte(3),correcbase,vistemp(nbasered),cointemp(nbasered)
real  (DP)       :: dep


! simplification d'ecriture
br = nbasered
limitej = 100

! facteurvis et coin traduisent le fait que les directions vis et coins initiales
! doivent etre des multiple fixes par la cristallo
if (HC) then
   facteurvis (:)  = (/3.0,3.0,3.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurcoin (:) = (/4.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurmixt (:) = (/1.0,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
elseif (BCC) then
   facteurvis (:)  = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurcoin (:) = (/2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurmixt (:) = (/0.334,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
elseif (ORT) then
   facteurvis (:)  = (/6.0,5.0,6.0,5.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurcoin (:) = (/10.0,10.0,5.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurmixt (:) = (/1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0/)
elseif(CFC) then
   facteurvis (:)  = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurcoin (:) = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurmixt (:) = (/0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
elseif(CS) then
   facteurvis (:)  = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurcoin (:) = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
   facteurmixt (:) = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
endif

! une suite qui indique la facon avec laquelle la base reduite est construite :
! ce sont les composate vis et coin des vecteur de la base reduite
! vis, vis+coin, coin, coin-vis .....
vistemp(:) =  (/1, 1, 0,-1,-1,-1, 0, 1/)
if (BCC) vistemp(:) =  (/1, 2, 0,-1,-1,-2, 0, 1/)
cointemp(:) = (/0, 1, 1, 1, 0,-1,-1,-1/)

! saisie des plan et direction de glissement de refernce (deje lus)

!print *, " "
!write(*,'(" Fichiers des Matrices = ",A50)') fichier


open (12,FILE=fichier,STATUS='OLD')

!*** Construction d'une base generale de nombre de systems de glissement * 8
!*
!*** Ordre choisi :
!*   vis,m1=coin+vis,coin,m2=coin-vis,-vis,m3=-m1,-coin,m4=-m2 et la meme
!*  chose avec le deuxieme vecteur coin associe au mme vecteur vis soit 16
!*  vecteurs par vecteur vis (8 vecteurs pour chacun des 12 systemes de
!*  glissement).
!*** Connectivite :
!*   Un segment peut etre connecte a un autre segment que si ce dernier le
!*  precede ou le suit a modulo 8 pres dans la liste des vecteurs du systeme
!*  de glissement
!*** veclin : vecteur ligne
!*   vecdep : vecteur deplacement associe

do indicesys = 1,nb_slip_types
    inverse = 1

! On verifie la conformite entre le fichier materiau et la base cristallo
    read (12,*) vecteur(1:3)
    if (etat(vecteur,slip(indicesys)%planes(1:3)) /= 3) &
    stop " Plans dans MATERIAU incompatible avec ROTATIONS"

    read (12,*) vecteur(1:3)
    if (etat(vecteur,slip(indicesys)%directions(1:3)) /= 3) &
    stop " Burgers dans MATERIAU incompatible avec ROTATIONS"

    read (12,*) vecteur(1)
    if (vecteur(1) /= slip(indicesys)%nsystemes) &
    stop " nombre de syst. dans MATERIAU incompatible avec ROTATIONS"

! generation des boucles associe a chaque systeme
    do i = 1, slip(indicesys)%nsystemes
        do itemp=1,3
            read (12,*) matrix(itemp,1:3)
        end do

! calcul des vecteurs plan, Burgers, et direction de la coin
        plan_ref(1:3) = Slip(indicesys)%planes(1:3)
        direction_ref(1:3) = Slip(indicesys)%directions(1:3)

        plan = matmul(matrix,plan_ref)
        plan = reduire_vec(plan)    ! le vecteur entier le plus petit

        if (HC .and. modulo(i,IDEUX) == 0) plan (:) = -1 * Plan (:)
        if (HC) inverse = 1

        vis_dir = inverse * matmul(matrix,direction_ref)
        vis_dir = reduire_vec(vis_dir)     ! le vecteur entier le plus petit

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Attension ! Convecntion : si la coin est definie par la formule suivante
!     alors la boucle de la base des vecteurs s'ettend avec une cission
!     dans le meme sens du vecteur VIS.! reper direct : Plan,Coin,Vis,Plan,Coin,Vis ..
! cela signifi egalement que le demi plan associe au Coin1 est // au plan
        coin_dir = reduire_vec(prodivec(vis_dir,plan))

!     Si on prend   coin_dir = - prodivec(vis_dir,plan), la boucle collapse
!     avec la meme contrainte
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (HC .or. ORT) vis_dir(:) = vis_dir(:) * facteurvis(indicesys)
        if (BCC .or. HC .or. ORT) coin_dir(:) = NINT(coin_dir(:) * facteurcoin(indicesys))
! le but etant de recuperer un vecteur Burgers et C entiers
!      On obtent alors : [ 0 3 3] a la place de 1/2[0 1 1]
!       et une translation selon C = [4 4 4] au lieu de 4/3[1 1 1] = coin_dir
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! nouvelle procedure pour optimiser les vecteurs de deplacement dex mixtes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pour la mixte m1 = vis + coin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        preci = norivect(vis_dir) * PRECISION

! signification des parametres:
! correcbase1: le facteur multiplicateur de la base entiere pour conformite avec la mixte1
! correcbase2: le facteur multiplicateur de la base entiere pour conformite avec la mixte2
! correcbase: le facteur multiplicateur de la base entiere pour les deux mixte
! correcdep1: facteur multiplicateur du deplacement de m1 pour comptabilite avec correcbase
! correcdep2: facteur multiplicateur du deplacement de m2 pour comptabilite avec correcbase
        correcbase1 = 1
! calul de du vecteur ligne de m1
        if (.not. BCC) then
           mixte(1:3) = facteurmixt(indicesys)*(coin_dir(1:3) - vis_dir(1:3))  ! vecteur mixte  2
        else
           mixte(1:3) = facteurmixt(indicesys)*(coin_dir(1:3) - 2*vis_dir(1:3))  ! vecteur mixte  2
        endif
! calul de du vecteur deplacement de m1
        depmixte(1:3) = reduire_vec(prodivec(plan(1:3),mixte(1:3))) ! deplacement minimal

        depsurprojecvis = abs(real(inorivect(depmixte(1:3)),DP)/  &
        real(iprodsca(vis_dir(1:3),depmixte(1:3)),DP))

        depsurprojeccoin = abs(real(inorivect(depmixte(1:3)),DP)/  &
        real(iprodsca(coin_dir(1:3),depmixte(1:3)),DP))

! test du cas ou la base n'est pas symetrique : BCC
        if (.not. egalite (depsurprojecvis,depsurprojeccoin)) then
           if(depsurprojecvis < depsurprojeccoin) then
              print *, depsurprojecvis , depsurprojeccoin
              stop" projection vis > projection coin: cas non prevu, desolé. BY"
           else
              correcvis = nint(depsurprojecvis/depsurprojeccoin)
              if (.not. egalite(depsurprojecvis,depsurprojeccoin*correcvis)) then
                 print *, depsurprojecvis , depsurprojeccoin*correcvis
                 stop " projection coin /= entier * projection vis. Trop complique. BY"
              else
                 depsurprojec = depsurprojeccoin
              endif
           endif
        else
           depsurprojec = depsurprojeccoin
        endif

! calul du rapport du deplacement de m1 sur la projection des voisins de m1
!        depsurprojec = abs(real(inorivect(depmixte(1:3)),DP)/real(iprodsca(vis_dir(1:3),depmixte(1:3)),DP))
        projecsurdep = 1.0/depsurprojec

!       print *, " sys =", indicesys
!        print *, " vis_dir(1:3) =", vis_dir(1:3)
!        print *, "  coin_dir(1:3) =",  coin_dir(1:3)
!       print *, " mixte(1:3) =", mixte(1:3)
!       print *, " depmixte(1:3) =", depmixte(1:3)
!        print *, " depsurprojec =", depsurprojec

! 1) le paradis : le deplacement = la projection des voisins : rien a faire
        if (egalite(depsurprojec,UN)) then
           correcdep1 = 1 ; correcbase1 = 1
        elseif (entier(depsurprojec)) then  ! le deplacement est multiple de la projection
           correcdep1 = 1 ; correcbase1 = nint(depsurprojec,DPI)
        elseif (entier(projecsurdep)) then  ! la projection est multiple du deplacement
           correcdep1 = nint(projecsurdep,DPI) ; correcbase1 = 1
! le rapport est un nombre rationnel, on calcul le facteur multiplicatif pour obtenir un entier
        else
           do j = 1,limitej  ! facteur d'homotethie de la base, plafonne a 20
!               print *,j, depsurprojec * j
               if (entier(depsurprojec * j)) then
                  correcdep1 = j ; correcbase1 = nint (depsurprojec * j)
                  exit
               endif
           enddo
           if (j == limitej ) stop " indice j depasse la boucle1"
        endif

        correcbase2 = 1
! calul de du vecteur ligne de m2
        mixte(1:3) = facteurmixt(indicesys)*(coin_dir(1:3) - vis_dir(1:3))  ! vecteur mixte  2
! calul de du vecteur deplacement de m2
        depmixte(1:3) = reduire_vec(prodivec(plan(1:3),mixte(1:3))) ! deplacement minimal
! calul du rapport du deplacement de m2 sur la projection des voisins de m2
        depsurprojecvis = abs(real(inorivect(depmixte(1:3)),DP)/  &
        real(iprodsca(vis_dir(1:3),depmixte(1:3)),DP))

        depsurprojeccoin = abs(real(inorivect(depmixte(1:3)),DP)/  &
        real(iprodsca(coin_dir(1:3),depmixte(1:3)),DP))

! test du cas ou la base n'est pas symetrique : BCC
        if (.not. egalite (depsurprojecvis,depsurprojeccoin)) then
           if(depsurprojecvis < depsurprojeccoin) then
              print *, depsurprojecvis , depsurprojeccoin
              stop" projection vis > projection coin: cas non prevu, desolé. BY"
           else
              correcvis = nint(depsurprojecvis/depsurprojeccoin)
              if (.not. egalite(depsurprojecvis,depsurprojeccoin*correcvis)) then
                 print *, depsurprojecvis , depsurprojeccoin*correcvis
                 stop " projection coin /= entier * projection vis. Trop complique. BY"
              else
                 depsurprojec = depsurprojeccoin
              endif
           endif
        else
           depsurprojec = depsurprojeccoin
        endif
        projecsurdep = 1.0/depsurprojec
!        print *, " sys =", indicesys
!        print *, " mixte2(1:3) =", mixte(1:3)
!        print *, " depmixte2(1:3) =", depmixte(1:3)
!        print *, " depsurprojecvis =", depsurprojecvis
!!        print *, " depsurprojeccoin =", depsurprojeccoin
!        print *, " depsurprojec =", depsurprojec

! 1) le paradis : le deplacement = la projection des voisins
        if (egalite(depsurprojec,UN)) then
           correcdep2 = 1 ; correcbase2 = 1
        elseif (entier(depsurprojec)) then  ! la de placement est multiple de la projection
           correcdep2 = 1 ; correcbase2 = nint(depsurprojec,DPI)
        elseif (entier(projecsurdep)) then  ! la projection est multiple du deplacement
           correcdep2 = nint(projecsurdep,DPI) ; correcbase2 = 1
        else
           do j = 1,limitej  ! facteur d'homotethie de la base, plafonne a 20
!               print *,j, depsurprojec * j
               if (entier(depsurprojec * j)) then
                  correcdep2 = j ; correcbase2 = nint (depsurprojec * j)
                  exit
               endif
           enddo
           if (j == limitej ) stop " indice j depasse la boucle1"
        endif
!        print *, " correcdep1", correcdep1;       print *, "correcbase1 ", correcbase1
!        print *, "correcdep2 ", correcdep2;       print *, "correcbase2 ", correcbase2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Traitement des conflits entre mixte1 et mixte2
! calcul du facteur de correction totale de la base et
! remise a jour du facteur de multiplication des vecteur deplacement

        if (Modulo(correcbase1,correcbase2) == 0) then
           correcbase = correcbase1
           correcdep2 = correcdep2 * correcbase1/correcbase2
        else if (Modulo(correcbase2,correcbase1) == 0) then
           correcbase = correcbase2
           correcdep1 = correcdep1 * correcbase2/correcbase1
        else
           print *, " Cas merdique : les correction des vis pour m1 et m2 ne sont pas identiques"
           print *, "correcbase1 =",correcbase1
           print *, "correcbase2 =", correcbase2
           print *, " attention il faut choisir normalement le plus petit diviseur commun descorrecbase"
           correcbase = correcbase1 * correcbase2
           correcdep1 = correcdep1 * correcbase2
           correcdep2 = correcdep2 * correcbase1
        endif


! homogeneisation des longueurs entre les differents type de systems de glissement
        if (HC .and. indicesys == 1) then
           correcbase = correcbase * 2
           correcdep1 = correcdep1 * 2
           correcdep2 = correcdep2 * 2
        endif

        if (HC .and. indicesys > 1) then
           correcbase = correcbase * 1
           correcdep1 = correcdep1 * 1
           correcdep2 = correcdep2 * 1
        endif

        if (ORT .and. indicesys == 1) then
           correcbase = correcbase * 61
           correcdep1 = correcdep1 * 61
           correcdep2 = correcdep2 * 61
        endif

        if (ORT .and. indicesys == 2) then
           correcbase = correcbase * 1037
           correcdep1 = correcdep1 * 1037
           correcdep2 = correcdep2 * 1037
        endif

        if (ORT .and. indicesys > 2) then
           correcbase = correcbase * 17
           correcdep1 = correcdep1 * 17
           correcdep2 = correcdep2 * 17
        endif

!        print*, "correcbase =", correcbase; print*, "correcdep1 =", correcdep1
!        print*, "correcdep2 =", correcdep2

        vis_dir (1:3)  = correcbase * vis_dir(1:3)
        coin_dir (1:3) = correcbase * coin_dir(1:3)
        facteur = 0
! generation des vecteurs de la base
        do j = 1, br
            itemp = ((indicesys - 1)*slip(1)%nsystemes + i-1) * br + j
! generation des vecteurs plans (bvecnor)
            bvecnor(1:3,itemp) = plan (1:3)

! generation des vecteurs ligne des segments et de leurs normes

            bveclin(1:3,itemp) = vistemp(j)*vis_dir(1:3) + cointemp(j)*coin_dir(1:3)
            if (modulo(j,IDEUX)== 0) bveclin(1:3,itemp)=bveclin(1:3,itemp)*facteurmixt(indicesys)
!            print*, "itemp =", itemp, "  bveclin =", bveclin(1:3,itemp)
            VecNorLin(1:3,itemp)=normaivec(BVECLIN(1:3,itemp))

! generation des vecteurs de deplacement des segments et de leurs normes

            vecteur(1:3) = reduire_vec(prodivec(bvecnor(1:3,itemp),bveclin(1:3,itemp)))
! Attention les vecteur deplacement ne sont pas encore corriges
            if (modulo(j,IQUATRE) == 0 ) then   ! cas des mixte2
               facteur = correcdep2
            elseif(modulo(j,IDEUX) == 0) then ! cas des mixtes 1
               facteur = correcdep1
            elseif (modulo(j,IQUATRE) == 1) then ! cas des vis
               facteur = facteurcoin(indicesys)*correcbase
            else
               facteur = facteurvis(indicesys)*correcbase
               if (BCC) facteur = facteur / 2
            endif
            bvecdep(1:3,itemp) = vecteur(1:3) * facteur
            VecNorDep(1:3,itemp)=normaivec(BVECDEP(1:3,Itemp))
! write (*,61) itemp, bvecnor(1:3,itemp), bveclin(1:3,itemp),bvecdep(1:3,itemp)
        enddo
        if (.not.BCC)inverse = -1  * inverse
    enddo
enddo

close(12)
!call write_basedevecteurs
!reste a faire une petite boucle de multiplication des vecteur deplacement
!En principe les vecteur de la base ont etet construits de tel facon que les vecteurs de
!deplacement sont ou bien egaux a la projection des voisin ou inferieurs et dans ce
!cas faut les rallonger pour etre egaux a la projection des voisins.
!Le rallongement est forcement entier , sinon c'est une erreuurrrrrrrrrrr

!###########################################################
! Mise a jours des norme de vecteurs de deplacement et
! verification de la compatibilite de la base des vecteurs
do j = 0, NTSG - 1
    do itemp = 1, br    ! 1-----4
        ivo = itemp - 1   ! voisin en o
        ive = itemp + 1   ! voisin en e
        if (itemp == 1) ivo = br ! correction
        if (itemp == br) ive = 1  ! correction
        i = j*br +itemp   ! vrais indice du segment 1---nbase
        if (HC .and. i > 24 .and. modulo(i,IDEUX) ==1) bvecdep(:,i) = bvecdep(:,i)/2
        ivo = ivo + j*br   ! vrais indice du voisin en o
        ive = ive + j*br   ! vrais indice du voisin en e

61      format(I5,' :  (',3(I3,1X),') [',3(I5,1X),']. Dep. : [',3(I5,1X),']')
        facteur1 = iprodsca(bveclin(:,ivo),bvecdep(:,i))
        facteur2 = iprodsca(bveclin(:,ive),bvecdep(:,i))


! 1) choisir le voisin dont la projection sur bvecdep est la plus petite
        if (abs(facteur1) < abs(facteur2)) then
           projection = abs(real(facteur1,DP))/norivect(bvecdep(:,i))
        else
           projection = abs(real(facteur2,DP))/norivect(bvecdep(:,i))
        endif

! 2) le deplacement de i ne doit depasser la projection de voisin
        dep = norivect(bvecdep(1:3,i))
        preci = dep * PRECISION

! 3) si la projection n'est pas egale au deplacment de I, on change bvecdep
        if (.not.egalite(projection,dep)) bvecdep(:,I) = nint(Bvecdep(:,I)*(projection/dep),DPI)

! 4) Verification definitive de projection = bvecdep
        projection = abs(real(iprodsca(bveclin(:,ivo),bvecdep(:,i)),DP)/norivect(bvecdep(:,i)))
        dep = norivect(bvecdep(1:3,i))
        preci = dep * PRECISION
        if (.not.egalite(real(nint(projection/dep),DP),projection/dep)) then
           write (*,61) ivo, bvecnor(1:3,ivo), bveclin(1:3,ivo),bvecdep(1:3,ivo)
           write (*,61) i, bvecnor(1:3,i), bveclin(1:3,i),bvecdep(1:3,i)
           write (*,61) ive, bvecnor(1:3,ive), bveclin(1:3,ive),bvecdep(1:3,ive)
           print *, " "
           print *, "i= ", i," projec =",projection," dep =", dep
           print *, "dep = ", dep," projec/dep =",nint(projection/dep,DPI)
           stop " Apres correction : projection de voisin /= bvecdep"
        endif


! 4) Si le plan de glissement n'est pas perpendiculaire au l ou dep stop
        if ( iprodsca (bvecnor(:,i),bvecdep(:,i)) /= 0) stop " n non perpendicualir sur dep"
        if ( iprodsca (bvecnor(:,i),bveclin(:,i)) /= 0) stop " n non perpendicualir sur l"
    enddo
enddo

end subroutine generer_basedevecteurs




!#############################################################################
!############           initialisation des variables globales       ##########
!#############################################################################
subroutine initialiser

implicit none


!**** PARANNOIAC SETTING ****
RAUDIS = ZERO
RAUDMO_I = ZERO
RAUDJONC_I = zero
VITMOY = ZERO
EPSDOT = ZERO
! Initialisation du tenseur de deformation plastique
EPSPLA(1:3,1:3) = ZERO
kk = 0
npar = izero
accutime = zero
airevis = zero
airecoin = zero
airemixte = zero
AireSys(1:NTSG) = zero
AireSysInst(1:NTSG) = zero
AireVisSys(1:NTSG)  = zero
AireCoinSys(1:NTSG) = zero
TrAppSys(1:NTSG) = zero
TrAppSysInst(1:NTSG) = zero
TrIntSys(1:NTSG) = zero
TrIntSysInst(1:NTSG) = zero
TrApp = zero
TrAppInst = zero
TrInt = zero
TrIntInst = zero
centrei(1:3) = -1
IDEP(1:NSEGMAX)=0
IBOITE(1:NSEGMAX)=0
IXCOORB(1:NSEGMAX)=0
IYCOORB(1:NSEGMAX)=0
IZCOORB(1:NSEGMAX)=0
QUIDEP(1:NSEGMAX)=0
RO(3,NSEGMAX)=0.
ForceTOT(NSEGMAX) = zero
ForceTL(NSEGMAX)=zero
ForceAP(NSEGMAX)=zero
DAIRE(NSEGMAX)=zero
SEG(1:NSEGMAX)%VECLIN   = nbase
SEG(1:NSEGMAX)%NORME	= 0
SEG(1:NSEGMAX)%VOISO	= nsegmax
SEG(1:NSEGMAX)%VOISE	= nsegmax
SEG(1:NSEGMAX)%VNNO	= nsegmax
SEG(1:NSEGMAX)%VNNE	= nsegmax
SEG(1:NSEGMAX)%IVNNO	= 1
SEG(1:NSEGMAX)%IVNNE	= 1
SEG(1:NSEGMAX)%IJONC	= nsegmax
SEG(1:NSEGMAX)%RESDEP   = zero
SEG(1:NSEGMAX)%JONC	= .false.
SEG(1:NSEGMAX)%GD	= izero
SEG(1:NSEGMAX)%DISEG	= .false.
SEG(NSEGMAX)%norme = nsegmax
Loi(1:NLV_max)%Arrhenius  = 0
Loi(1:NLV_max)%bf = zero
Loi(1:NLV_max)%tau0  = zero
Loi(1:NLV_max)%H0    = zero
Loi(1:NLV_max)%deltaG0eff       = zero
Loi(1:NLV_max)%tau0eff       = zero
Loi(1:NLV_max)%tauath     = zero
Loi(1:NLV_max)%coef_p        =zero
Loi(1:NLV_max)%coef_q        = zero
Solli_sys(1:NTSG_max) = UN
mufoisadivb = zero
RAUSYS(:)= zero
CRANSYS(:)= zero
GAMSYS(:)=zero
GAMMADOTSYS(:)=0.d0
kk_debut_jonc = 0  ! iteration a la quelle commence la premiere jonction
kkjonc = 0         ! nb d'iteration ou la jonction existe
Sigma = zero
RETURN

end subroutine initialiser


!#########################################################################
!# Procedure de tabulation des deplacements et autres tableau            #
!############################################################## 02/2001 ##
subroutine Tabulation

implicit none

integer(kind=DPI)  :: v1(3),v2(3),v3(3),d1(3),d2(3),d3(3),Jindice
integer(kind=DPI)  :: i,j,jj,ii,k,qq,nbcouple,bi(3), bj(3),li(3),lj(3),bjonc(3)
integer(kind=DPI)  :: buf,buf1,buf2,vec(1:3),indj,indk,lignei(3), lignej(3)
integer(kind=DPI)  :: bufaxe,pbufaxe,Libufaxe,Ljbufaxe,k1,k2,br,ni(3),nj(3)
real(DP)           :: depsurprojec, projecsurdep,x,rv(3),rv1(3)
logical :: CtiUneVis,condition

br = nbasered      ! pour faciliter l'ecriture

!*** TABULATION ANTI SPIRALE EN PT D'ANCRAGE :
!     RENVOIE LE CARACTERE REMPLACANT LES DEUX INITIAUX
!      (Juste la traduction en indice de la base réduite d'une somme vectorielle)
SC(:,:)=0
CoefSC(:,:)=0
!allocate (SC(nbase,nbase))
!allocate (CoefSC(nbase,nbase))
qq = 0
do i = 1,nb_slip_types
    do j = 1, slip(i)%nsystemes
        k = qq*nbasered
!        print *, i,j,k
        if(facteurmixt(i) == 1 .or. facteurmixt(i) == 0.334) then
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
           SC(k + 7,k + 4)= k + 5	! ATTENTION : il faut une base de type PM VCM1M2
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
        else
           stop " facteur mixt non pris en compte pour la tabulation de SC"
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
!              print *, "i,j ", i ,j
!              print *, "Sc, Coef ", sc(i,j),coefsc(i,j)
!              print *, v1
!              print *, v2
!              stop
           endif
        endif
    enddo
enddo

!*** Construction du tableau des vecteurs normalises
do I=1,NBASE

    if (bveclin(1,i) /= 0) then
       Inn(i) = 1
    elseif (bveclin(2,i) /= 0) then
       Inn(i) = 2
    elseif (bveclin(3,i) /= 0) then
       Inn(i) = 3
    else
       stop " Base de vecteur erronee"
    endif

    NormLin(I) = norivect(BVECLIN(1:3,I))
    VecNorLin(1:3,I) = normaivec(BVECLIN(1:3,I))
    NormDep(I) = norivect(BVECDEP(1:3,I))
    VecNorDep(1:3,I) = normaivec(BVECDEP(1:3,I))
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
    SISEG(I)   = sicose2v(bveclin(1:3,I),bveclin(1:3,assoc(I,1)))
    if (siseg(I).eq.0) SISEG(I)=siseg(I-1)
!    BISYSEG(I) = ((I-IUN)/(IDEUX*NBASERED))+IUN    ! bisysteme de glissement
    SYSEG(I)   = ((I-IUN)/(NBASERED))+IUN      !   systeme de glissement
    rv(:) = dble(bvecnor(:,I))
    x = cose2v(Z,rv)
    rv(:) = dble(bveclin(:,assoc(I,1)))
    SchmidSys(SYSEG(I)) = x*cose2v(Z,rv)
enddo

do I=1,NBASE
    PartVis(i) = Real(abs(iprodsca(bveclin(:,i), Bveclin(:,assoc(I,1)))),DP)
    PartVis(i) = PartVis(I)/normlin(i)/normlin(assoc(I,1))
!    print *, i, partvis (i)
    do J=1,2,1
!*** indice j
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
        if (J.eq.conec(I,1).or.J.eq.conec(I,2)) then
           NBROT(I,J,1)=0 !*** Segment directement connectable
           NBROT(I,J,2)=0
           NBROT(I,J,3)=0
        else
           if (I.eq.J) then
              NBROT(I,J,1)=1 !*** 1 rotule
              NBROT(I,J,2)=conec(I,1)
              NBROT(I,J,3)=conec(I,2)
           else
              if (conec(I,2).eq.conec(J,1)) then
                 NBROT(I,J,1)=1 !*** 1 rotule
                 NBROT(I,J,2)=conec(I,2)
                 NBROT(I,J,3)=conec(I,2)
              else
                 if (conec(I,1).eq.conec(J,2)) then
                    NBROT(I,J,1)=1 !*** 1 rotule
                    NBROT(I,J,2)=conec(I,1)
                    NBROT(I,J,3)=conec(I,1)
                 else
                    if (conec(conec(I,1),1).eq.conec(j,2).and.(nbasered.ge.8)) then
                       NBROT(I,J,1)=2 !*** 2 rotule
                       NBROT(I,J,2)=conec(I,1)
                       NBROT(I,J,3)=conec(J,2)
                    else
                       if (conec(conec(I,2),2).eq.conec(j,1).and.(nbasered.ge.8)) then
                          NBROT(I,J,1)=2 !*** 2 rotule
                          NBROT(I,J,2)=conec(I,2)
                          NBROT(I,J,3)=conec(J,1)
                       else
                          NBROT(I,J,1)=(nbasered/2)-1 !*** 3 rotules : recouvrement !!!
                          NBROT(I,J,2)=0
                          NBROT(I,J,3)=0
                       endif
                    endif
                 endif
              endif
           endif
        endif
    enddo
enddo

!**************** TABULATION POUR LE GLISSEMENT DEVIE ************

do I=1,NBASE
    do J=1,NBASE
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
                    if(sicose2v((bveclin(1:3,i)+bveclin(1:3,j)),&
                    bveclin(1:3,assoc(i,1))).gt.0) then
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
MODDEPLA(1:NBASE,1:NBASERED,1:NBASERED) = 0

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
        buf = iprodsca(v1,d2)
        buf1 = iabs(buf)

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
              print*,  " depsurprojec =", depsurprojec
              stop
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
               if(.not. egalite(prodsca(rv,rv1),ZERO)) then
                  print *, " j,i,k ", j,i,k
                  print *, "v1 =", v1,fac1dep(i,j)
                  print *, "v2 =", v2
                  print *, "v3 =", v3,fac1dep(i,k)
                  print *, " I apres dep=",rv
                  print *, "prod =", prodsca(rv,rv1)
                  print *, " fac2dep =", x
                  stop " fac1dep de vo et ve incompatible"
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
                     Stop "INIT:ERREUR 1: il faut changer la procedure de MODDEPLA"
                  endif
               endif

               x = dabs(fac1dep(i,k))    ! compatibilite avec le VE
               if (.not. entier(x)) then
                  if (x < 1) then
                     qq = nint(ppem(x)/x)
                  else
                     Stop "INIT:ERREUR 2: il faut changer la procedure de MODDEPLA"
                  endif
               endif

               x = dabs(fac2dep(i,j,k))    ! compatibilite avec I
               if (.not. entier(x)) then
                  k2 =  nint(ppem(x)/x)
               endif
               moddepla(i,j,k) = ppcm3(jj,int(qq,DPI),k2)
            endif

!            write (*,*) j,i,k," moddepla ",moddepla(i,j,k)
!           write (*,*)dabs(fac1dep(i,j)) ,dabs(fac1dep(i,k)),dabs(fac2dep(i,j,k))
        enddo
    enddo
enddo


!###########################################
!##### tabulation des axes de jonction #####
!###########################################

Devseg(1:NBASE,1:NBASE) = Izero
Devsegtemp(1:NBASE) = Izero
do i=1,nbase,1
    do j=1,nbase,1
        axejonci(i,j)=0         !   INITIALISATION
        axejoncj(i,j)=0         !   INITIALISATION
        SYSOBST(i,j)=999        !   INITIALISATION
        SYSCONNEC(i,j)=999   !POUR CONNECTER LES SEGMENTS ENTRE EUX
        SYSSTAT(i,j)=0   	!0 ANNIHIL ET COPLA 1 HIRTH 2 GLISSILE 3 LOMER
! SYSSTAT : permet de savoir quel type de jonctions on a pour un couple i j

        bi =  bveclin(:,assoc(i,1))       ! vecteur de burgers pour systeme i
        bj =  bveclin(:,assoc(j,1))       ! vecteur de burgers pour systeme j
        bjonc = bi(:) + bj(:)            ! vecteur de jonction brut (resultant du choix initial
! des vecteur de Burgers dans la base reduite
! Attention : on verifie plus tard le signe des vecteur lignes de la jonction pour i et pour j
! Cette modification est introduite pour faire un choix physique des vecteurs de jonction qui doit
! conduire a une jonction mixte de Hirth : (bi + bj) n'est pas perpendiculaire a l'axe de jonction
!
        li = bveclin(:,i)       ! vecteur de ligne i
        lj = bveclin(:,j)       ! vecteur de ligne j
        ni = bvecnor(:,i)       ! normal au plan i
        nj = bvecnor(:,j)

!*** TABULATION DES TYPES D'OBSTACLES ET DE CONNECTIVITES

        if (iabs(etat(ni,nj)) >= 2) then   ! meme plan de gliss
           if (iabs(etat(bi,bj)) >= 2) then     !meme systeme de gliss
              sysconnec(i,j) = 1
              sysOBST(i,j) = 1
           else
              sysconnec(i,j) = 0     ! systems differents mais coplanaires
              sysOBST(i,j) = 0
           endif
        else             ! plan differents
           if (iabs(etat(bi,bj)) >= 2) then       ! meme b
              sysconnec(i,j) = 2     ! vecteur de BUrgers : configuration de GD
! nouveau : 10/10 : devseg(i,J) renvoie l'indice du VL correspondant a I dans
! le systeme devie contenant J
              devseg(i,J) = assoc(J,1) + deassoc(i) -1
! on maintient la forme initial de DEVSEG pour pouvoir continuer a fonctioner pour les CFC
              devsegtemp(i) = devseg(i,J)

!              print *, " i j dev ", I ,J, devseg(i,j)

              if ((iabs(etat(li,bi)) >= 2) .and. (iabs(etat(lj,bj)) >= 2)) then
                 sysOBST(i,j) = 1   ! les segments i ET J correspondent a des vis
              else
                 sysOBST(i,j) = 6
              endif
           else
              sysconnec(i,j) = 3     ! Burgers differents et non coplanaires
              sysOBST(i,j) = 6
           endif
        endif

!*** TABULATION DES COUPLES DE LIGNE DES BINOMES DE JONCTIONS

! 1) d'abord : on ne travaille qu'avec des plans differents
        if (abs(etat(ni,nj)) > 1) cycle      ! ni,nj : non parallels

! 2) et il faut que les vecteurs de burgers soient differents
! la ligne suivante a ete enleve car on trait le cas des systemes devies comme jonction
!        if (abs(etat(bi,bj)) > 1) cycle      ! bi,bj : non parallels

! Cette  condition est etablie pour toutes les crystallo.

!***********************************************************************************
!***********************************************************************************
!***********************************************************************************

! ANNULEE ! critere energitique de Frank: jonction sss bi2 + bj2 >= somme2
!           if (inorivect(somme) <= (inoivect(bi) + inorivect(bj))) then

! identification des couples de jonction pour les deux systemes :
! il existent deux couples, chaque un est constitues de 2 vecteurs
! parallels et de sens oppose (forcement parallels a l'axe
! commun des 2 plans de glissement)

        do ii = assoc(i,1),assoc(i,br)
            do jj = assoc(j,1), assoc(j,br)

! eliminer les segments dans les deux boucles (de i et de j) qui ne sont
! pas parallels a l'axe de jonction .
                if (iabs(etat(bveclin(:,ii),bveclin(:,jj))) < 2) cycle
! Reste les 4 couples possibles

! cas ou bi et bj ne sont pas perpendiculaire
                if (etat(bi,bj) /= 0 ) then

! on introduit le critere de selection des couples de jonction en fonction de l'orientation
! initiale des bi et de bj
! critere etablie apres discussion ladislas-Benoit-Ronan_Ghiath 26/04/01
                   if (etat(bi,bj)*etat(bveclin(:,ii),bveclin(:,jj)) < 0) then

!   (ii - assoc(i,1))/4+1 = 1 pour ii : vis,m1,c,m2 ; et 2 pour -vis, -m1,-c et -m2
!                      CAJ(i,j,(ii-assoc(i,1))/4+1,1) = ii
!                      CAJ(i,j,(ii-assoc(i,1))/4+1,2) = jj
                      if(iprodsca(bi,bvecdep(1:3,ii)) * iprodsca(bj,bvecdep(1:3,jj)) == 0)then
                         SYSSTAT(i,j)=2			! Glissile
                      else
                         IF (abs(etat(bi,bj)) >= 2) THEN
                            SYSSTAT(i,j)=0			! Annihilation
                         ELSE
                            SYSSTAT(i,j)=3			! Lomer
                         ENDIF
                      endif
! inclure une restriction geometrique pour la formation de la
! jonction : les lignes de dislocation ne doivent pas faire
! un angle optu avec la direction de jonction .....
                      if ((etat (li,bveclin(:,ii)) >= 0) .AND. (etat(lj,bveclin(:,jj)) >= 0)) then
                         axejonci(i,j) = ii
                         axejoncj(i,j) = jj
                      endif
                   endif

! cas ou bi et bj sont perpendiculaire : jonction de Hirth
! chercher la configuration de la jonction de Hirth Mixte :
! b_jonction ( = bi + bj) ne doit pas etre perpendiculaire a l'axe de jonction

                else   ! cas ou bi et bj sont perpendiculaires

! 1) Premier Cas : si b_jonction est perpendiculaire a l'axe de jonction, alors
! il faut choisir les couples de jonction de sens oppose pour avoir bjon effective = bi-bj
                   if (etat(bjonc,bveclin(:,ii)) == 0) then
                      if (etat(bveclin(:,ii), bveclin(:,jj)) <= -2) then
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,1) = ii
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,2) = jj
                         if ((etat (li,bveclin(:,ii)) >= 0).AND.(etat(lj,bveclin(:,jj)) >= 0)) then
                            axejonci(i,j) = ii
                            axejoncj(i,j) = jj
                         endif
                      endif

! 2) deuxieme cas : b_jonction n'est pas perpendiculaire a l'axe de jonction, alors
! il faut choisir les couples de jonction de mem sens pour avoir bjonc effective = bi+bj
                   else
                      if (etat(bveclin(:,ii), bveclin(:,jj)) >= 2) then
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,1) = ii
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,2) = jj
                         if ((etat (li,bveclin(:,ii)) >= 0).AND.(etat(lj,bveclin(:,jj)) >= 0)) then
                            axejonci(i,j) = ii
                            axejoncj(i,j) = jj
                         endif
                      endif
                   endif
                endif
            enddo
        enddo

! dans le cas ou les axesjonc sont toujours nuls on assouplie le critere energitique
! et on refait la meme boucle pour choisir une des deux couples qui restent non favorabvle
! energitiquement

        if (axejonci(i,j)*axejoncj(i,j) /= 0) cycle
        do ii = assoc(i,1),assoc(i,br)
            do jj = assoc(j,1), assoc(j,br)

! eliminer les segments dans les deux boucles (de i et de j) qui ne sont
! pas parallels a l'axe de jonction .
                if (iabs(etat(bveclin(:,ii),bveclin(:,jj))) < 2) cycle
! Reste les 4 couples possibles

! cas ou bi et bj ne sont pas perpendiculaire
                if (etat(bi,bj) /= 0 ) then

! on elimine le critere etablie apres discussion ladislas-Benoit-Ronan_Ghiath 26/04/01
!                   if (etat(bi,bj)*etat(bveclin(:,ii),bveclin(:,jj)) < 0) then

!   (ii - assoc(i,1))/4+1 = 1 pour ii : vis,m1,c,m2 ; et 2 pour -vis, -m1,-c et -m2
!                      CAJ(i,j,(ii-assoc(i,1))/4+1,1) = ii
!                      CAJ(i,j,(ii-assoc(i,1))/4+1,2) = jj

! inclure une restriction geometrique pour la formation de la
! jonction : les lignes de dislocation ne doivent pas faire
! un angle optu avec la direction de jonction .....
                   if ((etat (li,bveclin(:,ii)) >= 0) .AND. (etat(lj,bveclin(:,jj)) >= 0)) then
                      axejonci(i,j) = ii
                      axejoncj(i,j) = jj
                   endif
!                   endif

! cas ou bi et bj sont perpendiculaire : jonction de Hirth
! chercher la configuration de la jonction de Hirth Mixte :
! b_jonction ( = bi + bj) ne doit pas etre perpendiculaire a l'axe de jonction

                else   ! cas ou bi et bj sont perpendiculaires

! 1) Premier Cas : b_jonction est perpendiculaire a l'axe de jonction, alors
! il faut choisir les couples de jonction de sens oppose pour avoir bjon effective = bi-bj
                   if (etat(bjonc,bveclin(:,ii)) == 0) then
                      if (etat(bveclin(:,ii), bveclin(:,jj)) <= -2) then
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,1) = ii
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,2) = jj
                         if ((etat (li,bveclin(:,ii)) >= 0).AND.(etat(lj,bveclin(:,jj)) >= 0)) then
                            axejonci(i,j) = ii
                            axejoncj(i,j) = jj
                         endif
                      endif

! 2) deuxieme cas : b_jonction n'est pas perpendiculaire a l'axe de jonction, alors
! il faut choisir les couples de jonction de mem sens pour avoir bjonc effective = bi+bj
                   else
                      if (etat(bveclin(:,ii), bveclin(:,jj)) >= 2) then
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,1) = ii
!                         CAJ(i,j,(ii-assoc(i,1))/4+1,2) = jj
                         if ((etat (li,bveclin(:,ii)) >= 0).AND.(etat(lj,bveclin(:,jj)) >= 0)) then
                            axejonci(i,j) = ii
                            axejoncj(i,j) = jj
                         endif
                      endif
                   endif
                endif
            enddo
        enddo
    enddo
enddo

! extension du tableau nbrot(...) au systemes devies
do i=1,nbase,1
    do j=1,nbase,1
        if(devseg(i,j) /= 0) then
           nbrot(i,j,1) = nbrot(devseg(i,j),j,1)
        endif
    enddo
enddo


! ICI on dit definit le numero de la loi de vitesse associee a chaque vecteur de la base

ii=0
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
!  Ceci est utilise seulement pour les HC
! =======================================================================
!  elles servent a detecter les probleme des sous reseau. le reper est constitue de :
if (.not. CFC) then
! ici on defini les trois translations elementaires pour la construction de
! l'espace, c-a-d le reseau principal.
! ces dierctions sont :
!         a1 : la direction a1 : 1/2 [ 0  1 -1 ] = [0 66 -66]
!         a2 : la direction a2 : 1/2 [-1  1  0 ] = [-66 66 0]
!         c : la direction normal : 2/3 [1 1 1]  = [88 88 88]
! les points de l'espace doivent etre accessibles apartir du point [000] et de ces
! translations uniquement.

! pour les HC, les trois direction entieres doivent etre :
! a1 = [0 66 -66]
! a2 = [-66 0 66]
! c = [-88 -88 -88]

   if (HC) then
      x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
      y_reseau(:) = bveclin (:,nbasered+1)   !========  vecteur vis 2eme systeme
      z_reseau(:) = bveclin (:,3)            !========  vecteur coin 1er systeme
   elseif (ORT) then
      x_reseau(:) = bveclin (:,nbasered+1)            !========  vecteur vis 2eme systeme
      y_reseau(:) = bveclin (:,3)   !========  vecteur coin 1er systeme
      z_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
   elseif (BCC) then
!  Attention: si les vecteurs x,y,z ne sont pas egaux
! aux bveclin 1,4,25 la procedure de decomposition en vecteurs elementaires
! est fausse
      x_reseau(:) = bveclin (:,1)            !========  vecteur vis 1er systeme
      y_reseau(:) = bveclin (:,4)     !========  vecteur mixte 1er systeme
      z_reseau(:) = bveclin (:,25)            !========  vecteur coin 1er systeme
   endif


! =======================================================================
!  rappel des unites de translation pour la definition du reseau spatial
! ====    elles sont definies dans le modules varglob    ==============
! =======================================================================

   write(*,'(/, "Les translation elementaires sont :")')
   write(*, '("      a1  = [ ",3(I6,X)," ]",/)') x_reseau(1:3)
   write(*, '("      a2  = [ ",3(I6,X)," ]",/)') y_reseau(1:3)
   write(*, '("      C   = [ ",3(I6,X)," ]",/)') z_reseau(1:3)

endif
write(*,'(/,"Longueur segemnts elementaire VIS   = ",F8.4," A")')  normlin(1)*avalue*1D10
write(*,'(/,"Longueur segemnts elementaire MIXT  = ",F8.4," A")')  normlin(2)*avalue*1D10
write(*,'(/,"Longueur segemnts elementaire COIN  = ",F8.4," A")')  normlin(3)*avalue*1D10
!#########################################################################################

! Trois vecteur elementaires pour la discretisation de l'espace

end subroutine Tabulation

!************************************************************************
! Ce sous-programme ecrit dans le fichier (BVD.cfc, cc,cs, etc.)
! les vecteurs de la base de simulation
!  Ghiath Monnet -  07/02/01
!************************************************************************
subroutine write_basedevecteurs

implicit none

integer (kind=DPI) :: i,j,itemp,k
character(len=60)  :: list

61 format(I3,': (',3(I5,1X),')[',3(I6,1X),']. Dep: [',3(I6,1X),']')
62 format(11I5)
63 format(I3,"  ",3(1X,I6),"  ",3(I6,1X),"      ",3(I3,1X))



list = '../out/BVD.'//crystal_structure(1:3)

open (14,FILE=list,STATUS='UNKNOWN')  ! fichier de sortie

! d'abord on ecrit la liste brut des vecteur de discretisation
! pour l'utilisation d'autre programmes


write (14,*) nbase , avalue         ! unique de programme : lu par camera ....
do i = 1,nbase
    write (14,63) i, bveclin(1:3,i),bvecdep(1:3,i),bvecnor(1:3,i)
enddo

write(14,*) "  "
write(14,*) "  "
if (HC) write(14,*) "Indice        Plan        ligne (/22)                 deplacement "
if (.not. HC) write(14,*) "Indice        plan       ligne                 deplacement "

! ensuite il y a l'ecriture formatee de la BVD pour une lecture agreable
do itemp = 1,NTSG
    do j = 1,8
        i = (itemp - 1) * nbasered + j
        if (HC) write (14,61) i, bvecnor(1:3,i), bveclin(1:3,i)/22,bvecdep(1:3,i)
        if (.not. HC) write (14,61) i, bvecnor(1:3,i), bveclin(1:3,i),bvecdep(1:3,i)
    enddo
    write (14,*) '    -----------------------------------------------------'
    write (14,*) ' '
enddo

close(14)

end subroutine write_basedevecteurs

!************************************************************************
! Ce sous-programme ecrit dans le fichier (facdep.cfc,bcc,cs, etc.)
! les facteur de deplacement
!  Ghiath Monnet -  07/02/03
!************************************************************************
subroutine write_facdep

implicit none

integer (kind=DPI) :: i,j,itemp,k
character(len=60)  :: list

64 format(" j,i,k =", 3I4," ; fac1(i,j)=",F6.2," ; fac1dep(i,k)=",F6.2," ; fac2dep=",F6.2," ; modep=",I3)



list = '../out/facdep.'//crystal_structure(1:3)

open (14,FILE=list,STATUS='UNKNOWN')  ! fichier de sortie

! d'abord on ecrit la liste brut des vecteur de discretisation
! pour l'utilisation d'autre programmes



! on ecrit pour information la liste des facteur de deplacement
do i = 1, nbase
    do j = 1,nbasered
        do k = 1, nbasered
            write (14,64)j,i,k, fac1dep (i,j), fac1dep(i,k), fac2dep(i,j,k), moddepla(i,j,k)
        enddo
    enddo
enddo


close(14)

end subroutine write_facdep

!##################################################################################
! lecture des sources de dislocation initiale mais il peut servire
! pour generer deux segments en configuration de jonction si nsegm = 0 dans lefichier
! segments. Les informations correspondantes sont le fichier JONCTION
!##################################################################################
subroutine lire_segments

implicit none

real (DP):: taille,taux1,taux2,boite,decal1, decal2,solli1,solli2
integer (kind=DPI) :: i,carac1,carac2,centre(3),sens, Oi(3),Oj(3),ivis1,ivis2,l(3) ,transla(3)
integer(kind=DPI) :: npas1,npas2


fichseg="../in/"//segments

open(3,file=fichseg,STATUS='OLD')

! lecture de la table de sollicitation des differents systemes
! le principe : tauapp(definitive) = tauapp(calcule) * solli_sys(sys)
! solli_sys(i) = 0 on impose tau app = 0  pour le systeme I
! solli_sys(i) = 1 on n'intervient pas sur tau app
read(3,*) Solli_Sys(1:NTSG)

! Nombre de segments dans la configuration initiale
read(3,*) Nsegm
print *, Nsegm
! Dimention des cotes de la boite de simulation
read(3,*) modur(1),modur(2),modur(3)
print *, modur(1),modur(2),modur(3)

! Finalement les segemnts
do I = 1, NSEGM
    read (3,*)	seg(I)%O(1:3),seg(I)%norme,SEG(I)%VECLIN,SEG(I)%Voiso,SEG(I)%Voise
    print *, I,seg(I)%O(1:3)

    if (SEG(I)%VECLIN > nbase .or. SEG(I)%VECLIN < 1 ) then
       print *, " Donnees erronnees pour le segment :", i
       stop
    endif
end do

if(cartograph == 1 .or. cartograph == 3) then
   read(3,*) axe1_carto,axe2_carto
   read(3,*) phi1_carto,phi2_carto
   read(3,*) fichier_carto
endif

close(3)

! Initialisation des couloire de GB
! Attention ce ligne devront etre efface et passe dans un fichier initiale
GB_101_SIZE = 50000
GB_121_SIZE = 50000
GB_111_SIZE = 300

if (nsegm /= IZERO) RETURN

open (3,FILE='../in/jonction',STATUS='OLD')

read(3,*) Plusseg  ! nombre de segemnts initiaux
Nsegm = Nsegm + Plusseg
Plusseg = Izero
read(3,*) boite ; ! taille de la boite den micron
read(3,*) carac1 ! coin du systeme 1
read(3,*) taux1  ! facteur boite / taille de segemnts
read(3,*) npas1
read(3,*) decal1
decal1 = decal1 / 100.0
read(3,*) solli1
solli_sys(syseg(carac1)) = solli1
read(3,*) carac2 ! vis du systeme 2
read(3,*) taux2 ! facteur boite / taille de segemnts
read(3,*) npas2
read(3,*) decal2
decal2 = decal2 / 100.0
read(3,*) solli2
solli_sys(syseg(carac2)) = solli2

!print *, solli_sys(syseg(carac1)),syseg(carac1)
!print *, solli_sys(syseg(carac2)),syseg(carac2)

close(3)

l(:) = (nint( (((boite * 1.D-6 ) / 2.0)/avalue),DPI) / facteur_boite) * facteur_boite

centre(:) = noeud(l,2)

write(*, '(/," Centre de la boite (",3(I10,1x),")")') centre(1:3)


Modur(:) = centre(:)*2

!print *, " centre =", centre
!print *, " modur =", modur

if (etat(centre,noeud(centre,2)) /= 3) stop "centre  mal calcule"

ivis1 = 1 + ((carac1-1)/NBASERED)*NBASERED ! indice du vecteur vis associe

SEG(1)%VECLIN = carac1


! taille du seg en bvd
taille =  real(modur(1),DP)/norivect(bveclin(1:3,carac1))
!print *, "taille =", taille, seg(1)%veclin
seg(1)%norme = ((nint(taille,DPI) * 2) / 2)/taux1 ! taille en vd

! le vecteur transla relie l'origine souhaite du segemnt au centre de la boite
if (modulo(carac1,IDEUX) == 0 .and. (carac1 - ivis1) < 4) then
   l (:) = bveclin(:,ivis1)
elseif (modulo(carac1,IDEUX) == 0 .and. (carac1 - ivis1) > 4) then
   l (:) = bveclin(:,ivis1+4)
else
   l(:) = bvecdep(1:3,carac1)
endif


l(:) = bvecdep(1:3,carac1)

transla (1:3) = npas1*l(1:3) + (1+decal1)*(seg(1)%norme/2)*bveclin(1:3,carac1)

l(1:3) = centre(1:3) - transla(1:3)
transla(:) = noeud (l(:),2)

!print *, "decalage 1/O =",sens1*npas1*l(1:3)
!print *, "transla =",transla (1:3)

if (etat(transla,noeud(transla,2)) /= 3) stop " mauvais transla1"
!if (iprodsca(transla(:),bvecnor(:,carac1)) /= 0 ) stop " mauvaise translation1 dans generer_segment"

seg(1)%O(:) = transla(1:3) ! origine en a


! calculer pour le segment 2 seg(1)%O

SEG(2)%VECLIN = carac2
! norme du seg 2 en bvd
taille =  real(modur(1),DP)/norivect(bveclin(1:3,carac2))
!print *, "taille =", taille, seg(1)%veclin
seg(2)%norme = ((nint(taille,DPI) * IDEUX) / IDEUX)/taux2 ! taille en vd
ivis2 = 1 + ((carac2-1)/NBASERED)*NBASERED ! vecteur vis associe

!if (iprodsca(bveclin(:,ivis1),bveclin(:,ivis2)) > 0) then
!   sens = 1
!else
!   sens = -1
!endif

! le vecteur transla relie l'origine souhaite du segemnt au centre de la boite
if (modulo(carac2,IDEUX) == 0 .and. (carac2 - ivis2) < 4) then
   l (:) = bveclin(:,ivis2)
elseif (modulo(carac2,IDEUX) == 0 .and. (carac2 - ivis2) > 4) then
   l (:) = bveclin(:,ivis2+4)
else
   l(:) = bvecdep(1:3,carac2)
endif

l(:) = bvecdep(1:3,carac2)
!if (modulo(carac2,Ideux) == 0) npas2 = (npas2/ 55) * 55
transla (1:3) = npas2*l(1:3) + (1+decal2)*(seg(2)%norme/2)*bveclin(1:3,carac2)
!print *, "decalage 1/O =",sens2*npas2*l(1:3)
!print *, "transla =",transla (1:3)

l(:) = centre(:) - transla(:)
transla(:)= noeud (l(:),2)

if (etat(transla,noeud(transla,2)) /= 3) stop " mauvais transla2"
!if (iprodsca(transla(:),bvecnor(:,carac2)) /= 0 ) stop " mauvaise translation2 dans generer_segment"
!print *, "2 ",sens2*npas, ivis2

seg(2)%O(:) = transla(:)

! definition des voin : points d'encrage
SEG(1)%Voiso = 0
SEG(1)%Voise = 0
SEG(2)%Voiso = 0
SEG(2)%Voise = 0

end subroutine lire_segments

!##################################################################################
! lecture des donnees particlues por les HC
!##################################################################################

subroutine lire_particules

implicit none

integer (kind=DPI) :: iini,i

open(3,file="../in/particules",STATUS='OLD')

! Nombre de particue
read(3,*) Npar
if(npar == izero) return

! S'il y a des particules, on deplace les vrais segments de nar cases dans le type derive
Do iini = 1, nsegm
    i = nsegm - iini + 1
    SEG(i + npar) = SEG(i)
enddo

Nsegm = Nsegm + npar

! Finalement les particlues
do I=1,Npar
    SEG(i) = SEG(NSEGMAX)
    seg(i)%norme = izero
    read (3,*)	seg(I)%O(1:3),seg(i)%resdep
    seg(i)%veclin = IZERO
    seg(i)%voiso = nsegmax ;   seg(i)%vnno = nsegmax
    seg(i)%voise = nsegmax ;   seg(i)%vnne = nsegmax
end do
close(3)
end subroutine lire_particules


!#############################################################################
!###determination des parametres relies a la configuration initile des segemnts
!#############################################################################


subroutine Configuration

Implicit NONE

INTEGER(DPI) :: I,O(3),Li(3),compt,dislo,isource
REAL(DP) :: X1,Lr(3)
! =======================================================================
! =======================================================================
! ============       verification de l'axe irrationnel axeirra   =========
! Cet axe est tres utile pour comparer les vecteurs parallels
! il ne doit pas etre ni parallele ni perpendiculaires au vecteurs de
! discretisation de la simulation
! on verifie que axeirra n'est pas sur le reseau principal. le test doit
! etre valid pour toutes les structures cristallographique
if (.not. CFC .and. etat(axeirra,noeud(axeirra,2)) == 3 ) then
   print *, "axeirra =",axeirra
   print *, " noeud(axeirra,2)=",noeud(axeirra,2)
   stop " Axeira : n'est pas si irrationnel que ca ?????????"
endif
write (*,'(/," Vitesse unitaire associee a VIS(1)   = ",F10.6," m/s")') &
normdep(1)*avalue/deltat
write (*,'(/," Vitesse unitaire associee a Mixte(2) = ",F10.6," m/s")') &
normdep(2)*avalue/deltat
write (*,'(/," Vitesse unitaire associee a COIN(3)  = ",F10.6," m/s")') &
normdep(3)*avalue/deltat


! les translations elementaires : possible de les reproduire a partire de ces translations
IF (HC .OR. BCC ) THEN
   Do i=1, nbase
       O(1:3) = bveclin(1:3,i)
       if (etat(O,noeud(O,1)) /= 3) then
          print *, " problem du vecteur ligne du segments :",I
          print *, " Veclin : ", bveclin(1:3,i)
          print *, " noed : ", noeud(O,1)
          stop
       endif
   enddo
ENDIF

! pour les test de segments individuels et jonction ......
! si le nombre de segments lu dans segs est nul, le programme lit
! le fichier "jonction" dans lequel on determine la dimension de la boite
! et les caracteristiques de un ou deux segmnets

!==========================================================================================
!============ remise en conformite des dimensions de la boite      ========================
!==========================================================================================

!Mise en conformite des dimension de la boite en function de la crystallo
Modur(:) = (modur(:) / Facteur_boite) * Facteur_boite

! remettre les dimensions de la boite sur les noeuds du reseau
! les calcul montrent que les coordonnees de la boite doivent etre des multiples
! du facteur  Facteur_boit , defini dans le module de constantes
!call segsinfo(" avnt modulo                   ")
if (.not. CFC) then
! si nsegm = 0, alors on appelle generer_segements pour creer une configuration
! pour l'etude des jonction


!==========================================================================================
!============= verification des dimensions e de conformite==================================
!==========================================================================================

   O(1) = modur(1)
   O(2) = 0
   O(3) = 0
   if (etat(O,noeud(O,2)) /= 3) stop " load:mauvais facteur de multiplication selon X"
   O(1) = 0
   O(2) = Modur(2)
   O(3) = 0
   if (etat(O,noeud(O,2)) /= 3) stop " load:mauvais facteur de multiplication selon Y"

   O(1) = 0
   O(2) = 0
   O(3) = Modur(3)
   if (etat(O,noeud(O,2)) /= 3) stop " load:mauvais facteur de multiplication selon Z"

! si le programme n'est pas arrete, c'est que les poits suivants :
! [ modur(1) ; 0 ; 0 ],  [ 0 ; modur(2) ; 0 ] et [ 0 ; 0 ; modur(3) ]
! sont sur le reseau de la simu

! teste ultime
   if (etat(modur,noeud(modur,2)) /= 3) stop " load:le point modur(1:3) n'est pas sur le reseau"

!==========================================================================================
!============ remise des originirne des dislo sur le reseau pricnipal=====
!==========================================================================================
   do i = npar+1,nsegm

       if (SEG(I)%Voiso /= 0) SEG(I)%Voiso = i-1
       if (SEG(I)%Voise /= 0) SEG(I)%Voise = i+1

       if (seg(i)%voiso == 0 .or. seg(i)%voiso == nsegmax) then
          if (HC.or.BCC) then
             O(:) = seg(i)%O(:)
             seg(i)%O(:) =  noeud(O,2)
          endif
       else
          seg(i)%O(:) = seg(seg(I)%voiso)%O(:) + &
          seg(seg(i)%voiso)%norme * bveclin(:,seg(seg(i)%voiso)%veclin)
       endif
   enddo
else
   do i=1,nsegm
       if (seg(i)%voiso == 0 .and. CFC) seg(I)%O(:) = (seg(I)%O(:) / 8)*8

       if (seg(i)%voiso /= 0) then
          seg(i)%O(:) = seg(seg(I)%voiso)%O(:) + &
          seg(seg(i)%voiso)%norme * bveclin(:,seg(seg(i)%voiso)%veclin)
       endif
   enddo
endif

If (HC) then
   do i = 1, Npar
       O(:) = seg(i)%O(:)
       seg(i)%O(:) =  noeud(O,2)
       seg(i)%resdep = seg(i)%resdep * 1.0D6
   enddo
endif




! Calul du nombre total de sources de Frank_Read
write(*, '(/,"Nombre de particules : ",I8," ; nombre de segemnts :",I8)') Npar, Nsegm- Npar

Nombre_FR = IZERO

do i= npar+1,nsegm
    if( seg(i)%voiso == 0) nombre_FR = Nombre_FR + 1
enddo

print *, " "
print *, "*******************************************************************"
write(*, '("   Nombre de source de Frank_Read  : ",I10)') Nombre_FR
print *, " "

isource = izero
do dislo = npar + 1, nsegm
    if(seg(dislo)%voiso /= 0) cycle
    isource = isource + 1
    i  = dislo
    Li(1:3) = Izero
    compt = izero
    do while (i /= 0)
    print *, i,seg(i)%voise
        compt = compt + 1
        Li(:) = Li(:) + seg(i)%norme * bveclin(:,seg(i)%veclin)
        i = seg(i)%voise
        if (compt > nsegm) then
           print *, " boulce infini dans init7: Source FR"
           stop
        endif
    enddo
    Lr(:)  = avalue * real(Li(:),DP)
    x1 = norvect(Lr) * 1D+6
    write(*, '(" Source : ", I4,", ; Systeme :", I6," ;  taille = ", F5.2," microns")') &
    isource, syseg(seg(dislo)%veclin), x1
enddo
print *, "*******************************************************************"

! verification de la compatibilité entre le
if (cartograph == 1 .and. nombre_fr /= 2) then
   print *, " impossible de faire la carto pour N source FR /= 2"
   stop
endif

!#######################################################
! remettre les segemnts dans la boite si ils sont d'hors
!#######################################################
seg(1:nsegm)%o(1)=modulo(seg(1:nsegm)%o(1),modur(1))
seg(1:nsegm)%o(2)=modulo(seg(1:nsegm)%o(2),modur(2))
seg(1:nsegm)%o(3)=modulo(seg(1:nsegm)%o(3),modur(3))
print *, " "
print *, "*******************************************************************"
write(*, '("Dimensions de la boite en a (",3(I10,1x),")")') MODUR(1:3)
write(*, '("Dimensions en microns       (",3(F10.3,1x),")")') MODUR(1:3)*avalue*1D6
print *, "*******************************************************************"
print *, " "


if(krc == 0)  then
   write (*,*) 'KRC == 0. Interdit , By'
   stop
endif
! Un petit teste de verification
if (mod(KISAUVE,KRC).ne.IZERO) then
   write (*,*) 'il faut que le nb de pas avant sauvegarde soit un multiple de krc'
   print *, kisauve, krc
   stop
endif


!==========================================================================================
!==========================================================================================
! =======================================================================
!   Mise en conformite des parametre de la siluation en fonction
!   du mode de deformation choisi
select case(mode_deformation)
case (1)        ! mode de deformation en Contrainte uniaxial imposee
   STRATE	= .false.     ! Clef de control en vitesse de def = cste
   METFUT	= .false.     ! Clef de control en metallofute
   FATIG		= .false. 	  ! Clef de control en fatigue
   FSIG = 0   ! Signe initial (fatigue)
   EPSMAX = 0 ! Amplitude maximale de deformation(fatigue)

case (2)        ! mode de deformation uniaxial imposee
   STRATE	= .True.     ! Clef de control en vitesse de def = cste
!    SigmaPoint    sert dans ce mode pendant l'interval de calculde moyenne d'epsilon
   METFUT	= .false.     ! Clef de control en metallofute
   FATIG		= .false. 	  ! Clef de control en fatigue
   EPSMAX = 0! Amplitude maximale de deformation(fatigue)
   FSIG   = 0 ! Signe initial (fatigue)

case (3)        ! mode de deformation en Fluage (Sigma0 initial uniaxiale)
   STRATE	= .false.     ! Clef de control en vitesse de def = cste
   SigmaPoint = 0.0
   METFUT	= .false.     ! Clef de control en metallofute
   FATIG		= .false. 	  ! Clef de control en fatigue
   DeltaEpsilon = 0.0      ! interval ou on calcul la moyenne de deformation
   EPSMAX = 0 ! Amplitude maximale de deformation(fatigue)
   FSIG   = 0 ! Signe initial (fatigue)

case (4)        ! mode de deformation en Fatigue
   STRATE	= .false.     ! Clef de control en vitesse de def = cste
   METFUT	= .false.     ! Clef de control en metallofute
   FATIG		= .TRUE. 	  ! Clef de control en fatigue
   DeltaEpsilon = 0     ! pour tester

case (5)        ! mode de deformation en MetalloFute
   STRATE	= .false.     ! Clef de control en vitesse de def = cste
   METFUT	= .TRUE.     ! Clef de control en metallofute
   FATIG		= .false. 	  ! Clef de control en fatigue
   EPSMAX = 0 ! Amplitude maximale de deformation(fatigue)
   FSIG   = 0 ! Signe initial (fatigue)
end select


!==========================================================================================
!==========================================================================================
!==========================================================================================

! Echelle spatiale
BDIVPA = BDIVA/PII

!#############################
! GEOMETRIE DE L'ECHANTILLON
!#############################
xloclmax = minval(ModuR(:))/Itrois	!*** LONGUEUR maximal d'un segmet


print *, " augmentation des longueurs de discretisation pour les system non actifs : "
print *, " remise à 0.5 micron"

compt = Izero
do I=1,NBASE
    LOMA(I) = INT(XLOMAX/NormLin(I),DPI) !*** PLUS RAPIDE
    LOCLMA(I) =INT(XLOCLMAX/NormLin(I),DPI) !*** PLUS RAPIDE
    JSLOMA(I) =INT(XLOMAX/(NormLin(I)*4.0),DPI) !*** PLUS RAPIDE
    if (solli_sys(syseg(i)) < PRECISION) then
!       print *, i, loma(i),INT(Xlomax_nact/NormLin(I),DPI)
       LOMA(I) = INT(XLOMAX_nact/NormLin(I),DPI) !*** PLUS RAPIDE
       compt = compt + IUN
    endif
enddo
if (compt /= Izero) then
   write(*,'(/," Discretisation de systemes non actifs pour :", I8, " systemes.")') compt/nbasered
endif

write(*, '(/,"Longueur de discretisation (en vis   elementaire) = ",I8)') loma(1)
write(*, '(/,"Longueur de discretisation (en mixte elementaire) = ",I8)') loma(2)
write(*, '(/,"Longueur de discretisation (en coin  elementaire) = ",I8)') loma(3)


if (xlomax >= xloclmax) &
stop " longueur de discretisation ne doit pas depasser la longueur maximale tolerer"

VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
BDIVAV = BDIVA/VOLUME
write(*, '(/,"Norme maximale des segemnts VIS = ",I8)')  loclma(1)
write(*, '(/,"Norme maximale des segemnts MIXTE = ",I8)')  loclma(2)
write(*, '(/,"Norme maximale des segemnts COIN = ",I8)')  loclma(3)

x1 = (dfloat(modur(1)))**2 + (dfloat(modur(2)))**2 + (dfloat(modur(3)))**2
x1 = Dsqrt(x1)
LDIAGOMAX = 100.0/x1

!*******************************
! Quelques initialisations utils
!*******************************
! Longueur de discretisation de la courbure en unite a ???
ICROSMAX = INT(XLOMAX/4.,DPI)

! Tailles caracteristiques de boucles de dislocation pour les procedures
! de netoyage, on defini deux longueur utils
MICROBOU=INT(XLOMAX/3.,DPI)
MINIBOU=INT(XLOMAX*2./3.,DPI)

!************************************
! MAILLAGE POUR LA METHODE DES BOITES
!************************************
IF(LINTEN) THEN
! la taille des boites et leur mise a jour est ici forcee au maximum
! pour eliminer du calcul des boites tous les segments
   write (*,*) 'MODE TENSION DE LIGNE SEULE'
   ISPLIT = 1	!*** NOMBRE DE BOITE
   KRC    = 1	!*** NOMBRE D'ITERATIONS SANS ACTUALISATION
ENDIF

VolMetBoi = Volume/isplit
x1 = 1.0/3.0
aBoi = VolMetBoi**(x1)
nbxyz(:)=real(ModuR(:),DP)/aboi

Inbxyz(:)=int(nbxyz(:),DPI)

!*** Nb de boites par dimension impair pour pouvoir renconstruire un config
!    autour d'une boite (ne semble pas vital en mode bord libre)

!--- BL ---
if ((Inbxyz(1)/2)*2.eq.Inbxyz(1).and.nbxyz(1).ge.Inbxyz(1)) Inbxyz(1)=Inbxyz(1)+1
if ((Inbxyz(1)/2)*2.eq.Inbxyz(1).and.nbxyz(1).lt.Inbxyz(1)) Inbxyz(1)=Inbxyz(1)-1
if ((Inbxyz(2)/2)*2.eq.Inbxyz(2).and.nbxyz(2).ge.Inbxyz(2)) Inbxyz(2)=Inbxyz(2)+1
if ((Inbxyz(2)/2)*2.eq.Inbxyz(2).and.nbxyz(2).lt.Inbxyz(2)) Inbxyz(2)=Inbxyz(2)-1
if ((Inbxyz(3)/2)*2.eq.Inbxyz(3).and.nbxyz(3).ge.Inbxyz(3)) Inbxyz(3)=Inbxyz(3)+1
if ((Inbxyz(3)/2)*2.eq.Inbxyz(3).and.nbxyz(3).lt.Inbxyz(3)) Inbxyz(3)=Inbxyz(3)-1
!--- BL ---

aboixyz(:)=real(ModuR(:),DP)/real(Inbxyz(:),DP)
isplit=inbxyz(1)*inbxyz(2)*inbxyz(3)

Modui(:) = Inbxyz(:)

!*** Calcul de quantites utile pour la methode des boites
ModuXLSPLIT(:) = DFLOAT(ModuR(:))/DFLOAT(ModuI(:))	!*	Taille des boites
ModuXLSMOIT(:) = ModuXLSPLIT(:)*HALF								!*	Centre des boites
ModuXLINV(:)   = UN/ModuXLSPLIT(:)									!*	Facteur pour determiner la boite


! =======================================================================
!  Calcule des parametres des lois de vitesses de la simulation
! =======================================================================

print *, " "
print *, " "
print *, "Nombre des lois de vitesse lues : NLV  =",NLV

do i = 1, NLV_Max

    Loi(i)%arrhenius = Loi_brute(i)%arrhenius
    if (Loi(i)%arrhenius == 2) then       ! modele de double decrochement
! exposant p de la loie de vitesse (pour les vis)
       Loi(i)%coef_p = Loi_brute(i)% coef_p
! exposant q de la loie de vitesse (pour les vis)
       Loi(i)%coef_q = Loi_brute(i)% coef_q
! cission critique a T = 0 K
       Loi(i)%Tau0eff =  Loi_brute(i)%tau0eff
! cission athermique
       Loi(i)%Tauath =  Loi_brute(i)%tauath
! calcul de la constante preexponentielle
       Loi(i)%H0 =  Loi_brute(i)%V0 / Loi_brute(i)%L0
! calcul de l'energie d'activation effective : normalise afin d'etre
! directement utilisable pour le calcule de la vitesse
       Loi(i)%DeltaG0eff = -1.0 *(Loi_brute(i)%deltag0/Boltzmann/TEMPERATURE)
!       print *, " loi =", i," DGefff = ", Loi(i)%DeltaG0eff
    elseif (Loi(i)%arrhenius == 1) then! activation de type obstacle locel
! exposant p de la loie de vitesse (pour les vis)
       Loi(i)%coef_p = Loi_brute(i)% coef_p
! exposant q de la loie de vitesse (pour les vis)
       Loi(i)%coef_q = Loi_brute(i)% coef_q
! cission critique a T = 0 K
       Loi(i)%Tau0eff =  Loi_brute(i)%tau0eff
! cission critique a T = 0 K
       Loi(i)%Tauath =  Loi_brute(i)%tauath

       Loi(i)%H0 =  Loi_brute(i)%V0/avalue
       if((HC .or. BCC) .and. temperature < 300.0) then
! 1.D6*Exp(-temperature/30.29) = 50 pour T = 300K
! 1.D6*Exp(-temperature/30.29) = 1000000 pour T = 1K
! autrement dit si Loi(i)%H0 est la meme pour vis et coin, le rapport de vitesse vaut
! 50 et 1000000 pour les differentes parametres
          Loi(i)%H0 = Loi(i)%H0 * 1.D6*Exp(-temperature/30.29)
          write(*, '(/,"Le rapport des vitesse a la CRSS = ",F10.2)') &
           Loi_brute(i)%V0/Loi_brute(1)%V0*1.D6*Exp(-temperature/30.29)
       endif

! calcul de la constante preexponentielle
! on considere que la vitesse n'est pas proportionnelle
! a la longueur (activation thermique autre que la friction de reseau)

! calcul de l'energie d'activation effective : normalise afin d'etre
! directement utilisable pour le calcule de la vitesse
       Loi(i)%DeltaG0eff = -1.0 *(Loi_brute(i)%deltag0/Boltzmann/TEMPERATURE)
!       print *, " loi =", i," DGefff = ", Loi(i)%DeltaG0eff

    elseif (Loi(i)%arrhenius == 0) then
       Loi(i)%BF   = Loi_brute(i) % Coef_visqueux*avalue
       Loi(I)%TAU0 = Loi_brute(i) % Max_friction
    else
       stop " Type de loi de vitesse inconnu !!!!"
    endif
! procedure teporaire: definition de loi de mobilite pour la relaxation
!    Loi(NLV_MAX)%bf = 5.D-4 *avalue
!    Loi(NLV_MAX)%TAU0 = zero
!    print *," loi_brute    ",loi_brute(i)
!    print *," loi    ",loi(i)
End do

If ( NLV > NLV_max) then
   print *, " nombre de lois de vitesse trop grand"
   stop
endif
call flush (6)

RETURN
end subroutine configuration




!************************************************************************
! Ce sous-programme ecrit dans le fichier (tableux.cfc,cc,cs, etc.)
! les tableaux : quelue tableaux utiles pour la simulation
!  Ghiath Monnet  09/02/01
!************************************************************************
subroutine write_tableaux

implicit none

integer (kind=DPI) :: i,j


60 format(9I5)

open (13,FILE=fichtab,STATUS='UNKNOWN')

write (13,*) ' les tableaux  : '
write (13,*) ' '
write (13,*)'    I    J   axei axej sysobst'
write (13,*) ' '

do i = 1,nbase
    do j=1,nbase
        write (13,60) i,j,axejonci(i,j),axejoncj(i,j),sysobst(i,j)
        if (modulo(j,IHUIT) == 0) then
           write (13,*) '    -----------------------------------------------------'
        endif
    enddo
    write (13,*) '    -----------------------------------------------------'
enddo

close(13)

end subroutine write_tableaux



!#########################################################################
!# Tabulation du calcul de la tension de ligne.                          #
!############################################################# 05/11/98 ##
subroutine TABTLIN

!*** A adapter au cas des mixtes !!!
!tabulation de la tension de ligne dans l approximation du gradiant de
!l energie de ligne. on supose les segments a 90 deg

implicit none

real(kind=DP)	:: FA2
integer (kind=DPI) 	:: I

!*** Constantes utiles au calcul des forces
FA	= (UN-DPOISS)
FACT	= -1*HALF/(UN-DPOISS)
FICTA	= PERCEN*QUART/PII/FA*BDIVA*BDIVA

!*** Tableau contenant la partie log du calcul de la tension de ligne
FA2	= (DEUX-DPOISS)
do I=1,10000

!*** Attention : en fonction de la formule commentee ci dessous, la
!*** solution est soit isotrope soit anisotrope

    DTABV(I)	= FA*(DLOG(DEUX*DFLOAT(I)*NormLin(3) &
    /(BDIVA))+UN)-FA2
! DTABV(I)	= (DLOG(DEUX*DFLOAT(I)*NormLin(3) &    ! Cas anisotrope
!             	  /(BDIVA))+UN)-FA2
    DTABC(I)	= FA*(DLOG(DEUX*DFLOAT(I)*NormLin(1) & ! Cas isotrope
    /(BDIVA))+UN)-FA2
enddo

end subroutine TABTLIN



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
deallocate (DEVSEG)
deallocate (DEVSEGtemp)
deallocate (CONEC)
deallocate (DECONEC)
deallocate (SC)
deallocate (CoefSC)
deallocate (AXEJONCi)
deallocate (AXEJONCj)
deallocate (SYSOBST)
deallocate (SYSCONNEC)
deallocate (LIEN)
deallocate (LIENS)
deallocate (NBROT)
deallocate (GDROT)
deallocate (INVCONEC)
deallocate (PartVis)
deallocate (SYSEG)
deallocate (LOMA)
deallocate (JSLOMA)
deallocate (LOCLMA)
deallocate (FAC1DEP)
deallocate (FAC2DEP)
deallocate (MODDEPLA)
deallocate (VECNORLIN)
deallocate (VECNORDEP)
deallocate (NORMLIN)
deallocate (NORMDEP)
deallocate (RAUSYS)
deallocate (CRANSYS)
deallocate (GAMSYS)
deallocate (SchmidSys)
deallocate (TrAppSys)
deallocate (TrIntSys)
deallocate (GAMMADOTSYS)
! nouveau

deallocate (Bf)
deallocate (Tau0)

end subroutine Desallocation


!#############################################################################
!# Le programme appele dans le main pour charger touts les donnees         #
!################################################################# 05/11/98 ##
subroutine LOAD

implicit none



!#########################################################################
!# Parametres de controle : chargement et reconstitution des donnees     #
!# physiques ou topologiques utile a la simulation                       #
!#########################################################################
call lire_donnees


if (norvect(Z) < 0.01) then
   UNIAXIALE = .FALSE.
   stop " procedure de chargement tensoriel n'est pas ajour avec les facteurs de schmid"
   call Lire_tenseur
else
   UNIAXIALE = .TRUE.
endif

call parametrage


!#############################################################################
! Initialisation des fichiers en fonction du mode de demarage de la simulation
!#############################################################################
call initialisation_fichiers

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
call generer_basedevecteurs


! Ecriture de la basee de vecteurs de discretisation format brut et soigne
call write_basedevecteurs

call initialiser
!#######################################################################
! tabulations annexes decrivant les correlations de la base des vecteurs
!#######################################################################
call tabulation

! Ecriture des tableaux de facteur de deplacement dans les simulation
call write_facdep

if(cartograph == 2) THEN
   call carto
   print *, " "
   print *, " "
   print *, " "
   print *, " "
   print *, " "
   Print *, "GENERATION DES FICHIERS DE CARTO DANS ../IN/CARTO/ TERMINEE"
   Print *, "SORTIE NORMALE DU PROGRAMME"
   PRINT*, "POUR EXECUTER DD, VEUILLEZ METTRE LA CLEF GENERER_CARTO SUR .FALSE."
   print *, " "
   print *, " "
   print *, " "
   print *, " "
   STOP
ENDIF

call lire_segments

if (HC) call lire_particules


call configuration
if(cartograph == 5) THEN
   call double_config
   Print *, "Dedoublement de la configuration initial"
ENDIF



! ecriture du fichier tableaux.cfc (.cc,.cs,...) contentant les tableaux :
! axejonci, axejonj, caj (crees par tabul_annexes)
call write_tableaux



end subroutine LOAD

end module INIT

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================


