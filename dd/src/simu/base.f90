
!===================================================================================================
!========================    DEBUT    MODULE   "BASE"  =============================================
!===================================================================================================

include '00copyright.f90'

!> \brief ???
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English

Program BASE

use constantes
use VARGLOB
use BRICAMAT

implicit none

real(kind=DP) :: facteurcoin(10),facteurvis(10),facteurmixt(10)
real(kind=DP) :: projection,depsurprojecvis,preci,echelle
real(kind=DP) :: projecsurdep,depsurprojeccoin,depsurprojec
real(kind=DP) :: dep,valeur

integer(DPI)  :: matrix(3,3),vecteur(3), plan(3),i,j,inverse,correcvis
integer(DPI)  :: itemp,vis_dir(3),coin_dir (3),plan_ref(3), direction_ref (3),ivo,ive,mixte(3)
integer(DPI)  :: facteur
integer(DPI)  :: correcbase1,correcdep1,correcbase2,br
integer(DPI)  :: indicesys,limitej,correcdep2
integer(DPI)  :: facteur1, facteur2,depmixte(3),correcbase,vistemp(nbasered),cointemp(nbasered)
integer(DPI)  :: ktemp,syst,Ivec(4)
character     :: carac*1,materiau*60, control*60,fichier*60,struc*3

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
close (1)

open(1,file="../in/"//materiau,STATUS='OLD')

!**   Module de cisaillement MU a T = 0 K
read(1,*) valeur

!**   pente de la courbe Module de cisaillement = fonction (Temperature)
read(1,*) valeur

!*** Module de POISSON
read (1,*) valeur

!**  tauIII = contrainte de debut du stade III (restauration)
read(1,*) valeur
read(1,*) valeur
read(1,*) valeur

!**  lecture de la structure cristallographique
read(1,*) struc(1:3)

!***  Coefficients de friction visqueux BF=5.E-5(Pa.s)/(mu=42GPa)
!**   Friction seche (Peierls) TAU0*b=b*1.0D0(MPa)/(mu=42GPA)
read (1,*) Nb_slip_types

syst=0
do itemp = 1, Nb_slip_types
    read (1,*) valeur        ! norme vecteur Buergers en Angstrom de la famille
    read (1,*) Slip(itemp)%planes(1:3)
    read (1,*) Slip(itemp)%directions(1:3)

    if (etat( Slip(itemp)%planes(1:3),Slip(itemp)%directions(1:3)) /= 0) &
      stop "Erreur de donne crystallo : vis n'est pas perpendiculaire au plan"

    read (1,*) Slip(itemp)%Nsystemes   ! nb de systemes de glissement

! Chaque systeme de glissement a sa norme vecteur de vecteur de Burgers
    do ktemp=(itemp-1)*syst+1,itemp*Slip(itemp)%Nsystemes
        Slip(ktemp)%VecBurgers=valeur
    enddo
    syst=Slip(itemp)%Nsystemes
    read (1,*)   Ivec(1:4)
ENDDO

! lecture du nombre de lois de vitesse a lire
! Attention par defaut la loi n 1 est celle de la relaxation

close (1)

!########################################################INIT
!# Chargement des parametres de control de la simulation
!########################################################
open(2,file="../in/"//control,STATUS='OLD')

read(2,*) Mode_deformation  ! mode de deformation dans la simulation
read(2,*) cartograph        ! clef pour la generation de carto de reaction locales

read (2,*) echelle     ! longueur de reference de la simulation exprime en b

close(2)

!==========================================================================================
! =======================================================================
! Definition de tous les parametres relies a la cristallographie
! FACTEUR_CRIS le rapport entre le modul du vecteur vis (ex. [440] et le module du
! vecteur de Burgers dans le repere cristallo(ex. 1/2 [110])

! Initilisation des clefs choix de la cristallo
DC = .false.; CS = .false.; BCC = .false.; CFC = .false.; HCP = .false.; ORT = .false.; MGO = .false.
If (struc(1:2) == 'CS') then
   CS = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_CS'
   facteur_cris = FACTEUR_CS
elseif (struc(1:3) == 'CFC')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_CFC'
   CFC = .true.
   facteur_cris = FACTEUR_CFC
elseIf (struc(1:3) == 'BCC')   then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_BCC'
   BCC = .true.
   facteur_cris =  FACTEUR_BCC
elseIf (struc(1:3) == 'HCP')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_HC'
   HCP = .true.
   facteur_cris =  FACTEUR_HC
elseif (struc(1:3) == 'ORT') then
   ORT = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_ORT'
   facteur_cris = FACTEUR_ORT
elseif (struc(1:3) == 'MGO') then
   MGO = .true.
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_MGO'
   facteur_cris = FACTEUR_MGO
elseif (struc(1:3) == 'DC')  then
   fichier = '../in/CRISTALLOGRAPHIE/ROTATIONS_DC'
   DC = .true.
   facteur_cris = FACTEUR_DC
else
   stop " structure crystallographique non identifiee ,  BY"
endif
if(.not. CS .and. .not. BCC .and. .not. HCP .and. .not. CFC .and. .not. ORT .and. .not. MGO.and. .not. DC) &
   stop " la constante angle_vis n'est pas definie ?, By !"
!==========================================================================================
!==========================================================================================
!==========================================================================================

write(*, '(" Le Facteur homogeneisation spatiale = ", I5)') facteur_cris
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

!################################################################# 05/11/98 ##
!                                 Allocation
!################################################################# 05/11/98 ##

allocate (BVECLIN(3,NBASE))
allocate (BVECDEP(3,NBASE))
allocate (BVECNOR(3,NBASE))
allocate (VECNORLIN(3,NBASE))
allocate (VECNORDEP(3,NBASE))

!################################################################# 05/11/98 ##

! simplification d'ecriture
br = nbasered
limitej = 100

! facteurvis et coin traduisent le fait que les directions vis et coins initiales
! doivent etre des multiple fixes par la cristallo
if (HCP) then
   facteurvis (:)  = (/TROIS,TROIS,TROIS,zero,zero,zero,zero,zero,zero,zero/)
   facteurcoin (:) = (/Quatre,un,un,zero,zero,zero,zero,zero,zero,zero/)
   facteurmixt (:) = (/un,half,half,zero,zero,zero,zero,zero,zero,zero/)
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
!*  chose avec le deuxieme vecteur coin associe au meme vecteur vis soit 16
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

        if (HCP .and. modulo(i,IDEUX) == 0) plan (:) = -1 * Plan (:)
        if (HCP) inverse = 1

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
        if (HCP .or. ORT) vis_dir(:) = int(float(vis_dir(:)) * facteurvis(indicesys))
        if (BCC .or. HCP .or. ORT) coin_dir(:) = NINT(coin_dir(:) * facteurcoin(indicesys))
! le but etant de recuperer un vecteur Burgers et C entiers
!      On obtent alors : [ 0 3 3] a la place de 1/2[0 1 1]
!       et une translation selon C = [4 4 4] au lieu de 4/3[1 1 1] = coin_dir
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! nouvelle procedure pour optimiser les vecteurs de deplacement dex mixtes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pour la mixte m1 = vis + coin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        preci = norivect(vis_dir) * numtols

! signification des parametres:
! correcbase1: le facteur multiplicateur de la base entiere pour conformite avec la mixte1
! correcbase2: le facteur multiplicateur de la base entiere pour conformite avec la mixte2
! correcbase: le facteur multiplicateur de la base entiere pour les deux mixte
! correcdep1: facteur multiplicateur du deplacement de m1 pour comptabilite avec correcbase
! correcdep2: facteur multiplicateur du deplacement de m2 pour comptabilite avec correcbase
        correcbase1 = 1
        correcdep1 = 1
        correcdep2 = 1

! calul de du vecteur ligne de m1
        if (.not. BCC) then
           mixte(1:3) = int(facteurmixt(indicesys)*float(coin_dir(1:3) + vis_dir(1:3)))  ! vecteur mixte  1
        else
           mixte(1:3) = int(facteurmixt(indicesys)*float(coin_dir(1:3) - 2*vis_dir(1:3)))  ! vecteur mixte  2
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
              print*," projection vis > projection coin: cas non prevu, desole. BY"
              stop
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
        mixte(1:3) = int(facteurmixt(indicesys)*float(coin_dir(1:3) - vis_dir(1:3)))  ! vecteur mixte  2
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
              print*," projection vis > projection coin: cas non prevu, desole. BY"
              stop
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
           print *, " correcbase1 =",correcbase1
           print *, " correcbase2 =", correcbase2
           print *, " attention il faut choisir normalement le plus petit diviseur commun descorrecbase"
           correcbase = correcbase1 * correcbase2
           correcdep1 = correcdep1 * correcbase2
           correcdep2 = correcdep2 * correcbase1
        endif

! homogeneisation des longueurs entre les differents type de systems de glissement
        if (HCP .and. indicesys == 1) then
           correcbase = correcbase * 2
           correcdep1 = correcdep1 * 2
           correcdep2 = correcdep2 * 2
        endif

        if (HCP .and. indicesys > 1) then
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

        if (MGO .and. indicesys == 1) then
           correcbase = correcbase * 2
           correcdep1 = correcdep1 * 2
           correcdep2 = correcdep2 * 2
        endif

        if (MGO .and. indicesys > 1) then
           correcbase = correcbase * 6
           correcdep1 = correcdep1 * 3
           correcdep2 = correcdep2 * 3
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
            if (modulo(j,IDEUX)== 0) bveclin(1:3,itemp) = int(float(bveclin(1:3,itemp))*facteurmixt(indicesys))
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
               facteur = int(facteurcoin(indicesys)*correcbase)
            else
               facteur = int(facteurvis(indicesys)*correcbase)
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

! boucle de verification que tous les indices de tous les vecteur ligne doiven,t etre paire
! sinon, il y a une erreur de cacul de centre des segemnts dans force

itemp = IUN
do i = 1,NTSG*br
    vecteur(:) = bveclin(:,i)
    if((modulo(vecteur(1),ideux) /= izero) .or. &
       (modulo(vecteur(2),ideux) /= izero) .or. &
       (modulo(vecteur(3),ideux) /= izero)) itemp = IDEUX
enddo
if(itemp > iun) then
   do i = 1,NTSG*br
       bveclin(:,i) = IDEUX *  bveclin(:,i)
   enddo
endif


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
        if (HCP .and. i > 24 .and. modulo(i,IDEUX) ==1) bvecdep(:,i) = bvecdep(:,i)/2
        ivo = ivo + j*br   ! vrais indice du voisin en o
        ive = ive + j*br   ! vrais indice du voisin en e

62      format(I5,' :  (',3(I3,1X),') [',3(I5,1X),']. Dep. : [',3(I5,1X),']')
        vecteur(:) = bveclin(:,i)
        if((modulo(vecteur(1),ideux) /= izero) .or. &
           (modulo(vecteur(2),ideux) /= izero) .or. &
           (modulo(vecteur(3),ideux) /= izero))     &
            stop 'Un des vecteur est impaire - interdit pour force'
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
        preci = dep * numtols

! 3) si la projection n'est pas egale au deplacment de I, on change bvecdep
        if (.not.egalite(projection,dep)) bvecdep(:,I) = nint(Bvecdep(:,I)*(projection/dep),DPI)

! 4) Verification definitive de projection = bvecdep
        projection = abs(real(iprodsca(bveclin(:,ivo),bvecdep(:,i)),DP)/norivect(bvecdep(:,i)))
        dep = norivect(bvecdep(1:3,i))
        preci = dep * numtols
        if (.not.egalite(real(nint(projection/dep),DP),projection/dep)) then
           write (*,62) ivo, bvecnor(1:3,ivo), bveclin(1:3,ivo),bvecdep(1:3,ivo)
           write (*,62) i, bvecnor(1:3,i), bveclin(1:3,i),bvecdep(1:3,i)
           write (*,62) ive, bvecnor(1:3,ive), bveclin(1:3,ive),bvecdep(1:3,ive)
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

61 format(I3,': (',3(I5,1X),')[',3(I6,1X),']. Dep: [',3(I6,1X),']')
63 format(I5,"  ",3(1X,I6),"  ",3(I6,1X),"      ",3(I6,1X))

materiau = '../in/CRISTALLOGRAPHIE/BVD.'//struc(1:3)

open (15,FILE=materiau,STATUS='UNKNOWN')  ! fichier de sortie
! d'abord on ecrit la liste brut des vecteur de discretisation
! pour l'utilisation d'autre programmes
write (15,*) nbase       ! unique de programme : lu par camera ....
write (15,*) "   "
do i = 1,nbase
    write (15,63) i, bvecnor(1:3,i), bveclin(1:3,i),bvecdep(1:3,i)
    if(modulo(i,ihuit) == izero) write (15,*) "   "
enddo
write(15,*) "  " ; write(14,*) "  " ; write(14,*) "  "
if (HCP) write(15,*)        "Indice        Plan        ligne (/22)                 deplacement "
if (BCC) write(15,*)        "Indice        Plan        ligne (/6)                 deplacement "
if (.not. HCP) write(15,*) "Indice        plan           ligne                    deplacement "

do itemp = 1,NTSG
    do j = 1,8
        i = (itemp - 1) * nbasered + j
        if (HCP) then
           write (15,61) i, bvecnor(1:3,i), bveclin(1:3,i)/22,bvecdep(1:3,i)
        elseif (BCC) then
           write (15,61) i, bvecnor(1:3,i), bveclin(1:3,i)/6,bvecdep(1:3,i)
        else
           write (15,61) i, bvecnor(1:3,i), bveclin(1:3,i),bvecdep(1:3,i)
        endif
    enddo
    write (15,*) "   "
enddo

close(15)


!#############################################################################
!# Desallocation des tableaux de la simulation (allacated by dimentionner_base)
!################################################################# 05/11/98 ##
deallocate (BVECLIN)
deallocate (BVECDEP)
deallocate (BVECNOR)
Deallocate (VECNORLIN)
Deallocate (VECNORDEP)

write (*,*) "========================================================= "
write (*,*) "========================================================= "
write (*,*) " Successful generation of the discretization vectors "
write (*,*) " Written in : ", materiau
write (*,*) "========================================================= "
write (*,*) "========================================================= "

end program BASE

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
