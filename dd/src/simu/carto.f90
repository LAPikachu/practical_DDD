
!===================================================================================================
!========================    DEBUT    MODULE   "CARTOGRAPHIE"  =====================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the junction maps.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module cartographie

use constantes
use bricamat
use varbase
use varglob


implicit none


Real(DP), parameter :: pisurdeg = PII/180.0
Real(DP), parameter :: degsurpi = 180.0/PII
integer (DPI) :: ivis1,vis1(3),ivis2,vis2(3), echan1,echan2,VecS1(3,100),VecS2(3,100)
integer(DPI) :: decalage,axe(3),Iaxe1,Iaxe2

integer(DPI), parameter :: LengthGD = 8             !< Variables used in the case includeGD, seg length
integer(DPI), parameter :: VeclinGD = 3             !< Variables used in the case includeGD, seg veclin

real(DP) :: phi1Axe,phi2axe,phi1(100),phi2(100)

contains

!##########################################################################
!# subroutine SELECAB : procedure fournissant les deux vecteursde base    #
!#            de part et d'autre de la direction de la source I.          #
!#            sys indiqur le system (1 ou 2) ,ia1 et ib1 sont les indices #
!#            dans la base des vecteurs selectionnes                      #
!###########################################################14/04/99#######
subroutine selecAB(I,sys,iA,iB)

implicit none

integer(DPI) :: I, Ia,ib,sys
integer(DPI) :: X(3), Va(3),Vb(3),v,k,plan(3),j,phi
if (sys == 1) then
   v = ivis1
   phi = int(phi1(I)*degsurpi,DPI)
   X(:) = reduire_vec(VecS1(:,I))
else
   phi = int(phi2(I)*degsurpi,DPI)
   v = ivis2
   X(:) = reduire_vec(VecS2(:,I))
endif

iA = izero ; ib = izero
plan(:) = bvecnor(:,v)

Do J = 0 , 7
    k = J + 1
    if(k == IHUIT ) k = izero
    Va(:) = reduire_vec(bveclin(:, v + K))
    Vb(:) = reduire_vec(bveclin(:, v + j ))
    ! si X est entre Va et Vb, alors, les produits vectoriels ci-dessous doivent etre dans la meme direction du plan
    ! de glissement
    if ((etat(prodivec(Va,X),plan) >= izero) .and. (etat(prodivec(X,Vb),plan) >= izero)) then
       iA = v + k ; ib = v + j
       exit
    endif
enddo
if(ia == izero .or. ib == izero) then
   print *, " system = ", sys, " source =", I
   print *, " The surch of parallel BVD vectors in the two sys failed, sorry"
   stop
endif

end subroutine selecAB

! ---------------------------------------------------------------------------------------

!##########################################################################
!# subroutine decompo : decompose le vecteur S en n* A + m*B. il renvoie  #
!#            donc n,m,l'erreur sur la longueur et sur l'angle            #
!###########################################################14/04/99#######
subroutine decompo(S,A,B,n,m,deltal,deltaphi)

implicit none

integer(DPI) :: S(3), A(3),B(3),n,m,res(3),cs
real(DP) :: normeS,normeA,normeB,deltal,deltaphi,temp,alpha,beta,proja,projb

n = 0
m = 0
deltal = 0.0
deltaphi = 0.0

!if (iprodsca(S,A) < 0 ) stop "decomp:erreur: angle optu entre A et S, pas normal"
!if (iprodsca(S,B) < 0 ) stop "decomp:erreur: angle optu entre B et S, pas normal"
normeS = norivect(S)
normeA = norivect(A)
normeB = norivect(B)

if(normeS < un) stop " vecteur S nul : decompo impossible"
if(normeA < un) stop " vecteur A nul : decompo impossible"
if(normeB < un) stop " vecteur B nul : decompo impossible"

if(normeS < normeA) stop " vecteur S < A : decompo impossible"
if(normeS < normeB) stop " vecteur S < B : decompo impossible"

cs = 0
if(S(1) /= 0) then
   cs = 1
elseif(S(2) /= 0) then
   cs = 2
else
   cs = 3
endif

! cas ou S est // a A
if (abs(etat(S,A)) > 1 ) then
   m = 0
   n = nint(dble(S(cS))/dble(A(cS)),DPI)
   ! cas ou S est // a B
elseif(abs(etat(S,B)) > 1 ) then
   n = 0
   m = nint(dble(S(cS))/dble(B(cS)),DPI)
else
   alpha = dacos(icose2v(S,A))
   beta = dacos(icose2v(S,B))
   proja = abs(dsin(alpha)*normeA)
   projb = abs(dsin(beta)*normeB)
   n = INT((normeS/normeA)/(dcos(alpha) + dsin(alpha)/dtan(beta)),DPI)
   m = INT((normeS/normeB)/(dcos(beta) + dsin(beta)/dtan(alpha)),DPI)
endif

Res(:) = n * A(:) + m * B(:)
temp = norivect(Res)
deltal = dabs(NormeS - temp) / NormeS
deltaphi = dacos(icose2v(S,Res))

end subroutine Decompo

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
!##########################################################################
!# subroutine decoupsource: discretise la source (n* A + m*B) en un       #
!#            nombre de segemnts (nseg) esgueur et sur l'enage            #
!#            donc n,m,l'erreur sur la longueur et sur l'enage            #
!###########################################################14/04/99#######
subroutine decoupsource (iA1,iB1,n,m,decoup,VLS,lonS,Nseg,includeGD)

implicit none

integer(DPI)  :: ia1,ib1,A1(3),B1(3),n,m,res(3),rn,rm,dn,dm,sys1,sys2,Nseg
integer(DPI)  :: half_res(3)          !< The half vector res
integer(DPI)  :: cum_vec(3)           !< The cummulated source length
integer(DPI)  :: VLS(1000),lons(1000),i,lomaa1,Lomab1,decoup
integer(DPI)  :: j
real(DP)      :: temp
real(DP)      :: ProjLongTest         ! The projected length of cum_vec
real(DP)      :: ProjLongRes
logical, intent(in) :: includeGD        !< The key parameter to activate the possibility to include a GD segment
                                        ! at the center of the source

VLs(:)  = 0
Lons(:) = 0
Nseg    = 0
sys1    = 0
sys2    = 0

! The two reference directions used to defined the segment direction
! At this stage the finale segment is discretised with only two segments
a1(:) = bveclin(:,ia1)
b1(:) = bveclin(:,ib1)

! The accumulated length of line
res(:)      = n *a1(:) + m * b1(:)
half_res(:) = res(:) / 2
ProjLongRes = (res(1)*res(1) + res(2)*res(2) + res(3)*res(3)) / 2

if (n == 0 .and. m == 0) stop " decoupsource impossible : n and m are nuls"

! The discretization reference length for the two reference direction
temp = norivect(res)/decoup
lomaa1 = nint(temp/norivect(a1),DPI)
Lomab1 = nint(temp/norivect(b1),DPI)


! The segment is parallel two one of the reference direction (n or m = 0)

if (n*m == 0) then

  ! The two trivial cases
  if (n /= 0) then
    Nseg        = 1
    VLs(Nseg)   = ia1
    Lons(Nseg)  = n
    sys1        = sys1 + n
  else
    Nseg        = 1
    VLS(Nseg)   = ib1
    LonS(Nseg)  = m
    sys2        = sys2 + m
  endif

else

  ! The non trivial cases

  ! We first calcul the length of the segments reconstructing the future segment length
  dn = n / lomaa1      ! The discrete length in direction a1
  rn = n - dn*lomaa1   ! The length of the last segment a1
  dm = m / Lomab1      ! The discrete length in direction b1
  rm = m - dm*Lomab1   ! The length of the last segment b1

  print *, "lomaa =",lomaa1,"  dn =",dn, "  rn =",rn
  print *, "lomab =", Lomab1,"  dm =",dm, "  rm =",rm

  ! On which side do we start building the line
  if (dn > dm ) then

    if (rm /= 0) then
      Nseg       = Nseg + 1
      VLS(Nseg)  = ib1
      LonS(Nseg) = rm
      sys2       = sys2 + rm
    endif

    if (dm == 0) then
      lomaa1  = 0
      dn      = 0
      rn      = n - dn*lomaa1
    else
      lomaa1  = n / dm
      dn      = n / lomaa1
      rn      = n - dn*lomaa1
    endif

    if (rn /= 0) then
      Nseg       = Nseg + 1
      VLS(Nseg)  = ia1
      LonS(Nseg) = rn
      sys1       = sys1 + rn
    endif

    ! The dm remaining segments can now be constructed
    do i = 1, dm
      ! The b1 segments
      Nseg        = Nseg + 1
      VLS(Nseg)   = ib1
      LONS(Nseg)  = Lomab1
      sys2        = sys2 + Lomab1
      ! The a1 segments
      Nseg        = Nseg + 1
      VLS(Nseg)   = ia1
      LONS(Nseg)  = lomaa1
      sys1        = sys1 + lomaa1
    enddo

  elseif (dm > dn) then       ! The opposite case

    if (rn /= 0) then
      Nseg        = Nseg + 1
      VLS(Nseg)   = ia1
      LonS(Nseg)  = rn
      sys1        = sys1 + rn
    endif

    if (dn == 0) then
      lomab1  = 0
      dm      = 0
      rm      = m - dm*lomab1
    else
      lomab1  = m / dn
      dm      = m / lomab1
      rm      = m - dm*lomab1
    endif

    if (rm /= 0) then
      Nseg        = Nseg + 1
      VLS(Nseg)   = ib1
      LonS(Nseg)  = rm
      sys2        = sys2 + rm
    endif

    ! The dn remaining segments can now be constructed
    do i = 1, dn
      ! The a1 segments
      Nseg = Nseg + 1
      VLS(Nseg) = ia1
      LONS(Nseg) = lomaa1
      sys1 = sys1 + lomaa1
      ! The b1 segments
      Nseg = Nseg + 1
      VLS(Nseg) = ib1
      LONS(Nseg) = Lomab1
      sys2 = sys2 + Lomab1
    enddo

  elseif (dm == dn) then          ! The symmetric case

    if (rn /= 0) then
      Nseg        = Nseg + 1
      VLS(Nseg)   = ia1
      LonS(Nseg)  = rn
      sys1        = sys1 + rn
    endif

    if (rm /= 0) then
      Nseg        = Nseg + 1
      VLS(Nseg)   = ib1
      LonS(Nseg)  = rm
      sys2        = sys2 + rm
    endif

    ! The dn remaining segments can now be constructed
    do i = 1, dn
      ! The a1 segments
      Nseg = Nseg + 1
      VLS(Nseg) = ia1
      LONS(Nseg) = lomaa1
      sys1 = sys1 + lomaa1
      ! The b1 segments
      Nseg = Nseg + 1
      VLS(Nseg) = ib1
      LONS(Nseg) = Lomab1
      sys2 = sys2 + Lomab1
    enddo

  endif

endif

! if needed we inter a edge segment in the GD system at the source center
if (includeGD) then

  cum_vec(:)  = 0

LOOPEXT:  do i = 1, Nseg

    ia1     = VLS(i)      ! the veclin of i
    a1(:)   = bveclin(:,ia1)
    lomaa1  = LONS(i)     ! the length of segment i

    do j = 1, lomaa1

      cum_vec(:) = cum_vec(:) + a1(:)
      ProjLongTest = (cum_vec(1)*res(1) + cum_vec(2)*res(2) + cum_vec(3)*res(3))

      if (ProjLongTest > ProjLongRes) then

        if (j == lomaa1) then

          ! The simplest case at the end of the segment
          ! No need to cut the segment length

          ! The following segment index is shifted by one
          VLS((i+2):(Nseg+1)) = VLS((i+1):Nseg)
          LONS((i+2):(Nseg+1)) = LONS((i+1):Nseg)

          VLS(i+1)  = VeclinGD         ! The GD character we want to introduce along the line
          LONS(i+1) = LengthGD        ! The GD segment length we want to introduce along the line

          Nseg = Nseg + 1     ! The list contain one additional segment

          exit LOOPEXT

        else

          ! The general case, the segment must be cut in two pieces

          ! The following segment index is shifted by one
          VLS((i+3):(Nseg+2)) = VLS((i+1):Nseg)
          LONS((i+3):(Nseg+2)) = LONS((i+1):Nseg)

          LONS(i) = j           ! The new length of i

          VLS(i+1)  = VeclinGD     ! The GD character we want to introduce along the line
          LONS(i+1) = LengthGD     ! The GD segment length we want to introduce along the line

          VLS(i+2)  = VLS(i)       ! The character of the second part of i
          LONS(i+2) = lomaa1 - j   ! The length of the second part of i

          Nseg = Nseg + 2   ! The list contain 2 additional segment

          exit LOOPEXT

        endif

      endif

    enddo

  enddo LOOPEXT

endif

print *, " after discreti Nseg =", Nseg

end subroutine decoupsource

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------

subroutine carto

implicit none

real(DP), parameter :: arrondeg = 0.1D0 ! arrondi en degre sur les angles de depart pour eviter
! les probleme des fonction trigo
integer(DPI) :: i,j,sys1,sys2,ivec(3),centre(3),ii,jj,A1(3),B1(3),A2(3),VLS1(1000),lons1(1000)
integer(DPI) :: coin1(3),coin2(3),B2(3),n,m,iA1,iA2,iB1,iB2,Nseg1,Iseg,segparsouI(200),segparsouJ(200)
integer(DPI) :: Nseg2,origine1(3),VLS2(1000),lons2(1000),origine2(3),decoup,jstart,indice,tmp1,tmp2

real(DP) :: Phi1_Start, Phi1_End,Phi2_Start,Phi2_End,taille_source,phi1ini(100),phi2ini(100)
real(DP) :: normeS,NormeAxe,deltal,deltaphi,temp1,Taille_boite,tolerance
real(DP) :: normevis1,normecoin1,normevis2,normecoin2,xxx,PIIar
character:: car1*2,car2*2,car*4,fich*9

logical  :: includeGD

open (60,FILE="../in/info.carto",STATUS='UNKNOWN')
write(60,*)"   I     j    phi1      phi2"
write(60,*)" "

!########################################################
!# Chargement des parametres de control de la simulation
!########################################################
open(10,file="../in/carto_def",STATUS='OLD')


read(10,*) taille_boite  !    taille de la boite de simulation
read(10,*) taille_source !    taille des sources initiales
read(10,*) decoup        !    longueur de dicretisation de reference pour les sources
read(10,*) decalage      !    decalage / centre de la boite
read(10,*) sys1          !    indice du system 1 (reference pour le choix de direction de jonction)
read(10,*) sys2          !    indice du system 2
read(10,*) Phi1_Start    !    borne inferieur de phi1 en degre
read(10,*) Phi1_End      !    borne superieur de phi1 en degre
read(10,*) Phi2_Start    !    borne inferieur de phi2 en degre
read(10,*) Phi2_End      !    borne superieur de phi2 en degre
read(10,*) tolerance     !    tolerance sur les angles en degre
read(10,*) echan1        !    nb d'echantillonnage pour phi1 (mettre un nombre impaire)
read(10,*) echan2        !    nb d'echantillonnage pour phi2 (mettre un nombre impaire)

close(10)

! ici on fait un arrondi de toutes les angles d'entree afin d'evite :/opt/intel/idbe/9.0/bin:/
PIIar = PII + arrondeg*pisurdeg

if (taille_boite < taille_source) &
     stop " The simulated volume linear dimension cannot be smaller than the source length"
if (decoup > 1000) stop " nb of segment discretization too large"
if (decoup < 1) stop " nb of segment discretization too small"
if (echan1 < 1) stop " nb of angle increments on sys 1 too small"
if (echan2 < 1) stop " nb of angle increments on sys 2 too small"
if (echan1 > 100) stop " nb of angle increments on sys 1 too large"
if (echan2 > 100) stop " nb of angle increments on sys 2 too small"

if (Phi1_Start < -180 .or. Phi1_Start > 180 .or.    &
     Phi1_end < -180 .or. Phi1_end > 180  .or.      &
     Phi2_start < -180 .or. Phi2_start > 180  .or.  &
     Phi2_end < -180 .or. Phi2_end > 180 )          &
   stop " The angle domain must be in between [-180 , 180] "

! arrondir les angles d'entree
Phi1_Start = Phi1_Start + arrondeg
Phi1_End = Phi1_End - arrondeg
Phi2_Start = Phi2_Start + arrondeg
Phi2_End = Phi2_End - arrondeg

! calcul de la position du centre de ma boite
ivec(:) = (nint((taille_boite * 0.5D-6 / avalue),DPI)/(2*facteur_boite))*2*facteur_boite
!centre(:) = noeud(ivec,2)
print *, "Facteur boite = ",facteur_boite
centre(:) = (ivec(:)/(ideux* facteur_boite)) * ideux * facteur_boite
Modur(:) = centre(:)*IDEUX

write(*, '("  Dimensions of the volume in AVALUE (",3(I10,1x),")")') MODUR(1:3)
write(*, '("  Coordinates of the volume center   (",3(I10,1x),")")') Centre(1:3)
write(*, '("  Dimensions of the volume in micron (",3(F10.3,1x),")")') MODUR(1:3)*avalue*1D6


NormeS = taille_source * 1.D-6 / avalue

ivis1 = (sys1-1)*nbasered + 1
ivis2 = (sys2-1)*nbasered + 1

! vecteurs vis et coins des deux systems
vis1(:) = bveclin(:,ivis1)
coin1(:)= bveclin(:,ivis1+6)
vis2(:) = bveclin(:,ivis2)
coin2(:)= bveclin(:,ivis2+6)

Normevis1 = norivect(vis1)
Normecoin1= norivect(coin1)
Normevis2 = norivect(vis2)
Normecoin2= norivect(coin2)

! transformation en Radian
Phi1_Start = Phi1_Start * pisurdeg
Phi1_End = Phi1_End * pisurdeg
Phi2_Start = Phi2_Start * pisurdeg
Phi2_End = Phi2_End * pisurdeg

! choix du vecteur de l'axe de jonction
Iaxe1 = 0
Iaxe2 = 0
BOU1:do i = ivis1,ivis1+nbasered-1
    BOU2: do j = ivis2,ivis2+nbasered-1
        if(etat(bveclin(:,i),bveclin(:,j)) > 1) then  ! I et J sont anti-parallels et meme direction
           ! teste de Frank fait un premier trie (on met les etats repulsif au centre de la carto)
           if ((etat(bveclin(:,ivis1),bveclin(:,i))) >=0) then
              if ((etat(bveclin(:,ivis2),bveclin(:,j))) >=0) then
                 if (etat(bveclin(:,ivis1),bveclin(:,ivis2)) > 0) then
                    Iaxe1 = i   ! indice de l'axe de jonction dans le system 1
                    Iaxe2 = j   ! indice de l'axe de jonction dans le system 2
                    axe(:) = bveclin(:,i)
                    print *, " The junction direction =", axe(:)
                    print *, " Iaxe1 =", Iaxe1, "   and Iaxe2 =", Iaxe2
                    exit BOU1
                 endif
              else
                 if (etat(bveclin(:,ivis1),bveclin(:,ivis2+4)) > 0) then
                    Iaxe1 = i   ! indice de l'axe de jonction dans le system 1
                    Iaxe2 = j   ! indice de l'axe de jonction dans le system 2
                    axe(:) = bveclin(:,i)
                    print *, " The junction direction =", axe(:)
                    print *, " Iaxe1 =", Iaxe1, "   and Iaxe2 =", Iaxe2
                    exit BOU1
                 endif
              endif
           else
              if ((etat(bveclin(:,ivis2),bveclin(:,j))) >=0) then
                 if (etat(bveclin(:,ivis1+4),bveclin(:,ivis2)) > 0) then
                    Iaxe1 = i   ! indice de l'axe de jonction dans le system 1
                    Iaxe2 = j   ! indice de l'axe de jonction dans le system 2
                    axe(:) = bveclin(:,i)
                    print *, " The junction direction =", axe(:)
                    print *, " Iaxe1 =", Iaxe1, "   and Iaxe2 =", Iaxe2
                    exit BOU1
                 endif
              else
                 if (etat(bveclin(:,ivis1+4),bveclin(:,ivis2+4)) > 0) then
                    Iaxe1 = i   ! indice de l'axe de jonction dans le system 1
                    Iaxe2 = j   ! indice de l'axe de jonction dans le system 2
                    axe(:) = bveclin(:,i)
                    print *, " The junction direction =", axe(:)
                    print *, " Iaxe1 =", Iaxe1, "   and Iaxe2 =", Iaxe2
                    exit BOU1
                 endif
              endif
           endif

           ! finalement : pour les cas douteux, parmi les deux possibilite on
           ! prend celui qui est dans le quart de plan entre bveclin(ivis1) et bveclin(ivis1+6)
           ! une fois choisit, on selection l'indice Iaxe2 de vecteur parallele dans le system 2
           if (((etat(bveclin(:,ivis1),bveclin(:,i)) >= 0).and.      &
                (etat(bveclin(:,ivis1+6),bveclin(:,i)) >= 0)) .or.   &
                (etat(bveclin(:,ivis1+5),bveclin(:,i)) == 3)       ) then
              Iaxe1 = i   ! indice de l'axe de jonction dans le system 1
              Iaxe2 = j   ! indice de l'axe de jonction dans le system 2
              axe(:) = bveclin(:,i)
              print *, " The junction direction =", axe(:)
              print *, " Iaxe1 =", Iaxe1, "   and Iaxe2 =", Iaxe2
              exit BOU1
           endif
        endif
    enddo BOU2
enddo BOU1

!Verification de l'axe choisit
if(Iaxe1*Iaxe2 == 0 ) stop " The junction direction could not be calculated"

if(etat(bveclin(:,Iaxe1),bveclin(:,Iaxe2)) /= 3) &
   stop " The vectors Iaxe1 and Iaxe2 are not parallel or have different length: erreur"

axe(:) = reduire_vec(axe)
NormeAxe = norivect(axe)

if(etat(axe(:),bvecnor(:,ivis1)) /= 0 ) stop " AXE not included in plan 1"
if(etat(axe(:),bvecnor(:,ivis2)) /= 0 ) stop " AXE not included in plan 2"

!if(etat(axe(:),bveclin(:,ivis1)) < 0 .or. etat(axe(:),bveclin(:,ivis1+6)) < 0) &
!stop " l'axe est mal oriente / vis1 et -coin1"

! l'angle entre l'axe de jonction et les deux directions vis des deux systems
phi1axe = datan2(dble(iprodsca(axe,coin1))/norivect(coin1),dble(iprodsca(axe,vis1))/norivect(vis1))
phi2axe = datan2(dble(iprodsca(axe,coin2))/norivect(coin2),dble(iprodsca(axe,vis2))/norivect(vis2))

print *, "Phi1axe =", phi1axe * degsurpi
print *, "Phi2axe =", phi2axe * degsurpi

! determination de l'angle de chaque future source

temp1 = (Phi1_End - Phi1_start) /dble(echan1) ! interval de variation de phi
print *," Phi1 increment =", temp1*degsurpi

if(echan1 > 1 .and. abs(tolerance) > abs(temp1*degsurpi)) then
   print *, " attention: tolerance > deltatphi in sys 1"
   read *
endif

! pour inclure le point de depart
if (echan1 > 1 ) echan1 = echan1 + 1
Do i = 1,echan1
    xxx =  phi1_start + (i-1) * temp1
    if(xxx > PIIAR) xxx = xxx - 2.0 * PII
    if(xxx < -1.0*PIIAR) xxx = 2.0 * PII + xxx
    phi1ini(I) = xxx          ! angle / axe
    xxx = phi1ini(I) + phi1axe ! angle / VIS1
    if(xxx > PIIAR) xxx = xxx - 2.0 * PII
    if(xxx < -1.0*PIIAR) xxx = 2.0 * PII + xxx
    phi1(I) = xxx
    print *, "source1  :",I
    print *, "Phiini1 = ", phi1ini(I)*degsurpi,"       Phi1 = ", phi1(I)*degsurpi
    IVec(:) = nint(normeS*dcos(phi1(i))/normevis1,DPI)*Vis1(:) + &
         nint(normeS*dsin(phi1(i))/normecoin1,DPI)*Coin1(:)

    xxx = datan2(dble(iprodsca(IVEC,coin1))/norivect(coin1),&
         dble(iprodsca(IVEC,vis1))/norivect(vis1))
    if(xxx > PIIAR) xxx = xxx - 2.0 * PII
    if(xxx < -1.0*PIIAR) xxx = 2.0 * PII + xxx
    if(abs(xxx-phi1(I))*degsurpi>tolerance .and. abs(xxx+phi1(I))*degsurpi>tolerance) then
       print *, " phi1eff = ", xxx*degsurpi
       print *, " vecs1 =", ivec
       print*, " identification of vecs1 failed"
       stop
    else
       VecS1(:,I) = Ivec(:)
       print *, "Phi1eff = ", xxx*degsurpi
    endif
enddo

temp1 = (Phi2_End - Phi2_start) /dble(echan2) ! interval de variation de phi
print *, " Phi1 increment =", temp1*degsurpi
if(echan2 > 1 .and. abs(tolerance) > abs(temp1*degsurpi)) then
   print*, " attention: tolerance > deltatphi in sys 2"
   read *
endif
! pour inclure le point de depart
if(echan2 > 1) echan2 = echan2 + 1
Do i = 1,echan2
    xxx = phi2_start + (i-1) * temp1          ! angle / axe
    if(xxx > PIIAR) xxx = xxx - 2.0 * PII
    if(xxx < -1.0*PIIAR) xxx = 2.0 * PII + xxx
    phi2ini(I) = xxx
    xxx = phi2ini(I) + phi2axe ! angle / VIS2
    if(xxx > PIIAR) xxx = xxx - 2.0 * PII
    if(xxx < -1.0*PIIAR) xxx = 2.0 * PII + xxx
    phi2(I) = xxx
    print *, "source2  :",I
    print *, "Phiini2 = ", phi2ini(I)*degsurpi, "   Phi2    = ", phi2(I)*degsurpi

    Ivec(:) = nint(normeS*dcos(phi2(i))/normevis2)*Vis2(:) + &
         nint(normeS*dsin(phi2(i))/normecoin2)*Coin2(:)
    xxx = datan2(dble(iprodsca(IVEC,coin2))/norivect(coin2),&
         dble(iprodsca(IVEC,vis2))/norivect(vis2))
    if(xxx > PIIAR+0.) xxx = xxx - 2.0 * PII
    if(xxx < -1.0*PIIAR) xxx = 2.0 * PII + xxx
    if(abs(xxx-phi2(I))*degsurpi>tolerance .and. abs(xxx+phi2(I))*degsurpi>tolerance) then
       print*, " Phi2eff = ", xxx *degsurpi
       print*, " vecs2 =", ivec
       print*, " identification of vecs2 failed"
       stop
    else
       VecS2(:,I) = Ivec(:)
       print*, "Phi2eff = ", xxx *degsurpi
    endif
enddo

! on sait mainteneant que chaque config est caracterise par :
!  phi1 : angle entre la ligne de la source 1 et la direction VIS1
!  phi2 : angle entre la ligne de la source 2 et la direction VIS2
!  phi1 et phi2 correspondent a phi1ini et phi2ini : les deux angle par
!  rapport a l'axe de jonction AXE

open (50,FILE="../out/stat.carto",STATUS='UNKNOWN')
if(abs(etat(vis1,vis2)) < 2 ) then
   write(50,'(" Phi1    Phi2Jonc    Phi2Cross   Phi2Rep   SpacingInit   SpacingEnd   %jonc   NTJ   Ljonc     file",/)')
else
   write(50,'(" Phi1    Phi2Annil    Phi2Cross   Phi2Rep    SpacingInit    SpacingEnd   %GD    file",/)')
endif
close(50)
open(30,FILE="../exec/carto",STATUS='UNKNOWN')
write(30,*)"#! /usr/local/bin/bash"
write(30,*)" "
write(30,*)"# script used to proceed with the carto computations "
write(30,*)"# Segments file are one by one copy from ../in/carto/segs.**** "
write(30,*)"# in ../in/segs (the segment file name which must be defined in input.dd "
write(30,*)"# After each copy mM is runned to compute data "
write(30,*)"# ............................................................"


open(40,FILE="../exec/gcarto",STATUS='UNKNOWN')
write(40,*)"#! /usr/local/bin/bash"
write(40,*)" "
write(40,*)"# script used to proceed with the carto computations "
write(40,*)"# Segments file are one by one copy from ../in/carto/segs.**** "
write(40,*)"# in ../in/segs (the segment file name which must be defined in input.dd "
write(40,*)"# After each copy mM is runned to compute data "
write(40,*)"# ............................................................"



do i = 1,echan1
    do j = 1,echan2
        write(60,'(I5,I5,F10.2,F10.2)')i,j,phi1ini(i)*degsurpi,phi2ini(j)*degsurpi
    enddo
enddo
write(60,'(///)')
write(60,'("-------------------------------------------------------------------")')


! boucle de generation des sources 1
ISEG = Izero
Do i = 1, echan1
    write(*,*)  "=========================================      Source 1 num =",i
    write(60,*) "=========================================      Source 1 num =",i
    ivec (:) = VecS1(:,I)
    !    write(*,'(" Vecteur source = (",I7,x,I7,x,I7,")")') ivec(1:3)
    write(60,'(" Vecteur source = (",I7,x,I7,x,I7,")")') ivec(1:3)
    A1(:) = (/0,0,0/)
    B1(:) = (/0,0,0/)
    call selecAB(I,IUN,iA1,iB1)
    a1 = bveclin(:,ia1)
    b1 = bveclin(:,ib1)
    ! verifcation que les vecteurs de discretisation font angle aigu
    if(etat(a1,b1) <= 0) then
       print *,"sourceI =",I, "; iA1,iB1", ia1,ib1
       print*, " Bad selection of a1 and b1"
       stop
    endif

    ! verifcation que la source est entre les deux vecteurs
    if(etat(ivec,a1) <= 0 .or. etat(ivec,b1) <= 0 ) then
       print *,"sourceI =",I, "; iA1,iB1", ia1,ib1
       print*, " Bad selection of a1 and b1"
       stop
    endif

    if(.not. egalite(dabs(dacos(icose2v(a1,b1))),                         &
         dabs(dacos(icose2v(IVec,a1))) + dabs(dacos(icose2v(IVec,b1)))))  &
       stop " VecS1 is not in between a1 and b1 "

    call decompo (VecS1(:,I),A1,B1,n,m,deltal,deltaphi)
    deltaphi = deltaphi * degsurpi
    write(*,'(" iA =",I3,"   nA =", I6,"  iB =",I3,"   nB =",I6)') ia1,n,ib1,m
    write(60,'(" iA =",I3,"   nA =", I6,"  iB =",I3,"   nB =",I6)') ia1,n,ib1,m
    write(*,'(" Delta(L) =", F5.2 ,"%   ---   delta(phi) = ", F5.2," Degre")') &
         deltal*100.0, deltaphi
    write(60,'(" Delta(L) =", F5.2 ,"%   ---   delta(phi) = ", F5.2," Degre")') &
         deltal*100.0, deltaphi

    if(n < 0) stop " n < 0 for 1"
    if(m < 0) stop " m < 0 for 1"

    if(deltaphi > tolerance) then
       print *, " Initial source =", vecS1(:,I)
       print *, " Calculated source ", n * a1(:) +  m * b1(:)
       print *, " Sorry, the angle increment could not be achieved"
       print *, " To solve that problem, please increase the angle increment or descrease "
       print *, " the simulation scale (ECHELLE)."
       stop
    endif

    includeGD = .false.   ! IncludeGD must be .true. to include a GD segment at the center of the source in subroutine decoupesource

    call decoupsource (iA1,iB1,n,m,decoup,VLS1,LONS1,Nseg1,includeGD)

    write(60,'(" Nb of seg building the source ", I20)') Nseg1

    segparsouI(i) = Nseg1
    ivec(:) = 0
    temp1 = 0

    do ii = 1, Nseg1
        ivec(:) = Ivec(:) + bveclin(:,VLS1(II))*LONS1(II)
    enddo

    !    print *, "S discretise =", ivec
    if (inorivect(ivec(:) - n*A1(:) - m*B1(:)) /= 0 ) stop "Bad discretization of S1"

    ! calcule de lorigine de la source 1
    decalage = decalage / IDEUX * IDEUX   ! parce que le pas minimal de decalage est 3*vis1(:) + 2*coin1(:)
    origine1(:) = centre(:) - ivec(:)/2 + decalage*Bvecdep(:,Iaxe1)

    print *, "dec1=", origine1(:)-centre(:)

    if (iprodsca(origine1-centre,bvecnor(:,ivis1)) /= 0) stop " Shift out of plan 1"

    ii = 1

    if (Nseg1 == 1) then

      ISEG              = ISEG + IUN
      seg(iseg)%O(:)    = origine1(1:3)
      seg(iseg)%norme   = Lons1(ii)
      seg(iseg)%veclin  = VLS1(II)
      seg(iseg)%voiso   = izero
      seg(iseg)%voise   = izero
      write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin,           &
                                      seg(iseg)%voiso,seg(iseg)%voise,iseg

    else

      ! The first segment
      ISEG              = ISEG + IUN
      seg(iseg)%O(:)    = origine1(1:3)
      seg(iseg)%norme   = Lons1(ii)
      seg(iseg)%veclin  = VLS1(II)
      seg(iseg)%voiso   = izero
      seg(iseg)%voise   = IUN
      ivec(:)           = origine1(1:3) + Lons1(ii)* bveclin(:,VLS1(II))
      write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin,           &
                                      seg(iseg)%voiso,seg(iseg)%voise,iseg

      ! The intermediate segments
      if(Nseg1 > 2) then

        do ii = 2 , Nseg1-1
          ISEG              = ISEG + IUN
          seg(iseg)%O(:)    = Ivec(1:3)
          seg(iseg)%norme   = Lons1(ii)
          seg(iseg)%veclin  = VLS1(II)
          seg(iseg)%voiso   = -1
          seg(iseg)%voise   = IUN
          ivec(:)           = ivec(:) + Lons1(ii)* bveclin(:,VLS1(II))
          write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin,       &
                                          seg(iseg)%voiso,seg(iseg)%voise,iseg
        enddo

      endif

      ! The last segment
      ii = Nseg1
      ISEG              = ISEG + IUN
      seg(iseg)%O(:)    = Ivec(1:3)
      seg(iseg)%norme   = Lons1(ii)
      seg(iseg)%veclin  = VLS1(II)
      seg(iseg)%voiso   = -1
      seg(iseg)%voise   = IZERO
      write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin,           &
                                      seg(iseg)%voiso,seg(iseg)%voise,iseg

    endif

enddo

! Beginning of the second segment building

Jstart = Iseg + IUN
Do j = 1, echan2
    write(*,*)  "*********************     System 2 num =",J
    write(60,*) "*********************     System 2 num =",J
    ivec(:) = vecS2(:,j)
    !    write(*,'(" Vecteur source = (",I7,x,I7,x,I7,")")') ivec(1:3)
    write(60,'(" Vecteur source = (",I7,x,I7,x,I7,")")') ivec(1:3)
    A2(:) = (/0,0,0/)
    B2(:) = (/0,0,0/)
    call selecAB(J,IDEUX,iA2,iB2)
    a2 = bveclin(:,ia2)
    b2 = bveclin(:,ib2)
    if(etat(a2,b2) <= 0) stop "Bad selection of a2 and b2"
    ! verifcation que la source est entre les deux vecteurs
    if(etat(ivec,a2) <= 0 .or. etat(ivec,b2) <= 0 ) then
       print *,"sourceJ =",J, "; iA1,iB1", ia2,ib2
       print*, "Bad selection of a2 and b2"
       stop
    endif
    if(.not. egalite(dabs(dacos(icose2v(a2,b2))),                         &
         dabs(dacos(icose2v(IVec,a2))) + dabs(dacos(icose2v(IVec,b2)))))  &
       stop " VecS2 is not in between : a2 et b2 "

    call decompo (VecS2(:,J),A2,B2,n,m,deltal,deltaphi)
    deltaphi = deltaphi * degsurpi
    write(*,'(" iA =",I3,"   nA =", I6,"  iB =",I3,"   nB =",I6)') ia2,n,ib2,m
    write(60,'(" iA =",I3,"   nA =", I6,"  iB =",I3,"   nB =",I6)') ia2,n,ib2,m
    write(*,'(" Delta(L) =", F5.2 ,"%   ---   delta(phi) = ", F5.2," Degree")') &
         deltal*100.0, deltaphi
    write(60,'(" Delta(L) =", F5.2 ,"%   ---   delta(phi) = ", F5.2," Degree")') &
         deltal*100.0, deltaphi
    !   print *," Le nom de fichier correspondant : segs.",car
    !   print *, "S decompo =", n*A2(:) + m*B2(:)

    if(n < 0) stop " n < 0 pour 2"
    if(m < 0) stop " m < 0 pour 2"

    if(deltaphi > tolerance) then
       print *, " Initial source =", vecS2(:,j)
       print *, " Calculated source ", n*a2(:) +  m * b2(:)
       print *, " Sorry, the angle increment could not be achieved"
       print *, " To solve that problem, please increase the angle increment or descrease "
       print *, " the simulation scale (ECHELLE)."
       stop
    endif

    includeGD = .true.   ! IncludeGD must be .true. to include a GD segment at the center of the source in subroutine decoupesource

    call decoupsource (iA2,iB2,n,m,decoup,VLS2,LONS2,Nseg2,includeGD)
    write(60,'(" Number of segments in the source ", I20)') Nseg2

    segparsouJ(J) = Nseg2
    ivec(:) = 0
    do jj = 1, Nseg2
        ivec(:) = Ivec(:) + bveclin(:,VLS2(JJ))*LONS2(JJ)
    enddo
    !        print *, "S2 discretise =", ivec

    if (.not. includeGD) then
      if(inorivect(ivec(:) - n*A2(:) - m*B2(:)) /= 0 ) stop "Bad discretization of S2"
    endif

    ! calcul de l origin de la source 2

    if (includeGD) then

      origine2(:) = centre(:)-ivec(:)/2 + decalage * Bvecdep(:,Iaxe2)

    else

      origine2(:) = centre(:)-ivec(:)/2 + decalage * Bvecdep(:,Iaxe2)

      if(iprodsca(origine2-centre,bvecnor(:,ivis2)) /= 0) stop " Outplane shift for plane 2"

    endif

    ii = 1
    if(Nseg2 == 1) then
       ISEG = ISEG + IUN
       seg(iseg)%O(:) = origine2(1:3) ; seg(iseg)%norme = Lons2(ii)
       seg(iseg)%veclin = VLS2(II) ; seg(iseg)%voiso = izero ; seg(iseg)%voise = izero
       write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin, &
            seg(iseg)%voiso,seg(iseg)%voise,iseg
    else
       ISEG = ISEG + IUN
       seg(iseg)%O(:) = origine2(1:3) ; seg(iseg)%norme = Lons2(ii)
       seg(iseg)%veclin = VLS2(II) ; seg(iseg)%voiso = izero ; seg(iseg)%voise = IUN

       ivec(:) = origine2(1:3) + Lons2(ii)* bveclin(:,VLS2(II))
       write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin, &
            seg(iseg)%voiso,seg(iseg)%voise,iseg
       if(Nseg2 > 2) then
          do ii = 2 , Nseg2-1
              ISEG = ISEG + IUN
              seg(iseg)%O(:) = Ivec(1:3) ; seg(iseg)%norme = Lons2(ii)
              seg(iseg)%veclin = VLS2(II) ; seg(iseg)%voiso = -1 ; seg(iseg)%voise = IUN
              ivec(:) = ivec(:) + Lons2(ii)* bveclin(:,VLS2(II))
              write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin, &
                   seg(iseg)%voiso,seg(iseg)%voise,iseg
          enddo
       endif
       ii = Nseg2
       ISEG = ISEG + IUN
       seg(iseg)%O(:) = Ivec(1:3) ; seg(iseg)%norme = Lons2(ii)
       seg(iseg)%veclin = VLS2(II) ; seg(iseg)%voiso = -1 ; seg(iseg)%voise = IZERO
       write(60,'(3I8,I6,I4,2I8,I5)') seg(iseg)%O(:),seg(iseg)%norme,seg(iseg)%veclin, &
            seg(iseg)%voiso,seg(iseg)%voise,iseg
    endif
enddo

do ii = 1, nseg1
enddo

! boucle double d'ecriture
ISEG = IZERO
ii = iUN

4 format(I5,3I10,I5,I7,4I7,L3,I7,L3,I8)

Do i = 1, echan1
    Nseg1 = segparsouI(I)
    jj = Jstart
    Do j = 1, echan2
        car1 = char(i/10 + 48)//char(i- (i/10)*10 + 48)
        car2 = char(j/10 + 48)//char(j- (j/10)*10 + 48)
        car = car1//car2
        write(30,*)" cp ../in/carto/segs.",car," ../in/SEGS"
        write(40,*)" cp ../in/carto/segs.",car," ../in/SEGS"
        write(30,*)" ../bin/mm"
        write(40,*)" ../bin/gmm"
        fich = "segs."//car
        indice = izero
        open(20,FILE="../in/carto/"//fich,STATUS='UNKNOWN')
        write(20,*)" 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
        Nseg2 = segparsouJ(J)
        write(20,*) Nseg1 + Nseg2
        write(20,*) modur

        do iseg = ii , ii + Nseg1-1
            indice = indice + IUN
            tmp1 = seg(iseg)%voiso  ;          tmp2 = seg(iseg)%voise
            seg(iseg)%O(1:3) = (seg(iseg)%O(1:3)/(facteur_boite))*facteur_boite
            if(tmp2 == 1) tmp2 = -1
!            write(20,4) seg(iseg)%O(1:3),seg(iseg)%norme,seg(iseg)%veclin,seg(iseg)%voiso,seg(iseg)%voise,iseg
            write (20,4) indice,seg(iseg)%O(1:3),seg(iseg)%norme,seg(iseg)%veclin,tmp1,tmp1,tmp2,tmp2,&
                 seg(Iseg)%JONC, izero,seg(iseg)%gd,iseg
        enddo

        do iseg = jj , jj + Nseg2-1
            indice = indice + IUN
            tmp1 = seg(iseg)%voiso  ;          tmp2 = seg(iseg)%voise
            seg(iseg)%O(1:3) = (seg(iseg)%O(1:3)/(facteur_boite))*facteur_boite
            if(tmp2 == 1) tmp2 = -1
!            write(20,4) seg(iseg)%O(1:3),seg(iseg)%norme,seg(iseg)%veclin,seg(iseg)%voiso,seg(iseg)%voise,iseg
            write (20,4) indice,seg(iseg)%O(1:3),seg(iseg)%norme,seg(iseg)%veclin,tmp1,tmp1,tmp2,tmp2,&
                 seg(I)%JONC, izero,seg(i)%gd,iseg
        enddo
        jj = jj + Nseg2
        Write(20,'(//,I10,I10,"   Iaxe1 and Iaxe2")')Iaxe1,Iaxe2
        Write(20,'(/,F8.2,5x,F8.2,"   phi1 and phi2",/)') phi1ini(I)*degsurpi,phi2ini(J)*degsurpi
        write(20,*) fich,"          Name of the input segment file"
        write(20,*) " "
        write(20,*) relax_TL, "         Number of line tension relaxation steps (relax_TL)"
        write(20,*) relax_INT,"         Number of relaxation steps before loading"
        write(20,'(////)')

        Write(20,*)" List of parameter used to define this intial configuration."
        Write(20,*)"   "
        Write(20,*)"   "
        Write(20,*)" Crystal type  =  ",crystal_structure(1:3)
        write(20,'(" Dimension of the simulated volume = ", F8.5," micron" )') taille_boite
        write(20,'(/," Simulation unit (avalue) = ", F8.5," A" )') avalue * 1D10
        write(20,'(/," Source segments must be made with less than ", I8, " segments" )') decoup
        write(20,'(/," Angle increment for SOURCE 1 / Junction axis = ", F8.2, " degres" )') &
             phi1ini(I) * degsurpi
        write(20,'(/," Angle increment for SOURCE 2 / Junction axis = ", F8.2, " degres" )') &
             phi2ini(j) * degsurpi

        write(20,'(/," Vecteur source = (",I7,x,I7,x,I7,")")') vecs1(1:3,I)
        write(20,'(/," iA1 =",I3,"  iB1 =",I3)') ia1,ib1
        !        write(*,'(" Vecteur source = (",I7,x,I7,x,I7,")")') VecS2(1:3,J)
        write(20,'(/," iA2 =",I3,"  iB2 =",I3)') ia2,ib2
        close(20)
    enddo
    ii = ii + Nseg1
enddo
print *, Jstart
close(30)
close(40)
close(60)

end subroutine carto

!*****************************************************
!*****************************************************

subroutine double_config
Integer(DPI) :: dislo,vl, O(3),s,vec(3),i

PLUSSEG = 0

Do dislo = 1, nsegm
    if (seg(dislo)%voise == 0) then
       vl = seg(dislo)%veclin
       O(:) = seg(dislo)%O(:)+ seg(dislo)%norme*bveclin(:,vl)

       ! composante de decalage parallele a vis1
       Vec(:) = Bvecdep(:,assoc(vl,3))*(Nint(Ldedou/avalue/normdep(assoc(vl,3)))/4)*4
       ! on ajoute une composante de decalage parallele a coin1 ou coin2

       O(:) = O(:) + Vec(:)
       s = dislo
       do while (S /= 0)
           vl = seg(s)%veclin
           PLUSSEG = PLUSSEG + 1
           i = NSEGM + PLUSSEG
           !           print *,"s =", s,"      i=", i
           seg(i)%O = O(:)
           seg(i)%norme = seg(s)%norme
           if(seg(s)%voise == 0) then
              seg(i)%voiso = 0
           else
              seg(i)%voiso = i-1
           endif
           if(seg(s)%voiso == 0) then
              seg(i)%voise = 0
           else
              seg(i)%voise = i+1
           endif
           if(vl+4 <= nbase .and. syseg(vl) == syseg(vl+4)) then
              seg(i)%veclin = vl + 4
           else
              seg(i)%veclin = vl - 4
           endif
           s = seg(s)%voiso
           O(:) = seg(i)%O(:)+ seg(i)%norme*bveclin(:,seg(i)%veclin)
       enddo
    endif
enddo

print *, " Initial number of segments =", nsegm
print *, " Number of segments added =", PLUSSEG
NSEGM = NSEGM + PLUSSEG

PLUSSEG = 0

end subroutine double_config

end module cartographie
