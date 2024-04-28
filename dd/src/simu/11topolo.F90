
!===================================================================================================
!========================    DEBUT    MODULE   "TOPOLO"  ===========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to topological problems taking places in the computations.
!> This include boundary condition limits.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module TOPOLO

use VARGLOB
use DEBUG
use CONNEC
use microstructure

implicit none

logical   :: Recouv(nsegmax),mathieugdo,mathieugde,condition
logical, allocatable, save        :: Key_Box_Border(:)  !< key is true when a box is on a border of the simulated volume
integer(DPI)                      :: NBoxBorder         !< total number of box on a border
integer(DPI), allocatable, save   :: BoxBorder(:)       !< tab of box that are at the border of the simulation volume
contains

!#################################################################################
! Definition of the Greengard domains in the simulated volume
!#################################################################################
subroutine Allocation_initiale

implicit none

integer(kind=DPI)     :: IX, IY, IZ
integer(kind=DPI)     :: I, J, K
integer(kind=DPI)     :: BX, BY, BZ
integer(kind=DPI)     :: tmp, nparmax_boite, compteur
integer(kind=DPI)     :: X(3)
logical               :: EmptyBoite

#ifdef PA
integer(kind=DPI)               :: ii,jj,IB,JB,KB
real(kind=DP)                   :: Boxcenter(3), highest_length
#endif

! if allocation_dynamique_boites = True, we try to calculate the optimum number of domains
if(allocation_dynamique_boites) then

  ! Equation defining the optimum number of domains with homogeneous dislocation distribution
  Nboites = nint(DSQRT(real(KRC*27*NSEGM,DP)))

  ! Special case corrections
  if(Nboites > nsegm) Nboites  = nsegm
  if(Nboites < 64) then
    L_boite = maxval(modur(:))              ! Number of domains smaller than 4x4x4 is considered as useless
  else
    L_boite = (Volume / NBoites)**(TIERS)   ! The optimum solution assuming a volume with homogenous shape
  endif

endif

! If the simulation box contains particles, domains cannot be smaller than the diameter of particles
if (particules) then
  if (L_boite < maxval(par(1:npar)%R)) then
    write(*,*) "--------------------------------------------------------------"
    write(*,'( "A particle is too large for optimum multi-pole domains!" )')
    write(*,*) "You must run the simulation by imposing a constant L_boite = ",  &
                L_boite,npar,(1.2*maxval(real(par(1:npar)%R,DP)))/avalue
    write(*,*) "--------------------------------------------------------------"
    stop
  endif
endif

! The minimum size of the Greengard domains must be at least the discretization length
if (L_boite < deux * Xlomax) then
  write(*,*) "--------------------------------------------------------------"
  write(*,'( "Discretization length is too large for optimum multi-pole domains!" )')
  write(*,*) "Xlomax =", Xlomax, "  optimum L_domains =" , l_boite
  L_boite =  deux * Xlomax
  write(*,*) "Multipole L_domains is increased to : ", L_boite, '==>', (L_boite * avalue*1.e6), 'micron'
  write(*,*) "--------------------------------------------------------------"
endif

! if the discretization length is to large we stop the Greengard algo
if ( .not. methode_boites .or.  (trois * L_boite) > real(minval(modur(:)))) then

  write (*,*) " Multipole domains cannot be used", nint(deux * L_boite),minval(modur(:))
  NBoites  = IUN
  NBoitesX = IUN
  NBoitesY = IUN
  NBoitesZ = IUN

  nsegmax_boite = nsegmax
  nparmax_boite = npar

else
  ! We look now for the optimum domain dimensions as function of the simulated volume shape

  ! as a starting point we evaluate the number of domains in each directions
  NBoitesX = int(Modur(1)/L_boite)
  NBoitesY = int(Modur(2)/L_boite)
  NBoitesZ = int(Modur(3)/L_boite)

  ! The maximum number of segments and particles we can have in a domain for the Greengard algo
  nsegINI = nsegm
  NBoites = NBoitesX*NBoitesY*NBoitesZ
  if (GB /= 0) then
    nsegmax_boite = (nsegmax / IDIX)
  else
    nsegmax_boite = (nsegmax / NBoites + 1) * ICENT
  endif
  if (nsegmax_boite < ListeSegRangeIni) nsegmax_boite = ListeSegRangeIni + icent
  nparmax_boite = (npar / NBoites + 1) * IDIX
  if(nparmax_boite < 10)  nparmax_boite = IDIX

  ! The particular case when the number of domains in one direction is less than 4
  ! Then the following coding does not work and the multipole algo must not be applied
  if (NBoitesX < 4 .or. NBoitesY < 4 .or. NBoitesZ < 4) then
    write (*,*) " One dim of the simulated volume is too small for Multipoles", NBoitesX, NBoitesY, NBoitesZ
    NBoites  = IUN
    NBoitesX = IUN
    NBoitesY = IUN
    NBoitesZ = IUN

    nsegmax_boite = nsegmax
    nparmax_boite = npar
  endif

endif

! Calculation of the final size of the boxes
tailleBoite(1) = real(Modur(1),DP) / NBoitesX
tailleBoite(2) = real(Modur(2),DP) / NBoitesY
tailleBoite(3) = real(Modur(3),DP) / NBoitesZ

write(*,'("Setting  : Initial number of Greengard domains",I6," (",I3,",",I3,",",I3,") &
          &  Size(",ES10.2E2,",",ES10.2E2,",",ES10.2E2,")")')                          &
     NBoites,NBoitesX,NBoitesY,NBoitesZ,TailleBoite(1:3)

if(Xlomax > int(minval(TailleBoite(:)))) then
  write(*,'(" After domains calculation : discretization length > 0.5 volume dimensions")')
  stop "STOP !!!"
endif

! Ultimate check
if((real(Modur(1),DP) - TailleBoite(1) * NBoitesX) > numtols .or. &
   (real(Modur(2),DP) - TailleBoite(2) * NBoitesY) > numtols .or. &
   (real(Modur(3),DP) - TailleBoite(3) * NBoitesZ) > numtols)     &
  stop " A problem exist in the greengard domain calculation !"

! allocation of domains for the multipoles algo
allocate (NsegBoite(NBoites))
allocate (IndexBoite(NBoites))
allocate (B3D_1D(NBoitesX,NBOITESY,NBoitesZ))
allocate (B1D_3D(NBoites,3))
allocate (SIGBOX(3,3,NBoites)) ! contrainte aucentre des boites
allocate (IndBoiteVois(NBoites, 27))

if (nboitesX < 4 .or. nboitesY < 4 .or. nboitesZ < 4) then
  nboxborder = 1
else
  nboxborder = 2*((NboitesX-2)*(NboitesY-2)+(NboitesX-2)*(NboitesZ-2)+(NboitesY-2)*(NboitesZ-2))+   &
               4*(NboitesX+NboitesY+NboitesZ-6) + 16 +                                              &
               2*((NboitesX-4)*(NboitesY-4)+(NboitesX-4)*(NboitesZ-4)+(NboitesY-4)*(NboitesZ-4))+   &
               4*(NboitesX+NboitesY+NboitesZ-12)
endif

allocate (boxborder(nboxborder))
allocate (NBox2Shift(NBoites))
allocate (Box2Shift(NBoites,26))
allocate (Key_Box_2Border(NBoites))
allocate (Key_Box_Border(NBoites))

! Gammabox allocation is placed here because the knowledge of NBoites is needed
if ((.not. allocation_dynamique_boites) .and. kk==0 .and. calculate_gammabox) then
  allocate (Gammabox(Nboites,NTSG))   ! In order to calculate the gamma on each box
  Gammabox(:,:) = 0 !Initialization
endif

#ifdef PA
! The boxes close to the boundary
allocate(list_boxes(NBoites))
#endif

! tabulation intermediaire utile pour l'utilisation de la methode de boites
if(NBoites /= IUN) then

  RCMAKE = .true.

  Do IX = 1,NBoitesX
    Do IY = 1,NBoitesY
      Do IZ = 1,NBoitesZ

        ! cette formule represente une relation entre Iboite et les trois indices des boites
        ! c'est utile surtout dans contact parce qu'il faut deplacer les segment par ordre predefini
        J = IX + (IY-1) * NBoitesX + (IZ-1)* (NBoitesX*NBoitesY)
        B3D_1D(IX,IY,IZ) = J
        B1D_3D(J,1) = IX
        B1D_3D(J,2) = IY
        B1D_3D(J,3) = IZ

        ! calcul du tableau logique de voisinage de sous_boites
        ! chaque boite doit avoir 26 boites + elle-meme comme voisines
        ! Ix,IY,IZ sont les indice de chaque boite
            compteur = iun
            IndBoiteVois(B3D_1D(IX,IY,IZ),compteur) = B3D_1D(IX,IY,IZ) ! premiere boite consideree voisine de BI est BI elle-meme

            do I = -IUN , IUN
              do J = -IUN , IUN
                do K =  -IUN , IUN

                  ! desole : les CPL shift ou tourne ne sont pas pris en compte
                  BX = modulo((IX + I - IUN),NBoitesX) + IUN
                  BY = modulo((IY + J - IUN),NBoitesY) + IUN
                  BZ = modulo((IZ + K - IUN),NBoitesZ) + IUN
                  ! print *, "-->", IX,IY,IZ,I,J,K,BX,BY,BZ
                  tmp = BX + (BY-1) * NBoitesX + (BZ-1)* (NBoitesX*NBoitesY)

                  if (I /= izero .or. J/= izero .or. K /= izero) then
                    compteur = compteur + iun
                    IndBoiteVois(B3D_1D(IX,IY,IZ),compteur) = tmp
                  endif

                  if (compteur > 27 ) stop 'erreur dans l indexation des boites voisines'

                enddo
              enddo
            enddo

      enddo
    enddo
  enddo

else

  RCMAKE = .false.   ! clef pour forcer le calcul de sigma INT CP
  B3D_1D(:,:,:) = IUN
  B1D_3D(IUN,:) = IUN
  IndBoiteVois(:,:) = iun

endif

!sq In what follows, we define the box which are on a border of the simulated volume
! this will be usefull for the obstacle detection. Because of the use of PBC, segments
! that are near a border of the simulation volume can be inside the volume itself or
! in its replica. This is why a shift procedure is applied defore contact detection to
! define the smallest distance between a segment i and a potential obstacle j.
!initialization
NBoxBorder = izero !total number of boxes that are at the border of the simulation volume
BoxBorder(:) = izero !tab of boxes that are at the border of the simulation volume
NBox2Shift(:) = izero ! total number of boxes to shift for a given box BI at the border of the simulation box
Box2Shift(:,26) =izero ! no. of boxes to shift for a given box BI at the border of the simulation box
Key_Box_Border(:) = .false. !key is true when a box is on a border of the simulated volume
Key_Box_2Border(:) = .false. !key is true when a box is on a border or at the second position near a border
                             ! of the simulated volume

! if (nboites > 64) then ! if there are less than 5 boxes per sides of the simulation volume, this procedure is useless
if (NBoitesX > 3 .and. NBoitesY > 3 .and. NBoitesZ > 3) then
  do J = 1, NBoites
    IX = B1D_3D(J,1)
    IY = B1D_3D(J,2)
    IZ = B1D_3D(J,3) ! no. of the considered box

    !#1 In a first step, we make a tab of the boxes that are near the border
    if (IX==1 .or. IX==nboitesX .or. iY==1 .or. IY==nboitesY .or. IZ==1 .or. IZ==nboitesZ) then
      Key_Box_Border(J)   = .true.
      Key_Box_2Border(J)  = .true.
      NBoxBorder = NBoxBorder + 1
      BoxBorder(NBoxBorder) = J
      ! print *, '#',NBoxBorder , J, ix, iy, iz
    elseif (IX==2 .or. IX==nboitesX-1 .or. iY==2 .or. IY==nboitesY-1 .or. IZ==2 .or. IZ==nboitesZ-1) then
      Key_Box_2Border(J) = .true.
      NBoxBorder = NBoxBorder + 1
      BoxBorder(NBoxBorder) = J
      ! print *, '#',NBoxBorder , J, ix, iy, iz
    endif
  enddo

  ! A test to ensure that the number of box near the borders of the simulated volume is correct
  compteur = 2*((NboitesX-2)*(NboitesY-2)+(NboitesX-2)*(NboitesZ-2)+&
             (NboitesY-2)*(NboitesZ-2))+ 4*(NboitesX+NboitesY+NboitesZ-6) + 16 + &
             2*((NboitesX-4)*(NboitesY-4)+(NboitesX-4)*(NboitesZ-4)+ &
             (NboitesY-4)*(NboitesZ-4))+ 4*(NboitesX+NboitesY+NboitesZ-12)
  if (NBoxBorder /= compteur) then
    print *, 'Oups! There is a problem with the initial number of Greengard domains!', NBoxBorder, compteur
    stop
  endif

  !#2 In a second step, for all the boxes that are near a border, we define the tab
  ! of boxes that have to be shifted
  do I = 1, NBoxBorder
    J = BoxBorder(I) ! N° of the box near a border
    IX = B1D_3D(J,1)
    IY = B1D_3D(J,2)
    IZ = B1D_3D(J,3)

    !initialization
    NBox2Shift(J)= 0
    Box2Shift(J, 1:26)= 0

    ! among all the neighboors of J, we retain the boxes that are AT the
    ! border of the simulted volume (J itself also)
    do compteur= 1, 27
      k = IndBoiteVois(J,compteur) ! K a neighboor of J
      if (Key_Box_Border(K)) then
        NBox2Shift(J) = NBox2Shift(J) + iun
        Box2Shift(J,NBox2Shift(J)) = K
      endif
    enddo

    if (NBox2Shift(J) < 7 .or.NBox2Shift(J) > 26) then
      print *, 'PB, nbr de Box2shift > 26' , NBox2Shift(J)
      stop
    endif

  enddo

else ! there are less than 4 boxes per sides of the simulation volume

  ! we gain nothing with this during shift procedure
  Key_Box_2Border(:) = .false.

endif

! The IndexBoite(i)%ListeSeg can now be allocated
! domains (boites) outside of the effective dislocated volume are set to dim=1 to minimize the total tab size
Do IX = 1,NBoitesX
  Do IY = 1,NBoitesY
    Do IZ = 1,NBoitesZ

      EmptyBoite = .false.
      I = B3D_1D(IX,IY,IZ)

      if (GB == 1) then

        ! Frontiers between boites and the effective simulation volume are tested
        if ( IX * tailleBoite(1) < dom_effective_deb(1)) EmptyBoite = .true.
        if ( IY * tailleBoite(2) < dom_effective_deb(2)) EmptyBoite = .true.
        if ( IZ * tailleBoite(3) < dom_effective_deb(3)) EmptyBoite = .true.
        if ( (IX-1) * tailleBoite(1) > dom_effective_end(1)) EmptyBoite = .true.
        if ( (IY-1) * tailleBoite(2) > dom_effective_end(2)) EmptyBoite = .true.
        if ( (IZ-1) * tailleBoite(3) > dom_effective_end(3)) EmptyBoite = .true.

      endif

      if (EmptyBoite) then
        allocate(IndexBoite(I)%ListeSeg(1))
      else
        allocate(IndexBoite(I)%ListeSeg(NSEGMAX_boite))
      endif

    enddo
  enddo
enddo

! initialization
SIGBOX(:,:,:) = izero
NsegBoite(:) = IZERO

do i = 1,NBoites
  IndexBoite(i)%ListeSeg(:) = izero
enddo

! allocation des tableaux pour la methode de boites appliquee aux particules
if (particules) then

  ! calcul du nombre moyens de particules par boite

  allocate (NparBoite(NBoites))
  allocate (listeparBoite(NBoites,nparmax_boite))
  allocate (Npartester(Nboites))
  allocate (listepartester(Nboites,10*nparmax_boite))

  NparBoite(:) = IZERO
  ListeparBoite(:,:) = izero
  IBOITEP(1:Npar) = IZERO

  do I =  1, Npar

      X(:) =  par(i)%C(:)
      if (Nboites /= iun) then
         ! le modulo permet de respecter les CPL pour les boite
         Bx = modulo(int(X(1) / tailleBoite(1),DPI),NboitesX) + IUN
         By = modulo(int(X(2) / tailleBoite(2),DPI),NboitesY) + IUN
         Bz = modulo(int(X(3) / tailleBoite(3),DPI),NboitesZ) + IUN
         J = B3D_1D(Bx,By,Bz)

         if(NparBoite(J) > nparmax_boite) then
            print *," nparmax_boite =",nparmax_boite, " npar  =",npar
            print *,"J,IB,JB,KB,",J,Bx,By,Bz
            print *,"NsegBoite(IB,JB,KB)",NparBoite(J)
            print *,"erreur REGION : debordement de tableau listeboite kk=", kk
            stop
         endif

         IBOITEP(I) = J
         NparBoite(J) =  NparBoite(J) + IUN
         Listeparboite(J,NparBoite(J)) = I
      else
         ! WHEN NO MULTIDOMAINE applied, all segments are tested every where
         NparBoite(IUN) = npar
         listeparBoite(IUN,i) = i
         IBOITEP(I) = IUN
      ENDIF

  enddo

  J= izero
  !==================================================================
  ! verification de la procedure
  do i = 1, Nboites
      J = J + Nparboite(I)
  enddo
  if (npar /=  J) stop  " probleme de methode de  boites avec les par"
  !==================================================================
  Npartester(:) = izero
  Listepartester(:,:) = izero
  if(nboites > iun) then
     do I = 1, NBOITES
         do J = 1, NBOITES
             ! et on stocke les particules visites dans le tableau listeparteste.
             if(.not. any(IndBoiteVois(I,:)-J .eq. IZERO)) CYCLE
             if(NparBoite(J) == Izero) CYCLE
             do k = 1, NparBoite(J)
                 tmp = listeparboite(J,k)
                 Npartester(I) = npartester(I) + IUN
                 listepartester(I,npartester(I)) = tmp
             enddo
         enddo
         if(Npartester(I) > 10*nparmax_boite)then
            print *, " Erreur: topolo: debordement de tableau listeparboite:",I
            stop
         endif
     enddo
  elseif(NBoites == iun) then
     Npartester = Npar
     !**** nouveau CHANGEMENT DE REPERE  et + rapide
     do tmp = 1 , Npar
         listepartester(IUN,tmp) = tmp
    enddo
  else
     stop " Undefined Nboite in TOPOLO"
  ENDIF

!   do i = 1, Nboites
!     write(*,'(" Boite : ",I5,"  par to be tested : ",I6," particules")') I,Npartester(I)
!   enddo
!   read(*,*)

endif

#ifdef PA
! Determination of the boxes close to the boundary (case Nb_phase = 2)
if (Nb_phase == 2) then

  ! Teste sur des joints type (100) et (110): OK
  ! Les conditions aux limites periodiques ne sont par contre pas appliquees dans la recherche de boite

  RBOX(:) = tailleBoite(:) * Half       ! The box center

  ! Only one boundary can be considered in the present form of the code
  if (NbPlan > 1 .and. Nb_Phase > 1) &
  stop 'The procedure used to determined the boxes close to the boundaries was tested for 1 boundary only!'

  do jj = 1,NbPlan


    ! Initializations
    list_boxes(:) = 0
    ii = 0
    normsurf(1:3)=normavec(Plane_MillerR(1:3,jj))

    ! Only the case with a boundary oriented along a <100> direction and as far from the simulated volume limits
    ! is clearly tested. For other cases, some mistakes could appera because the CLP are not taken into account in
    ! the calculation of list_box
    if (normsurf(1)*normsurf(2) /= 0 .or. normsurf(2)*normsurf(3) /= 0 .or. normsurf(1)*normsurf(3) /= 0) then
      write(*,*)
      write(*,*) "Only the case with a boundary normal parallel to a box axis (kind <100>) and at least"
      write(*,*) "far from L_Boite to the volume bordure is safely programmed!!"
      write(*,*)
      write(*,*) "With another kind of axis you may have errors in the determination of the boxes close to "
      write(*,*) "the boundary because the CLP are not taken into account on the calculation."
      write(*,*)
      stop
    elseif ((normsurf(1) * Plane_pos(jj) + &
             normsurf(2) * Plane_pos(jj) + &
             normsurf(3) * Plane_pos(jj) <= L_Boite) .or. &
            (modur(1) - abs(normsurf(1) *  Plane_pos(jj)) + &
             modur(2) - abs(normsurf(2) *  Plane_pos(jj)) + &
             modur(3) - abs(normsurf(3) *  Plane_pos(jj)) <= L_Boite)) then
      write(*,*)
      write(*,*) "Only the case with a boundary placed at least far from L_Boite to the volume bordure is"
      write(*,*) "safely programmed!!"
      write(*,*)
      write(*,*) "With another kind of axis you may have errors in the determination of the boxes close to "
      write(*,*) "the boundary because the CLP are not taken into account on the calculation."
      write(*,*)
      stop
    endif

    ! Highest length between the box center and the boundary along the boundary normal direction, for a
    ! box touching/crossing the boundary
    highest_length = abs(RBOX(1) * normsurf(1)) + abs(RBOX(2) * normsurf(2)) + abs(RBOX(3) * normsurf(3))

    ! We want to select the boxes containing the boundary and their neighbours
    ! --> highest possible length = 3 * highest_length (case of a box those the neighobouring box is
    ! tangent to the boundary)
    highest_length = 3 * highest_length

    ! Loop on the boxes
    do J = 1,Nboites

      IB = B1D_3D(J,1)
      JB = B1D_3D(J,2)
      KB = B1D_3D(J,3)

      ! Coordinate of the box centers calculation
      ! Be careful: the boxcenter could be outside the frame
      Boxcenter(1) = ((IB-1) + HALF) * tailleBoite(1)
      Boxcenter(2) = ((JB-1) + HALF) * tailleBoite(2)
      Boxcenter(3) = ((KB-1) + HALF) * tailleBoite(3)

      ! Test: (i,j,k) are the Miller indices of the plan
      ! d is the position such as i*X1 + j*X2 + k*X3 = d
      ! if i*boxcenter(1) + j*boxcenter(2) + k*boxcenter3 - d' <= the highest distance
      ! between the center and a point of the box,
      ! then the possibility that the box is in the plan neighbourhood exists
      ! d' = d * normsurf(i) in order
      if (normsurf(1)/=0) then
        if (abs(normsurf(1) * Boxcenter(1) +                                              &
                normsurf(2) * Boxcenter(2) +                                              &
                normsurf(3) * Boxcenter(3) -                                              &
            dsign(un,normsurf(1)) * normsurf(1) *  Plane_pos(jj)) <= highest_length) then

          ii = ii + 1
          list_boxes(ii) = J

        endif
      elseif (normsurf(2)/=0) then
        if (abs(normsurf(1) * Boxcenter(1) +                                              &
                normsurf(2) * Boxcenter(2) +                                              &
                normsurf(3) * Boxcenter(3) -                                              &
            dsign(un,normsurf(2)) * normsurf(2) *  Plane_pos(jj)) <= highest_length) then

          ii = ii + 1
          list_boxes(ii) = J

        endif
      else
          if (abs(normsurf(1) * Boxcenter(1) +                                            &
                  normsurf(2) * Boxcenter(2) +                                            &
                  normsurf(3) * Boxcenter(3) -                                            &
              dsign(un,normsurf(3)) * normsurf(3) *  Plane_pos(jj)) <= highest_length) then

            ii = ii + 1
            list_boxes(ii) = J

          endif
      endif

      NBoites_boundary = ii ! number of boxes close to the boundary for this phase

    enddo
  enddo


endif  ! fin allocation Nb_phase=2
#endif

deallocate(boxborder)
deallocate(key_box_border)

end subroutine Allocation_initiale

!######################################################################################
! The subroutine to compute the optimum number of Greengard domains during simulation #
!######################################################################################
subroutine Allocation_dynamique

implicit none

integer(kind=DPI)     :: IX, IY, IZ
integer(kind=DPI)     :: I, J, K
integer(kind=DPI)     :: BX, BY, BZ
integer(kind=DPI)     :: tmp, nparmax_boite, compteur, nvNBoites
integer(kind=DPI)     :: X(3)
logical               :: EmptyBoite

! New evaluation of nboite is made each time the number of segment is increased by 20%
if (abs(nsegm - nsegINI) < int(0.2 * float(nsegINI))) RETURN

! The last tested value of nsegm (number of segment used to evaluate Greengard volumes) is updated
nsegINI = nsegm

! The magic formula to calculate the optimum number of domains in the Greengard algo
nvnboites = nint(DSQRT(real(KRC*27*NSEGM,DP)))

! The two stupid cases
if(nvnboites >  nsegm) nvnboites  = nsegm   ! Not enough segments per domains to be efficient
if(nvnboites < 64) nvnboites = IUN          ! Not enough domains to be helpful

! The numbers of optimum domains must increase significantly to find a new tilling of the simulated volume
if (ICINQ * abs(nvnboites - nboites) < nboites) RETURN

! We should be able to find a solution, so we must allocate new tabs
deallocate (NsegBoite)
do i = 1, NBoites
  deallocate(IndexBoite(i)%ListeSeg)
enddo
deallocate (IndexBoite)
deallocate (B3D_1D)
deallocate (B1D_3D)
deallocate (IndBoiteVois)
deallocate (SIGBOX) ! contrainte aucentre des boites
deallocate (NBox2Shift)
deallocate (Box2Shift)
deallocate (Key_Box_2Border)

if (particules) then
  ! calcul du nombre moyens de particules par boite
  Deallocate (NparBoite)
  Deallocate (listeparBoite)
  Deallocate (Npartester)
  Deallocate (listepartester)
endif

! Calculation of the new optimum L_Boite
L_boite = (Volume / NVnBoites)**(TIERS)

! If the simulation box contains particles, domains cannot be smaller than the diameter of particles
if (particules) then
  if (L_boite < maxval(par(1:npar)%R)) then
    write(*,'("A particle is too large for optimum multipole domains!" )')
    print*,"You must run the simulation by imposing constant L_boite = ",(1.2*maxval(real(par(1:npar)%R,DP))*avalue)
    stop
  endif
endif

! The minimum size of the Greengard domains must be at least the discretization length
if (L_boite < deux * Xlomax) then
  write(*,*) "--------------------------------------------------------------"
  write(*,'( "Discretization length is too large for optimum multi-pole domains!" )')
  write(*,*) "Xlomax =", Xlomax, "  optimum L_domains =" , l_boite
  L_boite =  deux * Xlomax
  write(*,*) "Multipole L_domains is increased to : ", L_boite, '==>', (L_boite*avalue * 1.e6), 'micron'
  write(*,*) "--------------------------------------------------------------"
endif

!************************************
! MAILLAGE POUR LA METHODE DES BOITES
!************************************

! if the discretization length is to large we stop the Greengard algo
if ( .not. methode_boites .or.  (trois * L_boite) > real(minval(modur(:)))) then

  write (*,*) " Multipole domains are not used", (deux * L_boite),minval(modur(:))
  NBoites = IUN
  NBoitesX = IUN
  NBoitesY = IUN
  NBoitesZ = IUN

  nsegmax_boite = nsegmax
  nparmax_boite = npar
  RCMAKE = .false.

else

  ! nombre approximative des boites dans les 3D
  NBoitesX = int(Modur(1)/L_boite)
  NBoitesY = int(Modur(2)/L_boite)
  NBoitesZ = int(Modur(3)/L_boite)

  ! The new accepted NBoites value
  NBoites = NBoitesX*NBoitesY*NBoitesZ

  if (GB /= 0) then
    nsegmax_boite = (nsegmax / IDIX)
  else
    nsegmax_boite = (nsegmax / NBoites + 1) * ICENT
  endif
  if (nsegmax_boite < ListeSegRangeIni) nsegmax_boite = ListeSegRangeIni + icent
  nparmax_boite = ((npar_max / NBoites) + 1) * IDIX

  ! The particular case when the number of domains in one direction is less than 4
  ! Then the following coding does not work and the multipole algo must not be applied
  if (NBoitesX < 4 .or. NBoitesY < 4 .or. NBoitesZ < 4) then
    write (*,*) " One dim of the simulated volume is too small for Multipoles", NBoitesX, NBoitesY, NBoitesZ
    NBoites  = IUN
    NBoitesX = IUN
    NBoitesY = IUN
    NBoitesZ = IUN

    nsegmax_boite = nsegmax
    nparmax_boite = npar
  endif

endif

! Calculation of the exact size of the boxes.
tailleBoite(1) = real(Modur(1),DP) / NBoitesX
tailleBoite(2) = real(Modur(2),DP) / NBoitesY
tailleBoite(3) = real(Modur(3),DP) / NBoitesZ

write(*,'("Setting : Nb Greengard domains",I6," (",I2,",",I2,",",I2,") Size (",ES10.2E2,",",ES10.2E2,",",ES10.2E2,")")') &
     NBoites,NBoitesX,NBoitesY,NBoitesZ,TailleBoite(1:3)

if(Xlomax > int(minval(TailleBoite(:)))) then
   write(*,'("After correction : Discretization length > 0.5 Greengard domains")')
   print *, "This is a problem STOP!"
   stop
endif

! verification ultime
if( (real(Modur(1),DP) - TailleBoite(1)*NBoitesX) > numtols .or.  &
    (real(Modur(2),DP) - TailleBoite(2)*NBoitesY) > numtols .or.  &
    (real(Modur(3),DP) - TailleBoite(3)*NBoitesZ) > numtols )     &
   stop "Identification of the Greengard domains size failed STOP!"

! allocation des tableaux pour la methode de boites
allocate (NsegBoite(NBoites))
allocate (IndexBoite(NBoites))
allocate (B3D_1D(NBoitesX,NBOITESY,NBoitesZ))
allocate (B1D_3D(NBoites,3))
allocate (SIGBOX(3,3,NBoites)) ! contrainte aucentre des boites
allocate (IndBoiteVois(NBoites,27))

if (nboitesX < 4 .or. nboitesY < 4 .or. nboitesZ < 4) then
  nboxborder = 1
else
  nboxborder = 2*((NboitesX-2)*(NboitesY-2)+(NboitesX-2)*(NboitesZ-2)+(NboitesY-2)*(NboitesZ-2))+   &
               4*(NboitesX+NboitesY+NboitesZ-6) + 16 +                                              &
               2*((NboitesX-4)*(NboitesY-4)+(NboitesX-4)*(NboitesZ-4)+(NboitesY-4)*(NboitesZ-4))+   &
               4*(NboitesX+NboitesY+NboitesZ-12)
endif

allocate (boxborder(nboxborder))
allocate (NBox2Shift(NBoites))
allocate (Box2Shift(NBoites,26))
allocate (Key_Box_2Border(NBoites))
allocate (Key_Box_Border(NBoites))

if(NBoites /= IUN) then

  RCMAKE = .true.   ! clef pour forcer le calcul de sigma INT CP
  Do IX = 1,NBoitesX
    Do IY = 1,NBoitesY
      Do IZ = 1,NBoitesZ
        ! cette formule represente une bijonction entre Iboite et les trois indices des boites
        ! c est utile surtout dans contact parce qu'il faut deplacer les segemnt par ordre predefini
        J = IX + (IY-1) * NBoitesX + (IZ-1)* (NBoitesX*NBoitesY)
        B3D_1D(IX,IY,IZ) = J
        B1D_3D(J,1) = IX
        B1D_3D(J,2) = IY
        B1D_3D(J,3) = IZ

        ! calcul du tableau logique de voisinage de sous_boites
        ! chaque boite doit avoir 26 boites + elle-meme comme voisines
        ! Ix,IY,IZ sont les indice de chaque boite
        compteur = iun
        IndBoiteVois(B3D_1D(IX,IY,IZ),compteur) = B3D_1D(IX,IY,IZ) ! premiere boite consideree voisine de BI est BI elle-meme

        do I = -IUN , IUN
          do J = -IUN , IUN
            do K =  -IUN , IUN
              ! desole : les CPL shift ou tourne ne sont pas pris en compte
              BX = modulo((IX + I - IUN),NBoitesX) + IUN
              BY = modulo((IY + J - IUN),NBoitesY) + IUN
              BZ = modulo((IZ + K - IUN),NBoitesZ) + IUN
              ! print *, IX,IY,IZ,BX,BY,BZ
              tmp = BX + (BY-1) * NBoitesX + (BZ-1)* (NBoitesX*NBoitesY)
              if (I /= izero .or. J/= izero .or. K /= izero) then
                compteur = compteur + iun
                IndBoiteVois(B3D_1D(IX,IY,IZ),compteur) = tmp
              endif
              if (compteur > 27 ) stop 'Error in the IndBoiteVois tabulation STOP!'

            enddo
          enddo
        enddo

      enddo
    enddo
  enddo

else

  RCMAKE = .false.   ! clef pour forcer le calcul de sigma INT CP
  B3D_1D(:,:,:) = IUN
  B1D_3D(IUN,:) = IUN
  IndBoiteVois(:,:) = iun

endif

!sq In what follows, we define the box which are on a border of the simulated volume
! this will be useful for the obstacle detection. Because of the use of PBC, segments
! that are near a border of the simulation volume can be inside the volume itself or
! in its replica. This is why a shift procedure is applied before contact detection to
! define the smallest distance between a segment i and a potential obstacle j.

!initialization
NBoxBorder = izero !total number of boxes that are at the border of the simulation volume
BoxBorder(:) = izero !tab of boxes that are at the border of the simulation volume
NBox2Shift(:) = izero ! total number of boxes to shift for a given box BI at the border of the simulation box
Box2Shift(:,26) =izero ! no. of boxes to shift for a given box BI at the border of the simulation box
Key_Box_Border(:) = .false. !key is true when a box is on a border of the simulated volume
Key_Box_2Border(:) = .false. !key is true when a box is on a border or at the second position near a border
                             ! of the simulated volume

if (nboitesX > 3 .and. nboitesY > 3 .and. nboitesZ > 3) then ! if there are less than 5 boxes per sides of the simulation volume, this procedure is useless
  do J = 1, NBoites
    IX = B1D_3D(J,1)
    IY = B1D_3D(J,2)
    IZ = B1D_3D(J,3) ! no. of the considered box

    !#1 In a first step, we make a tab of the boxes that are near the border
    if (IX==1 .or. IX==nboitesX .or. iY==1 .or. IY==nboitesY .or. IZ==1 .or. IZ==nboitesZ) then
      Key_Box_Border(J)   = .true.
      Key_Box_2Border(J)  = .true.
      NBoxBorder = NBoxBorder + 1
      BoxBorder(NBoxBorder) = J
      ! print *, '#',NBoxBorder , J, ix, iy, iz
    elseif (IX==2 .or. IX==nboitesX-1 .or. iY==2 .or. IY==nboitesY-1 .or. IZ==2 .or. IZ==nboitesZ-1)  then
      Key_Box_2Border(J) = .true.
      NBoxBorder = NBoxBorder + 1
      BoxBorder(NBoxBorder) = J
      ! print *, '#',NBoxBorder , J, ix, iy, iz
    endif
  enddo

  ! A test to ensure that the number of box near the borders of the simulated volume is correct
  compteur = 2*((NboitesX-2)*(NboitesY-2)+(NboitesX-2)*(NboitesZ-2)+(NboitesY-2)*(NboitesZ-2))+   &
             4*(NboitesX+NboitesY+NboitesZ-6) + 16 +                                              &
             2*((NboitesX-4)*(NboitesY-4)+(NboitesX-4)*(NboitesZ-4)+(NboitesY-4)*(NboitesZ-4))+   &
             4*(NboitesX+NboitesY+NboitesZ-12)
  if (NBoxBorder /= compteur) then
    print *, 'There is a problem with the number of Greengard domains!', NBoxBorder, compteur
    stop
  endif

  !#2 In a second step, for all the boxes that are near a border, we define the tab
  ! of boxes that have to be shifted
  do I = 1, NBoxBorder
    J = BoxBorder(I) ! N° of the box near a border
    IX = B1D_3D(J,1)
    IY = B1D_3D(J,2)
    IZ = B1D_3D(J,3)

    !initialization
    NBox2Shift(J)= 0
    Box2Shift(J, 1:26)= 0

    ! among all the neighboors of J, we retain the boxes that are AT the
    ! border of the simulted volume (J itself also)
    do compteur = 1, 27
      k = IndBoiteVois(J,compteur) ! K a neighboor of J
      if (Key_Box_Border(K)) then
        NBox2Shift(J) = NBox2Shift(J) + iun
        Box2Shift(J,NBox2Shift(J)) = K
      endif
    enddo

    if (NBox2Shift(J) < 7 .or.NBox2Shift(J) > 26) then
      print *, 'PB, nbr of Box2shift > 26' , NBox2Shift(J)
      stop
    endif

  enddo

else ! there are less than 4 boxes per sides of the simulation volume
  ! we gain nothing with this during shift procedure
  Key_Box_2Border(:) = .false.

endif

! The IndexBoite(i)%ListeSeg can now be allocated
! domains (boites) outside of the effective dislocated volume are set to dim=1 to minimize the total tab size
Do IX = 1,NBoitesX
  Do IY = 1,NBoitesY
    Do IZ = 1,NBoitesZ

      EmptyBoite = .false.
      I = B3D_1D(IX,IY,IZ)

      if (GB == 1) then

        ! Frontiers between boites and the effective simulation volume are tested
        if ( IX * tailleBoite(1) < dom_effective_deb(1)) EmptyBoite = .true.
        if ( IY * tailleBoite(2) < dom_effective_deb(2)) EmptyBoite = .true.
        if ( IZ * tailleBoite(3) < dom_effective_deb(3)) EmptyBoite = .true.
        if ( (IX-1) * tailleBoite(1) > dom_effective_end(1)) EmptyBoite = .true.
        if ( (IY-1) * tailleBoite(2) > dom_effective_end(2)) EmptyBoite = .true.
        if ( (IZ-1) * tailleBoite(3) > dom_effective_end(3)) EmptyBoite = .true.

      endif

      if (EmptyBoite) then
        allocate(IndexBoite(I)%ListeSeg(1))
      else
        allocate(IndexBoite(I)%ListeSeg(NSEGMAX_boite))
      endif

    enddo
  enddo
enddo

! initialization
SIGBOX(:,:,:)             = izero
NsegBoite(:)              = izero
do i = 1,NBoites
  IndexBoite(i)%ListeSeg(:) = izero
enddo

! CLP shiftees pour methode des boites: tableau de vecteurs elementaires de shifts
! en unites de taille de sous-boite. l'arrondi ne devrait pas avoir d'effet
! sur le calcul des contraintes
! i=1 ---> sortie face a ;! i=2 ---> sortie face b ; i=3 ---> sortie face c
! sur le calcul des contraintes
! i=1 ---> sortie face a ; i=2 ---> sortie face b; i=3 ---> sortie face c
do i=1,3,1
  SHIFTBOITE(1,i) = int(SHIFT(1,i)/tailleboite(1))
  SHIFTBOITE(2,i) = int(SHIFT(2,i)/tailleboite(2))
  SHIFTBOITE(3,i) = int(SHIFT(3,i)/tailleboite(3))
enddo

write(*,'("Number of multipole domains ",I6," (",I2,",",I2,",",I2,") Nseg=",I6)')NBoites,NBoitesX,NBoitesY,NBoitesZ,nsegm
write(*,'("Size of multipole domains in AVALUE (",ES10.2E2,",",ES10.2E2,",",ES10.2E2,")")')tailleBoite(1:3)
Do I = 1, NBoitesX
  Do J = 1, NBoitesY
    Do K = 1, NBoitesZ
      if(NsegBoite(B3D_1D(i,j,k)) /= izero) &
           write(*,'(" Boite : (",I2,",",I2,",",I2,") contains",I6," particules")')&
           I,j,k,NsegBoite(B3D_1D(i,j,k))
    enddo
  enddo
enddo

! allocation des tableaux pour la methode de boites appliquee aux particules
if (particules) then

  ! calcul du nombre moyens de particules par boite

  allocate (NparBoite(NBoites))
  allocate (listeparBoite(NBoites,nparmax_boite))
  allocate (Npartester(Nboites))
  allocate (listepartester(Nboites,10*nparmax_boite))

  NparBoite(:)        = izero
  ListeparBoite(:,:)  = izero
  IBOITEP(1:Npar)     = izero

  do I =  1, Npar

    X(:) =  par(i)%C(:)

    if (Nboites /= iun) then
      ! le modulo permet de respecter les CPL pour les boite
      Bx = modulo(int(X(1) / tailleBoite(1),DPI),NboitesX) + IUN
      By = modulo(int(X(2) / tailleBoite(2),DPI),NboitesY) + IUN
      Bz = modulo(int(X(3) / tailleBoite(3),DPI),NboitesZ) + IUN
      J = B3D_1D(Bx,By,Bz)

      if(NparBoite(J) > nparmax_boite) then
        print *,"nparmax_boite =",nparmax_boite, " npar  =",npar
        print *,"J,IB,JB,KB,",J,Bx,By,Bz
        print *,"NsegBoite(IB,JB,KB)",NparBoite(J)
        print *,"Error in REGION : NparBoite is not correct, kk=", kk
        stop
      endif

      IBOITEP(I) = J                            ! Boite de la part
      NparBoite(J) =  NparBoite(J) + IUN        ! Nombre de part dans une boite
      Listeparboite(J,NparBoite(J)) = I         ! Liste de part dans une boite

    else

      ! WHEN NO MULTIDOMAINE applied, all segments are tested every where
      NparBoite(IUN) = npar
      listeparBoite(IUN,I) = I
      IBOITEP(I) = IUN

    endif

    ! write(*,'("part index:",I5," box index including part:",I5," part number in the box:",I5)') I,IBoiteP(I), NparBoite(IBOITEP(I))

  enddo

  !==================================================================
  ! verification de la procedure
  J= izero
  do i = 1, Nboites
      J = J + Nparboite(I)
  enddo

  if (npar /=  J) stop "A problem with the Greengard method with particles"
  !==================================================================

  Npartester(:) = izero
  Listepartester(:,:) = izero

  if(nboites > iun) then

    do I = 1, NBOITES

      do J = 1, NBOITES

        ! On stocke les particules visites dans le tableau listeparteste.
        if(.not. any(IndBoiteVois(I,:)-J .eq. IZERO)) CYCLE
        if(NparBoite(J) == Izero) CYCLE

        do k = 1, NparBoite(J)
          tmp = listeparboite(J,k)
          Npartester(I) = npartester(I) + IUN
          listepartester(I,npartester(I)) = tmp
        enddo

      enddo

      if(Npartester(I) > 10*nparmax_boite) then
          print *, " Erreur: topolo: debordement de tableau listeparboite:",I
          stop
      endif
    enddo

  elseif(NBoites == iun) then

    Npartester = Npar
    !**** nouveau CHANGEMENT DE REPERE  et + rapide
    do tmp = 1 , Npar
        listepartester(IUN,tmp) = tmp
    enddo
  else
    stop " Undefined Nboite in TOPOLO"
  endif

!   do i = 1, Nboites
!     write(*,'(" Boite : ",I5,"  par a tester : ",I6," particules")') I, Npartester(I)
!   enddo

endif

deallocate(boxborder)
deallocate(key_box_border)

end subroutine Allocation_dynamique

!##############################################################################
!> Discretization :
!>
!> Principle : When a segment is too long, we cut it.
!> As a matter of fact : We particularly look at the dislocation curve (i.e. the neighboring segments)
!> so that the discretized picture is not too far from the (continuous) reality.
!> In addition : Plenty of specific cases !
!>
!##############################################################################
subroutine DISCRETI

use CONNEC

implicit none

integer(kind=DPI) :: ILONG    !< Norme of segment i
integer(kind=DPI) :: JLONG    !< Norme of segment j
integer(kind=DPI) :: segi     !< Variable for segment indices
integer(kind=DPI) :: segj     !< Variable for segment indices
integer(kind=DPI) :: itemp    !< Temporary variable for segment indices
integer(kind=DPI) :: VLi      !< Vector number for segment i (VECLIN)
integer(kind=DPI) :: VLj      !< Vector number for segment j (VECLIN)
integer(kind=DPI) :: i1       !<
integer(kind=DPI) :: i2       !<
integer(kind=DPI) :: j1       !<
integer(kind=DPI) :: j2       !<
integer(kind=DPI) :: K        !<
integer(kind=DPI) :: IBV      !<
integer(kind=DPI) :: bvnno    !<
integer(kind=DPI) :: bvnne    !<
integer(kind=DPI) :: ve       !<
integer(kind=DPI) :: vo       !<
integer(kind=DPI) :: TOI(3)   !<
integer(kind=DPI) :: HandIo   !<
integer(kind=DPI) :: HandIe   !<
integer(kind=DPI) :: HandJo   !<
integer(kind=DPI) :: HandJe   !<
integer(kind=DPI) :: compteur !<
integer(kind=DPI) :: Iside    !<
integer(kind=DPI) :: RefEi    !<
integer(kind=DPI) :: RefOj    !<
integer(kind=DPI) :: RefEj    !<
integer(kind=DPI) :: CNN      !<
integer(kind=DPI) :: LIO      !<
integer(kind=DPI) :: LJO      !<
integer(kind=DPI) :: LIE      !<
integer(kind=DPI) :: LJE      !<
integer(kind=DPI) :: LOMAI    !<
integer(kind=DPI) :: lomaIR   !<
integer(kind=DPI) :: seggd    !<
integer(kind=DPI) :: lgd      !<
integer(kind=DPI) :: iloop    !< a loop index used in the microloop annihilation
integer(kind=DPI) :: NextSeg  !< a segment index used in the microloop annihilation

real (kind=DP)    :: Loop_length    !< The cumulated length of dislocation line in the microloop test

integer(kind=DPI) , dimension(3) :: Oi        !< Origin of segment i
integer(kind=DPI) , dimension(3) :: Oj        !< Origin of segment j
integer(kind=DPI) , dimension(3) :: Ei        !< End of segment i
integer(kind=DPI) , dimension(3) :: Ej        !< End of segment j
integer(kind=DPI) , dimension(3) :: DECAL     !<
integer(kind=DPI) , dimension(3) :: totoBUGv1 !<
integer(kind=DPI) , dimension(3) :: totoBUGv2 !<
integer(kind=DPI) , dimension(20):: ListA     !< List of segments building a microloop
integer(kind=DPI) , dimension(20):: ListB     !< List to test the inclusion of a segment to ListA

logical           :: Recouvrement !<
logical           :: condition    !<
logical           :: JunctLenMax  !<
logical           :: closed_loop  !< a microloop is found

! Initialization
PlusSeg             = izero
oi(:)               = izero
decal(:)            = izero
k                   = nsegm / idix
OUT(:)              = .false.
GD2made             = .false.
Seg(Nsegm+1:nsegm+K)%probadev = zero

!------------------------------------------------------------------------------------
! pour que la procedure de boite soit pertinante, il ne faut pas laisser les segments
! s'allonger sans limite...: ici on prend la limite 2 LOMA(I)
! on decoup par trois sur trois tours (parce qu'on ne peut pas changer les bornes
! d'une boucle de facon dynamique
do Segi = 1 , Nsegm

  ! We can pass these segments
  if (segi == seg(Segi)%voiso .and. segi == seg(segi)%voise) then
    out(segi) = .true.        ! a paranoid update
    cycle  ! Segment to eliminate
  endif

  if (seg(segi)%unload) cycle                                       ! segment part of an infinit line and pinned

  ! The wait algo is re-initialized on all segments if Period /= 0
  if (Period /= IZERO) seg(Segi)%wait = modulo(seg(Segi)%wait,Period)

  ilong = seg(SEGI)%norme

  if (seg(Segi)%jonc) then

    segj = seg(Segi)%Ijonc

    ! a paranoid test
    if (.not. seg(SEGJ)%jonc .or. seg(segj)%ijonc /= Segi) then
      call conf(Segi)
      call conf(segj)
      print *, " Discreti : Jonction assymetrique entre i=",Segi," et j :", seg(Segi)%ijonc
      stop
    endif

  else

    ! Nothing to do with pivot segments
    if (ilong == izero) cycle

    IBV   = seg(SEGI)%veclin
    lomaI = loma(IBV)
    vo    = seg(SEGI)%voiso
    bvnno = seg(SEGI)%vnno
    ve    = seg(SEGI)%voise
    bvnne = seg(SEGI)%vnne

    ! Nothing to do with very small segments
    if(ilong < ITROIS) cycle

    ! Critical length to introduce a new pivot segment at pinning points: O side
    if (2*ilong > LomaI .and. vo == nsegmax) then

      call rotuleO(Segi,3422)
      vo    = seg(SEGI)%voiso
      bvnno = seg(SEGI)%vnno
      ve    = seg(SEGI)%voise
      bvnne = seg(SEGI)%vnne

    endif

    ! Critical length to introduce a new pivot segment at pinning points: E side
    if (2*ilong > LomaI .and. ve == nsegmax) then

      call RotuleE(Segi,3421)
      vo    = seg(SEGI)%voiso
      bvnno = seg(SEGI)%vnno
      ve    = seg(SEGI)%voise
      bvnne = seg(SEGI)%vnne

    endif

    ! We force discretization if
    ! case 1) length I > 2 x  lomaI
    ! case 2) length I > lomaI and the segment at O or E of I is a pinning point
    ! case 3) seg(i)%diseg discretization is imperatively needed for that segment and long > 3
    ! case 4) I > 0.5 lomaI and the segment at O or E of I is immobile (jonc ou GD) or the segment touches the surface.
    ! case 5) Segment I is touching on both sides, O and E, two junctions
    ! Case 6) The segment is an infinit segment created with the PBC
    condition = .false.
    if (SEGI > iforet ) then
      condition =                                                                                             &
        (Ilong > 2*lomai                                                                                .or.  & ! Case 1
        ((vo == nsegmax .or. ve == nsegmax) .and. Ilong > lomai)                                        .or.  & ! Case 2
        (seg(Segi)%diseg .and. ((Ilong > ITROIS .and. seg(segi)%surface < IUN) .or. Ilong > ICINQ ))    .or.  & ! Case 3
        (((seg(bvnno)%gd + seg(bvnne)%gd) > izero .or. seg(bvnno)%jonc .or.                                   &
           seg(bvnne)%jonc .or. seg(segi)%surface > IZERO)                                              .and. &
         IDEUX*Ilong > lomai)                                                                           .or.  & ! Case 4
        (seg(bvnno)%jonc .and. seg(bvnne)%jonc)                                                         .or.  & ! Case 5
        (seg(Segi)%vnno == seg(Segi)%vnne .and. seg(Segi)%veclin == seg(seg(Segi)%vnno)%veclin))                ! Case 6
    endif

    if (.not. condition) cycle     ! Segment i must not be discretized

    ! Information relative to segment I
    Oi(:)          = seg(SEGI)%O(:)
    decal (:)      = bveclin(:,ibv)
    Ei(:)          = Oi(:) + Ilong*decal(:)
    Ilong          = Ilong / ITROIS


    ! The two new segments i1 and i2
    Plusseg  = plusseg + ideux
    i1       = Nsegm + PLUSSEG - IUN
    i2       = Nsegm + PLUSSEG

    if (kkdebug) then
      write(379,*) "***Discreti of",SEGI, " Long=",  Ilong, " diseg =",seg(Segi)%diseg
      write(379,*) "new segs : I1 =",I1," et I2 = ",I2, " plusseg",plusseg
    endif

    ! Definition of segment i1 and i2
    seg(i1)%dom      = seg(SEGI)%dom
    seg(i1)%norme    = Ilong
    seg(i2)%norme    = Ilong
    seg(SEGI)%norme  = seg(Segi)%norme - Ideux * Ilong
    seg(I1)%veclin   = IBV
    seg(I2)%veclin   = IBV
    seg(I1)%gd       = izero
    seg(I2)%gd       = izero
    seg(i1)%jonc     = .false.
    seg(i1)%Ijonc    = nsegmax
    seg(i2)%jonc     = .false.
    seg(i2)%Ijonc    = nsegmax
    seg(I1)%wait     = IZERO
    seg(I2)%wait     = IZERO
    seg(I2)%bloquer  = .false.     ! Better to initialize this quantity since we do not know
    seg(I2)%unload   = .false.     ! which part of the initial segment I was touching a barrier
    seg(I1)%unload   = .false.
    seg(Segi)%unload = .false.
    SEG(I1)%grain    = seg(SEGI)%grain
    seg(I2)%grain    = seg(SEGI)%grain
    seg(I1)%anglevis = seg(SegI)%anglevis
    seg(I2)%anglevis = seg(SegI)%anglevis
    singulier(i1)    = .false.
    singulier(segi)  = .false.
    singulier(i2)    = .false.

    !Which segment is touching free surface after discretization
    select case (seg(SEGI)%surface)

      ! the segment is not touching a surface
      case (0)
        seg(i1)%surface       = IZERO
        seg(i1)%VarFreePlan   = IZERO
        seg(SEGI)%surface     = IZERO
        seg(SEGI)%VarFreePlan = IZERO
        seg(i2)%surface       = IZERO
        seg(i2)%VarFreePlan   = IZERO

      ! the origin of segment is touching a surface
      case (1)
        seg(i1)%surface       = seg(SEGI)%surface
        seg(i1)%VarFreePlan   = seg(segi)%VarFreeplan
        seg(i1)%zerotl        = seg(segi)%zerotl
        singulier(i1)         =.true.
        seg(segi)%VarFreeplan = IZERO
        seg(segi)%zerotl      =.false.
        seg(SEGI)%surface     = IZERO
        seg(i2)%surface       = IZERO
        seg(i2)%VarFreePlan   = IZERO
        if (kkdebug) write(379,*) "Discretization of surface segment type 1",i1,seg(i1)%surface,SEGI,kk

      ! the end of segment is touching a surface
      case(2)
        seg(i2)%surface       = seg(SEGI)%surface
        seg(i2)%VarFreePlan   = seg(segi)%VarFreeplan
        seg(i2)%zerotl        = seg(segi)%zerotl
        singulier(i2)         = .true.
        seg(segi)%VarFreeplan = IZERO
        seg(segi)%zerotl      =.false.
        seg(SEGI)%surface     = IZERO
        seg(i1)%surface       = IZERO
        seg(i1)%VarFreePlan   = IZERO
        if (kkdebug) write(379,*) "Discretization of surface segment type 2",  &
                                  i2,seg(i2)%surface,SEGI,seg(SEGI)%surface,kk

      case DEFAULT
        print *,'Discretization part, problem with',segi,' seg(i)%surface at step',kk
        call seginfo (segi,'surface problem 1')
        stop  'surface problem 1'
    end select

    ! attention: si i est GD, alors il faut mettre les seg I1 et I2 (non GD) dans le systeme
    ! de glissement des bvnno, bvnne
    if(seg(Segi)%gd > izero .and. bvnno /= nsegmax) seg(I1)%veclin = segdevIJ(IBV,seg(bvnno)%veclin)
    if(seg(Segi)%gd > izero .and. bvnne /= nsegmax) seg(I2)%veclin = segdevIJ(IBV,seg(bvnne)%veclin)

    ! calcul des origines
    seg(i1)%O(:)   = Oi(:)
    seg(Segi)%O(:) = Oi(:) + Ilong * decal(:)
    seg(i2)%O(:)   = Oi(:) + (Ilong + seg(SEGI)%norme) * decal(:)

    ! traitement de connectivite entre :  Vo et I1
    seg(i1)%voiso = vo
    seg(i1)%voise = SegI
    seg(i1)%vnno  = bvnno
    seg(i1)%vnne  = SegI
    if (vo /= nsegmax) seg(vo)%voise = i1

    ! traitement de connectivite entre :  I1 et SegI
    seg(Segi)%voiso = I1
    seg(Segi)%voise = I2
    seg(Segi)%vnno  = I1
    seg(Segi)%vnne  = I2

    ! traitement de connectivite entre :  SegI et I2
    seg(i2)%voiso = SegI
    seg(i2)%voise = ve
    seg(I2)%vnno  = SegI
    seg(I2)%vnne  = bvnne

    ! traitement de connectivite entre :  I2 , Ve
    if (ve /= nsegmax) seg(ve)%voiso = i2

    seg(i1)%resdep = seg(Segi)%resdep
    seg(i2)%resdep = seg(Segi)%resdep
#ifdef MDC
    seg(i1)%SIGFE(:) = seg(Segi)%SIGFE(:)
    seg(i2)%SIGFE(:) = seg(Segi)%SIGFE(:)
#endif
    seg(i1)%jonc      = .false.
    seg(i2)%jonc      = .false.
    seg(i1)%gd        = izero
    seg(i2)%gd        = izero
    seg(I1)%wait      = izero
    seg(I2)%wait      = izero
    seg(I2)%bloquer   = .false.     ! exeptionnellement on reinitialise
    seg(I1)%bloquer   = .false.     ! pour la discretisation forcee
    seg(Segi)%bloquer = .false.     ! sinon on a des pb au barrieres
    seg(I2)%unload    = .false.
    seg(I1)%unload    = .false.
    seg(Segi)%unload  = .false.
    seg(I1)%grain     = seg(SEGI)%grain
    seg(I2)%grain     = seg(SEGI)%grain

    !*verifiaction pour le pavage
    if (DesorientGrain) then
      call CorrectifGrain(SEGI,3002)
      call CorrectifGrain(I1,3003)
      call CorrectifGrain(I2,3004)
    endif

    if (Nbcvxdom > IUN) then
      if ( seg(segi)%dom == seg(ve)%dom .and. seg(ve)%vnne/=nsegmax) then
        if(Any(Oi(:) /= modulo( Oi(:),modur(:) ) ) .or. &
           Any(Ei(:) /= modulo( Ei(:),modur(:) ) ) ) then
           call assign_domain(I2,30001)
           call assign_domain(SEGI,30002)
        else
          seg(I1)%dom=seg(SEGI)%dom
          seg(I2)%dom=seg(SEGI)%dom
        endif
      else
        call assign_domain(I2,30001)
        call assign_domain(SEGI,30002)
      endif ! %dom
    else
       seg(I2)%dom=IUN
       seg(SEGI)%dom=IUN
    endif

    call connecIJ(bvnno,i1,5241)
    call connecIJ(I1,SEGI,52441)
    call connecIJ(Segi,I2,5141)
    call connecIJ(i2,bvnne,54141)

  endif

enddo

! The new nsegm value
nsegm = nsegm + plusseg

if(nsegm > nsegmax) then
  write(*,*) 'Problem in (1) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif

PlusSeg = izero

!********************************************************************
!****   We introduce a kneecap in the vicinity of pinning points
!***    when the length of segment is larger than LomaI/10
!********************************************************************
do Segi = 1, Nsegm

    ! bypass tests
    if (seg(Segi)%norme == izero)                              cycle
    if (segi == seg(SegI)%voiso .and. segi == seg(segi)%voise) cycle
    lomaI = Loma(seg(SEGI)%veclin)
    if (seg(Segi)%norme < INVDIX*lomaI)                        cycle

    if (seg(Segi)%voiso == nsegmax .and. seg(Segi)%surface < IUN) then

      ! ajout d'une rotule entre j et NSEGMAX pour ne pas immobiliser le segement j
      PlusSeg = PlusSeg +1
      K = nsegm + PlusSeg

      call init_seg(K)

      if(KKdebug) write(379,*)  " Topolo:On ajoute rotule",K," entre Nsegmax et ",Segi
      IBV = seg(SEGI)%veclin

      seg(K)%veclin       = nbrot(IBV,IBV,2)
      seg(K)%norme        = 0
      seg(K)%O(:)         = seg(SEGI)%O(:)
      seg(K)%voise        = SEGI
      seg(K)%vnne         = SEGI
      SEG(K)%bloquer      = seg(SEGI)%bloquer
      seg(K)%grain        = seg(SEGI)%grain
      SEG(K)%unload       = seg(SEGI)%unload
      seg(K)%surface      = IZERO
      seg(K)%VarFreePlan  = IZERO
      seg(k)%dom          = seg(SEGI)%dom
      ! mise a jours de I
      seg(SEGI)%voiso     = K
      seg(SEGI)%vnno      = NSEGMAX


    elseif (seg(Segi)%voise == nsegmax .and. seg(Segi)%surface < IUN) then

      ! on ajoute un segemnt entre I et NSEGMAX
      PlusSeg = PlusSeg +1
      K = nsegm + PlusSeg

      if(KKdebug) write(379,*)   " Topolo:On ajoute rotule",K," entre ",Segi,"  et Nsegmax "

      call init_seg(K)
      IBV = seg(SEGI)%veclin

      seg(K)%veclin     = nbrot(IBV,IBV,2)
      seg(K)%norme      = izero
      seg(K)%O(:)       = seg(SEGI)%O(:)+ seg(Segi)%norme * bveclin(:,IBV)
      seg(K)%voiso      = SEGI
      seg(K)%vnno       = SEGI
      seg(K)%bloquer    = seg(SEGI)%bloquer
      seg(K)%grain      = seg(Segi)%grain
      seg(K)%unload     = seg(SEGI)%unload
      seg(K)%surface    = IZERO
      seg(K)%VarFreePlan= IZERO
      seg(k)%dom        = seg(SEGI)%dom
      ! mise a jours de I
      seg(SEGI)%voise = K
      seg(SEGI)%vnne = Nsegmax

    endif

enddo

!********************************************************************
! ici on verifie que le caractere GD ou jonc sont toujours necessaire
do Segi = 1, Nsegm

  if (segi == seg(SegI)%voiso .and. segi == seg(Segi)%voise) cycle

  IBV   = SEG(SEGI)%VECLIN
  lomaI = loma(IBV)
  ilong = SEG(SEGI)%NORME
  bvnno = seg(Segi)%vnno
  bvnne = seg(SEGI)%vnne

  ! dans la condition suivante on essaie de resoudre le problem de mobilite autour du GD
  ! on intervient lorsqu'il s'agit d'iun segemnt mixte=coin +vis avec un voisin GD de taille 1
  if (seg(SEGI)%gd < iun .and. Ilong < 3 .and. ilong /= 0) then

    ! le voisin a l'origine est GD et de taille 1 : probelem de mobilite
    condition = (modulo(IBV,IDEUX) == 0                                                                    .and. &
                 coefSc(assoc(IBV,1), assoc(IBV,3)) == 1                                                   .and. &
                 ((seg(bvnno)%gd > izero .and. seg(bvnno)%norme == 1 .and. seg(Segi)%voiso /= bvnno) .or.        &
                  (seg(bvnne)%gd > izero .and. seg(bvnne)%norme == 1 .and. seg(segi)%voise /= bvnne))              )

    if (condition) then

      PlusSeg = PlusSeg + 2
      i2      = nsegm + PlusSeg
      i1      = i2 - 1
      seg(i1) = seg(Segi)
      seg(i2) = seg(Segi)

      if (seg(bvnno)%gd > izero) then

        if (modulo(IBV,int(nbasered/2,DPI)) == 0) then
          seg(I1)%veclin = IBV - 1
          if (modulo(IBV,int(nbasered,DPI)) == 0) then
            seg(I2)%veclin = IBV + 1 - nbasered
          else
            seg(I2)%veclin = IBV + 1
          endif
        else
          seg(I1)%veclin = IBV + 1
          seg(I2)%veclin = IBV - 1
        endif

        seggd = bvnno

      else

        if (modulo(IBV,int(nbasered/2,DPI)) == 0) then
          seg(I2)%veclin = IBV - 1
          if (modulo(IBV,int(nbasered,dpi)) == 0) then
              seg(I1)%veclin = IBV + 1 - nbasered
          else
              seg(I1)%veclin = IBV + 1
          endif
        else
          seg(I2)%veclin = IBV + 1
          seg(I1)%veclin = IBV - 1
        endif
        seggd = bvnne

      endif

      lgd = seg(seggd)%veclin

      totoBUGv1(1:3) = bveclin(1:3,LGD)
      totoBUGv2(1:3) = bveclin(1:3,assoc(IBV,1))

      if(dot_product(totoBUGv1(1:3),totoBUGv2(1:3)) > 0) then
        LGD = assoc(IBV,1)
      else
        LGD = assoc(IBV,5)
      endif

      seg(seggd)%veclin = LGD

      ! on remplace i par i1 + i2
      seg(I1)%O(:)  = seg(Segi)%O(:)
      seg(I2)%O(:)  = seg(i1)%O(:) + Ilong * bveclin(:,seg(i1)%veclin)

      ! elimination de i
      print *, "Mobilite de point GD i ==> i1,i2 ",Segi ,i1,i2
      seg(seg(Segi)%voiso)%voise = i1
      seg(seg(Segi)%voise)%voiso = i2
      seg(i1)%voise = i2
      seg(i2)%voiso = i1
      print *, " bvnno,i1 =",bvnno,i1
      call ConnecIJ(bvnno,i1,1004)
      call ConnecIJ(i1,i2,1005)
      call ConnecIJ(i2,bvnne,1006)
      seg(Segi)%norme = 0
      seg(Segi)%voiso = Segi
      seg(Segi)%voise = Segi
      out(Segi)       = .true.

      if (KKdebug) call seginfo(i1,"apres Mobilite de gd         ")

      cycle

    endif

  endif

enddo

! The new nesm value
nsegm = nsegm + PlusSeg
if(nsegm > nsegmax) then
  write(*,*) 'Problem in (2) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif
PlusSeg = izero

!*** FIN DES BOUCLES RELATIVES A LA DISCRETISATION DES SEGMENTS LONGS *******

!****************************************************************************
!******** BEGINNING OF THE JUNCTION ANALYSIS ********************************
!****************************************************************************

if(kkdebug) write(379,*)  "================ TREATMENT OF JUNCTIONS ==============="

! Initialization
RauSysjonc(1:NTSG) = zero

! We test all the segments (this loop could be optimized!)
INI: do Segi = 1, Nsegm

  ! We look only for segments in junction configuration and not to eliminate
  if (.not. seg(Segi)%jonc)                                   cycle INI
  if (SegI == seg(Segi)%voise .and. SegI == seg(Segi)%voiso ) cycle INI

  SegJ = seg(Segi)%Ijonc

  ! A junction must not be treated two times
  if (segI >= SegJ)  cycle INI

  VLi     = seg(Segi)%veclin
  VLj     = seg(segj)%veclin
  CNN     = Inn(VLi)                  ! a non zero component of vector VLi
  lomaiR  = int(loma(VLi)*quart, DPI)
  bvnno   = seg(segi)%vnno
  bvnne   = seg(segi)%vnne
  ilong   = seg(segi)%norme
  jlong   = seg(segj)%norme

  ! initialization
  gd2elimin = izero  ! This quantity must be initialised for each segment I

  if(kkdebug) then
    write(379,*)  "*******  Treatment Of junction : I =", segI , " and J = ", segJ
    call check_connec ("   Before junction Treatment  ")
    call Seginfo (segi, "segi, before junction Treatment")
    call Seginfo (segj, "segj, before junction Treatment")
  endif

  ! paranoid test : We control that the two junction segments were made at the same time
  if (seg(segi)%tjonc /= seg(segJ)%tjonc) &
    stop "Discreti - erreur : tjonc different pour I et J"

  ! Segments are shifted to avoid any problems coming from PBC
  Ei(1:3)  = (ilong * Bveclin(1:3,VLi)) / ideux
  TOi(1:3) = (seg(Segi)%O(1:3) + Ei(1:3) - modur(1:3)/2)
  Oi(1:3)  = modulo(seg(segi)%o(1:3)-toi(1:3), modur(1:3))
  Oj(1:3)  = modulo(seg(segJ)%o(1:3)-toi(1:3), modur(1:3))
  Ei(1:3)  = modulo(seg(segi)%o(1:3) +                                     &
                    ilong*bveclin(1:3,VLi)-toi(1:3), modur(1:3))
  Ej(1:3)  = modulo(seg(segJ)%o(1:3) +                                     &
                    jlong*bveclin(1:3,VLj)-toi(1:3), modur(1:3))

  ! The problem of very long junctions
  JunctLenMax = .false.

  ! Junction length is tested, if this length is close to the dimensions of the simulated volume (45%)
  ! This junction is too long and may create problems with the PBC
  ! The junction is artificially eliminated and a flag is made to alarm users
  do Iside = 1, 3

    if ( (ilong * abs(bveclin(Iside,VLi)) - Int(0.45 * float(modur(Iside)),DPI)) > 0 ) JunctLenMax = .true.
    if ( (jlong * abs(bveclin(Iside,VLj)) - Int(0.45 * float(modur(Iside)),DPI)) > 0 ) JunctLenMax = .true.

  enddo

  if (JunctLenMax) then

    print *, "A very long junction was eliminated in the simulation !"
    print *, "The simulated volume or the initial dislocation density may be too small!"

    seg(segI)%jonc = .false.
    seg(segi)%tjonc = izero
    seg(segi)%Ijonc = nsegmax

    seg(segJ)%jonc = .false.
    seg(segJ)%tjonc = izero
    seg(segJ)%Ijonc = nsegmax
    cycle INI

  endif

  ! We test if the two segments are still aligned and in contact!
  decal(:) = Oj(:) - Oi(:)
  if (decal(cnn) /= Izero) then

    if (dot_product(bvecdep(:,VLi),decal(:)) /= izero .or. dot_product(bvecnor(:,VLi),decal(:)) /= izero) then

      if (sysobst(VLi,VLj) == iun) then

        ! The two segments are collinear and for some reason the two lines did not annihilate
        ! after the step 1 of junctions formation needed some time to create line length to annihilate
        ! the junction configuration can be destroy

        ! The i an j junction property elimination
        seg(Segi)%jonc  = .false.
        seg(Segi)%tjonc = izero
        seg(Segi)%Ijonc = nsegmax
        seg(Segj)%jonc  = .false.
        seg(Segj)%tjonc = izero
        seg(Segj)%Ijonc = nsegmax

        ! The segment connectivity is systematically reconstructed
        ! since the junction direction is maybe not needed any more
        if(ilong == Izero)  then
          bvnno = seg(segi)%vnno
          bvnne = seg(segi)%vnne
          call connecij (bvnno, bvnne, 5896)
        endif

        if(jlong == Izero) then
          bvnno = seg(segj)%vnno
          bvnne = seg(segj)%vnne
          call connecij (bvnno, bvnne, 5897)
        endif

        ! We can cycle to the next segment Segi
        cycle INI

      else  ! Two junction segment not aligned are found

        write(379,*) " kk = ", kk , " I = ",Segi," J=",SEGJ
        write(379,*) " decal(:) = ", decal(:)
        write(379,*) " bvecdep(:) = ", bvecdep(:,VLi)
        call conf(Segi)
        call conf(Segj)
        write(379,*) " Two junction segments not aligned!"
        stop  " Two junction segments not aligned!"

      endif

    endif

  endif

  ! definition of Ei, Oj and Ej with respect to Oi
  RefEi = Ei(cnn) - Oi(cnn)
  RefOj = Oj(cnn) - Oi(cnn)
  RefEj = Ej(cnn) - Oi(cnn)

  ! The particular case of at least one junction segment with zero length
  if( .not. ((Ilong /= Izero) .and. (jlong /= Izero)) ) then

    ! tjonc > tjoncmax and one junction segment have zero length = We could not zip the junction
    ! the vnn of a jonction segment with zero length is a jonction segment in the same direction
    ! Hence the junction must be eliminated
    if (seg(segi)%tjonc > tjoncmax                                                                               .or. &
        (ilong == izero .and. ((seg(seg(segi)%vnno)%jonc .and. seg(segi)%veclin == seg(seg(segi)%vnno)%veclin)        &
                          .or. (seg(seg(segi)%vnne)%jonc .and. seg(segi)%veclin == seg(seg(segi)%vnne)%veclin))) .or. &
        (jlong == izero .and. ((seg(seg(segj)%vnno)%jonc .and. seg(segj)%veclin == seg(seg(segj)%vnno)%veclin)        &
                          .or. (seg(seg(segj)%vnne)%jonc .and. seg(segj)%veclin == seg(seg(segj)%vnne)%veclin)))      &
        ) then

      if (kkdebug) write(379,*) "i or j have zero length and tjonc> tjoncmax ==> junction elimination"

      ! The i an j junction property elimination
      seg(Segi)%jonc  = .false.
      seg(Segi)%tjonc = izero
      seg(Segi)%Ijonc = nsegmax
      seg(Segj)%jonc  = .false.
      seg(Segj)%tjonc = izero
      seg(Segj)%Ijonc = nsegmax

      ! The segment connectivity is systematically reconstructed
      ! since the junction direction is maybe not needed any more
      if(ilong == Izero)  then
        bvnno = seg(segi)%vnno
        bvnne = seg(segi)%vnne
        call connecij (bvnno, bvnne, 5896)
      endif

      if(jlong == Izero) then
        bvnno = seg(segj)%vnno
        bvnne = seg(segj)%vnne
        call connecij (bvnno, bvnne, 5897)
      endif

    ! if tjonc < tjoncmax, the junction segments are still in the germination stage
    ! we must wait one step more to see the junction evolution
    else

      if (kkdebug) then
        write(379,*) "i or j have zero length and tjonc < tjoncmax"
        write(379,*) " ==> we must wait one step more to see the junction evolution"
      endif

      seg(Segi)%tjonc = seg(Segi)%tjonc + IUN
      seg(SegJ)%tjonc = seg(SegJ)%tjonc + IUN

    endif

    ! Some thing has been made we can cycle to the next segment Segi
    cycle INI

  endif

  ! Initialization of important geometrical variables
  LIO = izero   ! The length of I outside of the junction section on the O side
  LIE = izero   ! ''''''''''''' I '''''''''''''''''''''''''''''''''''''' E ''''
  LJO = izero   ! ''''''''''''' J '''''''''''''''''''''''''''''''''''''' O ''''
  LJE = izero   ! ''''''''''''' J '''''''''''''''''''''''''''''''''''''' E ''''

  ! Do we still have overlapping of segment I and J to form a junction
  Recouvrement = .false.

  if(kkdebug) then
    ! write(379,'(" RefEi=", I3, "; RefOj=", I3,"; RefEj=", I3)') RefEi, RefOj, RefEj
    call Seginfo (segi, "segi, before junction reaction")
    call Seginfo (segj, "segj, before junction reaction")
  endif

  ! Calculation of the geometrical variables defining the junction configuration
  ! Different cases must be taken into account depending the Oi, Oj, Ei and Ej relative positions
  if (RefEi > izero) then

    !1- On selectionne d'abord les segments i ayant Ei a droite de Oi
    if (RefOj > izero) then

      ! Cas ou Oj est a droite de Oi : soit Oj est entre Oi et Ei soit il est a droite de Ei
      if (RefOj < RefEi) then

        ! Cas ou Oj est entre Oi et Ei ==> recouvrement
        Recouvrement = .true.

        if(RefOj < RefEj) then

          ! Cas o Ej est a droite de Oj
          ! 3 cas suivant la position de Ej par rapport a Ei : Oi-Oj-Ej-Ei ou Oi-Oj-(Ei=Ej) ou Oi-Oj-Ei-Ej
          LIO = RefOj/ bveclin(cnn,VLi)

          if(RefEj > RefEi) then
            ! Cas : Oi-Oj-Ei-Ej
            LJE = (RefEj-RefEi)/bveclin(cnn,VLj)
          else
            ! Cas : Oi-Oj-(Ei=Ej) et Oi-Oj-Ej-Ei
            LIE = (RefEi-RefEj)/bveclin(cnn,VLi)
          endif

        else

          ! Cas o Ej est a gauche de Oj
          ! 3 cas suivant la position de Ej par rapport a Oi : Oi-Ej-Oj-Ei ou (Oi=Ej)-Oj-Ei ou Ej-Oi-Oj-Ei
          LIE = (RefEi-RefOj)/bveclin(cnn,VLi)

          if(RefEj > izero) then
            ! Cas : Oi-Ej-Oj-Ei
            LIO = RefEj/bveclin(cnn,VLi)
          else
            ! Cas : (Oi=Ej)-Oj-Ei et Ej-Oi-Oj-Ei
            LJE = RefEj/ bveclin(cnn,VLj)
          endif

        endif

      else

        ! Cas o Oj est a droite de Ei ou il est confondu avec Ei
        if(RefEj < RefEi) then

          ! Cas o Ej est situe a gauche de Ei ==> recouvrement
          Recouvrement = .true.

          if(RefOj > RefEi) then

            ! Cas o Oj est a droite de Ei
            ! 3 cas suivant la position de Ej par rapport a Oi : Oi-Ej-Ei-Oj ou (Oi=Ej)-Ei-Oj ou Ej-Oi-Ei-Oj
            LJO = (RefEi-RefOj)/bveclin(cnn,VLj)

            if(RefEj > izero) then
              ! Cas : Oi-Ej-Ei-Oj
              LIO = RefEj/ bveclin(cnn,VLi)
            else
              !Cas : (Oi=Ej)-Ei-Oj et Ej-Oi-Ei-Oj
              LJE = RefEj/ bveclin(cnn,VLj)
            endif

          else

            ! Cas o Ei et Oj sont confondus
            ! 3 cas suivant la position de Ej par rapport a Oi: Oi-Ej-(Ei=Oj) ou (Oi=Ej)-(Ei=Oj) ou Ej-Oi-(Ei=Oj)
            if(RefEj > izero) then
              ! Cas : Oi-Ej-(Ei=Oj)
              LIO = RefEj/ bveclin(cnn,VLi)
            else
              ! Cas : Ej-Oi-(Ei=Oj) et (Ej=Oi)-(Ei=Oj)
              LJE = RefEj/ bveclin(cnn,VLj)
            endif

          endif

        else
          ! Oj et Ej sont situe a droite de Ei, c-a-d qu'il n'y a pas de recouvrement non ponctuel entre i et j
          Recouvrement = .false.

        endif

      endif

    else

      ! Cas o Oj est soit confondu avec Oi soit il est a gauche de Oi
      if(RefEj  > izero) then

        ! Cas o Ej est a droite de Oi ==> recouvrement
        Recouvrement = .true.

        if(RefOj < izero)then

          ! Cas ou Oj est a gauche de Oi
          ! 3 Cas suivant la position de Ej par rapport a Ei : Oj-Oi-Ej-Ei ou Oj-Oi-(Ei=Ej) ou Oj-Oi-Ei-Ej
          LJO = -RefOj/ bveclin(cnn,VLj)

          if(RefEj > RefEi) then
            ! Cas :  Oj-Oi-Ei-Ej
            LJE = (RefEj-RefEi)/ bveclin(cnn,VLj)
          else
            ! Cas : Oj-Oi-Ej-Ei et Oj-Oi-(Ei=Ej)
            LIE = (RefEi-RefEj)/ bveclin(cnn,VLi)
          endif

        else

          ! Cas o Oj est confondu avec Oi
          ! 3 Cas suivant la position de Ej par rapport a Ei : (Oj=Oi)-Ej-Ei ou (Oj=Oi)-(Ei=Ej) ou (Oj=Oi)-Ei-Ej
          if(RefEj > RefEi) then
            ! Cas : (Oi=Oj)-Ei-Ej
            LJE = (RefEj-RefEi)/ bveclin(cnn,VLj)
          else
            ! Cas : (Oi=Oj)-(Ei=Ej) et (Oi=Oj)-Ej-Ei
            LIE = (RefEi-RefEj)/bveclin(cnn,VLi)
          endif

        endif

      else

          ! Oj et Ej sont situe a gauche du Oi, c-a-d qu'il n'y a pas de recouvrement non ponctuel entre i et j
          Recouvrement = .false.

      endif

    endif

  else

    !2- On selectionne ensuite les segments i ayant Ei a gauche de Oi
    if(RefOj < izero) then

     ! Cas o Oj est situe a gauche de Oi: soit Oj est confondu avec Ei, soit il est a droite ou a gauche de Ei
      if(RefOj > RefEi) then

        ! Cas o Oj est a droite de Ei, c-a-d qu'il appartient au segment i ==> recouvrement
        Recouvrement = .true.

        if(RefOj < RefEj) then

          ! Cas o Ej est a droite de Oj
          ! 3 Cas suivant la position de Ej par rapport a Oi: Ei-Oj-Ej-Oi ou Ei-Oj-(Oi=Ej) ou Ei-Oj-Oi-Ej
          LIE = (RefEi-RefOj)/bveclin(cnn,VLi)

          if(RefEj > izero) then
            ! Cas : Ei-Oj-Oi-Ej
            LJE = RefEj/bveclin(cnn,VLj)
          else
            ! Cas :  Ei-Oj-Ej-Oi et Ei-Oj-(Oi=Ej)
            LIO = RefEj/bveclin(cnn,VLi)
          endif

        else

          ! Cas o Ej est a gauche de Oj
          ! 3 Cas suivant la position de Ej par rapport a Ei: Ej-Ei-Oj-Oi ou (Ej=Ei)-Oj-Oi ou Ei-Ej-Oj-Oi
          LIO = RefOj/bveclin(cnn,VLi)

          if(RefEj > RefEi) then
            ! Cas : Ei-Ej-Oj-Oi
            LIE = (RefEi-RefEj)/bveclin(cnn,VLi)
          else
            ! Cas : Ej-Ei-Oj-Oi et (Ej=Ei)-Oj-Oi
            LJE = (RefEj-RefEi)/ bveclin(cnn,VLj)

          endif

        endif

      else

        ! Cas o Oj est situe a gauche de Ei ou il est confondu avec Ei
        if(RefEj > RefEi) then

         ! Cas o Ej est a droite de Ei ==> recouvrement
          Recouvrement = .true.

          if(RefOj < RefEi) then

            ! Cas o Oj est a gauche de Ei
            ! 3 cas suivant la position de Ej par rapport a Oi: Oj-Ei-Ej-Oi ou Oj-Ei-(Ej=Oi) ou Oj-Ei-Oi-Ej
            LJO = (RefEi-RefOj)/bveclin(cnn,VLj)

            if(RefEj > izero) then
              ! Cas : Oj-Ei-Oi-Ej
              LJE = RefEj/ bveclin(cnn,VLj)
            else
              ! Cas : Oj-Ei-Ej-Oi et Oj-Ei-(Ej=Oi)
              LIO = RefEj/bveclin(cnn,VLi)
            endif

          else

            ! Cas o Oj est confondu avec Ei
            ! 3 Cas suivant la position de Ej par rapport a Oi: (Oj=Ei)-Ej-Oi ou (Oj=Ei)-(Ej=Oi) ou (Oj=Ei)-Oi-Ej
            if(RefEj > izero) then
              ! Cas : (Oj=Ei)-Oi-Ej
              LJE = RefEj/ bveclin(cnn,VLj)
            else
              ! Cas : (Oj=Ei)-Ej-Oi et (Oj=Ei)-(Ej=Oi)
              LIO = RefEj/ bveclin(cnn,VLi)
            endif

          endif

        else

          ! Oj et Ej sont situe a gauche du Ei, c-ad qu'il n'y a pas de recouvrement non ponctuel entre i et j
          Recouvrement = .false.

        endif

      endif

    else

      ! Cas o Oj est a droite de Oi ou il est confondu avec Oi
      if(RefEj < izero) then

        ! Cas o Ej est situe a gauche de Oi ==> recouvrement
        Recouvrement = .true.

        if(RefOj /= izero)then
          ! Oj est non confondu avec Oi, c-ad que Oj est a droite de Oi
          ! 3 cas suivant la position de Ej par rapport a Ei : Ei-Ej-Oi-Oj ou (Ei=Ej)-Oi-Oj ou Ej-Ei-Oi-Oj
          LJO= -RefOj/bveclin(cnn,VLj)

          if(RefEj > RefEi) then
            ! Cas :  Ei-Ej-Oi-Oj
            LIE = (RefEi-RefEj) / bveclin(cnn,VLi)
          else
            ! Cas : (Ei=Ej)-Oi-Oj et Ej-Ei-Oi-Oj
            LJE = (RefEj-RefEi)/ bveclin(cnn,VLj)
          endif

        else

          ! Oj est confondu avec Oi
          ! 3 cas suivant la position de Ej par rapport a Ei : Ei-Ej-(Oi=Oj) ou (Ei=Ej)-(Oi=Oj) ou Ej-Ei-(Oi=Oj)
          if(RefEj > RefEi) then
            ! Cas : Ei-Ej-(Oi=Oj)
            LIE = (RefEi-RefEj)/ bveclin(cnn,VLi)
          else
            ! Cas : (Ei=Ej)-(Oi=Oj) et Ej-Ei-(Oi=Oj)
            LJE = (RefEj-RefEi)/bveclin(cnn,VLj)
          endif

        endif

      else

        ! Oj et Ej sont situe a droite de Oi, c-a-d qu'il n'y a pas de recouvrement non ponctuel entre i et j
        Recouvrement = .false.

      endif

    endif

  endif

  if (kkdebug) write(379,'(" Recouv=",L2,"; LomaiR=", I3, "; LIO=",                         &
                          & I3,"; LJO=", I3, "; LIE=", I3, "; LJE=", I3,", Connec=",I3 )')  &
                          Recouvrement, LomaiR, LIO, LJO, LIE, LJE,  sysconnec(VLi,VLj)

  !********************************************************************************
  !* Junction Annihilation must not be considered for segments touching a surface,
  !* it must proceed in that case from the surface elimination
  !********************************************************************************
  if (Recouvrement .and. seg(segI)%surface==izero .and. seg(segj)%surface==izero)  then

    itemp = Ilong - Lio - Lie

    if(Itemp /= Jlong - Ljo - Lje ) then
      print *, "Two segments building one junction are not overlapping !",segI,segJ,KK
      stop
    endif

    rausysjonc(syseg(Vli)) = rausysjonc(syseg(Vli)) + itemp * normlin(VLI)
    rausysjonc(syseg(VlJ)) = rausysjonc(syseg(VlJ)) + itemp * normlin(VLJ)

    ! Cas d'annihilation : recouvrement et meme vecteur de Burgers
    if(sysconnec(VLi,VLj) == ideux) then

      ! Determination of hands of old junction for the cutting of segments : iO side
      HandIo    = seg(segI)%vnno
      compteur  = izero

      do while (seg(handIo)%norme == izero .and. compteur < 100)
        HandIo    = seg(handIo)%vnno
        compteur  = compteur + iun
      enddo

      if (compteur > 99) then
        print * , " Discreti: erreur de voisinage O de i",segi
        stop
      endif

      if (HandIo == NSEGMAX) HandIo = seg(segI)%voiso

      ! Determination of hands of old junction for the cutting of segments : iE side
      HandIe = seg(segI)%vnne
      compteur  = izero

      do while (seg(handIe)%norme == izero .and. compteur < 100)
        HandIe    = seg(handIe)%vnne
        compteur  = compteur + iun
      enddo

      if (compteur > 99) then
        print * , " Discreti: erreur de voisinage E de i",segi
        stop
      endif

      if (HandIe == NSEGMAX) HandIe = seg(segI)%voise

      ! Determination of hands of old junction for the cutting of segments : jO side
      HandJo = seg(segJ)%vnno
      compteur  = izero

      do while (seg(handJo)%norme == izero .and. compteur < 100)
        HandJo    = seg(handJo)%vnno
        compteur  = compteur + iun
      enddo

      if (compteur > 99) then
          print * , " Discreti: erreur de voisinage O de j",segJ
          stop
      endif

      if (HandJo == NSEGMAX) HandJo = seg(segJ)%voiso

      ! Determination of hands of old junction for the cutting of segments : jE side
      HandJe = seg(segJ)%vnne
      compteur  = izero

      do while (seg(handJe)%norme == izero .and. compteur < 100)
          HandJe    = seg(handJe)%vnne
          compteur  = compteur + iun
      enddo

      if (compteur > 99) then
          print * , " Discreti: erreur de voisinage E de j",segJ
          stop
      endif

      if (HandJe == NSEGMAX) HandJe = seg(segJ)%voise

      ! A special case
      if ( HandJo == segi .or. HandIo == segj ) then
          if (kkdebug) write(379,*)  "Cas ignore !!! - Double recouvrement "
          cycle INI
      endif

      !special case to avoid small loop touching free surfaces!!!!
      if (seg(HandIe)%surface > IZERO .and. seg(handJo)%surface > IZERO) cycle
      if (seg(HandIo)%surface > IZERO .and. seg(handJe)%surface > IZERO) cycle

      ! Case of complete annihilation of I
      if (LIO == Izero .and. LIe == izero) then

        ! Suppresion of all segments between hands of junction
        k = seg(HandIo)%voise

        if (k == NSEGMAX) k = seg(segi)%voiso !the first neighbour is a pinning point

        compteur  = izero
        do while (K /= HandIe .and. compteur < 100)

          if (kkdebug) write(379,*) " Discreti - Annihilation of i : elimination of ",K,"; compteur= ",compteur

          if (seg(k)%jonc) then

            itemp             = seg(k)%Ijonc
            seg(K)%jonc       = .false.
            seg(K)%tjonc      = izero
            seg(K)%Ijonc      = nsegmax
            seg(itemp)%jonc   = .false.
            seg(itemp)%tjonc  = izero
            seg(itemp)%Ijonc  = nsegmax

            GD2made = .false.
            if(seg(itemp)%norme == izero) call voisinage(itemp,-87)

          endif

          itemp         = seg(k)%voise
          seg(k)%norme  = izero
          if (seg(k)%gd > iun) then
            gd2elimin = gd2elimin + 1    ! We keep the info about eliminated gd2
            if (kkdebug) write(379,*) "A GD2 segment will be eliminated in loop i, gd2elimin =",gd2elimin
          endif
          seg(k)%voiso  = k
          seg(k)%voise  = k
          out(k)        = .true.
          k             = itemp
          compteur      = compteur + iun

        enddo

        if (compteur > 99) then
          print * , " Discreti: error 1 ", segI, segj, kk
          stop
        endif

      endif

      ! case of complete annihilation of J
      if(LJO == Izero .and. LJe == izero) then

        k = seg(HandJo)%voise

        if (k == NSEGMAX) k = seg(Segj)%voiso !the first neighbour is a pinning point

        compteur  = izero
        do while (K /= HandJe .and. compteur < 100)

          if(kkdebug) write(379,*)  " Discreti - Annihilation of j: elimination of ",K,"; compteur= ",compteur

          if (seg(k)%jonc) then

            itemp             = seg(k)%Ijonc
            seg(K)%jonc       = .false.
            seg(K)%tjonc      = izero
            seg(K)%Ijonc      = nsegmax
            seg(itemp)%jonc   = .false.
            seg(itemp)%tjonc  = izero
            seg(itemp)%Ijonc  = nsegmax

            GD2made = .false.
            if (seg(itemp)%norme == izero) call voisinage(itemp,-88)

          endif

          itemp         = seg(k)%voise
          seg(k)%norme  = izero
          if (seg(k)%gd > iun) then
            gd2elimin = gd2elimin + 1    ! We keep the info about eliminated gd2
            if (kkdebug) write(379,*) "A GD2 segment will be eliminated in loop j, gd2elimin =",gd2elimin
          endif
          seg(k)%voiso  = k
          seg(k)%voise  = k
          out(k)        = .true.
          k             = itemp
          compteur      = compteur + iun

        enddo

        if (compteur > 99) then
          print * , " Discreti: error 2 ", segI, segj, kk
          stop
        endif

      endif

      ! Starting re-connection of hands
      i1 = seg(segi)%voiso
      i2 = seg(segi)%voise
      j1 = seg(segj)%voiso
      j2 = seg(segj)%voise

      ! case 1 : I and J disappear
      if(LIO == Izero .and. LIe == izero .and. LJO == Izero .and. LJe == izero) then

        if(kkdebug) write(379,*) " Suppression of i and j ...................................."

        if(handio /= nsegmax) seg(HandIo)%voise = HandJe
        if(handje /= nsegmax) seg(HandJe)%voiso = HandIo

        if (gd2elimin > izero) then
          gd2made = .false.
          call connecIJ(HandIo, HandJe, -354)
        else
          call connecIJ(HandIo, HandJe, 354)
        endif

        if(handjo /= nsegmax) seg(HandJo)%voise = HandIe
        if(HandIe /= nsegmax) seg(HandIe)%voiso = HandJo

        if (gd2elimin > izero) then
          gd2made = .false.
          call connecIJ(HandJo, HandIe, -324)
        else
          call connecIJ(HandJo, HandIe, 324)
        endif

        if(kkdebug) then
          call Seginfo (HandIo, "HandIo, after suppression of i and j  ")
          call Seginfo (HandIe, "HandIe, after suppression of i and j  ")
          call Seginfo (HandJo, "HandJo, after suppression of i and j  ")
          call Seginfo (HandJe, "HandJe, after suppression of i and j  ")
        endif

      ! case 2-  Only I disappears
      elseif (LIO == Izero .and. LIE == izero) then

        if (kkdebug) write(379,*)  " Suppression of i ........................................"

        ! case 2.1-  LJO /= 0 , Configurations: Ej-Oi-Ei-Oj and (Ej=Oi)-Ei-Oj
        if (LJO /= izero) then

          seg(segJ)%norme = LJO
          seg(segJ)%voise = HandIe

          if (Handie /= nsegmax) seg(Handie)%voiso = SegJ

          if (gd2elimin > izero .or.            &
             seg(seg(segi)%vnne)%gd > iun .or.  &
             seg(seg(segj)%vnno)%gd > iun .or.  &
             seg(seg(segj)%vnne)%gd > iun) then
             gd2made = .false.
            call connecIJ(SegJ, HandIe, -3244)
          else
            call connecIJ(SegJ, HandIe, 3244)
          endif

          if (HandJe == NSEGMAX .and. HandIo == NSEGMAX) then

            seg(j2)%voise       = j2
            seg(j2)%voiso       = j2
            out(j2)             = .true.
            seg(j2)%surface     = IZERO
            seg(j2)%VarFreePlan = IZERO

          ! case 2.1.1-  test on the end of J - Configuration:  Ej-Oi-Ei-Oj
          elseif (LJE /= izero) then

            seg(segi)%veclin  = VLj
            seg(segi)%norme   = LJE
            seg(segi)%voiso   = HandIo
            seg(segi)%voise   = j2
            out(segi)         = .false.

            if (HandIo /= nsegmax)  seg(HandIo)%voise = segi
            if (J2/= nsegmax)       seg(j2)%voiso = segi

            if (gd2elimin > izero .or.            &
               seg(seg(segi)%vnne)%gd > iun .or.  &
               seg(seg(segj)%vnno)%gd > iun) then
              gd2made = .false.
              call connecIJ(HandIo, segi, -78)
              gd2made = .false.
              call connecIJ(segi,HandJe, -44)
            else
              call connecIJ(HandIo, segi, 78)
              call connecIJ(segi,HandJe, 44)
            endif

            if(kkdebug) write(379,*)  " I Suppression - 2.1.1 ..............................."

          ! case 2.1.2-  Configuration: (Ej=Oi)-Ei-Oj
          else

            if (J2 /= nsegmax)      seg(j2)%voiso = handIo
            if (HandIo /= nsegmax)  seg(HandIo)%voise = j2

            if (gd2elimin > izero .or.            &
               seg(seg(segi)%vnne)%gd > iun .or.  &
               seg(seg(segj)%vnno)%gd > iun) then
               gd2made = .false.
              call connecIJ(HandIo, HandJe, -124)
            else
              call connecIJ(HandIo, HandJe, 124)
            endif

            if(kkdebug) write(379,*)  " I Suppression - 2.1.2 ..............................."

          endif

        ! case 2.2-  LJo = 0 , Configuration: Ej-Oi-(Ei=Oj)
        else

          if (J1 /= nsegmax) seg(j1)%voise = HandIe

          if (handIe /= nsegmax) seg(HandIe)%voiso = j1

          if (gd2elimin > izero .or.              &
              seg(seg(segi)%vnno)%gd > iun .or.   &
              seg(seg(segj)%vnne)%gd > iun) then
            gd2made = .false.
            call connecIJ(HandJo, HandIe, -14425)
          else
            call connecIJ(HandJo, HandIe, 14425)
          endif

          seg(segJ)%O(:) = seg(segj)%O(:) + (seg(segj)%norme - LJE) * Bveclin(:,VLj)

          if (Nbcvxdom > IUN)  then
            call assign_domain(segj,11000)
          else
            seg(segj)%dom=iun
          endif
          seg(segJ)%Norme = LJE

          if( HandIo /= nsegmax) seg(HandIo)%voise = SegJ

          seg(SegJ)%voiso = HandIo

          if (gd2elimin > izero) then
            gd2made = .false.
            call connecIJ(HandIo, SegJ, -14128)
          else
            call connecIJ(HandIo, SegJ, 14128)
          endif

          if(kkdebug) write(379,*)  " I Suppression - 2.2 .............................."

        endif

        if(kkdebug) then
          call Seginfo (HandIo, "HandIo, after suppression of i")
          call Seginfo (HandIe, "HandIe, after suppression of i")
          call Seginfo (HandJo, "HandJo, after suppression of i")
          call Seginfo (HandJe, "HandJe, after suppression of i")
        endif

      ! case 3-  Only J disappears
      elseif (LJO == Izero .and. LJe == izero) then

        if(kkdebug) write(379,*)  " Suppression of J ...................................."

        ! case 3.1-  LIO /= 0 , Configurations: Oi-Ej-Oj-Ei and Oi-Ej-(Oj=Ei)
        if (LIo /= izero) then

          seg(segi)%norme = LIO
          seg(segi)%voise = HandJe
          out(segi)       = .false.

          if (handJe /= nsegmax) seg(HandJe)%voiso = Segi

          if (gd2elimin > izero .or.              &
              seg(seg(segi)%vnno)%gd > iun .or.   &
              seg(seg(segi)%vnne)%gd > iun .or.   &
              seg(seg(segj)%vnne)%gd > iun) then
            gd2made = .false.
            call connecIJ(Segi, HandJe, -3245)
          else
            call connecIJ(Segi, HandJe, 3245)
          endif

          ! case 3.1.1-  test on the end of I - Configuration: Oi-Ej-Oj-Ei
          if (LIE /= izero) then

            seg(segj)%veclin = VLi
            seg(segj)%norme = LIE
            seg(segj)%voiso  = HandJo
            seg(segj)%voise  = i2
            out(segj)        = .false.

            if (handJo /= nsegmax)  seg(HandJo)%voise = segj
            if (i2 /= nsegmax)      seg(i2)%voiso = segj

            if (gd2elimin > izero .or.              &
                seg(seg(segi)%vnno)%gd > iun .or.   &
                seg(seg(segj)%vnne)%gd > iun) then
              gd2made = .false.
              call connecIJ(HandJo, segj, -79)
              gd2made = .false.
              call connecIJ(segj, HandIe, -4400)
            else
              call connecIJ(HandJo, segj, 79)
              call connecIJ(segj, HandIe, 4400)
            endif

            if(kkdebug) write(379,*)  " J Suppression - 3.1.1  .............................."

          ! case 3.1.2- Configuration: Oi-Ej-(Oj=Ei)
          else

            if (I2 /= nsegmax) seg(i2)%voiso = Handjo

            if (handJo /= nsegmax) seg(handJo)%voise = i2

            if (gd2elimin > izero .or.              &
                seg(seg(segi)%vnno)%gd > iun .or.   &
                seg(seg(segj)%vnne)%gd > iun) then
              gd2made = .false.
              call connecIJ(Handjo, handie, -127)
            else
              call connecIJ(Handjo, handie, 127)
            endif

            if(kkdebug) write(379,*)  " J Suppression - 3.1.2  .............................."

          endif

        ! case 3.2-  LIO = 0 , Configuration : (Oi=Ej)-Oj-Ei
        else

          if (I1 /= nsegmax) seg(i1)%voise = HandJe

          if (HandJe /= nsegmax) seg(HandJe)%voiso = i1

          if (gd2elimin > izero .or.            &
             seg(seg(segi)%vnne)%gd > iun .or.  &
             seg(seg(segj)%vnno)%gd > iun) then
            gd2made = .false.
            call connecIJ(HandIo, HandJe, -144)
          else
            call connecIJ(HandIo, HandJe, 144)
          endif

          seg(segi)%O(:)  =  seg(segj)%O(:)
          seg(segi)%dom   =  seg(segj)%dom
          seg(segi)%Norme = LiE

          if(handJo /= nsegmax) seg(HandJo)%voise = Segi

          seg(Segi)%voiso = HandJo
          out(segi)       = .false.

          if (gd2elimin > izero .or.                      &
             seg(seg(segi)%vnne)%gd > iun .or.            &
             seg(seg(seg(segi)%vnne)%vnne)%gd > iun .or.  &
             seg(seg(seg(segj)%vnne)%vnne)%gd > iun .or.  &
             seg(seg(segj)%vnno)%gd > iun) then
            gd2made = .false.
            call connecIJ(HandJo, Segi, -141)
          else
            call connecIJ(HandJo, Segi, 141)
          endif

          if(kkdebug) write(379,*)  " J Suppression - 3.2 .............................."

        endif

        if(kkdebug) then
          call Seginfo (HandIo, "HandIo, after suppression of i")
          call Seginfo (HandIe, "HandIe, after suppression of i")
          call Seginfo (HandJo, "HandJo, after suppression of i")
          call Seginfo (HandJe, "HandJe, after suppression of i")
        endif

      ! case 4-  Neither of I nor J disappear
      else

        ! case where LIO and LJO /= 0 and LIE = LJE = 0 - Configuration: Oi-Ej-Ei-Oj
        if(Lio /= izero) then

          if (kkdebug) write(379,*)  " Discreti, shortening of I and J  "

          seg(segi)%norme = LIO
          seg(segJ)%norme = LJO
          seg(segI)%voise = j2
          seg(segJ)%voise = i2
          out(segi)       = .false.

          if (J2 /= nsegmax) seg(j2)%voiso = segI
          if (I2 /= nsegmax) seg(i2)%voiso = segJ

          if (gd2elimin > izero .or.            &
             seg(seg(segi)%vnno)%gd > iun .or.  &
             seg(seg(segi)%vnne)%gd > iun .or.  &
             seg(seg(segj)%vnno)%gd > iun .or.  &
             seg(seg(segj)%vnne)%gd > iun) then
            gd2made = .false.
            call connecIJ(SegI, handJe, -122)
            gd2made = .false.
            call connecIJ(SegJ, handIe, -222)
          else
            call connecIJ(SegI, handJe, 122)
            call connecIJ(SegJ, handIe, 222)
          endif

          if (kkdebug) write(379,*)  "faute 2 .............................."

        ! case where LIE and LJE /= 0 and LIO = LJO = 0 - Configuration: Ej-Oi-Oj-Ei
        else

          if(kkdebug) write(379,*)  " Discreti, modify begining of I and J  ",segi, segj

          k = seg(segi)%vnno
          seg(segI)%O(:)  = seg(segj)%O(:)
          seg(segI)%dom   = seg(segj)%dom
          seg(segI)%norme = LIe
          seg(segI)%voiso = j1
          out(segI)       = .false.

          if( J1/= nsegmax)  seg(j1)%voise = segI

          if (gd2elimin > izero .or.              &
              seg(seg(segi)%vnno)%gd > iun .or.   &
              seg(seg(segj)%vnne)%gd > iun) then
            gd2made = .false.
            call connecIJ(HandJO, segI, -268)
          else
            if (seg(HandIo)%gd > iun) then
              gd2made = .false.
              call connecIJ(HandJO, segI, -267)   ! made of a GD2 must be forced here
            else
              call connecIJ(HandJO, segI, 267)
            endif
          endif

          seg(segJ)%O(:) =  seg(segj)%O(:) + ((seg(segj)%norme-LJE) * Bveclin(:,VLj))

          if (Nbcvxdom > IUN) then
            call assign_domain(segj,11001)
          else
            seg(segj)%dom = IUN
          endif
          seg(segJ)%norme = LJe
          seg(segJ)%voiso = I1

          if (I1 /= nsegmax) seg(i1)%voise = segJ

          if (gd2elimin > izero .or.            &
             seg(seg(segi)%vnne)%gd > iun .or.  &
             seg(seg(segj)%vnno)%gd > iun) then
            gd2made = .false.
            call connecIJ(k,SegJ, -133)
          else
            if (seg(HandJo)%gd > iun .or.       &
                seg(seg(HandJo)%vnno)%gd > iun)  then
              gd2made = .false.
              call connecIJ(k,SegJ, -132)   ! made of a GD2 must be forced here
            else
              call connecIJ(k,SegJ, 132)
            endif
          endif

          if(kkdebug) write(379,*)  "faute 4 .............................."

        endif

        if(kkdebug) then
          call Seginfo (segi, "segi, after annhilation - Neither of I nor J disappeared")
          call Seginfo (segj, "segj, after annhilation - Neither of I nor J disappeared")
          write(379,*)  "Now we look for the possible formation of a microloop after annihilation"
        endif

        !**********************
        ! Now it is a good idea to test if we did not form a microloop attached to segi which could be eliminated
        closed_loop = .false.
        ListA(1:20) = -iun
        NextSeg     = seg(segi)%voise
        Loop_length = zero

        do iloop = 1, 20

          ListA(iloop) = NextSeg      ! The index of segments listed in this loop is saved

          if (NextSeg == Nsegmax) exit      ! we touch a pining point, hence this is not a loop

          Loop_Length = Loop_Length + seg(NextSeg)%norme * NormLin(seg(NextSeg)%veclin)

          if (NextSeg == segi) then
            closed_loop = .true.
            exit
          endif

          NextSeg = seg(NextSeg)%voise

        enddo

        if (kkdebug) write(379,*) "The microloop test is made on the segi side made of ",  &
                                   iloop,closed_loop,int(loop_Length),microboucle

        ! A small loop can be eliminated
        if (closed_loop .and. int(loop_Length) < microboucle) then

          if(kkdebug) write(379,*)  "A microloop is found on the segi side made of ", iloop

          compteur = iloop          ! The number of segments to eliminate
          NextSeg  = segi

          do iloop = 1, compteur

            ! A segment in the loop take part to a junction
            if (seg(NextSeg)%jonc) then

              k = SEG(NextSeg)%ijonc

              ! We must test if k is a segment part of the micro-loop we want to eliminate
              ListB(1:20) = k
              if (any(ListA == ListB)) then

                seg(k)%JONC = .false.         ! Nothing to do since the segment will be eliminated
                seg(k)%IJONC = nsegmax
                seg(k)%tJONC = 0

              else          ! k is not part of the micro-loop

                !We liberate the complementary junction segment k
                if (seg(k)%norme /= 0) then

                  seg(k)%JONC = .false.       !simple case
                  seg(k)%IJONC = nsegmax
                  seg(k)%tJONC = 0

                else

                  call stopjonc(k,55506)

                  ! Attention, elimination of k can be requested by stopjonk
                  if(seg(k)%voiso == k .and. seg(k)%voise == k) then
                    out(k)=.true.
                  endif

                endif

              endif

              !We liberate the NextSeg junction segment
              seg(NextSeg)%jonc = .false.
              seg(NextSeg)%ijonc = nsegmax
              seg(NextSeg)%tjonc = 0

            endif

            if (seg(NextSeg)%gd > izero) then
!              if (seg(NextSeg)%gd == iun)   oldntgd1 = oldntgd1 + 1
              if (seg(NextSeg)%gd == ideux) then
!                oldntgd2 = oldntgd2 + 1
                if (kkdebug) write(379,*) 'In discreti - 1, oldntgd2 after icrement =', oldntgd2
              endif
            endif

            ! Segment part of the loop are eliminated
            k = seg(NextSeg)%voise
            seg(NextSeg)%voiso = NextSeg
            seg(NextSeg)%voise = NextSeg
            out(NextSeg)       = .true.

            if (kkdebug) write(379,*) "Segment part of the loop segi are eliminated",  &
                                      out(NextSeg), seg(NextSeg)%voiso, seg(NextSeg)%voise

            NextSeg = k         ! The next segment to eliminate

          enddo

        endif                   ! End of the segi microloop elimination

        !**********************
        ! Now it is a good idea to test if we did not form a microloop attached to segj which could be eliminated

        ! We first check that Segj is not a segment eliminated in the above Segi micro-loop test
        if ((segJ /= seg(segj)%voise .and. segJ /= seg(segj)%voiso) .and. .not. out(segJ)) then

          closed_loop = .false.
          ListA(1:20) = -iun
          NextSeg     = seg(segj)%voise
          Loop_length = zero

          do iloop = 1, 20

            ListA(iloop) = NextSeg      ! The index of segments listed in this loop is saved

            if (NextSeg == Nsegmax) exit      ! we touch a pining point, hence this is not a loop

            Loop_Length = Loop_Length + seg(NextSeg)%norme * NormLin(seg(NextSeg)%veclin)

            if (NextSeg == segj) then
              closed_loop = .true.
              exit
            endif

            NextSeg = seg(NextSeg)%voise

          enddo

          if (kkdebug) write(379,*) "The microloop test is made on the segj side made of ", &
                                      iloop,closed_loop,int(loop_Length),microboucle

          ! A small loop should be eliminated
          if (closed_loop .and. int(loop_Length) < microboucle) then

            if(kkdebug) write(379,*)  "A microloop is found on the segj side made of ", iloop

            compteur = iloop          ! The number of segments to eliminate
            NextSeg  = segj

            do iloop = 1, compteur

              ! A segment in the loop take part to a junction
              if (seg(NextSeg)%jonc) then

                k = SEG(NextSeg)%ijonc

                ! We must test if k is a segment part of the micro-loop we want to eliminate
                ListB(1:20) = k
                if (any(ListA == ListB)) then

                  seg(k)%JONC = .false.       ! Nothing to do since the segment will be eliminated
                  seg(k)%IJONC = nsegmax
                  seg(k)%tJONC = 0

                else          ! k is not part of the micro-loop

                  !We liberate the complementary junction segment k
                  if (seg(k)%norme /= 0) then

                    seg(k)%JONC = .false.       !simple case
                    seg(k)%IJONC = nsegmax
                    seg(k)%tJONC = 0

                  else

                    call stopjonc(k,55996)

                    ! Attention, elimination of k can be requested by stopjonk
                    if(seg(k)%voiso == k .and. seg(k)%voise == k) out(k)=.true.

                  endif

                endif

                !We liberate the NextSeg junction segment
                seg(NextSeg)%jonc = .false.
                seg(NextSeg)%ijonc = nsegmax
                seg(NextSeg)%tjonc = 0

              endif

              if (seg(NextSeg)%gd > izero) then
!                if (seg(NextSeg)%gd == iun)   oldntgd1 = oldntgd1 + 1
                if (seg(NextSeg)%gd == ideux) then
!                  oldntgd2 = oldntgd2 + 1
                  if (kkdebug) write(379,*) 'In discreti - 2, oldntgd2 after icrement =', oldntgd2
                endif
              endif

              ! Segment part of the loop are eliminated
              k = seg(NextSeg)%voise

              seg(NextSeg)%voiso = NextSeg
              seg(NextSeg)%voise = NextSeg
              out(NextSeg)=.true.

              if (kkdebug) write(379,*) "Segment part of the loop segj are eliminated",       &
                                        out(NextSeg), seg(NextSeg)%voiso, seg(NextSeg)%voise

              NextSeg = k         ! The next segment to eliminate

            enddo

          endif                   ! End of the segj microloop elimination

        endif

      endif

    ! Traitement des cas ou on a recouvrement mais on a pas le meme vecteur de Burgers
    ! ==> On incremente tjonc et on discretise LIO, LJO, LIE et LJE s'ils sont plus grands que lomaiR
    else

      seg(segi)%tjonc = seg(segi)%tjonc + iun
      seg(segj)%tjonc = seg(segi)%tjonc

      if (LIO >  LomaiR) then
          if (kkdebug) write(379,*)  "-------------------------  LIO >  LomaiR"
          Plusseg = Plusseg + iun
          K = nsegm + Plusseg
          seg(k)        = seg(segi)
          seg(k)%O(:)   = seg(segi)%O(:)
          seg(K)%norme  = LIO
          seg(k)%veclin = VLi
          seg(K)%jonc   = .false.
          seg(K)%tjonc  = izero
          seg(K)%Ijonc  = nsegmax
          if (seg(segi)%surface==IUN) then
            seg(k)%surface        = IUN
            seg(k)%VarFreePlan    = seg(segi)%VarFreePlan
            seg(segi)%surface     = IZERO
            seg(segi)%VarFreePlan = IZERO
          else
            seg(k)%surface      = IZERO
            seg(k)%VarFreePlan  = IZERO
          endif
          seg(segi)%O(:)  = seg(segi)%O(:) + LIO * Bveclin(:,VLi)
          seg(segi)%norme = ilong - LIO
          seg(segi)%voiso = K
          seg(K)%voise    = segi
          seg(segi)%vnno  = K
          seg(K)%vnne     = segi
          if (Nbcvxdom > IUN) then
            call assign_domain(k,11002)
            call assign_domain(segi,11003)
          else
            seg(k)%dom    = IUN
            seg(segi)%dom = IUN
          endif
          itemp = seg(K)%voiso
          if(itemp /= nsegmax) then
            seg(itemp)%voise  = K
            seg(itemp)%vnne   = K
          endif
          call connecIJ(k,segi,241)
          if (kkdebug) write(379,*) "LIO ---------------------------------"
          call voisinage(K,-242)
          if (sysconnec(VLi,seg(seg(segi)%voise)%veclin) /= 1 .and.  &
             (seg(segi)%voise/=nsegmax) .and. seg(segi)%gd > izero) then
            seg(K)%gd = seg(segi)%gd
          else
            seg(K)%gd = izero
          endif
          seg(K)%DISEG    = .false.
          seg(K)%BLOQUER  = .false.
          seg(K)%unload   = .false.
          seg(k)%Nsysdev  = izero
          seg(k)%probadev = izero
          seg(k)%taudev   = izero
          if (kkdebug) call seginfo(SegI, " Apres traitement de LIO                  ")

      endif

      if ( LIE > lomaiR) then
          if (kkdebug) write(379,*)  "-------------------------      LIE > lomaiR     "
          Plusseg = Plusseg + iun
          K = nsegm + Plusseg
          seg(K)        = seg(segi)
          seg(K)%norme  = LIE
          seg(K)%veclin = VLi
          seg(K)%O(:)   = seg(segi)%O(1:3) + (seg(segi)%norme - LIE) * Bveclin(:,VLi)
          if (seg(segi)%surface == IDEUX) then
            seg(k)%surface        = IDEUX
            seg(k)%VarFreePlan    = seg(segi)%VarFreePlan
            seg(segi)%surface     = IZERO
            seg(segi)%VarFreePlan = IZERO
          else
            seg(k)%surface      = IZERO
            seg(k)%VarFreePlan  = IZERO
          endif
          seg(K)%jonc     = .false.
          seg(K)%tjonc    = izero
          seg(K)%Ijonc    = nsegmax
          seg(segi)%voise = K
          seg(K)%voiso    = segi
          seg(segi)%vnne  = K
          seg(K)%vnno     = segi
          seg(segi)%norme = seg(segi)%norme - LIE

          itemp = seg(K)%voise
          if(itemp /= nsegmax) then
            seg(itemp)%voiso  = K
            seg(itemp)%vnno   = K
          endif
          if (Nbcvxdom > IUN) then
            call assign_domain(k,11004)
          else
            seg(k)%dom  = IUN
          endif
          !print*, k, seg(K)%norme, seg(K)%O(:)
          !print*, segi, seg(segi)%norme, seg(segi)%O(:)
          call connecIJ(segi,K,243)
          if (kkdebug) write(379,*) "LIE ----------------------------"
          call voisinage(K,-244)
          if (sysconnec(VLi,seg(seg(segi)%voise)%veclin) /= 1 .and. &
              (seg(segi)%voise/=nsegmax) .and. seg(segi)%gd > izero) then
            seg(K)%gd = seg(segi)%gd
          else
            seg(K)%gd = izero
          endif
          seg(K)%DISEG    = .false.
          seg(K)%BLOQUER  = .false.
          seg(K)%unload   = .false.
          seg(k)%Nsysdev  = izero
          seg(k)%probadev = izero
          seg(k)%taudev   = izero
          if (kkdebug) call seginfo(SegI, " Apres traitement de LIE                   ")
      endif

      if ( LJO > lomaiR) then
          if (kkdebug) write(379,*)  "-------------------------    LJO > lomaiR      "
          Plusseg = plusseg + iun
          K = nsegm + Plusseg
          seg(k) = seg(segj)
          seg(k)%O(:) = seg(segj)%O(:)
          seg(K)%norme = LJO
          seg(k)%veclin = VLj
          seg(K)%jonc = .false. ; seg(K)%tjonc = izero ; seg(K)%Ijonc = nsegmax
          if (seg(segj)%surface==IUN) then
              seg(k)%surface=seg(segj)%surface
              seg(k)%VarFreePlan=seg(segj)%VarFreePlan
              seg(segj)%surface=IZERO
              seg(segj)%VarFreePlan=IZERO
          else
            seg(K)%surface=IZERO
            seg(K)%VarFreePlan=IZERO
          endif
          seg(segj)%O(:) = seg(segj)%O(:) + LJO * Bveclin(:,VLj)
          seg(segj)%norme = jlong - LJO
          seg(segj)%voiso  = K
          seg(K)%voise = segj
          seg(segj)%vnno  = K
          seg(K)%vnne = segj
          itemp = seg(K)%voiso
          if(itemp /= nsegmax) then
            seg(itemp)%voise = K
            seg(itemp)%vnne = K
          endif
          if (Nbcvxdom > IUN) then
            call assign_domain(k,11005)
            call assign_domain(segj,11006)
          else
            seg(k)%dom=IUN
            seg(segj)%dom=IUN
          endif
          call connecIJ(k,segj,245)
          if (kkdebug) write(379,*) "LIO --------------------------------"
          call voisinage(K,-246)
          if (sysconnec(VLi,seg(seg(segj)%voise)%veclin) /= 1 .and. &
              (seg(segj)%voise/=nsegmax) .and. seg(segj)%gd > izero) then
            seg(K)%GD = seg(segj)%gd
          else
            seg(K)%GD = izero
          endif
          seg(K)%DISEG = .false.
          seg(K)%BLOQUER = .false.
          seg(K)%unload = .false.
          seg(k)%Nsysdev = izero
          seg(k)%probadev = izero
          seg(k)%taudev = izero
          !CALL seginfo(Segj, " apres traitement de LIO                      ")
          if (kkdebug) call seginfo(Segj, " Apres traitement de LJO                   ")
      endif

      if ( LJE > lomaiR) then
          if (kkdebug) write(379,*)  "-------------------------     LJE > lomaiR      "
          Plusseg = plusseg + iun
          K = nsegm + Plusseg
          seg(K) = seg(segj)
          seg(K)%O(:) = seg(segj)%O(:) + (seg(segj)%norme - LJE) * Bveclin(:,VLj)
          seg(K)%norme = LJE
          seg(K)%veclin = VLj
          seg(segj)%voise  = K
          seg(K)%voiso = segj
          seg(segj)%vnne  = K
          seg(K)%vnno = segj
          seg(segj)%norme = seg(segj)%norme - LJE
          seg(K)%jonc = .false.
          seg(K)%tjonc = izero
          seg(K)%Ijonc = nsegmax
          if (seg(segj)%surface==IDEUX) then
              seg(k)%surface=seg(segj)%surface
              seg(k)%VarFreePlan=seg(segj)%VarFreePlan
              seg(segj)%surface=IZERO
              seg(segj)%VarFreePlan=IZERO
          else
            seg(K)%surface=IZERO
            seg(K)%VarFreePlan=IZERO
          endif
          itemp = seg(K)%voise

          if(itemp /= nsegmax) then
            seg(itemp)%voiso = K
            seg(itemp)%vnno = K
          endif
          if (Nbcvxdom > IUN) then
            call assign_domain(k,11007)
          else
            seg(k)%dom=IUN
          endif

          call connecIJ(segj,K,247)
          call voisinage(K,-248)

          if (sysconnec(VLi,seg(seg(segj)%voise)%veclin) /= 1 .and. &
              (seg(segj)%voise/=nsegmax) .and. seg(segj)%gd > izero) then
            seg(K)%GD = seg(segj)%gd
          else
            seg(K)%GD = izero
          endif
          seg(K)%DISEG = .false.
          seg(K)%BLOQUER = .false.
          seg(K)%unload = .false.
          seg(k)%Nsysdev = izero
          seg(k)%probadev = izero
          seg(k)%taudev = izero
          if (kkdebug) call seginfo(SegJ, " Apres traitement de LJE                 ")
      endif
    endif

  ! Cas ou on a pas de recouvrement ==> 2 cas:
  ! 1- Si tjonc > tjoncmax, on casse la jonction    2- Sinon, on incremente tjonc
  else

    if (seg(segi)%tjonc > tjoncmax) then

      if(kkdebug) write(379,*)  "Pas de recouvrement et tjonc >>, On CASSE LA JONCTION "
      k =seg(segI)%Ijonc
      seg(segI)%jonc = .false.
      seg(segi)%tjonc = izero
      seg(segi)%Ijonc = nsegmax
      seg(k)%jonc = .false.
      seg(k)%tjonc = izero
      seg(k)%Ijonc = nsegmax

      if(kkdebug) call Seginfo (segi, "I, apres destruction de la jonction ")
      if(kkdebug) call Seginfo (segj, "J, apres destruction de la jonction ")

    else

      if(kkdebug) write(379,*)  "Pas de recouvrement et tjonc << , ON INCREMENTE TJONC"
      seg(segi)%tjonc = seg(segi)%tjonc + IUN
      seg(segJ)%tjonc = seg(segJ)%tjonc + IUN

    endif

  endif

enddo INI

nsegm = nsegm + plusseg

if(nsegm > nsegmax) then
  write(*,*) 'Problem in (3) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif

PlusSeg = 0

! !***********************************************************
! !A loop to check if annihilation has been made correctly
! !***********************************************************
! B10: do segi = 1,NSEGM
!
!   if (out(segi)) cycle B10
!
!   if (seg(segi)%jonc) then
!
!     if (seg(segi)%tjonc <= tjoncmax+1) cycle B10
!     if (tyseg(seg(segi)%veclin) /= 1) cycle B10
!     if (syseg(seg(seg(segi)%ijonc)%veclin) /= sysdev(syseg(seg(segi)%veclin),1)) cycle B10
!     if (seg(segi)%norme /= 0 .and. seg(seg(segi)%ijonc)%norme /= 0) then
!
!       SegJ = seg(Segi)%Ijonc
!
!       VLi     = seg(Segi)%veclin
!       VLj     = seg(segj)%veclin
!       ilong   = seg(segi)%norme
!       jlong   = seg(segj)%norme
!
!       !Segments are shifted to avoid any problems coming from PBC
!       Ei(1:3)  = (ilong * Bveclin(1:3,VLi)) / ideux
!       TOi(1:3) = (seg(Segi)%O(1:3) + Ei(1:3) - modur(1:3)/2)
!       Oi(1:3)  = modulo(seg(segi)%o(1:3)-toi(1:3), modur(1:3))
!       Oj(1:3)  = modulo(seg(segJ)%o(1:3)-toi(1:3), modur(1:3))
!       Ei(1:3)  = modulo(seg(segi)%o(1:3) +                                     &
!                     ilong*bveclin(1:3,VLi)-toi(1:3), modur(1:3))
!       Ej(1:3)  = modulo(seg(segJ)%o(1:3) +                                     &
!                     jlong*bveclin(1:3,VLj)-toi(1:3), modur(1:3))
!
!       !The two segments are connected only with one point, so annihilation is not possible
!       if (sum(Oi(1:3)-Ej(1:3)) == 0 .or. sum(Oj(1:3)-Ei(1:3)) == 0) cycle B10
!
!       print*,'step = ',kk
!       print*, segi, seg(segi)%veclin,seg(segi)%norme, seg(segi)%gd, seg(segi)%jonc, seg(segi)%ijonc
!       print*, seg(segi)%ijonc, seg(seg(segi)%ijonc)%veclin,seg(seg(segi)%ijonc)%norme, &
!           seg(seg(segi)%ijonc)%gd, seg(seg(segi)%ijonc)%jonc, seg(seg(segi)%ijonc)%ijonc
!
!       call seginfo(segi, 'Critical junction problem i!!!')
!       call seginfo(seg(segi)%ijonc, 'Critical junction problem ijonc!!!')
!       stop 'Critical junction problem !!!'
!
!     endif
!
!   endif
!
! enddo B10            ! End if i segment do loop


end subroutine DISCRETI

!##############################################################
!> Procedure NET: The cleaning procedure
!>   A series of loops are used to test segments configuration and
!>   eliminate useless segments before the end of a time step.
!!##############################################################
subroutine NET

implicit none

integer(kind=DPI),parameter     :: NBOUCLE = 300
integer(kind=DPI)               :: IIVOIS(Nsegmax)
integer(kind=DPI)               :: I           !< Segment indice
integer(kind=DPI)               :: ib
integer(kind=DPI)               :: SIB
integer(kind=DPI)               :: Jindex
integer(kind=DPI)               :: JindexVeclin
integer(kind=DPI)               :: ibtest
integer(kind=DPI)               :: sibtest
integer(kind=DPI)               :: ILONBOU
integer(kind=DPI)               :: ILONBOUVIS
integer(kind=DPI)               :: J
integer(kind=DPI)               :: JJ
integer(kind=DPI)               :: ICONT2
integer(kind=DPI)               :: INDICE
integer(kind=DPI)               :: IAV
integer(kind=DPI)               :: IAP
integer(kind=DPI)               :: K
integer(kind=DPI)               :: PIECE
integer(kind=DPI)               :: Li          !< Segment line direction
integer(kind=DPI)               :: nsegm_lim   !< An uper limit for the nsegm tabs

integer(kind=DPI),DIMENSION(3)  :: Oi   !< Segment Origin
integer(kind=DPI),DIMENSION(3)  :: Ei   !< Segment End

real(kind=DP)                   :: SegLength         !< Cumulated not signed length of segments of loops
real(kind=DP)                   :: ProjScrewLength   !< Cumulated signed screw length of loops
real(kind=DP)                   :: ProjEdgeLength    !< Cumulated signed edge length of loops

real(kind=DP),DIMENSION(3)      :: Intersec
real(kind=DP),DIMENSION(3)      :: VecO   !< Vector segment origin / surface intersection point
real(kind=DP),DIMENSION(3)      :: VecE   !< Vector segment end / surface intersection point

logical                         :: ZAPONS(Nsegmax)
logical                         :: ibflag
logical                         :: nojunction  !< A dislocation line is free of any junction
logical,DIMENSION(nsegmax)      :: unload_old

integer(kind=DPI)               :: vnno, vnne
integer(kind=DPI)               :: seg_vnno, seg_vnne !< intermediate var for optimization
integer(kind=DPI)               :: sid_memo   !< The id of the segment we pin
real(kind=DP)                   :: JIndex_l   !< A length used in the unload procedure
real(kind=DP)                   :: sid_lmax   !< A length used in the unload procedure
real(kind=DP)                   :: sid_length !< A length used in the unload procedure
real(kind=DP)                   :: sid_ratio  !< A length ratio used in the infinit loop pinning procedure

logical,dimension(NbCvxDom)     :: InsideDomains  !< True if segment touching extremity is inside the domain i
real (kind=DP)                  :: DotProdO       !< Scalar Product

!*****************
! Initialization
!*****************
nsegm_lim             = nsegm * ideux
if (nsegm_lim > nsegmax) nsegm_lim = nsegmax

ZAPONS(1:nsegm_lim)   = .false.
PlusSeg               = izero

!*****************************************************************
! Very small loops accumulated as a result of cross-slip must be
! eliminated. They should dynamically disappear
! but the simulation time step is too large to reproduce that
!*****************************************************************
B1: do I = 1,NSEGM

    ! Segment I is part of a loop we already tested
    if (ZAPONS(I)) cycle B1

    ZAPONS(I) = .true. ! The segment I is noted as tested

    ! The segments we already know that they must be eliminated
    if (out(I) .or. (seg(I)%voiso == I .and. seg(I)%voise == I)) then
      seg(I)%probadev = zero
      seg(I)%wait = izero
      OUT(I) = .true.
      if(kkdebug) write(379,fmt='(" Loop B1 ", I9,": will be eliminated!" )') I
      cycle B1
    endif

    ! Segment I is a pinning point, we can skip this segment since it is not part of a loop
    if (SEG(I)%VOISO == NSEGMAX .or. SEG(I)%VOISE == NSEGMAX) cycle B1

    ! Initialization
    IIVOIS(1)   = I
    ILONBOU     = SEG(I)%NORME * nint(NormLin(SEG(I)%VECLIN),DPI)
    ILONBOUVIS  = ILONBOU
    if (TYSEG(SEG(I)%VECLIN) /= 1) ILONBOUVIS=IZERO

    ! Check if segment has a non zero vnno or vnne (paranoid test)
    seg_vnno = seg(i)%vnno
    seg_vnne = seg(i)%vnne
    if (((seg_vnno /= nsegmax) .and. ( seg(seg_vnno)%norme == izero )         &
        .and. seg(seg_vnno)%gd < iun .and. .not. seg(seg_vnno)%jonc) .or.     &
        ((seg_vnne /= nsegmax) .and. ( seg(seg_vnne)%norme == izero )         &
        .and. seg(seg_vnne)%gd < iun .and. .not. SEG(seg_vnne)%jonc)) then
      call seginfo(I,"VNNO or VNNE norm is zero (1)")
      call FLUSH(379)
      stop "VNNO or VNNE norm is zero (1)"
    endif

    !************************************
    ! Do loop on segments J connected to I
    ! J can be as big as nsegm + 1 when only one dislocation exist
    B2: do J = 2, nsegm + 1

      IIVOIS(J) = SEG(IIVOIS(J-1))%VOISO
      Jindex    = IIVOIS(J)

      if (Jindex == NSEGMAX) exit B2   ! A pining point is found we can exit

      ! The segments we already know that they must be eliminated
      if (out(jindex) .or. (seg(jindex)%voiso == jindex .and. seg(jindex)%voise == jindex)) then
        stop 'Something strange in loop B2 of NET'
      endif

      ! Check if segment has a non zero vnno or vnne
      seg_vnno = seg(Jindex)%vnno
      seg_vnne = seg(Jindex)%vnne
      if (((seg_vnno /= nsegmax) .and. ( seg(seg_vnno)%norme == izero )         &
          .and. seg(seg_vnno)%gd < iun .and. .not. seg(seg_vnno)%jonc) .or.     &
          ((seg_vnne /= nsegmax) .and. ( seg(seg_vnne)%norme == izero )         &
          .and. seg(seg_vnne)%gd < iun .and. .not. SEG(seg_vnne)%jonc)) then
        call seginfo(Jindex,"VNNO or VNNE norm is zero  (2)")
        call FLUSH(379)
        stop "VNNO or VNNE norm is zero  (2)"
      endif

      ! Calculation of the new line length increment
      PIECE   = SEG(Jindex)%NORME * nint(NormLin(SEG(Jindex)%VECLIN),DPI)
      ILONBOU = ILONBOU + PIECE                                              ! The total loop length
      if (TYSEG(SEG(Jindex)%VECLIN) == 1) ILONBOUVIS = ILONBOUVIS + PIECE    ! The screw section length

      if (ILONBOU < (IDEUX * MICROBOUCLE)) then     ! 1 tests on the loop size

        if (ILONBOU < MICROBOUCLE) then             ! 2 tests on the loop size (smaller)

          if(Jindex == I) then                      ! The loop is closed

!             ! --------------------------------------------------------------------------------
!             ! This section is a beta test not yet validated (DB)
!
!             ! We want to test if the micro_loop is an interstitial or vacancy type of loop
!             ! To do such test, we test the loop rotation sign with respect to the Burgers vector
!             ! Only Edge segments are considered for this test to speed up and simplify the calculation
!             vecnormX(:) = VecNorLin(:,seg(i)%veclin)   ! The test start with segment i (whatever the segment character)
!             Burgerveci(:) = VecNorLin(1:3,ASSOC(seg(i)%veclin,1))   ! The slip and cross-slip system Burgers vector
!             testloopsign = zero
!
!             Do ib = 2, J
!
!               sib = IIVOIS(ib)    ! the next segment to test
!
!               if (TYSEG(SEG(sib)%VECLIN) /= 3) cycle            ! We keep only the edge segments
!
!               vecnormY(:) = VecNorLin(:,seg(sib)%veclin)
!
!               if (abs(prodsca(vecnormX,vecnormY)) > 0.999999) cycle  ! The two vectors are idem
!               vech(:) = normavec(prodvec(vecnormX,vecnormY))    ! The normal to the two tested edge segments
!               signvech = prodsca(vech(:),Burgerveci(:))         ! Projection in the Burgers vector direction
!               testloopsign = testloopsign + signvech            ! We sum-up all contribution to determine the average loop direction
!               vecnormX(:) = vecnormY(:)
!
!             enddo
!
!             ! The size of vacancy loops to eliminate is decreased by a factor 4 for one sign of loops
!             if (ILONBOU > (MICROBOUCLE*QUART) .and. testloopsign > Un ) exit B2
!
!             ! --------------------------------------------------------------------------------

            if(kkdebug) call seginfo(IIVOIS(1),"Loop info 2nd step")

            ! All the segments part of this loop must be eliminated (J-1 segments)
            Do ib = 1, (J-1)

              sib = IIVOIS(ib)

              ! A segment in the loop take part to a junction
              if (seg(sib)%JONC) then

                k = SEG(sib)%IJONC

                !We liberate the complementary segment k
                if (seg(k)%norme /= 0) then
                  seg(k)%JONC = .false.       !simple case
                  seg(k)%IJONC = nsegmax
                  seg(k)%tJONC = 0
                else
                  ! Very particular cases are here considered.
                  ! A junction can be made between two segments part of the same loop
                  ! as an effect of collinear annihilation. In this case the STOPJONC
                  ! procedure cannot be used!
                  ibflag = .false.
                  do ibtest = ib + 1, J
                    sibtest = IIVOIS(ibtest)
                    ! We look for two segments forming one junction in a same loop
                    if (sibtest == k) then
                      ibflag = .true.       ! We flag this case
                      exit
                    endif
                  enddo

                  if (ibflag) then
                    ! The flagged cases must be eliminated without STOPJONC
                    seg(k)%JONC  = .false.
                    seg(k)%IJONC = nsegmax
                    seg(k)%tJONC = 0
                  else
                    call stopjonc(k,56306)
                    ! Attention, elimination of k can be requested by stopjonk
                    if(SEG(k)%voiso.eq.k .and. SEG(k)%voise.eq.k) OUT(k)=.true.
                  endif

                endif

                !We liberate the sib junction segment
                seg(sib)%JONC = .false.
                seg(sib)%IJONC = nsegmax
                seg(sib)%tJONC = 0

              endif

              ! Info regarding the GD segment eliminated in NET must be keeped
              if (seg(sib)%gd > iun) then
                oldntgd2 = oldntgd2 + 1
                if (kkdebug) write(379,*) 'In Net - 1, oldntgd2 after icrement =', oldntgd2
              elseif (seg(sib)%gd > izero) then
                oldntgd1 = oldntgd1 + 1
              endif

              OUT(sib) = .true.     ! Segment part of the loop are eliminated

              if(kkdebug) write(379,fmt='("2nd step ", I9,": loop annihilation ")') sib

            enddo

            ! All the segments in this loop have been noted to be eliminated, we can cycle to the next I
            exit B2

          endif

        else

          ! Those larger loops are eliminated if they are mostly of screw character
          if(Jindex == I) then        ! Closing test

            ! The loop is dominated by screw length
            if ((ILONBOU - ILONBOUVIS) < ILONBOUVIS) then

              if(kkdebug) call seginfo(IIVOIS(1)," info boucle etape 3")

              ! All the segments part of this loop must be eliminated (J-1 segments)
              Do ib = 1, (J-1)

                sib = IIVOIS(ib)

                ! A segment in the loop take part to a junction
                if (seg(sib)%JONC) then

                  k = SEG(sib)%IJONC

                  !We liberate the complementary segment junction
                  if (seg(k)%norme /= 0) then
                    seg(k)%JONC = .false.       !simple case
                    seg(k)%IJONC = nsegmax
                    seg(k)%tJONC = 0
                  else
                  ! Very particular cases are here considered.
                  ! A junction can be made between two segments part of the same loop
                  ! as an effect of collinear annihilation. In this case the STOPJONC
                  ! procedure cannot be used!
                    ibflag = .false.
                    do ibtest = ib+1, J
                      sibtest = IIVOIS(ibtest)
                      ! We look for two segments forming one junction in a same loop
                      if (sibtest == k) then
                        ! Oui le binome de jonc est aussi sur la boucle
                        ibflag = .true.
                        exit
                      endif
                    enddo

                    if (ibflag) then
                      ! The flagged cases must be eliminated without STOPJONC
                      seg(k)%JONC = .false.
                      seg(k)%IJONC = nsegmax
                      seg(k)%tJONC = 0
                    else
                      call stopjonc(k,66306)
                      ! Attention, elimination of k can be requested by stopjonk
                      if(SEG(k)%voiso.eq.k .and. SEG(k)%voise.eq.k) OUT(k)=.true.
                    endif
                  endif

                  !We liberate the sib junction segment
                  seg(sib)%JONC = .false.
                  seg(sib)%IJONC = nsegmax
                  seg(sib)%tJONC = 0

                endif

                ! Info regarding the GD segment eliminated in NET must be keeped
                if (seg(sib)%gd > iun) then
                  oldntgd2 = oldntgd2 + 1
                  if (kkdebug) write(379,*) 'In Net - 2, oldntgd2 after icrement =', oldntgd2
                elseif (seg(sib)%gd > izero) then
                  oldntgd1 = oldntgd1 + 1
                endif

                OUT(sib) = .true.     ! Segment part of the loop are eliminated

                if(kkdebug) write(379,*) " etape 3,annihilation boucle de I", sib

              enddo

              ! All segments in the loop are noted to be eliminated, we can cycle to the next I segment
              exit B2

            else

              ! This loop must not be eliminated, we can cycle to the next I segment
              exit B2

            endif   ! Test on long loops dominated by screw segments

          endif     ! Test on closed long loops

        endif       ! Test on short loops

      endif         ! Test to eliminate loop

      if (ZAPONS(Jindex)) exit B2      ! We endup on a section of dislocation already tested, no need to continue
      ZAPONS(Jindex) = .true.           ! We remember that this segment J is tested

    enddo B2        ! The J segment do loop

enddo B1            ! The I segment do loop

!***********************************************************
! A loop to check if all the gd segments are still usefull
!***********************************************************
!   B0: do i = 1,NSEGM
!
!     if (out(i)) cycle B0
!
!     if (seg(i)%gd > izero ) then
!       if (syseg(seg(seg(i)%vnno)%veclin) == syseg(seg(seg(i)%vnne)%veclin)) then
!
!         if (seg(i)%norme /= izero) then
!           Stop 'more work is needed here !!!!'
!         else
!           ! We keep the information about the gd eliminated
!           if (seg(i)%gd == iun) oldntgd1 =  oldntgd1 + iun
!           if (seg(i)%gd == ideux) oldntgd2 =  oldntgd2 + iun
!           seg(i)%gd = izero
!           ! The gd segment is replaced by a pivot segment compatible with the vnno and vnne
!           if (seg(seg(seg(i)%vnno)%voiso)%veclin /= nsegmax) then
!             seg(i)%veclin = seg(seg(seg(i)%vnno)%voiso)%veclin
!           elseif (seg(seg(seg(i)%vnne)%voise)%veclin /= nsegmax) then
!             seg(i)%veclin = seg(seg(seg(i)%vnne)%voise)%veclin
!           else
!             Stop 'No pivot segment could be defined in loop B0 of Net !'
!           endif
!
!           !print*, i, seg(i)%veclin,seg(i)%norme, seg(i)%gd
!           !print*, seg(i)%vnno,seg(seg(i)%vnno)%veclin,syseg(seg(seg(i)%vnno)%veclin)
!           !print*, seg(i)%vnne,seg(seg(i)%vnne)%veclin,syseg(seg(seg(i)%vnne)%veclin)
!           !print*,SEGdev(seg(i)%veclin,syseg(seg(seg(i)%vnno)%veclin)),SEGdevIJ(seg(i)%veclin,seg(seg(i)%vnno)%veclin)
!           !read *
!         endif
!
!         stop ' A problem is found in loop B0 in NET'
!
!       endif
!     endif
!
!   enddo B0            ! End if i segment do loop


!*******************************************************************************************
! In this loop we eliminate the small segments that may touch two free surfaces at both ends
! This loop is justified only if GB /= 0
if (GB /= izero) then
  B3: do I = 1,NSEGM
    if(out(i) .or. seg(i)%voiso == seg(i)%voise) cycle
    if (Nbcvxdom > IUN) then
      if((seg(i)%vnno == NSEGMAX .and. seg(i)%surface /= IUN   .and. seg(i)%voiso /= NSEGMAX) .or. &
       (seg(i)%vnne == NSEGMAX .and. seg(i)%surface /= IDEUX .and. seg(i)%voise /= NSEGMAX)) then

        call seginfo(i, 'surface problem 2')
        stop 'surface problem 2'

      endif
    endif

    if (seg(I)%surface > IUN .and. seg(I)%surface < ITROIS) then

      IIVOIS(1) = I

      ! Initialisation du compteur de longueur des micro-boucles
      ILONBOU = SEG(I)%NORME*nint(NormLin(SEG(I)%VECLIN),DPI)

      B4: do J = 2,NBOUCLE

        IIVOIS(J) = SEG(IIVOIS(J-1))%VOISO
        Jindex    = IIVOIS(J)
        PIECE=SEG(Jindex)%NORME*nint(NormLin(SEG(Jindex)%VECLIN),DPI)
        ILONBOU = ILONBOU+PIECE

        ! The next segment is in the volume and is not a free ending segment
        if(seg(Jindex)%surface < IUN .and. seg(Jindex)%voiso /= nsegmax ) cycle B4

        if(seg(Jindex)%surface == IUN .or. seg(Jindex)%surface == ITROIS .or. seg(Jindex)%voiso == nsegmax) then

          if (ILONBOU < MICROBOUCLE) then

            ! This dislocation is touching free surface at both ends
            ! and it is small. It must be eliminated
            Do ib = 1, J

              sib = IIVOIS(ib)
              if (seg(sib)%JONC) then
                k = SEG(sib)%IJONC
                !We liberate the complementary segment junction
                if (seg(k)%norme /= 0) then
                  seg(k)%JONC = .false.       !simple case
                  seg(k)%IJONC = nsegmax
                  seg(k)%tJONC = 0
                else
                  ! si les deux segments jonc sont sur la microboucle
                  ! (c est possible avec les annihilations colli) il
                  ! ne faut pas faire a stopjonc qui risque d
                  ! introduire de nouveau segments
                  ibflag = .false.
                  do ibtest = ib+1, J
                      sibtest = IIVOIS(ibtest)
                      if (sibtest == k) then
                        ! Oui le binome de jonc est aussi sur la boucle
                        ibflag = .true.
                        exit
                      endif
                  enddo
                  if (ibflag) then
                    ! cas particulier des jonc sur une meme boucle
                    seg(k)%JONC = .false.
                    seg(k)%IJONC = nsegmax
                    seg(k)%tJONC = 0
                  else
                    call stopjonc(k,56306)            !non trivial case
                    ! Attention, elimination of k can be requested by stopjonk
                    if(SEG(k)%voiso.eq.k .and. SEG(k)%voise.eq.k) OUT(k)=.true.
                  endif
                endif

                !We liberate the junction
                seg(sib)%JONC = .false.
                seg(sib)%IJONC = nsegmax
                seg(sib)%tJONC = 0
              endif

              OUT(sib) = .true.   !! To be eliminated !!
              seg(sib)%surface = IZERO
              seg(sib)%VarFreePlan = IZERO
              if(kkdebug) then
                write(379,*) " etape 3,annihilation of dislocation touching free surfaces on both sides", sib
                call seginfo(sib," info boucle etape 3")
              endif

            enddo

            OUT(seg(IIVOIS(1))%voise) = .true.
            seg(seg(IIVOIS(1))%voise)%surface = IZERO
            seg(seg(IIVOIS(1))%voise)%VarFreePlan = IZERO
            OUT(seg(Jindex)%voiso) = .true.
            seg(seg(Jindex)%voiso)%surface = IZERO
            seg(seg(Jindex)%voiso)%VarFreePlan = IZERO

            exit B4   !exit

          endif

          exit B4
        endif

      enddo B4

    endif

  enddo B3
endif

! initialisation
ICONT2 = -1
! A problem exist here. Out() cannot be initialized here since segments in small
! loops to eliminate was defined just above and can be in between nsegm and
! nsegm+plusseg
!OUT((NSEGM+1):(NSEGM+PlusSeg)) = .false.

nsegm = nsegm + PlusSeg
if(nsegm > nsegmax) then
  write(*,*) 'Problem in (4) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif
PlusSeg = 0

!************************************************************************
! Test on the length and coordinates of segments touching a free surface
!************************************************************************
if(kkdebug) write(379,*)  "  > Test on the length and coordinates of segments touching a free surface"

if (GB /= 0) then

  B5: do I = 1,NSEGM ! Loop on all segments

    !Select segments which are touching a free surface and which are not meant to be deleted yet
    if (seg(i)%surface > Izero .and. .not. OUT(i)) then

      if(kkdebug) then
        write(379,fmt='(">>> In NET:B5 - looking at segment ",I7,": touching the FS ",I4," dom ",I4)') &
        i,seg(i)%VarFreePlan,seg(i)%dom
      endif

      ! We look for segment not correctly defined with respect to free surface
      if (seg(i)%surface > ITROIS) then !< Key to identify segments touching free surfaces (0=no,1=origin,2=end)
        print *, "Found a segment surface key > 3 in subroutine NET",i,seg(i)%surface
        call seginfo(i,"surf prob")
        stop
      endif

      ! If we still have a kneecap touching a free surface, the work has not correctly be done
      if (seg(i)%surface > IZERO .and. seg(i)%norme==IZERO .and. .not. OUT(i)) then
        print *,"surface pinning point problem", i ,"kk",kk
        call seginfo(i,"pin")
        stop
      endif

      ! Usefull informations on segment I
      Oi(1:3) = seg(i)%o(1:3) ! Origin coordinates
      Li = seg(i)%veclin      ! Index of line vector and displacement
      J  = seg(i)%VarFreePlan ! Index of free surface touched by the segment

      ! Point of intersection between the direction of the segment and the free surface J
      Intersec(1:3) = real(InterPlanSeg(Plane_MillerR(1:3,j),Plane_pos(j),real(Oi,DP),BVECLIN(:,Li)))

      Ei(1:3) = Oi(1:3) + seg(i)%norme * BVECLIN(1:3,Li) !End of I from its origin, direction and norm

      InsideDomains(:) = .true.

      if (seg(i)%surface > IUN) then !If segment touch the free surface with its end (case 2)

        ! We check that the End of segments I is touching a free surface is correctly set
        Ei(1:3) = Oi(1:3) + seg(i)%norme * BVECLIN(1:3,Li)

        VecE(1:3) = real(Ei(1:3) - Intersec(1:3),DP) ! Vector from the intersection point on FS and segment I end

        ! We check that the segment E is outside - Meaning dot product (intersect,E) with plane normale positive
        if (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) < ZERO .and. kk > IZERO .and. &
           norvect(intersec) > 1e-12) then !and segment end and intersection point are not the same (recall that norvect return a float)
          print*, "Stop! in NET, a segment with surface index IUN has its end inside the volume  ","I : ",i,"KK : ",KK
          print*, "Line Vector : ", BVECLIN(1:3,Li)
          print*, "Segment origin :",Oi(1:3)
          print*, "Segment end : ", Ei(1:3)
          print*, "Segment norme :", seg(i)%norme
          print*, "Free surface :", seg(i)%VarFreePlan
          print*, "Intersection point",intersec(1:3)
          print*, "Norme of Intersection point",norvect(intersec)
          print*, "Dot product between vector (Intersect - Ei) and Plane normal :",dot_product(VecE(1:3),Plane_MillerR(1:3,j))
          call seginfo(i,"Segment with surface IDEUX has its end inside the volume")
         stop
        endif

        !Testing that the extremity point minus one veclin is in a domain

        !VecE is now a point
        VecE(1:3) = real(Ei(1:3),DP) - real(BVECLIN(1:3,Li),DP)


        Do jj= 1,NbPlanMax ! Loop over surfaces
          If ( .Not. InsideDomains(Plane_dom(jj)) ) cycle !If we already know that segment extremity is outside the domain of the plane, cycle

          normsurf(:) = Plane_MillerR(1:3,jj)
          DotProdO=VecE(1)*normsurf(1)+VecE(2)*normsurf(2)+VecE(3)*normsurf(3)-Plane_pos(jj)

          if (kkdebug) write (379,fmt='("VFP",x,I2,x,"dom",x,I2,x," DotProdO",F15.3)') &
                                        & jj,Plane_dom(jj), DotProdO

          if (DotProdO > numtol_dotp) then
            InsideDomains(Plane_dom(jj)) =.False.
          endif

        Enddo

        if (.not. Any(InsideDomains(:)) ) then
            write (379,fmt='("!> Connection to free surface error")')
            write (379,fmt='("!> The segment is attached to far from the free surface - case IDEUX")')
            call seginfo(i,"The segment is attached to far from the free surface - case IDEUX")
            stop "Connection to free surface error case IDEUX - check debug file"
        endif



      else ! If segment I is touching the free surface with its origin

        ! We check that the connection to the free surface is correctly set

        VecO(1:3) = real(Oi(1:3) - Intersec(1:3),DP) ! Vector from the intersection point on FS and segment I origin

        ! We check that the segment O is outside - Meaning dot product (intersect,0) with plane normale positive
        if (dot_product(VecO(1:3),Plane_MillerR(1:3,j)) < ZERO .and. kk > IZERO .and. &
           norvect(intersec) > 1e-12) then !and segment end and intersection point are not the same (recall that norvect return a float)
          print*, "Stop! in NET, a segment with surface index 2 or 3 has its end inside the volume  ","I : ",i,"KK : ",KK
          print*, "Line Vector : ", BVECLIN(1:3,Li)
          print*, "Segment origin :",Oi(1:3)
          print*, "Segment end : ", Ei(1:3)
          print*, "Segment norme :", seg(i)%norme
          print*, "Free surface :", seg(i)%VarFreePlan
          print*, "Intersection point",intersec(1:3)
          print*, "Norme of Intersection point",norvect(intersec)
          print*, "Dot product between vector (Intersect - Ei) and Plane normal :",dot_product(VecO(1:3),Plane_MillerR(1:3,j))
          call seginfo(i,"Segment with surface index 2 has its origin inside the volume")
         stop
        endif

        !Testing that the extremity point minus one veclin is in a domain

        !VecO is now a point
        VecO(1:3) = real(Oi(1:3),DP) + real(BVECLIN(1:3,Li),DP)


        Do jj= 1,NbPlanMax ! Loop over surfaces
          If ( .Not. InsideDomains(Plane_dom(jj)) ) cycle !If we already know that segment extremity is outside the domain of the plane, cycle

          normsurf(:) = Plane_MillerR(1:3,jj)
          DotProdO=VecO(1)*normsurf(1)+VecO(2)*normsurf(2)+VecO(3)*normsurf(3)-Plane_pos(jj)

          if (kkdebug) write (379,fmt='("VFP",x,I2,x,"dom",x,I2,x," DotProdO",F15.3)') &
                                        & jj,Plane_dom(jj), DotProdO

          if (DotProdO > numtol_dotp) then
            InsideDomains(Plane_dom(jj)) =.False.
          endif

        Enddo

        if (.not. Any(InsideDomains(:)) ) then
            write (379,fmt='("!> Connection to free surface error")')
            write (379,fmt='("!> The segment is attached to far from the free surface - case IUN")')
            call seginfo(i,"The segment is attached to far from the free surface - case IUN")
            stop "Connection to free surface error case IUN - check debug file"
        endif

      endif

      ! Test if there is no unconsistency between free surface index and connection type
      if ((seg(i)%surface /= IZERO .or. seg(i)%varfreeplan /= IZERO) &
         .and. (seg(i)%surface*seg(i)%varfreeplan == IZERO)) then
        print*, "found seg(i)%surface and seg(i)%varfreplan not consistent in subroutine NET",i, kk
        call seginfo(i,"seg surf consistency")
        stop
      endif

    endif ! End of the tests on segments touching a free surface and not set to deletion
    !We verify that domain attribution has been done correctly
    if (kk>IZERO .and. .not. out(i)) call check_domain(i,321)
  enddo B5
endif !GB /= ZERO

! Initialization
ZAPONS(1:nsegm_lim)   = .false.
unload_old(1:nsegm)   = seg(1:Nsegm)%unload ! The previous list of unload segments is saved
seg(1:nsegm_lim)%unload = .false.
frankloop             = zero
frankline             = zero
frankfreeline         = zero
LoopNumb              = izero
FreeLoopNumb          = izero

!************************************************************************
! Statistique on dislocation loops and pinning of free infinite loops
!************************************************************************
B7: do I = 1,NSEGM

    ! Segment I is part of a loop we already tested
    if (ZAPONS(I)) cycle B7

    ZAPONS(I) = .true. ! The segment I part of a new dislocation is now tested

    ! The segments we already know that they must be eliminated
    if (out(I) .or. (seg(I)%voiso == I .and. seg(I)%voise == I)) cycle B7

    ! Segment I is a pinning point, we can skip this segment since it is not part of a loop
    if (SEG(I)%VOISO == NSEGMAX .or. SEG(I)%VOISE == NSEGMAX) cycle B7

    ! Initialization
    nojunction      = .true.   ! Dislocations loops are initially supposed not pined by junctions
    IIVOIS(1)       = I
    JIndex_l        = real(seg(I)%norme,DP)
    SegLength       = zero
    ProjScrewLength = zero
    ProjEdgeLength  = zero

    ! at least one effective junction is pining this line section
    if (seg(I)%jonc .and. (seg(I)%tjonc > (tjoncmax+1))) nojunction = .false.

    !************************************
    ! Do loop on segments J connected to I
    ! J can be as big as nsegm when only one dislocation exist
    B8: do J = 2, nsegm + 1

      IIVOIS(J) = SEG(IIVOIS(J-1))%VOISO
      Jindex    = IIVOIS(J)

      if (Jindex == NSEGMAX) exit B8     ! A pining point is found we must exit B8

      ! Segment length resolved in the Burgers direction is summed to test the Frank rule for the loops
      JIndex_l        = real(seg(Jindex)%norme,DP)
      JindexVeclin    = seg(Jindex)%veclin
      if (seg(Jindex)%jonc) then
        SegLength     = SegLength + JIndex_l * normlin(JindexVeclin) * half  ! junction density must be devided by two
      else
        SegLength     = SegLength + JIndex_l * normlin(JindexVeclin)
      endif
      ProjScrewLength = ProjScrewLength + JIndex_l * ProjVis(JindexVeclin)
      ProjEdgeLength  = ProjEdgeLength  + JIndex_l * Projcoin(JindexVeclin)

      ! At least one junction is pining this line section
      if (seg(Jindex)%jonc .and. (seg(Jindex)%tjonc > (tjoncmax+1))) nojunction = .false.

      ! The case of loops that may be infinit loop artefact made by the PBC
      if(Jindex == I) then

        FrankLoop = FrankLoop + SegLength  ! Loop length are summed
        LoopNumb  = LoopNumb + 1           ! The number of loop found is incremented

        if (abs(ProjScrewLength) > 0.1 .or. abs(ProjEdgeLength) > 0.1) then

          FrankLine     = FrankLine + SegLength     !Cumulated length of Inf lines
          FrankLoop     = FrankLoop - SegLength     ! This is not a real loop, but an inf line
          FreeLoopNumb  = FreeLoopNumb + 1          ! Number of free inf lines
          LoopNumb      = LoopNumb - 1              ! This is not a real loop, but an inf line

          ! The special case of inf line not pined by any junction
          if (nojunction) then

            FrankFreeLine = FrankFreeLine + SegLength   ! length of free Inf lines

            if (key_infline) then

              ! This dislocation line must be unloaded (pinned) since nothing limits its displacements
              sid_memo   = -999
              sid_lmax   = -999.
              sid_ratio  = 0.5

              ! Ieration with smaller and smaller length of segment to pin
              do
                ! Loop on the segment part of the infint loop
                do ib = 1, J

                  sib = iivois(ib)

                  ! If a segment was 'unloaded' at the previous step we always use it again
                  if (unload_old(sib)) then
                    sid_memo   = sib
                    exit
                  endif

                  ! We want to pin a not too small segment and as far as possible from the point (0,0,0)
                  if (real(seg(sib)%norme)/real(loma(seg(sib)%veclin)) > sid_ratio) then

                    sid_length = sqrt( real(modulo(seg(sib)%o(1),modur(1))**2    + &
                                            modulo(seg(sib)%o(2),modur(2))**2    + &
                                            modulo(seg(sib)%o(3),modur(3))**2) )

                    if (sid_length > sid_lmax) then

                      sid_lmax = sid_length
                      sid_memo = sib
                      sid_ratio = real(seg(sib)%norme)/real(loma(seg(sib)%veclin))
                      !print*,'pinned ib=', ib, J, sid_lmax, sid_memo, sid_ratio
                      !read *

                    endif

                  endif

                enddo

                ! The segment at the border of the box is pinned
                if (sid_memo < izero) then

                  ! The test is made again, but with a smaller value for sid_ratio
                  sid_ratio =  sid_ratio - 0.1

                  if (kkdebug)  write(379,*) " Iteration on sid_ratio are needed to find a segment to pin!", sid_ratio

                  ! too many iterations. We must stop
                  if (sid_ratio < 0.01) then
                    write(379,*) "A problem is found in the infinite segment procedure at step = ", KK
                    write(379,*) "sid_ratio = ",sid_ratio
                    write(379,*) "Number of segments part of the loop = ",J
                    Do ib = 1, J
                      sib = iivois(ib)
                      vnno = seg(sib)%vnno
                      vnne = seg(sib)%vnne
                      write(379,*) sib, vnno, vnne, seg(sib)%O(:), seg(sib)%norme
                    enddo
                    stop 'Pb in the infinit segments procedure, see the debug files'
                  endif

                else

                  ! We found the segment to pin
                  exit

                endif

              enddo

              ! This segment must be unload to pin locally this infinit free line
              seg(sid_memo)%unload = .true.
              seg(sid_memo)%wait = izero     ! Wait must be initialized to visite the force calculation, hence pinning
              seg(sid_memo)%bloquer = .true.
              seg(sid_memo)%resdep  = zero
              if (kkdebug)  write(379,*) " The segment = ", sid_memo, " is noted to unload!"

            endif

          endif

        endif

        exit B8

      endif

      if (ZAPONS(Jindex)) exit B8          ! We endup on a section already tested, no need to continue
      ZAPONS(Jindex) = .true.               ! We remember that this segment J has been tested

    enddo B8        ! The J segment do loop

enddo B7          ! The I segment do loop

!******************************************************
! Elimination of unneeded segments (OUT(I) = true)
!******************************************************
do I = 1,nsegm

  if ((I+(ICONT2+1)) > NSEGM) exit ! test de sortie du tri

  if (OUT(I)) then                  ! Ce segment doit etre elimine

    seg(i)%surface      = IZERO  ! to avoid any problem variables of segment touching free surfaces
    seg(i)%varfreeplan  = IZERO  ! are reinitialized

69  CONTINUE

    ICONT2 = ICONT2 + 1

    if((I+(ICONT2+1)) > NSEGM) exit ! test de sortie du tri

    INDICE = NSEGM - ICONT2

    if(OUT(INDICE)) then
      if (seg(indice)%jonc) then
        print *,"NET1: I = ",i, " a supprimer et jonc ??"
        k = SEG(indice)%IJONC
        seg(k)%jonc = .false.
        seg(k)%Ijonc = nsegmax
        seg(Indice)%jonc = .false.
        seg(Indice)%Ijonc = nsegmax
      endif
      goto 69
    endif
    ! Cleanup operation

    ! Si un segment qui participait a une jonction est elimine, sont binome doit etre relaxe.
    if (seg(i)%jonc) then
      print *,"NET2: I = ",i, " a supprimer et jonc ??"
      k = SEG(I)%IJONC
      seg(k)%jonc = .false.
      seg(k)%Ijonc = nsegmax
      seg(I)%jonc = .false.
      seg(I)%Ijonc = nsegmax
    endif

    if (kkdebug) write(379,*) "NET: le seg = ",i, " sera ecrase par : ",INDICE,out(indice)
    out(I)    = .false.
    seg(I)    = seg(INDICE)
    LoiSEG(I) = LoiSEG(INDICE)

    ! Si le nv segment participait a une jonction il faut remmetre a jour son binome
    if (seg(INDICE)%jonc) then
      k = SEG(Indice)%IJONC
      if(.not. seg(k)%jonc) then
        print *,"NET: I = ",i, " jonction assymetrique ??, KK = ",KK
        seg(I)%jonc = .false.
        seg(I)%Ijonc = nsegmax
      else
        SEG(K)%IJONC = I
      endif
    endif

    IDEP(I) = IDEP(INDICE)
    DAIRE(I) = DAIRE(INDICE)
    ! Gestion des vosins de I
    IAV = SEG(INDICE)%VOISO
    IAP = SEG(INDICE)%VOISE
    ! On modifie le type derive pour les segments en contact
    if(IAV.ne.NSEGMAX) SEG(IAV)%VOISE = I
    if(IAP.ne.NSEGMAX) SEG(IAP)%VOISO = I
    call voisinage(i,100)
    !if(Indice == 20806) call seginfo(i," I apres  voisinage")

    seg(INDICE)%SURFACE=IZERO
    seg(INDICE)%VARFREEPLAN=IZERO

  endif

enddo

NSEGM = NSEGM-(ICONT2+1)

if (nsegm == IZERO) stop ' We have lost all the segments'

! debugage  : verification de la symetrie des sourcess de frank-read
if (nsegm < 1) then
  I = 1
  do while (seg(I)%voiso /= nsegmax)
      I = I + 1
  enddo

  J = I
  I = 1
  do while (seg(I)%voise /= nsegmax)
      I = I + 1
  enddo

  K = I
  do while (J /= K .and. seg(J)%voise /= K)
      if(kkdebug) write(379,*)  "J,k",J,K
      if(seg(J)%norme /= seg(K)%norme .or. abs(seg(j)%resdep-seg(k)%resdep) > 0.0001) then
        print *, " probelme entre ",J,K
        print *, seg(j)%resdep,seg(k)%resdep
        call disloinfo("source non symetrique dans net ")
        stop
      endif
      j = seg(j)%voise ; k = seg(k)%voiso
  enddo
endif

!reinitialitation of the segments just eliminated
seg(nsegm + 1: nsegm + ICONT2+1)%resdep   = zero
seg(nsegm + 1: nsegm + ICONT2+1)%depinst  = zero
seg(nsegm + 1: nsegm + ICONT2+1)%probadev = zero
seg(nsegm + 1: nsegm + ICONT2+1)%taudev   = zero
seg(nsegm + 1: nsegm + ICONT2+1)%anglevis = -un

seg(nsegm + 1: nsegm + ICONT2+1)%O(1)     = izero
seg(nsegm + 1: nsegm + ICONT2+1)%O(2)     = izero
seg(nsegm + 1: nsegm + ICONT2+1)%O(3)     = izero
seg(nsegm + 1: nsegm + ICONT2+1)%veclin   = izero
seg(nsegm + 1: nsegm + ICONT2+1)%norme    = izero
seg(nsegm + 1: nsegm + ICONT2+1)%voiso    = nsegmax
seg(nsegm + 1: nsegm + ICONT2+1)%voise    = nsegmax
seg(nsegm + 1: nsegm + ICONT2+1)%vnno     = nsegmax
seg(nsegm + 1: nsegm + ICONT2+1)%vnne     = nsegmax
seg(nsegm + 1: nsegm + ICONT2+1)%ijonc    = nsegmax
seg(nsegm + 1: nsegm + ICONT2+1)%tjonc    = izero
seg(nsegm + 1: nsegm + ICONT2+1)%wait     = izero
seg(nsegm + 1: nsegm + ICONT2+1)%grain    = izero
seg(nsegm + 1: nsegm + ICONT2+1)%dom      = izero
seg(nsegm + 1: nsegm + ICONT2+1)%surface  = izero
seg(nsegm + 1: nsegm + ICONT2+1)%VarFreePlan = izero
seg(nsegm + 1: nsegm + ICONT2+1)%nsysdev  = izero
seg(nsegm + 1: nsegm + ICONT2+1)%NPhase   = izero

seg(nsegm + 1: nsegm + ICONT2+1)%jonc     = .false.
seg(nsegm + 1: nsegm + ICONT2+1)%gd       = izero
seg(nsegm + 1: nsegm + ICONT2+1)%diseg    = .false.
seg(nsegm + 1: nsegm + ICONT2+1)%bloquer  = .false.
seg(nsegm + 1: nsegm + ICONT2+1)%unload   = .false.
seg(nsegm + 1: nsegm + ICONT2+1)%zerotl   = .false.

#ifdef MDC
seg(nsegm + 1: nsegm + ICONT2+1)%SIGFE(1) = zero
seg(nsegm + 1: nsegm + ICONT2+1)%SIGFE(2) = zero
seg(nsegm + 1: nsegm + ICONT2+1)%SIGFE(3) = zero
seg(nsegm + 1: nsegm + ICONT2+1)%SIGFE(4) = zero
seg(nsegm + 1: nsegm + ICONT2+1)%SIGFE(5) = zero
seg(nsegm + 1: nsegm + ICONT2+1)%SIGFE(6) = zero
#endif

! The out variable must be reinitialize from 1 since it is used for the free surface part
! (defining the segment touching for the first time a free surface) in module microstructure.
OUT(1 : nsegm + ICONT2+1)=.false.

end subroutine NET

end module TOPOLO

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
