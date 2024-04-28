
!===================================================================================================
!========================    DEBUT    MODULE   "RESUL"  ============================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the dislocation dynamics mathematical analysis to
!> extract pertinent mechanical and physical properties.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module RESUL

use VARGLOB

use BRICAMAT
use DEBUG
use Microstructure

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Subroutine where the main outputs of mm simulations are calculated !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine STATIS

implicit none

real(kind=DP) :: UNDELTA      !<
real(kind=DP) :: XLONG        !<
real(kind=DP) :: InvVol       !< 1 / Volume (a unit)
real(kind=DP) :: const2       !<
real(kind=DP) :: const3       !<
real(kind=DP) :: temp         !<
real(kind=DP) :: factspheric  !<

real(kind=DP) :: SONseuil       !<
real(kind=DP) :: raudisact      !<
real(kind=DP) :: RAU_GROSEG     !< Total density of long segments at step kk
real(kind=DP) :: Px,Py,Pz       ! Used to write the gammabox file

real(kind=DP), DIMENSION(NTSG) :: rauforet !<
real(kind=DP), DIMENSION(NTSG) :: somme    !< To check the calculation of gammabox

integer(kind=DPI)   :: I          !<
integer(kind=DPI)   :: Li         !<
integer(kind=DPI)   :: ISYST      !<
integer(kind=DPI)   :: nloi       !<
integer(kind=DPI)   :: JK         !<
integer(kind=DPI)   :: sigsursys  !<
integer(kind=DPI)   :: ns         !<
integer(kind=DPI)   :: ii         !<
integer(kind=DPI)   :: j          !<
integer(kind=DPI)   :: jj         !<
integer(kind=DPI)   :: ll         !<
integer(kind=DPI)   :: k          !<
integer(kind=DPI)   :: sysact     !<

integer(kind=DPI)   :: ni,isum,jsum   ! Used to write the gammabox file

integer,parameter   :: nsmax = 12       ! Maximum number of slip systems considered in the output files
real(kind=DP)       :: toto(nsmax)      ! Intermediate tabs needed for simple writing
integer(kind=DPI)   :: itoto(nsmax)

#ifdef PA
real(kind=DP)       :: Eps_pl_mphase, &! Total plastic deformation in the loading direction in case of a multi-phase simulation
                       temp_phase
#endif

!*************************************************************************
!*** Reinitialisation des grandeurs instantannees des calcul de STATIS ***
!*************************************************************************

! RAUDMO,VITMOY,EPSDOT,EPS_PL,airevis,airecoin,airemixte

RAUDIS = 0.d0        !*** Densite totale des disloctions
RAUDIS_matrice = 0.d0      !*** Densite par systeme dans la matrice
RAUDIS_precipite = 0.d0    !*** Densite par systeme dans le precipite
RAUDMO_I = 0.d0      !*** Densite des dislocation mobile moyennee sur 4 iterations (m-2)
RAUSYS(:) = ZERO     !*** Densite par systeme
RAU_GROSEG = ZERO   ! Density of segments not treated with the multipole approach
RAUSYS_precipite(:) = ZERO     !*** Densite par systeme
NSEGMJONC = 0        !*** Compteur dislo mobiles
VITA = 0.0D0         !*** Accumulateur pour la vitesse moyenne des dislo
GAMMADOTSYS(:)=0.D0
ns = NTSG
InvVol = 1.0/(VOLUME*AVALUE*AVALUE)
const3 = zero
factspheric = zero

UNDELTA = UN/DELTAT*AVALUE
SONseuil = SONDMS*0.95
TrAppInst = zero
TrIntInst = zero

! Array initialization for DesorientGrain
if (desorientgrain) then
  RAUGRAIN(:)=0.d0              ! Density in the grain
  RAUGRAINSYS(:,:)=0.d0         ! Density per slip systems in the grains
  TensRotGrain(:,:,:)=0.d0
  gamsysGRAIN(:,:)=0.d0
  const3 = InvVol*nbgrains      ! Volume of the grains (supposed to be regular)
endif

if (GB==IDEUX) factspheric =  VOLUME*AVALUE**3/(QUATRE/TROIS*PII*(GB2_size*1.D-6)**3)

!*******************************
!*** Boucle sur les segments ***
!*******************************
do I = 1,NSEGM
    Li = SEG(I)%VECLIN
    XLONG = SEG(I)%NORME*NormLin(Li)
    ISYST = SYSEG(Li)
    !*** Cumule des longeurs pour la densite totale de dislocation
    RAUDIS = RAUDIS+XLONG

    !*** Calculation of the total density of long segments   ***
    if(MOD(KK,KSTAT).eq.IZERO) then
      if(seg(I)%norme * normlin(seg(I)%veclin)  > lengthlimit) then
          RAU_GROSEG = RAU_GROSEG + Xlong
      endif
    endif


    if (seg(i)%NPhase == IUN) RAUDIS_precipite = RAUDIS_precipite+XLONG

    !***Cumul pour chaque grain
    if (DesorientGrain) then
      call CorrectifGrain(I,3005) !verif necessaire les initialisations ne sont pas tres stables
      RAUGRAIN(seg(i)%grain)=RAUGRAIN(seg(i)%grain)+xlong
    endif
    !*** 08/11/00 : STAT SUR LES JONCTIONS
    if(seg(i)%jonc) then
      NSEGMJONC = NSEGMJONC+1
    endif
    !*** Statistique sur la densite par systeme de glissement
    RAUSYS(ISYST) = RAUSYS(ISYST)+XLONG

    if (seg(I)%NPhase ==iun) RAUSYS_precipite(ISYST) = RAUSYS_precipite(ISYST)+XLONG

    !*** Statistique sur la densite par systeme de glissement pour chaque
    !*** grain
    if (DesorientGrain) RAUGRAINSYS(seg(i)%grain,ISYST) = RAUGRAINSYS(seg(i)%grain,ISYST)+XLONG
enddo

!************************************************************************************
!************************************************************************************
! Determination du deplacement maximal d'ecranter les deplacements excessives :

Do nloi = 1,NLV

  !************************************************************************************
  !* calcul des deplacement moyens : vitesse initiale (brute), vitesse ecrantee,VITESSE vraie
  !************************************************************************************

  ! calculation of a mean displacement per unit length (avalue)
  if (Lseg_V(nloi) > zero) D_vraie(nloi) = D_vraie(nloi) / LSEG_V(nloi)

  if (KK < Relax_TL) then
    Dlimite(nloi) = Dlimite_max
  else
    Dlimite(nloi) =  D_vraie(nloi) * facteur_depmax
    if (Dlimite(nloi) < Modep_max)   Dlimite(nloi) = Dlimite_max
    if (Dlimite(nloi) > Dlimite_max) Dlimite(nloi) = Dlimite_max
  endif

  if (kk == iterinfo) &
      write(379,'("nloi =", I4,"  Lseg=",F7.3,"  Di =",F7.3,"  Dv =",F7.3," Dlim=",F7.2)') &
      nloi, Lseg(nloi)*avalue*1.D6,D_brute(nloi),D_vraie(nloi),dlimite(nloi)

enddo
!************************************************************************************


!*** Vitesse moyenne en valeur absolue des dislocations (m/s)
if(RAUDMO_I.ne.ZERO) then
  VITMOY = (VITMOY*TROIS+VITA/RAUDMO_I)*QUART
endif

!*** Densite par loi de mobilite en m-2
rau_screw = rau_screw * InvVol
rau_edge = rau_edge * InvVol

!*** Densite des dislocation mobile moyennee sur 4 iterations (m-2)
RAUDMO = (RAUDMO*TROIS+RAUDMO_I * InvVol)*QUART

Eps_el =  sigma * half / (un + Dpoiss)
Eps_pl = zero
Eps_pl_inst = zero
raudisact = zero
raudisjonc = zero
raudisjonc_precipite = zero

if (desorientgrain) then
  raugrainjonc(:)=zero
  rausysjoncgrain(:,:)=zero
endif

do i = 1, NTSG

    ! We start to calculate plastic deformation only after the relaxation steps
    if (kk <= relax_reac) AireSys(i) = zero

    if (kkdebug) write(379,*) "Junction length for Sys =",I," is ", RauSysjonc(I)
    RauSys(I) = RauSys(I)* InvVol
    RauSys_precipite(I) = RauSys_precipite(I)* InvVol
    RauSysjonc(I) = RauSysjonc(I)* InvVol
    raudisjonc = raudisjonc + RauSysjonc(I)

    ! Calculation for the multi-grain simulation
    if (DesorientGrain) then
      RauSysjoncgrain(:,I) = RauSysjoncgrain(:,I)* CONST3
      raugrainjonc(:) = raugrainjonc(:) + RauSysjoncgrain(:,I)
    endif

    if(solli_sys(i) > 0.99) raudisact = raudisact +  RauSys(I)

    ! The calculation of airevis and airecoin must be signed to be correct
    aireVis   = aireVis + AireVisSys(i)
    AireCoin  = AireCoin + AireCoinSys(i)
    AireMixte = AireMixte + AireMixtSys(i)

    if (DesorientGrain) then

      do j = 1,nbgrains

        ! strain in the tensile axis direction
        Eps_pl = Eps_pl + AireSysGrain(j,i) * SchmidSysGrain(j,i) * BdivA / Volume
        ! Strain rate as no meaning in multi-grains simulation
        Eps_pl_inst = 0

      enddo

    else

      ! strain in the tensile axis direction
      Eps_pl = Eps_pl + AireSys(i) * SchmidSys(I) * BdivA / Volume
      ! strain increment in the tensile axis direction
      Eps_pl_inst = Eps_pl_inst + AireSysInst(i) * SchmidSys(I) * BdivA / Volume

    endif

    ! calcul du travail applique cumule par system
    TrAppSys(i) = TrAppSys(i) + TrAppSysInst(i)

    ! calcul du travail applique total instanane
    TrAppInst = TrAppInst + TrAppSysInst(i)

    ! calcul du travail interne cumule par system
    TrIntSys(i) = TrIntSysInst(i) + TrIntSysInst(i)

    ! calcul du travail applique total instanane
    TrIntInst = TrIntInst + TrIntSysInst(i)
enddo

#ifdef PA
!*** Calculation of the total plastic strain in the tensil direction for a multi-phase simulation ***
if (Nb_phase == 2) then
  ! The total deformation is due to the displacement of the dislocations contained in all the phases
  ! The strain contribution of each phase is so added
  ! The aim of the calculation is to provide a correct control on the applied strain rate

  ! Initialization
  Eps_pl_mphase = zero ! Total plastic deformation in the tensil direction in case of a multi-phase simulation


  ! SUM on all the phases
  ! a MPI_ALLREDUCE will be used, rather than a combination MPI_REDUCE between the header proc plus MPI_BCAST, because
  ! the second solution may need a MPI_BARRIER to work.
  ! As a consequence each proc will add its contribution to the sum. As a preliminary, we need so to divide the value
  ! of this contribution by the number of procs of one phase
  temp_phase = Eps_pl
  temp_phase = temp_phase / TAILLE_PAIRIMPAIR

  CALL MPI_ALLREDUCE(temp_phase,Eps_pl_mphase,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

endif
#endif


!*** Densite totale des disloctions en m-2
RAUDIS = RAUDIS * InvVol - half * raudisjonc
RAUDIS_precipite = RAUDIS_precipite * InvVol - half * raudisjonc_precipite
RAUDIS_matrice = RAUDIS - RAUDIS_precipite
sysact = IUN    ! index of the most active system

do i = 1, NTSG
!    if (schmidsys(I) < Schmid .and. rausys(i) > UN) print *, "changt 1 de sysact de ",i,"a",sysact
    if (schmidsys(I) < Schmid .and. rausys(i) > UN) sysact = I
enddo

do i = 1, NTSG
!    if ( abs(AireSys(i)) >  abs(AireSys(sysact)) .and. schmidsys(i) >= schmidsys(sysact)) &
!        print *, "changt 1 de sysact de ",i,"a",sysact
    if ( abs(AireSys(i)) >  abs(AireSys(sysact)) .and. schmidsys(i) >= schmidsys(sysact)) sysact = I
enddo

!*** meme chose pour chaque grain par syst de glissement
if (DesorientGrain) RAUGRAINSYS(:,:) = RAUGRAINSYS(:,:) * CONST3 - half * rausysjoncgrain(:,:)


!***cas des barrieres spheriques, il faut modifier la densite
if (GB==IDEUX) RAUDIS=RAUDIS*factspheric

! pour le HCP avec particules, on normalise la vitesse de deformation pour que les simulations
! fonctionnent a 10 12 m-2 de densite de dislocation
! la variable temp donne la valeur correspondant a la vitesse de deformation conseillee
! il s'agit de calculer la vitesse resultante du mouvement au moyenne d'un pas par segment
! et par iteration

!if(npar > izero .and. HCP ) then
if(relative_rate) then
   if(rau_ini < zero) then
      rau_ini = raudis
   else
      epsilonpoint = strainRate_ini * raudis / rau_ini
      if(modulo (kk-1, 10*kpredraw) == izero) &
       write(*,'("SR_ini=",F8.4,"s-1; Rau_ini =",F8.4,"E+12 ====> SR new =", F8.2," s-1; rau =",F8.2,"E+12;")') &
      strainRate_ini, rau_ini*1E-12,epsilonpoint,raudis*1E-12
   endif
endif

!write(*,'(" La vitesse de deformation conseillee = ", F15.10," s-1")')  temp
! la deformation total = eps elastique + eps plastique(eps_pl)
! le modeule de young egal E = 2(1+nu)G donc
!if (.not. HCP) then
!   call mes("Attention :  EPSO, desormais inclut la partie elastique. !!")
!endif
!print *,  sigma, Eps_pl ,Eps_el

!if (npar > izero) then
EPSO =  Eps_pl
!else
!   EPSO =  Eps_pl + Eps_el
!endif
! facteur a rajouter pour le calcul de deformation de la sphere ca n est
! pas le meme volume
if (GB==IDEUX) EPSO=EPSO*factspheric
! attention: le signe est arbitraire mais la valeur absolu est correcte

! Plastic fields are recalculated at each steps
TensRot(1:3,1:3)=zero
TensEps(1:3,1:3)=zero
TensDist(1:3,1:3)=zero

Do i = 1, NTSG
    ! Calculation of gamma for each slip system
    gamSys(i) = aireSys(i) * BdivA / Volume
    ! Total rotation tensor of the simulation box
    do j=1,3
        do k=1,3
            TensRot(j,k) = TensRot(j,k) + gamSys(i)*TensRotSys(i,j,k)
            TensEps(j,k) = TensEps(j,k) + gamSys(i)*TensEpsSys(i,j,k)
            TensDist(j,k) = TensDist(j,k) + gamSys(i)*TensDistSys(i,j,k)
        enddo
    enddo
enddo

! if sysact is not defined
!*************************************
! The multi-grains simulation additional calculations
if (DesorientGrain) then
  do jj=1,nbgrains
      do i=1,NTSG
          gamsysGRAIN(jj,i)= aireSysGRAIN(jj,i) * BdivA / Volume * nbgrains ! volume d un grain
          do j=1,3
              do k=1,3
                  TensRotGrain(jj,j,k) = TensRotGrain(jj,j,k)+gamsysgrain(jj,i)*TensRotSysGrain(jj,i,j,k)
              enddo
          enddo
      enddo
  enddo
endif
!*************************************

! calcul du nombre moyen de crans formes
Do i = 1, NTSG
    ! calcul de la densite de dislocation de foret pour le system I
    temp = zero
    if (RauSys(I) > zero) then
      rauForet(I) = RauDis - RauSys(I)
      !     calcul du nombre de crans forme/unite de longueur de dislocation /iteration
      ! il s'agit de diviser le nombre total de crans forme sur chaque systeme:
      !                  airesysinst(i) * avalue**2 * rauforet
      ! par la longueur totale de dislocations sur le systeme en question:
      !                  RauSys(i) * VOLUME * avalue**3
      ! ce qui fait :
      temp = (AireSysInst(I) * RauForet(I) )/(RauSys(I)*VOLUME*avalue)
      ! on ajoute temp aux crans deja forme
      if (kkdebug) then
        write(379,*) "Sys =", I, " cran/uni  =", temp
        write(379,*) "RauForet =",  RauForet(I)*1D-10, " RauSys=", RauSys(I)*1D-10
        write(379,*) "Swept area at this step : ", AireSysInst(I)
        write(379,*) "Sys length  =", rausys(i)* VOLUME * avalue**3
      endif
      CranSys(I) = CranSys(I) + temp
    endif
enddo
! if Sysact is not loaded to the maximum Schmid factor, correct it.
! vitesse de cisaillement par systeme
GAMMADOTSYS(1:NTSG) = GamSys(1:NTSG) / deltat

! calcul du travail applique total cumule
TrApp = TrApp + TrAppInst

! calcul du travail interne total cumule
TrInt = TrInt + TrIntInst

EPSOLD(1:(NstatControl-1)) = EPSOLD(2:NstatControl)

if (Nb_phase /= 2) then
  EPSOLD(NstatControl) = EPS_pl / deltat
else
#ifdef PA
  ! Recall: Nb_phase =2 is only possible in parallel mode
  ! Multi-phase simulation: the control is made on the total plastic strain in the loading axis direction
  EPSOLD(NstatControl) = EPS_pl_mphase / deltat
#endif
endif

!---------------------------------------------------------------------------
! Calculation of the plastic deformation at each step is too fluctuating and
! smoothing procedure are needed to correctly monitor the simulation

! ! Calculation made of two simple moving average
! EPSDOT = (SUM(EPSOLD((NstatControl/2+1):NstatControl))-             &
!           SUM(EPSOLD(1:(NstatControl/2))))/(NstatControl**2*0.25D0)

! WMA = Espdot is calculate with a standard weighted moving average
EPSDOT = zero
do i = 1, NStatControl-1
  EPSDOT = EPSDOT + i * (EPSOLD(i+1) - EPSOLD(i))
enddo
EPSDOT = EPSDOT / (half*NStatControl*(NStatControl-1))   ! Average is made with (NStatControl-1) elements


do sigsursys = 1, ns, 1

!   ! Calculation made of two simple moving average
!   GAMMADOTSYS(sigsursys) = (SUM(GAMMAOLDSYS(sigsursys,(NstatControl/2+1):NstatControl))-             &
!                             SUM(GAMMAOLDSYS(sigsursys,1:(NstatControl/2))))/(NstatControl**2*0.25D0)

  ! WMA = GAMMADOTSYS is calculate with a standard weighted moving average
  GAMMADOTSYS(sigsursys) = zero
  do i = 1, NStatControl-1
    GAMMADOTSYS(sigsursys) = GAMMADOTSYS(sigsursys) + i * (GAMMAOLDSYS(sigsursys,i+1) - GAMMAOLDSYS(sigsursys,i))
  enddo
  GAMMADOTSYS(sigsursys) = GAMMADOTSYS(sigsursys) / (half*NStatControl*(NStatControl-1))       ! Average is made with (NStatControl-1) elements

enddo

!*** -EF- *** Sauvegarde pour les EF tous les KKEF
! if(MOD(KK,KKEF).eq.IZERO) then
!  call WRITEGAUSS(82)
! endif
!*** -EF- ***

!*** Sauvegarde des grandeures utiles tous les "Kstat" iterations
!****** DEBUT DE STAT OUT *******
ll = izero
Do i = 1 , nsegm
    if(seg(i)%norme == izero .and. .not. seg(i)%jonc .and. seg(i)%gd < iun) ll = ll + 1
enddo
ll = int(real(ll,DP) / real(nsegm,DP) * 100)

if (NTSG > nsmax) then
  print *, "STOP,kk=",kk, " stat.txt can be written only if NTSG < 13"
  stop
endif

if(MOD(KK,KSTAT).eq.IZERO) then

!*** Total density of long segments
RAU_GROSEG = RAU_GROSEG *InvVol
!*** end of the calculation of the long segment density ***

  !   call DATE_AND_TIME(DATE=DAT,TIME=TIM)

#ifdef PA
  if ((Mon_Rang_PairImpair + IUN) >  IUN) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  else
#endif

    if(rauforet(sysact) < UN) rauforet(sysact) = IUN
    const2 = sigma * xmu * Schmid   ! maximum shear stress
    ! const3 = (sigma*xmu*schmidsys(sysact)- loi(2)%friction)/xmu/vecBurgers/DSQRT(RauForet(sysact))
    const3 = (const2 - loi(2)%friction)/xmu/vecBurgers/DSQRT(Raudis)
    ! 11 colonnes

27  format(I8,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,F11.7,1x,F11.7,1x,F11.7,1x,F11.7,1x,F20.15)

    open( 8,FILE=file_graph,STATUS='OLD',POSITION='APPEND')
          write(8,27) kk, &            ! number of simulation step
          (Eps_el+Eps_pl)*100, &  ! total deformation
          Eps_pl*100, &           ! plastic deformation
          sigma * xmu * 1D-6 , &  !  the stress in Pa
           Eps_pl*100/Schmid, &   !Gamma
          sigma * xmu * 1D-6 * Schmid , & ! maximum resolved shear stress (in Pa)
          epsdot  , &             ! strain rate
          raudis  , &             ! dislocation density
          (sigma * xmu * Schmid - (loi(2)%friction*sign(un,sigma)))/xmu/vecBurgers/DSQRT(Raudis),& !alpha
          ZPrime(1:3),accutime    ! the present tensile axis and total time of simulation
    close (8)

    open(16,FILE=file_eps,STATUS='OLD',POSITION='APPEND')
    write(16,fmt='(I14 ,12E14.3)') kk,TensEps(1,1),TensEps(2,2),TensEps(3,3),TensEps(2,3),TensEps(1,3),TensEps(1,2),  &
                                   SigApp(1,1),SigApp(2,2),SigApp(3,3),SigApp(2,3),SigApp(1,3),SigApp(1,2)
    close (16)

216 format(2X,E13.5,2x,I5,2x,I5,2x,I5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,   &
            2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,I5,2x,I5,2x,I5,2x,I5,2x,I5,2x,I5)

    open( 99,FILE=file_stat,STATUS='OLD',POSITION='APPEND')
    write(99,216) Eps_pl*100,         & ! plastic deformation
        nbjonction , ntgd ,           & ! total NUMBER of junctions and gd
        ll ,                          & ! percentage of kneecaps  over total number of segment
        airevis,airecoin,airemixte,   & ! swept areas by screw and edge dislocation
        rau_screw, rau_edge,          & ! dislocation densities for mobility laws 1 and 2
        L_average_screw,              & ! average "effective" length of screw segments
        L_average_edge,               & ! average length of non-screw segments
        tau_average_screw,            & ! average shear stress on screw segments
        tau_average_edge,             & ! average shear stress on non-screw segments
        frankloop*InvVol,             & !< Density of loop
        frankline*InvVol,             & !< Density of infinit lines
        frankfreeline*InvVol,         & !< Density of free infinit lines (not pinned)
        LoopNumb,                     & !< Number of loops
        FreeLoopNumb,                 & !< Number of infinit free lines
        ntgd1,                        & !< Number of GD segments of type 1 = made by colli annihilation
        ntgd2,                        & !< Number of GD segments of type 2 = made by cross-slip
        oldntgd1,                     & !< Number of GD segments of type 1 = made by colli annihilation and eliminated
        oldntgd2                        !< Number of GD segments of type 2 = made by cross-slip and eliminated

    close (99)

    open(25,FILE=file_gammap,STATUS='OLD',POSITION='APPEND')
    open(28,FILE=file_rau,STATUS='OLD',POSITION='APPEND')
    open(29,FILE=file_raujonc,STATUS='OLD',POSITION='APPEND')
    open(38,FILE=file_gamma,STATUS='OLD',POSITION='APPEND')
    open(48,FILE=file_travapp,STATUS='OLD',POSITION='APPEND')
    open(49,FILE=file_travint,STATUS='OLD',POSITION='APPEND')

    if ((.not. allocation_dynamique_boites) .and. calculate_gammabox) then
      open(61,FILE=file_gammabox,FORM='UNFORMATTED',STATUS='OLD',POSITION='APPEND')

      ! The grid dimension is redefined and is used as header for the gammabox file
      if (KK == KSTAT) then
        write(61)  int(NBoites,DPI_S)
        write(61)  int(NBoitesX,DPI_S),tailleboite(1)
        write(61)  int(NBoitesY,DPI_S),tailleboite(2)
        write(61)  int(NBoitesZ,DPI_S),tailleboite(3)
        write(61)  int(NTSG,DPI_S)
      endif
    endif

    !* multicristals outfiles
    if (DesorientGrain) then
      open(47,FILE='../out/polydens.txt',STATUS='OLD',POSITION='APPEND')
      open(46,FILE='../out/polyrotation.txt',STATUS='OLD',POSITION='APPEND')
      open(45,FILE='../out/polygamma.txt',STATUS='OLD',POSITION='APPEND')
    endif

    ! write(*,'("Vv", F10.5, " m/s SRc ", F10.5 "106 /s,  L = ",F10.5," nm")') &
    !     eps_pl_inst/deltat/schmid/vecburgers/raudis,55*vecburgers*raudis*schmid/1000000,&
    !     raudis*volume*1000000000*avalue**3

    ! 13 colonnes
37  format(1X,13(E11.4,1X))

    toto(1:nsmax) = zero  ! This initialization is needed to set at zero columns between ns+ 1 and nsmax
    itoto(1:nsmax) = izero

    toto(1:ns) = abs(GAMMADOTSYS(1:ns))
    write(25,37) Epso*100, (toto(JK),JK=1,nsmax)

    toto(1:ns) = RAUSYS(1:ns)
    write(28,'(1X,22(E11.4,1X))') Epso*100, RAU_GROSEG, (toto(JK),JK=1,nsmax)

    ! 13 colonnes
337 format(1X,13(E11.4,1X),12I5)

    toto(1:ns) = RAUSYSjonc(1:ns)
    itoto(1:ns) = njoncsys(1:ns)
    write(29,337) Epso*100, (toto(JK),JK=1,nsmax), (itoto(JK),JK=1,nsmax)

    ! 13 colonnes
    toto(1:ns) = GAMSYS(1:ns)
    write(38,37) Epso*100, (toto(JK),JK=1,nsmax)

    ! 13 colonnes
    toto(1:ns) = TrAppSys(1:ns)
    write(48,37) Epso*100, (toto(JK),JK=1,nsmax)

    ! 13 colonnes
    toto(1:ns) = TrIntSys(1:ns)
    write(49,37) Epso*100, (toto(JK),JK=1,nsmax)

    ! if needed, we save the gammabox file
    if ((.not. allocation_dynamique_boites) .and. calculate_gammabox) then

      ! Ecriture des gamma par systeme dans chaque boite, 13 colonnes par boites
      ! Reecriture du header a chaque pas de sauvegarde
      if (MOD(KK,KSTAT) .eq. zero) then
        ! The simulation time step
        write(61) int(kk,DPI_S)

        ! In order to check the gammabox calculation, the total gamma is calculated
        somme(:)=0
        do jsum = 1,NTSG
          do isum = 1,NBoites
            somme(jsum) = somme(jsum) + Gammabox(isum,jsum)
          enddo
        enddo
        write(61) (somme(JK) * volboite / volume, JK=1,NTSG)

        ! The gammabox field is saved
        do ni = 1,NBoites
          Px = (real(B1D_3D(ni,1),DP) - half) * tailleBoite(1)
          Py = (real(B1D_3D(ni,2),DP) - half) * tailleBoite(2)
          Pz = (real(B1D_3D(ni,3),DP) - half) * tailleBoite(3)
          write(61) Px,Py,Pz,(Gammabox(ni,JK), JK=1,NTSG)
        enddo

        ! This mark is used to end the step info
        write(61) '&'
      endif

    endif

    ! Calculation for the multi-grains simulation
36  format(1X,17(E11.4,1X))
    if (DesorientGrain) then

      !**
      do ii=1,nbgrains
          write(47,36) Epso*100,raugrainsys(ii,:)
      enddo
      write(47,*)

      !**
      do ii=1,nbgrains
          do ll=1,3
              write(46,*) TensRotGrain(ii,ll,:)
          enddo
          write(46,*)
      enddo
      write(46,*)

      !**
      do ii=1,nbgrains
          write(45,36) sigma*xmu*1.d-6,gamsysgrain(ii,:)
      enddo
      write(45,*)

    endif

    close(25)
    close(28)
    close(29)
    close(38)
    close(48)
    close(49)
    close(47)
    close(46)
    close(45)
    close(44)
    if ((.not. allocation_dynamique_boites) .and. calculate_gammabox) close(61)

    if (nsegm+1000 > nsegmax) then
      write (*,*) "NSEGMAX is too small : ",nsegm,'/',nsegmax
    endif

#ifdef PA
    ! Le proc zero a fini d ecrire, on libere tout le monde
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
#endif

endif
!****** FIN STAT OUT ********

end subroutine STATIS

!###########################################################################
!#  Calculation of the number of segments junction and gd                  #
!##############################################################22/12/98#####

subroutine reacloc

implicit none

integer(DPI)  :: I,J,ISYST

nbjonction    = izero
ntgd          = izero
ntgd1         = izero
ntgd2         = izero
NTFS          = izero
njoncsys(:)   = izero

! Calculation of the number of GD and junction segments
do i = 1,nsegm
    if (seg(i)%voiso /= i) then
      if (seg(i)%gd > izero) then
        ntgd = ntgd + 1
        if (seg(i)%gd > iun) then
          ! The number of segment GD made by cross-slip
          ntgd2 = ntgd2 + 1
          if (kkdebug) write(379,*) 'In reacloc - 1, ntgd2 after icrement =', ntgd2,'with seg i=',I
        else
          ! The number of segment GD made by colli annihilation
          ntgd1 = ntgd1 + 1
        endif
      endif

      if (seg(i)%jonc) then
        j = seg(i)%Ijonc
        if(seg(j)%jonc) then
          nbjonction = nbjonction + 1
          isyst=syseg(seg(i)%veclin)
          njoncsys(isyst)=njoncsys(isyst)+1
        else
          print *, " Reacloc: Jonc assyme entre i=",i," et j :", j, "kk =",kk
          seg(i)%jonc = .false.
          seg(i)%Ijonc = NSEGMAX
          seg(j)%Ijonc = NSEGMAX
          stop
        endif
      endif
      if (seg(i)%surface > IZERO) then
        NTFS = NTFS +1
      endif
    endif
enddo

!! A paranoid test
!if (GB == izero .and. modulo((ntgd2+oldntgd2),2) /= 0) then
!  write(*,*) 'At step ', KK, 'ntgd2 + oldntgd2 is found assymetric', ntgd2 + oldntgd2, ntgd2, oldntgd2
!  stop
!endif

if (nbjonction /= 0) then
  if (kk_debut_jonc == 0) kk_debut_jonc = kk
  KKjonc = KKjonc + 1
endif

if (ntgd /= 0) then
  if (kk_debut_GD == 0) then
    kk_debut_GD = kk
    kkGD = 1
  else
    KKGD = KKGD + 1
  endif
endif

if (modulo(nbjonction,ideux) /= 0) print *, " reacloc :  Impair number of junction segments ?"

end subroutine reacloc

!###########################################################################
!#  Calcul de la surface entre la ligne de dislo et l'axe de jonction      #
!##############################################################22/12/98#####
subroutine stat_carto

implicit none

integer(DPI)  :: i,ii,J,jj,dislo1, Oi(3),vec(3),dislo2,vecl(3),facdiv
integer(DPI)  :: Ojonc(3),Ejonc(3),njonc,ljonc
real(DP)      :: temp,distance_carto,d(4)
logical :: ouvert,test,nulle,colin


facdiv = 1
if(HCP) facdiv = 22

if (cartograph == 0) print *, "STOP,kk= ",kk, " cartograph /= 1 : stat_carto not defined"
test = .false.
colin = .false.
if(sysconnec(axe1_carto,axe2_carto) == 2) colin = .true.

if (kk == Nstep) test = .true.

dislo1 = 0  ! debut de la dislo 1
dislo2 = 0  ! debut de la dislo 2
do i = 1,nsegm
    if(seg(i)%voiso == NSEGMAX) then                   !!***
      if (dislo1 == 0 .and. dislo2 == 0) then
        dislo1 = I
      elseif(dislo1 /= 0 .and. dislo2 == 0) then
        dislo2 = I
        exit
      endif
    endif
enddo

if(dislo1*dislo2 == 0) print *, "STOP,kk=",kk, " dislo1 or dislo2 is not defined"

! calcul de la distance minimale entre les deux dislocation : calcul lourd sur n**2
! l'idee est de ballayer point par point la premiere dislo et pour chaque point calculer
! la distance minimal avec l'autre dislo

i = dislo1
distance_carto = modur(1)
ouvert = .false.
njonc = 0
Ojonc(:) = -1 ! debut de la jonction
Ejonc(:) = -1 ! fin de la jonction

BI : Do while (i /= NSEGMAX)
    BII : Do ii = 0, seg(i)%norme-1
        Oi(:) = seg(i)%O(:) + bveclin(:,seg(I)%veclin)*ii
        j = dislo2
        nulle = .false.
        BJ : Do while(j /= Nsegmax)
            BJJ : Do JJ = 0, seg(J)%norme-1
                vec(:) =  (Oi(:)-(seg(j)%O(:)+bveclin(:,seg(J)%veclin)*jj))
                temp = norivect(vec)
                if(temp < distance_carto) distance_carto = temp
                if(NTGD == 0 .and. (abs(vec(1)) + abs(vec(2)) + abs(vec(3)) == 0)) then
                  if (Ojonc(1) == -1) then
                    Ojonc(:) = Oi(:)
                    Ejonc(:) = Oi(:)
                  else
                    Ejonc(:) = Oi(:)
                  endif
                endif
            enddo BJJ
            j = seg(j)%voise
        enddo BJ

    enddo BII
    i = seg(i)%voise
enddo BI

! calcul de longueur de recouvrement des deux dislocation

ljonc = 0
if(.not. colin) then
  vec(:) = (Ejonc(:) - Ojonc(:))/facdiv
  ljonc = nint(facdiv*norivect(vec)/normlin(axe1_carto),DPI)
endif


if(NTGD /= 0) then
  dislo1 = 0
  dislo2 = 0
  do i = 1,nsegm
      if(seg(i)%gd > izero) then
        if (dislo1 == 0 .and. dislo2 == 0) then
          dislo1 = I
        elseif(dislo1 /= 0 .and. dislo2 == 0) then
          dislo2 = I
          exit
        endif
      endif
  enddo
  d(1) = norivect(seg(dislo2)%O(:) - seg(dislo1)%O(:))
  vec(:) =  seg(dislo2)%O(:)+ seg(dislo2)%norme*bveclin(:,seg(dislo2)%veclin)
  vecl(:) =  seg(dislo1)%O(:)+ seg(dislo1)%norme*bveclin(:,seg(dislo1)%veclin)
  d(2) = norivect(vec(:) - seg(dislo1)%O(:))
  d(3) = norivect(seg(dislo2)%O(:) - vecl(:))
  d(4) = norivect(vec(:) - vecl(:))
  distance_carto = 1000000
  do i = 1,4
      if(d(i) < distance_carto) distance_carto = d(i)
  enddo
endif

!On normalise par la taille du vecteur de l'axe de jonction
distance_carto = distance_carto / normlin(axe1_carto)

if (kk == relax_TL .and. cartograph == Ideux) then
  dist_carto_ini = distance_carto
  print *, " Initial spacing between sources before reaction =", dist_carto_ini
endif

10 format (F10.2,F10.2,"  -1000.00  -1000.00",   F10.1,F10.1,F10.2,"% ",A12) ! pour jonction
20 format (F10.2,"  -1000.00",F10.2,"  -1000.00",F10.1,F10.1,F10.2,"% ",A12) ! etat croise
30 format (F10.2,"  -1000.00  -1000.00",F10.2,   F10.1,F10.1,F10.2,"% ",A12) ! repulsif

1 format (F10.2,F10.2,"  -1000.00  -1000.00",    F10.1,F10.1,F10.2,"% ",I5,I6,A12) ! pour jonction
2 format (F10.2,"  -1000.00",F10.2,"  -1000.00",F10.1,F10.1,F10.2,"% ",I5,I6,A12) ! etat croise
3 format (F10.2,"  -1000.00  -1000.00",F10.2,   F10.1,F10.1,F10.2,"% ",I5,I6,A12) ! repulsif

if(kk == NSTEP) then

  open(9,FILE='../out/stat.carto',STATUS='OLD',POSITION='APPEND')

  if(colin) then
    temp = 100.0*dble(kkGD)/dble(kk - kk_debut_GD)
    print *, kkGD,kk_debut_GD," KK......"
    if(ntgd > 0 .and. distance_carto >= ideux) then
      write(9,10) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,fichier_carto
      write(*,10) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,fichier_carto
    elseif(ntgd == 0 .and. distance_carto > dist_carto_ini) then
      write(9,30) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,fichier_carto
      write(*,30) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,fichier_carto
    else
      write(9,20) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,fichier_carto
      write(*,20) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,fichier_carto
    endif
    !   on considere une jonction quand la taille de jonction depasse 10
  else
    temp = 100.0*dble(kkJONC)/dble(kk - kk_debut_JONC)
    if(ljonc >= ideux) then
      write(9,1) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,nbjonction,ljonc,fichier_carto
      write(*,1) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,nbjonction,ljonc,fichier_carto
      !   on considere etat repulssif quand dist_carto_ini < distance final
    elseif(distance_carto > dist_carto_ini) then
      write(9,3) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,nbjonction,ljonc,fichier_carto
      write(*,3) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,nbjonction, ljonc,fichier_carto
    else
      write(9,2) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,nbjonction, ljonc,fichier_carto
      write(*,2) phi1_carto,phi2_carto,dist_carto_ini,distance_carto,temp,nbjonction, ljonc,fichier_carto
    endif
  endif
  close(9)
endif

end subroutine stat_carto


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Subroutine where the main outputs of MDC simulations are calculated !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MDC
subroutine STATISMDC

implicit none

real(kind=DP) :: UNDELTA      !<
real(kind=DP) :: XLONG        !<
real(kind=DP) :: InvVol       !< 1 / Volume (a unit)
real(kind=DP) :: const2       !<
real(kind=DP) :: const3       !<
real(kind=DP) :: temp         !<
real(kind=DP) :: factspheric  !<

real(kind=DP) :: SONseuil , rauforet(NTSG), raudisact
real(kind=DP) :: somme(NTSG)    !< To check the calculation of gammabox
real(kind=DP) :: RAU_GROSEG     !< Total density of long segments at step kk
real(kind=DP) :: Px,Py,Pz       !< Used to write the gammabox file

integer(kind=DPI)   :: I,Li,ISYST,nloi,JK,sigsursys,ns,ii,j,jj,ll,k,sysact
integer(kind=DPI)   :: ni,isum,jsum   ! Used to write the gammabox file

integer,parameter   :: nsmax = 12       ! Maximum number of slip systems considered in the output files
real(kind=DP)       :: toto(nsmax)      ! Intermediate tabs needed for simple writing
integer(kind=DPI)   :: itoto(nsmax)

#ifdef PA
real(kind=DP)       :: Eps_pl_mphase, &! Total plastic deformation in the loading direction in case of a multi-phase simulation
                       temp_phase
#endif

!*************************************************************************
!*** Reinitialisation des grandeurs instantannees des calcul de STATIS ***
!*************************************************************************
!Schmid = 0
!SchmidSys(:) = 0

! RAUDMO,VITMOY,EPSDOT,EPS_PL,airevis,airecoin,airemixte

RAUDIS = 0.d0        !*** Densite totale des disloctions
RAUDIS_matrice = 0.d0      !*** Densite par systeme dans la matrice
RAUDIS_precipite = 0.d0    !*** Densite par systeme dans le precipite
RAUDMO_I = 0.d0      !*** Densite des dislocation mobile moyennee sur 4 iterations (m-2)
RAUSYS(:) = ZERO     !*** Densite par systeme
RAUSYS_precipite(:) = ZERO     !*** Densite par systeme
RAU_GROSEG = ZERO   ! Density of segments not treated with the multipole approach
NSEGMJONC = 0        !*** Compteur dislo mobiles
VITA = 0.0D0         !*** Accumulateur pour la vitesse moyenne des dislo
GAMMADOTSYS(:)=0.D0
ns = NTSG
InvVol = 1.0/(VOLUME*AVALUE*AVALUE)
const3 = zero
factspheric = zero

UNDELTA = UN/DELTAT*AVALUE
SONseuil = SONDMS*0.95
TrAppInst = zero
TrIntInst = zero

! The multi-grains extra calculations
if (desorientgrain) then
  RAUGRAIN(:)=0.d0          ! Dislocation density in the grains
  RAUGRAINSYS(:,:)=0.d0     ! Slip system density in the grains
  TensRotGrain(:,:,:)=0.d0
  gamsysGRAIN(:,:)=0.d0
  const3 = InvVol*nbgrains  ! Volume of the grains (supposed to be regular)
endif

if (GB==IDEUX) factspheric =  VOLUME*AVALUE**3/(QUATRE/TROIS*PII*(GB2_size*1.D-6)**3)

!*******************************
!*** Boucle sur les segments ***
!*******************************
do I = 1,NSEGM
    Li = SEG(I)%VECLIN
    XLONG = SEG(I)%NORME*NormLin(Li)
    ISYST = SYSEG(Li)
    !*** Cumule des longeurs pour la densite totale de dislocation
    RAUDIS = RAUDIS+XLONG

    !*** Calculation of the total density of long segments   ***
    if(MOD(KK,KSTAT).eq.IZERO) then
      if(seg(I)%norme * normlin(seg(I)%veclin)  > lengthlimit) then
          RAU_GROSEG = RAU_GROSEG + Xlong
      endif
    endif

    if (seg(i)%NPhase == IUN) RAUDIS_precipite = RAUDIS_precipite+XLONG

    !***Cumul pour chaque grain
    if (DesorientGrain) then
      call CorrectifGrain(I,3005) !verif necessaire les initialisations ne sont pas tres stables
      RAUGRAIN(seg(i)%grain)=RAUGRAIN(seg(i)%grain)+xlong
    endif
    !*** 08/11/00 : STAT SUR LES JONCTIONS
    if(seg(i)%jonc) then
      NSEGMJONC = NSEGMJONC+1
    endif
    !*** Statistique sur la densite par systeme de glissement
    RAUSYS(ISYST) = RAUSYS(ISYST)+XLONG

    if (seg(I)%NPhase ==iun) RAUSYS_precipite(ISYST) = RAUSYS_precipite(ISYST)+XLONG

    !*** Statistique sur la densite par systeme de glissement pour chaque
    !*** grain
    if (DesorientGrain) RAUGRAINSYS(seg(i)%grain,ISYST) = RAUGRAINSYS(seg(i)%grain,ISYST)+XLONG
enddo


!!************************************************************************************
!!************************************************************************************
!! Determination du deplacement maximal d'ecranter les deplacements excessives :
!
Do nloi = 1,NLV

  !************************************************************************************
  !* calcul des deplacement moyens : vitesse initiale (brute), vitesse ecrantee,VITESSE vraie
  !************************************************************************************

  ! calculation of a mean displacement per unit length (avalue)
  if (Lseg_V(nloi) > zero) D_vraie(nloi) = D_vraie(nloi) / LSEG_V(nloi)

  if (KK < Relax_TL) then
    Dlimite(nloi) = Dlimite_max
  else
    Dlimite(nloi) =  D_vraie(nloi) * facteur_depmax
    if (Dlimite(nloi) < Modep_max)   Dlimite(nloi) = Dlimite_max
    if (Dlimite(nloi) > Dlimite_max) Dlimite(nloi) = Dlimite_max
  endif

  if (kk == iterinfo) &
    write(397,'("nloi =", I4,"  Lseg=",F7.3,"  Lseg_V=",F7.3,"  Di =",F7.3,"  Dv =",F7.3," Dlim=",F7.2)') &
       nloi, Lseg(nloi)*avalue*1.D6,Lseg_V(nloi)*avalue*1.D6,D_brute(nloi),D_vraie(nloi),dlimite(nloi)

enddo

!************************************************************************************
!*** Vitesse moyenne en valeur absolue des dislocations (m/s)
if(RAUDMO_I.ne.ZERO) then
  VITMOY = (VITMOY*TROIS+VITA/RAUDMO_I)*QUART
endif

!*** Densite par loi de mobilite en m-2
rau_screw = rau_screw * InvVol
rau_edge = rau_edge * InvVol

!*** Densite des dislocation mobile moyennee sur 4 iterations (m-2)
RAUDMO = (RAUDMO*TROIS+RAUDMO_I * InvVol)*QUART

Eps_el =  sigma * half / (un + Dpoiss)
Eps_pl = zero
Eps_pl_inst = zero
raudisact = zero
raudisjonc = zero
raudisjonc_precipite = zero
if (desorientgrain) then
  raugrainjonc(:)=zero
  rausysjoncgrain(:,:)=zero
endif

do i = 1, NTSG

    ! We start to calculate plastic deformation only after the relaxation steps
    if (kk <= relax_reac) AireSys(i) = zero

    if (kkdebug) write(379,*) "Junction length for Sys =",I," is ", RauSysjonc(I)
    RauSys(I) = RauSys(I)* InvVol
    RauSys_precipite(I) = RauSys_precipite(I)* InvVol
    RauSysjonc(I) = RauSysjonc(I)* InvVol
    raudisjonc = raudisjonc + RauSysjonc(I)
    !******************************
    !** meme chose pour chaque grain
    if (DesorientGrain) then
      RauSysjoncgrain(:,I) = RauSysjoncgrain(:,I)* CONST3
      raugrainjonc(:) = raugrainjonc(:) + RauSysjoncgrain(:,I)
    endif
    !******************************
    if(solli_sys(i) > 0.99) raudisact = raudisact +  RauSys(I)

    aireVis  = aireVis + AireVisSys(i)
    AireCoin  = AireCoin + AireCoinSys(i)
    AireMixte  = AireMixte + AireMixtSys(i)

    if (DesorientGrain) Stop 'DesorientGrain mode impossible in MDC'

    ! strain in the tensil axis direction
    Eps_pl = Eps_pl + AireSys(i) * SchmidSys(I) * BdivA / Volume

    ! strain increment in the tensile axis direction
    Eps_pl_inst = Eps_pl_inst + AireSysInst(i) * SchmidSys(I) * BdivA / Volume

    ! calcul du travail applique cumule par system
    TrAppSys(i) = TrAppSys(i) + TrAppSysInst(i)

    ! calcul du travail applique total instanane
    TrAppInst = TrAppInst + TrAppSysInst(i)

    ! calcul du travail interne cumule par system
    TrIntSys(i) = TrIntSysInst(i) + TrIntSysInst(i)

    ! calcul du travail applique total instanane
    TrIntInst = TrIntInst + TrIntSysInst(i)
enddo

#ifdef PA
!*** Calculation of the total plastic strain in the tensil direction for a multi-phase simulation ***
if (Nb_phase == 2) then
  ! The total deformation is due to the displacement of the dislocations contained in all the phases
  ! The strain contribution of each phase is so added
  ! The aim of the calculation is to provide a correct control on the applied strain rate

  ! Initialization
  Eps_pl_mphase = zero ! Total plastic deformation in the tensil direction in case of a multi-phase simulation


  ! SUM on all the phases
  ! a MPI_ALLREDUCE will be used, rather than a combination MPI_REDUCE between the header proc plus MPI_BCAST, because
  ! the second solution may need a MPI_BARRIER to work.
  ! As a consequence each proc will add its contribution to the sum. As a preliminary, we need so to divide the value
  ! of this contribution by the number of procs of one phase
  temp_phase = Eps_pl
  temp_phase = temp_phase / TAILLE_PAIRIMPAIR

  CALL MPI_ALLREDUCE(temp_phase,Eps_pl_mphase,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

endif
#endif


!*** Densite totale des disloctions en m-2
RAUDIS = RAUDIS * InvVol - half * raudisjonc
RAUDIS_precipite = RAUDIS_precipite * InvVol - half * raudisjonc_precipite
RAUDIS_matrice = RAUDIS - RAUDIS_precipite
sysact = IUN    ! index of the most active system

do i = 1, NTSG
!    if (schmidsys(I) < Schmid .and. rausys(i) > UN) print *, "changt 1 de sysact de ",i,"a",sysact
    if (schmidsys(I) < Schmid .and. rausys(i) > UN) sysact = I
enddo

do i = 1, NTSG
!    if ( abs(AireSys(i)) >  abs(AireSys(sysact)) .and. schmidsys(i) >= schmidsys(sysact)) &
!        print *, "changt 1 de sysact de ",i,"a",sysact
    if ( abs(AireSys(i)) >  abs(AireSys(sysact)) .and. schmidsys(i) >= schmidsys(sysact)) sysact = I
enddo

!*** meme chose pour chaque grain par syst de glissement
if (desorientgrain) RAUGRAINSYS(:,:) = RAUGRAINSYS(:,:) * CONST3 - half * rausysjoncgrain(:,:)


!***cas des barrieres spheriques, il faut modifier la densite
if (GB==IDEUX) RAUDIS=RAUDIS*factspheric

! pour le HCP avec particules, on normalise la vitesse de deformation pour que les simulations
! fonctionnent a 10 12 m-2 de densite de dislocation
! la variable temp donne la valeur correspondant a la vitesse de deformation conseillee
! il s'agit de calculer la vitesse resultante du mouvement au moyenne d'un pas par segment
! et par iteration
!
!if(npar > izero .and. HCP ) then
if(relative_rate) then
   if(rau_ini < zero) then
      rau_ini = raudis
   else
      epsilonpoint = strainRate_ini * raudis / rau_ini
      if(modulo (kk-1, 10*kpredraw) == izero) &
       write(*,'("SR_ini=",F8.4,"s-1; Rau_ini =",F8.4,"E+12 ====> SR new =", F8.2," s-1; rau =",F8.2,"E+12;")') &
      strainRate_ini, rau_ini*1E-12,epsilonpoint,raudis*1E-12
   endif
endif
!
!write(*,'(" La vitesse de deformation conseillee = ", F15.10," s-1")')  temp
! la deformation total = eps elastique + eps plastique(eps_pl)
! le modeule de young egal E = 2(1+nu)G donc
!if (.not. HCP) then
!   call mes("Attention :  EPSO, desormais inclut la partie elastique. !!")
!endif
!print *,  sigma, Eps_pl ,Eps_el

!if (npar > izero) then
EPSO =  Eps_pl
!else
!   EPSO =  Eps_pl + Eps_el
!endif
! facteur a rajouter pour le calcul de deformation de la sphere ca n est
! pas le meme volume
if (GB==IDEUX) EPSO=EPSO*factspheric
! attention: le signe est arbitraire mais la valeur absolu est correcte
!
! Plastic fields are recalculated at each steps
TensRot(1:3,1:3)=zero
TensEps(1:3,1:3)=zero
TensDist(1:3,1:3)=zero
!
Do i = 1, NTSG
    ! Calculation of gamma for each slip system
    gamSys(i) = aireSys(i) * BdivA / Volume
    ! Total rotation tensor of the simulation box
    do j=1,3
        do k=1,3
            TensRot(j,k) = TensRot(j,k) + gamSys(i)*TensRotSys(i,j,k)
            TensEps(j,k) = TensEps(j,k) + gamSys(i)*TensEpsSys(i,j,k)
            TensDist(j,k) = TensDist(j,k) + gamSys(i)*TensDistSys(i,j,k)
        enddo
    enddo
enddo
! if sysact is not defined
!*************************************
!* calcul du gamma par grain (DesorientGrain) et calcul du tenseur de rotation
!* correspondant
if (DesorientGrain) then
  do jj=1,nbgrains
      do i=1,NTSG
          gamsysGRAIN(jj,i)= aireSysGRAIN(jj,i) * BdivA / Volume * nbgrains ! volume d un grain
          do j=1,3
              do k=1,3
                  TensRotGrain(jj,j,k) = TensRotGrain(jj,j,k)+gamsysgrain(jj,i)*TensRotSysGrain(jj,i,j,k)
              enddo
          enddo
      enddo
  enddo
endif
!*************************************
!
! calcul du nombre moyen de crans formes
Do i = 1, NTSG
!    ! calcul de la densite de dislocation de foret pour le system I
    temp = zero
    if (RauSys(I) > zero) then
      rauForet(I) = RauDis - RauSys(I)
!      !     calcul du nombre de crans forme/unite de longueur de dislocation /iteration
!      ! il s'agit de diviser le nombre total de crans forme sur chaque systeme:
!      !                  airesysinst(i) * avalue**2 * rauforet
!      ! par la longueur totale de dislocations sur le systeme en question:
!      !                  RauSys(i) * VOLUME * avalue**3
!      ! ce qui fait :
      temp = (AireSysInst(I) * RauForet(I) )/(RauSys(I)*VOLUME*avalue)
!      ! on ajoute temp aux crans deja forme
      if (kkdebug) then
        write(379,*) "Sys =", I, " cran/uni  =", temp
        write(379,*) "RauForet =",  RauForet(I)*1D-10, " RauSys=", RauSys(I)*1D-10
        write(379,*) "Swept area at this step : ", AireSysInst(I)
        write(379,*) "Sys length  =", rausys(i)* VOLUME * avalue**3
      endif
      CranSys(I) = CranSys(I) + temp
    endif
enddo
! if Sysact is not loaded to the maximum Schmid factor, correct it.
! vitesse de cisaillement par systeme
GAMMADOTSYS(1:NTSG) = GamSys(1:NTSG)/deltat
!
! calcul du travail applique total cumule
TrApp = TrApp + TrAppInst
!
! calcul du travail interne total cumule
TrInt = TrInt + TrIntInst
if(kkdebug) then
  write(379,*) "Epsold =",EPSold(:)
  write(379,*) "Epso =",EPSo
endif

if (mode_deformation_key .and. mode_deformation >= 0 .and. mode_deformation < 8) then

  EPSOLD(1:(NstatControl-1)) = EPSOLD(2:NstatControl)
  !
  !if (Nb_phase /= 2) then
    EPSOLD(NstatControl) = EPS_pl / deltat
  !else
  !#ifdef PA
  !  ! Recall: Nb_phase =2 is only possible in parallel mode
  !  ! Multi-phase simulation: the control is made on the total plastic strain in the loading axis direction
  !  EPSOLD(NstatControl) = EPS_pl_mphase / deltat
  !#endif
  !endif

  !---------------------------------------------------------------------------
  ! Calculation of the plastic deformation at each step is too fluctuating and
  ! smoothing procedure are needed to correctly monitor the simulation

  ! ! Calculation made of two simple moving average
  ! EPSDOT = (SUM(EPSOLD((NstatControl/2+1):NstatControl))-             &
  !           SUM(EPSOLD(1:(NstatControl/2))))/(NstatControl**2*0.25D0)

  ! WMA = Espdot is calculate with a standard weighted moving average
  EPSDOT = zero
  do i = 1, NStatControl-1
    EPSDOT = EPSDOT + i * (EPSOLD(i+1) - EPSOLD(i))
  enddo
  EPSDOT = EPSDOT / (half*NStatControl*(NStatControl-1))       ! Average is made with (NStatControl-1) elements

endif

do sigsursys = 1, ns, 1

!   ! Calculation made of two simple moving average
!   GAMMADOTSYS(sigsursys) = (SUM(GAMMAOLDSYS(sigsursys,(NstatControl/2+1):NstatControl))-             &
!                             SUM(GAMMAOLDSYS(sigsursys,1:(NstatControl/2))))/(NstatControl**2*0.25D0)

  ! WMA = GAMMADOTSYS is calculate with a standard weighted moving average
  GAMMADOTSYS(sigsursys) = zero
  do i = 1, NStatControl-1
    GAMMADOTSYS(sigsursys) = GAMMADOTSYS(sigsursys) + i * (GAMMAOLDSYS(sigsursys,i+1) - GAMMAOLDSYS(sigsursys,i))
  enddo
  GAMMADOTSYS(sigsursys) = GAMMADOTSYS(sigsursys) / (half*NStatControl*(NStatControl-1))       ! Average is made with (NStatControl-1) elements

enddo

!
!*** -EF- *** Sauvegarde pour les EF tous les KKEF
! if(MOD(KK,KKEF).eq.IZERO) then
!  call WRITEGAUSS(82)
! endif
!*** -EF- ***
!
!*** Sauvegarde des grandeures utiles tous les "Kstat" iterations
!****** DEBUT DE STAT OUT *******
ll = izero
Do i = 1 , nsegm
    if(seg(i)%norme == izero .and. .not. seg(i)%jonc .and. seg(i)%gd < iun) ll = ll + 1
enddo
ll = int(real(ll,DP) / real(nsegm,DP) * 100)

if (NTSG > nsmax) then
  print *, "STOP,kk=",kk, " stat.txt can be written only if NTSG < 13"
  stop
endif

if(MOD(KK,KSTAT).eq.IZERO) then

!*** Total density of long segments
RAU_GROSEG = RAU_GROSEG *InvVol
!*** end of the calculation of the long segment density ***

  !   call DATE_AND_TIME(DATE=DAT,TIME=TIM)
#ifdef PA
  if ((Mon_Rang_PairImpair + IUN) >  IUN) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  else
#endif

!    if(rauforet(sysact) < UN) rauforet(sysact) = IUN
     const2 = sigma * xmu * Schmid   ! maximum shear stress
     ! const3 = (sigma*xmu*schmidsys(sysact)- loi(2)%friction)/xmu/vecBurgers/DSQRT(RauForet(sysact))
     const3 = (const2 - loi(2)%friction)/xmu/vecBurgers/DSQRT(Raudis)
     ! 11 colonnes


341  format(I8,4x,I8,4x,I8)
     open(67,FILE='../out/debug/MDCstepstat.txt',STATUS='OLD',POSITION='APPEND')
     write(67,341) kk,NBDEP,NSEGM
     close(67)

27   format(I8,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,E13.5,1x,F11.7,1x,F11.7,1x,F11.7,1x,F11.7,1x,F20.15)

     open( 8,FILE=file_graph,STATUS='OLD',POSITION='APPEND')
          write(8,27) kk, &            ! number of simulation step
          (Eps_el+Eps_pl)*100, &  ! total deformation
          Eps_pl*100, &           ! plastic deformation
          sigma * xmu * 1D-6 , &  !  the stress in Pa
          Eps_pl*100, &   !Gamma
          sigma * xmu * 1D-6 * Schmid , & ! maximum resolved shear stress (in Pa)
          epsdot  , &             ! strain rate
          raudis  , &             ! dislocation density
          (sigma * xmu * Schmid - (loi(2)%friction*sign(un,sigma)))/xmu/vecBurgers/DSQRT(Raudis),& !alpha
          ZPrime(1:3),accutime    ! the present tensile axis and total time of simulation
     close (8)

216  format(2X,E13.5,2x,I5,2x,I5,2x,I5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,   &
            2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,E13.5,2x,I5,2x,I5,2x,I5,2x,I5,2x,I5,2x,I5)

     open( 99,FILE=file_stat,STATUS='OLD',POSITION='APPEND')
     write(99,216) Eps_pl*100,        & !< plastic deformation
          nbjonction , ntgd,          & !< total NUMBER of junctions and gd
          ll,                         & !< percentage of kneecaps  over total number of segment
          airevis,airecoin,airemixte, & !< swept areas by screw and edge dislocation
          rau_screw, rau_edge,        & !< dislocation densities for mobility laws 1 and 2
          L_average_screw,            & !< average "effective" length of screw segments
          L_average_edge,             & !< average length of non-screw segments
          tau_average_screw,          & !< average shear stress on screw segments
          tau_average_edge,           & !< average shear stress on non-screw segments
          frankloop*InvVol,           & !< Density of loop
          frankline*InvVol,           & !< Density of infinit lines
          frankfreeline*InvVol,       & !< Density of free infinit lines (not pinned)
          LoopNumb,                   & !< Number of loops
          FreeLoopNumb,               & !< Number of infinit free lines
          ntgd1,                      & !< Number of GD segments of type 1 = made by colli annihilation
          ntgd2,                      & !< Number of GD segments of type 2 = made by cross-slip
          oldntgd1,                   & !< Number of GD segments of type 1 = made by colli annihilation and eliminated
          oldntgd2                      !< Number of GD segments of type 2 = made by cross-slip and eliminated

     close (99)

    open(25,FILE=file_gammap,STATUS='OLD',POSITION='APPEND')
    open(28,FILE=file_rau,STATUS='OLD',POSITION='APPEND')
    open(29,FILE=file_raujonc,STATUS='OLD',POSITION='APPEND')
    open(38,FILE=file_gamma,STATUS='OLD',POSITION='APPEND')
    open(48,FILE=file_travapp,STATUS='OLD',POSITION='APPEND')
    open(49,FILE=file_travint,STATUS='OLD',POSITION='APPEND')

    if ((allocation_dynamique_boites .eqv. .false.) .and. calculate_gammabox) then
      open(61,FILE=file_gammabox,FORM='UNFORMATTED',STATUS='OLD',POSITION='APPEND')

      ! The grid dimension is redefined and is used as header for the gammabox file
      if (KK == KSTAT) then
        write(61)  int(NBoites,DPI_S)
        write(61)  int(NBoitesX,DPI_S),tailleboite(1)
        write(61)  int(NBoitesY,DPI_S),tailleboite(2)
        write(61)  int(NBoitesZ,DPI_S),tailleboite(3)
        write(61)  int(NTSG,DPI_S)
      endif
    endif

    !* multicristals outfiles
    if (DesorientGrain) then
      open(47,FILE='../out/polydens.txt',STATUS='OLD',POSITION='APPEND')
      open(46,FILE='../out/polyrotation.txt',STATUS='OLD',POSITION='APPEND')
      open(45,FILE='../out/polygamma.txt',STATUS='OLD',POSITION='APPEND')
    endif

    ! write(*,'("Vv", F10.5, " m/s SRc ", F10.5 "106 /s,  L = ",F10.5," nm")') &
    !     eps_pl_inst/deltat/schmid/vecburgers/raudis,55*vecburgers*raudis*schmid/1000000,&
    !     raudis*volume*1000000000*avalue**3

    ! 13 colonnes
37  format(1X,13(E11.4,1X))

    toto(1:nsmax) = zero  ! This initialization is needed to set at zero columns between ns+ 1 and nsmax
    itoto(1:nsmax) = izero

    toto(1:ns) = abs(GAMMADOTSYS(1:ns))
    write(25,37) Epso*100, (toto(JK),JK=1,nsmax)

    toto(1:ns) = RAUSYS(1:ns)
    write(28,'(1X,22(E11.4,1X))') Epso*100, RAU_GROSEG, (toto(JK),JK=1,nsmax)

    ! 13 colonnes
337 format(1X,13(E11.4,1X),12I5)

    toto(1:ns) = RAUSYSjonc(1:ns)
    itoto(1:ns) = njoncsys(1:ns)
    write(29,337) Epso*100, (toto(JK),JK=1,nsmax), (itoto(JK),JK=1,nsmax)

    ! 13 colonnes
    toto(1:ns) = GAMSYS(1:ns)
    write(38,37) Epso*100, (toto(JK),JK=1,nsmax)

    ! 13 colonnes
    toto(1:ns) = TrAppSys(1:ns)
    write(48,37) Epso*100, (toto(JK),JK=1,nsmax)

    ! 13 colonnes
    toto(1:ns) = TrIntSys(1:ns)
    write(49,37) Epso*100, (toto(JK),JK=1,nsmax)

    if ((allocation_dynamique_boites .eqv. .false.) .and. calculate_gammabox) then

      ! Ecriture des gamma par systeme dans chaque boite, 13 colonnes par boites
      ! Reecriture du header a chaque pas de sauvegarde
      if (MOD(KK,KSTAT) .eq. zero) then
        ! The simulation time step
        write(61) int(kk,DPI)

        ! In order to check the gammabox calculation, the total gamma is calculated
        somme(:)=0
        do jsum = 1,NTSG
          do isum = 1,NBoites
            somme(jsum) = somme(jsum) + Gammabox(isum,jsum)
          enddo
        enddo
        write(61) (somme(JK) * volboite / volume, JK=1,NTSG)

        ! The gammabox field is saved
        do ni = 1,NBoites
          Px = (real(B1D_3D(ni,1),DP) - half) * tailleBoite(1)
          Py = (real(B1D_3D(ni,2),DP) - half) * tailleBoite(2)
          Pz = (real(B1D_3D(ni,3),DP) - half) * tailleBoite(3)
          write(61) Px,Py,Pz,(Gammabox(ni,JK), JK=1,NTSG)
        enddo

        ! This mark is used to end the step info
        write(61) '&'
      endif

    endif

    !***** donnees specifique au polycristal
    !**************************************
36  format(1X,17(E11.4,1X))
    if (DesorientGrain) then
      !** 17 colonnes
      do ii=1,nbgrains
          write(47,36) Epso*100,raugrainsys(ii,:)
      enddo
      write(47,*)

      !**
      do ii=1,nbgrains
          do ll=1,3
              write(46,*) TensRotGrain(ii,ll,:)
          enddo
          write(46,*)
      enddo
      write(46,*)

      !**
      do ii=1,nbgrains
          write(45,36) sigma*xmu*1.d-6,gamsysgrain(ii,:)
      enddo
      write(45,*)

    endif

    close(25)
    close(28)
    close(29)
    close(38)
    close(48)
    close(49)
    close(47)
    close(46)
    close(45)
    close(44)
    if ((allocation_dynamique_boites .eqv. .false.) .and. calculate_gammabox) close(61)

    if (nsegm+1000 > nsegmax) then
      write (*,*) "NSEGMAX is too small : ",nsegm,'/',nsegmax
    endif

#ifdef PA
    ! Le proc zero a fini d ecrire, on libere tout le monde
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
#endif

endif
!!****** FIN STAT OUT ********

end subroutine STATISMDC
#endif


end module RESUL

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
