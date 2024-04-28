
!===================================================================================================
!========================    DEBUT    MODULE  "ELASTI"   ===========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to calculations directly related to the dislocation
!> elastic theory such as the P-K Forces.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module ELASTI

use VARGLOB
use BRICAMAT
use DEBUG
use CONNEC
use microstructure
#ifdef MDC
use mdctools
#endif

implicit none

real(kind=DP)                             :: SIG_CP(3,3,NSEGMax)  !<
real(kind=DP)                             :: TauPar(Nsegmax)      !<
real(kind=DP)                             :: mufoisadivb          !<
real(kind=DP)                             :: unkt_force           !<

real(kind=DP),dimension(3,NSEGMAX)        :: T      ,&      !< Line vector
                                             B      ,&      !< Burger vector
                                             BPVT           !< Vector normal to B and T vectors
real(kind=DP),dimension(3,3)              :: SIGTMP

integer (kind=DPI),dimension(3,NSEGMAX)   :: RO,         &  !< Coordinates where the stress is calculated
                                             A_COOR,     &  !< Coordinates at segments origin
                                             P1,         &  !< first point used to calculate the radius of curvature
                                             P2,         &  !< Second point used to calculate the radius of curvature
                                             P3             !< Third point used to calculate the radius of curvature
integer (kind=DPI),dimension(NSEGMAX)     :: AB_LENGTH,  &  !< Segment length
                                             VLseg,      &  !< Segment character
                                             tabGD,      &  !< Segment GD type
                                             SegDom         !< Segment domain

Logical,dimension(NSEGMAX)                :: tabjonc, & !<
                                             approxLT

#ifdef MDC
real (kind=DP),parameter                  :: halfsqrt2 = 0.707106781185d0 !< \f$\sqrt{2}\f$

real(kind=DP)                             :: SIG_EF(3,3,NSEGMax)          !<
real(kind=DP)                             :: SIG_NS(3,3,NSEGMax)          !<
real(kind=DP)                             :: unsurxmu                     !<

integer(kind=DPI),dimension(NSEGMAX)      :: listsegmdcstress   !<

integer (kind=DPI)                        :: nsegmdcstress  !<
#endif

#if defined (MDC) && defined(PA)
 real(kind=DP) :: SIGFETEMP(6*NSEGMAX)  !<
#endif


#ifdef PA
! For the multi-phase calculations
real(kind=DP)                             :: SIGMA_JOINT(3,3,NSEGMAX) !<

integer (kind=DPI),allocatable            :: ID_RO(:)   !<
#endif

!========================================================================

contains

!##################################################################
!>  Initialization of loading mode : mainly NstatControl ...      #
!##################################################################
subroutine Init_Solli

implicit none

integer(kind=DPI)  :: i, j, sys
integer(kind=DPI)  :: IVLt        !< Segment line direction index
integer(kind=DPI)  :: IVLd        !< Segment line direction index in the cross slip plan
integer(kind=DPI)  :: SegSlipSys  !< Segment Slip System index

real(kind=DP)      :: y1, y2
real(kind=DP)      :: raudis_tot        !< Total dislocation density
real(kind=DP)      :: YTemp(ntsg)       !< A temp variable appearing in the schmid factor calculation
real(kind=DP)      :: YTempGrain(nbgrains,ntsg)  !< A temp variable appearing in the schmid factor calculation in multi-grain simulation

! Initialization
Ytemp(1:NTSG)                 = -un
YTempGrain(1:NbGrains,1:ntsg) = -un
SIGMA0                        = SIGMA0*1.D6/XMU
sigmapoint                    = sigmapoint*1.D6  ! conversion en Pa
sigmapoint                    = sigmapoint / Xmu
Mufoisadivb                   = XMU * avalue / vecburgers
bdivpa                        = VecBurgers /PII/avalue

write(*,*)
write(*,*) "Loading parameters"
write(*,'("  The shear modulus  = ",F6.1," GPa")') xmu*1D-9
write(*,'("  The Poisson ratio  = ",F6.4," GPa")') DPOISS

if(uniaxiale) then
  write(*,'("  The tensil axis    =  [",3F7.3," ]")') Z (1:3)
  call  ORIENTDIEG(un,Z,SIGAPP,.TRUE.)
  write(*,*)
  write(*,'(" The applied stress tensor")')
  write(*,'(" [",3F7.3," ]")') SIGAPP (:,1)
  write(*,'(" [",3F7.3," ]")') SIGAPP (:,2)
  write(*,'(" [",3F7.3," ]")') SIGAPP (:,3)
else
  write(*,'(" The strain gauge axis    =  [",3F7.3," ]")') Z (1:3)
  write(*,*)
  write(*,'(" The applied stress tensor")')
  write(*,'(" [",3F7.3," ]")') Tensapp (:,1)
  write(*,'(" [",3F7.3," ]")') tensapp (:,2)
  write(*,'(" [",3F7.3," ]")') tensapp (:,3)
endif

! Calculation of the Schmid factor on each segments existin in the initial configuration
!  and the maximum value to calculate the loading CRSS
schmid = zero
if(nsegm == izero) stop " in subroutine Ini_solli, NSEGM = zero ??"

! Identification of the nominal maxi Schmid factor
do i = 1, nsegm

  IVLt        = seg(I)%veclin
  SegSlipSys  = syseg(ivlt)

  ! The multi-grains simulation additional calculations
  if (DesorientGrain) then

    y1 = dabs(SchmidSysGrain(seg(I)%grain,SegSlipSys)) * solli_sys(SegSlipSys)
    YTempGrain(seg(I)%grain,SegSlipSys) = y1

    ! The Schmid factors on cross slip systems of I
    if (GLDEV) then
      do sys = 1, nbsysdev
        IVLd = segdev(IVLt, sys)
        y1 = dabs(SchmidSysGrain(seg(I)%grain,syseg(ivld))) * solli_sys(syseg(ivld))
        YTempGrain(seg(I)%grain,syseg(IVLd)) = y1
      enddo
    endif

  endif

  y1 = dabs(schmidsys(SegSlipSys)) * solli_sys(SegSlipSys)
  ytemp(SegSlipSys) = y1
  if (y1 > Schmid) schmid = y1

  ! The Schmid factors on cross slip systems of I
  if (GLDEV) then
    do sys = 1, nbsysdev
      IVLd = segdev(IVLt, sys)
      y1 = dabs(schmidsys(syseg(ivld))) * solli_sys(syseg(ivld))
      ytemp(syseg(IVLd)) = y1
      if (y1 > Schmid) schmid = y1
    enddo
  endif

enddo

if (schmid == zero) then
  print *, " "
  print *, "Any slip system present in the initial configuration have a Schmid factor /= 0"
  print *, "The code variable -Schmid- is artificially set = 1 !"
  print *, " "
  schmid = un
endif

print *, " "
print *, "Existing slip systems, their Schmid factor and loading status"

! Screen display of the Schmid factors for the activated slip systems
if (DesorientGrain) then

  ! The multi-grains simulation case
  do j = 1, Nbgrains
    write(*, '(" ")')
    write(*, '("   Grain number : ", I2)') j
    write(*, '(" ")')
    do i = 1, NTSG

        if(YTempGrain(j,i) >= -0.01) then
          write(*, '("   System :", I4,"  Schmid =", F6.3,"  Solli =",F6.3)')&
                i,SchmidSysGrain(j,i),solli_sys(i)
        else
          YTempGrain(j,i) = zero
        endif

    enddo
  enddo

else

  ! Only one grain in the volume
  do i = 1, NTSG

      if(ytemp(i) >= -0.01) then
        write(*, '("   System :", I4,"  Schmid =", F6.3,"  Solli =",F6.3)')&
              i,SchmidSys(i),solli_sys(i)
      else
        write(*, '("   Systeme :", I4,"  Schmid =", F6.3,"  Solli =",F6.3, "  Not present")')  &
             i,SchmidSys(i),solli_sys(i)
        ytemp(i) = zero
      endif

  enddo

endif

write(*,*)
write(*, '(" The maximum Schmid factor = ",F6.3,/)') schmid

if (schmid < 0.01 .and. cartograph == 0 .and. .not. key_nucleation  &
    .and. .not. key_crack .and. Mode_deformation /= 8 )    &
  stop "All the Schmid factor equal zero!!!"

! if the Key Shear is true, the value of Sigma0 should be increased in order
! to increase the shear stress = sigma0*Schmid to the initial value sigma0
if (shear) Sigma0 = sigma0 / Schmid
write(*,'(" The initial stress       = ",F6.1," MPa")') sigma0*1D-6*xmu
write(*,'(" The initial shear stress = ",F6.1," MPa",/)') sigma0*1D-6*xmu* Schmid

! Calculation of the total dislocation denity
RauDis_tot = zero
do I = 1 , nsegm
  Raudis_tot = Raudis_tot + SEG(i)%NORME * NormLin(SEG(i)%VECLIN)
enddo
! computing dislocation density
Raudis_tot = raudis_tot / (volume * avalue * avalue)
write(*, '(" Total dislocation density  : ", F10.5," E+12 m-2")') raudis_tot * 1.D-12
write(*, '(" Source dislocation density : ", F10.5," E+12 m-2")') raudis * 1.D-12
write(*, '(" OTHER  dislocation density : ", F10.5," E+12 m-2")') (raudis_tot - raudis) * 1.D-12
write(*, '("===  Average dislocation velocilty (one populated system in shear mode): ",F10.4,1x," m/s")') &
     EpsilonPoint/Raudis/VecBurgers
raudis = raudis_tot

! Important initialization
! Beta = Facteur d'ajustement
! XALPHA = definition du ratio entre taug et taud
! BETA = parametre d ajustement de l'experience pour la proba de CS
! tauIII = contrainte de debut du stade III (restauration)
TAUIII = TAUIII/XMU
Ytemp(1:ntsg) = 350*BDIVA**3*AVALUE**3    !Volume d activation du Cu: V=350*b**3
UNKT_force = XMU*Ytemp(1)*1.0D0/(1.38D-23*TEMPERATURE)
ARGEXP = 1.0D15*BETA*DELTAT*avalue
Z(1:3) = Z(1:3)/DSQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))

! The number of forest segments (immobile segments)
! These segments must be at the begining of the input segment file
IFORET = 0

! Initialization of strain if we are not restarting
if(sideja /= IUN) then
  EPSOLD(:)         = 0.D0
  GAMMAOLDSYS(:,:)  = 0.D0
  Eps_Pl            = 0.D0
endif

! Warning if NSTATCONTROL is not properly defined
if (Nstatcontrol < 10) then
  print*," "
  print *, " Control is made using:", Nstatcontrol, " iterations "
  print *, " It is too small, you should select > 10 !!!"
  stop
elseif (Nstatcontrol > 2000 )then
  print *, " Control is made using:", Nstatcontrol, " iterations "
  print *, " It is too big, you should select < 2000 !!!"
  stop
else
  print *, "Control is made using:", Nstatcontrol, " iterations "
endif

select case(mode_deformation)

  case (0:7)
    ! interesting information to print
    y1 = RAID * (EpsilonPoint) * DELTAT * xmu * schmid * 1.D-6
    y2 = SigmaPoint * DELTAT * xmu * schmid * 1.D-6

    if (shear) then
      write (*,*) "The control is performed on the SHEAR stress  "
      write (*,*) "   "
    endif

    write  (*,'(" The maximum shear stress increment per step = ", F20.15," MPa")') y1
    write  (*,'(" (mode_deformation=3,5,6,7) The maximum shear stress increment per step = ", F20.15," MPa")') y2
    if (y1 > numtold) write  (*,'(" Number of steps for an increase of 1 MPa : = ", F20.2)') 1/y1
    if (y2 > numtold) write  (*,'(" (mode_deformation=3,5,6,7) Number of steps for an increase of 1 MPa : = ", F20.2)') 1/y2

end select

end subroutine Init_Solli

!##############################################################
!> In this subroutine we calculate the simulation loading conditions at each steps
!##############################################################
subroutine SOLLI

implicit none

real (kind=DP)  :: epsilonpointtemp
real (kind=DP)  :: delsig1
real (kind=DP)  :: delsig2
! real (kind=DP)  :: Torque(3,3)
integer         :: i
integer         :: j

!**************************************
!> Calculation of the shift of the tensile axis Z due to plastic rotation
!**************************************
select case (RigidBody)

case(1)
   ZPrime=Z

case(2)
   do i=1,3
      ZPrime(i)=zero
      do j=1,3
        !  ZPrime(i)=Zprime(i)+(MatIdent(i,j)+TensRot(i,j))*Z(j)        ! Plastic rotation is accounted for
         ZPrime(i)=Zprime(i)+(MatIdent(i,j)+TensDist(i,j))*Z(j)       ! Plastic distortion is accounted for
      enddo
   enddo

end select

! During very large simulation a smooth numerical deviation of TensRot symmetry can be
! observed, hence better to normalize ZPrime at each step.
ZPrime=normavec(ZPrime)

!*********************************************************************
!> Calculation of the loading stress tensor as function of the loading axis
!! and the crystal orientation
!*********************************************************************

! The applied stress is imposed only at steps kk > Relax_Reac
if (kk <= Relax_Reac+1) then
  sigma       = zero
  sigapp(:,:) = zero
  if (kk == Relax_Reac+1) then
     sigma = sigma0
     if (uniaxiale) then
        select case (RigidBody) ! no homogen body rotation for desoriented multicrystals
        case(1)
           call  ORIENTDIEG(SIGMA,Z,SIGAPP,.TRUE.)      !Without the solid body rotation
        case(2)
           call  ORIENTDIEG(SIGMA,ZPrime,SIGAPP,.TRUE.)  !With the solid body rotation
        end select
     else
        call  TensIncr(SIGMA,tensapp,SIGAPP,.TRUE.)
     endif
  endif
endif

!> \note The different loading control options \{

!> - if the fatigue loading condition
if(FATIG) then

  if (kk <= Nstatcontrol) then
    DELSIG = zero
  else
    ! The criteria to reverse the loading
    if(REVERS .and. DABS(EPSO) > EPSMAX) then                                 ! Without elastic strain
!     if(REVERS .and. DABS(EPSO + SIGMA/(2.*(1.+Dpoiss))) > EPSMAX) then      ! with the elastic strain
      FSIG = -FSIG
    endif

    ! Definition of the domains where the loading sign can be reversed
    if(DABS(EPSO) > EPSMAX .or. (FSIG*EPSO) < ZERO) then                ! Without elastic strain
!     if(DABS(EPSO + SIGMA/(2.*(1.+Dpoiss))) > EPSMAX .or.     &          ! with the elastic strain
!      (FSIG*(EPSO + SIGMA/(2.*(1.+Dpoiss)) )) < ZERO) then
      REVERS = .false.
    else
      REVERS = .true.
    endif

    ! The default loading mode in fatigue is metalofute
    if(FSIG > ZERO) then
      DELSIG = RAID * (EpsilonPoint - EPSDOT) * DELTAT
      if (DELSIG > zero .and. sigma > zero) DELSIG = DELSIG * dix   ! A trick to speed up load
      if (DELSIG < zero .and. sigma < zero) DELSIG = zero           ! A trick to speed up unload
    else
      DELSIG = RAID * (EpsilonPoint*FSIG - EPSDOT) * DELTAT
      if (DELSIG < zero .and. sigma < zero) DELSIG = DELSIG * dix   ! A trick to speed up load
      if (DELSIG > zero .and. sigma > zero) DELSIG = zero           ! A trick to speed up unload
    endif
  endif

else        ! the non fatigue loading conditions

 ! During the relaxation steps delsig = 0
  if(KK <= Relax_reac) then

    delsig = zero

  else

    ! The default solution is constant stress increment
    delsig = sigmapoint * deltat

    !> - The constant strain rate control
    if(STRATE) then

      if (DesorientGrain) Stop 'Strain rate control cannot be used in multi-grains simulation'

      ! 1 - Stress increment is calculate with instantaneous values of eps_dot
      DELSIG1 = RAID * (EpsilonPoint-EPSDOT) * DELTAT
      ! 2 - Stress increment is purely elastic and calculate with eps_pl
      DELSIG2 = RAID * (EpsilonPoint) * DELTAT * dsign(un,(EpsilonPoint*(KK-Relax_reac)*deltat-Eps_pl))
      ! Stress increment must be delayed, if the rate control (2) have the wrong sign

      if (delsig1*delsig2 < 0) then
        ! zero is prefered to delsig2 in that case
        if (delsig1 < 0) then
          delsig = DELSIG1
        else
          delsig = zero
        endif
      else
        delsig = DELSIG1
      endif

      !*   Les NstatControl first steps must be made with delsig = 0
      if (kk < StepControl) delsig = zero

      ! We boost DELSIG if we are much lower than the friction stress
      if (abs(sigma)*schmid*mufoisadivb < Loi(2)%tau0 * 0.8) DELSIG = DELSIG*dix

    endif

    !> - The 'MetalloFute' control loading mode
    if(METFUT) then

      if(FSIG > ZERO) then
        DELSIG = ZERO
        ! The METFUT control is now based on instantaneous value of epspoint as well as on the average value of
        ! epspoint. The two tests are needed to avoid any delay in the stress increment.
        if (((Sigma*FSIG) <= zero) .or.                                       &
            (epsdot.LT.epsilonpoint .and.                                     &
             (epsold(nstatcontrol)-epsold(nstatcontrol-1)).LT.epsilonpoint)) then
          DELSIG = SigmaPoint *DELTAT  !* &
        endif
      else
        DELSIG = ZERO
        if (((Sigma*FSIG) <= zero) .or.                                       &
            (epsdot > (epsilonpoint*FSIG) .and.                               &
             (epsold(nstatcontrol)-epsold(nstatcontrol-1)) > (epsilonpoint*FSIG))) then
          DELSIG = SigmaPoint*DELTAT*FSIG
        endif
      endif

      if (kk <= Nstatcontrol) DELSIG = zero

    endif

    !> - The 'MetalloFute2' control loading mode (for the materials with high lattice friction)
    if(METFUT2) then

      DELSIG = ZERO
      epsilonpointtemp = densitemob * VecBurgers * VITMAX
      !         print*,"vitmax=",vitmax
      if (epsilonpointtemp < EPSILONPOINT .and. epsilonpointtemp > 0) then
        print*,"epsilonpointtemp = ",epsilonpointtemp,"epsdot=",epsdot
        DELSIG = RAID*(EpsilonPointtemp-EPSDOT)*DELTAT
      else
        DELSIG = RAID*(EpsilonPoint-EPSDOT)*DELTAT
      endif

    endif

    !> - The 'StrainRateMode8' control mode. The FE_Sig loading is adjusted to copy cst strain rate control
    if(StrainRateMode8) then

      if (DesorientGrain) Stop 'Strain rate control cannot be used with the loading mode 8'

      ! 1 - Stress increment is calculate with instantaneous values of eps_dot
      DELSIG = (EpsilonPoint- EPSDOT) * DELTAT / ref_deform

      !*   Les NstatControl first steps must be made with delsig = 0
      if (kk >= StepControl) then

        Stress_fact = Stress_fact + DelSig

        if (Stress_fact > stress_fmax) Stress_fact = stress_fmax
        if (Stress_fact < stress_fmin) Stress_fact = stress_fmin

      else

        Stress_fact = stress_fmin

      endif

      if (modulo(kk,kpredraw) == 0)   write(*,*) "Stress_fact, DeltaStrainRate", Stress_fact,(EpsilonPoint-EPSDOT)

    endif

    !> - The creep control mode  \ingroup creep
    !! \remarks \{ The condition to activate a climb 'timebreak' in the simulation creep mode \n
    !! 1) the plastic strain rate is smaller than the reference imposed strain rate epsilonpoint \n
    !! 2) A sufficient number of glide steps have been made since the last climb timebreak \n
    !! 3) At least one segment is blocked at the boundary of a particle.\}
    if (creep) then

      ! The conditions were dislocation dynamics with only glide is not sufficient and must be "boosted"
      ! by switching to a creep "timebreak". During each timebreak, the next possible event of particle by-pass with climb
      ! is identified and the corresponding segment is allowed to shear its blocking particle.

      ! Initialization
      DELSIG = ZERO
      timebreak = .false.

      ! To permit a new climb 'timebreak' we must have
      ! - the plastic strain rate is smaller than the reference imposed strain rate epsilonpoint
      ! - A sufficient number of glide steps have been made since the last climb timebreak
      ! - At least one segment is blocked at the boundary of a particle.
      if (epsdot*epspointsign < epsilonpoint*epspointsign .and. kk > StepNextTimeBreak .and. npar_touch /= izero) then

        ! A new climb timebreak is required
        timebreak = .true.
        ! The step ending the delay before a new climb timebreak becomes possible is calculated
        ! each climb timebreak event must be separated by at least "nstatcontrol" steps
        StepNextTimeBreak = KK + nstatcontrol

        print*,'timebreak event at step',kk, '   npar_touch=', npar_touch,'   npar_cut=',npar_cut

      endif

    endif

  endif

endif !> \}

! During the relaxation steps no sigma increment are considered
if (kk <= Relax_Reac) then
   DELSIG = zero
endif

!*** Calcul de l'increment sur le tenseur des contraintes
if (uniaxiale) then

  if(shear) DELSIG = DELSIG / Schmid
  SIGMA = SIGMA + DELSIG

!   call  ORIENTDIEG(SIGMA,Z,SIGAPP,.TRUE.)  !Withou the solid body rotation
  call  ORIENTDIEG(SIGMA,ZPrime,SIGAPP,.TRUE.)  !With the solid body rotation

  ! Calculation of the updated value of SchmidSys accounting for the
  ! tensil axis rotation. The sigapp tensor was updated at the previous kk step
  if (sigma /= zero) then
    call GlidCompPK(sigapp, nbase, ntsg, VecNorLin, VecNorNor, SchmidSys)
    SchmidSys = SchmidSys / sigma
  endif

  ! An artificial elastique Torque field is added to compensate the rigid body rotations parallel to the tensile axis.
  ! This contribution must not be taken into account in the Schmid factor analysis
!     call  PBCCompStrain(TORQUE)
!     SIGAPP(:,:) = SIGAPP(:,:) + TORQUE(:,:)      ! The torque force is added to the tensile loading

else

  if(StrainRateMode8) then



  else

    SIGMA = SIGMA + DELSIG

    ! In this simple solution solid body rotation is not accounted for
    call  TensIncr(DELSIG,tensapp,SIGPLUS,.TRUE.)

    SIGAPP(:,:) = SIGAPP(:,:) + SIGPLUS(:,:)
    !   print *, " delsig =", delsig, "sigplus  =", sigplus (1,1)

    ! Calculation of the updated value of SchmidSys accounting for the
    ! tensil axis rotation. The sigapp tensor was updated at the previous kk step
    if (sigma /= zero) then
      call GlidCompPK(sigapp, nbase, ntsg, VecNorLin, VecNorNor, SchmidSys)
      SchmidSys = SchmidSys / sigma
    endif

  endif

endif

#ifdef MDC
if (zcall) then
open(321, file='../out/SIGAPP4Zeb.dat',STATUS='OLD')
! Factor 2 for non-diagonal components using impose_elset_dof(_reaction) in
! Zebulon
! See "Simulation et homogeneisation de muicrostructure periodique - Justin
! Dirrenberger, Samuel Forrest for more details 2010
write(321,*) SIGAPP(1,1)/1.D6*XMU*VOLUME, &
SIGAPP(1,2)/1.D6*XMU*VOLUME*2., &
SIGAPP(1,3)/1.D6*XMU*VOLUME*2., &
SIGAPP(2,2)/1.D6*XMU*VOLUME   , &
SIGAPP(2,3)/1.D6*XMU*VOLUME*2., &
SIGAPP(3,3)/1.D6*XMU*VOLUME
close(321)
endif
#endif

end subroutine SOLLI

!###########################################################################
!#  Cette procedure permet de calculer la Tau resolue qui s exerce sur   #
!# chaque segment. Cette Tau est le fruit du calcul :                    #
!# -1) des contraintes internes au centre du segment "i"                   #
!#    ces contraintes sont la somme de deux contributions. Une             #
!#    contribution issue des segments j tres eloignes et approximee par la #
!#    methode des boites (cette contribution n'est pas necessairement      #
!#    recalculee a chaque pas de simulation). Une contribution des plus    #
!#    proches voisins qui est calculee de facon explicite.                 #
!# -2) La tension de ligne sur les segments (2 caculs possibles)           #
!# -3) La contrainte de Peierls                                            #
!############################################################### 10/12/98 ##
subroutine FORCE

implicit none

! Pour la Methode des Boites
integer(kind=DPI) :: P1i(3),P2i(3),P3i(3),compteur,jj
integer(kind=DPI) :: VLID,INDICE,K,IB,I,J,IVA  !*** DIFFERENTS
integer(kind=DPI) :: nSysDev,Long, VFP
integer(kind=DPI) :: CARAC,I1O,I1E, NBDIM(3),DECAL(3)

real(kind=DP)  :: angleini(nsegm),tmp,adivb,lnquatre,Linedir(3)
real(kind=DP)  :: G_INT(3),G_APP(3),SIGint(3,3),alpha,W,sigprim2(3,3)
real(kind=DP)  :: FORAPP,FORINT,ACTHERM,PROBA,TIR,taueff,taueffdev
real(kind=DP)  :: RAYON2,ANGLE,TFACT,TSIGN,rlong,tmpTAUMON,LogTerm,TAUCRIT
#ifdef MDC
real(kind=DP)  :: SigTFE(3,3)
#endif
#ifdef PA
real(kind=DP), dimension(NSEGM) :: TAUTLSUM,TAUTLCOM
#endif

real(kind=DP), dimension(NSEGM+100) :: DISIG

! Anisotropic line tension variables
real(kind=DP)     :: AngleDeg,DisdAng1,DisdAng2,Tension1,Tension2,LinedirNorm
integer(kind=DPI) :: AngleInd,VLI

real(kind=DP), dimension(NbSysDev+IUN) :: tmpTAUapp, tmpTAUINT

real(kind=DP), dimension(3) :: VUMONT, VCOURB
real(kind=DP), dimension(3) :: ZGrain               !< The tensile axis direction in the grains for multi-grains simulation
logical                     :: condition,condition_cs,FlagPbLT,testbl1,testbl2,acceptLT

! Variables used in the local line tension regularization
!integer(kind=DPI), dimension(3) :: trAB,Tr
!integer(kind=DPI)               :: ivli,ii,jj
!integer(kind=DPI)               :: listP1(11),listP3(11),ROOP1(3,11),ROOP3(3,11),P1temp,P3temp,compteurP1,compteurP3

real(kind=DP), dimension(3)     :: P3P2,P1P2
real(kind=DP)                   :: arclength
!real(kind=DP)                   :: normP3P2,normP1P2,cosP1P2vP3P2,invcosP1P2vP3P2

! Variables used in mode 8 load (imported constant heterogeous stress field)
integer(kind=DPI)    :: ROx,ROy,ROz

REAL(kind=DP) :: Tau_limite   !< Tau_limite is the maxi stress which must be considered in the simulation.
                              !! If a stress value is larger than Tau-limite it means that some computation was wrong
                              !! Tau_limite is defined as 10 time the maximum internal stress defined in the simulation
                              !! control file and used to defined singulier segment.

!*******************
!* Initialization *
!*******************

seg(1:nsegm)%diseg  = .false.
singulier(1:nsegm)  = .false.
TFACT               = un / (4.D0*Pii*(UN-dpoiss))
PLUSSeg             = izero
Loiseg(1:nsegm)     = izero
condition           = .false.
FlagPbLT            = .false.

!*** Symetrisation du tenseur des contraintes appliquees
SIGAPP(2,1) = SIGAPP(1,2)
SIGAPP(3,1) = SIGAPP(1,3)
SIGAPP(3,2) = SIGAPP(2,3)

!*** table intermediaire utile pour l'optimisation de different calculs
#ifdef MDC
SIG_EF(:,:, 1:NSEGM) = ZERO
listsegmdcstress(1:NSEGM) = IZERO
SIG_NS(:,:,1:NSEGM)=ZERO
if (zcall) then
  seg(1:nsegm)%SIGFE(1) = zero
  seg(1:nsegm)%SIGFE(2) = zero
  seg(1:nsegm)%SIGFE(3) = zero
  seg(1:nsegm)%SIGFE(4) = zero
  seg(1:nsegm)%SIGFE(5) = zero
  seg(1:nsegm)%SIGFE(6) = zero
endif
#endif

#if defined (MDC) && defined(PA)
SIGFETEMP(1:6*nsegmax)=zero
#endif

#ifdef PA
TAUTLSUM(1:NSEGM) = ZERO
TAUTLCOM(1:NSEGM) = ZERO
#endif

SIG_CP(:,:,1:NSEGM) = ZERO
TAUAPP(1:NSEGM)   = ZERO
TAUINT(1:NSEGM)   = ZERO
TAUMON(1:NSEGM)   = 1.D-20  ! The climb stress is not initialized exactly to zero to avoid problems in the vclimb computation
TauTOT(1:NSEGM)   = ZERO
TauTL(1:NSEGM)    = ZERO
tabGD(1:nsegm)    = seg(1:nsegm)%gd
tabjonc(1:nsegm)  = seg(1:nsegm)%jonc
VLseg(1:nsegm)    = seg(1:nsegm)%veclin
segDom(1:nsegm)   = seg(1:nsegm)%dom
VL_CS(1:nsegm)    = Izero
AdivB             = UN / BdivA
lnquatre          = DLOG(quatre)
approxLT(1:NSEGM) = .false.


do K = 1,NSEGM

  I = VLseg(K)

  if(K == seg(K)%voiso .and. K == seg(K)%voise) CYCLE

  if(seg(K)%probadev > zero .and. tyseg(I) /= iun)then
      print *, " erreur de CS debut elasti : non vis, I=",K, " KK =",kk
      print *, K,I, tyseg(I),seg(K)%probadev
      call seginfo(K, " avant erreur          ")
      stop
  endif

  DISIG(K) = real(SiSEG(I),DP)
  T(1:3,K) = VecNorLin(1:3,I)
  B(1:3,K) = VecNorLin(1:3,ASSOC(I,1))
  !*** B produit vectoriel T
  BPVT(1,K) = B(2,K)*T(3,K)-B(3,K)*T(2,K)
  BPVT(2,K) = B(3,K)*T(1,K)-B(1,K)*T(3,K)
  BPVT(3,K) = B(1,K)*T(2,K)-B(2,K)*T(1,K)

enddo

! Calculation of the grain's applied stress tensor in multi-grains simulation with disorientation
if (GB /= IZERO .and. DesorientGrain) then
  do I=1,Nbgrains
    ! print*,'I=',I,'MatRotGrain=',MatRotGrain(I,:,:)

    ! Calculation of the grain tensile axis
    do j = 1,3
      ZGrain(j) = MatRotGrain(i,j,1)*Z(1) +   &
                  MatRotGrain(i,j,2)*Z(2) +   &
                  MatRotGrain(i,j,3)*Z(3)
    enddo

    ! Calculation of the stress tensor in the grain
    call  ORIENTDIEG(SIGMA,ZGrain,SigPrim2,.TRUE.)
    ! print*,'I=',I,'SigPrim2=',SigPrim2(:,:) ; read(*,*)

    SIGAPPoriente(I,:,:) = SIGprim2(:,:)

  enddo
endif

!*** Redefinition en reel double precision
AB_LENGTH(1:NSEGM) = SEG(1:NSEGM)%NORME

! Tau_limite is the maxi stress which must be considered in the simulation.
! If a stress value is larger than Tau-limite it means that some computation was wrong
! Tau_limite is defined as 10 times the maximum internal stress used to defined singulier segment
! (usually the line tension is more or less ten times larger than the internal stress)
! defined in the simulation control file
if (.not. key_crack) then
  Tau_limite = TauINT_limite * 10.
else
  ! In the crack loading problem, tau_limite = TauINT_limite
  ! due to the applied  singular stress field around the crack tip.
  Tau_limite = TauINT_limite
endif

!*****************************************************************************
!* Domain decomposition for the multipoles algo and for the obstacles search *
!*****************************************************************************
if (kkdebug) write(379,*) "Before  call region"

Call REGION   ! outputs : Nsegboite and IndexBoite%ListeSeg

#ifdef MDC

#if defined(MDC) && defined(PA)
if (Mon_Rang == IZERO) then
#endif

  ! Send stress locations to Zebulon
  if (zcall) call send_stresspoints_location
  ! Long distance contribution is given by the FEs

  !call recv_stresspoints_stress(mdclrstress, nsegm)

#if defined(PA) && defined(MDC)
endif
#endif

! For short-range interactions, the contribution is given by the analytical solution
! Still TODO - Anisotropic solution - parallelization with monseg(i)

! Attention en parallel il faut intialiser a tout les steps sigint
! puisque la totalite du tableau n est plus calcule par chaque proc
! et que les element pris en compt par les proc varient a chaque step

Call SIGMA_INT_CP

#else

!**********************************************
!*   The short distance stress contribution   *
!**********************************************
if (kkdebug) write(379,*) "Before  Call SIGMA_INT_CP"

Call SIGMA_INT_CP

if (kkdebug) write(379,*) "After  Call SIGMA_INT_CP"

!*******************************************************************
!* Calcul des contraintes a longues distances au centre des boites *
!* le calcul n'est fait que si on a realloue les boite (RCMAKE)    *
!* ou si kk = frequense de raffraichissement KRC                   *
!* SIGMA_INT_LP produit le tableau SigBox(ix,iy,iz,i                *
!* contrainte a longue porte) qui est la contrainte                 *
!* au centre de la boite induite par toutes les boites non_voisines *
!*******************************************************************
if (kkdebug) write(379,*) "Before  Call SIGMA_INT_LP"

if (NBoites > IUN  .and. kk > relax_TL .and. (RCMAKE .or. (modulo(KK-1,KRC) == IZERO))) call SIGMA_INT_LP

if (kkdebug) write(379,*) "After  Call SIGMA_INT_LP"

RCMAKE = .FALSE.   ! Par defaut au prochain pas inutile de recalculer les LP
#endif

! desole : les CPL shift ou tourne ne sont pas pris en compte
if(decalageCPL) stop " Methode de boite incompatoible avec shift"

!****************************************************************************
!* Calcul de la Force de Peach-Koeler resolue dans certaines directions:     *
!*      - vis : intersection plan de glissement, plan de glissement devie    *
!*      - coin : plan de glissement, plan de monte                           *
!* Variables de sortie (en unite reduite : Mu.a^2) :                         *
!* TAUAP    : Tau appliquee dans la direction de glissement                  *
!* TAUAPdev : Tau appliquee dans la direction de glissement devie (sys J)    *
!* TAUEF    : Tau effective dans la direction de glissement                  *
!* TAUEFdev : Tau effective dans la direction de glissement devie (sys J)    *
!* TAUINT   : Tau interne dans le direction de glissement                    *
!* TAUINdev : Tau interne dans le direction de glissement devie   (sys J)    *
!* TAUMON   : The climb resolved shear stress for the non-screw segments     *
!****************************************************************************
#ifdef MDC

  if (zcall) call recv_stresspoints_stress

#endif

! Loop on the domains
do IB = 1,Nboites

  if(NsegBoite(IB) == izero) CYCLE

  ! Loop on the segments in the domains
  bI: do Indice = 1,NsegBoite(IB)

    ! Initialization
    VLI = -100         ! VLI Should never be used without be redefined in the loop
    PROBA = ZERO
    NSYSDEV = IZERO
    I = IndexBoite(IB)%ListeSeg(Indice)

    ! This segment is to be eliminated
    if (i == seg(I)%voiso .and. i == seg(i)%voise) cycle bI

    ! This segment is a pivot segment
    if (AB_length(I) == IZERO) cycle bI

    ! This segment is a GD or JONC segment
    if (tabjonc(i) .or. tabgd(i) > izero) cycle bI

    ! This segment is defined as part of the immobile forest
    if (I.le.IFORET) cycle bI

    ! The segment is not touching a pinning point
    if (SEG(I)%VOISO /= NSEGMAX .and. SEG(I)%VOISE /= NSEGMAX) then

        ! Initialization
        tmpTAUapp(1:nbsysdev+iun) = zero
        tmpTAUINT(1:nbsysdev+iun) = zero
        tmpTAUMON                 = 1.D-20  ! The climb stress is not initialized exactly to zero
                                            ! to avoid problems in the vclimb computation
        VLI                       = VLseg(I)
        long                      = ab_length(I)
        i1o                       = seg(i)%vnno
        i1e                       = seg(i)%vnne

        ! The status of singular segment is defined. A segment is singular if it is a first
        ! neighbour of a junction, or a GD segment, or it is touching a free surface
        condition = (tabjonc(i1o)                       .or. &
                     tabjonc(i1e)                       .or. &
                     (tabgd(i1o) + tabgd(i1e)) > izero  .or. &
                     seg(I)%surface > izero)
        if(condition) then
          singulier(i) = .true.
          if (kkdebug) then
            if (tabjonc(i1o) .or. tabjonc(i1e))       &
              write(379,*)  i, ": is made singulier segment because it touch a GD "
            if ((tabgd(i1o) + tabgd(i1e)) > izero)    &
              write(379,*)  i, ": is made singulier segment because it touch a Junction"
            if (seg(I)%surface > izero)               &
              write(379,*)  i, ": is made singulier segment because it touch a Surface"
          endif
        endif

#ifdef MDC
        if((seg(i)%wait < ITROIS) .or. seg(i)%unload .or. zcall) then

#else
        ! The following is useless if wait > 2
        if(seg(i)%wait < ITROIS .or. seg(i)%unload) then
#endif
          sigint(:,:)=zero
          ! Be aware that carac must be defined outside the parallel part since it is
          ! used latter in the subroutine in a non parallel part
          carac = TYSEG(VLI)

#ifdef PA
          if (monseg(i))then
#endif
            seg(i)%taudev = zero

#ifdef MDC
            !* the long distance contribution + applied stress is given by the FEs
            if (mode_deformation == IHUIT) then

              ROX =int(real((modulo(Ro(1,I),modur(1))-FEgrid_or(1))/(real(FEgrid_size(1))/real(dimFE_Sig(1),DP)),DP))+IUN
              ROY =int(real((modulo(Ro(2,I),modur(2))-FEgrid_or(2))/(real(FEgrid_size(2))/real(dimFE_Sig(2),DP)),DP))+IUN
              ROZ =int(real((modulo(Ro(3,I),modur(3))-FEgrid_or(3))/(real(FEgrid_size(3))/real(dimFE_Sig(3),DP)),DP))+IUN

              SigappS(1,1) = FE_SigappS(ROX,ROY,ROZ,1)
              SigappS(2,2) = FE_SigappS(ROX,ROY,ROZ,2)
              SigappS(3,3) = FE_SigappS(ROX,ROY,ROZ,3)
              SigappS(2,3) = FE_SigappS(ROX,ROY,ROZ,4)
              SigappS(1,3) = FE_SigappS(ROX,ROY,ROZ,5)
              SigappS(1,2) = FE_SigappS(ROX,ROY,ROZ,6)
              SigappS(3,2) = SigappS(2,3)
              SigappS(3,1) = SigappS(1,3)
              SigappS(2,1) = SigappS(1,2)

              Sigapps(:,:) = Sigapps(:,:) * stress_fact

              !sigint(:,:)  = SIG_CP(1:3,1:3,i)+SIGtmp(:,:)* Solli_sys(syseg(VLI))

              if (zcall) then

                sigtFE(:,:)=SIG_EF(1:3,1:3,i)-SIG_NS(1:3,1:3,I)

#if defined (MDC) && defined (PA)
                SIGFETEMP(6*(i-1)+1)=SigtFE(1,1)
                SIGFETEMP(6*(i-1)+2)=SigtFE(2,2)
                SIGFETEMP(6*(i-1)+3)=SigtFE(3,3)
                SIGFETEMP(6*(i-1)+4)=SigtFE(2,3)
                SIGFETEMP(6*(i-1)+5)=SigtFE(1,3)
                SIGFETEMP(6*(i-1)+6)=SigtFE(1,2)
#else
                seg(i)%SIGFE(1)=SigtFE(1,1)
                seg(i)%SIGFE(2)=SigtFE(2,2)
                seg(i)%SIGFE(3)=SigtFE(3,3)
                seg(i)%SIGFE(4)=SigtFE(2,3)
                seg(i)%SIGFE(5)=SigtFE(1,3)
                seg(i)%SIGFE(6)=SigtFE(1,2)
#endif

                Sigint(:,:) = SIG_CP(1:3,1:3,I)-SIG_NS(1:3,1:3,I)+SIG_EF(1:3,1:3,I)

              else

                Sigint(1,1) = seg(i)%SIGFE(1)+SIG_CP(1,1,I)
                Sigint(2,2) = seg(i)%SIGFE(2)+SIG_CP(2,2,I)
                Sigint(3,3) = seg(i)%SIGFE(3)+SIG_CP(3,3,I)
                Sigint(2,3) = seg(i)%SIGFE(4)+SIG_CP(2,3,I)
                Sigint(1,3) = seg(i)%SIGFE(5)+SIG_CP(1,3,I)
                Sigint(1,2) = seg(i)%SIGFE(6)+SIG_CP(1,2,I)
                Sigint(3,2) = Sigint(2,3)
                Sigint(3,1) = Sigint(1,3)
                Sigint(2,1) = Sigint(1,2)

              endif

            else
            ! The applied stress + the long distance contribution

              if (zcall) then

                SigappS(:,:) = SIG_EF(1:3,1:3,i)-SIG_NS(1:3,1:3,I)

#if defined (MDC) && defined (PA)
                SIGFETEMP(6*(i-1)+1)=Sigapps(1,1)
                SIGFETEMP(6*(i-1)+2)=Sigapps(2,2)
                SIGFETEMP(6*(i-1)+3)=Sigapps(3,3)
                SIGFETEMP(6*(i-1)+4)=Sigapps(2,3)
                SIGFETEMP(6*(i-1)+5)=Sigapps(1,3)
                SIGFETEMP(6*(i-1)+6)=Sigapps(1,2)
#else
                seg(i)%SIGFE(1)=SigappS(1,1)
                seg(i)%SIGFE(2)=SigappS(2,2)
                seg(i)%SIGFE(3)=SigappS(3,3)
                seg(i)%SIGFE(4)=SigappS(2,3)
                seg(i)%SIGFE(5)=SigappS(1,3)
                seg(i)%SIGFE(6)=SigappS(1,2)
#endif

              else

                SigappS(1,1) = seg(i)%SIGFE(1)
                SigappS(2,2) = seg(i)%SIGFE(2)
                SigappS(3,3) = seg(i)%SIGFE(3)
                SigappS(2,3) = seg(i)%SIGFE(4)
                SigappS(1,3) = seg(i)%SIGFE(5)
                SigappS(1,2) = seg(i)%SIGFE(6)
                SigappS(3,2) = SigappS(2,3)
                SigappS(3,1) = SigappS(1,3)
                SigappS(2,1) = SigappS(1,2)

              endif

              ! The short-range contribution inside the inclusion

              sigint(:,:)  = SIG_CP(1:3,1:3,i)

            endif
#else
            if (Nb_phase /= 2) then

              Sigint(:,:) = SIGBOX(1:3,1:3,IB) + SIG_CP(1:3,1:3,I)

            else
#ifdef PA
              ! recall: Nb_phase == 2  is posible only in parallel
              Sigint(:,:) = SIGBOX(1:3,1:3,IB) + SIG_CP(1:3,1:3,I) + SIGMA_JOINT(1:3,1:3,I)
#endif
            endif

            ! If the applied stress is provided by an external calculation, e.g.FEM
            ! We simply define the valeue of SigappS at Ro(i)
            If (mode_deformation == IHUIT) then

              ROX =int(real((modulo(Ro(1,I),modur(1))-FEgrid_or(1))/(real(FEgrid_size(1),DP)/real(dimFE_Sig(1),DP)),DP))+IUN
              ROY =int(real((modulo(Ro(2,I),modur(2))-FEgrid_or(2))/(real(FEgrid_size(2),DP)/real(dimFE_Sig(2),DP)),DP))+IUN
              ROZ =int(real((modulo(Ro(3,I),modur(3))-FEgrid_or(3))/(real(FEgrid_size(3),DP)/real(dimFE_Sig(3),DP)),DP))+IUN

              ! print*, 'ROx,Roy,Roz,ro(:,i)', ROx,Roy,Roz,ro(:,i)
              ! print*, 'FEgrid_or,FEgrid_or+FEgrid_size', FEgrid_or,FEgrid_or+FEgrid_size

              SigappS(1,1) = FE_SigappS(ROX,ROY,ROZ,1)
              SigappS(2,2) = FE_SigappS(ROX,ROY,ROZ,2)
              SigappS(3,3) = FE_SigappS(ROX,ROY,ROZ,3)
              SigappS(2,3) = FE_SigappS(ROX,ROY,ROZ,4)
              SigappS(1,3) = FE_SigappS(ROX,ROY,ROZ,5)
              SigappS(1,2) = FE_SigappS(ROX,ROY,ROZ,6)
              SigappS(3,2) = SigappS(2,3)
              SigappS(3,1) = SigappS(1,3)
              SigappS(2,1) = SigappS(1,2)

              Sigapps(:,:) = Sigapps(:,:) * stress_fact

              !read(*,*)

              SigappS(:,:) = stress_fact * SigappS(:,:)

            endif

            ! Solution for a loading stress gradient along z axis in the simulation box (bending test)
            ! SigappS(:,:) = SigappS(:,:) * real((RO(3,I)-(Modur(3)/IDEUX)),DP)/real((Modur(3)/IDEUX),DP)

            !## Calculate stress tensor for a sigle-edge cracked specimen under tension along z direction
            if ( key_crack .and. ( mode_deformation <= IQUATRE .or. mode_deformation >= ITROIS ) ) then

              call calculate_crackfield(RO(1:3,I))

              SigappS(1,1) = Sigcrack(1) !+ Sigapp(1,1)
              SigappS(2,2) = Sigcrack(2) !+ Sigapp(2,2)
              SigappS(3,3) = Sigcrack(3) !+ Sigapp(3,3)

              SigappS(1,2) = Sigcrack(4) !+ Sigapp(1,2)
              SigappS(2,3) = Sigcrack(5) !+ Sigapp(2,3)
              SigappS(1,3) = Sigcrack(6) !+ Sigapp(1,3)

              SigappS(3,2) = SigappS(2,3)
              SigappS(3,1) = SigappS(1,3)
              SigappS(2,1) = SigappS(1,2)

            endif

            ! sigtmp est le tenseur applique avec poids de systems impose
            ! Attention : ICI on intervient pour changer artificiellement la composante des
            ! contraintes appliquees et , de fait, changer tautot et tauef
            ! toute l'informatin est contenu dans Solli_sys(((IVL(I)-1)/nbasered)+1)
            if (mode_deformation /= IHUIT .and. .not. key_crack) then
              if (DesorientGrain) then
                SigappS(:,:) = SIGAPPoriente(seg(I)%grain,:,:)*Solli_sys(syseg(VLI))
              else
                SigappS(:,:) = SIGapp(:,:)* Solli_sys(syseg(VLI))
              endif
            endif
#endif

            if (seg(i)%wait < ITROIS .or. seg(i)%unload) then

            !*** Calcul de G = SIGEFF.B
            !*   G_INT =  vecteur interne
            !*   G_app =  vecteur appliquee
            G_int(:) = zero
            G_app(:) = zero

            do J = 1,3
              G_int(1) = G_int(1)  + SIGint(1,J)*B(J,I)
              G_int(2) = G_int(2)  + SIGint(2,J)*B(J,I)
              G_int(3) = G_int(3)  + SIGint(3,J)*B(J,I)

              G_app(1) = G_app(1) + SIGappS(1,J)*B(J,I)
              G_app(2) = G_app(2) + SIGappS(2,J)*B(J,I)
              G_app(3) = G_app(3) + SIGappS(3,J)*B(J,I)
              !print *,'J=',j,'SigAPPs(1,J)',SIGappS(1,J),'B(J,I)',B(J,I)
            enddo
            !Print *,'G_APP=',G_app

            !*   FORCE = G.produit_vectoriel.T
            !*   TAU = FORCE.Direction de glissement normalisee mais non signee
            !*   La direction 1 est la direction de glissement
            !*   Les directions > 1 sont direction de glissement deviee
            if(carac == 1) then
              !**  Cas des vis,
              do IVA = 1,3
                J = IPERM(IVA)  !*** SUR LES INDICES JE SUPPOSE
                K = IPERM(IVA+1)

                FORAPP = G_APP(J) * T(K,I) - G_APP(K) * T(J,I)
                FORINT = G_int(J) * T(K,I) - G_INT(K) * T(J,I)

                tmp =  VecNorDep(IVA,VLI) ! vecteur dep du seg de vecteur ligne VLI

                tmpTAUapp(iun) = tmpTAUapp(iun) + FORAPP * tmp
                tmpTAUINT(iun) = tmpTAUINT(iun) + FORINT * tmp

                do compteur = IDEUX , NbSysDev + iun

                    ! print *, " ===============",VecNorDep(IVA,IVLD(I,compteur - iUN))
                    ! print *, VLI,compteur, IVLD(I,compteur - iUN)

                    ! vecteur dep dans le sys system devie numero J
                    VLID = Segdev(VLI,compteur-IUN)
                    tmp = VecNorDep(IVA,VLID)

                    tmpTAUapp(compteur) = tmpTAUapp(compteur) + FORAPP * tmp
                    tmpTAUINT(compteur) = tmpTAUINT(compteur) + FORINT * tmp

                enddo

              enddo

            else

              !**  Cas des coins et mixtes,
              do IVA=1,3

                J = IPERM(IVA)
                K = IPERM(IVA+1)
                VUMONT(IVA) = (VecNorDep(J,VLI)*T(K,I) - VecNorDep(K,VLI)*T(J,I)) ! Vecteur unitaire de montee
                FORAPP = G_APP(J) * T(K,I) - G_APP(K) * T(J,I)
                FORINT = G_int(J) * T(K,I) - G_int(K) * T(J,I)

                tmp =  VecNorDep(IVA,VLI) ! vecteur dep dans le sys system devie numero J

                tmpTAUapp(iun) = tmpTAUapp(iun) + FORAPP * tmp
                tmpTAUINT(iun) = tmpTAUINT(iun) + FORINT * tmp
                tmpTAUMON      = tmpTAUMON + (ForApp + forint) * VUMONT(IVA)

              enddo

            endif

! #ifdef PA
!                 endif
!                 ! The computations bypassed in parallel mode must be stop
!                 ! here since important quantity like anglevis are determined
!                 ! in the line tension calculation. This aspect could be
!                 ! modified to optimized the simulation.
! #endif

            taueff = tmpTAUapp(iUN) + tmpTAUint(iUN)
            !print *,'taueff',tmptauapp,seg(i)%veclin
            !read(*,*)
            taueffdev = izero
            SEG(I)%nsysdev  = -IUN

            if(carac == 1) then
              ! determination of the cross system with maximum RSS
              do iva = ideux , nbsysdev + iun
                tmp = tmpTAUapp(iVA) + tmpTAUint(iva)
                if(abs(tmp) >= abs(taueffdev)) then
                  taueffdev = tmp
                  nsysdev = iva - IUN
                endif
              enddo
              seg(i)%Nsysdev = nsysdev

              if(CFC) then
                ! Here, to facilitate anihilation between close
                ! segments, only TAUint is considered at the
                ! deviation step !!!
                seg(i)%taudev = tmpTAUint(nsysdev+IUN)
              else
                seg(i)%taudev = taueffdev
              endif

              !-----------------------------------------------------------------!
              !           write (*,'(I10, E20.5,10x,E20.5)') I, taueff,  taueffdev
            endif

          !********************************************************************************
          !** line tension ************************************************************
          !********************************************************************************

          !-------------------------------------------------------------------!
          !     CETTE PARTIE CALCUL LA TENSION DE LIGNE A PARTIR DU RAYON DE  !
          !     COURBURE D UN CERCLE PASSANT PAR LE CENTRE DE TROIS SEGMENTS  !
          !     MITOYENS. (UNITES : FORCE EN MU.A^2)       !
          !-------------------------------------------------------------------!

          ! anglevis is initialized with a stupid vale to check if the calculation is correctly made in the following
          seg(i)%anglevis = -1.0

          ! This quantities are defined exactly in the sigma_int_cp subroutine
          P1i(:) = P1(:,i)
          P2i(:) = P2(:,i)
          P3i(:) = P3(:,i)

          P1P2(1:3) = P1i(1:3) - P2i(1:3)
          P3P2(1:3) = P3i(1:3) - P2i(1:3)

          testbl1=.false.
          testbl2=.false.

          ! The segment is touching a surface or interface
          if (seg(i)%surface > IZERO) then

            if (seg(i)%surface > ITROIS) then
              print *,'Problem with segment touching free surface (line tension) at step',KK,I,seg(i)%surface
              stop
            endif

            ! If segment on free surface and changed from domain, Line tension is not calculate the next step
            if (seg(i)%zerotl) then

              Tautl(I)      = ZERO
              seg(i)%zerotl = .False.
              seg(I)%anglevis = (tyseg(seg(I)%veclin)-1)*dble(Pii)*quart   ! we approximate anglevis

            else

              if (seg(i)%surface > IUN) then

                do jj = nbplanDom + 1,nbplanMax - NbFreePlan !test on P1i: if it belongs to fixed boundary

                  if (P1i(1)*Plane_MillerI(1,jj) + &
                      P1i(2)*Plane_MillerI(2,jj) + &
                      P1i(3)*Plane_MillerI(3,jj) == Plane_pos(jj)) testbl1=.true.
                enddo

                Linedir(1:3)=real((P2i(1:3) - P1i(1:3)),DP)

              else

                do jj = nbplanDom + 1,nbplanMax - NbFreePlan

                  if (P3i(1)*Plane_MillerI(1,jj) + &
                      P3i(2)*Plane_MillerI(2,jj) + &
                      P3i(3)*Plane_MillerI(3,jj) == Plane_pos(jj)) testbl2=.true.
                enddo

                Linedir(1:3)=real((P2i(1:3) - P3i(1:3)),DP)

              endif

              Linedirnorm = (Linedir(1)*Linedir(1) + Linedir(2)*Linedir(2) + Linedir(3)*Linedir(3))

              if (Linedirnorm == 0.) then
                write(379,*) 'Segment ', I
                write(379,*) 'Norm ' , seg(i)%norme
                write(379,*) 'P1i ', P1i
                write(379,*) 'P2i ', P2i
                write(379,*) 'P3i ', P3i
                write(379,*) 'LineDir ', Linedir
                write(379,*) 'Linedirnorm = ZERO Elasti - Surface line tension '

                call seginfo(i, 'ZERO Elasti - Surface line tension')
                stop "Linedirnorm = ZERO Elasti - Surface line tension"
              endif

              if (testbl1 .or. testbl2) then

                Tautl(i)=ZERO

              else

                Linedir(:) = Linedir(:) / dsqrt(Linedirnorm)

                ANGLE=Linedir(1)*vecnorlin(1,assoc(VLi,1))+   &
                      Linedir(2)*vecnorlin(2,assoc(VLi,1))+   &
                      Linedir(3)*vecnorlin(3,assoc(VLi,1))

                if (dabs(angle) > UN) angle = dsign(UN,angle)
                seg(i)%anglevis = dacos(ANGLE)
                ANGLE       = dacos(ANGLE)
                angleini(i) = angle
                TauAPP(I)   = tmpTAUapp(Iun)
                VFP = seg(i)%VARFREEPLAN

                acceptLT = .false.    ! The default solution is not LT calculation

                if (seg(i)%surface==IUN) then

                  acceptLT = (P3P2(1)*Plane_MillerR(1,VFP) +         &
                              P3P2(2)*Plane_MillerR(2,VFP) +         &
                              P3P2(3)*Plane_MillerR(3,VFP) < ZERO)

                else if (seg(i)%surface==IDEUX) then

                  acceptLT = (P1P2(1)*Plane_MillerR(1,VFP)+         &
                              P1P2(2)*Plane_MillerR(2,VFP)+         &
                              P1P2(3)*Plane_MillerR(3,VFP) < ZERO)

                endif

                !Line tension is calculated for segment touching the surface ONLY if segment points toward the surface

                if (acceptLT) then

                  TauTl(i)=surfLT(i,vli,angle,Linedir)

                  if (dabs(TauTL(I)) > Tau_limite) then
                    TauTL(I)  = Tau_limite * Dsign(un,TauTL(I))
                    singulier(i) = .true.
                  endif

                else

                  TauTL(i)=IZERO

                endif


              endif

            endif

          else

            ! si les segments sont alignes,prodvec = zero, TL = zero, sans calcul
            Decal(1:3) = P1i(1:3) - P2i(1:3)
            NBdim(1:3) = P3i(1:3) - P2i(1:3)

            condition = (((Decal(2)*NBdim(3) /= Decal(3)*NBdim(2))  .or.             &
                          (Decal(3)*NBdim(1) /= Decal(1)*NBdim(3))  .or.             &
                          (Decal(1)*NBdim(2) /= Decal(2)*NBdim(1))) .or. approxLt(i))

            ! the following test is true only when P1i,P2i,P3i are aligned
            if (condition) then

              ! Calculation of curvature radius vector
              if (approxLt(i)) then
                VCOURB(1:3) = -(P1P2(1:3)+P3P2(1:3))
                alpha = UN
              else
                VCOURB(1:3) = VECCOURBCRAMER(P1P2,P3P2,alpha)
              endif

              ! The square norm of VCOURB is the curvature radius
              RAYON2 = VCOURB(1)*VCOURB(1) + VCOURB(2)*VCOURB(2) + VCOURB(3)*VCOURB(3)

              ! If P1i, P2i and P3i are aligned, the value of the line tension is equal to zero
              if (RAYON2 /= ZERO) then

                ! Angle de la tangente locale a la direction vis
                ! (calcul utile seulement dans le cas anisotrope)
                ANGLE = VCOURB(1)*vecnorlin(1,assoc(VLi,1)) +   &
                        VCOURB(2)*vecnorlin(2,assoc(VLi,1)) +   &
                        VCOURB(3)*vecnorlin(3,assoc(VLi,1))

                ANGLE = ANGLE / dsqrt(rayon2)

                ! angle between the normal to I and b
                if (dabs(angle) > UN) angle = dsign(UN,angle)

                seg(i)%anglevis = dabs(dble(PII) * HALF - dacos(ANGLE))
                ANGLE       = dble(PII)      * HALF - dacos(ANGLE)

                angleini(i) = angle

                ! TSIGN est la projection de VCOURB  dans la direction de deplacement
                TSIGN = (VECnorDEP(1,VLI)*VCOURB(1) + &
                         VECnorDEP(2,VLI)*VCOURB(2) + &
                         VECnorDEP(3,VLI)*VCOURB(3))

                RAYON2 = dsqrt(rayon2)

                arclength = alpha*RAYON2


                TSIGN  = SIGN(1.D0,TSIGN)

                ! different line tension approximation
                select case(linten)

                !-------------------------------------------------------------------------
                case(Izero)     !Line tension in the constante energy approximation: Friedel Formula

                  TAUTL(I) = TSIGN * TFACT * BdivA / RAYON2

                !-------------------------------------------------------------------------
                case(Iun)       ! formula of Bacon [1967] based on DeWitt

                  LogTerm = DLOG(arclength*AdivB)

                  if (abs(LogTerm) < un) then
                    LogTerm = 1.
                  endif

                  TAUTL(I)= TSIGN * TFACT * BdivA /(RAYON2)*                   &
                          ((UN + dpoiss - trois * dpoiss * DSIN(ANGLE)**2)*    &
                          LogTerm - (deux*dpoiss*COS(2.0D0*ANGLE)))

                !-------------------------------------------------------------------------
                case(Ideux)     !Isotropic line tension: Foreman formula

                  LogTerm = DLOG(arclength*AdivB)

                  if (abs(LogTerm) < un) then
                    FlagPbLT = .true.
                    LogTerm = 1.
                  endif

                  TAUTL(I)= TSIGN * TFACT * BdivA /(RAYON2)*                   &
                          ((UN + dpoiss - trois * dpoiss * DSIN(ANGLE)**2)*    &
                          LogTerm - (dpoiss*COS(2.0D0*ANGLE)))

                !-------------------------------------------------------------------------
                case(Itrois)    !Anisotropic line tension solution; a Tabulated formula

                  AngleDeg = abs(angle*180./PII) ! angle in now defined with degree

                  !Identification of the closest values of the tabulated line tension values
                  if (AngleDeg == 180.) then
                    AngleInd = 90
                  else
                    AngleInd=Int(AngleDeg*DisdInc)+1
                    if (AngleInd < 1 .and. AngleInd > 90) &
                      stop  'There is a probleme with the inisotropic line tension interpolation'
                  endif

                  DisdAng1 = Disd(1,Angleind)
                  DisdAng2 = Disd(1,Angleind+1)
                  Tension1 = Disd(3,Angleind)
                  Tension2 = Disd(3,Angleind+1)

                  TAUTL(I)= TSIGN*(Tension1+(AngleDeg-DisdAng1)*(Tension2-Tension1)/    &
                          (DisdAng2-DisdAng1))/Rayon2

                !-------------------------------------------------------------------------
                case(Iquatre)      !Isotropic line tension: Mohles formula
                    ! new modification:Ghiath: the proportionality with log R/b is only valid
                    ! for R >2b. It not, the enrgy goes to zero almost proportionally to R
                    if(arclength > deux*Bdiva) then
                      LogTerm = DLOG(deux*arclength*AdivB)
                    else
                      LogTerm = (half * arclength *AdivB) * Lnquatre
                    endif
                    TAUTL(I)= TSIGN * TFACT * BdivA /(RAYON2) *               &
                          (UN + dpoiss - Trois* dpoiss * DSIN(ANGLE)**2) *     &
                          LogTerm

                endselect

                if (dabs(TauTL(I)) > Tau_limite) then
                  TauTL(I)  = Tau_limite * Dsign(un,TauTL(I))
                  singulier(i) = .true.
                endif

              else

                ! The line tension cannot be calculated for this segment (very strange connection of segment)
                TAUTL(i)        = ZERO
                seg(i)%anglevis = (tyseg(seg(I)%veclin)-1)*dble(Pii)*quart   ! We approximate anglevis as much as possible

              endif

            else

              ! The line tension equal zero since the line section is flat
              TAUTL(i)        = ZERO
              seg(i)%anglevis = (tyseg(seg(I)%veclin)-1)*dble(Pii)*quart   ! We approximate anglevis as much as possible

            endif

          endif

          ! The special case of segments with anglevis not very well defined
          ! Those segments are usually oscillating and therefore are polluting
          ! the statistic on sweep areas.
          ! all those segment are set of mixte character
          if (SUM(P1P2(:)*P1P2(:)) < rlocal2 .or. SUM(P3P2(:)*P3P2(:)) < rlocal2) seg(i)%anglevis = un


          !**************************************
          !* Tau final = fonction de mode de relaxation : Tau_INT+Tension de ligne+Tau_TL
          !**************************************
          ! stocking values of force components in globale arrays
          !           TAUAPP(I) = tmpTAUapp(IUN)
          !           TAUINT(I) = tmpTAUINT(IUN)
          !           TAUMON(I) = tmpTAUMON

          if (I.GT.IFORET) then
              if(kk > Relax_Reac) then
                ! Correction de FTL et Fint si probleme de procedure de calcul: remise a Tau_limite"
                ! print *,Tau_limite,1D12*rlong(I)/mufoisadivb(syseg(ivl(i)))
                if (dabs(tmpTAUapp(Iun)) > Tau_limite)then
                    tmpTAUapp(Iun) = Tau_limite * Dsign(un,tmpTAUapp(Iun))
                    singulier(i) = .true.
                endif
                TauAPP(I) = tmpTAUapp(Iun)
              endif

              if (dabs(tmpTauINT(Iun)) > TauINT_limite) then
                singulier(i) = .true.
                if (dabs(tmpTauINT(Iun)) > Tau_limite) then
                  tmpTauINT(IUN)  = Tau_limite * Dsign(un,tmpTauINT(Iun))
                endif
              endif
              ! in cases where tau int is too small it should be truncated
              !                print *, I,tmptauINT(IUN)
              if(abs(tmptauINT(IUN)) > 1E-7) then
                TauINT(I)  = tmpTauINT(Iun)
              else
                TauINT(I)  = zero
              endif

              !> The climb stress on segment I is saved in tab 'TauMon' \ingroup creep
              if (creep) TauMon(I) = tmpTauMon
          endif

          ! important: if the segment is of effective screw character or
          ! it is screw and enough long ,or.
          ! it is the neighbor of a screw segment
          if (i1o /= Nsegmax .and. i1e /= Nsegmax) then

            if (seg(i)%anglevis < angle_vis                                                     .or.    &
                (carac == 1 .and. (tyseg(seg(i1o)%veclin)==1 .or. tyseg(seg(i1e)%veclin)==1))       ) then

              LoiSeg(i) = int(numero_loi(assoc(VLi,1)))

            elseif (carac == 1) then

              ! important: if the segment is not effective screw it can be the prolongation
              ! of an effective screw segment, we attribute to it the effective screw character
!                 if(tyseg(seg(i1o)%veclin) == 1
              ! important: if the segment is not effective screw but
              ! it is screw, we attribute to it the mobility law defined for the 2d character
              LoiSeg(i) = int(numero_loi(assoc(VLi,2)))

            else

              LoiSeg(i) = int(numero_loi(VLi))

            endif

          else

            ! This segment is touching a pinning point
            LoiSeg(i) = int(numero_loi(VLi))

          endif

          !**********************
          !** Glissement devie **
          !**********************

          rlong           = normlin(VLI) * long
          SEG(I)%PROBADEV = ZERO
          proba           = zero

!#ifdef PA
!         ! The proba calculation must be bypassed since the stress
!         ! was calculated only for the monseg(i) list
!         if (monseg(i))then
!#endif
          if (GLDEV .and. KK > Relax_reac) then
            !print*,'debut GLDEV for segment i=',i

            !*1) Only screw segments can cross-slip
            if (carac == 1) then
              !print *,'carac == 1'

              !*2) To cross-slip, the segment length must be long enough
              !*    to facilitate a bowout in the cross-slip plane
              LogTerm = DLOG(RLONG*AdivB)
              if (RLONG < bdivA) LogTerm = 1.d0

              ! The bowout stress (10% of Foreman critical stress)
              TAUCRIT = 0.1 * (BdivA / (DEUX*PII*RLONG) * LogTerm)

              ! If the segment is too small we do not want to test CS whatever the stress
              if (RLONG < XLOMAX*INVDIX) then
                !print*, 'segment i is too small to cross-slip',i,TAUEFFDEV,TAUCRIT
                TauCrit = TauEffDev
              endif

              if(abs(TAUEFFDEV) > abs(TAUCRIT))then
                !print *,'TAUCRIT passed' ,abs(TAUEFFDEV),abs(TAUCRIT),abs(TauEFF),RLONG

                !*3) Segments pinned at one end cannot cross-slip
                if (SEG(I)%Vnno /= NSEGMAX  .AND. SEG(I)%Vnne /= NSEGMAX) then
                  !print *,'pinned segment passed'

                  !*4) The segment character must be close to the screw direction
                  if(seg(i)%anglevis <= angle_vis) then
                    !print *,'screw direction passed'

                    !*5) Segments blocked at obstacle interface cannot cross-slip
                    if((.not.seg(i)%bloquer).and.(.not.seg(seg(i)%voiso)%bloquer)                 &
                        .and.(.not.seg(seg(i)%vnno)%bloquer).and.(.not.seg(seg(i)%voise)%bloquer) &
                        .and.(.not.seg(seg(i)%vnne)%bloquer))  then
                      !print *,'bloquer passed'

                      !*6) Additional rules for FCC crystals
                      if(CFC) then

                        !*7) Balance between stress in glide and CS plane
                        if(XALPHA*DABS(TAUEFFdev) >= DABS(TAUEFF)) then
                        !print *, 'passed the stress balance',I, TAUEFFdev,TAUEFF
                          !** The cross-slip activation energy computation

                          ! A simple criteria to identify quasi immobile segment
                          if (DABS(Taueff) <= 0.2 * DABS(tmpTauAPP(iUN))) then

                            ! expression for quasi immobile segments
                            ! The old Kubin-Devincre formula
                            ! W = ((DABS(tmpTAUint(iUN))) - TAUIII)*UNKT_force
                            ! The formulla modified by Sansal-Devincre to promote cs on segments
                            ! scillating around an equilibrium position within a strong stress gradient
                            W = (((DABS(tmpTAUapp(iUN))+DABS(tmpTAUint(iUN)))*half)-TAUIII)*UNKT_force
                            !W = (((DABS(tmpTAUapp(iUN)))*half)-TAUIII)*UNKT_force

                          else

                            ! Expression for mobile segments, strongly attracted in the cs plane
                            ! This expression is based on Brown model for screw annihilation
                            W = (DABS(TAUEFFdev) - TAUIII)*UNKT_force

                          endif

                          ACTHERM = 1.4D22      !********* 1.2 au depart

                          ! A cut-off is used to avoid problems with the exponantiel
                          if(W < 51.D0) ACTHERM = DEXP(W)

                          ! Cross-slip criteria for FCC segments
                          PROBA = ARGEXP * ACTHERM * RLONG

                        endif

                      elseif(HCP .or. BCC .or. DC) then
                        ! Cross-slip criteria for HCP and BCC segments

                        if(dabs(TAUEFFdev) >=DABS(TAUEFF))proba = UN
                        !if(dabs(TAUEFFdev) >=DABS(TAUEFF)) then
                        !   print *,'taueff passed',taueffdev,taueff,seg(i)%veclin
                        !   read(*,*)
                        !endif

                      else

                        stop "No rules defined for CS in this crystallographic structure "

                      endif

                      seg(i)%probadev = PROBA   ! The cross-slip probability
                      ! To diffuse the nsysdev information in parallel computation
                      ! when wait=2 the VL_CS tab is used. A negative value is set
                      ! to save the informatioin and say that the segment did not
                      ! cross slip at that step
                      if (seg(i)%wait == IDEUX) VL_CS(i) = -seg(i)%nsysdev

                    endif
                  endif
                endif
              endif
            endif
          endif
          endif
#ifdef PA
        endif !monseg(i)
#endif
      endif
    endif

    ! For the WAIT Algo we must save some information
    if(SEG(I)%WAIT > IDEUX) then
      Proba = SEG(I)%PROBADEV
      Nsysdev = SEG(I)%nsysdev
    endif

#ifdef PA
    if (monseg(i))then
#endif

    ! The Metropolis criteria : We test the proba of cross-slip
    if(PROBA > 1.0D-5) then
        TIR = taus88()
!            TIR = 0.5         ! To test CS in parallel runing it is necessary to set a constant value
        do while (TIR < zero .or. TIR >= un)
            TIR = taus88()
        enddo

        ! The metropolis test
        if(TIR < PROBA) then

          ! the following test is to prevent CS when a close segment at origin or extremity
          ! is of Jonc or GD type
          condition_CS = .true.
          i1o = seg(i)%voiso
          i1e = seg(i)%voise
          iva = Loma(VLseg(I))
          compteur = Izero
          do while ((i1o /= nsegmax) .and. (compteur < IVA) .and. (condition_cs))
              compteur = compteur + AB_length(i1o)
              if (tabgd(i1o) > izero .or. tabjonc(i1o)) condition_cs = .false.
              i1o = seg(i1o)%voiso
          enddo
          compteur = Izero
          do while ((i1e /= nsegmax) .and. (compteur < IVA) .and. (condition_cs))
              compteur = compteur + ab_length(i1e)
              if (tabgd(i1e) > izero .or. tabjonc(i1e)) condition_cs = .false.
              i1e = seg(i1e)%voise
          enddo

          IF(condition_cs) THEN

              !*** TOUT EST LA :
              !*** TOUTES LES CONDITIONS SONT REUNIES : GLISSEMENT DEVIE
              ! At this point we know if a segment is supose to cross-slip
              ! but the modification of the connectivity and glide plane is
              ! made latter in the depredic subroutine
              VL_CS(I) = segdev(VLI,nsysdev)     ! the new veclin index
              singulier(i) = .true.              ! CS segment are always initially singulier
              TauTOT(I) = seg(i)%taudev
              ! print *,i,kk,'cross slip condition true'
          endif
        endif
    endif

#ifdef PA
    endif
#endif

  enddo bI      ! End of the loop on the segments in domains

enddo           ! End of the loop on domains

#if defined (MDC) && defined (PA)
if (zcall)  call sig_FEsincro
#endif

#ifdef PA
TAUTLCOM(1:NSEGM)=TAUTL(1:NSEGM)
CALL MPI_ALLREDUCE(TAUTLCOM,TAUTLSUM,NSEGM,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_PAIRIMPAIR,IERR)
TAUTL(1:NSEGM) = TAUTLSUM(1:NSEGM)
#endif

! A warning is made, to inform that the Log expression in the line tension formula was modified at this step
!if (FlagPbLT) write(*,*) 'limits of line tension definition are reached at step!',KK

b3: do I = 1,NSEGM
#ifdef PA
   if (monseg(i))then
#endif

  if (i == seg(I)%voiso .and. i == seg(i)%voise) cycle

  if(LoiSeg(i) == IZERO) THEN
    LoiSeg(i) = int(numero_loi(VLseg(I)))
  endif

  if ((VL_CS(I) <= izero .and. SEG(i)%wait < itrois) .or. seg(i)%unload) then

    ! modification of the line tension tension value :
    ! when I touchs a pinning point the curvature of I is ill defined. Its LT is replaced
    ! by that of itsd neighbor
    if (seg(i)%surface == IZERO) then
      if(seg(i)%vnno == nsegmax) tauTL(I) = tauTL(seg(i)%vnne)
      if(seg(i)%vnne == nsegmax) tauTL(I) = tauTL(seg(i)%vnno)
    endif


    !*** Imposed stress field inside a cubic region used for calculation made in collaboration with P. Geantil
    !   If (mode_deformation == IHUIT) then
    !       ROX = int(Ro(1,I)/(modur(1)/dimFE_Sig(1)))+IUN
    !       ROY = int(Ro(2,I)/(modur(2)/dimFE_Sig(2)))+IUN
    !       ROZ = int(Ro(3,I)/(modur(3)/dimFE_Sig(3)))+IUN
    !       if (ROX>=2 .and. ROX<=6) then
    !         if (ROY>=2 .and. ROY<=6) then
    !           if (ROZ>=2 .and. ROZ<=6) then
    !
    !       TauTOT(I) = TauTL(I)+TauAPP(I)
    !
    !           else
    !             TauTOT(I) =TauTL(I)+TauAPP(I)+TauINT(I)
    !           endif
    !         else
    !           TauTOT(I) =TauTL(I)+TauAPP(I)+TauINT(I)
    !         endif
    !       else
    !         TauTOT(I) =TauTL(I)+TauAPP(I)+TauINT(I)
    !       endif
    !   else
    !***

    !=================================================================================================
    ! Calculation of the total (effective) stress
    if (seg(i)%unload) then

      ! The particular case of segments part of an infinite line and free of junctions
      ! Thus line is forced to been pinned with one immobile segment par of the free inf line
      TauTOT(i) = 0

    else

      TauTOT(I) = TauTL(I) + TauAPP(I) + TauINT(I)      ! The effective stress
!       TauTOT(I) = TauTL(I) + TauAPP(I)                  ! The spaghetti stress mode
!       TauTOT(I) = TauAPP(I)							  ! The ghost stress mode
    endif
!     if (syseg(seg(i)%veclin) /= 1) TauTOT(I) = 0.     ! all segments except segments in sys 1 are immobile
    !=================================================================================================

    if (TauTOT(i) > Tau_limite) then
      singulier(i)  = .true.
      TauTOT(i)     = Tau_limite * Dsign(un,TauTOT(i))
    endif

    ! The special case of line tension relaxation
    if (kk <= relax_tl) TauTOT(I) =TauTL(I)

    !***
    !endif
    !***

    ! Check if the stress screening works
    if (abs(TAUAPP(I)) > TAU_LIMITE .or.    &
        abs(TAUINT(I)) > TAU_LIMITE .or.    &
        abs( TAUTL(I)) > TAU_LIMITE) then
      write(*,*) "Critical Shear Stresses problem at step:", KK
      write(*,fmt='("--> Calculation of shear stresses for Segment i= ", I7)') i
      write(*,fmt='("TauTOT     : ",F20.15)') TauTOT(i)
      write(*,fmt='("TauAPP     : ",F20.15)') TauAPP(i)
      write(*,fmt='("TauTL      : ",F20.15)') TauTL(i)
      write(*,fmt='("TauINT     : ",F20.15)') TauINT(i)
      write(*,fmt='("TauMon     : ",F20.15)') TauMon(i)
      write(*,fmt='("TauLimite  : ",F20.15)') Tau_limite
      stop
    endif

  endif

#ifdef PA
  endif
#endif

ENDDO b3

NSEGM=NSEGM+PLUSSeg
PLUSSeg = IZERO


if (kkdebug) then
  write(379,*)  " Resultats de la procedure force  ......................"
  write(379,*)  " Matrice de contraintes : " ;  print *, " "
  write(379,'(3(F10.7,2x))') sigapp(1,1:3)
  write(379,'(3(F10.7,2x))') sigapp(2,1:3)
  write(379,'(3(F10.7,2x))') sigapp(3,1:3)
  do i= 1,nsegm
    if (AB_length(i) /= 0) then
      carac = VLseg(i)
      tmp = 1.D-6*Xmu
      write(379,'("S:",I5,"; Tap= ",F13.2,";Tint=",F12.2,"; &
            & TL=",F12.2," ; T= ",F13.2," ; S= ",L1," W=",I2," Loi=",I2," AngVis=",F7.2)')  &
            i, TAUAPP(I)*tmp,tauINT(I)*tmp, TAUTL(I)*tmp,                                   &
            TAUTOT(I)*tmp, Singulier(I), seg(i)%wait, LoiSeg(i), seg(i)%anglevis
    endif
  enddo
endif

end subroutine FORCE

!###########################################################################
!#  Region is the subroutine proceeding with the domain decomposition used #
!#  in the multipoles algo and the obstacle search. Outputs of REGION are: #
!#    NsegBoite(:) nb of seg in each domain                                #
!#    IndexBoite(:)%listeseg(:) : the list of segs in each domains         #
!#    RO(I) coordinates of the segs middle                                 #
!###########################################################################
subroutine REGION

implicit none

integer(kind=DPI)         :: i, j, ib, jb, kb, wait
integer(kind=DPI)         :: ListeSeg_I(nsegmax)
integer                   :: ListeSeg_size
integer                   :: ListeSegTest
integer                   :: NSegBoite_I    !< Number of segment in each ListeSeg tab elements
integer                   :: ListeSegRange  !< The range value used in the ListeSeg tab allocation
integer                   :: MinRange       !< The size limit used to reduce the ListeSeg tab dim
integer                   :: MaxRange       !< The upper size limit used to increase the ListeSeg tab dim
integer(kind=DPI)         :: X(3)

real(kind=DP)             :: invtailleboite1, invtailleboite2, invtailleboite3

#ifdef PA
! For the biphase calculation
integer(kind=DPI)                 :: J1,J2,J3,L
integer(kind=DPI),allocatable     :: liste_RO_temp(:,:)
#endif

! Initialisation
NsegBoite(:) = izero
invtailleboite1 = un/tailleBoite(1)
invtailleboite2 = un/tailleBoite(2)
invtailleboite3 = un/tailleBoite(3)

! Attention : ici on remet les O des seg dans la cellule de ref pour les tests de localisation
seg(1:nsegm)%O(1) = modulo(seg(1:nsegm)%O(1),modur(1))
seg(1:nsegm)%O(2) = modulo(seg(1:nsegm)%O(2),modur(2))
seg(1:nsegm)%O(3) = modulo(seg(1:nsegm)%O(3),modur(3))

! The definition of LenghtLimit must be redefined at every step since the size of Greengard domains can change
lengthlimit = izero          ! LegnhtLimit must not be used without the Greengard algorithm
if(Nboites /= iun) lengthlimit = ideux * int (minval(Tailleboite(1:3)))

! test on the shape of the shape of Greengard domains
if (maxval(Tailleboite(1:3)) > (2 * minval(Tailleboite(1:3))) .and. Nboites /= 1) then
  stop "Problem ! Greengard domains are assumed to have approximately cubic shape"
endif

! definition of NbGroseg and Gorseg(:)
Nb_groseg = izero
groseg(:) = izero
LongSeg(:) = .false.

#ifdef MDC
nsegmdcstress=IZERO
#endif

! Loop on the segment to define their distribution in the domains
do I =  1, NSEGM
  wait = IZERO
  ! Such segments have been eliminated in the first discretisation step and must not be include in the Greengard domains
  if (i == seg(I)%voiso .and. i == seg(i)%voise) CYCLE

  VLseg(I)    = seg(I)%veclin
  A_COOR(:,I) = SEG(I)%O(:)
  X(:) = (AB_LENGTH(I)*bveclin(:,VLseg(I))) / ideux
  X(:) =  A_COOR(:,I) + X(:)
  RO(:,I) =  X(:)
  if (AB_LENGTH(I) == IUN) then
    if (seg(i)%surface == IUN) then
      X(:) =  A_COOR(:,I) + bveclin(:,VLseg(I))
      RO(:,I) =  X(:)
    endif
    if (seg(i)%surface == IDEUX) then
      X(:) = A_COOR(:,I)
      RO(:,I) = X(:)
    endif
  endif

!#ifdef MDC
!  if (AB_LENGTH(I) /= IZERO) then
!    nsegmdcstress = nsegmdcstress + IUN
!    listsegmdcstress(nsegmdcstress) = i
!  endif
!#endif

  ! When building the region list, we update the info on quasi-immobile segments not included in the force calculation
  if (Period /= IZERO) then
      wait = modulo(seg(i)%wait,Period)
      seg(i)%wait = wait
  endif

#ifdef MDC
  if(zcall) then
    if (AB_LENGTH(I) /= IZERO &
        .and. .not. tabjonc(I) .and. tabgd(I) < iun) then
      nsegmdcstress = nsegmdcstress + IUN
      listsegmdcstress(nsegmdcstress) = i
    endif
  endif
#endif

  if (Nboites /= iun) then
    ! The definition of LongSeg = segments which cannot be taken into account with the Greengard algorithm
    if( (seg(I)%norme * normlin(VLseg(I)))  > lengthlimit) then
      NB_Groseg         = NB_Groseg + IUN
      Groseg(NB_groseg) = I
      LongSeg(I)        = .true.
!       ! We check that groseg are only junction segments
!       if (.not. seg(i)%jonc) print*,NB_Groseg,(seg(I)%norme * normlin(VLseg(I)))*avalue,lengthlimit*avalue
    endif

    ! Be aware that the function floor must be prefered to function int
    ! in the modulo operations
    IB = modulo(floor(X(1)*invtailleBoite1),int(NboitesX))+iun
    JB = modulo(floor(X(2)*invtailleBoite2),int(NboitesY))+iun
    KB = modulo(floor(X(3)*invtailleBoite3),int(NboitesZ))+iun
    J = B3D_1D(IB,JB,KB)

    if( NsegBoite(J) + IUN > nsegmax_boite) then
      print *," nsegmax_boite =",nsegmax_boite, " nsegm =",nsegm
      print *,"J,IB,JB,KB,",J,IB,JB,KB
      print *,"NsegBoite(IB,JB,KB)",NsegBoite(J)
      print *,"erreur REGION : debordement de tableau listeboite kk=", kk
      stop
    endif

    IBOITE(I) = J
    NsegBoite(J) =  NsegBoite(J) + IUN
    IndexBoite(J)%ListeSeg(NsegBoite(J)) = I

  else

    ! WHEN NO MULTIDOMAINE applied, all segments are tested every where
    NsegBoite(IUN) = NsegBoite(IUN) + IUN
    IndexBoite(IUN)%ListeSeg(NsegBoite(IUN)) = i
    IBOITE(I) = IUN

  endif

enddo

if(kkdebug .or. Nb_groseg*10 > nsegm) write(379,*)  " % of long segments :", (NB_groseg*100 / Nsegm)

if(KKdebug) then
  write(379,*)  "NB  long segments :", groseg(1:nb_groseg)
  write(379,*)  "NB  long segment % jonc :", tabjonc(groseg(1:nb_groseg))
  write(379,*)  " info sur les boites dans REGION"
  write(379,'("Number of multipole domains ",I6," (",I2,",",I2,",",I2,")  &
        & Nseg=",I6)')NBoites,NBoitesX,NBoitesY,NBoitesZ,nsegm
  write(379,'("Size of multipole domains in AVALUE (",F9.3,",",F9.3,",",F9.3,")")') TailleBoite(1:3)
endif

# ifdef PA

! Nb_phase=2: the segments close to a boundary are memorized in liste_RO with its
! corresponding box number. Those segments index is saved in ID_RO
!
! The structure of list_RO is as follow:
! liste_RO(1,I)   : index of the box containing segment I
! liste_RO(2:4,I) : center coordinates of segment I
! I is comprised between 1 and Nb_seg ( the number of segments in boxes close to the boundary)
if (Nb_phase == 2) then

  ! allocation of temporary lists with length proprotionnal to NSegm and
  ! NBoites because the final lists need the knowledge of the number of full boxes
  ! close to the boundary and the number of segments contained by these boxes.
  allocate(liste_RO_temp(4,Nsegm))
  allocate(ID_RO(Nsegm))

  ! Initialization
  liste_RO_temp(:,:)  = IZERO
  ID_RO(:)            = IZERO

  I = IZERO   ! Total number of segments (in all boxes) close to a boundary

  do J1 = 1,NBoites_boundary

    j2 = list_boxes(J1)     ! Index of the boxes close to a boundary

    if (NsegBoite(J2)/=0) then

      L = IZERO         ! L is the number of non zero segment in a box

      do J3 = 1, NsegBoite(J2)

        J = IndexBoite(J2)%ListeSeg(J3)

        if (AB_length(J) /= 0) then

          I = I + IUN
          L = L + IUN

          ! List of effective segments to consider
          liste_RO_temp(1,I) = J2           ! Index of the boxes close to a boundary

          liste_RO_temp(2,I) = RO(1,J)      ! Center of the segment
          liste_RO_temp(3,I) = RO(2,J)
          liste_RO_temp(4,I) = RO(3,J)

          ID_RO(I) = J                      ! Index of the segments used for SIGMA_JOINT calculation

        endif

      enddo

    endif

  enddo

  Nb_seg = I          ! Total number of segments close to a boundary of non zero length segments

  ! allocation of minimum lists with length proprotionnal to Nb_seg
  allocate(liste_RO(4,Nb_seg))

  liste_RO(:,1:Nb_seg) = liste_RO_temp(:,1:Nb_seg)

  deallocate(liste_RO_temp)

endif
# endif

if (Nboites /= iun) then
  ! To optimize memory size and memory accesses, the dimension of the huge array IndexBoite(:)%ListeSeg(:)
  ! is re-allocate to fit as much as possible the real number of segments in each domain
  ! We first check that the ListeSegRange value is still correctly evaluated for
  ! the runnning simulation
  ListeSegTest = 0

  do i = 1, NBoites

    ListeSeg_size = size(IndexBoite(i)%ListeSeg)

    if (ListeSeg_size /= 1) then   ! ListeSeg_size = 1 for domains outside of the effective dislocated volume

      NsegBoite_I   = int(NsegBoite(i))

      ListeSegRange = ListeSeg_size    ! The initial dim of the ListeSeg tab element
      MinRange = ListeSeg_size / 3     ! When possible we decrease ListeSeg if its
                                       ! dim is more than 3 time the number of segments
      MaxRange = NsegBoite_I * 2       ! When possible we increase ListeSeg if its
                                       ! dim is less than 2 time the number of segments

      ! After new allocation of the boites, the optimum dimension of the boite must be set again
      if (ListeSegRange == Nsegmax_boite) then

        ListeSeg_size = ListeSegRangeIni
        if (ListeSeg_size < MaxRange) ListeSeg_size = int(real(NsegBoite_I) * (2. + ListeSegRangeFac))

        ! We save the elements of ListeSeg
        ListeSeg_I(1:NsegBoite_I) = IndexBoite(i)%ListeSeg(1:NsegBoite_I)
        ! We deallocate the ListeSeg tab
        deallocate(IndexBoite(i)%ListeSeg)
        ! The new dim of the ListeSeg is calculated (50 is added to avoid problems in domains with very few segments)
        allocate(IndexBoite(i)%ListeSeg(1:ListeSeg_size))
        ! The elements of ListSeg are recovered
        IndexBoite(i)%ListeSeg(1:NsegBoite_I)                = ListeSeg_I(1:NsegBoite_I)
        IndexBoite(i)%ListeSeg(NsegBoite_I+1:ListeSeg_size)  = izero

        ! print*,kk,i,'initialization',ListeSegRange,ListeSeg_size,NsegBoite_I
        ! read *
      endif

      ! print*, '---',i,MinRange,MaxRange,NsegBoite_I,ListeSegRange,   &
      !     int(real(NsegBoite_I) * (2. + ListeSegRangeFac)),ListeSegRangeFac

      ! We compare the number of segments in the domain (boite) and the dim of the ListeSeg tab
      ! The factor ListeSegRangeFac is used to avoid many repeated allocation
      ! If needed we decrease the ListeSeg tab dim
      if (NsegBoite_I < MinRange .and. ListeSegRange > ListeSegRangeIni) then

        ListeSeg_size = int(real(NsegBoite_I) * (3. - ListeSegRangeFac))        !  The default solution
        if (ListeSeg_size < ListeSegRangeIni) ListeSeg_size = ListeSegRangeIni  !  The minimum value

        ! We save the elements of ListeSeg
        ListeSeg_I(1:NsegBoite_I) = IndexBoite(i)%ListeSeg(1:NsegBoite_I)
        ! We deallocate the ListeSeg tab
        deallocate(IndexBoite(i)%ListeSeg)
        ! The new dim of the ListeSeg is calculated (50 is added to avoid problems in domains with very few segments)
        allocate(IndexBoite(i)%ListeSeg(1:ListeSeg_size))
        ! The elements of ListSeg are recovered
        IndexBoite(i)%ListeSeg(1:NsegBoite_I)                = ListeSeg_I(1:NsegBoite_I)
        IndexBoite(i)%ListeSeg(NsegBoite_I+1:ListeSeg_size)  = izero

        ! print*,kk,i,'decrease',ListeSegRange,ListeSeg_size,NsegBoite_I
        ! read *

      ! if the number of segments in the boite is too close to the boite dimension, the boite dimension
      ! must be increased
      elseif (ListeSegRange < MaxRange) then

        ListeSeg_size = int(real(NsegBoite_I) * (2. + ListeSegRangeFac))

        ! We save the elements of ListeSeg
        ListeSeg_I(1:NsegBoite_I) = IndexBoite(i)%ListeSeg(1:NsegBoite_I)
        ! We deallocate the ListeSeg tab
        deallocate(IndexBoite(i)%ListeSeg)
        ! The new dim of the ListeSeg is calculated (50 is added to avoid problems in domains with very few segments)
        allocate(IndexBoite(i)%ListeSeg(1:ListeSeg_size))
        ! The elements of ListSeg are recovered
        IndexBoite(i)%ListeSeg(1:NsegBoite_I)                = ListeSeg_I(1:NsegBoite_I)
        IndexBoite(i)%ListeSeg(NsegBoite_I+1:ListeSeg_size)  = izero

        ! print*,kk,i,'increase',ListeSegRange,ListeSeg_size,NsegBoite_I
        ! read *

      endif

    endif

  !   if (ListeSeg_size > ListeSegTest) ListeSegTest = ListeSeg_size

  enddo

endif

!procedure to reallocate sweptsurfdata array!!!

#ifdef MDC

#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
#endif

    if (zcall) then

       if (nsegmdcstress > nsegm_INI) then
         write(379,*) 'increasing size of the sweptsurfdata array'
         write(379,*) 'nsegm = ', nsegm
         write(379,*) 'nsegm_INI old= ', nsegm_INI
         deallocate(sweptsurfdata)
         nsegm_INI = int(nsegmdcstress*1.6,DPI)
         write(379,*) 'nsegm_INI new= ', nsegm_INI
         if (nsegm_INI < NSEGMAX) then
           allocate(sweptsurfdata(nsegm_INI*mdc_timestep_ratio*16+1))
           write(379,*) 'allocation using NSEGM_INI', nsegm_INI
         else
           allocate(sweptsurfdata(nsegmax*mdc_timestep_ratio*16+1))
           nsegm_INI=nsegmax
           write(379,*) 'allocation using NSEGMAX', nsegmax
           write(379,*) 'NSEGM_INI == NSEGMAX'
         endif
         sweptsurfdata(:)=0
       endif
    endif

#if defined(PA) && defined(MDC)
  endif
#endif

#endif


! print*,'---->', maxval(NsegBoite(1:NBoites)), ListeSegTest, Nsegmax_boite
! read *

end subroutine REGION

!###########################################################################
!# Subroutine SIGMA_INT_CP : Calcul de la somme des champ des segemnts     #
!# voisins: dans la meme boite et dans les boites voisines                 #
!###########################################################################
subroutine SIGMA_INT_CP

implicit none

integer (kind=DPI)                      :: I,IB,JB,indiceI,indiceJ,iini,ivli
integer (kind=DPI)                      :: I1O,I1E,NsegTester,listeSegTester(NSEGM),DCFJ_norm
integer (kind=DPI),dimension(3)         :: R,Tr

real    (kind=DP),dimension(3)          :: RA,RB,ROOmod,TrABmod
real    (kind=DP),dimension(3)          :: Tbis,Bbis,BPVTbis
real    (kind=DP),dimension(6)          :: SIGINTERN
real    (kind=DP)                       :: realnorme
logical                                 :: condition

!variables for determine P1,P2,P3 used in Line Tension regularisation
integer (kind=DPI),dimension(3)         :: P1Iini,P2Iini,P3Iini,Oi,Ei
integer(kind=DPI), dimension(3)         :: trAB,P1temp,P3temp
integer(kind=DPI)                       :: ii,jj,compteurP1,compteurP3,vnn,halflength
integer(kind=DPI)                       :: listP1Iini(11),listP3Iini(11),ROOP1Iini(3,11),ROOP3Iini(3,11)

real(kind=DP), dimension(3)             :: P3P2,P1P2
real(kind=DP)                           :: normP3P2,normP1P2,cosP1P2vP3P2,invcosP1P2vP3P2

! variable for the calculation of non-singular dislocation stress field (NSDSF)
real(kind=DP) :: d2,s1,s2,a2_p_d2,inv_a2_p_d2
real(kind=DP) :: s1Ra3inv,Ra1_3inv,Ra1_inv
real(kind=DP) :: s2Ra3inv,Ra2_3inv,Ra2_inv
real(kind=DP) :: s1_03,s1_13,s1_05,s1_15,s1_25
real(kind=DP) :: s2_03,s2_13,s2_05,s2_15,s2_25
real(kind=DP) :: s_03,s_13,s_05,s_15,s_25
real(kind=DP) :: m4pn,mn4pn,a2m8p
real(kind=DP) :: d_pv_b_ps_t,commun,R_ps_t,Ra1,Ra2
real(kind=DP),dimension(3) :: cai_d,d_pv_b,t_pv_b
real(kind=DP),dimension(6) :: t_pt_t
real(kind=DP),dimension(6) :: d_pt_d,t_pt_d,t_pt_t_pv_b,d_pt_t_pv_b,t_pt_d_pv_b
real(kind=DP),dimension(6) :: I_03,I_13,I_05,I_15,I_25

! Optimisation vriables
real(kind=DP)      :: aX            !<
real(kind=DP)      :: aY            !<
real(kind=DP)      :: aZ            !<
real(kind=DP)      :: modurX        !<
real(kind=DP)      :: modurY        !<
real(kind=DP)      :: modurZ        !<
real(kind=DP)      :: HalfModurX    !<
real(kind=DP)      :: HalfModurY    !<
real(kind=DP)      :: HalfModurZ    !<
real(kind=DP)      :: invmodurX     !<
real(kind=DP)      :: invmodurY     !<
real(kind=DP)      :: invmodurZ     !<

! variable intermediare pour optimisation
integer(kind=DPI)                       :: IX,IY,IZ,BX,BY,BZ,I1,J1,K1

#ifdef PA
! Variables pour la parallelisation
REAL(kind=DP),allocatable     :: Vitesse_p(:)
INTEGER(kind=DPI),allocatable :: COUTSEG(:)
INTEGER(kind=DPI),allocatable :: Charge(:)
INTEGER(kind=DPI),allocatable :: NSEGBP(:)
INTEGER(kind=DPI),allocatable :: BDEB(:)
INTEGER(kind=DPI),allocatable :: SEGDEB(:)
INTEGER(kind=DPI),allocatable :: NSEG(:)

REAL(KIND=DP)        :: VTOTAL,INVVTOTAL
REAL(KIND=DP)        :: prep,prep0,prep1,calc
INTEGER(kind=DPI)    :: CHRG,NCOUPLESEG,INTP,INTP2
INTEGER(kind=DPI)    :: BOITED,SEGD,BOITEF,SEGF

! Used for the calculation of SIGJOINT
INTEGER(kind=DPI)           :: RANG_DEST,SEGD_JOINT,SEGF_JOINT,IS
INTEGER(kind=DPI)           :: Nb_seg_stored ! number of segments in the total box list close to the boundary
REAL(kind=DP)               :: inv_homothetie
real(kind=DP),allocatable   :: SIGJOINT(:,:), SIGJOINT_Sum(:,:), SIGJOINT_Recv(:,:)
#endif

#ifdef MDC
real(kind=DP), dimension(11)        :: buf
real(kind=DP)                       :: dt2,dn2,dt
real(kind=DP)                       :: dist2,ExtR2
logical                             :: ied
real(kind=DP) :: a2_p_d2_MDC,inv_a2_p_d2_MDC
real(kind=DP) :: s1Ra3inv_MDC,Ra1_3inv_MDC,Ra1_inv_MDC
real(kind=DP) :: s2Ra3inv_MDC,Ra2_3inv_MDC,Ra2_inv_MDC
real(kind=DP) :: s1_03_MDC,s1_13_MDC,s1_05_MDC,s1_15_MDC,s1_25_MDC
real(kind=DP) :: s2_03_MDC,s2_13_MDC,s2_05_MDC,s2_15_MDC,s2_25_MDC
real(kind=DP) :: s_03_MDC,s_13_MDC,s_05_MDC,s_15_MDC,s_25_MDC
real(kind=DP) :: a2m8p_MDC,Ra1_MDC,Ra2_MDC
real(kind=DP),dimension(6) :: I_05_MDC,I_15_MDC,SIGINTNS
real(kind=DP)       :: invtailleboite1,invtailleboite2,invtailleboite3
integer(kind=DPI) :: IDEBX,IDEBY,IDEBZ
#endif

! Maybe we should change IDEB name for the different uses
#ifndef MDC
integer(kind=DPI)                       :: IDEB
#endif

#if defined(MDC) && defined(PA)
integer(kind=DPI)                       :: IDEB
#endif


!initialization of NSDSF parameters
m4pn=0.25d0/(1.0d0-dpoiss)
mn4pn=m4pn*dpoiss
a2m8p=bspread2*0.125d0

! Initialization of optimization variables
modurX      = real(modur(1),DP)
modurY      = real(modur(2),DP)
modurZ      = real(modur(3),DP)
HalfModurX  = HalfModur(1)
HalfModurY  = HalfModur(2)
HalfModurZ  = HalfModur(3)
invmodurX = UN/modurX
invmodurY = UN/modurY
invmodurZ = UN/modurZ

#ifdef MDC
a2m8p_MDC=halfthickness2*0.125d0
if (Nboites /= IUN) then

  invtailleboite1 = un/tailleBoite(1)
  invtailleboite2 = un/tailleBoite(2)
  invtailleboite3 = un/tailleBoite(3)

  IDEBX = nint((SR_cutoff*invtailleboite1),DPI)
  IDEBY = nint((SR_cutoff*invtailleboite2),DPI)
  IDEBZ = nint((SR_cutoff*invtailleboite3),DPI)

  if ((2*IDEBX+1)>=NboitesX .or. (2*IDEBY+1)>=NboitesY .or. (2*IDEBZ+1)>=NboitesZ ) &
    stop 'Inclusion halfthickness > half simulation box'

else

  IDEBX = izero
  IDEBY = izero
  IDEBZ = izero

endif
#endif

#ifdef PA

prep0 = real(MPI_WTIME(),DP)   ! La clock

!**************!
!  Allocation  !
!**************!
allocate(Vitesse_p(TAILLE_PAIRIMPAIR))
allocate(NSEG(TAILLE_PAIRIMPAIR))
allocate(Charge(TAILLE_PAIRIMPAIR))
allocate(BDEB(TAILLE_PAIRIMPAIR))
allocate(SEGDEB(TAILLE_PAIRIMPAIR))
allocate(COUTSEG(NBOITES))
allocate(NSEGBP(NBOITES))

! Take care !!!! For parallel computation SIG_CP must be initialized at each steps and for each procs
SIG_CP(1:3,1:3,1:NSEGM)=zero

!***********************************************************!
! The reference velocity of the cores in parallel computing
!***********************************************************!
Vitesse_p(1:TAILLE_PAIRIMPAIR)=UN/real(TAILLE_PAIRIMPAIR,DP)   ! We guess that the speed of all the procs is identical
VTotal = UN

!***************************************!
! calculation dispatch between the procs
!***************************************!
 COUTSEG(1:NBOITES)=0

! The number of operations to be made for each segments is function of its domain
! with the multi-poles algorithm
if(nboites > un) then

  do IB = 1, NBOITES
    IX=B1D_3D(IB,1)-IUN
    IY=B1D_3D(IB,2)-IUN
    IZ=B1D_3D(IB,3)-IUN
    INTP2=0
    do I1 = -IUN , IUN
      do J1 = -IUN , IUN
        do K1 =  -IUN , IUN
          BX = modulo(IX+I1,NBoitesX)+IUN
          BY = modulo(IY+J1,NBoitesY)+IUN
          BZ = modulo(IZ+K1,NBoitesZ)+IUN
          INTP  = BX + ((BY-1)+ (BZ-1)*NBoitesY)*NBoitesX
          INTP2 = INTP2 + NSEGBOITE(INTP)
        enddo
      enddo
    enddo

    COUTSEG(IB)=INTP2   ! nombre de segments voisin vu par un segment dans la boite IB

  enddo

else

  COUTSEG(IUN)=NSEGM   ! if only on domain

endif

! Total number of operation to proceed
NCOUPLESEG = DOT_PRODUCT(COUTSEG(1:NBOITES),NsegBoite(1:NBOITES))

! The unit for parallel calculation
INVVTOTAL = NCOUPLESEG / VTOTAL

 Charge(1:TAILLE_PAIRIMPAIR) = nint(INVVTOTAL * Vitesse_p(1:TAILLE_PAIRIMPAIR),DPI)
 Charge(TAILLE_PAIRIMPAIR) = NCOUPLESEG - SUM(Charge(1:(TAILLE_PAIRIMPAIR - 1)))

! Dispatch of the work to do by each process
IB=1
INTP=0
NSEG(:)=0
NSEGBP(:)=NSEGBoite(:)
BDEB(1)=1
SEGDEB(1)=1

do I = 1, (TAILLE_PAIRIMPAIR - 1)

  CHRG=Charge(I)

  do while(CHRG > 0 .AND. CHRG >= (NSEGBP(IB)*COUTSEG(IB)))
    CHRG = CHRG-NSEGBP(IB)*COUTSEG(IB)
    NSEG(I) = NSEG(I)+NSEGBP(IB)
    NSEGBP(IB) = 0
    IB = IB + 1
  enddo

  if (COUTSEG(IB)/=0) then
    INTP = modulo(CHRG,COUTSEG(IB)) ! le reste de charge ne faisant pas le calcul d un segment entier
    INTP2 = (CHRG-INTP)/COUTSEG(IB) ! nombre de segment de IB encore possible de traiter
  else
    INTP = 0
    INTP2 = 0
  endif

  Charge(I) = Charge(I) - INTP
  Charge(I+1) = Charge(I+1) + INTP
  NSEG(I) = NSEG(I) + INTP2
  NSEGBP(IB) = NSEGBP(IB) - INTP2
  BDEB(I+1) = IB

  if (BDEB(I) /= BDEB(I+1)) then
    SEGDEB(I+1) = INTP2 + 1
  else
    SEGDEB(I+1) = SEGDEB(I) + INTP2 + 1
  endif

enddo

!******************!
!   Initialization
!******************!
MONSEG(1:NSEGM) = .FALSE.
BOITED = BDEB(MON_RANG_PAIRIMPAIR + 1)                    ! the starting box for the procs
SEGD   = SEGDEB(MON_RANG_PAIRIMPAIR + 1)                  ! the starting segment for the procs

if ((MON_RANG_PAIRIMPAIR + 1) /= TAILLE_PAIRIMPAIR) then
  BOITEF = BDEB(MON_RANG_PAIRIMPAIR + 2)                  ! the ending box for the procs
  SEGF   = SEGDEB(MON_RANG_PAIRIMPAIR + 2) - 1            ! the ending segment for the procs
else
  BOITEF = Nboites
  SEGF   = NsegBoite(Nboites) + 1
endif

prep1 = real(MPI_WTIME(),DP)       ! Time at the end of the parallel preparation
prep = (prep1-prep0)            ! Time dedicated to prepare the parallel computation

#endif

! boucle sur les boites
LIini:do IB = 1,Nboites

    if (NsegBoite(IB) == izero) cycle LIini

#ifdef PA
    ! debut du domaine ou le process va travailler
    if (IB < BOITED) cycle
#endif

    ! boucle sur les seg dans la boite IX,IY,IZ
    bIini: do IndiceI = 1,NsegBoite(IB)

        Iini = INdexBoite(IB)%ListeSeg(IndiceI)
        IVLI = VLseg(iini)

#ifdef PA
        ! debut et fin exact du domain de travail du process
        IF ( IB >  BOITEF)                          exit LIini
        IF ((IB == BOITEF) .AND. (IndiceI > SEGF))  exit LIini
        IF ((IB == BOITED) .AND. (IndiceI < SEGD))  cycle
        MONSEG(Iini)=.TRUE.
#endif

        !****************************** On ne calule pas la contrainte si :

        ! for waiting segments no stress calculation
        ! exept for wait=0,1 2
        if(seg(Iini)%wait> IDEUX) cycle bIini

        !  Le segment est de longueur nulle,
        if(AB_length(Iini) == IZERO) cycle bIini

        ! le segment est JONC ET GD rotule
        if(tabjonc(iini) .or. tabgd(iini) > izero) cycle bIini

        !  Le segment fait partie de la foret,
        if(Iini.le.IFORET) cycle bIini

        !  L'une des 2 extremites du segment est ancre

        if((SEG(Iini)%VOISO /= NSEGMAX) .and. (SEG(Iini)%VOISE /= NSEGMAX)) then

          !!!!!!!! CALCulation of P1 P2 and P3 usefull for line tension regularisation !!!!!!!

          i1o = seg(iini)%vnno
          i1e = seg(iini)%vnne

          ! DETERMINATION of the points needed for the calculation of the local radius of curvature:
          ! P2 is always the center of I
          P2Iini(1:3) = A_coor(1:3,Iini) + (AB_length(Iini) * bveclin(1:3,IVLI)) / ideux

          ROOP1Iini(:,:) = IZERO
          listP1Iini(:)  = IZERO
          ROOP3Iini(:,:) = IZERO
          listP3Iini(:)  = IZERO

          ! P1Iini is :
          ! - the origin of I if I1o is pinning point
          ! - the center of I1o if I1o is a real segment (not a part of a long segment beeing descretised)
          ! - the center of the effective neighbour before discretisation

          if (i1o  == nsegmax .or. tabjonc(i1o)) then

              P1Iini(1:3) = A_coor(1:3,Iini)
              compteurP1  = iun
              ROOP1Iini(1:3,compteurP1) = A_coor(1:3,Iini) - P2Iini(1:3)
              listP1Iini(compteurP1) = Iini

          else

              ! Oi : is the beginning of the effectiv I1o
              ! Ei : is the end of the effectiv I1o
              Ei(1:3) = A_coor(1:3,Iini)

              ! to avoid PBC problem all coordinate are calculated relatively to the origine of I
              Tr(1:3) =   AB_length(I1o) * bveclin(1:3,VLseg(i1o)) / IDEUX
              Oi(1:3) = Ei(1:3) - Tr(1:3) * IDEUX

              compteurP1 = IUN

              if (RLOCAL2 < ZERO) then

                ! P1Iini is the center of the vnn segment
                P1Iini(1:3) = Oi(1:3) + Tr(1:3)
                ROOP1Iini(1:3,compteurP1) = (Oi(1:3) + Tr(1:3)) - P2Iini(1:3)
                listP1Iini(compteurP1) = i1o

              else

                ! The vector connecting P2Iini and the center of i1o
                ROOP1Iini(1:3,compteurP1) = (Oi(1:3) + Tr(1:3)) - P2Iini(1:3)

                ! Do we need to look for more distant segment i10 to regularize the line tension curvature ?
                if (SUM(real(ROOP1Iini(:,compteurP1),DP)*real(ROOP1Iini(:,compteurP1),DP)) < Rlocal2) then

                  !The half length of i1o segment
                  TrAB(1:3) = Tr(1:3)

                  ! Ra and Rb are vectors between P2Iini and the extremities of i1o
                  RA(1:3)=real(ROOP1Iini(1:3,compteurP1)+TrAB(1:3),DP)
                  RB(1:3)=real(ROOP1Iini(1:3,compteurP1)-TrAB(1:3),DP)

                  vnn = i1o
                  listP1Iini(compteurP1) = i1o

                  ! exact calculation of the line tension impose that the local radius of curvature is calculated
                  ! with points P1 and P3 not very close from P2Iini. A minimum distance of Rlocal is imposed!
                  do while (seg(vnn)%vnno /= nsegmax .and.          &
                            .not. tabjonc(vnn)       .and.          &
                            tabgd(vnn) < iun         .and.          &
                            SUM(RA*RA) < Rlocal2     .and.          &
                            SUM(RB*RB) < Rlocal2)

                    ! The next i1o segment along the dislocation line is tested
                    if(compteurP1 > idix) exit
                    vnn = seg(vnn)%vnno
                    compteurP1 = compteurP1 + IUN

                    ! The vector connecting P2Iini and the center of previous vnn
                    Tr(1:3) = TrAB(1:3)

                    !The half length of the new vnn segment
                    TrAB(1:3) = (AB_length(vnn)*bveclin(1:3,VLseg(vnn))) / ideux

                    !The center of the new vnn segment
                    ROOP1Iini(1:3,compteurP1) = ROOP1Iini(1:3,CompteurP1-IUN) - (Tr(1:3) + TrAB(1:3))
                    listP1Iini(compteurP1) = vnn
                    ! Ra and Rb are vectors between P2Iini and the extremities of new i1o
                    RA(:)=real(ROOP1Iini(:,compteurP1)+TrAB(:),DP)
                    RB(:)=real(ROOP1Iini(:,compteurP1)-TrAB(:),DP)

                  enddo

                  P1Iini(1:3) = P2Iini(1:3) + ROOP1Iini(1:3,compteurP1) ! P1 is the center of the vnn segment

                  ! A special case for segments touching the surface. Since we don t do the obstacle research
                  ! segments overlap may happen and induce wrong ROO calculation
                  if (SUM(ROOP1Iini(1:3,compteurP1)*ROOP1Iini(1:3,compteurP1)) == 0 &
                     .and. seg(Iini)%surface > 0) then

                    P1Iini(1:3) = A_coor(1:3,Iini)
                    compteurP1=IUN
                    ROOP1Iini(1:3,compteurP1) = A_coor(1:3,Iini) - P2Iini(1:3)
                    listP1Iini(compteurP1) = Iini

                  endif

                endif

              endif

              ! Now we look for the problem of long straight lines cut in many segments,
              ! straight lines section must be treated as one long segment, whatever rlocal
              if (compteurP1 == IUN) then

                ! The initial vnn is redefined
                vnn = i1o
                listP1Iini(compteurP1) = i1o

                do while(seg(vnn)%vnno /= nsegmax .and. VLseg(seg(vnn)%vnno) == VLseg(i1o))

                  if(compteurP1 > idix) exit
                  vnn = seg(vnn)%vnno
                  compteurP1 = compteurP1 + IUN

                  !The half length of the new vnn segment
                  Oi(1:3) =  Oi(1:3) - AB_length(vnn) * bveclin(1:3,VLseg(vnn))

                  !The center of the new vnn segment
                  ROOP1Iini(1:3,compteurP1) = (Oi(1:3)+Ei(1:3))/IDEUX - P2Iini(1:3)
                  listP1Iini(compteurP1) = vnn

                enddo

                P1Iini(1:3) = P2Iini(1:3) + ROOP1Iini(1:3,compteurP1)

              endif

          endif  ! End of the tests on the P1 side


          ! P3Iini is :
          ! - the end     of I if I1e is pinning point
          ! - the center  of I1e if I1e is a real segment (not a part of a long segment beeing descretised)
          ! - the center of the effective neighbor before discretisation

          if (i1e  == nsegmax .or. tabjonc(i1e)) then

              P3Iini(1:3) = A_coor(1:3,Iini) + AB_length(Iini) * bveclin(1:3,IVLI)  ! end of I
              compteurP3=IUN
              ROOP3Iini(1:3,compteurP3) = A_coor(1:3,Iini) - P2Iini(1:3)
              listP3Iini(compteurP3) = Iini

          else

              ! begining of effective I1e = end of I

              Oi(1:3) = A_coor(1:3,Iini) + AB_length(Iini) * bveclin(1:3,IVLI)

              ! starting end of I1e

              Tr(:) =   AB_length(i1e) * bveclin(1:3,VLseg(i1e)) / IDEUX
              Ei(1:3) = Oi(1:3) + Tr(:) * IDEUX

              compteurP3 = IUN

              if (Rlocal2 < zero) then

                P3Iini(1:3) = Oi(1:3) + Tr(:) ! P3 is the center of the vnn segment
                ROOP3Iini(1:3,compteurP3) = (Oi(1:3) + Tr(1:3)) - P2Iini(1:3)
                listP3Iini(compteurP3) = i1e

              else

                ROOP3Iini(1:3,compteurP3)= (Oi(1:3) + Tr(1:3)) - P2Iini(1:3) ! The vector connecting P2Iini and the center of i1e

                ! Do we need to look for more distant segment i1e to regularize the line tension curvature ?

                if (SUM(real(ROOP3Iini(:,compteurP3),DP)*real(ROOP3Iini(:,compteurP3),DP)) < Rlocal2) then

                  !The half length of i1e segment
                  TrAB(:) = Tr(1:3)

                  ! Ra and Rb are vectors between P2Iini and the extremities of i1o
                  RA(:)=real(ROOP3Iini(:,CompteurP3)+TrAB(:),DP)
                  RB(:)=real(ROOP3Iini(:,CompteurP3)-TrAB(:),DP)

                  vnn = i1e
                  listP3Iini(compteurP3) = i1e

                  ! exact calculation of the line tension impose that the local radius of curvature is calculated
                  ! with points P1Iini and P3 not very close from P2Iini. A minimum distance of Rlocal is imposed!

                  do while (seg(vnn)%vnne/= nsegmax         .and.     &
                            .not.tabjonc(vnn)               .and.     &
                            tabgd(vnn) < iun                .and.     &
                            SUM(RA*RA) < Rlocal2            .and.     &
                            SUM(RB*RB) < Rlocal2)

                    if(compteurP3 > idix) exit
                    ! The next i1e segment along the dislocation line is tested
                    vnn = seg(vnn)%vnne
                    compteurP3 = compteurP3 + IUN

                    ! The vector connecting P2Iini and the center of previous vnn
                    Tr(1:3) = TrAB(1:3)

                    !The half length of the new vnn segment
                    TrAB(1:3) = (AB_length(vnn)*bveclin(1:3,VLseg(vnn))) / ideux

                    ! The vector connecting P2Iini and the center of vnn
                    ROOP3Iini(1:3,compteurP3) = ROOP3Iini(1:3,compteurP3-IUN) + (Tr(1:3) + TrAB(1:3))
                    listP3Iini(compteurP3)=vnn

                    ! Ra et Rb distances entre le centre de i et les extremites de vnn
                    RA(:) = real(ROOP3Iini(:,compteurP3) - TrAB(:),DP)
                    RB(:) = real(ROOP3Iini(:,compteurP3) + TrAB(:),DP)

                  enddo

                  P3Iini(1:3) = P2Iini(1:3) + ROOP3Iini(1:3,compteurP3) ! P3 is the center of the vnn segment

                  ! A special case for segments touching the surface. Since we don t do the obstacle research
                  ! segments overlap may happen and induce wrong ROO calculation
                  if (SUM(ROOP3Iini(1:3,compteurP3)*ROOP3Iini(1:3,compteurP3)) == 0 &
                     .and. seg(Iini)%surface > 0) then
                     P3Iini(1:3) = A_coor(1:3,Iini) + AB_length(Iini)* bveclin(1:3,IVLI)
                     compteurP3=IUN
                     ROOP3Iini(1:3,compteurP3) = P3Iini(1:3) - P2Iini(1:3)
                     listP3Iini(compteurP3) = Iini
                  endif

                endif

              endif

              ! Now we look for the problem of long straight lines cut in many segments
              ! straight lines section must be treated as one long segment, whatever rlocal
              if (compteurP3 == IUN) then
                 ! The initial vnn is redefined
                 vnn = i1e
                 listP3Iini(compteurP3) = vnn

                 do while(seg(vnn)%vnne /= nsegmax .and. VLseg(seg(vnn)%vnne) == vlseg(i1e))
                      if(compteurP3 > idix) exit
                      vnn=seg(vnn)%vnne
                      compteurP3 = compteurP3 + IUN

                      !The half length of the new vnn segment

                      Ei(1:3) = Ei(1:3) + AB_length(vnn) * bveclin(1:3,VLseg(vnn))
                      !The center of the new vnn segment
                      ROOP3Iini(1:3,compteurP3) = (Oi(1:3)+Ei(1:3))/IDEUX -P2Iini(1:3)
                      listP3Iini(compteurP3) = vnn
                 enddo

                 P3Iini(1:3) = P2Iini(1:3) + ROOP3Iini(1:3,compteurP3)

              endif

          endif ! End of the P3 tests

          !P3P2 and P1P2 are the two chords used to compute the radius of curvature.
          P3P2(1:3)=real(P3Iini(1:3)-P2Iini(1:3),DP)
          P1P2(1:3)=real(P1Iini(1:3)-P2Iini(1:3),DP)

          ! length of the two chordes
          normP3P2=dsqrt(P3P2(1)*P3P2(1)+P3P2(2)*P3P2(2)+P3P2(3)*P3P2(3))
          normP1P2=dsqrt(P1P2(1)*P1P2(1)+P1P2(2)*P1P2(2)+P1P2(3)*P1P2(3))

          ! angle between the two cordes
          if (normP1P2 /= zero .and. normP3P2 /= zero) then

            cosP1P2vP3P2=(P3P2(1)*P1P2(1)+P3P2(2)*P1P2(2)+P3P2(3)*P1P2(3))/(normP3P2*normP1P2)

          else

            cosP1P2vP3P2 = un
            approxLT(Iini)=.true.
            !      print*, approxLT(iini), '09elasti line 2499'
            !                  write(379,*)  'P1', P1Iini
            !                  write(379,*)  'P2', P2Iini
            !                  write(379,*)  'P3', P3Iini
            !                  call seginfo(iini, 'linetension')
            !                  read(*,*)

          endif

          !two cases:
          !1) the angle between P1P2 and P2P3 is larger then pi/2 cosP1P2vP3P2 < 0:
          !nothing to do, we can always define a circumference passing from P1 P2 and P3
          !and the total length of the arc is the P1P2 arc + P3P2 arc
          !2) the angle between P1P2 and P2P3 is smaller then pi/2  cosP1P2vP3P2 > 0:
          ! We change the length of P2P1 or P3P2 in order to define a circumference passing
          ! from P1 P2 and P3 in which the total length of the arc is the sum of P1P2 arc and the P3P2 arc
          !

          if (cosP1P2vP3P2 > zero .and. normP1P2 /= ZERO .and. normP3P2 /= ZERO) then

            invcosP1P2vP3P2=UN/cosP1P2vP3P2

            ! the P3P2 arc length is reduced
            if ((normP1P2*invcosP1P2vP3P2 - normP3P2) < 1.e-12) then

              do ii = 1,CompteurP3 !P3 point already computed are exploited

                P3temp=P2Iini+ROOP3Iini(:,ii)
                P3P2(1:3)=real(P3temp(1:3)-P2Iini(1:3),DP)
                normP3P2=dsqrt(P3P2(1)*P3P2(1)+P3P2(2)*P3P2(2)+P3P2(3)*P3P2(3))

                if (normP1P2 /= zero .and. normP3P2 /= zero) then
                 cosP1P2vP3P2=(P3P2(1)*P1P2(1)+P3P2(2)*P1P2(2)+P3P2(3)*P1P2(3))/(normP3P2*normP1P2)
                else
                 cosP1P2vP3P2 = un
                endif

                if (cosP1P2vP3P2 <= zero) then
                  P3Iini=p3temp
                  cycle
                endif

                invcosP1P2vP3P2 = un / cosP1P2vP3P2

                if ((normP1P2*invcosP1P2vP3P2 - normP3P2) > 1e-12 ) then
                  P3Iini=P3temp
                  cycle
                else
                  compteurP3 = ii-1
                  exit
                endif

              enddo

              if(ii == iun) then ! the first neighbour is reached

                compteurP3=iun
                halflength=AB_length(i1e)

                do jj = 1,halflength ! we change the definition of P3

                  P3Iini=P2Iini+ROOP3Iini(:,ii) - jj * bveclin(1:3,VLseg(i1e))/IDEUX
                  P3P2(1:3)=real(P3Iini(1:3)-P2Iini(1:3),DP)
                  normP3P2=dsqrt(P3P2(1)*P3P2(1)+P3P2(2)*P3P2(2)+P3P2(3)*P3P2(3))

                  if (normP1P2 /= zero .and. normP3P2 /= zero) then
                    cosP1P2vP3P2=(P3P2(1)*P1P2(1)+P3P2(2)*P1P2(2)+P3P2(3)*P1P2(3))/(normP3P2*normP1P2)
                  else
                    cosP1P2vP3P2 = un
                  endif

                  if (cosP1P2vP3P2 <=zero) exit

                  invcosP1P2vP3P2=UN/cosP1P2vP3P2

                  if (normP1P2*invcosP1P2vP3P2 - normP3P2 > 1e-12 ) exit

                  if (jj==halflength) then
                    approxLT(Iini)=.true.
                    P3Iini=P2Iini+ROOP3Iini(:,ii)
                    if(kkdebug) write(379,*) "problem P3", i, kk , approxLT(Iini)
                  endif

                enddo

              endif

            ! the P1P2 arc length is reduced
            elseif ((normP1P2 - normP3P2*invcosP1P2vP3P2) > -1.e-12) then

              do ii = 1,CompteurP1 !P1 point already computed are exploited

                P1temp=P2Iini+ROOP1Iini(:,ii)
                P1P2(1:3)=real(P1temp(1:3)-P2Iini(1:3),DP)
                normP1P2=dsqrt(P1P2(1)*P1P2(1)+P1P2(2)*P1P2(2)+P1P2(3)*P1P2(3))

                if (normP1P2 /= zero .and. normP3P2 /= zero) then

                  cosP1P2vP3P2=(P3P2(1)*P1P2(1)+P3P2(2)*P1P2(2)+P3P2(3)*P1P2(3))/(normP3P2*normP1P2)

                else

                  cosP1P2vP3P2 = un

                endif

                if (cosP1P2vP3P2 <= zero) then
                  P1Iini=p1temp
                  cycle
                endif

                invcosP1P2vP3P2=UN/cosP1P2vP3P2

                if ((normP1P2 - normP3P2*invcosP1P2vP3P2) < -1.e-12) then
                  P1Iini=P1temp
                  cycle
                else
                  compteurP1 = ii-1
                  exit
                endif

              enddo

              if(ii == iun) then ! the first neighbour is reached

                compteurP1 = IUN
                halflength = AB_length(i1o)

                do jj = 1,halflength! we change the definition of P3

                  P1Iini=P2Iini+ROOP1Iini(:,ii) + jj * bveclin(1:3,VLseg(i1o))/IDEUX
                  P1P2(1:3)=real(P1Iini(1:3)-P2Iini(1:3),DP)
                  normP1P2=dsqrt(P1P2(1)*P1P2(1)+P1P2(2)*P1P2(2)+P1P2(3)*P1P2(3))

                  if (normP1P2 /= zero .and. normP3P2 /= zero) then

                   cosP1P2vP3P2=(P3P2(1)*P1P2(1)+P3P2(2)*P1P2(2)+P3P2(3)*P1P2(3))/(normP3P2*normP1P2)

                  else

                   cosP1P2vP3P2 = un

                  endif

                  if (cosP1P2vP3P2 <=zero) exit
                  invcosP1P2vP3P2=UN/cosP1P2vP3P2

                  if (normP1P2 - normP3P2*invcosP1P2vP3P2 < -1.e-12) exit

                  if (jj==halflength) then
                   approxLT(Iini)=.true.
                   P1Iini=P2Iini+ROOP1Iini(:,ii)
                   if(kkdebug) write(379,*) "problem P1",i, kk , approxLT(Iini)
                  endif

                enddo

              endif

            endif

          endif

          ! The three points used to calculate the locale radius of curvature are defined
          P1(1:3,Iini) = P1Iini(1:3)
          P2(1:3,Iini) = P2Iini(1:3)
          P3(1:3,Iini) = P3Iini(1:3)

          !!!!!!!!!!!!!!!!!!!!!

          DCFJ_norm = INT(DCFJ / Normlin(IVLI))
          ! determination du centre de calcul pour I
          !*** Calcul de la contrainte au centre des segments "i"
          !*** Attention au test logique ou il fait n'importe quoi et on s'en
          !*** rend difficlement compte !!!
          ! Ghiath / pour simlifier la lecture : la condtion sous laquelle on deplace le
          ! point de calcul de force a une distance DCFJ du noeud est evaluee separement

          condition = ((IDEUX * DCFJ_norm < AB_length(iini)) .and. &
                       (tabjonc(i1o) .or. tabjonc(i1e) .or. (tabgd(i1o) + tabgd(i1e)) > izero))

          if (condition) then
            ! cas ou le voisin en O est gd ou jonc : calcul de la force a DCFJ du noeud
            ! print *, " Li =", IVL(I),"    DCFJeff =", DCFJeff
            if (tabjonc(i1o) .or. tabgd(i1o) > izero) then
               ! = calculer la force a DCFJ pas de la jonction si norme > x et si
               ! un des deux voisin est un bras de jonction
               R(1:3) = A_COOR(1:3,Iini)+ BVecLin(1:3,IVLi)* DCFJ_norm !
            else
               R(1:3) = A_COOR(1:3,I1E) - BVecLin(1:3,IVLi)* DCFJ_norm !dcalj
               ! decoupeavannihil = .false.     print *, " A_COOR(1:3,I1E)A",A_COOR(1:3,I1E)
            endif
          else
            ! Cas normal : calcul classique de la force au centre de I
            R(1:3) = RO(:,Iini)
          endif

          Tr(:) =  Int(HalfModur(1:3)  - R(1:3),DPI) ! vecteur de translation applique sur tous les segments
          ! pour les CPL

          !******************************************
          !* ballayage sur tous les segemnts voisins
          !******************************************
          SIGINTERN(:) = zero
#ifdef MDC
          SIGINTNS(:) = zero
#endif

           ! ============================================================================================
           ! ====  ici on retient la boite et les 26 boite 1ere voisins  ================================
           ! ============================================================================================

           ! pour optimiser la recherche des boites a CP
           ! il vaut mieux une triple boucle explicite, plutot qu un
           ! tableau indexe et parcourant toute les boites
#ifdef MDC
           Nsegtester = izero
           IX=B1D_3D(IB,1)-IUN
           IY=B1D_3D(IB,2)-IUN
           IZ=B1D_3D(IB,3)-IUN


           BOITE_I1: do I1 = -IDEBX , IDEBX
               bx = modulo(IX+I1,NBoitesX)+IUN
               BOITE_J1: do J1 = -IDEBY , IDEBY
                   by = modulo(IY+J1,NBoitesY)+IUN
                   BOITE_K1: do K1 =  -IDEBZ , IDEBZ
                       bz = modulo(IZ+K1,NBoitesZ)+IUN


#else
           IDEB=IUN
           IX=B1D_3D(IB,1)-IUN
           IY=B1D_3D(IB,2)-IUN
           IZ=B1D_3D(IB,3)-IUN
           if (NBOITES==1) IDEB=0  ! Cas ou l on a qu une boite
           Nsegtester = izero

           BOITE_I1: do I1 = -IDEB , IDEB
               bx = modulo(IX+I1,NBoitesX)+IUN
               BOITE_J1: do J1 = -IDEB , IDEB
                   by = modulo(IY+J1,NBoitesY)+IUN
                   BOITE_K1: do K1 =  -IDEB , IDEB
                       bz = modulo(IZ+K1,NBoitesZ)+IUN
#endif
                       JB= BX + ((BY-1)+ (BZ-1)*NBoitesY)*NBoitesX

                       ! boucle sur les segments dans chaque boites
                       BI : do IndiceJ = 1, NSEGBoite(JB)
                           i = IndexBoite(JB)%ListeSeg(indiceJ)

                           ! On laisse de cote certains segments qui sont caches
                           if (ab_length(i) /= izero) then
                              if (I==I1O .and. .not. tabjonc(I) .and. tabgd(I) < iun)     CYCLE BI
                              if (I==I1E .and. .not. tabjonc(I) .and. tabgd(I) < iun)     CYCLE BI
                              if (i == iini)                                              CYCLE BI
                              ! new ! long segments are added further
                              if (longseg(i))                                             CYCLE BI
                              if (NbCVXDom > iun .and. Domconnec(segDom(i),segDom(iini))) CYCLE BI
                              if(out(i))               CYCLE BI
                              Nsegtester = nsegtester + IUN
                              listesegtester(nsegtester) = I
                           endif
                       enddo BI
                   enddo BOITE_K1
               enddo BOITE_J1
           enddo BOITE_I1

           ! we add long segments to the list
           DO I1 = 1, Nb_groseg
              I = GROSEG(I1)
              if (I==I1O .and. .not.tabjonc(I) .and. tabgd(I) < iun) CYCLE
              if (I==I1E .and. .not.tabjonc(I) .and. tabgd(I) < iun) CYCLE
              if (i == iini)                                         CYCLE
              Nsegtester = nsegtester + IUN
              listesegtester(nsegtester) = I
           enddo



           do IndiceJ = 1, Nsegtester

            I = listesegtester(IndiceJ)                           ! Index of the segment to test


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Segments included in the regularized local line tension must be excluded
            ! i and iini are very close and may be along the same line
            ! We must check that the interaction between both segment is not already
            ! accounted for in the regularized local line tension

            if(any(listP1Iini(1:compteurP1)-I .eq. IZERO)) cycle  ! segments used in the LT regularisation are eliminated (P1 side)

            if(any(listP3Iini(1:compteurP3)-I .eq. IZERO)) cycle  ! segments used in the LT regularisation are eliminated (P3 side)

            !if (out(i)) cycle    ! Such segment does not need to be considered since it is to eliminate

              ! The vector connecting the I segment center and the point coordinate (with modulo translation)
              ! where the stress calculation is needed :
              ! ROOmod(:)= HalfModur(:) - modulo((ro(:,i) + Tr(:)),ModuR(:))
              aX = real(ro(1,i) + Tr(1),DP)
              aY = real(ro(2,i) + Tr(2),DP)
              aZ = real(ro(3,i) + Tr(3),DP)
              ROOmod(1)= HalfModurX - (ax - floor(ax*invmodurx)*modurx) !modulo(ax,ModuRX)
              ROOmod(2)= HalfModurY - (ay - floor(ay*invmodury)*modury) !modulo(ay,ModuRY)
              ROOmod(3)= HalfModurZ - (az - floor(az*invmodurz)*modurz) !modulo(az,ModuRZ)
              IVLI = VLseg(I)
              Tbis(:) = bveclin(1:3,IVLI)
              realnorme = real(AB_LENGTH(i),DP)
#ifdef MDC
              dist2=ROOmod(1)*ROOmod(1)+ROOmod(2)*ROOmod(2)+ROOmod(3)*ROOmod(3)
              ExtR2=SR_cutoff + realnorme*half
              ExtR2=ExtR2*ExtR2
              !if the distance between segment centers <= than the cut-off radius we already know we must compute the short-range correction.
              !On the other hand, if the distance > we need a test to verify if the center of the segment (R) in which we want to compute the force is
              !inside the homogeneization region of the segment i
              if (dist2 > ExtR2) cycle
              if (dist2 > SR_cutoff2) then
                buf(4)    =ROOmod(1)*tbis(1)+ROOmod(2)*tbis(2)+ROOmod(3)*tbis(3) ! v.t
                buf(5)    =buf(4)/(tbis(1)*tbis(1)+tbis(2)*tbis(2)+tbis(3)*tbis(3))! (v.t)/(t.t)
                buf(6:8)  =buf(5)*tbis(1:3) ! vt=(v.t)/(t.t)*t
                buf(9:11) =ROOmod(1:3)-buf(6:8) ! vn=v-vt
                dn2=buf(9)*buf(9)+buf(10)*buf(10)+buf(11)*buf(11)
                dt2=buf(6)*buf(6)+buf(7)*buf(7)+buf(8)*buf(8)
                dt=dsqrt(dt2)
                ied=((dt-realnorme*half)**2+dn2 < SR_cutoff2) ! test if gauss point is inside a end sphere
                if (.not.(ied .or. (dn2 < SR_cutoff2 .and. dt<realnorme*half))) cycle! R is outside the homogeneizetion area of segment I
              endif
#endif

              TrABmod(:) = realnorme *Tbis(1:3) * half


              ! Ra and Rb are vectors pointing the origin and the end of segment I
              RA(:)=ROOmod(:)+TrABmod(:)
              RB(:)=ROOmod(:)-TrABmod(:)


              !! Pour eliminer les CP entre segments d un meme systeme
              !! e.g. le systeme 4 dans les manipes modeles
              !! Attention RA et RB ne sont pas normalises
              !if ((syseg(seg(iini)%veclin) == 4) .and.  &  ! avec la self contrainte
              !    (syseg(seg(i)%veclin) == 4) .and.     &  ! avec la self contrainte
              !    (abs(dot_product(RA(:),real(bvecnor(:,seg(i)%veclin),DP))) .gt. 0.1)) cycle BI

              ! Les caracteristiques du segment utiles
              Tbis(:)=T(:,I)
              Bbis(:)=B(:,I)
              BPVTbis(:)=BPVT(:,I)

              R_ps_t=RA(1)*Tbis(1)+RA(2)*Tbis(2)+RA(3)*Tbis(3)
              s1=-R_ps_t
              cai_d=RA(:)-R_ps_t*Tbis(:)
              d2=cai_d(1)*cai_d(1)+cai_d(2)*cai_d(2)+cai_d(3)*cai_d(3) ! distance au segment
              R_ps_t=RB(1)*Tbis(1)+RB(2)*Tbis(2)+RB(3)*Tbis(3)
              s2=-R_ps_t

              t_pv_b=-BPVTbis(:)

              !!d_pv_b=prodvec(cai_d,B)
              d_pv_b(1)=cai_d(2)*Bbis(3)-cai_d(3)*Bbis(2)
              d_pv_b(2)=cai_d(3)*Bbis(1)-cai_d(1)*Bbis(3)
              d_pv_b(3)=cai_d(1)*Bbis(2)-cai_d(2)*Bbis(1)

              d_pv_b_ps_t=d_pv_b(1)*Tbis(1)+d_pv_b(2)*Tbis(2)+d_pv_b(3)*Tbis(3)


              !!t_pt_d=SymTens6Vec(localT,cai_d)
              t_pt_d(1) = (Tbis(1)*cai_d(1))*2.d0
              t_pt_d(2) = (Tbis(2)*cai_d(2))*2.d0
              t_pt_d(3) = (Tbis(3)*cai_d(3))*2.d0
              t_pt_d(4) = Tbis(2)*cai_d(3)+cai_d(2)*Tbis(3)
              t_pt_d(5) = Tbis(1)*cai_d(3)+cai_d(1)*Tbis(3)
              t_pt_d(6) = Tbis(1)*cai_d(2)+cai_d(1)*Tbis(2)


              !!t_pt_t_pv_b=SymTens6Vec(T,t_pv_b)
              t_pt_t_pv_b(1) = (Tbis(1)*t_pv_b(1))*2.d0
              t_pt_t_pv_b(2) = (Tbis(2)*t_pv_b(2))*2.d0
              t_pt_t_pv_b(3) = (Tbis(3)*t_pv_b(3))*2.d0
              t_pt_t_pv_b(4) = Tbis(2)*t_pv_b(3)+t_pv_b(2)*Tbis(3)
              t_pt_t_pv_b(5) = Tbis(1)*t_pv_b(3)+t_pv_b(1)*Tbis(3)
              t_pt_t_pv_b(6) = Tbis(1)*t_pv_b(2)+t_pv_b(1)*Tbis(2)
              !!d_pt_t_pv_b=SymTens6Vec(cai_d,t_pv_b)
              d_pt_t_pv_b(1) = (cai_d(1)*t_pv_b(1))*2.d0
              d_pt_t_pv_b(2) = (cai_d(2)*t_pv_b(2))*2.d0
              d_pt_t_pv_b(3) = (cai_d(3)*t_pv_b(3))*2.d0
              d_pt_t_pv_b(4) = cai_d(2)*t_pv_b(3)+t_pv_b(2)*cai_d(3)
              d_pt_t_pv_b(5) = cai_d(1)*t_pv_b(3)+t_pv_b(1)*cai_d(3)
              d_pt_t_pv_b(6) = cai_d(1)*t_pv_b(2)+t_pv_b(1)*cai_d(2)

              !!t_pt_d_pv_b=SymTens6Vec(localT,d_pv_b)
              t_pt_d_pv_b(1) = (Tbis(1)*d_pv_b(1))*2.d0
              t_pt_d_pv_b(2) = (Tbis(2)*d_pv_b(2))*2.d0
              t_pt_d_pv_b(3) = (Tbis(3)*d_pv_b(3))*2.d0
              t_pt_d_pv_b(4) = Tbis(2)*d_pv_b(3)+d_pv_b(2)*Tbis(3)
              t_pt_d_pv_b(5) = Tbis(1)*d_pv_b(3)+d_pv_b(1)*Tbis(3)
              t_pt_d_pv_b(6) = Tbis(1)*d_pv_b(2)+d_pv_b(1)*Tbis(2)


              !!I_03=m4pn*(d_pv_b_ps_t*TENS6ID+d_pt_t_pv_b)-m4p*t_pt_d_pv_b
              I_03(1)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(1))-0.25d0*t_pt_d_pv_b(1)
              I_03(2)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(2))-0.25d0*t_pt_d_pv_b(2)
              I_03(3)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(3))-0.25d0*t_pt_d_pv_b(3)
              I_03(4)=m4pn*d_pt_t_pv_b(4)-0.25d0*t_pt_d_pv_b(4)
              I_03(5)=m4pn*d_pt_t_pv_b(5)-0.25d0*t_pt_d_pv_b(5)
              I_03(6)=m4pn*d_pt_t_pv_b(6)-0.25d0*t_pt_d_pv_b(6)


              !a2_d2 = a2+d2;
              !temp=1./a2_d2;
              !a2_p_d2=halfthickness2+d2
              a2_p_d2=bspread2+d2
              inv_a2_p_d2=1.d0/a2_p_d2

              !Ra    = sqrt( [a2_d2 a2_d2] + s.*s);
              !Rainv=1./Ra;
              !Ra3inv=Rainv.*Rainv.*Rainv;
              Ra1=a2_p_d2+s1*s1
              Ra1_inv=1.d0/sqrt(Ra1)
              Ra1_3inv=Ra1_inv/Ra1 ! 1/ + cher que 2* ? OUI
              s1Ra3inv=s1*Ra1_3inv

              !Ra1_3inv=Ra1_inv*Ra1_inv*Ra1_inv
              Ra2=a2_p_d2+s2*s2
              Ra2_inv=1.d0/sqrt(Ra2)
              Ra2_3inv=Ra2_inv/Ra2 ! 1/ + cher que 2* ? ENFIN PAREIL (dans la barre d'erreur)
              !Ra2_3inv=Ra2_inv*Ra2_inv*Ra2_inv
              s2Ra3inv=s2*Ra2_3inv

              s1_03=s1*Ra1_inv*inv_a2_p_d2
              s2_03=s2*Ra2_inv*inv_a2_p_d2
              s_03=s2_03-s1_03



#ifdef MDC
              if (zcall) then
                 a2_p_d2_MDC=halfthickness2+d2
                 inv_a2_p_d2_MDC=1.d0/a2_p_d2_MDC

                 !Ra    = sqrt( [a2_d2 a2_d2] + s.*s);
                 !Rainv=1./Ra;
                 !Ra3inv=Rainv.*Rainv.*Rainv;
                 Ra1_MDC=a2_p_d2_MDC+s1*s1
                 Ra1_inv_MDC=1.d0/sqrt(Ra1_MDC)
                 Ra1_3inv_MDC=Ra1_inv_MDC/Ra1_MDC ! 1/ + cher que 2* ? OUI
                 s1Ra3inv_MDC=s1*Ra1_3inv_MDC

                 !Ra1_3inv=Ra1_inv*Ra1_inv*Ra1_inv
                 Ra2_MDC=a2_p_d2_MDC+s2*s2
                 Ra2_inv_MDC=1.d0/sqrt(Ra2_MDC)
                 Ra2_3inv_MDC=Ra2_inv_MDC/Ra2_MDC ! 1/ + cher que 2* ? ENFIN PAREIL (dans la barre d'erreur)
                 !Ra2_3inv=Ra2_inv*Ra2_inv*Ra2_inv
                 s2Ra3inv_MDC=s2*Ra2_3inv_MDC

                 s1_03_MDC=s1*Ra1_inv*inv_a2_p_d2_MDC
                 s2_03_MDC=s2*Ra2_inv*inv_a2_p_d2_MDC
                 s_03_MDC=s2_03_MDC-s1_03_MDC
              endif

#endif

              ! Si les segments sont vis pures on reduit le calcul
              if(TYSEG(IVLI) == 1) then

                 SIGINTERN(1) = SIGINTERN(1) + I_03(1) *s_03
                 SIGINTERN(2) = SIGINTERN(2) + I_03(2) *s_03
                 SIGINTERN(3) = SIGINTERN(3) + I_03(3) *s_03
                 SIGINTERN(4) = SIGINTERN(4) + I_03(4) *s_03
                 SIGINTERN(5) = SIGINTERN(5) + I_03(5) *s_03
                 SIGINTERN(6) = SIGINTERN(6) + I_03(6) *s_03

#ifdef MDC
                 if (zcall) then
                    SIGINTNS(1)=SIGINTNS(1) + I_03(1)*s_03_MDC
                    SIGINTNS(2)=SIGINTNS(2) + I_03(2)*s_03_MDC
                    SIGINTNS(3)=SIGINTNS(3) + I_03(3)*s_03_MDC
                    SIGINTNS(4)=SIGINTNS(4) + I_03(4)*s_03_MDC
                    SIGINTNS(5)=SIGINTNS(5) + I_03(5)*s_03_MDC
                    SIGINTNS(6)=SIGINTNS(6) + I_03(6)*s_03_MDC
                 endif
#endif

                 cycle
              endif

              !!t_pt_t=Tens6SelfVec(T) ! tabulable
              t_pt_t(1) = Tbis(1)*Tbis(1)
              t_pt_t(2) = Tbis(2)*Tbis(2)
              t_pt_t(3) = Tbis(3)*Tbis(3)
              t_pt_t(4) = Tbis(2)*Tbis(3)
              t_pt_t(5) = Tbis(1)*Tbis(3)
              t_pt_t(6) = Tbis(1)*Tbis(2)



              !!d_pt_d=Tens6SelfVec(cai_d)
              d_pt_d(1) = cai_d(1)*cai_d(1)
              d_pt_d(2) = cai_d(2)*cai_d(2)
              d_pt_d(3) = cai_d(3)*cai_d(3)
              d_pt_d(4) = cai_d(2)*cai_d(3)
              d_pt_d(5) = cai_d(1)*cai_d(3)
              d_pt_d(6) = cai_d(1)*cai_d(2)

              commun=m4pn*d_pv_b_ps_t

              I_13=-mn4pn*t_pt_t_pv_b(:)

              I_25=commun*t_pt_t


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




              !!I_15=a2m8p*t_pt_t_pv_b-commun*t_pt_d
              I_15(1)=a2m8p*t_pt_t_pv_b(1)-commun*t_pt_d(1)
              I_15(2)=a2m8p*t_pt_t_pv_b(2)-commun*t_pt_d(2)
              I_15(3)=a2m8p*t_pt_t_pv_b(3)-commun*t_pt_d(3)
              I_15(4)=a2m8p*t_pt_t_pv_b(4)-commun*t_pt_d(4)
              I_15(5)=a2m8p*t_pt_t_pv_b(5)-commun*t_pt_d(5)
              I_15(6)=a2m8p*t_pt_t_pv_b(6)-commun*t_pt_d(6)



              !!I_05=commun*(core2*TENS6ID+d_pt_d)-a2m8p*t_pt_d_pv_b
              I_05(1)=commun*(bspread2+d_pt_d(1))-a2m8p*t_pt_d_pv_b(1)
              I_05(2)=commun*(bspread2+d_pt_d(2))-a2m8p*t_pt_d_pv_b(2)
              I_05(3)=commun*(bspread2+d_pt_d(3))-a2m8p*t_pt_d_pv_b(3)
              I_05(4)=commun*d_pt_d(4)-a2m8p*t_pt_d_pv_b(4)
              I_05(5)=commun*d_pt_d(5)-a2m8p*t_pt_d_pv_b(5)
              I_05(6)=commun*d_pt_d(6)-a2m8p*t_pt_d_pv_b(6)


              SIGINTERN(1) = SIGINTERN(1) + (I_03(1)*s_03+I_13(1)*s_13+I_05(1)*s_05+I_15(1)*s_15+I_25(1)*s_25)
              SIGINTERN(2) = SIGINTERN(2) + (I_03(2)*s_03+I_13(2)*s_13+I_05(2)*s_05+I_15(2)*s_15+I_25(2)*s_25)
              SIGINTERN(3) = SIGINTERN(3) + (I_03(3)*s_03+I_13(3)*s_13+I_05(3)*s_05+I_15(3)*s_15+I_25(3)*s_25)
              SIGINTERN(4) = SIGINTERN(4) + (I_03(4)*s_03+I_13(4)*s_13+I_05(4)*s_05+I_15(4)*s_15+I_25(4)*s_25)
              SIGINTERN(5) = SIGINTERN(5) + (I_03(5)*s_03+I_13(5)*s_13+I_05(5)*s_05+I_15(5)*s_15+I_25(5)*s_25)
              SIGINTERN(6) = SIGINTERN(6) + (I_03(6)*s_03+I_13(6)*s_13+I_05(6)*s_05+I_15(6)*s_15+I_25(6)*s_25)

#ifdef MDC
              if (zcall) then
                 s1_13_MDC=-Ra1_inv_MDC
                 s2_13_MDC=-Ra2_inv_MDC
                 s_13_MDC=s2_13_MDC-s1_13_MDC

                 s1_05_MDC=(2.d0*s1_03_MDC+s1Ra3inv_MDC)*inv_a2_p_d2_MDC
                 s2_05_MDC=(2.d0*s2_03_MDC+s2Ra3inv_MDC)*inv_a2_p_d2_MDC
                 s_05_MDC=s2_05_MDC-s1_05_MDC

                 s1_15_MDC=-Ra1_3inv_MDC
                 s2_15_MDC=-Ra2_3inv_MDC
                 s_15_MDC=s2_15_MDC-s1_15_MDC

                 s1_25_MDC=s1_03_MDC-s1Ra3inv_MDC
                 s2_25_MDC=s2_03_MDC-s2Ra3inv_MDC
                 s_25_MDC=s2_25_MDC-s1_25_MDC


                 !!I_15=a2m8p*t_pt_t_pv_b-commun*t_pt_d
                 I_15_MDC(1)=a2m8p_MDC*t_pt_t_pv_b(1)-commun*t_pt_d(1)
                 I_15_MDC(2)=a2m8p_MDC*t_pt_t_pv_b(2)-commun*t_pt_d(2)
                 I_15_MDC(3)=a2m8p_MDC*t_pt_t_pv_b(3)-commun*t_pt_d(3)
                 I_15_MDC(4)=a2m8p_MDC*t_pt_t_pv_b(4)-commun*t_pt_d(4)
                 I_15_MDC(5)=a2m8p_MDC*t_pt_t_pv_b(5)-commun*t_pt_d(5)
                 I_15_MDC(6)=a2m8p_MDC*t_pt_t_pv_b(6)-commun*t_pt_d(6)


                 !!I_05=commun*(core2*TENS6ID+d_pt_d)-a2m8p*t_pt_d_pv_b
                 I_05_MDC(1)=commun*(halfthickness2+d_pt_d(1))-a2m8p_MDC*t_pt_d_pv_b(1)
                 I_05_MDC(2)=commun*(halfthickness2+d_pt_d(2))-a2m8p_MDC*t_pt_d_pv_b(2)
                 I_05_MDC(3)=commun*(halfthickness2+d_pt_d(3))-a2m8p_MDC*t_pt_d_pv_b(3)
                 I_05_MDC(4)=commun*d_pt_d(4)-a2m8p_MDC*t_pt_d_pv_b(4)
                 I_05_MDC(5)=commun*d_pt_d(5)-a2m8p_MDC*t_pt_d_pv_b(5)
                 I_05_MDC(6)=commun*d_pt_d(6)-a2m8p_MDC*t_pt_d_pv_b(6)

                 SIGINTNS(1)=SIGINTNS(1) + &
                             (I_03(1)*s_03_MDC+I_13(1)*s_13_MDC+I_05_MDC(1)*s_05_MDC+I_15_MDC(1)*s_15_MDC+I_25(1)*s_25_MDC)
                 SIGINTNS(2)=SIGINTNS(2) + &
                             (I_03(2)*s_03_MDC+I_13(2)*s_13_MDC+I_05_MDC(2)*s_05_MDC+I_15_MDC(2)*s_15_MDC+I_25(2)*s_25_MDC)
                 SIGINTNS(3)=SIGINTNS(3) + &
                             (I_03(3)*s_03_MDC+I_13(3)*s_13_MDC+I_05_MDC(3)*s_05_MDC+I_15_MDC(3)*s_15_MDC+I_25(3)*s_25_MDC)
                 SIGINTNS(4)=SIGINTNS(4) + &
                             (I_03(4)*s_03_MDC+I_13(4)*s_13_MDC+I_05_MDC(4)*s_05_MDC+I_15_MDC(4)*s_15_MDC+I_25(4)*s_25_MDC)
                 SIGINTNS(5)=SIGINTNS(5) + &
                             (I_03(5)*s_03_MDC+I_13(5)*s_13_MDC+I_05_MDC(5)*s_05_MDC+I_15_MDC(5)*s_15_MDC+I_25(5)*s_25_MDC)
                 SIGINTNS(6)=SIGINTNS(6) + &
                             (I_03(6)*s_03_MDC+I_13(6)*s_13_MDC+I_05_MDC(6)*s_05_MDC+I_15_MDC(6)*s_15_MDC+I_25(6)*s_25_MDC)
              endif
#endif

           enddo

           SIGINTERN(:) = SIGINTERN(:) * bdivpa

           SIG_CP(2,1,Iini)=SIGINTERN(6) ; SIG_CP(3,2,Iini)=SIGINTERN(4) ; SIG_CP(3,1,Iini)=SIGINTERN(5)
           SIG_CP(1,1,Iini)=SIGINTERN(1) ; SIG_CP(2,2,Iini)=SIGINTERN(2) ; SIG_CP(3,3,Iini)=SIGINTERN(3)
           SIG_CP(1,2,Iini)=SIGINTERN(6) ; SIG_CP(2,3,Iini)=SIGINTERN(4) ; SIG_CP(1,3,Iini)=SIGINTERN(5)

#ifdef MDC
           if (zcall) then
              SIGINTNS(:)=SIGINTNS(:)*bdivpa

              SIG_NS(2,1,Iini)=SIGINTNS(6) ; SIG_NS(3,2,Iini)=SIGINTNS(4) ; SIG_NS(3,1,Iini)=SIGINTNS(5)
              SIG_NS(1,1,Iini)=SIGINTNS(1) ; SIG_NS(2,2,Iini)=SIGINTNS(2) ; SIG_NS(3,3,Iini)=SIGINTNS(3)
              SIG_NS(1,2,Iini)=SIGINTNS(6) ; SIG_NS(2,3,Iini)=SIGINTNS(4) ; SIG_NS(1,3,Iini)=SIGINTNS(5)
           endif
#endif
           !if(kkdebug) then
!              write(379,*) "============================================================"
!              write(379,*) SIG_CP(:,:,Iini)
!              write(379,*) Nsegtester
!              write(379,*) "============================================================"
!           endif

        endif

    enddo BIini
enddo LIini

# ifdef PA
! For the Bi-Phased material simulations
if (Nb_phase == 2) then

  inv_homothetie = UN / real(homothetie,DP)

  ! Coordinates of segments is re-defined in a based common to the two phases (free from the homothetie factor)
  do I = 1,Nb_seg
    liste_RO(2:4,I) = nint(real(liste_RO(2:4,I),DP) * inv_homothetie)
  enddo

  ! if (MON_RANG_PAIRIMPAIR==0) print *,'mon_rang',MON_RANG,'liste_RO apres mult:',liste_RO(:,1:Nb_seg)

  ! Stockage valeur Nb_seg
  Nb_seg_stored = Nb_seg

  ! Point2point exchange of the liste_RO between the two phases
  if (Ma_couleur == 0) then
    RANG_DEST = MON_RANG + 1
  else
    RANG_DEST = MON_RANG - 1
  endif

  ! We first send the number of segments to exchange, this is to know the size of list to distribute in the next step
  CALL MPI_SENDRECV_REPLACE(Nb_seg,1,MPI_DOUBLE_PRECISION,RANG_DEST,123,                                                &
                            RANG_DEST,123,MPI_COMM_WORLD,MPI_STATUS_IGNORE,IERR)


  ! allocation of the needed lists
  allocate(liste_RO_Recv(4,Nb_seg))
  allocate(SIGJOINT(6,1:Nb_seg))
  allocate(SIGJOINT_Sum(6,1:Nb_seg))
  allocate(SIGJOINT_Recv(6,1:Nb_seg_stored))

  ! Initialization
  SIGJOINT(:,1:Nb_seg)              = ZERO
  SIGJOINT_Sum(:,1:Nb_seg)          = ZERO
  SIGJOINT_Recv(:,1:Nb_seg_stored)  = ZERO

  ! Transfert of liste_RO
  CALL MPI_SENDRECV(liste_RO,4*Nb_seg_stored,MPI_DOUBLE_PRECISION,RANG_DEST,121,                                        &
                    liste_RO_Recv,4*Nb_seg,MPI_DOUBLE_PRECISION,RANG_DEST,121,MPI_COMM_WORLD,MPI_STATUS_IGNORE,IERR)

  ! We do not need any more this list
  deallocate(liste_RO)


  ! The segment coordinate are re-defined in the proc corresponding phase
  do I = 1,Nb_seg
      liste_RO_Recv(2:4,I) = liste_RO_Recv(2:4,I) * homothetie
  enddo

  ! Computation of the stress on the segment center defined in the list_RO and induced by the segment associated to
  ! the proc phase and at short distance (closed boxes). Result of this stress computation is later transmitted back
  ! to the procs where the seg of list_RO are defined.

  if (Nb_Seg > IZERO) then

    ! Definition of a beginning and an ending box on this list for each proc
    SEGD_JOINT = nint(MON_RANG_PAIRIMPAIR * (real(Nb_seg) / real(TAILLE_PAIRIMPAIR))) + 1
    SEGF_JOINT = nint((MON_RANG_PAIRIMPAIR + 1) * (real(Nb_seg) / real(TAILLE_PAIRIMPAIR)))

    ! loop on the segments for the current proc
    LJIini: do IS = SEGD_JOINT,SEGF_JOINT

      ! The coordinate where the stress must be computed
      IB     = liste_RO_Recv(1,IS)
      R(1:3) = liste_RO_Recv(2:4,IS)

      ! The CLP translation vector is defined
      Tr(:) =  int(HalfModur(1:3),DPI)  - R(1:3)

      !******************************************
      !* ballayage sur tous les segments voisins
      !******************************************
      SIGINTERN(:) = zero

      ! ========================================================================
      ! ====  ici on retient la boite et les 26 boite 1ere voisins  ============
      ! ========================================================================

      ! pour optimiser la recherche des boites a CP
      ! il vaut mieux une triple boucle explicite, plutot qu un
      ! tableau indexe et parcourant toute les boites
      IDEB=IUN
      IX=B1D_3D(IB,1)-IUN
      IY=B1D_3D(IB,2)-IUN
      IZ=B1D_3D(IB,3)-IUN
      if (NBOITES==1) IDEB=0  ! Cas ou l on a qu une boite
      Nsegtester = izero
      listesegtester(:) = izero

      BOITE_J_I1: do I1 = -IDEB , IDEB
        bx = modulo(IX+I1,NBoitesX)+IUN
        BOITE_J_J1: do J1 = -IDEB , IDEB
          by = modulo(IY+J1,NBoitesY)+IUN
          BOITE_J_K1: do K1 =  -IDEB , IDEB
            bz = modulo(IZ+K1,NBoitesZ)+IUN
            JB= BX + ((BY-1) + (BZ-1) * NBoitesY) * NBoitesX
            BI_J : do IndiceJ = 1, NSEGBoite(JB) ! loop on segments in the box

              i = IndexBoite(JB)%ListeSeg(indiceJ)

              ! Some segments can be forget in this computation
              if (ab_length(i) /= izero) then
                if (longseg(i))                                   CYCLE BI_J
                if(i == seg(I)%voiso .and. i == seg(i)%voise)     CYCLE
                Nsegtester = nsegtester + IUN
                listesegtester(nsegtester) = I
              endif

            enddo BI_J
          enddo BOITE_J_K1
        enddo BOITE_J_J1
      enddo BOITE_J_I1

      ! we add long segments to the liste
      do I1 = 1, Nb_groseg
          I = GROSEG(I1)
          Nsegtester = nsegtester + IUN
          listesegtester(nsegtester) = I
      enddo

      do IndiceJ = 1, Nsegtester

        I = listesegtester(IndiceJ)
        if (i == seg(I)%voiso .and. i == seg(i)%voise) CYCLE

        ROOmod(:)= HalfModur(:) - modulo((ro(:,i) + Tr(:)),ModuR(:))
        IVLI = VLseg(I)
        TrABmod(:) = (AB_length(i)*bveclin(1:3,IVLI)) / ideux

        ! Ra et Rb sont l'origine et l'extremite de I avec le modulo
        RA(:)=real(ROOmod(:)+TrABmod(:),DP)
        RB(:)=real(ROOmod(:)-TrABmod(:),DP)

        ! Les caracteristiques du segment utiles
        Tbis(:)=T(:,I)
        Bbis(:)=B(:,I)
        BPVTbis(:)=BPVT(:,I)
        R_ps_t=RA(1)*Tbis(1)+RA(2)*Tbis(2)+RA(3)*Tbis(3)
        s1=-R_ps_t
        cai_d=RA(:)-R_ps_t*Tbis(:)
        d2=cai_d(1)*cai_d(1)+cai_d(2)*cai_d(2)+cai_d(3)*cai_d(3) ! distance au segment
        R_ps_t=RB(1)*Tbis(1)+RB(2)*Tbis(2)+RB(3)*Tbis(3)
        s2=-R_ps_t

        t_pv_b=-BPVTbis(:)

        !!d_pv_b=prodvec(cai_d,B)
        d_pv_b(1)=cai_d(2)*Bbis(3)-cai_d(3)*Bbis(2)
        d_pv_b(2)=cai_d(3)*Bbis(1)-cai_d(1)*Bbis(3)
        d_pv_b(3)=cai_d(1)*Bbis(2)-cai_d(2)*Bbis(1)

        d_pv_b_ps_t=d_pv_b(1)*Tbis(1)+d_pv_b(2)*Tbis(2)+d_pv_b(3)*Tbis(3)


        !!t_pt_d=SymTens6Vec(localT,cai_d)
        t_pt_d(1) = (Tbis(1)*cai_d(1))*2.d0
        t_pt_d(2) = (Tbis(2)*cai_d(2))*2.d0
        t_pt_d(3) = (Tbis(3)*cai_d(3))*2.d0
        t_pt_d(4) = Tbis(2)*cai_d(3)+cai_d(2)*Tbis(3)
        t_pt_d(5) = Tbis(1)*cai_d(3)+cai_d(1)*Tbis(3)
        t_pt_d(6) = Tbis(1)*cai_d(2)+cai_d(1)*Tbis(2)


        !!t_pt_t_pv_b=SymTens6Vec(T,t_pv_b)
        t_pt_t_pv_b(1) = (Tbis(1)*t_pv_b(1))*2.d0
        t_pt_t_pv_b(2) = (Tbis(2)*t_pv_b(2))*2.d0
        t_pt_t_pv_b(3) = (Tbis(3)*t_pv_b(3))*2.d0
        t_pt_t_pv_b(4) = Tbis(2)*t_pv_b(3)+t_pv_b(2)*Tbis(3)
        t_pt_t_pv_b(5) = Tbis(1)*t_pv_b(3)+t_pv_b(1)*Tbis(3)
        t_pt_t_pv_b(6) = Tbis(1)*t_pv_b(2)+t_pv_b(1)*Tbis(2)
        !!d_pt_t_pv_b=SymTens6Vec(cai_d,t_pv_b)
        d_pt_t_pv_b(1) = (cai_d(1)*t_pv_b(1))*2.d0
        d_pt_t_pv_b(2) = (cai_d(2)*t_pv_b(2))*2.d0
        d_pt_t_pv_b(3) = (cai_d(3)*t_pv_b(3))*2.d0
        d_pt_t_pv_b(4) = cai_d(2)*t_pv_b(3)+t_pv_b(2)*cai_d(3)
        d_pt_t_pv_b(5) = cai_d(1)*t_pv_b(3)+t_pv_b(1)*cai_d(3)
        d_pt_t_pv_b(6) = cai_d(1)*t_pv_b(2)+t_pv_b(1)*cai_d(2)

        !!t_pt_d_pv_b=SymTens6Vec(localT,d_pv_b)
        t_pt_d_pv_b(1) = (Tbis(1)*d_pv_b(1))*2.d0
        t_pt_d_pv_b(2) = (Tbis(2)*d_pv_b(2))*2.d0
        t_pt_d_pv_b(3) = (Tbis(3)*d_pv_b(3))*2.d0
        t_pt_d_pv_b(4) = Tbis(2)*d_pv_b(3)+d_pv_b(2)*Tbis(3)
        t_pt_d_pv_b(5) = Tbis(1)*d_pv_b(3)+d_pv_b(1)*Tbis(3)
        t_pt_d_pv_b(6) = Tbis(1)*d_pv_b(2)+d_pv_b(1)*Tbis(2)


        !!I_03=m4pn*(d_pv_b_ps_t*TENS6ID+d_pt_t_pv_b)-m4p*t_pt_d_pv_b
        I_03(1)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(1))-0.25d0*t_pt_d_pv_b(1)
        I_03(2)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(2))-0.25d0*t_pt_d_pv_b(2)
        I_03(3)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(3))-0.25d0*t_pt_d_pv_b(3)
        I_03(4)=m4pn*d_pt_t_pv_b(4)-0.25d0*t_pt_d_pv_b(4)
        I_03(5)=m4pn*d_pt_t_pv_b(5)-0.25d0*t_pt_d_pv_b(5)
        I_03(6)=m4pn*d_pt_t_pv_b(6)-0.25d0*t_pt_d_pv_b(6)


        !a2_d2 = a2+d2;
        !temp=1./a2_d2;
        !a2_p_d2=halfthickness2+d2
        a2_p_d2=bspread2+d2
        inv_a2_p_d2=1.d0/a2_p_d2

        !Ra    = sqrt( [a2_d2 a2_d2] + s.*s);
        !Rainv=1./Ra;
        !Ra3inv=Rainv.*Rainv.*Rainv;
        Ra1=a2_p_d2+s1*s1
        Ra1_inv=1.d0/sqrt(Ra1)
        Ra1_3inv=Ra1_inv/Ra1 ! 1/ + cher que 2* ? OUI
        s1Ra3inv=s1*Ra1_3inv

        !Ra1_3inv=Ra1_inv*Ra1_inv*Ra1_inv
        Ra2=a2_p_d2+s2*s2
        Ra2_inv=1.d0/sqrt(Ra2)
        Ra2_3inv=Ra2_inv/Ra2 ! 1/ + cher que 2* ? ENFIN PAREIL (dans la barre d'erreur)
        !Ra2_3inv=Ra2_inv*Ra2_inv*Ra2_inv
        s2Ra3inv=s2*Ra2_3inv

        s1_03=s1*Ra1_inv*inv_a2_p_d2
        s2_03=s2*Ra2_inv*inv_a2_p_d2
        s_03=s2_03-s1_03


        ! Si les segments sont vis pures on reduit le calcul
        if(TYSEG(IVLI) == 1) then

            SIGINTERN(1) = SIGINTERN(1) + I_03(1) * s_03
            SIGINTERN(2) = SIGINTERN(2) + I_03(2) * s_03
            SIGINTERN(3) = SIGINTERN(3) + I_03(3) * s_03
            SIGINTERN(4) = SIGINTERN(4) + I_03(4) * s_03
            SIGINTERN(5) = SIGINTERN(5) + I_03(5) * s_03
            SIGINTERN(6) = SIGINTERN(6) + I_03(6) * s_03

            cycle

        endif


        !!t_pt_t=Tens6SelfVec(T) ! tabulable
        t_pt_t(1) = Tbis(1)*Tbis(1)
        t_pt_t(2) = Tbis(2)*Tbis(2)
        t_pt_t(3) = Tbis(3)*Tbis(3)
        t_pt_t(4) = Tbis(2)*Tbis(3)
        t_pt_t(5) = Tbis(1)*Tbis(3)
        t_pt_t(6) = Tbis(1)*Tbis(2)

        !!d_pt_d=Tens6SelfVec(cai_d)
        d_pt_d(1) = cai_d(1)*cai_d(1)
        d_pt_d(2) = cai_d(2)*cai_d(2)
        d_pt_d(3) = cai_d(3)*cai_d(3)
        d_pt_d(4) = cai_d(2)*cai_d(3)
        d_pt_d(5) = cai_d(1)*cai_d(3)
        d_pt_d(6) = cai_d(1)*cai_d(2)

        commun=m4pn*d_pv_b_ps_t

        I_13=-mn4pn*t_pt_t_pv_b(:)

        I_25=commun*t_pt_t


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



        !!I_15=a2m8p*t_pt_t_pv_b-commun*t_pt_d
        I_15(1)=a2m8p*t_pt_t_pv_b(1)-commun*t_pt_d(1)
        I_15(2)=a2m8p*t_pt_t_pv_b(2)-commun*t_pt_d(2)
        I_15(3)=a2m8p*t_pt_t_pv_b(3)-commun*t_pt_d(3)
        I_15(4)=a2m8p*t_pt_t_pv_b(4)-commun*t_pt_d(4)
        I_15(5)=a2m8p*t_pt_t_pv_b(5)-commun*t_pt_d(5)
        I_15(6)=a2m8p*t_pt_t_pv_b(6)-commun*t_pt_d(6)



        !!I_05=commun*(core2*TENS6ID+d_pt_d)-a2m8p*t_pt_d_pv_b
        I_05(1)=commun*(bspread2+d_pt_d(1))-a2m8p*t_pt_d_pv_b(1)
        I_05(2)=commun*(bspread2+d_pt_d(2))-a2m8p*t_pt_d_pv_b(2)
        I_05(3)=commun*(bspread2+d_pt_d(3))-a2m8p*t_pt_d_pv_b(3)
        I_05(4)=commun*d_pt_d(4)-a2m8p*t_pt_d_pv_b(4)
        I_05(5)=commun*d_pt_d(5)-a2m8p*t_pt_d_pv_b(5)
        I_05(6)=commun*d_pt_d(6)-a2m8p*t_pt_d_pv_b(6)

        SIGINTERN(1) = SIGINTERN(1) + (I_03(1)*s_03 + I_13(1)*s_13 + I_05(1)*s_05 + I_15(1)*s_15 + I_25(1)*s_25)
        SIGINTERN(2) = SIGINTERN(2) + (I_03(2)*s_03 + I_13(2)*s_13 + I_05(2)*s_05 + I_15(2)*s_15 + I_25(2)*s_25)
        SIGINTERN(3) = SIGINTERN(3) + (I_03(3)*s_03 + I_13(3)*s_13 + I_05(3)*s_05 + I_15(3)*s_15 + I_25(3)*s_25)
        SIGINTERN(4) = SIGINTERN(4) + (I_03(4)*s_03 + I_13(4)*s_13 + I_05(4)*s_05 + I_15(4)*s_15 + I_25(4)*s_25)
        SIGINTERN(5) = SIGINTERN(5) + (I_03(5)*s_03 + I_13(5)*s_13 + I_05(5)*s_05 + I_15(5)*s_15 + I_25(5)*s_25)
        SIGINTERN(6) = SIGINTERN(6) + (I_03(6)*s_03 + I_13(6)*s_13 + I_05(6)*s_05 + I_15(6)*s_15 + I_25(6)*s_25)

      enddo

      ! Stockage des valeurs (Nb_seg*6 composantes)
      SIGJOINT(:,IS) = SIGINTERN(:)*bdivpa

    enddo LJIini

  endif

  ! on ne doit plus faire du point a point car tous les procs d'une phase ont besoin
  ! de connaitre SIGJOINT pour tous les segments
  ! en effet, le proc de rang pairimpair 0 dans une phase ne gere peut etre pas
  ! le segment gere par le proc de rang pairimpair similaire dans l'autre phase

  ! on va donc sommer la contribution de chaque proc a SIGJOINT sur COMM_PAIRIMPAIR pour chaque
  ! segment (sachant que les procs ne s'etant pas occupes d'un segment voit SIGJOINT=0 pour ce
  ! segment).
  CALL MPI_ALLREDUCE(SIGJOINT,SIGJOINT_Sum,6*Nb_seg,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_PAIRIMPAIR,IERR)

  ! Puis on va faire un send receive point a point pour reechanger SIGJOINT entre les 2 phases.
  CALL MPI_SENDRECV(SIGJOINT_Sum,6*Nb_seg,MPI_DOUBLE_PRECISION,RANG_DEST,121,                                   &
                    SIGJOINT_Recv,6*Nb_seg_stored,MPI_DOUBLE_PRECISION,RANG_DEST,121,MPI_COMM_WORLD,MPI_STATUS_IGNORE,IERR)

  ! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  ! THE RESULT OF THE SEGMENT CLOSE TO A BOUNDARY IS SAVED
  SIGMA_JOINT(:,:,1:NSEGM) = zero
  do IS = 1 ,Nb_seg_stored
    SIGMA_JOINT(1,1,ID_RO(IS)) = SIGJOINT_Recv(1,IS)
    SIGMA_JOINT(2,1,ID_RO(IS)) = SIGJOINT_Recv(6,IS)
    SIGMA_JOINT(3,2,ID_RO(IS)) = SIGJOINT_Recv(4,IS)
    SIGMA_JOINT(3,1,ID_RO(IS)) = SIGJOINT_Recv(5,IS)
    SIGMA_JOINT(1,1,ID_RO(IS)) = SIGJOINT_Recv(1,IS)
    SIGMA_JOINT(2,2,ID_RO(IS)) = SIGJOINT_Recv(2,IS)
    SIGMA_JOINT(3,3,ID_RO(IS)) = SIGJOINT_Recv(3,IS)
    SIGMA_JOINT(1,2,ID_RO(IS)) = SIGJOINT_Recv(6,IS)
    SIGMA_JOINT(2,3,ID_RO(IS)) = SIGJOINT_Recv(4,IS)
    SIGMA_JOINT(1,3,ID_RO(IS)) = SIGJOINT_Recv(5,IS)
  enddo

  ! We do not need any more the following list
  deallocate(liste_RO_Recv)
  deallocate(ID_RO)
  deallocate(SIGJOINT)
  deallocate(SIGJOINT_Sum)
  deallocate(SIGJOINT_Recv)

endif ! end of the case Nb_phase==2

# endif

#ifdef PA

prep0 = real(MPI_WTIME(),DP)       ! Time at the end of the stress calculation

!**************!
! deallocation
!**************!
deallocate(NSEGBP)
deallocate(COUTSEG)
deallocate(Charge)
deallocate(SEGDEB)
deallocate(BDEB)
deallocate(NSEG)
deallocate(Vitesse_p)

!Time dedicated to the stress calculation
calc = (prep0-prep1)
if (calc <= zero) calc= 1.

prep1 = real(MPI_WTIME(),DP)
prep = prep + (prep1-prep0)     ! Cumulated time of preparation for the parallel calculation

10 format(A21,I2,A6,F12.5,A6,F12.5)
if (modulo(KK,KISAUVE)==0) then
  write (*,10) "MPI: SIGMA_CP;  PID: ", MON_RANG_PAIRIMPAIR, " Cal: ", calc, " Vit: ", prep/calc
endif

#endif

end subroutine SIGMA_INT_CP

!###########################################################################
!# Subroutine SIGMA_INT_LP calculate the long-range stress field           #
!# contribution with the "box" algorithm. All calculation done we know     #
!# SIGBOX(:,:,IX,IY,IZ) the stress at the center of box domains.           #
!###########################################################################
subroutine SIGMA_INT_LP

implicit none

! Parameters

integer (kind=DPI)                  :: I,IB,JB,carac,indice,VLI,PBCIx,PBCIy,PBCIz
integer (kind=DPI), dimension(3)    :: PBCTr
real(kind=DP), dimension(3)         :: RBOX
real(kind=DP), dimension(3)         :: R        !< The center of the domain where the stress is calculated
real(kind=DP), dimension(3)         :: Tr       !< The translation vector for the modulo operations
real(kind=DP), dimension(3)         :: TrABmod  !<
real(kind=DP), dimension(3)         :: ROOmod   !< Coordinate of the domain center after modulo translation
real(kind=DP), dimension(3)         :: XModuR   !< The modur dimension in real format
real(kind=DP), dimension(3)         :: RA,RB
real(kind=DP), dimension(3)         :: Tbis,Bbis,BPVTbis
real(kind=DP), dimension(6)         :: SIGINTERN

! variable for the calculation of non-singular dislocation stress field (NSDSF)
real(kind=DP) :: d2,s1,s2,a2_p_d2,inv_a2_p_d2
real(kind=DP) :: s1Ra3inv,Ra1_3inv,Ra1_inv
real(kind=DP) :: s2Ra3inv,Ra2_3inv,Ra2_inv
real(kind=DP) :: s1_03,s1_13,s1_05,s1_15,s1_25
real(kind=DP) :: s2_03,s2_13,s2_05,s2_15,s2_25
real(kind=DP) :: s_03,s_13,s_05,s_15,s_25
real(kind=DP) :: m4pn,mn4pn,a2m8p
real(kind=DP) :: d_pv_b_ps_t,commun,R_ps_t,Ra1,Ra2
real(kind=DP),dimension(3) :: cai_d,d_pv_b,t_pv_b
real(kind=DP),dimension(6) :: t_pt_t
real(kind=DP),dimension(6) :: d_pt_d,t_pt_d,t_pt_t_pv_b,d_pt_t_pv_b,t_pt_d_pv_b
real(kind=DP),dimension(6) :: I_03,I_13,I_05,I_15,I_25

!PrintFieldVTK tabs
real(kind=DPI),allocatable      :: XCooR(:),YCoor(:),ZCoor(:)   !< Box Center Coordinates
real(kind=DPI),allocatable      :: FieldAmp(:)                  !< The field amplitude we want to trace
character(len=26)               :: file2print                   !< The name and place of the file to print

#ifdef PA

REAL(kind=DP),allocatable     :: Vitesse_p(:)
REAL(kind=DP),allocatable     :: PkSIGSend(:)
REAL(kind=DP),allocatable     :: PkSIGRecv(:)
INTEGER(KIND=DPI),allocatable :: COUTBOX(:)
INTEGER(KIND=DPI),allocatable :: Charge(:)

REAL(KIND=DP)       :: VTOTAL,INVVTOTAL
REAL(KIND=DP)       :: prep,prep0,prep1,calc
INTEGER(KIND=DPI)   :: CHRGD,CHRG,NCOUPLE,INTP,INTP2
INTEGER(KIND=DPI)   :: BFIN,IX,IY,IZ,BX,BY,BZ,I1,J1,K1

#endif

!initialization of NSDSF parameters
m4pn=0.25d0/(1.0d0-dpoiss)
mn4pn=m4pn*dpoiss
a2m8p=bspread2*0.125d0

! Initialization
SIGBOX(:,:,:) = ZERO
RBOX(:)       = tailleBoite(:) * Half       ! The box center
XModuR(:)     = ModuR(:)

#ifdef PA

prep0 = real(MPI_WTIME(),DP)       ! time a the beginning of the parallel preparation

!***************************!
! List of allocatable tables
!***************************!
allocate(Vitesse_p(TAILLE_PAIRIMPAIR))
allocate(Charge(TAILLE_PAIRIMPAIR))
allocate(COUTBOX(NBOITES))
allocate(PkSIGSend(6*NBOITES))
allocate(PkSIGRecv(6*NBOITES))

!***********************************************************!
! The reference velocity of the cores in parallel computing
!***********************************************************!
Vitesse_p(1:TAILLE_PAIRIMPAIR)=UN/real(TAILLE_PAIRIMPAIR,DP)   ! We guess that the speed of all the procs is identical
VTotal = UN

!***************************************!
! calculation dispatch between the procs
!***************************************!
 COUTBOX(1:NBOITES) = NSEGM

! The number of operations to be made is function of segments distribution in domain
! with the multi-poles algorithm
do IB = 1, NBOITES

  IX=B1D_3D(IB,1)-IUN
  IY=B1D_3D(IB,2)-IUN
  IZ=B1D_3D(IB,3)-IUN
  INTP2=0

  do I1 = -IUN , IUN
    do J1 = -IUN , IUN
      do K1 =  -IUN , IUN
          bx = modulo(int(IX+I1),int(NBoitesX))+IUN
          by = modulo(int(IY+J1),int(NBoitesY))+IUN
          bz = modulo(int(IZ+K1),int(NBoitesZ))+IUN
          INTP  = BX + ((BY-1)+ (BZ-1)*NBoitesY)*NBoitesX
          INTP2 = INTP2 + NSEGBOITE(INTP)
      enddo
    enddo
  enddo

  if (NSEGBOITE(IB) > IUN) COUTBOX(IB) = INTP2

enddo

! We add one virtual segment in each box to account for the empty domains
 COUTBOX(1:NBOITES)= NSEGM + 1 - COUTBOX(1:NBOITES)

! The total number of operations
NCOUPLE = SUM(COUTBOX(1:NBOITES))

! Calculation of the amount of work to be made by each procs
InvVTotal = NCouple / VTotal

Charge(1:(TAILLE_PAIRIMPAIR-1)) = nint(InvVTotal * Vitesse_p(1:(TAILLE_PAIRIMPAIR-1)),DPI)
Charge(TAILLE_PAIRIMPAIR) = NCOUPLE - SUM(Charge(1:(TAILLE_PAIRIMPAIR - 1)))  ! The last core finish the work

IB=1

do I = 1, (TAILLE_PAIRIMPAIR-1)
  CHRG = Charge(I)

  do while(CHRG.GT.0 .AND. CHRG.GE.COUTBOX(IB))
    CHRG = CHRG - COUTBOX(IB)
    IB = IB+1
  enddo

  Charge(I) = Charge(I) - CHRG
  Charge(I+1) = Charge(I+1) + CHRG
enddo

! Calculation of the amount of work to be made by the procs (MON_RANG_PAIRIMPAIR + 1)
CHRGD = SUM(Charge(1:MON_RANG_PAIRIMPAIR))
CHRG = Charge(MON_RANG_PAIRIMPAIR + 1)
PkSIGSend(:)=0

prep1 = real(MPI_WTIME(),DP)             ! time after preparation
prep = (prep1-prep0)            ! new parallel calculation extra cost

#endif

!***********************************************************************************!
!* Loop on the boxes where the long-range stress field contribution must be defined !
!***********************************************************************************!
Boite : do IB = 1,Nboites

#ifdef PA
    IF (CHRG==0) EXIT Boite
    IF (CHRGD.GE.COUTBOX(IB)) THEN
       CHRGD=CHRGD-COUTBOX(IB)
       CYCLE
    ENDIF
    CHRG = CHRG-COUTBOX(IB)

    ! for bi-phase calculations the internal stress must be defined in
    ! all the multipole domains
    if (Nb_phase /= 2) then
#endif

    ! We do not compute the stress at the center of an empty box
    if(NsegBoite(IB) < IUN) cycle BOITE

#ifdef PA
    endif
#endif

    SIGINTERN(:) = zero

    ! Coordinates of the box center where the stress must be computed
    R(1) = real(B1D_3D(IB,1),DP) * tailleBoite(1) - RBOX(1)
    R(2) = real(B1D_3D(IB,2),DP) * tailleBoite(2) - RBOX(2)
    R(3) = real(B1D_3D(IB,3),DP) * tailleBoite(3) - RBOX(3)

    if (kkdebug) write(379,*) "Index of boite in the LR calculation",IB

    ! The translation vector associates to PBC for the box where the stress is calculated
    Tr(:) =  int(HalfModur(1:3)  - R(1:3),DPI)

    !***************************************************************!
    !* The list of boxes to consider as a source of internal stress !
    !***************************************************************!
    S_boite: do JB = 1,Nboites

      ! Initially it was proposed to pass calculation in the empty domain to speed up calculation
      ! This approximation is very dangerous as dislocation may come in a domain between 2 Sigma_int_LP call
      ! if(NsegBoite(JB) < IUN) Cycle S_Boite

      ! The PBC images contributions to the LR calculation
      ! The default solution is to reconstruct a more or less isotropic MAXMODUR volume of the simulated volume
      Do PBCIx = -PBCIxDIM,PBCIxDIM
        PBCTr(1) = PBCIx * Modur(1)
      Do PBCIy = -PBCIyDIM,PBCIyDIM
        PBCTr(2) = PBCIy * Modur(2)
      Do PBCIz = -PBCIzDIM,PBCIzDIM
        PBCTr(3) = PBCIz * Modur(3)

        ! Domains to eliminate
        if(any(IndBoiteVois(IB,:) == JB) .and.                        & ! Side domains to JB + JB
           ((abs(PBCTr(1))+abs(PBCTr(2))+abs(PBCTr(3))) == Izero))    & ! Domains in the reference volume
           cycle

        !************************************************************!
        !* The loop on segments in the box source of internal stress !
        !************************************************************!
        BI: do Indice =  1, NsegBoite(JB)

          I = IndexBoite(JB)%ListeSeg(indice)

          ! if (i /= seg(I)%voiso .or. i /= seg(i)%voise) then
          if (.not. out(i)) then     ! Segment to eliminate must not be considered

            ! Pivot segments and very long segment must be eliminated from the calculation
            if(ab_length(I) /= izero .AND. .not. longseg(I)) then
               VLI = VLseg(i)
               carac = TYSEG(VLI)

               ! Modulo are applied to the center of segment coordinates (to account for PBC)
               ROOmod(:)= HalfModur(:) - (modulo((ro(:,i) + Tr(:)),XModuR(:)) + real(PBCTr(:),DP))

               TrABmod(:) = real((AB_length(i)*bveclin(1:3,VLI)), DP) * half

               ! The A to R and the B to R vectors
               RA(:)=ROOmod(:)+TrABmod
               RB(:)=ROOmod(:)-TrABmod

               ! Tabulated information of the segment I
               Tbis(:)=T(:,I)
               Bbis(:)=B(:,I)
               BPVTbis(:)=BPVT(:,I)

               ! Les caracteristiques du segment utiles
               R_ps_t=RA(1)*Tbis(1)+RA(2)*Tbis(2)+RA(3)*Tbis(3)
               s1=-R_ps_t
               cai_d=RA(:)-R_ps_t*Tbis(:)
               d2=cai_d(1)*cai_d(1)+cai_d(2)*cai_d(2)+cai_d(3)*cai_d(3) ! distance au segment
               R_ps_t=RB(1)*Tbis(1)+RB(2)*Tbis(2)+RB(3)*Tbis(3)
               s2=-R_ps_t

               t_pv_b=-BPVTbis(:)

               !!d_pv_b=prodvec(cai_d,B)
               d_pv_b(1)=cai_d(2)*Bbis(3)-cai_d(3)*Bbis(2)
               d_pv_b(2)=cai_d(3)*Bbis(1)-cai_d(1)*Bbis(3)
               d_pv_b(3)=cai_d(1)*Bbis(2)-cai_d(2)*Bbis(1)

               d_pv_b_ps_t=d_pv_b(1)*Tbis(1)+d_pv_b(2)*Tbis(2)+d_pv_b(3)*Tbis(3)

               !!t_pt_d=SymTens6Vec(localT,cai_d)
               t_pt_d(1) = (Tbis(1)*cai_d(1))*2.d0
               t_pt_d(2) = (Tbis(2)*cai_d(2))*2.d0
               t_pt_d(3) = (Tbis(3)*cai_d(3))*2.d0
               t_pt_d(4) = Tbis(2)*cai_d(3)+cai_d(2)*Tbis(3)
               t_pt_d(5) = Tbis(1)*cai_d(3)+cai_d(1)*Tbis(3)
               t_pt_d(6) = Tbis(1)*cai_d(2)+cai_d(1)*Tbis(2)

               !!t_pt_t_pv_b=SymTens6Vec(T,t_pv_b)
               t_pt_t_pv_b(1) = (Tbis(1)*t_pv_b(1))*2.d0
               t_pt_t_pv_b(2) = (Tbis(2)*t_pv_b(2))*2.d0
               t_pt_t_pv_b(3) = (Tbis(3)*t_pv_b(3))*2.d0
               t_pt_t_pv_b(4) = Tbis(2)*t_pv_b(3)+t_pv_b(2)*Tbis(3)
               t_pt_t_pv_b(5) = Tbis(1)*t_pv_b(3)+t_pv_b(1)*Tbis(3)
               t_pt_t_pv_b(6) = Tbis(1)*t_pv_b(2)+t_pv_b(1)*Tbis(2)
               !!d_pt_t_pv_b=SymTens6Vec(cai_d,t_pv_b)
               d_pt_t_pv_b(1) = (cai_d(1)*t_pv_b(1))*2.d0
               d_pt_t_pv_b(2) = (cai_d(2)*t_pv_b(2))*2.d0
               d_pt_t_pv_b(3) = (cai_d(3)*t_pv_b(3))*2.d0
               d_pt_t_pv_b(4) = cai_d(2)*t_pv_b(3)+t_pv_b(2)*cai_d(3)
               d_pt_t_pv_b(5) = cai_d(1)*t_pv_b(3)+t_pv_b(1)*cai_d(3)
               d_pt_t_pv_b(6) = cai_d(1)*t_pv_b(2)+t_pv_b(1)*cai_d(2)
               !!t_pt_d_pv_b=SymTens6Vec(localT,d_pv_b)
               t_pt_d_pv_b(1) = (Tbis(1)*d_pv_b(1))*2.d0
               t_pt_d_pv_b(2) = (Tbis(2)*d_pv_b(2))*2.d0
               t_pt_d_pv_b(3) = (Tbis(3)*d_pv_b(3))*2.d0
               t_pt_d_pv_b(4) = Tbis(2)*d_pv_b(3)+d_pv_b(2)*Tbis(3)
               t_pt_d_pv_b(5) = Tbis(1)*d_pv_b(3)+d_pv_b(1)*Tbis(3)
               t_pt_d_pv_b(6) = Tbis(1)*d_pv_b(2)+d_pv_b(1)*Tbis(2)

               !!I_03=m4pn*(d_pv_b_ps_t*TENS6ID+d_pt_t_pv_b)-m4p*t_pt_d_pv_b
               I_03(1)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(1))-0.25d0*t_pt_d_pv_b(1)
               I_03(2)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(2))-0.25d0*t_pt_d_pv_b(2)
               I_03(3)=m4pn*(d_pv_b_ps_t+d_pt_t_pv_b(3))-0.25d0*t_pt_d_pv_b(3)
               I_03(4)=m4pn*d_pt_t_pv_b(4)-0.25d0*t_pt_d_pv_b(4)
               I_03(5)=m4pn*d_pt_t_pv_b(5)-0.25d0*t_pt_d_pv_b(5)
               I_03(6)=m4pn*d_pt_t_pv_b(6)-0.25d0*t_pt_d_pv_b(6)

               !a2_d2 = a2+d2;
               !temp=1./a2_d2;
               !a2_p_d2=halfthickness2+d2
               a2_p_d2=bspread2+d2
               inv_a2_p_d2=1.d0/a2_p_d2
               !Ra    = sqrt( [a2_d2 a2_d2] + s.*s);
               !Rainv=1./Ra;
               !Ra3inv=Rainv.*Rainv.*Rainv;
               Ra1=a2_p_d2+s1*s1
               Ra1_inv=1.d0/sqrt(Ra1)
               Ra1_3inv=Ra1_inv/Ra1 ! 1/ + cher que 2* ? OUI
               s1Ra3inv=s1*Ra1_3inv

               !Ra1_3inv=Ra1_inv*Ra1_inv*Ra1_inv
               Ra2=a2_p_d2+s2*s2
               Ra2_inv=1.d0/sqrt(Ra2)
               Ra2_3inv=Ra2_inv/Ra2 ! 1/ + cher que 2* ? ENFIN PAREIL (dans la barre d'erreur)
               !Ra2_3inv=Ra2_inv*Ra2_inv*Ra2_inv
               s2Ra3inv=s2*Ra2_3inv
               s1_03=s1*Ra1_inv*inv_a2_p_d2
               s2_03=s2*Ra2_inv*inv_a2_p_d2
               s_03=s2_03-s1_03

               ! Si les segments sont vis pures on reduit le calcul
               if(carac == 1) then

                 SIGINTERN(1) = SIGINTERN(1) + I_03(1)*s_03
                 SIGINTERN(2) = SIGINTERN(2) + I_03(2)*s_03
                 SIGINTERN(3) = SIGINTERN(3) + I_03(3)*s_03
                 SIGINTERN(4) = SIGINTERN(4) + I_03(4)*s_03
                 SIGINTERN(5) = SIGINTERN(5) + I_03(5)*s_03
                 SIGINTERN(6) = SIGINTERN(6) + I_03(6)*s_03
                 cycle

               endif

               !!t_pt_t=Tens6SelfVec(T) ! tabulable
               t_pt_t(1) = Tbis(1)*Tbis(1)
               t_pt_t(2) = Tbis(2)*Tbis(2)
               t_pt_t(3) = Tbis(3)*Tbis(3)
               t_pt_t(4) = Tbis(2)*Tbis(3)
               t_pt_t(5) = Tbis(1)*Tbis(3)
               t_pt_t(6) = Tbis(1)*Tbis(2)

               !!d_pt_d=Tens6SelfVec(cai_d)
               d_pt_d(1) = cai_d(1)*cai_d(1)
               d_pt_d(2) = cai_d(2)*cai_d(2)
               d_pt_d(3) = cai_d(3)*cai_d(3)
               d_pt_d(4) = cai_d(2)*cai_d(3)
               d_pt_d(5) = cai_d(1)*cai_d(3)
               d_pt_d(6) = cai_d(1)*cai_d(2)
               commun=m4pn*d_pv_b_ps_t

               I_13=-mn4pn*t_pt_t_pv_b(:)
               I_25=commun*t_pt_t

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

               !!I_15=a2m8p*t_pt_t_pv_b-commun*t_pt_d
               I_15(1)=a2m8p*t_pt_t_pv_b(1)-commun*t_pt_d(1)
               I_15(2)=a2m8p*t_pt_t_pv_b(2)-commun*t_pt_d(2)
               I_15(3)=a2m8p*t_pt_t_pv_b(3)-commun*t_pt_d(3)
               I_15(4)=a2m8p*t_pt_t_pv_b(4)-commun*t_pt_d(4)
               I_15(5)=a2m8p*t_pt_t_pv_b(5)-commun*t_pt_d(5)
               I_15(6)=a2m8p*t_pt_t_pv_b(6)-commun*t_pt_d(6)

               !!I_05=commun*(core2*TENS6ID+d_pt_d)-a2m8p*t_pt_d_pv_b
               I_05(1)=commun*(bspread2+d_pt_d(1))-a2m8p*t_pt_d_pv_b(1)
               I_05(2)=commun*(bspread2+d_pt_d(2))-a2m8p*t_pt_d_pv_b(2)
               I_05(3)=commun*(bspread2+d_pt_d(3))-a2m8p*t_pt_d_pv_b(3)
               I_05(4)=commun*d_pt_d(4)-a2m8p*t_pt_d_pv_b(4)
               I_05(5)=commun*d_pt_d(5)-a2m8p*t_pt_d_pv_b(5)
               I_05(6)=commun*d_pt_d(6)-a2m8p*t_pt_d_pv_b(6)

               SIGINTERN(1) = SIGINTERN(1) + (I_03(1)*s_03 + I_13(1)*s_13 + I_05(1)*s_05 + I_15(1)*s_15 + I_25(1)*s_25)
               SIGINTERN(2) = SIGINTERN(2) + (I_03(2)*s_03 + I_13(2)*s_13 + I_05(2)*s_05 + I_15(2)*s_15 + I_25(2)*s_25)
               SIGINTERN(3) = SIGINTERN(3) + (I_03(3)*s_03 + I_13(3)*s_13 + I_05(3)*s_05 + I_15(3)*s_15 + I_25(3)*s_25)
               SIGINTERN(4) = SIGINTERN(4) + (I_03(4)*s_03 + I_13(4)*s_13 + I_05(4)*s_05 + I_15(4)*s_15 + I_25(4)*s_25)
               SIGINTERN(5) = SIGINTERN(5) + (I_03(5)*s_03 + I_13(5)*s_13 + I_05(5)*s_05 + I_15(5)*s_15 + I_25(5)*s_25)
               SIGINTERN(6) = SIGINTERN(6) + (I_03(6)*s_03 + I_13(6)*s_13 + I_05(6)*s_05 + I_15(6)*s_15 + I_25(6)*s_25)

            endif

          endif

        enddo BI

      enddo   ! The PBC Long-range contribution loop
      enddo
      enddo

    enddo S_boite

    ! A bdiva renormalization is needed here with the W Cai formula
    SIGINTERN(:) = SIGINTERN(:) * bdivpa

#ifdef PA
    IF (IZERO .NE. IUN) THEN

      ! pour un calcul parallele il faut toujours cette solution
      BFIN=6*IB
      PkSIGSend(BFIN-5:BFIN)=SIGINTERN(:)

    ELSE
#endif

      SIGTMP(2,1)=SIGINTERN(6) ; SIGTMP(3,2)=SIGINTERN(4) ; SIGTMP(3,1)=SIGINTERN(5)
      SIGTMP(1,1)=SIGINTERN(1) ; SIGTMP(2,2)=SIGINTERN(2) ; SIGTMP(3,3)=SIGINTERN(3)
      SIGTMP(1,2)=SIGINTERN(6) ; SIGTMP(2,3)=SIGINTERN(4) ; SIGTMP(1,3)=SIGINTERN(5)

      SIGBOX(:,:,IB) = SIGTMP(:,:)

#ifdef PA
    ENDIF
#endif

enddo Boite

#ifdef PA

prep0 = real(MPI_WTIME(),DP)           ! time at the end of the stress calculation
calc = (prep0 - prep1)              ! time associate to the stress calculation
if (calc <= zero) calc= 1.

! distribution of the stress result calculation between procs
CALL MPI_ALLREDUCE(PkSIGSend,PkSIGRecv,6*NBOITES,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

INTP=0

do IB=1,NBOITES
  SIGBOX(2,1,IB)=PkSIGRecv(INTP+6)
  SIGBOX(1,2,IB)=PkSIGRecv(INTP+6)
  SIGBOX(3,2,IB)=PkSIGRecv(INTP+4)
  SIGBOX(2,3,IB)=PkSIGRecv(INTP+4)
  SIGBOX(1,3,IB)=PkSIGRecv(INTP+5)
  SIGBOX(3,1,IB)=PkSIGRecv(INTP+5)
  SIGBOX(1,1,IB)=PkSIGRecv(INTP+1)
  SIGBOX(2,2,IB)=PkSIGRecv(INTP+2)
  SIGBOX(3,3,IB)=PkSIGRecv(INTP+3)
  INTP=INTP+6
enddo

! extra tabs needed for the parallel calculation can be deallocated
deallocate(Vitesse_p)
deallocate(Charge)
deallocate(COUTBOX)
deallocate(PkSIGSend)
deallocate(PkSIGRecv)

prep1 = real(MPI_WTIME(),DP)               ! time clock to evaluate extra calculation associated to the parallel calculation
prep = prep + (prep1-prep0)       ! sum of the extra parallel calculation

! Summary information is now ploted !!!
10 format(A21,I2,A6,F12.5,A6,F12.5)
if (mod(KK-1,KISAUVE)==0) then
  write (*,10) "MPI: SIGMA_LP;  PID: ", MON_RANG_PAIRIMPAIR,"; Cal:",calc,"; Vit:", prep/calc
endif

#endif

! When debuging we may wants to visualize some components of the stress field in the box grid
if (kkdebug) then

  ! The box field tabs are defined
  allocate(XCooR(NBOITES))
  allocate(YCooR(NBOITES))
  allocate(ZCooR(NBOITES))
  allocate(FieldAmp(NBOITES))

  do IB = 1, NBOITES
    ! Coordinates of the box center where the stress must be computed
    XCooR(IB) = real(B1D_3D(IB,1),DP) * tailleBoite(1) - RBOX(1)
    YCooR(IB) = real(B1D_3D(IB,2),DP) * tailleBoite(2) - RBOX(2)
    ZCooR(IB) = real(B1D_3D(IB,3),DP) * tailleBoite(3) - RBOX(3)
    FieldAmp(IB) = SIGBOX(3,3,IB)                             ! The field we wants to plot
  enddo

  file2print = '../out/debug/DataField.vtk'
  write(379,*) "Information on the SIGMA_INT_LP subroutine is saved in the debug directory"
  call PrintFieldVTK(NBOITES,XCooR,YCoor,ZCoor,FieldAmp,file2print)

  deallocate(XCooR)
  deallocate(YCooR)
  deallocate(ZCooR)
  deallocate(FieldAmp)

endif

end subroutine SIGMA_INT_LP

!###########################################################################
!# To annihilate plastic rotation around the tensile axis, we impose a torque field
!# opposed to the rotations generated by slip activity
!# Isotropic elasticity and Lame constants (Lambda and Mu) are used
!###########################################################################
subroutine PBCCompStrain(Torque)

implicit none

Logical           :: SENS
Integer           :: I,J
real(kind=DP)     :: NewTensRot(3,3),Torque(3,3),NewTorque(3,3),Sig6(6),Eps6(6)
real(kind=DP)     :: TorqueMax,SigAppMax,TSRatio

Torque(:,:) = zero

! Initialization
Sig6(:) = zero

! We go from the crystal system to the new reference frame with the Z axis parallel
! to the tensile axis
SENS=.true.

!Call TurnTens(ZPrime,TensRot,NewTensRot,SENS)
Call TurnTens(ZPrime,TensDist,NewTensRot,SENS)

! if(MOD(KK,KSTAT).eq.IZERO) then
!    print*,'NewTensRot'
!    print*,NewTensRot(1,1),NewTensRot(1,2),NewTensRot(1,3)
!    print*,NewTensRot(2,1),NewTensRot(2,2),NewTensRot(2,3)
!    print*,NewTensRot(3,1),NewTensRot(3,2),NewTensRot(3,3)
! endif

! For simplicity Einstein notation is used
Eps6(1) = -NewTensRot(1,1)
Eps6(2) = -NewTensRot(2,2)
Eps6(3) = -NewTensRot(3,3)
Eps6(4) = -NewTensRot(2,3)
Eps6(5) = -NewTensRot(1,3)
Eps6(6) = -NewTensRot(1,2)

! We defined a virtual stress tensor opposite to the rigid body rotation
Do I = 1,6
   Do J = 1,6
      Sig6(I) = Sig6(I) + Stiffness(I,J) * Eps6(J)
   Enddo
Enddo

! Diagonal coeficients must equal zero in the torque tensor
NewTorque(1,1) = zero
NewTorque(2,2) = zero
NewTorque(3,3) = zero

! A reverse torque field is defined in two direction normal to the tensile axis
! The stress is in Mu reduce unit

! NewTorque(2,3) = zero
! NewTorque(1,3) = zero
NewTorque(1,2) = zero
NewTorque(2,3) = Sig6(4) * InvXMu
NewTorque(1,3) = - Sig6(5) * InvXMu
! NewTorque(1,2) = Sig6(6) * InvXMu

NewTorque(3,2) = NewTorque(2,3)
NewTorque(3,1) = NewTorque(1,3)
NewTorque(2,1) = NewTorque(1,2)

! We go back from the new reference frame to the crystal system
SENS=.false.
Call TurnTens(Zprime,NewTorque,Torque,SENS)

if(KKdebug) then
  write(379,*) 'The torque tensor associate to zprime rotations'
  write(379,*) Torque(1,1),Torque(1,2),Torque(1,3)
  write(379,*) Torque(2,1),Torque(2,2),Torque(2,3)
  write(379,*) Torque(3,1),Torque(3,2),Torque(3,3)
endif

! A saturation value is imposed at 0.10 max of sigapp
TorqueMax = max(abs(Torque(1,1)),abs(Torque(2,2)),abs(Torque(3,3)),   &
               &abs(Torque(2,3)),abs(Torque(1,3)),abs(Torque(1,2)))
SigAppMax = max(abs(SigApp(1,1)),abs(SigApp(2,2)),abs(SigApp(3,3)),   &
               &abs(SigApp(2,3)),abs(SigApp(1,3)),abs(SigApp(1,2)))
if (SigAppMax /= zero) then
  TSRatio = TorqueMax/SigAppMax
  ! The stress correction can only be a small fraction of sigapp
  if (TSRatio > 0.1) then
    Torque(:,:) = 0.1 * Torque(:,:) / TSRatio
  endif
else
  Torque(:,:) = zero
endif

end subroutine PBCCompStrain

!###########################################################################
!# Computation of the local line tension existing on segments touching a   #
!# freesurface. The solution used here is the same one as K. Schwartz and  #
!# referenced in the Hirth and Lothe.                                      #
!###########################################################################
function surfLT(i,vli,anglelb,Linedir)

implicit none

integer(kind=DPI)  :: I   ,&  !< Segment indice
                      J   ,&  !< Varfreeplan indice
                      vli     !< Segment line vector

real(kind=DP)      :: lambda(3)     ,&
                      surfLTvec(3)  ,&
                      n1(3)         ,&
                      n2(3)         ,&
                      Linedir(3)    ,&
                      nglpl(3)      ,&
                      surfLT        ,&
                      ANGLElb       ,&
                      dcosB         ,& !< Cosinus of B (degree)
                      dsinB         ,& !< Sinus of B (degree)
                      sintheta      ,&
                      costheta      ,&
                      tantheta      ,&
                      BurgersOfI(3) ,& !< Burger line vector of segment I
                      surfnorm(3)      !< Normal of the surface

J=seg(i)%VarFreePlan
BurgersOfI(1:3) = BvecLin(1:3,assoc(vli,1))

if (seg(i)%norme > DEUX*DCFJ) then
  lambda = DCFJ
else
  lambda = seg(i)%norme*half
endif

surfnorm=Plane_MillerR(1:3,j) !normal to the surface j

surfnorm=normavec(surfnorm)

nglpl(1:3) = BVECNOR(1:3,SEG(i)%VECLIN) !normal to the glide plane

 costheta = cose2v(linedir,surfnorm) !angle between the normal to the surface and the real dislocation line

if (costheta ==ZERO) then
  costheta=0.00001
endif

if (costheta < ZERO) then
    costheta = -costheta
    linedir = -linedir
endif

sintheta = dsqrt(UN - costheta*costheta)
tantheta = sintheta / costheta

n1 = prodvec(surfnorm,linedir(1:3)) !the normal to the surface and the real dislocation line lie on this plane
n1 = prodvec(linedir(1:3),n1(1:3)) !
dcosB = dcos(anglelb)
dsinB = dsqrt(UN - dcosB*dcosB)
n2 = prodvec(nglpl,linedir)
n2 = n2(1:3) / dsqrt(n2(1)*n2(1) + n2(2)*n2(2) + n2(3)*n2(3))
if (dot_product(n2,BurgersOfI(1:3))<IZERO) then
 n2 = -n2
endif

if (dot_product(n1,surfnorm)<-1e-15) then
   print *, 'wrong sign'
   read(*,*)
endif
if (dot_product(n2,BurgersOfI(1:3))<-1e-15) then
   print *, 'wrong sign n2'
   read(*,*)
endif



surfLTvec = BdivA/(4*PII*(UN-dpoiss)*lambda)*(n1*abs((UN-dpoiss*dcosB*dcosB)*tantheta)+n2*(DEUX*dpoiss*dsinB*dcosB))
surfLT = (VECnorDEP(1,VLI)* surfLTvec(1)+ &
          VECnorDEP(2,VLI)* surfLTvec(2)+ &
          VECnorDEP(3,VLI)* surfLTvec(3))


end function surfLT

#ifdef MDC
!###############################################################################
!# Send to Zebulon the location of the points where the stress has to be known #
!###############################################################################
subroutine send_stresspoints_location

implicit none

integer(kind=DPI)  :: i

call zebmpipack(int(nsegmdcstress))

do i=1,nsegmdcstress
     call zebmpipack(int(ro(1:3,listsegmdcstress(i))),3)
enddo

call zebmpipacksend

end subroutine send_stresspoints_location

!###############################################################################
!# Get from Zebulon the long range stresses                                    #
!###############################################################################
subroutine recv_stresspoints_stress

implicit none

integer(kind=DPI) :: iseg,i
real(kind=DP), dimension(6*nsegmdcstress)  :: mdclrstress


mdclrstress(:) = ZERO
#if defined(MDC) && defined(PA)

if (Mon_Rang ==  IZERO) then
#endif

  call zebmpipackrecv

  call zebmpiunpack_doublearray(mdclrstress, int(6*nsegmdcstress))

  mdclrstress(1:(6*nsegmdcstress))=mdclrstress(1:(6*nsegmdcstress))*1.e6*unsurxmu;

#if defined(PA) && defined(MDC)
endif

CALL MPI_BCAST(mdclrstress,6*nsegmdcstress,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

#endif

!!! The stress tensor received by Zebulon is (sxx, syy, szz, sqrt(2)sxy,sqrt(2)*syx,sqrt(2)sxz)
!!! to get the right value in microMegas we MUST divide the off diagonal components by sqrt(2)

do i=1,nsegmdcstress
  iseg=listsegmdcstress(i)

  SIG_EF(1,1,iseg)=mdclrstress(6*(i-1)+1)
  SIG_EF(2,2,iseg)=mdclrstress(6*(i-1)+2)
  SIG_EF(3,3,iseg)=mdclrstress(6*(i-1)+3)
  SIG_EF(1,2,iseg)=mdclrstress(6*(i-1)+4)*halfsqrt2
  SIG_EF(2,3,iseg)=mdclrstress(6*(i-1)+5)*halfsqrt2
  SIG_EF(3,1,iseg)=mdclrstress(6*(i-1)+6)*halfsqrt2
  SIG_EF(2,1,iseg)=SIG_EF(1,2,iseg)
  SIG_EF(3,2,iseg)=SIG_EF(2,3,iseg)
  SIG_EF(1,3,iseg)=SIG_EF(3,1,iseg)


enddo

end subroutine recv_stresspoints_stress


subroutine initialize_MDC_var

implicit none

real(kind=DP)      :: ytemp(ntsg)

endcalcul=0
nbswept=0

if (sideja == 0) then

  nstep = 0  ! Initialise the number of steps to zero
!  raudis = zero
  raudmo = zero
  vitmoy = zero
  epsdot = zero

endif

halfthickness2  = real((halfthickness*halfthickness),DP)
bdivpa = VecBurgers /PII/avalue
SR_cutoff=real(2.5d0*halfthickness,DP)
SR_cutoff2=SR_cutoff*SR_cutoff

TAUIII = TAUIII/XMU
Ytemp(1:ntsg) = 350*BDIVA**3*AVALUE**3    !Volume d activation du Cu: V=350*b**3
UNKT_force = XMU*Ytemp(1)*1.0D0/(1.38D-23*TEMPERATURE)
ARGEXP = 1.0D15*BETA*DELTAT*avalue
unsurxmu = un/xmu


end subroutine initialize_MDC_var
#endif

#if defined(MDC) && defined(PA)
subroutine sig_FEsincro

implicit none

real(kind=DP) :: SIGFETEMP_S(1:6*nsegmax)
integer(kind=DPI) :: iseg,i

SIGFETEMP_S(1:6*nsegmax)=ZERO

Call MPI_ALLREDUCE(SIGFETEMP,SIGFETEMP_S,6*nsegmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

do i=1,nsegmdcstress
  iseg=listsegmdcstress(i)
  seg(iseg)%SIGFE(1)=SIGFETEMP_S(6*(iseg-1)+1)
  seg(iseg)%SIGFE(2)=SIGFETEMP_S(6*(iseg-1)+2)
  seg(iseg)%SIGFE(3)=SIGFETEMP_S(6*(iseg-1)+3)
  seg(iseg)%SIGFE(4)=SIGFETEMP_S(6*(iseg-1)+4)
  seg(iseg)%SIGFE(5)=SIGFETEMP_S(6*(iseg-1)+5)
  seg(iseg)%SIGFE(6)=SIGFETEMP_S(6*(iseg-1)+6)
enddo
end subroutine sig_FEsincro
#endif

!################################################################################
!# Calculation of an elastic, 2D , stress field generated by a crack
!# For more details about the stress intensity factor, please see Appendix A.1 in
!# 'Fatigue of materials', 2nd Edition, author : S. Suresh
!# The crack shape used is the standard shape (R,theta) of a crack in a isotropic material
!################################################################################
subroutine calculate_crackfield(Rtemp)

implicit none

real(kind=DPI) :: radiusVect(3)   !< Vector from crack tip to considered segment center

real(kind=DP) :: radiusNorm       !< Distance of the segment center from crack tip
real(kind=DP) :: theta            !< Angle between radiusVect and (0,1,0)
real(kind=DP) :: SIF           !< Stress intensity factor (\f$MPa*\sqrt(mm)$\f)
real(kind=DP) :: SIF_reduced   !< 'Effective' Stress intensity factor (\f$SIFact/sqrt(2*\pi*r)$\f)

integer (kind=DPI),dimension(3)   :: Rtemp  !< Coordinates where the stress is calculated

integer(kind=DPI)    :: ROx,ROy,ROz         !< Segment center coordinates

!real(kind=DP) :: Sigvm  !< Von Mises Stress

!Segment Position
ROX = Rtemp(1)
ROY = Rtemp(2)
ROZ = Rtemp(3)

cracktip(1)=ROX                    !So that field does not depend on segment x-coordinate
radiusVect(1)=real(ROX,8)-cracktip(1) !Zero obviously
radiusVect(2)=real(ROY,8)-cracktip(2)
radiusVect(3)=real(ROZ,8)-cracktip(3)

!Polar Coordinates from Cartesian coordinates
Call RPOLARS( radiusVect(2), radiusVect(3), radiusNorm, THETA)

SIF = SIGMA*sqrt(acrack*PII)*SIF_fact

SIF_reduced = SIF/sqrt(DEUX*PII*radiusNorm)

Sigcrack(2) =  SIF_reduced * cos(theta/DEUX) * (1 - sin(theta/DEUX) * sin(TROIS*theta/DEUX))
Sigcrack(3) =  SIF_reduced * cos(theta/DEUX) * (1 + sin(theta/DEUX) * sin(TROIS*theta/DEUX))
Sigcrack(5) =  SIF_reduced * cos(theta/DEUX) * (    sin(theta/DEUX) * sin(TROIS*theta/DEUX))
Sigcrack(1) =  DPOISS *(Sigcrack(2)+Sigcrack(3)) !Plane Strain
!Sigcrack(1) =  ZERO                             !Plane Stress

Sigcrack(4) = ZERO
Sigcrack(6) = ZERO

!Von mises equivalent stress (for information)
!Sigvm = sqrt(HALF*((Sigcrack(1,1)-Sigcrack(2,2))**2 + (Sigcrack(2,2)-Sigcrack(3,3))**2 &
!+(Sigcrack(3,3)-Sigcrack(1,1))**2 + 6.0d0*(Sigcrack(1,2)**2 + Sigcrack(2,3)**2 + Sigcrack(1,3)**2)))

end subroutine calculate_crackfield

end module ELASTI

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
