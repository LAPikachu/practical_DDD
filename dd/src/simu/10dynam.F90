!===================================================================================================
!========================    DEBUT   MODULE   "DYNAM"    ===========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the dislocation mobility rules and laws.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module DYNAM

use CONSTANTES
use DEBUG       !*** Encore indispensable......
use VARGLOB
use CONNEC

implicit none

real(kind=DP)      :: RDEP_PREDI(1:NSEGMax)

contains

!###########################################################
!# This function provide the result of mobility laws for   #
!# materials as a function of slip systems and dislocation #
!# caracteres.                                             #
!###########################################################
function vitesse(taueff,longueur,nloi,sys,NPhase)

real(kind = DP)       :: longueur, vitesse, taueff, facteur_cran, tmp, work_app
integer(kind = DPI)   :: nloi, sys, Nphase

!initialization
vitesse = zero

if(nloi > NLV .or. nloi < 1 ) then
  print *, kk, " erreur in the mobility law : undefined mobility law"
  return
endif

if(.not. effet_longueur .and. longueur > ModurMax) then
  print *, kk, " erreur in the mobility law : segment length > volume dimensions "
elseif(effet_longueur .and. longueur > ModurMax) then
  longueur = ModurMax
endif

if(longueur < un) then
  print *, kk, " erreur in the mobility law : segment length = 0  ", longueur
  return
endif

if(Loi(nloi)%arrhenius > 2 .OR. Loi(nloi)%arrhenius < 0) then
  print *, kk," type de loi de vitesse inconnu de VITESSE "
  return
endif

!=====================================================================================
!* During relaxation steps a standard damping mobility law is applied (mobility law 1)
!=====================================================================================
if (KK <= relax_reac .and. Loi(IUN)%arrhenius == 0) then
  ! Be aware that dry friction is applied during the relaxation steps
  taueff = taueff - Loi(IUN)%friction
  if(taueff > zero) then
    vitesse = taueff*VecBurgers/loi(IUN)%BF
  endif
  return
endif

!====================================================================
!  Dry friction stress contribution is taken into account
!====================================================================
taueff = taueff - Loi(nloi)%friction

!====================================================================
! If there is no enough driving force the velocity is simply zero
!====================================================================
if(taueff <= zero) return

!====================================================================
! Selection of the mobility law to consider
!====================================================================
if(Loi(nloi)%arrhenius == izero) then

  ! The viscous law case

  if(NPhase /= izero) then
    vitesse = taueff*VecBurgers/(2.*loi(nloi)%bf) ! I do not anderstand this line ??? BD
    print*,'strange mobility law !'
  else
    vitesse = taueff*VecBurgers/loi(nloi)%bf
  endif

else

  ! The thermally activated law case

  ! The work_app calculation
  if(taueff >= loi(nloi)%tau0) then ! case of error in the value of taueff

    work_app = Loi(nLoi)%deltaG0
    !      print *, "work_app ",work_app, taueff,loi(nloi)%tau0

  else

    ! determination of the work provided by the effective stress : it is equal to deltaG0 - deltaG
    ! This energy should be determined separately because it involves a sinh function
    work_app = Loi(nLoi)%deltaG0*(1-(1-(Taueff/loi(nloi)%tau0)**Loi(nloi)%coef_p)**Loi(nLoi)%coef_q)
    !      print *, "work_app nor",work_app

  endif

  if(Loi(nloi)%arrhenius == IUN) then

    ! case of thermally activated law, but length - independent mobility
    vitesse = Loi(nloi)%V0 * (exp(-UnkT * (Loi(nLoi)%deltaG0 - work_app)) -  &
                              exp(-UnkT * (Loi(nLoi)%deltaG0 + work_app)))

  elseif (Loi(nloi)%arrhenius == IDEUX) then

    ! calcul du nombre de cran deja formes sur longueur
    tmp = abs(longueur * avalue * CranSys(sys))
    ! si efficacite = 1, la longueur est divise en consequence
    Facteur_cran = un / (un + Effica_crans * tmp)
    ! facteur_ecran vaux 1 si efficacite ecran = 0
    if (kkdebug .and. tmp >= un) write(379,*) "nc =",tmp, " Fac", Facteur_cran
    longueur = longueur * facteur_cran

    vitesse = Loi(nloi)%H * longueur * (exp(-UnkT * (Loi(nLoi)%deltaG0 - work_app)) -  &
                                        exp(-UnkT * (Loi(nLoi)%deltaG0 + work_app)))

!     print *, " Loi(nloi)%H ", Loi(nloi)%H, "work_app",work_app
!     print *,"---", taueff,longueur*avalue*1E6,nloi
!     print *, " Loi(nloi)%pre ", Loi(nloi)%H*longueur*avalue, " exp ", &
!     exp(-Loi(nLoi)%deltaG0 * UnkT) * DEUX * DSINH(work_app* UnkT),vitesse

  else

    write(*,*) " Unkown type of mobility law number : ", nloi
    stop

  endif

endif

end function vitesse

!######################################################################
!# Calcul du deplacement des segments de dislocation a partir d'une   #
!# loi des vitesses de type frotement visqueux sans tenir compte      #
!# d'eventuelle reactions locales entre dislocation.                  #
!#                                                                    #
!# IDEP(I) : deplacement entier exprime en vecteur deplacement        #
!# SEG(I)%RDEP_PREDI : reste reel du deplacement                      #
!# SEG(I)VIT : vitesse instantannee du segment en metre par seconde   #
!########################################################## 15/11/98 ##
!> Calculation of the segments velocity and displacement in free flight condition
subroutine DEPPREDIC

implicit none

integer (kind=DPI) :: I,J,K,i1,i2,ii1,compteur,m,ibv,nloi,list(1000)
integer            :: Irdep_predi(nsegmax), Irdep_Index(nsegmax)                          !< Integer intermediate tables needed for optimization
integer (kind=DPI) :: nscrew, nedge                                                       !< number of screw and edge segments
real(kind=DP)      :: VITES,RESD,fac,Lcumul,DepInst,mufoisadivb
real(kind=DP)      :: taueff,taueffcumul,unsurb,longueur,signe,temp,lboite,vitmin_nloi1
real(kind=DP)      :: l_i1,vites_i1_ele(1000)
real(kind=DP)      :: length3seg                !< Length of dislocation section made of 3 connected segments
logical            :: condition1,condition2,condition3,zip(nsegm)

#ifdef PA

real(kind=DP),allocatable   :: DEPICOM(:)
real(kind=DP),allocatable   :: DEPICOMR(:)
real(kind=DP),allocatable   :: TAUTOTCOM(:)
real(kind=DP),allocatable   :: TAUTOTCOMR(:)
real(kind=DP),allocatable   :: ANGLEVISCOM(:)
real(kind=DP),allocatable   :: ANGLEVISCOMR(:)
real,allocatable            :: PROBAGDCOM(:)
real,allocatable            :: PROBAGDCOMR(:)
integer,allocatable         :: GDVECLCOM(:)
integer,allocatable         :: GDVECLCOMR(:)
integer,allocatable         :: LOISEGCOM(:)
integer,allocatable         :: LOISEGCOMR(:)

! Communication tabs needed for the creep mode
real(kind=DP),allocatable   :: TauMonCOM(:)
real(kind=DP),allocatable   :: TauMonCOMR(:)


! Initialization
allocate (DEPICOM(NSEGM))
allocate (DEPICOMR(NSEGM))
DEPICOM(1:NSEGM)  = ZERO
DEPICOMR(1:NSEGM) = ZERO

allocate (TAUTOTCOM(NSEGM))
allocate (TAUTOTCOMR(NSEGM))
TAUTOTCOM(1:NSEGM)  = ZERO
TAUTOTCOMR(1:NSEGM) = ZERO

allocate (ANGLEVISCOM(NSEGM))
allocate (ANGLEVISCOMR(NSEGM))
ANGLEVISCOM(1:NSEGM) = ZERO
ANGLEVISCOMR(1:NSEGM) =ZERO

allocate (LOISEGCOM(NSEGM))
allocate (LOISEGCOMR(NSEGM))
LOISEGCOM(1:NSEGM)  = IZERO
LOISEGCOMR(1:NSEGM) = IZERO


allocate (PROBAGDCOM(NSEGM))
allocate (PROBAGDCOMR(NSEGM))
PROBAGDCOM(1:NSEGM)=ZERO
PROBAGDCOMR(1:NSEGM)=ZERO

allocate (GDVECLCOM(NSEGM))
allocate (GDVECLCOMR(NSEGM))
GDVECLCOM(1:NSEGM)=IZERO
GDVECLCOMR(1:NSEGM)=IZERO

allocate (TauMonCOM(NSEGM))
allocate (TauMonCOMR(NSEGM))
TauMonCOM(1:NSEGM)=ZERO
TauMonCOMR(1:NSEGM)=ZERO

#endif

! Initialization
mufoisadivb = XMU * avalue / vecburgers
IDEP(1:NSEGMAX) = IZERO                     ! new segments may be added during the step, so
RDEP_PREDI(1:NSEGMAX) = ZERO                ! nsegmax must be considered for the initialization.
fac = facteur_cris/deltat
Lseg(1:NLV) = ZERO
D_brute(1:NLV) = ZERO                       ! Non corrected average displacement (velocity) per mobility law
unsurb = 1.0 / VecBurgers
vitmin_nloi1 = 0.5 * float(modep_max) / Deltat

! call distribution_crans

if (kkdebug) write(379,*)  " Dlimite =  ",dlimite(1:NLV)

nscrew = izero
nedge = izero
L_average_screw = zero
L_average_edge = zero
rau_screw = zero
rau_edge = zero
tau_average_screw  = zero
tau_average_edge = zero


!*** First loop on the segments

B1:  do I = 1, NSEGM

#ifdef PA
  if (monseg(I)) then

  LoiSegCom(I) = LoiSeg(I)      ! the loiseg status must be loaded for all monseg(I) segments
#endif

  ! Those segments are not supposed to move
  if (seg(i)%norme == 0 .or. seg(i)%jonc .or. seg(i)%gd > izero .or. out(i) .or. &
        seg(i)%voiso  == nsegmax .or. seg(i)%voise == nsegmax) then
    cycle B1
  endif

  ! Initializations
  vites = zero
  IBV = seg(I)%veclin
  longueur = real(SEG(I)%NORME)  ! The default value of the i segment length

  ! We check that a velocity law exist for this segment
  if (Loiseg(I) < 1 .or.Loiseg(I) > NLV) then
    call conf(i)
    print *, " LoiSeg(i =",i,")  not correctly initialized :", Loiseg(I)
    Loiseg(I) = int(numero_loi(IBV))
  endif
  ! The number of the velocity law of I is defined in subroutine Force
  nloi = LoiSeg(I)

    ! This segment is expected to move at this step
    if(seg(i)%wait < itrois .or. VL_CS(I) > izero) then

      ! The sign of the effective stress is saved
      signe = Dsign(UN,tauTOT(I))
      ! The abs value of the effective stress in real units
      taueff = dabs(tauTOT(I)) * Xmu

      ! A segment can move only if Taueff > zero
      if(taueff > ZERO) then

        ! The screw direction associate to the slip system of i
        k = assoc(IBV,1)

        ! For some mobility laws, the effective length in special character direction must be defined
        if (effet_longueur .and. seg(i)%anglevis < angle_vis) then
          ! Segment i is part of a screw section, only its screw length parallel to b must be considered
          if (PartVis(IBV) /= zero) then
            Longueur = longueur * normlin(k) * PartVis(IBV)
          else
            ! In this case the segment is assumed of elementary length! Since it is part of a screw line, it real
            ! velocity is defined latter (see loop B3)
            Longueur = normlin(k)
          endif
        else
          ! The length of the segment i is considered as the real length
          longueur = longueur * normlin(IBV)
        endif

        m = Syseg(IBV)

        ! Segments defined as "singulier" are moved independently of their mobility laws
        if (.not. singulier(I)) then

          !> The segment velocity evaluation (in ABS value)
          vites = vitesse(taueff,longueur,nloi,m,seg(i)%NPhase)

          ! For the materials with a dislocation mobility thermally activated. The mobility law applied to segments
          ! not parallel to the direction with the large lattice friction is usually artificial (see Loi(nloi)%arrhenius == 1)
          ! At low effective stress, the exp form used in this law is sometime source of artifacts. It decrease too much
          ! the segment mobility. To avoid that problem, a minimum segment velocity is imposed in the dynamics.
          if (Loi(nloi)%arrhenius == 1) then
            if (vites < vitmin_nloi1) vites = vitmin_nloi1
          endif

          !> The displacement evaluation
          vites = vites * signe
          RDEP_PREDI(I) = VITES*DELTAT/NormDep(IBV)

          if(kkdebug) write(379,*) "ini I, rdep ini, Tau",i,";", RDEP_PREDI(I),signe*taueff*1D-6

          ! if rdep_predi is too large, the resulting integer (Idep) can be any value
          ! the magnitude is thus reduced to less than the double precision upper limit value 3 10E9
          if(abs(RDEP_PREDI(I)) > 1E5) RDEP_PREDI(I) = 1.0E5 * dsign(UN,vites)

        else

          ! The case of singular "singulier" segments.
          ! Be aware that in this case, we add 2.D9 to rdep_predi in order to keep the information
          ! of singularity in the parallel computation

          ! A special treatment is made for the segments touching GDs just after cross-slip event
          ! to help the bow out of dislocation in the cross-slip plane
          if ((seg(seg(i)%vnno)%gd + seg(seg(i)%vnne)%gd) > izero &
              .and. abs(seg(i)%anglevis) < 2.*angle_vis) then

            ! When the applied stress is very large, its contribution to cross-slip bow out cannot be neglected
            if (abs(tauAPP(i)) > 0.01) then
              ! At large stress, the total effective stress can be considered
              temp = (tauAPP(i) + tauINT(I) + tauTL(I)) * xmu
            else
              ! At low stress, only the internal stress is taken into account to help initial cross-slip bow out
              temp = tauINT(I) * xmu
            endif

          else

            if (seg(i)%surface > IZERO) then
              temp = (tauAPP(i) + tauINT(I) + tauTL(I)) * xmu
            else
              ! Generale case:
              ! When segments are closed to a singular field (very large stress value), they must be close to an other
              ! dislocation line and then a contact reaction is expected. In the singular field "region", the contribution
              ! of the applied stress is not considered to help the contact reaction to proceed rapidly
              temp = (tauINT(I) + tauTL(I)) * xmu
            endif

          endif

          !> The displacement of singular segments
          if(abs(temp) > 1.D0) then

            ! Here the displacement is made huge and proportional to the stress. This is to move singular segments
            ! in a decreasing velocity order and then to break symmetry in the dynamics of line section where several
            ! consecutive segments are defined as singular.
            if(dabs(temp) < 1.D9) then
              RDEP_PREDI(I) = (dabs(temp)*1.D-9 + modep_max) * dsign(UN,temp)
            else
              ! For those segments it is useless to defined any logic in the segment displacement
              RDEP_PREDI(I) = deux * modep_max * dsign(UN,temp)
            endif

          else

            ! This segment is defined as singular, but its stress is tiny. Hence the segment displacement is set to zero
            RDEP_PREDI(I) = zero

          endif

        endif

      endif

    endif

#ifdef PA

    DEPICOM(I) = RDEP_PREDI(I)          !Diffuse the DepInst tab between nodes outside the do loop
    TAUTOTCOM(I) = TAUTOT(I)            !Diffuse the effective stress tab between nodes outside the do loop
    ANGLEVISCOM(I) = seg(i)%anglevis    !Diffuse the anglevis information

    if (GLDEV) then
      GDVECLCOM(I) = int(VL_CS(I),4)           !Diffuse the GDTABVEC tab between nodes outside the do loop
      PROBAGDCOM(I)= real(seg(i)%probadev)     !Diffuse the CS proba between nodes
    endif

    if (creep) then
      TauMonCOM(i) = Taumon(i)          !Diffuse the TauMon tab between nodes outside the do loop
    endif

    ! With the following trick, two information are simultaneously diffused between nodes: the Loiseg type and
    ! the singulier character. Use is made of the property that loiseg can only be > 0
    if (singulier(i)) then
      LOISEGCOM(I) = -LoiSeg(I)         ! The segment i is singular
    endif

  endif     ! end of monseg(i) selection

#endif

enddo B1

#ifdef PA

CALL MPI_ALLREDUCE(DEPICOM,DEPICOMR,NSEGM,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_PAIRIMPAIR,IERR)
CALL MPI_ALLREDUCE(TAUTOTCOM,TAUTOTCOMR,NSEGM,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_PAIRIMPAIR,IERR)
CALL MPI_ALLREDUCE(ANGLEVISCOM,ANGLEVISCOMR,NSEGM,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_PAIRIMPAIR,IERR)
CALL MPI_ALLREDUCE(LOISEGCOM,LOISEGCOMR,NSEGM,MPI_INTEGER,MPI_SUM,COMM_PAIRIMPAIR,IERR)

if (GLDEV) then
  CALL MPI_ALLREDUCE(GDVECLCOM,GDVECLCOMR,NSEGM,MPI_INTEGER,MPI_SUM,COMM_PAIRIMPAIR,IERR)
  CALL MPI_ALLREDUCE(PROBAGDCOM,PROBAGDCOMR,NSEGM,MPI_REAL,MPI_SUM,COMM_PAIRIMPAIR,IERR)
endif

if (creep) then
  !> Communication of the climb stress between procs for parallel computation \ingroup creep
  CALL MPI_ALLREDUCE(TauMonCOM,TauMonCOMR,NSEGM,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_PAIRIMPAIR,IERR)
endif

RDEP_PREDI(1:NSEGM) = DEPICOMR(1:NSEGM)
TAUTOT(1:NSEGM) = TAUTOTCOMR(1:NSEGM)
seg(1:nsegm)%anglevis = ANGLEVISCOMR(1:NSEGM)

deallocate (DEPICOM)
deallocate (DEPICOMR)
deallocate (TAUTOTCOM)
deallocate (TAUTOTCOMR)
deallocate (ANGLEVISCOM)
deallocate (ANGLEVISCOMR)

if (GLDEV) then
  VL_CS(1:nsegm) = GDVECLCOMR(1:NSEGM)
  seg(1:nsegm)%probadev = PROBAGDCOMR(1:nsegm)
endif

deallocate (PROBAGDCOM)
deallocate (PROBAGDCOMR)
deallocate (GDVECLCOM)
deallocate (GDVECLCOMR)

if (creep) then
  TauMon(1:nsegm) = TauMonCOMR(1:NSEGM)
endif

deallocate (TauMonCOM)
deallocate (TauMonCOMR)

! We extract the information on the loiseg type
LoiSeg(1:nsegm) = iabs(LOISEGCOMR(1:nsegm))
! We extract the information on the singulier character
do i=1,nsegm
  singulier(i) = .false.
  if (LOISEGCOMR(i) < izero) singulier(i) = .true.
enddo

deallocate (LOISEGCOM)
deallocate (LOISEGCOMR)

#endif

! This part is dedicated to the calculation of an effective stress and effective length on segments
! with a mobility law thermally activated. This part is not // since segment neighbors at long distance
! can be involved
if (effet_longueur) then

  ! Initialization
  zip(1:nsegm) = .false.

  ! The effective length of dislocation cannot be larger than the simulated volume dimension
  lboite = real(sqrt3 * ModurMax) / normlin(1)

  B3: Do i = 1 , nsegm

    ! This segment as already been tested
    if (zip(i)) cycle B3

    ! We first cycle segments that does not need to be tested
    if (seg(i)%norme == 0 .or.                                  &
        seg(i)%jonc .or. seg(i)%gd > izero .or.                 &
        seg(i)%voiso == nsegmax .or. seg(i)%voise == nsegmax )    cycle B3

    ! The segment is now tested
    zip(i) = .true.

    ! Segment useful info
    IBV = seg(i)%veclin
    lcumul = seg(i)%norme

    ! This segment must be considered
    if (seg(i)%anglevis < angle_vis .and. .not. singulier(i)) then

      lcumul = lcumul * partVis(IBV)
      taueffcumul = tauTOT(i)
      ii1 = 0
      list(1:1000) = -9999
      i1 = seg(i)%voiso
      vites_i1_ele(1:1000) = zero

      k = assoc(IBV,1)
      m = Syseg(IBV)
      nloi = numero_loi(k)

      ! We search in a do loop all neighbors on the o side of i with the same character
      do while (Lcumul < lboite     .and.         &
                i1 /=  nsegmax      .and.         &
                i1 /= i             .and.         &
                .not. seg(i1)%jonc  .and.         &
                seg(i1)%gd < iun )

        ! We reach the end of the tested line section
        if (seg(i1)%norme /= izero .and. (seg(i1)%anglevis > Angle_vis) .or. zip(i1) .or. singulier(i1)) exit

        if(seg(i1)%norme > 0) then

          ii1 = ii1 + 1
          list(ii1) = i1

          ! The length and stress of segment is sum up
          lcumul = lcumul + seg(i1)%norme * partVis(seg(i1)%veclin)
          taueffcumul = taueffcumul + tauTOT(i1)

          ! The mobility of an elementary section of line one this segment i1 is calculated
          taueff = dabs(tauTOT(i1))  * Xmu
          l_i1 = seg(i1)%norme * partVis(seg(i1)%veclin) * normlin(k)
          if (l_i1 /= zero) vites_i1_ele(ii1) = vitesse(taueff,l_i1,nloi,m,seg(i)%Nphase)

        endif

        zip(i1) = .true.
        i1 = seg(i1)%voiso

!         if (ii1 > 1000) then
!           call seginfo(I,"infinit loop 1 in dynam - lcomu")
!           exit
!         endif

      enddo

      ! We search in a do loop all neighbors on the e side of i with the same character
      i1 = seg(i)%voise
      do while (Lcumul< lboite        .and.       &
                i1 /=  nsegmax        .and.       &
                i1 /= i               .and.       &
                .not. seg(i1)%jonc    .and.       &
                seg(i1)%gd < iun)

        ! We reach the end of the tested line section
        if (seg(i1)%norme > 0 .and. (seg(i1)%anglevis > Angle_vis) .or. zip(i1) .or. singulier(i1)) exit

        if(seg(i1)%norme > 0) then

          ii1 = ii1 + 1
          list(ii1) = i1

          ! The length and stress of segment is sum up
          lcumul = lcumul + seg(i1)%norme * partVis(seg(i1)%veclin)
          taueffcumul = taueffcumul + tauTOT(i1)

          ! The mobility of an elementary section of line one this segment i1 is calculated
          taueff = dabs(tauTOT(i1))  * Xmu
          l_i1 = seg(i1)%norme * partVis(seg(i1)%veclin) * normlin(k)
          if (l_i1 /= zero) vites_i1_ele(ii1) = vitesse(taueff,l_i1,nloi,m,seg(i)%NPhase)

        endif

        zip(i1) = .true.
        i1 = seg(i1)%voise

!         if (ii1 > 1000) then
!             call seginfo(i,"infinit loop 2 in dynam - lcomu")
!             exit
!         endif

      enddo

      if (ii1 > 1000) then
        call seginfo(i,"infinit loop in dynam - lcomu")
        exit
      endif

      ! The total number of segment in the effective screw length is ii1 + 1
      lcumul = lcumul * normlin(k)
      ! The new LoiSeg type for segment I is set
      Loiseg(i) = int(nloi)

      ! we restrict calculation of averages to the most active slip systems, i.e.
      ! for those Schmid factor is > 95% maximum Schmid factor.
      if (desorientGrain) then
        if (SchmidSysGrain(seg(i)%grain,syseg(IBV)) > 0.95 * Schmid) then
          nscrew = nscrew + IUN
          L_average_screw = L_average_screw + lcumul
          rau_screw = rau_screw + lcumul
          tau_average_screw = tau_average_screw + tauTOT(i) * Xmu
        endif
      else
        if (SchmidSys(syseg(IBV)) > 0.95 * Schmid) then
          nscrew = nscrew + IUN
          L_average_screw = L_average_screw + lcumul
          rau_screw = rau_screw + lcumul
          tau_average_screw = tau_average_screw + tauTOT(i) * Xmu
        endif
      endif

      ! The average effective resolved stress is calculated
      taueff = dabs(taueffcumul / (ii1 +1))  * Xmu
      ! From this value we compute a new average displacement per segments
      signe = Dsign(un,taueff)

      ! The particular case of one isolate segment with partVis(IBV) = 0
      ! This segment is in average in the screw direction, but is too small to be considered as screw
      if(lcumul == zero) then
        lcumul  = seg(i)%norme
        nloi    = numero_loi(IBV)
      endif

      ! Two possibilities
      ! 1: The screw section glide at the average stress velocity
          vites = vitesse(taueff,lcumul,nloi,m,seg(i)%NPhase)

      ! print*,"1",vites
      ! 2: The screw section glide at a velocity accounting for the variations of the kink-pair nucleation along the line
      !    l_i1 = seg(i)%norme * partVis(seg(i)%veclin) * normlin(k)
      !    vites_i1_ele(ii1+1) = vitesse((dabs(tauTOT(i))*Xmu),l_i1,nloi,m)
      ! The velocity of the total line section made of ii1+1 segments
      !    vites = sum(vites_i1_ele(1:ii1+1))
      ! print*,"2",vites,vitesse((dabs(tauTOT(i))*Xmu),l_i1,nloi,m)/l_i1*lcumul
      ! read(*,*)

      ! The average displacement evaluation
      vites = vites * signe
      rdep_predi(i) = VITES*DELTAT/NormDep(IBV)

      ! Now in the list() there is ii1 segment sharing the same screw direction as i
      ! therefore, we apply the same displacement to all these segments.
      do i1 = 1 , ii1
        rdep_predi(list(i1)) = rdep_predi(i)
      enddo

    else

      ! Statistics on the non screw mobile segments
      if (SchmidSys(syseg(IBV)) > 0.95 * Schmid) then
        lcumul = lcumul * normlin(IBV)
        nedge = nedge + IUN
        L_average_edge = L_average_edge + lcumul
        rau_edge = rau_edge + lcumul
        tau_average_edge = tau_average_edge + tauTOT(I) * Xmu
      endif

    endif

  enddo B3

  ! Calculation of the average length and stresses on screw and non srew segments
  if(nscrew == izero) nscrew = IUN
  if(nedge == izero) nedge = IUN
  L_average_screw = L_average_screw / nscrew * avalue * 1.E6  ! in microns
  L_average_edge = L_average_edge / nedge * avalue * 1.E6  ! in microns
  tau_average_screw = tau_average_screw / nscrew * 1.E-6   ! in MPa
  tau_average_edge = tau_average_edge / nedge * 1.E-6   ! in MPa

endif

if (GLDEV) then
    ! The following strange loop is needed to recover the information regarding nsysdev,
    ! indeed, in parallel computation and when the wait function is active, the nsysdev
    ! information is diffused when wait=2 betwwen nodes thanks to the VL_CS table.
    DO I=1,NSEGM
        if (VL_CS(I) < IZERO) then
            seg(i)%nsysdev = -VL_CS(I)
            !print*,'ici',I,mon_rang,VL_CS(I),seg(i)%nsysdev
            VL_CS(I) = IZERO
        endif
    ENDDO
endif

compteur = 0

DO I=1,NSEGM

   ! here, we get back the information about singularity
    if(rdep_predi(I) > 1.D9) then
       rdep_predi(I) = rdep_predi(I) - 2.D9
       compteur = compteur + 1
       seg(i)%resdep = zero
    endif

    if (seg(i)%norme == 0 .or. seg(i)%jonc .or. seg(i)%gd > izero .or. out(i) .or. &
         seg(i)%voiso  == nsegmax .or. seg(i)%voise == nsegmax ) cycle

    ! the length of line per type of velocity law is accumulated for the calculation of D_brute
    temp  = SEG(I)%NORME * NormLin(seg(I)%veclin)
    nloi =  Loiseg(I)

    Depinst = RDEP_PREDI(I)
    if (kkdebug) write(379,*) "i,rdep,res",  I,';',rdep_predi(I),seg(I)%resdep
    Resd =  seg(I)%resdep + rdep_predi(I)
    RDEP_PREDI(I) = abs(Resd)
    IDEP(I)   = nint(Resd,DPI)

    ! The predicted displacement of each segments, per type of velocity law is accumulated.
    ! Displacements are multiply per the segments length for the calculation of a mean value per unit length
    if (.not. seg(i)%bloquer .and. IDEP(I) /= izero) then
      i1 = seg(i)%vnno
      i2 = seg(i)%vnne

      ! The length of vnno and vnne must be included otherwise the segment length is reduced after displacement
      if (i1 /= nsegmax .and. i2 /= nsegmax) then
        length3seg = temp + SEG(i1)%NORME*NormLin(seg(i1)%veclin) + SEG(i2)%NORME*NormLin(seg(i2)%veclin)
      elseif (i1 /= nsegmax) then
        length3seg = temp + SEG(i1)%NORME*NormLin(seg(i1)%veclin)
      elseif (i2 /= nsegmax) then
        length3seg = temp + SEG(i2)%NORME*NormLin(seg(i2)%veclin)
      else
        length3seg = temp
      endif
      Lseg(nloi) = Lseg(nloi) + length3seg
      D_brute(nloi) = D_brute(nloi) + (length3seg * abs(real(IDEP(I),DP)))
    endif

    ! ----------    definition of resdep   ----------

    seg(I)%resdep = Resd - real(IDEP(I),DP)
    ! ----------    definition of resdep   ----------

    ! A paranoiac test is made here, to check that no large resdep cumul is made
    if (abs(seg(I)%resdep) > 2*modep_max) then
        print *, "error :debug:depassement en resdep kk,i=", kk,i,seg(I)%resdep
        stop
    endif

    ! We must free the singulier segments from the wait procedure
    if(singulier(i)) seg(i)%wait = izero

    !The particular cases whose Depinst was not calculated in B1 loop
    if (seg(i)%wait > IDEUX) then
       if (seg(i)%bloquer) then
          seg(i)%resdep = zero
          IDEP(I) =izero
       else
          Resd = seg(i)%resdep+seg(I)%DepInst
          IDEP(I) = nint(Resd,DPI)
          seg(I)%resdep = Resd - real(IDEP(I),DP)
       endif
    endif

    !********************************************
    !The waiting segment optimization algorithm
    if (Period /= IZERO) then
       if (kk > RELAX_reac) then
          !Different conditions to increment Wait
          !* condition1 beginnig conditions from wait=0 to wait=2
          condition1 = (seg(i)%wait < ITROIS .and. DepInst < (UN/Period))
          !* condition2 wait incrementation from wait=3 to wait=Period
          condition2 = (seg(I)%wait > IDEUX)
          !* condition3 for segments bloced in barriers
          condition3 =(seg(i)%wait < ITROIS .and. seg(i)%bloquer)
          if (condition1.or.condition2.or.condition3) then
             seg(i)%wait=seg(i)%wait+IUN
             !The DepInst is saved at step 3 of Wait and will serve as increment for the
             !following (period-3) steps
             if (seg(i)%wait<IQUATRE) then
                seg(i)%DepInst = DepInst
                if (seg(i)%wait<ITROIS) then
                  seg(i)%DepInst = zero
                endif
             endif
          else
             seg(i)%wait=IZERO    !The seg I is re-mobilized
          endif

          ! The wait algorithm must systematically be reinitialized when a segment move
          if (IDEP(I) /= IZERO) then
            seg(i)%wait=IZERO    !The seg I is re-mobilized
          endif
       endif
    endif
    !********************************************

ENDDO

! calculation of the mean predicted displacement per unit length (avalue)
do nloi = 1 ,nlv
  If(Lseg(nloi) /= izero) D_brute(nloi) = D_brute(nloi) / Lseg(nloi)
enddo

if (kkdebug) then
   if(Numero_loi(1) /= Numero_loi(3)) then
        write(379,'("Pour IVL=1, D_brute =",F10.5,", pour IVL=2, Dbrute =",F20.5)') &
        D_brute(numero_loi(1)),D_brute(numero_loi(3))
   else
        write(379,'("Pour IVL=1, D_brute =" ,F10.5)') D_brute(numero_loi(1))
   endif
endif

!******************************************
! Troisieme boucle : calcul du deplacement entiers des segemnts
! Calcul du deplacement entier resultant
!******************************************
compteur = 0
B2: do I = 1, nsegm

  if (Idep(I) == izero) cycle
  !**********************************************************************
  ! reduction de la vitesse : cas general
  !* Si les segments ont une vitesse elevee (i.e. donnant lieu a une
  !deplacement trop grand, on reduit l amplitude de son deplacement
  nloi = Loiseg(I)
  if(ABS(rdep_predi(I)) > Dlimite(nloi)) then
     compteur = compteur + 1
  !      print *,I,ABS(rdep_predi(I)),Dlimite(nloi)
     Idep(I) = nint(Dlimite(nloi),DPI) * SIGN(IUN,idep(i))
     seg(I)%resdep = zero
  endif

enddo B2

segments_fous = IZERO
if (nbdep /= 0) then
   segments_fous = nint(100.0 * real(compteur)/real(nbdep),DPI)
endif

!***********************************************************************
! The moving segments are listed in decreasing velocity order
! This ordering is made with integer since we want a final list
! machine and compiler independent!
!***********************************************************************
! A smaller list, with only moving segment, is made
nbdep = izero
do i = 1, nsegm

  if (IDEP(i) /= izero) then

    nbdep = nbdep + 1
    Irdep_index(nbdep) = int(i)
    Irdep_predi(nbdep) = int(abs(Idep(I)))

  endif

enddo

! Sorting is made on the reduce list with the comb sort algo
call combsort(Irdep_index(1:nbdep),Irdep_predi(1:nbdep))

! The final quidep descending list is made
j = nbdep
do i = 1, nbdep
  quidep(i) = Irdep_index(j)
  j = j-1
enddo

! ! If we want to check the ordered list
! do j = 1, nbdep
!   print*,'-->',j, quidep(j),rdep_predi(quidep(j)),seg(quidep(j))%norme
! enddo
! read(*,*)


! LAST LOOP :  REDO CONNECTION with cross_slipped segments
!*****************************************
!**** the cross slip actualization  ******
!*****************************************
if (GLDEV .and. kk > relax_reac) then
    do i=1,nsegm
        if (VL_CS(I) /= izero) then
            ! The two following information must be reset globally before making the
            ! connection around cross-sliping segments for parallel computing
            seg(i)%veclin = VL_CS(I)
            singulier(i) = .true.
            !print*,'la',I,mon_rang,VL_CS(I)
            !call seginfo(I, " avant erreur          ")
            ! We check that cross-sliping segment are of screw type
            if(tyseg(seg(i)%veclin) /= iun .or. tyseg(VL_CS(I)) /= iun) then
                 print *, " erreur de CS: non vis, I=",I, " KK =",kk
                 print*,seg(i)%veclin,tyseg(seg(i)%veclin),VL_CS(I),tyseg(VL_CS(I))
                 call seginfo(I, " avant erreur          ")
                 stop
            endif
        endif
    enddo
    do i=1,nsegm
        if (VL_CS(I) /= izero) then
            K = Iboite(I)
            nloi = nsegm + plusseg + IUN
            i1 = seg(i)%vnno
            i2 = seg(i)%vnne
            call connecIJ(i1,I,1)  ! The connectivity management
            call connecIJ(I,i2,2)  ! The connectivity management
            do J = nloi , Nsegm + PlusSeg
              IBOITE(J) = K
              NsegBoite(K) =  NsegBoite(K) + IUN
              IndexBoite(K)%ListeSeg(NsegBoite(K)) = J
              seg(J)%wait     = IZERO
              seg(J)%bloquer  = seg(I)%bloquer
              seg(J)%unload   = seg(I)%unload
              seg(J)%grain    = seg(I)%grain
              ! We keep the information that the new gd segments was mad by cross-slip
              if (seg(J)%gd > izero) seg(J)%gd = ideux
              ! if (seg(J)%gd > izero) print*,'1 gd2 segment was introduced', J, seg(J)%gd
            enddo
            ! print*,'VL_CS',mon_rang,I,VL_CS(I),plusseg
        endif
    enddo
endif

Idep(nsegm + iun : nsegm + plusseg) = izero
Nsegm=nsegm+plusseg      ! New segments are accounted for
plusseg=izero
if(kkdebug) write(379,*)  "test 6"
if (segments_fous > 20   .and. compteur > 10 .and. nsegm > 1000) then
  write(*,'("Fraction of too fast segments:",I4," %. Dlim(2:3) =",2F6.1," pas/iter","    compteur=",I7,"    nsegm=",I7)') &
  segments_fous,Dlimite(2:3),compteur,nsegm
endif

if (kkdebug)  then
  lcumul = avalue * 1E9
  write(379,'("nsegm: ", i6, ",  nbdep: ", i6)') nsegm,nbdep
  do m=1,nbdep
    i = quidep(m)
    write(379,343) m,i,seg(i)%norme*normlin(seg(i)%veclin)*lcumul,loiseg(i),idep(I),     &
        (idep(I)+ seg(i)%resdep)*normdep(seg(i)%veclin)*avalue/deltat*1E6,         &
        normdep(seg(i)%veclin)*avalue/deltat, seg(i)%resdep, seg(i)%wait
    if(abs(seg(i)%resdep) > IUN) then
      write(379,*)  kk,"Fatal: dynam:deppredic: For seg I=",I," Resdep ",seg(i)%resdep
      stop
    endif
  enddo

343 format("ordre:", I6, " seg:", I6, " L =", F10.2, " nm; loi:", I2, " dep:", I5,    &
         " V=", ES10.3, " micron/s ; Ve:", ES10.3, " m/s;resd:", F6.2, ", wait",I2)
endif

end subroutine DEPPREDIC


!############################################################
!# ???                                                      #
!############################################################
subroutine DISTRIBUTION_CRANS

integer (DPI) :: disO, DisE, i,dislo,NC,NS,reste,NCd,m,ndislo,compt
logical :: zipseg(nsegmax)
zipseg(1:nsegm+PlusSeg) = .false.
ndislo = 0
if(kkdebug) write(379,*)  " debut distribution .........."
! Les dislocation encrees
do dislo = 1, nsegm + PlusSeg
if (seg(dislo)%voiso == nsegmax .and. seg(dislo)%voise == nsegmax) zipseg(dislo) =.true.
if(zipseg(dislo)) cycle
m = syseg(seg(dislo)%veclin)
if(seg(dislo)%norme == 0 .and. seg(dislo)%voiso == dislo .and. seg(dislo)%voise == dislo) then
   zipseg(dislo) = .true.
elseif (seg(dislo)%voiso == nsegmax) then
   ndislo = ndislo + 1
   DisO = dislo
   i = dislo
   NC = izero
   NCd = izero
   NS = izero
   COMPT = 0
   Do while (i /= nsegmax .and. syseg(seg(i)%veclin) == m)
       if (zipseg(i)) then
          print *, " DisloInfo : double voisinage 1 "
          exit
       elseif(i == seg(i)%voise) then
          print*, " Distrib cran: seg a eliminer ????"
          exit
       elseif(seg(i)%norme /= 0)then
          !              NC = NC + seg(i)%ncran
          !              NS = NS + 1
       endif
       zipseg(i) = .true.
       !           if (test)print *, " dislo =", dislo, "  seg =", i , " ncran = ", seg(i)%ncran
       i = seg(i)%voise
       COMPT = COMPT + 1
       if (compt > 10000) PRINT *, "STOP kk=",kk, " CYCLE INFINIE 1 dans DIstibution_crans "
   enddo
   DisE = I
   if (NS == 0) exit
   NCd = NC / NS
   reste = NC - NCd*NS
   if(kkdebug) write(379,*)  " ENCREE =", dislo, "NC  = ", Nc,"NS = ", NS
   If (NCd > 0) then
      if (kkdebug) write(379,*)  NCd, "cran  sur :",NS, " .Reste :", reste
      i = DisO
      compt = 0
      Do while(i /= nsegmax)
          !              if (seg(i)%norme /= 0) seg(i)%ncran = NCd
          if (i == DisE) then
             !                 seg(i)%ncran = seg(i)%ncran + reste
             exit
          endif
          i = seg(i)%voise
          COMPT = COMPT + 1
          if (compt > 10000) PRINT *, "STOP kk=",kk, " CYCLE INFINIE 3 dans DIstibution_crans "
      enddo
   endif
elseif (syseg(seg(seg(dislo)%voiso)%veclin) /= m) then
   ndislo = ndislo + 1
   DisO = dislo
   i = dislo
   NC = izero
   NCd = izero
   NS = izero
   compt = 0
   Do while (i /= nsegmax .and. syseg(seg(i)%veclin) == m)
       if (zipseg(i)) then
          print *, " DisloInfo : double voisinage 1 "
          exit
       elseif(i == seg(i)%voise) then
          print*, " Distrib cran: seg a eliminer ????"
          exit
       elseif(seg(i)%norme /= 0)then
          !              NC = NC + seg(i)%ncran
          !              NS = NS + 1
       endif
       !           if (test)print *, " dislo =", dislo, "  seg =", i , " ncran = ", seg(i)%ncran
       zipseg(i) = .true.
       i = seg(i)%voise
       COMPT = COMPT + 1
       if (compt > 10000) PRINT *, "STOP kk=",kk, " CYCLE INFINIE 2 dans DIstibution_crans "
   enddo
   DisE = I
   if (NS == 0) exit
   NCd = NC / NS
   reste = NC - NCd*NS
   if(kkdebug) write(379,*)  " ENCREE =", dislo, "NC  = ", Nc,"NS = ", NS
   If (NCd > 0) then
      if (kkdebug) write(379,*)  NCd, "cran  sur :",NS, " .Reste :", reste
      i = DisO
      compt = 0
      Do while(i /= nsegmax)
          !              if (seg(i)%norme /= 0) seg(i)%ncran = NCd
          if (i == DisE) then
             !                 seg(i)%ncran = seg(i)%ncran + reste
             exit
          endif
          i = seg(i)%voise
          COMPT = COMPT + 1
          if (compt > 10000) PRINT *, "STOP kk=",kk, " CYCLE INFINIE 3 dans DIstibution_crans "
      enddo
   endif
endif
enddo

! Les boucles
do dislo = 1, nsegm + PlusSeg
if(zipseg(dislo)) cycle
if(kkdebug) write(379,*)  " dislo =", dislo, "  BOUCLE "
ndislo = ndislo + 1
NC = izero
NCd = izero
NS = izero
i = dislo
compt = 0
Do while (i /= nsegmax)
    if (zipseg(i)) then
       print *, " DisloInfo : double voisinage 1 "
       exit
    elseif(i == seg(i)%voise) then
       print*, " Distrib cran: seg a eliminer ????"
       exit
    elseif(seg(i)%norme /= 0)then
       !           NC = NC + seg(i)%ncran
       NS = NS + 1
    endif
    !         if (test)print *, " dislo =", dislo, "  seg =", i , " ncran = ", seg(i)%ncran
    zipseg(i) = .true.
    i = seg(i)%voise
    COMPT = COMPT + 1
    if (compt > 10000) PRINT *, "STOP kk=",kk, " CYCLE INFINIE 4 dans DIstibution_crans "
    if (i == dislo) exit
enddo
if (NS == 0) exit
NCd = NC / NS
reste = NC - NCd*NS
if(kkdebug) write(379,*)  " BOUCLE =", dislo, "NC  = ", Nc,"NS = ", NS
If (NCd > 0) then
   if (kkdebug) write(379,*)  NCd, "cran  sur :",NS, " .Reste :", reste
   i = Dislo
   compt = 0
   Do while(i /= nsegmax)
       !           if (seg(i)%norme /= 0) seg(i)%ncran = NCd
       i = seg(i)%voise
       COMPT = COMPT + 1
       if (compt > 10000) PRINT *, "STOP kk=",kk, " CYCLE INFINIE 5 dans DIstibution_crans "
       if (i == dislo) then
          !              seg(i)%ncran = seg(i)%ncran + reste
          exit
       endif
   enddo
endif
enddo
end subroutine Distribution_crans

end module DYNAM

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
