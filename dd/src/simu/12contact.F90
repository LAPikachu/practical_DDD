
!===================================================================================================
!========================    DEBUT    MODULE   "CONTACT"  ==========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to segment displacement rules and contact reaction rules
!> such as annihilation, junction formation, etc ...
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module CONTACT

!*** Ce module est la source des plus grands problemes de MMPM
!*** Ces derniers sont souvent lies a un des points suivants :
!*** Accumulation de rotules
!*** Description 2D a 4 vecteurs (soit 2 de trop)
!*** Problemes de sous reseau de 3 type <> : VDi VLi et VLj
!*** Les vecteurs ne sont plus necessairement // ou perpendicuaires donc
!*** les procedures CV meme modifies ne sont pas toujours a la hauteur

use BRICAMAT
use VARGLOB
use CONNEC
use DEBUG
use Microstructure
#ifdef MDC
use mdctools
#endif

implicit none

integer(kind=DPI) ,dimension(NSEGMAX) :: listeSegTester     !<
integer(kind=DPI) ,dimension(NSEGMAX) :: listetmp           !<
integer(kind=DPI) ,dimension(3)       :: Oi                 !<
integer(kind=DPI) ,dimension(3)       :: Ei                 !<
integer(kind=DPI) ,dimension(3)       :: VDi                !<
integer(kind=DPI) ,dimension(3)       :: VLi                !<
integer(kind=DPI) ,dimension(3)       :: VDi1o_or           !<
integer(kind=DPI) ,dimension(3)       :: VDi1e_or           !<
integer(kind=DPI) ,dimension(3)       :: NVec               !<

integer(kind=DPI) :: i                !<
integer(kind=DPI) :: Li               !<
integer(kind=DPI) :: LiO              !<
integer(kind=DPI) :: i1o              !< First neighboring segment from segment origin (Approx. the same as VOISO)
integer(kind=DPI) :: i1e              !< First neighboring segment from segment extremity (Approx. the same as VOISE)
integer(kind=DPI) :: i2e              !< Second neighboring segment from segment origin
integer(kind=DPI) :: i2o              !< Second neighboring segment from segment extremity
integer(kind=DPI) :: limite_discreti  !<
integer(kind=DPI) :: limite_bloquage  !<
integer(kind=DPI) :: modep            !<
integer(kind=DPI) :: normsqrVDi       !<
integer(kind=DPI) :: vnn1o            !< First non-zero neighboring segment from segment origin
integer(kind=DPI) :: vnn1e            !< First non-zero neighboring segment from segment extremity
integer(kind=DPI) :: vnn2o            !< Second non-zero neighboring segment from segment origin
integer(kind=DPI) :: vnn2e            !< Second non-zero neighboring segment from segment extremity
integer(kind=DPI) :: LongSegi         !<

real,parameter    :: NsegSafty = 1.1   !< used to stop the simulation when Nsegm is 90% of NSEGMAX

real(kind=DP)     :: distiobs   !<

type obstacle
  integer (kind=DPI) ::              numseg   !< indice du segment obstacle
  integer (kind=DPI) ::              react    !< nature de la reaction 1=interaction m SG 2=annihilation 3=jonction
                                              !< 4=pb ss reseau 5=obs deja jonction
  logical            ::              plani    !< key is true if obstacle belong to the i plan
  integer (kind=DPI),dimension(3) :: ptobs    !< Coordinates of the interaction point
                                              !! ptobs is negative if the interaction is a segment
end type obstacle

type (obstacle), allocatable :: OBSTAB(:)     !< tableau des obstacles

! If we want to do the DETECT_LOOP in Two Pass
logical                           :: TwoPass        !<
integer ,allocatable              :: MiniListeSeg(:)!<
integer                           :: MiniNsegTester !<
integer                           :: dim_alloc      !< Size of the allocatable obstacle lists

#ifdef PA
! Variables and arrays need in the parallel DETECT_LOOP calculation
integer ,allocatable    :: MiniListe_Send(:)!<
integer ,allocatable    :: MiniListe_Recv(:)!<
integer ,allocatable    :: MiniNseg_Send(:) !<
integer ,allocatable    :: MiniNseg_Recv(:) !<
integer                 :: array_deb        !<
integer                 :: array_fin        !<
#endif

contains

!######################################################################################
!> update is the procedure dealing with obstacles and displacement.                   \n
!! It is made of the following calculation steps:                                     \n
!! - elimination of the immobile segments (junction, gd, etc ..)                      \n
!! - listing of potential obstacles during displacement                               \n
!! - management of problems of connectivity along the line before moving              \n
!! - accounting of local rules which may affect the displacement:                     \n
!!    -> topological hardening (line recombination)                                   \n
!!    -> fixed boundaries                                                             \n
!!    -> particles                                                                    \n
!!    -> forest obstacles (list of obstacles, definition of obstacles and selection)  \n
!! - swept area calculation                                                           \n
!! - the effective displacement:                                                      \n
!!    -> displacement with recombination of segments                                  \n
!!    -> simple displacement                                                          \n
!! - Management of the reaction mechanism after displacement                          \n
!!    -> junctions                                                                    \n
!!    -> annihilation                                                                 \n
!! - segments are centered to account for the PBC                                     \n
!! - the management of the two lists par_touch() and par_cut() in creep mode          \n
!######################################################################################
subroutine update

implicit none

integer(kind=DPI)  :: NsegTester      !< The number of obstacle to test
integer(kind=DPI)  :: nbobs           !< The reduce number of obstacle found in the loop

integer(kind=DPI)  :: j           !<
integer(kind=DPI)  :: itemp       !<
integer(kind=DPI)  :: ktemp       !<
integer(kind=DPI)  :: i2e         !<
integer(kind=DPI)  :: i2o         !<
integer(kind=DPI)  :: Ivnn1o      !<
integer(kind=DPI)  :: Ii1o        !<
integer(kind=DPI)  :: VoisI1o     !<
integer(kind=DPI)  :: Ivnn1e      !<
integer(kind=DPI)  :: lomai       !<
integer(kind=DPI)  :: LGD         !<
integer(kind=DPI)  :: ie          !<
integer(kind=DPI)  :: io          !<
integer(kind=DPI)  :: oe          !<
integer(kind=DPI)  :: oo          !<
integer(kind=DPI)  :: Ii1e        !<
integer(kind=DPI)  :: VoisI1e     !<
integer(kind=DPI)  :: ovnn1o      !<
integer(kind=DPI)  :: ovnn1e      !<
integer(kind=DPI)  :: bufvoise    !<
integer(kind=DPI)  :: AbsDepMax   !<
integer(kind=DPI)  :: AbsDeptest  !<
integer(kind=DPI)  :: absdepini   !<
integer(kind=DPI)  :: absdeprail  !<
integer(kind=DPI)  :: SIGdep      !< Sign of the displacement direction
integer(kind=DPI)  :: ABSdep      !< Amplitude of the displacement
integer(kind=DPI)  :: LONGSEGj    !<
integer(kind=DPI)  :: BURGi       !<
integer(kind=DPI)  :: BJ          !<
integer(kind=DPI)  :: bufnorme    !<
integer(kind=DPI)  :: nvnormeI    !<
integer(kind=DPI)  :: nvnormeO    !<
integer(kind=DPI)  :: DELLONG     !<
integer(kind=DPI)  :: LiE         !<
integer(kind=DPI)  :: LONAD       !<
integer(kind=DPI)  :: LONVOAD     !<
integer(kind=DPI)  :: LONVEAD     !<
integer(kind=DPI)  :: londiff     !<
integer(kind=DPI)  :: LONADtest   !<
integer(kind=DPI)  :: LONVOADtest !<
integer(kind=DPI)  :: LONVEADtest !<
integer(kind=DPI)  :: bufveclin   !<
integer(kind=DPI)  :: I2          !<
integer(kind=DPI)  :: I1          !<
integer(kind=DPI)  :: codeobst    !<
integer(kind=DPI)  :: DomOold     !< Store I domain at Origin
integer(kind=DPI)  :: DomEold     !< Store I Domain at End
integer(kind=DPI)  :: Lji         !<
integer(kind=DPI)  :: Ljj         !<
integer(kind=DPI)  :: interi      !<
integer(kind=DPI)  :: interj      !<
integer(kind=DPI)  :: HANDIO      !<
integer(kind=DPI)  :: HANDIE      !<
integer(kind=DPI)  :: HANDOO      !<
integer(kind=DPI)  :: HANDOE      !<
integer(kind=DPI)  :: HANDIOold   !<
integer(kind=DPI)  :: HANDIEold   !<
integer(kind=DPI)  :: HANDOOold   !<
integer(kind=DPI)  :: HANDOEold   !<
integer(kind=DPI)  :: GDHANDIO    !< A var to remember if we pass a GD before reaching HandIo
integer(kind=DPI)  :: GDHANDIE    !< A var to remember if we pass a GD before reaching HandIe
integer(kind=DPI)  :: GDHANDOO    !< A var to remember if we pass a GD before reaching HandOo
integer(kind=DPI)  :: GDHANDOE    !< A var to remember if we pass a GD before reaching HandOe
integer(kind=DPI)  :: DVNNo       !<
integer(kind=DPI)  :: DVNNI       !<
integer(kind=DPI)  :: DVNNE       !<
integer(kind=DPI)  :: jonc_cote_i !<
integer(kind=DPI)  :: jonc_cote_j !<
integer(kind=DPI)  :: ntmp        !<
integer(kind=DPI)  :: itemp_loop  !<
integer(kind=DPI)  :: NdomO       !<
integer(kind=DPI)  :: NdomE       !<
integer(kind=DPI)  :: surfR       !< surface reaction type
integer(kind=DPI)  :: SFP         !< index of free surface where reaction takes place
integer(kind=DPI)  :: temp1       !<
integer(kind=DPI)  :: temp2       !<
integer(kind=DPI)  :: InnI        !<
integer(kind=DPI)  :: CInnI       !<
integer(kind=DPI)  :: InnJ        !<
integer(kind=DPI)  :: liobr       !<
integer(kind=DPI)  :: liebr       !<
integer(kind=DPI)  :: nloi        !<
integer(kind=DPI)  :: debut       !<
integer(kind=DPI)  :: fin         !<
integer(kind=DPI)  :: Lj          !<
integer(kind=DPI)  :: counter_o   !< Counter of segment between vnno and i
integer(kind=DPI)  :: counter_e   !< Counter of segment between i and vnne
integer(kind=DPI)  :: ik          !<

integer(kind=DPI) ,dimension(3) :: OiP      !<
integer(kind=DPI) ,dimension(3) :: centrei  !<
integer(kind=DPI) ,dimension(3) :: EiP      !<
integer(kind=DPI) ,dimension(3) :: Oobs     !<
integer(kind=DPI) ,dimension(3) :: Eobs     !<
integer(kind=DPI) ,dimension(3) :: ptobst   !<
integer(kind=DPI) ,dimension(3) :: NBDIM    !<
integer(kind=DPI) ,dimension(3) :: DECAL    !<
integer(kind=DPI) ,dimension(3) :: OiPtest  !<
integer(kind=DPI) ,dimension(3) :: EiPtest  !<
integer(kind=DPI) ,dimension(3) :: Oitest   !<
integer(kind=DPI) ,dimension(3) :: Eitest   !<

real(kind=DP)      :: temp        !<
real(kind=DP)      :: radius      !<
real(kind=DP)      :: mudivV      !<
real(kind=DP)      :: length3seg  !< length of dislocation made of 3 connected segments

real(kind=DP) ,dimension(3)      :: rvdvoiso  !<
real(kind=DP) ,dimension(3)      :: rvdvoise  !<
real(kind=DP) ,dimension(3)      :: VDEi      !<
real(kind=DP) ,dimension(3)      :: VDOi      !<

!< Variables for particule detection
integer(kind=DPI) ,dimension(3)     :: VDO,VDE,OIT,EiT,Ci             !<
integer(kind=DPI)                   :: selected_par                   !< The particule we keep at the end
logical                             :: pinned                         !<
real(kind=DP)                       :: normvdo,normvde,normsegi       !<

!< Variables for boundaries treatment
logical :: bloquer      !<
logical :: notannili2e  !<
logical :: notannili2o  !<
logical :: itemp_mark   !<
logical :: Ltemp        !<
logical :: JUMP         !<
logical :: condition    !<
logical :: JoncOI       !<
logical :: JoncEi       !<
logical :: JoncOJ       !<
logical :: JoncEJ       !<

logical ,dimension(2)  :: i1orotule !<
logical ,dimension(2)  :: i1erotule !<

! Variables for the particle detection and creep procedure
integer(kind=DPI) ,dimension(3) :: Oj   !<
integer(kind=DPI) ,dimension(3) :: VLj  !<
integer(kind=DPI) ,dimension(3) :: OiOj !<
integer(kind=DPI) ,dimension(3) :: EiOj !<

real(kind=DP)       :: distij       !<
real(kind=DP)       :: normnvec     !<

integer (kind=DPI)  :: absdep_par !<
integer (kind=DPI)  :: ipar_cut   !<
integer (kind=DPI)  :: ipar_touch !<
integer (kind=DPI)  :: icont2     !<
integer (kind=DPI)  :: iswitch    !<
integer (kind=DPI)  :: touch2cut  !<
integer (kind=DPI)  :: vo         !<
integer (kind=DPI)  :: ve         !<
integer (kind=DPI)  :: ll         !<

real(kind=DP)       :: normVDI        !<
real(kind=DP)       :: OiOjprojN      !<
real(kind=DP)       :: taueff         !<
real(kind=DP)       :: tauclimb       !<
real(kind=DP)       :: theta          !<

real(kind=DP), allocatable   :: time_climb(:)     !<
real(kind=DP), allocatable   :: vitesse_climb(:)  !<

! Variables for SELEC subroutine
integer (DPI) :: obsta        !<
integer (DPI) :: jtemp        !<
integer (DPI) :: voitmp       !<
integer (DPI) :: jselec       !<
integer (DPI) :: lvoitmp      !<
integer (DPI) :: Absdepobsta  !<
integer (DPI) :: Ljselec      !<
integer (DPI) :: nvoi         !<

logical       :: pb3b         !< is true if there is a 3 body problem
logical       :: kneecap      !< is true if the obstacle is
logical       :: jnct         !<
logical       :: anseg        !<
logical       :: condition1   !<
logical       :: Condition2   !<
logical       :: FirstTime    !< Logical parameter noting first and second call of DETECT_LOOP
logical       :: simpleDisp   !< True if simple displacement
logical       :: recomb       !< True if recombination
logical       :: gdobst       !<
logical       :: i_passed     !< used in a loop to know if we passed i

integer(kind=DPI)   :: itmp             !< Temporary variable for storing segment indice in case we have problems while connecting free surfaces
integer(kind=DPI)   :: thecounter       !< Counter for seg indices

integer(kind=DPI)   :: domchange        !< if we change from domain in barriere_conc
integer(kind=DPI)   :: maxdep           !< The largest integer displacement during the step

logical      :: surf_correct

! Variables for the modulo operations
integer(kind=DPI) ,dimension(3) :: Tr         !< The shift vector between i center and the simulation volume center
integer(kind=DPI)               :: Trx        !<
integer(kind=DPI)               :: Try        !<
integer(kind=DPI)               :: Trz        !<
integer(kind=DPI)               :: ax         !<
integer(kind=DPI)               :: ay         !<
integer(kind=DPI)               :: az         !<
integer(kind=DPI)               :: BI         !< The Greengard domain associated to the initial segment i (center)
integer(kind=DPI)               :: Bio        !< The Greengard domain associated to Oi
integer(kind=DPI)               :: Bie        !< The Greengard domain associated to Ei
integer(kind=DPI)               :: BiPo       !< The Greengard domain associated to OiP
integer(kind=DPI)               :: BiPe       !< The Greengard domain associated to EiP
real                            :: invmodurx  !<
real                            :: invmodury  !<
real                            :: invmodurz  !<
real                            :: modurx     !<
real                            :: modury     !<
real                            :: modurz     !<
logical                         :: allshift   !< All the segment must be shifted for this i segment

real(kind=DP)                   :: invtailleboite1, invtailleboite2, invtailleboite3
real(kind=DP)                   :: daireBV      !< The segments elementary swept area

!-----------------------------------

! Allocation of the obstacle list size. The default value is 100
dim_alloc = 100

#ifdef PA
! // simulations with many procs are simulation with many segments,
! hence the dim_alloc value must be increased
if (TAILLE_PAIRIMPAIR * 20 > dim_alloc) dim_alloc= TAILLE_PAIRIMPAIR * 20
allocate (OBSTAB(dim_alloc))
allocate (MiniListeSeg(dim_alloc))
allocate (MiniListe_Send(dim_alloc))
allocate (MiniListe_Recv(dim_alloc))
#else
allocate (OBSTAB(dim_alloc))
allocate (MiniListeSeg(dim_alloc))
#endif

#ifdef OUTPUTSURF
write(999,'(e23.15)') AVALUE
write(999,'(3i8)') MODUR(1), MODUR(2), MODUR(3)
#endif

if (kkdebug) write(379,fmt='("**Enter Update subroutine**")')

! Test to stop the simulation if Nsegm is too closed to Nsegmax
if ( int(NsegSafty*float(Nsegm),DPI) > nsegmax ) then
  write (*, fmt='("The number of segments has reached : Nsegm = ",I7,"Maximum number of segments allowed = ",I7,"-",F10.3)') &
  nsegm,nsegmax,NsegSafty
  stop
endif

! Initialization
mudivV                = xmu/VOLUME
PlusSeg               = izero       ! Number of added segments
dAIRE(1:Nsegm)        = zero        ! Area swept by segment
seg(1:nsegm)%NPhase   = izero
AireSysInst(1:NTSG)   = zero
AirevisSys(1:NTSG)    = zero
AireCoinSys(1:NTSG)   = zero
AireMixtSys(1:NTSG)   = zero
D_Vraie(1:NLV)        = zero        !< The average mobile segment displacement
LSEG_V(1:NLV)         = zero        !< The total length of mobile segments for each mobility law
TrAppSysInst(1:NTSG)  = zero
TrIntSysInst(1:NTSG)  = zero
NSEGMO                = izero       ! Mobile dislocation counter
jonc_cote_i           = izero
jonc_cote_j           = izero
HANDIOold             = izero
HANDIEold             = izero
HANDOOold             = izero
HANDOEold             = izero

invtailleboite1       = un/tailleBoite(1)
invtailleboite2       = un/tailleBoite(2)
invtailleboite3       = un/tailleBoite(3)
volboite              = tailleBoite(1) * tailleBoite(2) * tailleBoite(3)  ! Box volume (a3 unit) for the Greengard method
BdivAboite            = BdivA / volboite                                  ! Used to calculate gammabox
allshift              = .false.

if (CREEP) then
  par_cut(:,4) = 0
  par_touch(:,4) = 0
endif

#ifdef PA
  if (TAILLE_PAIRIMPAIR > 1) then
    TwoPass = .true.              ! If parallel, We want to do the DETECT_LOOP in two pass
    allocate(MiniNseg_Send(TAILLE_PAIRIMPAIR))
    allocate(MiniNseg_Recv(TAILLE_PAIRIMPAIR))
  else
    TwoPass = .false.             ! The default solution, We want to do the DETECT_LOOP in one pass
  endif
#else
  TwoPass = .false.               ! The default solution, We want to do the DETECT_LOOP in one pass
#endif

#ifdef OUTPUTSURF
write(999,'(e23.15)') AVALUE
write(999,'(3i8)') MODUR(1), MODUR(2), MODUR(3)
#endif

#ifdef MDC
#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
#endif
sweptsurfdata(1)=-1
#if defined(MDC) && defined(PA)
 endif
#endif

#endif

if (NbDep /= 0) then
  maxdep = abs(IDEP(QUIDEP(1)))
  allocate(DEP_REAC(5,maxdep))
endif

!###############################################################################
! The loop for the displacement of segments is not made on all the segments
! but only on segment with a predicted displacement and in a predefined order
! fixed with the tab QUIDEP
!###############################################################################
SURI : do ik = 1,NBDEP

  ! initialization
  gdobst  = .false.
  jnct    = .false.    !< Junction key
  anseg   = .false.    !< Annihilation key
  kneecap = .false.
  surf_correct = .false.

  ! The previous i segment displacement required an "all shift" operation, hence we must
  ! put again all the segments in the reference volume. Such operation is supposed to be very
  ! rare since mobile segments are supposed to be discretized.
  if(allshift .or. NBOITES < IDEUX) then
    seg(1:nsegm)%o(1)=modulo(seg(1:nsegm)%o(1),modur(1))
    seg(1:nsegm)%o(2)=modulo(seg(1:nsegm)%o(2),modur(2))
    seg(1:nsegm)%o(3)=modulo(seg(1:nsegm)%o(3),modur(3))
    allshift = .false.
  endif

  ! Segments are moved thanks to the list quidep which is a list of moving segments from fastest to slowest
  i = QUIDEP(ik)

  ! If during the step I was modified to a segment to eliminate, hence no need to try to move it
  if ((i == seg(i)%voiso .and. i == seg(i)%voise) .or. out(i)) cycle SURI

  ! I during the step was modified to be part of an elementary closed loop, hence no need to try to move it
  ! This test must exclude segment source with two pinning points and infinit line made of of two segments
  ! connected together with the PBC
  if ( seg(i)%vnno == seg(i)%vnne .and. seg(i)%vnno /= nsegmax .and.  &
      (seg(i)%veclin /= seg(seg(i)%vnno)%veclin)) then
    if (kkdebug) call seginfo(I,"A small closed loop eliminated starting from I")
    call desbou(i)   ! This loop is eliminated
    cycle
  endif

  if (kkdebug) then
    write(379,*)  "************************************************************************** "
    write(379,fmt='(" >>>>> We move segment = ",I7, ":   with IDEP = ",I7)') i,IDEP(i)
    write(379,*)  "************************************************************************** "
  endif

  ABSdepini = IDEP(i)
  ABSDEPmax = IZERO
  longsegi  = seg(i)%norme

  ! The segments which must not be moved : junction segment + cross-slip segment + kneecap segment
  ! and segments with no prediction of displacement
  if (seg(i)%jonc .or. seg(i)%gd > izero .or. absdepini == IZERO .or. LONGSEGi == IZERO) then
    if(kkdebug) write(379,*)  "> Cycle over I : Jonc or GD or zero idep or zero length"
    cycle SURI
  endif

  ! The segments which must not be moved : Segment quasi immobile and managed with the wait algorithm
  if (seg(i)%wait >= ITROIS .and. (.not.seg(i)%bloquer)) then
      if(kkdebug) write(379,*)  "> Cycle over I : wait > 3 and no bloquer"
      cycle suri
  endif

  ! As an effect of segment displacements during the step, strange configurations with many GD segments may be formed
  ! We try first to simplify such strange configurations
  if (seg(i)%vnno == Nsegmax) then
    itmp = seg(i)%voiso             ! The special case of pinning points
  else
    itmp = seg(i)%vnno
  endif

  ! Number of segments between the vnno and the vnne is calculated
  thecounter = izero      ! Total counter
  counter_o  = izero      ! Counter on the O side
  counter_e  = izero      ! Counter on the E side
  i_passed   = .false.
  do while (itmp /= seg(seg(i)%vnne)%voise)

    if (thecounter > 15) then
      print *, "A very very strange configuration is found for i at the beginning of update, kk= ",kk
      call seginfo(i,"Too many segments to be honest!")
      stop
    endif

    ! Total number of segments
    thecounter  = thecounter + iun

    ! Number of segments on the O and E sides is incremented
    if (itmp == i) then
      i_passed = .true.
      counter_o  = counter_o + iun
      counter_e  = counter_e + iun
    elseif(i_passed) then
      counter_e  = counter_e + iun
    else
      counter_o  = counter_o + iun
    endif

    itmp                = seg(itmp)%voise         !Update itmp

  enddo

  ! We can simplify the i2o i i2e connectivity
  if (thecounter > 9) then

    if (kkdebug) call seginfo(i,'Configuration of i before simplification')

    if (counter_o > 5) then
      vnn1o = seg(i)%vnno
      GD2made = .false.
      call ConnecIJ(vnn1o,i,72)
    endif
    if (counter_e > 5) then
      vnn1e = seg(i)%vnne
      GD2made = .false.
      call ConnecIJ(i,vnn1e,73)
    endif

    if (kkdebug) call seginfo(i,'i2o i i2e connectivity was simplified before displacement')
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Identification of the 1st and 2nd vnno and vvne of segment i
  vnn1o = seg(i)%vnno
  vnn2o = seg(vnn1o)%vnno
  vnn1e = seg(i)%vnne
  vnn2e = seg(vnn1e)%vnne

  ! A small loop of 3 connected segments is not supposed to move
  if (vnn2o == vnn1e .and. idix*LONGSEGi < loma(seg(i)%veclin)) then
    if(kkdebug) write(379,*)  "> Cycle over I : I is part of a small loop"
    cycle suri
  endif

  ! A small loop to eliminate latter made of two segments interconnected with two GD segments
  if (vnn2o == vnn2e .and. (seg(vnn1o)%gd * seg(vnn1e)%gd) > izero) then
    if(kkdebug) write(379,*)  "> Cycle over I : I is part of a small loop made of two GD segments"
    cycle suri
  endif


  if(kkdebug) call seginfo(I, "Segment information before displacement (SURI Loop) ")

  ! The Number of segment is now increased to Nsegm + Pluseg
  dAIRE(1+NSEGM:NSEGM+PLUSSEG) = zero           ! The swept area for the new segments is set to zero
  Nsegm   = Nsegm + PlusSeg                     ! The new Nsegm value
  PlusSeg = izero
  obsta   = izero

  ! Initialization of possible new segments
  i1 = Nsegm  + 1
  i2 = Nsegm  + 20
  SEG(i1 : i2)%JONC         = .FALSE.
  SEG(I1 : I2)%IJONC        = nsegmax
  SEG(I1 : I2)%tJONC        = 0
  SEG(I1 : I2)%GD           = izero
  SEG(I1 : I2)%WAIT         = izero
  SEG(I1 : I2)%bloquer      = .false.
  SEG(I1 : I2)%unload       = .false.
  SEG(I1 : I2)%NORME        = izero
  SEG(I1 : I2)%probadev     = zero
  SEG(I1 : I2)%surface      = izero
  SEG(I1 : I2)%VarFreePlan  = izero
  SEG(I1 : I2)%dom          = izero
  SEG(I1 : I2)%anglevis     = -un

  ! test on the new Nsegm value
  if(nsegm > nsegmax) then
    write(*,*) '! In SURI loop, the number of segments exceed Nsegmax !'
    write(*,*) "Number of segments :",nsegm, "Number of segment to add :",PlusSeg
    stop
  endif

  ! Local variables to speed up computations
  ABSdep  = int(abs(absdepini),DPI)
  SIGdep  = int(SIGN(IUN,absdepini),dpi)
  Oi(:)   = SEG(i)%O(:)
  Li      = SEG(i)%VECLIN
  VLi(:)  = Bveclin(:,Li)
  VDi(:)  = Bvecdep(:,Li)
  nloi    = LoiSeg(I)
  InnI    = Inn(Li)
  CInnI   = bveclin(InnI,Li)
  LomaI   = loma(Li)
  BI      = IBoite(I)           ! The multipole domain associated to the segment i

  ! Topological quantities
  ! A minimum length for discretization is defined. Such length is important for some particular
  ! segments, like segment with a pinning point at one end.
  limite_discreti = int(float(Lomai) / 10.0,DPI)
  if(particules) limite_discreti = int(Lomai / trois,DPI)
  limite_bloquage = limite_discreti ! limite de taille de segment pour le bloquage a cause d'interaction impossible
                                    ! (sous reseau, ..)

  ! The new center for the PBC shift is defined
  Centrei(1:3) = Oi(1:3)+LONGSEGi/2*VLi(1:3)
  Centrei(1:3) = modulo(Centrei(1:3), modur(1:3)) ! shift procedure restricted to border boxes
                                                  ! requires modulo operations on i
  Tr(:) = modur(1:3)/2 - Centrei(1:3) ! vecteur de traslation applique sur tous les segments
  ! intermediate variables for optimization
  Trx = modur(1)/2 - Centrei(1) ! The shift vector from segment i center and the simulated volume center
  Try = modur(2)/2 - Centrei(2)
  Trz = modur(3)/2 - Centrei(3)
  modurx = real(modur(1))
  modury = real(modur(2))
  modurz = real(modur(3))
  invmodurx = 1./real(modur(1))
  invmodury = 1./real(modur(2))
  invmodurz = 1./real(modur(3))

  !********************************************************************************
  ! The list of segments to be tested for reactions with I is made (listesegtester)
  !********************************************************************************

  NSegTester  = izero
  NbObs       = izero

  ! Use is made of the multipole domains to reduce the number of segment to consider
  if (NBOITES > IUN) then

    ! In particular case when i is a very long seg, every segment inside the volume is a possible obstacle
    ! This case should exist only during the first steps of the simulation, before discretization
    if (longseg(i)) then

      ! for long seg, the following procedure is deliberately not optimized
      ! since long segments are not supposed to move
      if(kkdebug) write(379,*)  ">>> Modulo shift is applied to border LONG segments"

      ! In this particular case all j segment must be shifted latter in detect_loop
      allshift = .true.

      Nsegtester = izero
      do ktemp = 1, nsegm
        if (longseg(ktemp) .or. seg(ktemp)%surface  > IZERO) cycle     ! The longseg are added latter
        Nsegtester = Nsegtester + iun
        listesegtester(Nsegtester) = ktemp
      enddo

    else

      ntmp = izero

      do ktemp = iun, 27 ! Only the segments in the BI domain and the close domains are considered

        BJ = IndBoiteVois(BI, ktemp)

        ! Segments in the BJ domains must be considered as possible obstacles
        if (NSegBoite(BJ) /= Izero) then

          debut                 = ntmp + iun
          fin                   = ntmp + NsegBoite(BJ)
          listetmp(debut : fin) = IndexBoite(BJ)%ListeSeg(1:NsegBoite(BJ))
          ntmp                  = ntmp + NsegBoite(BJ)

        endif

      enddo

      ! The list of segments to test as possible obstacles is saved
      ! Note that the list of longseg is going to be added latter
      Nsegtester = izero
      do ktemp = 1, ntmp
        if (longseg(listetmp(ktemp)) .or. seg(ktemp)%surface > IZERO) cycle     ! The longseg are added latter
        Nsegtester = Nsegtester + iun
        listesegtester(nsegtester) = listetmp(ktemp)
      enddo


    endif

    ! The long segments, not yet taken into account, are systematically added
    do itemp = 1, Nb_groseg
      nsegtester = nsegtester + iun                         ! the final number of segments to test
      listesegtester(nsegtester) = GROSEG(itemp)
    enddo

  ! When no multipole domains are defined, we must test possible contact with all the segments
  else

    Nsegtester = Nsegm + plusseg
    do itemp = 1, Nsegtester
      listesegtester(itemp) = itemp
    enddo

  endif

  ! Some useful information to remember
  seg(i)%o(:) = modulo(seg(i)%o(:) + Tr(:),modur(:)) - Tr(:)
  Oi(:)         = seg(i)%O(:)
  Centrei(1:3)  = oi(1:3) + (LONGSEGi/2 * VLi(1:3))

  !**********************************************************************
  ! Beginning of the geometrical tests defining the connected segment to i
  ! and how they affect the displacement of i
  !**********************************************************************
  ! The connectivity of segment can be reconstruct before we move the segment I in special cases
  ! We test such reconstruction on the o and e side

  ! The screw segment direction associate to segment i is identified
  itemp  = assoc(Li,1)
  VLj(:) = bveclin(:,itemp)

  !***************
  !* 1st neighbor at the origin of segment i
  !***************

  LGD = seg(vnn1o)%veclin
  ! Problem to move a segment if it is close to two GD segments
  if (seg(vnn1o)%GD > izero) then
    ! a particular case: if seg(vnn2o)%GD > izero then connecij cannot be applied
    ! and the only possibility is to stop the segment displacement at this step
    if (seg(vnn2o)%GD > izero) then

      cycle suri

    else

      ! The vnn1o is a GD with non zero length, better to force this GD to be on the same slip system as i
      if(seg(vnn1o)%norme /= izero .and. sysconnec(LGD,Li) == 2) then

        seg(vnn1o)%veclin = SegDevIJ(LGD,Li)

      ! The vnn1o is a GD with zero length
      elseif (seg(vnn1o)%norme == 0) then

        ! if i is a screw segment, the GD segment is switched to the same slip system as i
        if(Tyseg(Li) == IUN) then

          seg(vnn1o)%veclin = Li

        ! reste une ambiguite : quand le GD est une rotule mal orientee
        ! pour le deplacement de I (composante negative de VL(GD) sur la direction de deplacement)
        else
          if (sigdep * dot_product(VDI,VLJ) > 0) then
            seg(vnn1o)%veclin = ITEMP
          else
            seg(vnn1o)%veclin = assoc(Li,5)
          endif
        endif
      endif

      ! si on a change le caractere de vnn1o, il faut refaire la connection
      if (seg(vnn1o)%veclin /= LGD) then
        if (kkdebug)  write(379,*)  ">>> Vnn1o vectlin had been changed - reconnecting segments"
        GD2made = .false.
        call ConnecIJ(vnn2o,vnn1o,20)
        call ConnecIJ(vnn1o,i,20340)
        if (kkdebug)  call check_connec("Vnn1o character (Veclin) had been changed")
      endif
    endif
  endif

  !***************
  !* 1er neighbor at the end of segment i
  !***************

  LGD = seg(vnn1e)%veclin
  ! Problem to move a segment if it is close to two GD segments
  if (seg(vnn1e)%GD > izero) then
    ! a particular case: if seg(vnn2e)%GD > izero then connecij cannot be applied
    ! and the only possibility is to stop the segment displacement at this step
    if (seg(vnn2e)%GD > izero) then
      cycle suri
    else

      ! si la norme du gd /= 0, la seule chose qu'on peut faire est de le remettre dans le
      ! sys de gliss. de I
      if(seg(vnn1e)%norme /= izero .and. sysconnec(LGD,Li) == 2) then
        seg(vnn1e)%veclin = SegDevIJ(LGD,Li)

        ! le GD est une rotule !
      elseif (seg(vnn1e)%norme == 0) then

        ! si i est une vis on change le caractere de la rotule GD eviter pb de connection
        if(Tyseg(Li) == IUN) then
            seg(vnn1e)%veclin = Li

            ! reste une ambiguite : quand le GD est une rotule mal orientee
            ! pour le deplacement de I (composante negative de VL(GD) sur la direction de deplacement)
        else
            if (sigdep * dot_product(VDI,VLJ) < 0) then
              seg(vnn1e)%veclin = ITEMP
            else
              seg(vnn1e)%veclin = assoc(Li,5)
            endif
        endif
      endif

      ! si on a change le caractere de vnn1e, il faut refaire la connection
      if (seg(vnn1e)%veclin /= LGD) then
        if (kkdebug)  write(379,*)  ">>> Vnn1e vectlin had been changed - reconnecting segments"
        GD2made = .false.
        call ConnecIJ(I,vnn1e,21)
        call ConnecIJ(vnn1e,vnn2e,21340)
        if (kkdebug)   call check_connec("Vnn1e character (Veclin) had been changed")
        !         if (kkdebug)   call seginfo(I," apres chaggment de character en E      ")
        !             call ConnecIJ(vnn1e,vnn2e,2165)
      endif
    endif
  endif


  ! The new list of IVNNO, IVNNE, VNNE and VNNO segments is defined
  vnn1o = SEG(i)%voiso     !*** Premier voisin non nul
  Ivnn1o = IUN
  do while(vnn1o /= nsegmax .and. seg(vnn1o)%norme ==0 .and. seg(vnn1o)%gd < iun .and. .not. seg(vnn1o)%jonc)
      Ivnn1o = Ivnn1o + IUN
      vnn1o = seg(vnn1o)%voiso
  enddo
  seg(i)%vnno = vnn1o

  vnn2o = SEG(vnn1o)%voiso     !*** 2eme voisin non nul
  do while(vnn2o /= nsegmax .and. seg(vnn2o)%norme ==0 .and. seg(vnn2o)%gd < iun .and. .not. seg(vnn2o)%jonc)
      vnn2o = seg(vnn2o)%voiso
  enddo
  seg(vnn1o)%vnno = vnn2o

  vnn1e = SEG(i)%voise     !*** Premier voisin non nul
  Ivnn1e = IUN
  do while(vnn1e /= nsegmax .and. seg(vnn1e)%norme ==0 .and. seg(vnn1e)%gd < iun .and. .not. seg(vnn1e)%jonc)
      Ivnn1e = Ivnn1e + IUN
      vnn1e = seg(vnn1e)%voise
  enddo
  seg(i)%vnne = vnn1e

  vnn2e = SEG(vnn1e)%voise     !*** Premier voisin non nul
  do while(vnn2e /= nsegmax .and. seg(vnn2e)%norme ==0 .and. seg(vnn2e)%gd < iun .and. .not. seg(vnn2e)%jonc)
      vnn2e = seg(vnn2e)%voise
  enddo
  seg(vnn1e)%vnne = vnn2e

  if (kkdebug) then
    write (379,fmt='("--Segment Neighboring--")')
    write (379,fmt='("Segment :          ",I7," vnno : ",I7," vnne : ",I7)') i,seg(i)%vnno,seg(i)%vnne
    write (379,fmt='("Segment at origin: ",I7," vnno : ",I7)') vnn1o,seg(vnn1o)%vnno
    write (379,fmt='("Segment at end:    ",I7," vnne : ",I7)') vnn1e,seg(vnn1e)%vnne
  endif


  ! traitement des problem de mobilite autour des rotule jonction prevu pour gd
  ! This test is not clear to me (BD) and should be further commented (see G. Monnet)

  i1o = SEG(i)%VOISO      !*** Premiers voisins servant de rail lors du
  !    deplacement initialise comme %vois i.e.
  !    1er voisin nul ou pas
  Ii1o=1      !*** Positions de i1o
  VoisI1O = SEG(i1o)%VOISO !*** 1er voisin nul ou pas de i1
  ! la condition pour laquelle on verifie la courbure

  if (voisi1o /= NSEGMAX) then
    condition = (seg(i1o)%norme == 0                  .and. &
                .not. seg(i1o)%jonc                   .and. &
                seg(i1o)%gd < iun                     .and. &
                (Li == SEG(VoisI1o)%veclin .or.             &
                 Li == SegDevIJ(seg(voisi1o)%veclin,Li))    &
                )
  else
    condition = (seg(i1o)%norme == 0                  .and. &
                .not. seg(i1o)%jonc                   .and. &
                seg(i1o)%gd < iun)
  endif

  if (condition) then
      if (kkdebug) write(379,*) "> I1o veclin change from ",SEG(i1o)%veclin
      SEG(i1o)%veclin = courbuO(Li,((sigdep+1)/2)+1)
      Lio=seg(i1o)%veclin
      seg(i1o)%probadev = zero                      ! Indeed, We change the segment character
      if (kkdebug) write(379,*)  " to ",Lio
      IDEP(i1o) = IZERO
      SEG(i1o)%RESDEP = ZERO
  else
      LiO = SEG(i1o)%VECLIN
  endif


  i1e = SEG(i)%VOISE
  Ii1e=1
  vnn1e = SEG(i)%VNNE
  VoisI1E = SEG(i1e)%VOISE

  If (voisi1e /= Nsegmax) then
    condition = (seg(i1e)%norme == 0                    .and. &
                .not. seg(i1e)%jonc                     .and. &
                seg(i1e)%gd < iun                       .and. &
                (Li == SEG(VoisI1e)%veclin .or.               &
                 Li == SegDevIJ(seg(voisi1e)%veclin,Li))      &
                )
  else
    condition = (seg(i1e)%norme == 0                    .and. &
                .not.SEG(i1e)%JONC                      .and. &
                seg(i1e)%gd < iun)
  endif

  if (condition) then
      if (kkdebug) write(379,*) "> I1e veclin change from ",SEG(i1e)%veclin
      SEG(i1e)%veclin = courbuE(Li,((sigdep+1)/2)+1)
      Lie=SEG(i1e)%veclin
      seg(i1e)%probadev = zero                      ! Indeed, We change the segment character
      if (kkdebug) write(379,*) " to ",LiE
      IDEP(i1e) = IZERO
      SEG(i1e)%RESDEP = ZERO
  else
      LiE = SEG(i1e)%VECLIN
  endif

  !******************
  !* 2nd neighbors
  !******************
  i2o=seg(i1o)%vnno
  i2e=seg(i1e)%vnne
  LiOBR = deassoc(LiO)
  LiEBR = deassoc(LiE)


  ! in some case (defined in the condition variable) when a segment touches the surface has length ==1
  ! the kneecap  is changed to handle the predicted displacement.
  if (seg(i)%surface > IZERO .and. longsegi==IUN) then

    if (seg(i)%surface == iun) then
      do j=1,nbasered
        condition = ((sigdep > IZERO .and. fac2dep(Li,j,LiEBR) == -IUN .and. &
                     fac1dep(Li,LiEBR) == -IUN .and. fac1dep(Li,j) == IUN) .or. &
                     (sigdep < IZERO .and. fac2dep(Li,j,LiEBR) == IUN .and. &
                     fac1dep(Li,LiEBR) == IUN .and. fac1dep(Li,j) == -IUN) .and. &
                     moddepla(Li,j,LiEBR) == IUN)
        if (condition) then
          LiOBR=j
          seg(i1o)%veclin = assoc(Li,j)
          LiO = seg(i1o)%veclin
          if (kkdebug) write(379,*) "Liobr changed to ,",Lio," to move a surface segment"
        endif
      enddo
    endif

    if(seg(i)%surface == Ideux) then
      do j=1,nbasered
        condition = ((sigdep > IZERO .and. fac2dep(Li,LiOBR,j) == -IUN .and. &
                     fac1dep(Li,LiOBR) == IUN .and. fac1dep(Li,j)== -IUN) .or. &
                     (sigdep < IZERO .and. fac2dep(Li,LiOBR,j) == IUN .and. &
                     fac1dep(Li,LiOBR) == -IUN .and. fac1dep(Li,j)== IUN) .and. &
                     moddepla(Li,LiObr,j) == IUN)
        if (condition) then
          LiEBR = j
          seg(i1e)%veclin = assoc(Li,j)
          LiE = seg(i1e)%veclin
          if (kkdebug) write(379,*) "LiEbr changed to ,",LiE," to move a surface segment"
        endif
      enddo
    endif

  endif
  !!!!
  ! some special cases where the simple displacement of i cannot easily proceed as a result of segment connectivity
  ! The caracter of the neighboring segment is incorrect
  if (modulo(LiOBR-Li,iquatre) == izero .or. modulo(LiEBR-Li,iquatre) == izero )then

      if (kkdebug) write(379,*) "! #1 Incompatibles rails after Vectlin modifications - Cycle (SURI loop) !"
      CYCLE surI

  endif

  lonad   = LONGSEGi+ nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep)
  lonvoad = SEG(i1o)%NORME+ nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep)
  lonvead = SEG(i1e)%NORME- nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep)

  if(kkdebug) then
    write(379,*) " **********"
    write (379,fmt='("Absdep : ",I7,"| SIGdep ",I7,"| modep",I7)') ABSdep,Sigdep,modep
    write (379,fmt='("-------------------------------")')
    write (379,fmt='("Length O-Neighbor After Displ. : ",I7," | fac1dep(Li,LiOBR)       : ",F7.3)') lonvoad,fac1dep(Li,LiOBR)
    write (379,fmt='("Length E-Neighbor After Displ. : ",I7," | fac1dep(Li,LiEBR)       : ",F7.3)') lonvead,fac1dep(Li,LiEBR)
    write (379,fmt='("Length Segment I  After Displ. : ",I7," | fac2dep(Li,LiOBR,LiEBR) : ",F7.3)') lonad,fac2dep(Li,LiOBR,LiEBR)
  endif

  ! The segment reducing the displacement amplitude is shift from v1o to vnn1o
  if(lonvoad < IZERO .and. SEG(i1o)%NORME == IZERO) then
      !* Non : deplacement special a l origine *
      if (vnn1o /= nsegmax .and. seg(vnn1o)%norme /= 0 .and. (seg(i1o)%gd < iun .or. seg(vnn1o)%veclin == Lio)) then
        !*** Deplacement des rotules avec le segment
        i1o = vnn1o
        LiO = SEG(i1o)%VECLIN
        LiOBR = deassoc(LiO)
        Ii1o = Ivnn1o
        VoisI1O = SEG(i1o)%VOISO
        lonad   = LONGSEGi+ nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep)
        lonvoad = SEG(i1o)%NORME+ nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep)
        lonvead = SEG(i1e)%NORME- nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep)
      endif
  endif

  ! The segment reducing the displacement amplitude is shift from v1e to vnn1e
  if ((lonvead < IZERO.and.SEG(i1e)%NORME == IZERO)) then
      !*** Meme chose mais en E
      !* Non : deplacement special a l extremite*
      if (vnn1e /= nsegmax .and. seg(vnn1e)%norme /= 0 .and. (seg(i1e)%gd < iun .or. seg(vnn1e)%veclin == Lie)) then
        i1e = vnn1e
        LiE = SEG(i1e)%VECLIN
        LiEBR = deassoc(LiE)
        Ii1e = Ivnn1e
        VoisI1E = SEG(i1e)%VOISE
        lonad   = LONGSEGi+ nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep)
        lonvoad = SEG(i1o)%NORME+ nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep)
        lonvead = SEG(i1e)%NORME- nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep)
      endif
  endif


  if (modulo(LiOBR-Li,iquatre) == izero .or. modulo(LiEBR-Li,iquatre) == izero )then

      if (kkdebug) write(379,*) "! #2 Incompatibles rails after Vectlin modifications - Cycle (SURI loop) !"

      CYCLE surI

  endif


  ! initialization of intermediate variable needed to move the segment
  VLi(1:3) = BVECLIN(1:3,Li)          !*** Vecteur ligne "signe"
  VDi(1:3) = BVECDEP(1:3,Li)*SIGDEP   !*** Vect deplacement avec signe theo
  !*** Vecteur deplacement "signe" des points Oi et Ei
  VDOi(1:3) = fac1dep(Li,LiOBR)*SIGDEP*real(BVECLIN(1:3,LiO),DP) !***
  VDEi(1:3) = fac1dep(Li,LiEBR)*SIGDEP*real(BVECLIN(1:3,LiE),DP) !*GD
  !*** Vecteur deplacement oriente vers l'interieur de la surface de i1o et i1e
  rVDvoiso = - fac1dep(LiO,deassoc(Li))*BVECDEP(1:3,ASSOC(Li,LiOBR))
  rVDvoise = fac1dep(LiE,deassoc(Li))*BVECDEP(1:3,ASSOC(Li,LiEBR))
  !* Vecteur normal au plan de glissement (signe donne par le vecteur ligne)
  Nvec(1:3) = BVECNOR(1:3,Li)
  !* Vecteur de BURGERS
  BURGi = ASSOC(Li,1)

  !*** Initialisations pour la mise a jour des voisins non nuls *************
  if (SEG(i1o)%NORME == 0 .and. seg(i1o)%gd < iun) then !*** 0 : longueur initialement nulle
    DVNNO=0      !*   1 : longueur initialement non nulle
  else
    DVNNO=1
  endif
  DVNNI=1
  if (SEG(i1e)%NORME == 0 .and. seg(i1e)%gd < iun) then !*** 0 : longueur initialement nulle
    DVNNE=0      !*   1 : longueur initialement non nulle
  else
    DVNNE=1
  endif

  ! The end of the segment i
  Ei(1:3) = Oi(1:3)+LONGSEGi*VLi(1:3)

  !******************************************************************************************
  !* End of the part defining the geometrical information needed before we move the segment i
  !******************************************************************************************

  !************************************************************************
  ! Topological hardening: the displacement of some segments is going to be
  ! reduce to account for the line connectivity topology. Better to reduce
  ! the amplitude of displacement to account for those problems before we
  ! surch for obtacles in the glide area.
  !************************************************************************
  ! L'idee est que c est inutile de s'amuser a chercher des obstacle
  ! en dehors du trapeze que l'on va balayer...

  ! ghiath : afin de generaliser la procedure mixt/mixt, on fait le test sur les segmnt i
  ! et ces voisins. L'idee est de garantir un longueur entier des segments en modulant
  ! le deplacement : absdep doit etre multiple de MODEP, qui est une variable integer
  ! calculee en fonction de caractersitiques geometriques des segmnts, ex cas des mixt/mixt
  ! dans les cfc MODEP=2

  ! initialisation
  temp1 = 1
  temp2 = 1
  modep = moddepla (Li,LiOBR,LiEBR)
  if (modep > modep_max ) then
      print *, "I=",I,"    Li =", Li
      call seginfo(i," Mauvais choix des railles : probleme de modep  ")
      call conf(i)
      print *, "I1O = ", I1O,"       LiOBR =", LiOBR
      !       call conf (i1o) ;      call conf (i1e)
      print *, "I1E =", I1E, "       LiEBR =", LiEBR
      print *, " fac1dep(Li,LiOBR)) =", fac1dep(Li,LiOBR)," fac1dep(Li,LiEBR)) =", fac1dep(Li,LiEBR)
      print *," tempi et modep", temp1, temp2, modep
      cycle SURI
  endif

  ABSdep = (ABSdep/modep)*modep !*** L'histoire du modulo 2 mixte?mixte

  if (abs(seg(I)%resdep) > 2 * modep_max) then
      print *, "depassement en resdep =", IDEP(I), sigdep * absdep
  endif

  lonad = -1
  ABSdep = ABSdep + modep

 !   if(kkdebug) then
 !     write(379,*) " **** Before Displacement"
 !     write (379,fmt='("Absdep : ",I7,"| SIGdep ",I7)') ABSdep,Sigdep
 !     write (379,fmt='("-------------------------------")')
 !     write (379,fmt='("Length O-Neighbor After Displ. : ",I7," | fac1dep(Li,LiOBR)       : ",F7.3)') lonvoad,fac1dep(Li,LiOBR)
 !     write (379,fmt='("Length E-Neighbor After Displ. : ",I7," | fac1dep(Li,LiEBR)       : ",F7.3)') lonvead,fac1dep(Li,LiEBR)
 !     write (379,fmt='("Length Segment I  After Displ. : ",I7," | fac2dep(Li,LiOBR,LiEBR) : ",F7.3)') &
 !           & lonad,fac2dep(Li,LiOBR,LiEBR)
 !   endif

  do while((lonad < IZERO .or. lonvoad < ZERO .or. lonvead < ZERO) .and. absdep >= 0)
      ABSdep = ABSdep - modep
      lonad   = LONGSEGi+nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep)
      lonvoad = SEG(i1o)%NORME+nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep)
      lonvead = SEG(i1e)%NORME-nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep)
  enddo

  if(kkdebug) then
    write(379,*) " **** After Displacement"
    write (379,fmt='("Absdep : ",I7,"| SIGdep ",I7)') ABSdep,Sigdep
    write (379,fmt='("-------------------------------")')
    write (379,fmt='("I1O :             ",I7,"| I : ",I7 ," | I1E : ",I7)') i1o,i,i1e
    write (379,fmt='("VFP :             ",I7,"| I : ",I7 ," | I1E : ",I7)') &
          & seg(i1o)%VarFreePlan,seg(i)%VarFreePlan,seg(i1e)%VarFreePlan
    write (379,fmt='("Surf connection : ",I7,"| I : ",I7 ," | I1E : ",I7)') seg(i1o)%surface,seg(i)%surface,seg(i1e)%surface
    write (379,fmt='("Length O-Neighbor After Displ. : ",I7," | fac1dep(Li,LiOBR)       : ",F7.3)') lonvoad,fac1dep(Li,LiOBR)
    write (379,fmt='("Length E-Neighbor After Displ. : ",I7," | fac1dep(Li,LiEBR)       : ",F7.3)') lonvead,fac1dep(Li,LiEBR)
    write (379,fmt='("Length Segment I  After Displ. : ",I7," | fac2dep(Li,LiOBR,LiEBR) : ",F7.3)') &
          & lonad,fac2dep(Li,LiOBR,LiEBR)
  endif

  condition1 = .false.
  if (lonvoad==IZERO .and. seg(i1o)%surface /= IZERO) then
    if (VLI(1)*Plane_MillerI(1,seg(i1o)%VarFreePlan)               &
        +VLI(2)*Plane_MillerI(2,seg(i1o)%VarFreePlan)              &
        +VLI(3)*Plane_MillerI(3,seg(i1o)%VarFreePlan) >= IZERO) condition1 = .true.
  endif
  condition2 = .false.
  if (lonvead==IZERO .and. seg(i1e)%surface /= IZERO) then
    if (VLI(1)*Plane_MillerI(1,seg(i1e)%VarFreePlan)              &
        +VLI(2)*Plane_MillerI(2,seg(i1e)%VarFreePlan)               &
        +VLI(3)*Plane_MillerI(3,seg(i1e)%VarFreePlan) <= IZERO) condition2 = .true.
  endif

  if (ABSdep /= IZERO .and. seg(i)%surface < IUN .and. (condition1 .or. condition2)) then

    do while((lonvead == IZERO .or. lonvoad == IZERO) .and. absdep > IZERO)
      ABSdep=ABSdep-modep! ABSdep=ABSdep-IUN
      if (ABSdep<=IZERO) Absdep=IZERO
      lonad   = LONGSEGi+ nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep)
      lonvoad = SEG(i1o)%NORME+ nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep)
      lonvead = SEG(i1e)%NORME- nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep)

      if(kkdebug) then
        write(379,*) " **********"
        write (379,fmt='("Absdep : ",I7,"| SIGdep ",I7)') ABSdep,Sigdep
        write (379,fmt='("-------------------------------")')
        write (379,fmt='("Length O-Neighbor After Displ. : ",I7," | fac1dep(Li,LiOBR)       : ",F7.3)') lonvoad,fac1dep(Li,LiOBR)
        write (379,fmt='("Length E-Neighbor After Displ. : ",I7," | fac1dep(Li,LiEBR)       : ",F7.3)') lonvead,fac1dep(Li,LiEBR)
        write (379,fmt='("Length Segment I  After Displ. : ",I7," | fac2dep(Li,LiOBR,LiEBR) : ",F7.3)') &
              & lonad,fac2dep(Li,LiOBR,LiEBR)
      endif

    enddo

    if (kkdebug) write(379,*) "Displacement reduced because of free surface incompatibility"
    seg(i)%diseg=.true.

  endif

  ! mettre absdep en conformite avec les railles : multiple de modep
  if(kkdebug) then
    write(379,*) " **** After Displacement"
    write (379,fmt='("Absdep : ",I7," | SIGdep ",I7)') ABSdep,Sigdep
    write (379,fmt='("AbsdepINI : ",I7," | Modep ",I7," | After rails ",I7)') absdepini, modep, absdep*sigdep
    write (379,fmt='("-------------------------------")')
    write (379,fmt='("I1O :             ",I7," | I : ",I7 ," | I1E : ",I7)') i1o,i,i1e
    write (379,fmt='("VFP :             ",I7," | I : ",I7 ," | I1E : ",I7)') &
          & seg(i1o)%VarFreePlan,seg(i)%VarFreePlan,seg(i1e)%VarFreePlan
    write (379,fmt='("Surf connection : ",I7," | I : ",I7 ," | I1E : ",I7)') seg(i1o)%surface,seg(i)%surface,seg(i1e)%surface
    write (379,fmt='("Length O-Neighbor After Displ. : ",I7," | fac1dep(Li,LiOBR)       : ",F7.3)') lonvoad,fac1dep(Li,LiOBR)
    write (379,fmt='("Length E-Neighbor After Displ. : ",I7," | fac1dep(Li,LiEBR)       : ",F7.3)') lonvead,fac1dep(Li,LiEBR)
    write (379,fmt='("Length Segment I  After Displ. : ",I7," | fac2dep(Li,LiOBR,LiEBR) : ",F7.3)') &
          & lonad,fac2dep(Li,LiOBR,LiEBR)
  endif

  ! modification de absdep a cause de ces voisins - la derniers avant une autre modification
  ! a cause des obstacle
  !    l'amplitude restant a disposition apres le
  !    mvt a venir

  ! The segment could not proceed with all its integer displacement during this step
  ! We keep the information about the displacement we lost in resdep
  seg(I)%resdep =  seg(I)%resdep  + (IDEP(I) - sigdep * absdep)
  ! It is useless to cumulate more than modep_max resdep. If the segment did not proceed this
  ! displacement it is topologically bloc by its rails.
  if (abs(seg(I)%resdep) > modep_max) then
    seg(I)%resdep = modep_max *sign(un,seg(I)%resdep)
    ! print*,'resdep correction is applyed',kk,i,modep_max,seg(I)%resdep
  endif

  !  Si absdep est negative c'est que le mouvement du segment est impossible
  !  compte tenue de ses railles
  if (absdep <= 0) then
      seg(i)%diseg = (longsegI > Limite_discreti .and. abs(idep(I)) > modep_max)
      ! remettre absdep a l'echelle apres reduction eventuel de deplacement
      ! a acuse des obstacles
      if (kkdebug) then
        write (379,*) " -----------------------------------------------------------------"
        write (379,fmt='( "Exiting SURI loop for segment =",I7," because of absdep = 0 after taking rail into account ")') i
        write (379,*) " -----------------------------------------------------------------"
        call seginfo(i,"Exiting SURI loop")
      endif
      cycle surI
  endif

  !*** Que fait-on du reste d'amplitude ?
  !*** INUTILE DE PERDRE SON TEMPS SI BLOCAGE ON
  !*** VERRA A l'ITERATION SUIVANTE (si on bouge pas
  !*** il n'y a aucune chance de changer la topologie)

  !*** DEBOGUONS
  !*** Si encore un pb c'est que la configuration n'etait pas propre !!!
  OiP(1:3) = Oi(1:3)+nint((ABSdep*VDOi(1:3)),DPI)
  EiP(1:3) = Ei(1:3)+nint((ABSdep*VDEi(1:3)),DPI) ! nouveaux fait par /Ghiath 23/05/01
  Oj(:) = BVECDEP(:,Lio)! vecteur deplacement du rail en Oi
  VDi1o_or(:) = 999999
  if (dot_product(VLi,Oj) > izero)then
      VDi1o_or(:) = Oj(:) !VD du rail en Oi oriente vers l'aire balaye par i
  else
      VDi1o_or(:) = -Oj(:)
  endif
  ! recherche du VD du rail en e oriente vers l'aire balayee de i
  Oj(:) = BVECDEP(:,Lie)! vecteur deplacement du rail en Ei
  VDi1e_or(:) = 999999
  if (dot_product(VLi,Oj) > izero) then
      VDi1e_or(:) = -Oj(:) !VD du rail en Oi oriente vers l'aire balaye par i
  elseif (dot_product(VLi,Oj) < izero) then
      VDi1e_or(:) = Oj(:)
  endif
  if (kkdebug) then
      write(379,*) " OiP(1:3) =", OiP(1:3)
      write(379,*) "  EiP(1:3) =",  EiP(1:3)
  endif

  !************************************************************************
  ! End of the topological hardening test: the displacement of segments is
  ! is now compatible with the line connectivity.
  !************************************************************************

  !**************************************************************************
  ! Influence of the fixed boundaries on the segments displacement
  ! Absdep is reduced in such way that segments are stop in front of boundaries
  !**************************************************************************

  ABSDEPrail             = ABSDEP
  DEP_REAC(1:5,1:ABSdep) = IZERO
  domchange              = IZERO

  if ((GB/=IZERO)) then
      !*ATTENTION AUX MODULOS,ils ne sont retablis qu'en fin de module, il faut les remettre
      !*pour traiter les barrieres
      Oitest(:)  = modulo(Oi(:),modur(:))
      Eitest(:)  = modulo(Ei(:),modur(:))
      OiPtest(:) = Oitest(1:3)+nint((ABSdep*VDOi(1:3)),DPI)
      EiPtest(:) = Eitest(1:3)+nint((ABSdep*VDEi(1:3)),DPI)
      domchange=IZERO
      DomOold = seg(i)%dom
      DomEold = seg(i1e)%dom

      if (NbCvxdom > IUN) then

        ListCvxdomsegI(1:NbCvxdom) = ListCvxdom(1:NbCvxdom)
        ListCvxdomsegI(1) = DomOold

        ListCvxdomsegI(DomOold) = 1
        if (DomOold /= DomEold) then
          do ll = 1,Nbcvxdom
            if(ListCvxdomsegI(ll)==DomEold) exit
          enddo
          ve = ListCvxdomsegI(2)
          ListCvxdomsegI(2) = DomEold
          ListCvxdomsegI(ll) = ve
          if (kkdebug) write(379,*)  ListCVXdomsegi(:),domEold
        endif
      else
        ListCvxdomsegI(1) = IUN
      endif
      ! On test la position par rapport aux barrieres
      select case(GB)
      case (IUN)
         call domain_barriers(i,i1o,i1e,OiP,EiP,bloquer,NdomO,NdomE,lonad,lonvoad,lonvead,domchange,surfR,SFP)
         if (.not. bloquer) DEP_REAC(1:5,absdep) = (/NdomO,NdomE,surfR,SFP,domchange/)
      case (IDEUX)
          call barriere_spherique(OiPtest,EiPtest,bloquer)
      end select

      !* bloquer : true : the segment can not perfomr all predicted displacement
      if ((bloquer .or. domchange>IZERO)) then

        AbsDepMax=AbsDep

        !* voir commentaire 10 lignes plus  bas
        if ( ABSdep /= IUN .and. ABSDEP /= IZERO ) then

          domchange=IZERO

          ! On recherche la premiere position qui depasse les barrieres
          if (kkdebug) write(379,*) "INITIAL DISPLACEMENT MUST be REDUCED"

          do AbsDeptest=modep,AbsDepMax,modep
          if (kkdebug) write(379,*) "TESTING MAXIMUN POSSIBLE DISPLACEMENT, ACTual disp=",ABsdeptest
              ! test sur les barrieres
              OiPtest(:)=Oitest(1:3)+nint((ABSdeptest*VDOi(1:3)),DPI)
              EiPtest(:)=Eitest(1:3)+nint((ABSdeptest*VDEi(1:3)),DPI)
              ! On test la position par rapport aux barrieres
              select case(GB)
                case (IUN)
                    lonadtest   = LONGSEGi+nint(fac2dep(Li,LiOBR,LiEBR)*ABSdeptest*SIGdep,DPI)
                    lonvoadtest = SEG(i1o)%NORME+nint(fac1dep(Li,LiOBR)*ABSdeptest*SIGdep,DPI)
                    lonveadtest = SEG(i1e)%NORME-nint(fac1dep(Li,LiEBR)*ABSdeptest*SIGdep,DPI)
                    call domain_barriers(i,i1o,i1e,OiPtest,EiPtest,bloquer,NdomO,NdomE,&
                                        lonadtest,lonvoadtest,lonveadtest,domchange,surfR,SFP)

                    if (.not. bloquer) then
                       DEP_REAC(1:5,absdeptest) = (/NdomO,NdomE,surfR,SFP,domchange/)
                    endif
                    if (domchange > IZERO .and. .not. bloquer) then
                       AbsDep=AbsDeptest
                       exit
                    endif
                case (IDEUX)
                    call barriere_spherique(OiPtest,EiPtest,bloquer)
              end select


              if (bloquer) then
                  if (kkdebug) write(379,*) "Maximum displacement reached : ", (AbsDeptest-IUN/modep)*modep
                  ! On est alle une case trop loin
                  ! On garde la derniere valeur de AbsDep admissible
                  !if (seg(i)%surface>ITROIS) then
                  AbsDep=AbsDeptest-modep
                  absdep=(absdep/modep)*modep ! Make sure that displacement is compatible with sub-network

                  exit
              endif
          enddo
          !* Pour IABSdep=1 on saute la boucle et on met directement a zero (gain de temps)

        else
          if (domchange < IUN) then
            ABSdep = IZERO
            seg(I)%resdep =  zero
            IDEP(I)= sigdep*absdep
          endif
        endif
        !Attention Absdep a ete recalcule
        absdep=(absdep/modep)*modep
        OiP(1:3) = Oi(1:3)+nint((ABSdep*VDOi(1:3)),DPI)
        EiP(1:3) = Ei(1:3)+nint((ABSdep*VDEi(1:3)),DPI)
        lonad   = LONGSEGi+nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep,DPI)
        lonvoad = SEG(i1o)%NORME+nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep,DPI)
        lonvead = SEG(i1e)%NORME-nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep,DPI)
        if (domchange < IUN) then
          if (surfR == IZERO .and. bloquer) seg(i)%bloquer = .true.
        else
          if (seg(I)%surface > IZERO)  seg(i)%zerotl = .True.
        endif
      else
         DEP_REAC(1:5,absdep) = (/NdomO,NdomE,surfR,SFP,domchange/)
      endif

  endif

  if (absdep <= IZERO) then
      seg(i)%diseg = (Icinq * longsegI > Lomai)
      if (kkdebug) then
        write (379,*) " -----------------------------------------------------------------"
        write (379,fmt='( "Exiting SURI loop for segment =",I7," because of absdep = 0 after taking GB into account ")') i
        write (379,*) " -----------------------------------------------------------------"
      endif
      if (seg(i)%surface == 1) seg(i1o)%veclin=conec(seg(i)%veclin,1)
      if (seg(i)%surface == 2) seg(i1e)%veclin=conec(seg(i)%veclin,2)
      cycle surI
  elseif (out(i)) then
      if (kkdebug) then
        write (379,*) " -----------------------------------------------------------------"
        write (379,fmt='( "Exiting SURI loop for segment =",I7," because of out(i) = True after taking GB into account ")') i
        write (379,*) " -----------------------------------------------------------------"
      endif
      cycle surI
  endif

  !*************************************************
  ! End of the boundaries correction on displacement
  !*************************************************

  !If the segment i is going to cross a free surface will be eliminated
  !we don't need to look for obstacles
  if (seg(i)%surface < ITROIS) then

    !**************************************************************************
    ! Influence of the particules on the segment displacement
    ! Absdep is reduced in such way that segments are stop in front of boundaries
    !**************************************************************************
    if (particules .and. kk > relax_int) then

      ! We look for particles in box BI
      if (Npartester(BI) > izero ) then

        if (kkdebug)  then
            write (379,*) " "
            write (379,*) "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
            write (379,*) "&&&&&&&&&&&&&& Starting detection of intersecting particules  &&&&&&&&&&&&&&&&&&&&&&&&"
            endif
            pinned = .false.
            absdep_par = absdep
            normvdi = normdep(Li)
            normsegI = longsegI * normlin(Li)  ! the length of segment I in a
            ! we select the normal to the slip plane such that (OiEi, OiOip, nvec) is right-hand coordinate system
            if(iprodsca (nvec,prodivec(Ei-Oi,Oip-Oi)) < zero) nvec(:) = -IUN * nvec(:)
            ! VDO: vector of displacement of neighbor in O: the outside vector normal to OiOip
            VDO = prodivec(Nvec,Oip-Oi)
            ! VDE: vector of displacement of neighbor in E: the outside vector normal to EiEip
            VDE = prodivec(Eip-Ei,Nvec)
            normnvec  = dsqrt(real(Nvec(1)*Nvec(1) + Nvec(2)*Nvec(2) + Nvec(3)*Nvec(3),DP))
            normVDO  = dsqrt(real(VDO(1)*Vdo(1) + Vdo(2)*Vdo(2) + Vdo(3)*Vdo(3),DP))
            normVDE  = dsqrt(real(VdE(1)*VdE(1) + VdE(2)*VdE(2) + VdE(3)*VdE(3),DP))

            detect_par : do itemp = 1, Npartester(BI)

                ktemp = Listepartester(BI,itemp)
                ! The particle center
                Oj(:) = par(ktemp)%C(:)
                Oj(:) = modulo(Oj(:) + Tr(:) , modur(:)) - Tr(:)

                ! The particle radius
                LONGSEGj = par(ktemp)%R
                if(LONGSEGj <= izero) stop "taille particule null"
                ! calculation of the radius of the particle in the segment glide plane
                ! we first calculate the normal distance between the particle center and the seg plan
                OiOj(:) = Oj(:)-Oi(:)
                OiOjprojN = real(dot_product(OiOj,Nvec),DP)/normNvec

                if(kkdebug) write(379,*) "par:",ktemp,"H/B:",OiOjprojN,"size",LONGSEGj
                ! since there is no easy procedure to confirm that the particle toughing the trapezium,
                ! we proceed by ignoring all the other particles following deferrent test

                !(1) The segment glide plane does not intersect the particle
                if (abs(OiOjprojN) >= LONGSEGj ) cycle DETECT_par
                ! The particle radius of intersection in the seg glide plan
                radius = DSQRT(LONGSEGj**2 - OiOjprojN**2)
                ! And all the calculations are now in the segment glide plane
                ! computing of the center of the intersection of the particule with the slip plane
                Oj(:) = Oj(:) - nint(Nvec(:)* OiOjprojN/normNvec )
                ! We consider the trapezium of motion formed by Oi,Ei,Eip,OiP
                ! (2) CYCLE if the circle is to the right of the trapezium
                EiOj(:) = Oj(:)-Ei(:)
                temp = prodisca(EiOj,VDE)/normvdE
                if(kkdebug) write (379,*) " Distance from I1E ",temp
                if(temp >= radius)  cycle DETECT_par
                OiOj(:) = Oj(:)-Oi(:)
                ! (3) CYCLE if the circle is to the left of the trapezium
                temp = prodisca(OiOj,VDO)/normvdO ! the minus sign because VDI pointing forwards
                if(kkdebug) write (379,*) " Distance from I1O ",temp
                if(temp >= radius)  cycle DETECT_par
                ! (4) CYCLE if the circle is BEHINDE the trapezium, ie behinde OiEi
                distIJ = prodisca(OiOj,VDI)/normvdi
                if(kkdebug) write (379,*) " Distance from I in the VDI direction ",distIJ
                ! the minus sign because VDI pointing forwards
                if(distIJ <= -radius)  cycle DETECT_par
                ! (3) CYCLE if the circle is in front of the trapezium, ie.in front of OipEip
                if(distIJ >= radius + absdep*normvdi)  cycle DETECT_par
                if(kkdebug) write (379,*) " particle ",ktemp, " touchsvery close to the trapezium, distIJ = ", distIJ
                ! Now, we know that the particle is very close, but not sure it cuts the trapezium
                ! because of the angle effects
                ! We start by testing if the segment is already touching the particle on its
                ! ends or between ends
                Oit = Oi ; Eit = Ei

                do I1 = 0 , absdep
                    if(kkdebug) write (379,*) " test intersection after absdep = ",I1
                    Ci(:) = OiOj(:) - int(distIJ/normvdi) * VDI(:) ! Ci is the projection of OiOj on OiEi
                    if(kkdebug) write (379,*) " OiOj(:)",OiOj(:)
                    if(kkdebug) write (379,*) " proj ", int(distIJ/normvdi) * VDI(:)

                    temp = real(dot_product(Ci,Eit-Oit),DP)/normsegI
                    if(kkdebug) then
                        write (379,*) "distance from Oi = ",norivect(OiOj), " from Ei = ", norivect(EiOj)
                        write (379,*) " projection on OiEi = ", temp
                    endif
                    if (norivect(OiOj) < radius .or. norivect(EiOj) < radius .or. &
                        (temp > zero .and. temp < normsegI)) then
                        if(I1 == izero) then
                            ! the segment touchs the particles before moving
                            pinned = .true.
                            selected_par = ktemp
                            if(kkdebug) write (379,*) " particle ",ktemp, " touchs segI before motion"
                            exit detect_par
                        else
                            ! SegI touchs the parI after displacment of I1
                            if(I1 < absdep_par) then
                              selected_par = ktemp
                              absdep_par = I1
                              if(kkdebug) write (379,*) " particle ",ktemp, " touchs segI after absdep = ",I1
                              cycle DETECT_par
                            endif
                        endif
                    else
                        ! if after displacment of I1 no intersection, we move SegI by 1 step
                        Oit(:) = Oit(:) + int(VDOi(1:3)) ! the position Oi as if SegI moved by 1 step
                        Eit(:) = Eit(:) + int(VDEi(1:3))! the position Ei as if SegI moved by 1 step
                        OiOj(:) = Oj(:) - Oit(:)
                        EiOj(:) = Oj(:) - Eit(:)
                        distIJ = real(prodisca(OiOj,VDI),DP)/normvdi
                        normsegI = norivect(Eit-Oit)
                    endif
                enddo
                ! we compute Oi and Ei because their values changed in detection of particles
            enddo detect_par

            if (kkdebug)  then
                write (379,*) "&&&&&&&&&&&&   END of detection of intersecting particules  &&&&&&&&&&&&&&&&&&&&&&&&&&"
                write (379,*) "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
                write (379,*) " "
            endif

            if (.not. pinned .and. absdep_par == absdep ) then
                if(kkdebug) write (379,*) " No particles are obstacles for the motion"
            elseif(pinned) then
                if (par(selected_par)%tau < zero) then
                   ! this is the convention of Orowan particles. No penetration
                   Jump = .false.
                else
                   taueff = tautot(I)*xmu
                   jump = JUMP_Cr(taueff,normsegi,selected_par,radius)
                endif

                if(jump) then
                    if (kkdebug) then
                      write (379,*) "KK",kk, "I =",I,"absdep = ", absdep, "absdep_par =",absdep_par
                      write (379,*) " SegI is activated and allowed to cross the particle ",selected_par
                      endif
                    else
                      if(normsegi > deux * radius) then
                          !write (379,*)"Y" ,I, longsegI*normlin(Li),radius,kk
                          seg(i)%diseg = .true.
                          !read(*,*)
                      endif
                      if (kkdebug) then
                          write (379,*) " segI intersect particle but NOT activated "
                          write (379,*) " seg(I)%diseg is set true "

                          write (379,*) " -----------------------------------------------------------------"
                          write (379,*) " sortie SURI pour absdep = 0 apres PARTICULE: contact non-active  "
                          write (379,*) " -----------------------------------------------------------------"
                      endif
                      cycle surI
                    endif
                elseif(absdep_par < absdep) then
                    if(absdep_par == izero) stop "absdep_par = zero and not pinned ? kk, I = "
                    absdep = absdep_par
                    if (kkdebug) then
                       write (379,*) " reduction of absedep after detection of particles      -------"
                       write (379,*) " The new absdep is =  ",absdep
                       write (379,*) " -----------------------------------------------------------------"
                    endif
                    ! if absdep has been changed, one must update Oip EiP, and lengths after disp
                    absdep=(absdep/modep)*modep
                    OiP(1:3) = Oi(1:3)+nint((ABSdep*VDOi(1:3)),DPI)
                    EiP(1:3) = Ei(1:3)+nint((ABSdep*VDEi(1:3)),DPI)
                    lonad   = LONGSEGi+nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep,DPI)
                    lonvoad = SEG(i1o)%NORME+nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep,DPI)
                    lonvead = SEG(i1e)%NORME-nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep,DPI)
                endif
            endif
        endif
        !*************************************************************
        !End of the section accounting for the influence of particules
        !*************************************************************

        !**************************************************************************
        ! Beginning of the forest obstacle procedure
        !**************************************************************************
        if (kkdebug) then
          write(379,*) '>>>>>>>>>>>>>>>>>>>>>Debut DETECT<<<<<<<<<<<<<<<<<<<<<<'
          write(379, '("i=",I6," Iboite =",I5," n seg a tester =",I6)')  i,iboite(i),NsegTester
        endif

        !----------------------------------------------------------------------------------------
        ! If we are in the relaxation steps or if a segment is touching a surface,
        ! the forest obstacle detection is not applied.
        !----------------------------------------------------------------------------------------
        if (kk > relax_int .and. NsegTester > izero .and. seg(i)%surface == izero) then
        !if (NsegTester > izero .and. seg(i)%surface == izero) then
        !***Avoid formation of junction during initial steps. Used for P. Geantil calculations.
        !if (kk > relax_int .and. NsegTester > izero .and. seg(i)%surface == izero &
        !  .and. kk > 3000) then
        !***

        !*******************************************************************
        ! Beginning of the forest obstacle detection loop.
        ! At the end the list of obstacles and their specificity are defined
        !*******************************************************************

        !initialization
        distiobs = absdep

        if (kkdebug) write(379,*) "==============nb segtester ", Nsegtester, "=== Nb groseg", Nb_groseg

        ! recherche du VD du rail en o oriente vers l'aire balayee de i
        normsqrVDi = VDi(1)*VDi(1) + VDi(2)*VDi(2) + VDi(3)*VDi(3)  !carre de la norme de VDi

        !*********************   DETECT_LOOP  ***************************
        ! This subroutine is a loop on all possible obstacles segments  *
        ! the most important obstacles are listed in arrays OBSTAB      *
        ! In parralel calculation is is favorable two repeat the call   *
        ! to reduce the number of data exchanged in the MPI Comm        *
        ! With the first pass we only buildup a reduce listeSegTester.  *
        !****************************************************************
#ifdef PA
    if (TwoPass .and. ((TAILLE_PAIRIMPAIR * 2) < NsegTester) .and. (TAILLE_PAIRIMPAIR /= 1)) then

      if (NsegTester /= 0)    then

        ! The cases where all segments must be shifted with a modulo
        if (allshift .or. NBOITES < IDEUX .or. key_box_2border(BI)) then

          do ktemp = 1, nsegtester
            j = listeSegTester(ktemp)

            ax = seg(j)%o(1) + Trx
            ay = seg(j)%o(2) + Try
            az = seg(j)%o(3) + Trz
            seg(j)%o(1) = (ax - int(floor(real(ax)*invmodurx)*modurx,DPI)) - Trx
            seg(j)%o(2) = (ay - int(floor(real(ay)*invmodury)*modury,DPI)) - Try
            seg(j)%o(3) = (az - int(floor(real(az)*invmodurz)*modurz,DPI)) - Trz
          enddo

        ! The particular case of long segments with the Greengard algo
        elseif(.not. key_box_2border(BI) .and. NBOITES /= iun) then

          do ktemp = (nsegtester-Nb_Groseg), nsegtester
            j = listeSegTester(ktemp)

            if (longseg(j)) then
              ax = seg(j)%o(1) + Trx
              ay = seg(j)%o(2) + Try
              az = seg(j)%o(3) + Trz
              seg(j)%o(1) = (ax - int(floor(real(ax)*invmodurx)*modurx,DPI)) - Trx
              seg(j)%o(2) = (ay - int(floor(real(ay)*invmodury)*modury,DPI)) - Try
              seg(j)%o(3) = (az - int(floor(real(az)*invmodurz)*modurz,DPI)) - Trz
            endif
          enddo

        endif

        ! First pass
        FirstTime = .true.              ! Is it the first time the subroutine DETECT_LOOP is called
        call Detect_loop(NsegTester,nbobs,FirstTime)

        ! A second pass with a tiny ListSegTester array is prepared
        NsegTester = MiniNsegTester
        listeSegTester(1:NsegTester) = MiniListeSeg(1:NsegTester)
        distiobs = absdep
        FirstTime = .false.

        ! Second pass
        if (NsegTester /= 0)    call Detect_loop(NsegTester,nbobs,FirstTime)

      endif

    else
#endif
      ! Only one call to DETECT_LOOP is made
      FirstTime = .true.      ! Is it the first time the subroutine DETECT_LOOP is called
      if (NsegTester /= 0)  then

        ! The cases where all segments must be shifted with a modulo
        if (allshift .or. NBOITES < IDEUX .or. key_box_2border(BI)) then

          do ktemp = 1, nsegtester
            j = listeSegTester(ktemp)

            ax = seg(j)%o(1) + Trx
            ay = seg(j)%o(2) + Try
            az = seg(j)%o(3) + Trz
            seg(j)%o(1) = (ax - int(floor(real(ax)*invmodurx)*modurx,DPI)) - Trx
            seg(j)%o(2) = (ay - int(floor(real(ay)*invmodury)*modury,DPI)) - Try
            seg(j)%o(3) = (az - int(floor(real(az)*invmodurz)*modurz,DPI)) - Trz
          enddo

        ! The cases of long segments with the Greengard algo
        elseif(.not. key_box_2border(BI) .and. NBOITES /= iun) then

          do ktemp = (nsegtester-Nb_Groseg), nsegtester
            j = listeSegTester(ktemp)

            if (longseg(j)) then
              ax = seg(j)%o(1) + Trx
              ay = seg(j)%o(2) + Try
              az = seg(j)%o(3) + Trz
              seg(j)%o(1) = (ax - int(floor(real(ax)*invmodurx)*modurx,DPI)) - Trx
              seg(j)%o(2) = (ay - int(floor(real(ay)*invmodury)*modury,DPI)) - Try
              seg(j)%o(3) = (az - int(floor(real(az)*invmodurz)*modurz,DPI)) - Trz
            endif
          enddo

        endif

        call Detect_loop(NsegTester,nbobs,FirstTime)

      endif
#ifdef PA
    endif
#endif

      ! If the amplitude of the segment displacement is less than one elementary
      ! displacement on the simulation network, We discretize the segment in order
      ! bypass the close obstacle
      if (Distiobs < modep) then
        seg(i)%diseg = (longsegI > Limite_discreti)
        ! remettre absdep a l'echelle apres reduction eventuel de deplacement
        ! a acuse des obstacles

        if (kkdebug) then
            write(379,*) " -----------------------------------------------------------------"
            write(379,*) " sortie du SURI pour absdep = 0 apres detection des obstacles"
            write(379,*) " ----------------------------------------------------------------"
        endif
        cycle surI
      endif

      !***************************************
      ! Identification of the type of reaction
      !***************************************
      do KTEMP = 1, nbobs

        if (OBSTAB(KTEMP)%react/=iquatre)then    !pas de pb de ss reseau

          j  = OBSTAB(KTEMP)%numseg
          Lj = seg(j)%veclin

          if (sysobst(Li,Lj) == iun) then

            ! i and j belong to the same or collinear slip system
            VLj(:) = bveclin(:,Lj)

            if ((VLj(1) + VLi(1) == izero) .and.    &
                (VLj(2) + VLi(2) == izero) .and.    &
                (VLj(3) + VLi(3) == izero)) then

              ! annihilation is possible if ptobs is (-1,-1,-1) with covering
              if (obstab(ktemp)%ptobs(1) /= -iun  .or.   &
                  obstab(ktemp)%ptobs(2) /= -iun  .or.   &
                  obstab(ktemp)%ptobs(3) /= -iun) then
                if (OBSTAB(ktemp)%react/=icinq) OBSTAB(ktemp)%react = itrois ! without covering, the reaction is treated as a junction
              else
                OBSTAB(ktemp)%react =  ideux !cas annihilation (j m SG ou vis SDevie)
              endif

            elseif (OBSTAB(ktemp)%react/=icinq) then
              OBSTAB(ktemp)%react = iun !i et j m SG mais pas annihilation
            endif

          ! i and j belong to different slip
          elseif (OBSTAB(ktemp)%react /= icinq .and. sysobst(Li,Lj) == isix)then    !sysobst = 6 pour systemes /= non coplanaires

              OBSTAB(ktemp)%react = itrois ! jonction possible

          endif

          if (SEG(j)%gd > izero .and. OBSTAB(ktemp)%react /= ideux) OBSTAB(ktemp)%react = isix

        endif

      enddo
      !***********************************
      ! End of the reaction identification
      !***********************************
      ! the out put of the last loop is to define OBSTAB(jselec)%react:
      ! 1 : same system without annihilation
      ! 2 : annihilation
      ! 3 : jonction
      ! 4 : probelem of sub lattice
      ! 5 : obsta is already junction
      ! 6 : obsta is already GD

      !*************************************
      ! End of the forest obstacle detection
      !*************************************

      if (kkdebug) write(379,*) '>>>>>>>>>>>>>>>>>>>>> entre DETECT et SELEC<<<<<<<<<<<<<<<<<<<<<<'

      if (nbobs > 0) then
        if (kkdebug) then
            write(379,*) 'Nbobstot:', Nbobs
            do Ktemp = 1, Nbobs
                write(379,*) ' ----------------------------  seg:', obstab(ktemp)%numseg
                write(379,*) ' ptobs:', obstab(ktemp)%ptobs(1:3)
                write(379,*) ' dist:', distiobs, ' plani:', obstab(ktemp)%plani, ' react:', obstab(ktemp)%react
                call conf(obstab(ktemp)%numseg)
            enddo
            write(379,*) "Fin des caracteristiques d'obstacless >>>>>>>>>>>>>>>>>>"
        endif
      else
        if (kkdebug) write(379,*) " il n'y pas d'obstacles detecte -----------------------------"
      endif

      !***********************************************************
      ! Selection of the obstacle segment to consider for reaction
      !***********************************************************
      ! the aim is to select the obstacle with which i will react. there are two cases :
      ! - there is only one obstacle => easy
      ! - there are some obstacles => we have to detect if there is a probleme of three body or if they are just
      !   neighboors

      ! variable initialization
      ! ------------------------------------------------------------------
      ! distiobs : distance from the closest obstacle
      ! obsta : distance from the selected obstacle when the reaction is possible
      absdepobsta = int(distiobs,DPI)
      ! ------------------------------------------------------------------
      codeobst = izero
      ptobst(1:3)= (/0,0,0/)      !* position de l'obstacle

      !* si il y a un obstacle le plus proche, plusobs = 0, si il y a plus d'obstacles aussi proches
      ! plusobst doit etre > 0
      if (nbobs /= izero) then
        Lvoitmp = izero ! line vector of obstacle jtmp
        jselec = -IUN
        !==================================================================================================
        ! the next IF conditions are used to detect three-body problems
        !==================================================================================================
        ! if there are more than 3obstacles at the same distance, it's too difficult to study theses cases
        ! even if they are neighboors
        ! In the same time we select the obstacle in the three following case :
        ! nbobs = 1  ; 2 obstacles forming the smae junction with possible annihilation ; nbobs = 3
        select case (Nbobs)
        case (1)
            ! if there is only one obstacle, it is the selected one
            jselec = iun
            pb3b = .false.

        ! case of annihilation with one of the two segments participating to a junction
        case (2)

            if (seg(obstab(iun)%numseg)%ijonc == obstab(ideux)%numseg) then
              ! this specific case is considered as 3 body problem but it may be easily treated
              ! when an annihilation is predicted with one of the segmemnt
              pb3b = .true.
            if (obstab(1)%react /= IDEUX .and. obstab(2)%react /= IDEUX)then
                Jselec = izero
              else
                if (obstab(1)%react==2) jselec = IUN
                if (obstab(2)%react==2) jselec = IDEUX
              endif

            ! the obstacles are not a junction binome
            else

              itemp = obstab(IUN)%numseg

              oo = seg(itemp)%vnno
              nvoi = iun
              do while (seg(oo)%norme == zero)
                oo = seg(oo)%vnno
                nvoi = nvoi + IUN
                if (nvoi > 19) then
                  print *, " Erreur 124 : contact kk,I =",KK,I
                  call seginfo(i,"Pb Erreur 124 i")
                  call seginfo(itemp,"Pb Erreur 124 i")
                  stop
                endif
              enddo

              oe = seg(itemp)%vnne
              nvoi = iun
              do while (seg(oe)%norme == zero)
                oe = seg(oe)%vnne
                nvoi = nvoi + IUN
                if (nvoi > 19) then
                  print *, " Erreur 125 : contact kk,I =",KK,I
                  call seginfo(i,"Pb Erreur 125 i")
                  call seginfo(itemp,"Pb Erreur 125 i")
                  stop
                endif
              enddo

              ! the two obstacles are not neighboors = > 3 body problem
              if (oo /= obstab(ideux)%numseg .AND. oe /= obstab(ideux)%numseg) THEN
                  pb3b = .true.
                  JSELEC = IZERO
                  ! the two obstacles are neighboors
              else
                  pb3b = .false.
                  ! if an obstacle is GD, only annihilation can be made
                  if ((SEG(obstab(1)%numseg)%gd + SEG(obstab(2)%numseg)%gd) > izero .and. &
                      (obstab(1)%react /= 2 .and. obstab(2)%react /= 2)) then
                    if (kkdebug) write(379,*) "pas d'obstacle selectione car un segment GD"
                    jselec = izero
                  else
                    ! coplanar case ==========================================================
                    if ( OBSTAB(1)%plani .and. OBSTAB(2)%plani ) then
                        if (kkdebug) write(379,*) 'SELEC: obstacle colineaire'
                        ! the 2 obstacles belong to the plane of i, this case cooresponds to the interaction
                        ! with the same system. the possible interaction is annihilation.
                        ! If this is not possible sgmt i stays in front of the obstacle
                        jselec = izero
                        if ( obstab(1)%react == 2 ) jselec = IUN
                        if ( obstab(2)%react == 2 ) jselec = IDEUX

                        ! neither coplanar with I
                        !===================================================================================
                    elseif (.not. OBSTAB(1)%plani .and. .not. OBSTAB(2)%plani) then
                        ! no segment is in the plane of i, it corresponds to junction even if the system is collinear
                        ! because no covering between segments
                        if (kkdebug) write(379,*) 'SELEC: obstacle non colineaire'
                        ! we should check if ptobst for the two segments is the same
                        if(obstab(IUN)%ptobs(1) /= obstab(Ideux)%ptobs(1) .or. &
                            obstab(IUN)%ptobs(2) /= obstab(Ideux)%ptobs(2) .or. &
                            obstab(IUN)%ptobs(3) /= obstab(Ideux)%ptobs(3)) then
                          jselec = izero
                          ! the same point for obstacle
                        else
                          ! it is not necessary to search a kneecap if there is a sub lattice problem
                          if (obstab(iun)%react  /= itrois .and. obstab(ideux)%react /= itrois) then
                              jselec = izero
                          else
                              temp1 = SEG(obstab(1)%numseg)%veclin
                              temp2 = SEG(obstab(2)%numseg)%veclin
                              ! to be able to perfrom the interaction the two obstacles must belong
                              ! to the same system
                              if (etatsys(temp1,temp2) /= izero) THEN
                                JSELEC = IZERO
                              ELSE
                                ! IF  oO /= obstab(ideux)%numseg, THEN  obstab(ideux)%numse is at the end of IUN
                                ! temp2 is the neighbor nn at the end of temp1 : temp1 -> temp2
                                if (oo /= obstab(ideux)%numseg) THEN
                                    temp1 =  obstab(iun)%numseg
                                    temp2 = obstab(ideux)%numseg
                                    Ltemp = .true.
                                else
                                    temp2 =  obstab(iun)%numseg
                                    temp1 =  obstab(ideux)%numseg
                                    Ltemp = .false.
                                endif

                                voitmp = temp1 ! neighboor of obstacle 1
                                nvoi = iun
                                ! new procedure: it is possible to find a kneecap between the two obstacles having alreading
                                ! the junction character
                                do while (SEG(voitmp)%norme == zero .and. Lvoitmp /= Ljj) ! peut on avoir plus de 4
                                                                                          ! rotules sans changer de sens?
                                    ! search of the extremity which belong to the plane of i
                                    voitmp = SEG(voitmp)%voiso ! neighboor of obstacle 1
                                    Lvoitmp = SEG(voitmp)%veclin
                                    Ljj = axejoncj(Li,Lvoitmp)  ! corresponds to the junction line between
                                                                ! the systems of Li and Ljselec
                                    Lji = axejonci(Li,Lvoitmp)
                                    nvoi = nvoi + iun
                                    if (nvoi > 19) exit
                                enddo
                                ! if nvoi > 10 then no kneecap was found convenient for the junction
                                if (nvoi > 19 .or. SEG(voitmp)%norme /= izero)  then
                                    if (kkdebug) write(379,*) 'pas de rotule bien orientee'
                                    Jselec = iUN
                                    if(seg(temp1)%veclin /=  seg(temp2)%veclin) then
                                      ! we have to choose between temp1 and temp2
                                      ! the creterion is the solution minimising the total number of kneecaps
                                      ! we start by testing temp1
                                      jtemp = seg(temp1)%veclin
                                      Ktemp = seg(temp2)%veclin
                                      Ljj = axejoncj(Li,jtemp)  ! corresponds to the junction line between
                                                                ! the systems of Li and Ljselec
                                      Lji = axejonci(Li,jtemp)

                                      nvoi = nbrot(Li,Lji,1) + nbrot(Lji,Li,1) + nbrot(jtemp,Ljj,1) + nbrot(Ljj,ktemp,1)
                                      ! calculation of the total number of kneecaps to introduce
                                      ! when temp2 is selected for the junction
                                      Ljj = axejoncj(Li,ktemp)  ! corresponds to the junction line between
                                                                ! the systems of Li and Ljselec
                                      Lji = axejonci(Li,ktemp)
                                      ! if the total num of kneecaps temp2 is lesse, we modify the default value fo jselec
                                      if((nbrot(Li,Lji,1)+nbrot(Lji,Li,1)+nbrot(jtemp,Ljj,1)+nbrot(Ljj,ktemp,1)) < nvoi) then
                                          if (Ltemp) Jselec = ideux
                                          JoncOj = .true.
                                      else
                                          if (.not. Ltemp) Jselec = ideux
                                          JoncEj = .true.
                                      endif
                                    else
                                      ! obstacles have the same veclin then no preference.
                                      ! the default: one choose the obstacle that ptobs is its origine
                                      if (Ltemp) Jselec = ideux
                                      Lvoitmp = SEG(jselec)%veclin
                                      nvoi = SEG(jselec)%veclin
                                      Ljj = axejoncj(Li,Lvoitmp)  ! corresponds to the junction line between
                                                                  ! the systems of Li and Ljselec
                                      Lji = axejonci(Li,Lvoitmp)
                                      JoncOj = .true.
                                    endif
                                else
                                    if (kkdebug) write(379,*) 'rotule qui va bien :', voitmp
                                    jselec = itrois
                                    obstab(jselec)%numseg = voitmp
                                    obstab(jselec)%ptobs(:) = obstab(IUN)%ptobs(:)
                                    obstab(jselec)%react = itrois
                                    obstab(jselec)%plani = .true.
                                    kneecap = .true.
                                endif
                              endif
                          endif
                        endif
                        ! systems:  annhilitation with screw segment of a collinear system or junction
                        ! => search of one obstacle well oriented for the reaction : one belonging to the plan of i

                        ! case ofone obstacle which belongs to the plane of i, the other is not .
                        ! we select the coplanar
                    else
                        Jselec = IUN
                        if (OBSTAB(IDEUX)%plani) Jselec = IDEUX
                        if (kkdebug) write(379,*) 'Contact: SELEC: 1 obstacle copla et 1 obstacle non copla'
                    endif
                  endif
              endif
            endif
        case (3)
            ! for 2 or 3 obstacles detected, we have to check if they are connected or if they are 2 different obstacles
            if ((SEG(obstab(1)%numseg)%gd + SEG(obstab(2)%numseg)%gd + SEG(obstab(3)%numseg)%gd) > izero .and. &
                (obstab(1)%react /= 2 .and. obstab(2)%react /= 2 .and. obstab(3)%react /= 2)) then
              jselec = izero
              if (kkdebug) write(379,*) "pas d'obstacle selectione car un segment GD"
            Else
              jselec = izero
              do jtemp = IUN, nbobs
                  itemp = obstab(jtemp)%numseg
                  oo = seg(itemp)%vnno
                  nvoi = iun
                  do while (seg(oo)%norme == zero)
                      oo = seg(oo)%vnno
                      nvoi = nvoi + IUN
                      if (nvoi > 19) then
                        print *, " Erreur 224 : contact kk,I =",KK,I
                        call seginfo(i,"Pb Erreur 224 i")
                        call seginfo(itemp,"Pb Erreur 224 i")
                        stop
                      endif
                  enddo
                  oe = seg(itemp)%vnne
                  nvoi = iun
                  do while (seg(oe)%norme == zero)
                      oe = seg(oe)%vnne
                      nvoi = nvoi + IUN
                      if (nvoi > 19) then
                        print *, " Erreur 225 : contact kk,I =",KK,I
                        call seginfo(i,"Pb Erreur 225 i")
                        call seginfo(itemp,"Pb Erreur 225 i")
                        stop
                      endif
                  enddo
                  ! oo and oe are neighbors with non zero negth of the segment itemp
                  if(jtemp == IUN) then
                      io = ideux
                      ie = itrois
                  elseif(jtemp == Ideux) then
                      io = iun
                      ie = itrois
                  else
                      io = iun
                      ie = ideux
                  endif
                  ! itemp: number of segment corresponding to the JTEMP saved obstacle
                  ! the following condition is fullfilled when the tow other obstacles are the neighbors of Itemp
                  if ((oo == obstab(io)%numseg .or. oo == obstab(ie)%numseg) .and. &
                        (oe == obstab(io)%numseg .or. oe == obstab(ie)%numseg)) then

                      ! the following condition is true when obstacle num Jtemp is parralell to I
                      ! but not the other two obstacles : ONLY THIS CASE CAN BE TREATED
                      ! all other cases considered as 3 body problem
                      if ((obstab(jtemp)%ptobs(1) == -IUN)  .and. &
                          (obstab(io)%ptobs(1) /= -IUN)  .and. &
                          (obstab(ie)%ptobs(1) /= -IUN)) then
                        jselec = jtemp
                        pb3b = .false.
                        exit
                      endif
                  endif
              enddo
            endif
        case(4:100)
            jselec = izero
            pb3b = .true.
        end select
        !==================================================================================================
        ! obstacle selection between real neighbors (or binome of junction) when no 3-body problem         !
        !==================================================================================================
        ! for nbobs = 1, 3 , the obstacle was selected above
        ! When Nbobs = 2, we do the following

        if (kkdebug) write(379,*) 'pb a 3 corps:', pb3b
        !  Here the obstacle is already choosen and known

        !==================================================================================================
        ! reaction determination                                                                          !
        !==================================================================================================
        ! this last part,defines : the two variables jnct and anseg,
        ! absdep taking into account obstacles : absdepobsta
        ! obsta : the definitive segment used after select
        ! codeobst : corresponds to the type of interaction

        ! the final obstacle is one of the detected segment
        !========================================================================================================
        OBSTA = izero                 !*** Pas d'obstacle
        if (jselec /= izero .and. jselec /= -iun)then

            select case (OBSTAB(jselec)%react)

            case(iun,iquatre,icinq,isix)
              ! - I: m SG; 4: sous reseau;  5: already junction
              ! - is yet a junction
              ! - make a subarray problem
              ! - if there is one segment CS and no annihilation
              ! there is no contact reaction, i stay in front of the obstacle
              !---------------------------------------------------------------------------
              if (LongSegI > limite_discreti) SEG(i)%diseg = .true.
              !---------------------------------------------------------------------------
              if (distiobs > modep) then
                  absdepobsta = int(distiobs,DPI) - IUN
              else
                  absdepobsta = izero
              endif
              codeobst = obstab(jselec)%react

            case (ideux)
              ! annihilation
              anseg = .true.
              obsta = obstab(jselec)%numseg
              codeobst = obstab(jselec)%react


            case (itrois)
              ! junction
              !  we have to check the franck's criteria to make reaction (include in axejoncj)
              ! if not we stay in front of the obstacle
              Ljselec = SEG(obstab(jselec)%numseg)%veclin
              Ljj = axejoncj(Li,Ljselec) ! corresponds to the junction line between the systems of Li and Ljselec
              Lji = axejonci(Li,Ljselec)
              if (kkdebug) write(379,*) 'Ljselec',Ljselec,Ljj
              ! jselec is in the intersection direction but is it in the good sens ?
              if (obstab(jselec)%plani) then
                  if (Ljselec/=Ljj) then
                    ! junction cannot be made, we stay in front of jselec
                    if (distiobs > modep) then
                        Absdepobsta = int(distiobs,DPI) - IUN
                    else
                        absdepobsta = izero
                    endif
                    if (kkdebug) write(379,*) 'SELEC : jonction impossible non energetiquement favorable'
                    !---------------------------------------------------------------------------
                    !    if (seg(i)%norme > LOMAi*half) SEG(i)%diseg = .true.
                    !---------------------------------------------------------------------------
                  else ! junction can be made
                    jnct = .true. ! key for junction
                    obsta = obstab(jselec)%numseg
                    ptobst = obstab(jselec)%ptobs
                    codeobst = obstab(jselec)%react
                  endif
              else
                  if (Ljj==0.or.Lji==0) then
                    ! junction cannot be made, we stay in front of jselec
                    if (distiobs > modep) then
                        Absdepobsta = int(distiobs,DPI) - IUN
                    else
                        absdepobsta = izero
                    endif
                    if (kkdebug) write(379,*) 'SELEC : jonction impossible non energetiquement favorable'
                    !---------------------------------------------------------------------------
                    ! if (seg(i)%norme > LOMAi*half) SEG(i)%diseg = .true.
                    !---------------------------------------------------------------------------
                  else ! junction can be made
                    jnct = .true.
                    obsta = obstab(jselec)%numseg
                    ptobst = obstab(jselec)%ptobs
                    codeobst = obstab(jselec)%react
                  endif
              endif

            end select


            ! there is no final obstacle or is not a detected one : a kneecap
            !========================================================================================================

        elseif (jselec == izero) then
            ! there is no reaction for the 3body problem and if there are two neighboors of the same
            ! glide system
            obsta = izero
            if (distiobs > modep) then
              absdepobsta = int(distiobs,DPI) - IUN
            else
              absdepobsta = izero
            endif

        elseif (jselec == -iun) then
            print *, '!!!PB : configuration d obstacles telle que l on ne peut pas selectionner!!!'
            print *, 'nbobs:', nbobs
            print *, 'kk:',kk,'i:',i,'nbobs:', nbobs
            stop
        endif

        if(obsta == izero) then
            !---------------------------------------------------------------------------
            if (LongSegI > limite_discreti) SEG(i)%diseg = .true.
            if (kkdebug) write(379,*) 'Forced discretisation: can not do interaction  %diseg ', SEG(i)%diseg
            !---------------------------------------------------------------------------
        endif

        if (kkdebug) then
            write(379,*) "apres DETECT ABSDEP = ",absdep," absdepobsta=", absdepobsta, " modep =",modep
            write(379,*) 'kk:',kk,'i:',i,'Nbobs:',Nbobs
            write(379,*) '- SELEC '
            write(379,*) 'jselec, kneecap',jselec,kneecap
            write(379,*) ' obsta:',obsta, ' jnct', jnct, ' anseg', anseg, ' pb3b', pb3b
            write(379,*) ' absdepobsta',absdepobsta,  'ptobst', ptobst
            write(379,*) 'apres selec l obstacle i:', i,' est j:',obsta
        endif

        !====================Attention : modification of ABSDEP   =========================================
        !==================================================================================================
        absdepobsta = (absdepobsta/modep) * modep            !nom a changer
        absdep = absdepobsta !nom a changer
        !==================================================================================================
        !==================================================================================================

        if (absdep <= izero) then
            seg(i)%diseg = (longsegI > limite_discreti)
            ! remettre absdep a l'echelle apres reduction eventuel de deplacement
            ! a acuse des obstacles
            if (kkdebug) then
              write(379,*) " -----------------------------------------------------------------"
              write(379,*) " sortie du SURI pour absdep = 0 apres SELECTION des obstacles"
              write(379,*) " ----------------------------------------------------------------"
            endif
            cycle surI
        endif

        !#######################################################################################################
        lonad   = LONGSEGi + nint(fac2dep(Li,LiOBR,LiEBR)*ABSdep*SIGdep)
        lonvoad = SEG(i1o)%NORME + nint(fac1dep(Li,LiOBR)*ABSdep*SIGdep)
        lonvead = SEG(i1e)%NORME - nint(fac1dep(Li,LiEBR)*ABSdep*SIGdep)

      endif     !  (IF nbobs /= zero)
      !******************************************************************
      ! End of selection of the obstacle segment to consider for reaction
      ! The main modification to segment connection are now defined.
      !******************************************************************

    ENDIF
  !****************************************************************
  ! End of the forest obstacle procedure     if(KK > relax_int ....
  !****************************************************************

  endif !!!end of if test to skip the obstacle detection

  if (GB/=IZERO .and. ABSDEP /=IZERO) then
    if (ABSDEPmax == izero .and. ABSDEP /= ABSDEPrail) then
       if (kkdebug) write(379,*) "calling barriere detection with ABSDEP = ",absdep

       OiP(1:3) = Oi(1:3)+nint((ABSdep*VDOi(1:3)),DPI)
       EiP(1:3) = Ei(1:3)+nint((ABSdep*VDEi(1:3)),DPI)

       call domain_barriers(i,i1o,i1e,OiP,EiP,bloquer,NdomO,NdomE,&
                            lonad,lonvoad,lonvead,domchange,surfR,SFP)
       seg(i)%dom = NdomO
       seg(i1e)%dom = NdomE
       seg(i)%surface= surfR
       seg(i)%VarfreePlan= SFP
       if (kkdebug) then
        write(379,*) "Domain assignement recalculated"
        write(379,*) "Dom i", NDomO
        write(379,*) "Dom i1e", NDomE
        write(379,*) "free surface reaction type", surfR
        write(379,*) "free surface", SFP
        write(379,*) "domchange", domchange
       endif
    else
      NdomO = DEP_REAC(1,absdep)
      NdomE = DEP_REAC(2,absdep)
      seg(i)%dom =  NdomO
      seg(i1e)%dom =  NdomE
      seg(i)%surface=  DEP_REAC(3,absdep)
      seg(i)%VarfreePlan=  DEP_REAC(4,absdep)
      domchange= DEP_REAC(5,absdep)

      if (kkdebug) then
        write(379,*) "Domain assignement using Barrier calculation"
        write(379,*) "Dom i", DEP_REAC(1,absdep)
        write(379,*) "Dom i1e", DEP_REAC(2,absdep)
        write(379,*) "free surface reaction type", DEP_REAC(3,absdep)
        write(379,*) "free surface", DEP_REAC(4,absdep)
        write(379,*) "domchange", DEP_REAC(5,absdep)
      endif
    endif

    if (seg(i)%voiso /= i1o) then
       ll=seg(i)%voiso
       do while (ll /= i1o)
          vo=seg(ll)%voiso
          seg(ll)%dom=NDomO
          ll=vo
       enddo
     endif

    if (seg(i)%voise /= i1e) then
        ll=seg(i)%voise
        do while (ll /= i1e)
          ve=seg(ll)%voise
          seg(ll)%dom = NDomE
          ll=ve
        enddo
    endif

  endif

  if (kkdebug) then
      write(379,*) " -----------------------------------------------------------------"
      write(379,*) " absdep definitive = ",absdep
      write(379,*) " longueurs o i e:",lonvoad,lonad,lonvead
      write(379,*) " longueurs Oip:",Oip,SEG(i1o)%O(1:3) + BVECLIN(1:3,SEG(i1o)%VECLIN)*lonvoad
      write(379,*) " longueurs Eip:",Eip,SEG(i1o)%O(1:3) + BVECLIN(1:3,SEG(i1o)%VECLIN)*lonvoad &
                                     + lonad * BVECLIN(1:3,seg(i)%veclin)
      write(379,*) "i1o,i,i1e", i1o,i,i1e
      write(379,*) " ----------------------------------------------------------------"
  endif

  ! If during the obstacle search I was modified and is now a segment to eliminate, hence no need to move it
  if ((i == seg(i)%voiso .and. i == seg(i)%voise) .or. out(i)) then
    if (kkdebug) write(379,*) "The segment i= ",i," is finally not moved and we cycle to the next !!!"
    cycle SURI
  endif

  !**************************************************
  !* Area swept by the segment and derived quantities
  !**************************************************
  ! Calculation of the area swept by the segment during its displacement
  ! As a convention: Positive area mean displacement in the direction imposed
  ! by the applied stress

#ifdef MDC

  if (SEG(i1o)%NORME==IZERO .and. lonvoad /= IZERO) seg(i1o)%SIGFE=half*(seg(i)%SIGFE(:)+seg(seg(i1o)%vnno)%SIGFE(:))
  if (SEG(i1e)%NORME==IZERO .and. lonvead /= IZERO) seg(i1e)%SIGFE=half*(seg(i)%SIGFE(:)+seg(seg(i1e)%vnne)%SIGFE(:))

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
if (Mon_Rang ==  IZERO) then
#endif

  if (ABSDEP /= IZERO ) then
    if (seg(i)%surface < IUN    .or.  &
        seg(i)%surface > IDEUX  .or.  &
       (seg(i)%surface > IZERO  .and. &
        seg(i)%surface < ITROIS .and. &
        domchange > IUN) ) then
      ! first element is reserved for storing nbswept at the end of the routine
      sweptsurfdata(1+16*nbswept+1)=int(ABSdep)
      sweptsurfdata(1+16*nbswept+2:1+16*nbswept+4)=int(Oi(1:3))
      sweptsurfdata(1+16*nbswept+5:1+16*nbswept+7)=int(Ei(1:3))
      sweptsurfdata(1+16*nbswept+8:1+16*nbswept+10)=int(VDOi(1:3)*ABSdep)
      sweptsurfdata(1+16*nbswept+11:1+16*nbswept+13)=int(VDEi(1:3)*ABSdep)
      !sweptsurfdata(1+16*nbswept+8:1+16*nbswept+10)=int(VDOi(1:3))
      !sweptsurfdata(1+16*nbswept+11:1+16*nbswept+13)=int(VDEi(1:3))
      sweptsurfdata(1+16*nbswept+14:1+16*nbswept+16)=int(bveclin(1:3, assoc(Li,1)))

!!!!!!!!!!!!
! Tests to be sure that we send a correct sweptsurface
!!!!!!!!!!!!
!      if ( ALL(Oi(1:3).EQ.Ei(1:3)) ) then
!        write(379,*) kk, " Oi == Ei in swepsurfdata", Oi,Ei
!        write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!        stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!      elseif (ALL(VDOi(1:3).EQ.(/ IZERO, IZERO, IZERO /))) then
!        write(379,*) kk, " VDOI is null vector", VDOi
!        write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!        stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!      elseif (ALL(VDEi(1:3).EQ.(/ IZERO, IZERO, IZERO /))) then
!        write(379,*) kk, " VDEI is null vector", VDEi
!        write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!        stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!      elseif (dot_product(VDOi,Sigdep*BVECDEP(1:3,Li)) < ZERO) then
!        write(379,*) kk, " Dot prod (VDOI,BVECDEP) is negative", dot_product(VDOi,BVECDEP(1:3,Li))
!        write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!        stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!      elseif (dot_product(VDEi,Sigdep*BVECDEP(1:3,Li)) < ZERO) then
!        write(379,*) kk, " Dot prod (VDEI,BVECDEP) is negative", dot_product(VDEi,BVECDEP(1:3,Li))
!        write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!        stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!      endif

      nbswept=nbswept+1
    endif
  endif

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon

endif
#endif

#endif

! TAG OJ
#ifdef OUTPUTSURF
  write(999,'(14i8)') ABSdep * sigdep , Li, Oi(:), Ei(:), nint(VDOi(:)), nint(VDEi(:))
#endif
! ENDTAG OJ

  ! Calculation of the area swept in unit a2
  daire(i) = 0.5 * float((lonad+longsegi) * ABSdep * SIGdep) * NORMLIN(Li) * NORMDEP(Li)

  ! In case of segments touching free surface some corrections are needed
  ! The deformation made by such segments is neglected
  if (seg(i)%surface < IUN .or. seg(i)%surface > ITROIS ) then

    ! This quantity is attribute to the slip systems activity
    itemp               = syseg(Li)
    AireSys(itemp)      = AireSys(itemp) + daire(I)
    AireSysInst(itemp)  = AireSysInst(itemp) + daire(I)

    ! The gammabox calculation
    if ((.not. allocation_dynamique_boites) .and. calculate_gammabox) then

      ! Box index for Oi, Ei, OiP, EiP
      BiO     = B3D_1D(modulo(floor(Oi(1)*invtailleBoite1),int(NboitesX))+iun,  &
                       modulo(floor(Oi(2)*invtailleBoite2),int(NboitesY))+iun,  &
                       modulo(floor(Oi(3)*invtailleBoite3),int(NboitesZ))+iun)
      BiE     = B3D_1D(modulo(floor(Ei(1)*invtailleBoite1),int(NboitesX))+iun,  &
                       modulo(floor(Ei(2)*invtailleBoite2),int(NboitesY))+iun,  &
                       modulo(floor(Ei(3)*invtailleBoite3),int(NboitesZ))+iun)
      BiPO    = B3D_1D(modulo(floor(OiP(1)*invtailleBoite1),int(NboitesX))+iun,  &
                       modulo(floor(OiP(2)*invtailleBoite2),int(NboitesY))+iun,  &
                       modulo(floor(OiP(3)*invtailleBoite3),int(NboitesZ))+iun)
      BiPE    = B3D_1D(modulo(floor(EiP(1)*invtailleBoite1),int(NboitesX))+iun,  &
                       modulo(floor(EiP(2)*invtailleBoite2),int(NboitesY))+iun,  &
                       modulo(floor(EiP(3)*invtailleBoite3),int(NboitesZ))+iun)

      ! The shear strain for each slip system is distributed in the gammabox tabs
      ! The swept area of oscillating segments must be shared in the gammabox(s) associated to the O and E side
      ! Also the swept area of oscillating segments must be shared between starting and ending coordinates
      daireBV = quart * daire(I) * BdivAboite

      Gammabox(BiO,itemp)  = Gammabox(BiO,itemp) + daireBV
      Gammabox(BiE,itemp)  = Gammabox(BiE,itemp) + daireBV
      Gammabox(BiPO,itemp) = Gammabox(BiPO,itemp) + daireBV
      Gammabox(BiPE,itemp) = Gammabox(BiPE,itemp) + daireBV

    endif

    ! Calculation for the polycrystal case
    if (desorientGrain) then
      AireSysGRAIN(seg(I)%grain,itemp) = AireSysGRAIN(seg(i)%grain,itemp) + daire(I)
    endif

    ! mobile dislocation density calculation
    NSEGMO      = NSEGMO+1
    RAUDMO_I    = RAUDMO_I + lonad*NORMLIN(Li)
    densitemob  = RAUDMO_I        ! utile pour le control en metallofute2

    ! Calcul of the mechanical work per step associate to the applied stress = tauApp*gamma
    ! attension: dans les simulation la force = tau*b*longueur(en a)/(mu*a)
    temp                = mudivV* daire(i)
    TrAppSysInst(itemp) = TrAppSysInst(itemp) + tauApP(i)*temp

    ! Calcul of the mechanical work per step associate to the internal stress
    TrIntSysInst(itemp) = TrIntSysInst(itemp) + (tauInt(I)+tauTL(I))*temp

    ! Usually, if a segment is improperly moved, his mobility law was not defined before
    if (nloi == 0) then
      write(*,*) "STOP,KK=",kk, " NLOI n'est pas defini dans contact ???"
      stop
    endif

    ! The effective displacement of each segments, per type of velocity law is accumulated.
    ! Displacements are multiply per the segments length for the calculation of a mean value per unit length
    if (.not. seg(i)%bloquer) then
      i1 = seg(i)%vnno
      i2 = seg(i)%vnne
      ! The length of vnno and vnne must be included otherwise the segment length is reduced after displacement
      if (i1 == nsegmax .and. i2 == nsegmax) then
        length3seg    = NORMLIN(Li)*seg(i)%norme
      elseif (i1 == nsegmax) then
        length3seg    = NORMLIN(Li)*seg(i)%norme + SEG(i2)%NORME*NormLin(seg(i2)%veclin)
      elseif (i2 == nsegmax) then
        length3seg    = NORMLIN(Li)*seg(i)%norme + SEG(i1)%NORME*NormLin(seg(i1)%veclin)
      else
        length3seg    = NORMLIN(Li)*seg(i)%norme + SEG(i1)%NORME*NormLin(seg(i1)%veclin) + &
                        SEG(i2)%NORME*NormLin(seg(i2)%veclin)
      endif

      Lseg_V(nloi)  = Lseg_V(nloi) + length3seg
      D_Vraie(nloi) = D_Vraie(nloi) + (absdep * length3seg)
    endif

    ! A paranoid test to check that anglevis is always defined
    if (seg(i)%anglevis < zero) then
      print*,'A problem is found in the calculation of anglevis: 12contact',kk,i,seg(i)%anglevis
      call seginfo(i,"Pb anglevis 1")
      stop
    endif

    ! Calculation of the area swept by screw, edge and mixed dislocations
    if (desorientGrain) then
      if (abs(seg(i)%anglevis) < Critic_angle) then
        airevisSys(itemp) = airevisSys(itemp) + daire(i)*sign(1.D0,SchmidSysGrain(seg(i)%grain,itemp))
      elseif (abs(halfpii-seg(i)%anglevis) < Critic_angle) then
        airecoinSys(itemp) = airecoinSys(itemp) + daire(i)*sign(1.D0,SchmidSysGrain(seg(i)%grain,itemp))
      else
        AireMixtSys(itemp) = AireMixtSys(itemp) + daire(i)*sign(1.D0,SchmidSysGrain(seg(i)%grain,itemp))
      endif
    else
      if (abs(seg(i)%anglevis) < Critic_angle) then
!         print*,'vis',seg(i)%anglevis,daire(i)*sign(1.D0,SchmidSys(itemp))
        airevisSys(itemp) = airevisSys(itemp) + daire(i)*sign(1.D0,SchmidSys(itemp))
      elseif (abs(halfpii-seg(i)%anglevis) < Critic_angle) then
!         print*,'coin',seg(i)%anglevis,daire(i)*sign(1.D0,SchmidSys(itemp))
        airecoinSys(itemp) = airecoinSys(itemp) + daire(i)*sign(1.D0,SchmidSys(itemp))
      else
!         print*,'mixt',seg(i)%anglevis,daire(i)*sign(1.D0,SchmidSys(itemp))
        AireMixtSys(itemp) = AireMixtSys(itemp) + daire(i)*sign(1.D0,SchmidSys(itemp))
      endif
    endif

  endif

  !*************************************
  ! End of the swept area calculation  *
  !*************************************

  !****************************************************************
  ! Effective displacement of the segment with or without reactions
  !****************************************************************

  !***********************************************************************************
  ! The recombination cases, i.e. the segment is recombined with one of its neighbors
  ! after displacement. The is the opposite of the discretisation procedure
  !***********************************************************************************

  ! No need to recombine segments if some annihilation is coming latter
  notannili2o = .not. (anseg .and. obsta == i2o)
  notannili2e = .not. (anseg .and. obsta == i2e)

  !***************************************
  !*** RECOMBINAISON i2o i1o I i1e i2e ***
  !***************************************

  ! We first store values of I neighboring in case something is going wrong with the segment connection
  ! Segments are listed from the vnno to the vnne

  ! The starting segment
  if (seg(i)%vnno == Nsegmax) then
    itmp = seg(i)%voiso             ! The special case of pinning points
  else
    itmp = seg(i)%vnno
  endif

  temp1 = SEG(i2O)%NORME + lonad + SEG(i2E)%NORME   ! TOTAL LENGTH after combination

  ! test for the i2e and i2o alignement
  If (i2o == NSEGMAX .or. i2e == NSEGMAX .or. seg(i)%surface > ideux) then
  condition=.False.
  else
    condition = (lonvoad == zero                                    .and. &
                lonvead == zero                                     .and. &
                Li == SEG(i2o)%VECLIN                               .and. &
                Li == SEG(i2e)%VECLIN                               .and. &
                .not.(SEG(i2o)%JONC .and. SEG(i2e)%JONC)            .and. &
                .not.(SEG(i2o)%JONC .and. anseg)                    .and. &
                .not.(SEG(i2e)%JONC .and. anseg)                    .and. &
                .not. SEG(i1o)%JONC                                 .and. &
                .not. SEG(i1e)%JONC                                 .and. &
                .not.(SEG(i2o)%JONC .and. SEG(i2e)%gd > izero)      .and. &
                .not.(SEG(i2e)%JONC .and. SEG(i2o)%gd > izero)      .and. &
                i2o /= i2e                                          .and. &
                seg(i2e)%vnne /= i2o                                .and. &
                temp1 < Lomai)
   endif

   recomb=.FALSE.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test on segments recombination after displacement !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (condition) then       ! i2o-i-i2e can be recombine to form one segment

    if (kkdebug) then
      write(379,*) "--------------------------------------------------------------------------"
      write(379,*) " Extended recombination i2o --> i2e : ",i2o,i2e
      write(379,*) " the length of i2o is extended to the end of i2e "
      write(379,*) "--------------------------------------------------------------------------"
      call conf(i)
      call seginfo(i,"Extended recombination i2o --> i2e")
    endif

    !* Gestion des mouvements ulterieurs
    idep(i2o)         = izero  ! If one segment was participating to a junction idep must be zero
    seg(i2O)%RESDEP   = zero
#ifdef MDC
    seg(i2o)%SIGFE(:) = seg(I)%SIGFE(:)
#endif


    ! if needed we exchange junction between i2e and i2o
    if (seg(i2E)%JONC) then

      temp2 = seg(i2E)%iJONC

      seg(i2E)%JONC     = .false.
      seg(i2E)%IJONC    = nsegmax
      seg(i2E)%tJONC    = 0

      seg(i2O)%JONC     = .true.
      seg(i2O)%IJONC    = temp2
      seg(i2O)%tJONC    = seg(temp2)%tJONC

      seg(temp2)%ijonc  = i2o

    endif

    ! The new i2o length
    seg(i2O)%norme = seg(i2O)%norme + lonad + seg(i2E)%norme

    ! if needed we exchange surface information between i2e and i2o
    if (seg(i2e)%surface == IDEUX)  then

      seg(i2o)%surface = IDEUX
      seg(i2o)%VarFreePlan = seg(i2e)%VarFreePlan
      seg(i2e)%surface = IZERO
      seg(i2e)%VarFreePlan = IZERO

    endif

    ! i1o i i1e i2e must be eliminated
    seg(i1o)%norme  = izero
    seg(I)%norme    = izero
    seg(i1e)%norme  = izero
    seg(i2e)%norme  = izero
    IDEP(i2e)       = izero

    itemp = seg(i2o)%vnno
    ! in all cases, we reconnecte itemp to i2e
    GD2made = .false.
    call ConnecIJ(itemp,i2o,22)

    itemp = seg(i2e)%vnne
    ! in all cases, we reconnecte itemp to i2e
    GD2made = .false.
    call ConnecIJ(i2o,itemp,23)

    !* Gestion des Voisins Non Nuls
    !* i devient i2o pour reaction eventuelle juste apres
    Loiseg(I2o) = Loiseg(I)
    ! attention: the indice of I is modified for future reaction
    i = i2o

    LONGSEGi = SEG(i)%NORME

    recomb=.TRUE.

  else

    !***********************************************************
    !* simple displacement or recombination on the o or E side *
    !***********************************************************

    simpleDisp=.true.

    ! For recombination involving segment I touching a free surface: lonad must be correct
    if ( seg(i)%surface > IZERO .and. seg(i)%surface < ITROIS ) then

      simpleDisp=.false.

      if (domchange < IDEUX ) then

        call resolve_freeSurface_hooking_up(i,i1o,i1e,lonad,lonvoad,lonvead,SIGdep,absdep,3670)
        if (domchange == IUN) then
           if (seg(i)%surface == IUN .and. (lonvead == IZERO .or. DVNNE == IZERO)) call voisinage(i1e,1021)
           if (seg(i)%surface == IDEUX .and. (lonvoad == IZERO .or. DVNNO == IZERO)) call voisinage(i1o,1021)
        endif
#ifdef MDC
        if (ABSDEP /= IZERO ) then

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
if (Mon_Rang ==  IZERO) then
#endif

          if (kkdebug) then
            write(379,*) "Swept surface correction for segment touching a surface after free_surface_hooking_up "
          endif
          ! first element is reserved for storing nbswept at the end of the routine
          Oip(:) = modulo(seg(i)%o(:) + Tr(:),modur(:)) - Tr(:)
          Eip(:) = Oip(:) + lonad*BVECLIN(1:3,Li)
          sweptsurfdata(1+16*nbswept+1)=int(ABSdep)
          sweptsurfdata(1+16*nbswept+2:1+16*nbswept+4)=int(Oi(1:3))
          sweptsurfdata(1+16*nbswept+5:1+16*nbswept+7)=int(Ei(1:3))
          sweptsurfdata(1+16*nbswept+8:1+16*nbswept+10)=int(Oip(:) - Oi(:))
          sweptsurfdata(1+16*nbswept+11:1+16*nbswept+13)=int(Eip(:) - Ei(:))
          !sweptsurfdata(1+16*nbswept+8:1+16*nbswept+10)=int(VDOi(1:3)*ABSdep)
          !sweptsurfdata(1+16*nbswept+11:1+16*nbswept+13)=int(VDEi(1:3)*ABSdep)
          sweptsurfdata(1+16*nbswept+14:1+16*nbswept+16)=int(bveclin(1:3, assoc(Li,1)))

!!!!!!!!!!!!
! Tests to be sure that we send a correct sweptsurface
!!!!!!!!!!!!
!
!          if ( ALL(Oi(1:3).EQ.Ei(1:3)) ) then
!            write(379,*) kk, " Oi == Ei in swepsurfdata", Oi,Ei
!            write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!            stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!          elseif (ALL((Oip(1:3) - Oi(1:3)).EQ.(/ IZERO, IZERO, IZERO /))) then
!            write(379,*) kk, " VecOI is null vector", (Oip(1:3) - Oi(1:3))
!            write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!            stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!          elseif (ALL((Eip(1:3) - Ei(1:3)).EQ.(/ IZERO, IZERO, IZERO /))) then
!            write(379,*) kk, " VecEI is null vector", (Eip(1:3) - Ei(1:3))
!            write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!            stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!          elseif (dot_product((Oip(1:3) - Oi(1:3)),Sigdep*BVECDEP(1:3,Li)) < ZERO) then
!            write(379,*) kk, " Dot prod (VDOI,SIGNBVECDEP) is negative", dot_product((Oip(1:3) - Oi(1:3)),Sigdep*BVECDEP(1:3,Li))
!            write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!            stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!          elseif (dot_product((Eip(1:3) - Ei(1:3)),Sigdep*BVECDEP(1:3,Li)) < ZERO) then
!            write(379,*) kk, " Dot prod (VDEI,SIGNBVECDEP) is negative", dot_product((Eip(1:3) - Ei(1:3)),Sigdep*BVECDEP(1:3,Li))
!            write(379,*) kk, nbswept, sweptsurfdata(1+16*nbswept+1:1+16*nbswept+16)
!            stop "Problem in sweptsurface sent to Zebulon : see debug file"
!
!
!          endif

          nbswept=nbswept+1
        endif

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon

endif
#endif

#endif


      endif

      call correct_freeSurface_connection(i,i1o,i1e,lonad,lonvoad,lonvead,domchange,absdep)

      if (seg(i)%surface == 1) seg(i1o)%veclin=conec(seg(i)%veclin,1)
      if (seg(i)%surface == 2) seg(i1e)%veclin=conec(seg(i)%veclin,2)

      ! in case that DOMCHANGE is larger than ZERO, we cycle. The displacement
      ! has been applied and the connectivity has been updated in the correct_freeSurface_connection subroutine
      if (DOMCHANGE > IZERO) cycle surI



    endif

    !*******************************
    !*** Recombination i2o i1o I ***
    !*******************************

    temp2 = seg(i2o)%norme
    temp1 = lonad + temp2

    condition=.False.
    If (i2O /= nsegmax .and. domchange < IUN) then
      condition = ( lonvoad == zero                     .and. &
                  Li == SEG(i2O)%VECLIN                 .and. &
                  i2e /= i2o                            .and. &
                  .not. seg(i1o)%JONC                   .and. &
                  SEG(i1o)%gd < iun                     .and. &
                  .not.(seg(i2O)%jonc .and. (anseg .or. seg(i)%surface > IZERO)) .and. &
                  seg(i2o)%gd < iun                     .and. &
                  .not. seg(i)%JONC                     .and. &
                  seg(i2o)%vnno /= i2e                  .and. &
                  seg(i)%surface < itrois               .and. &
                  seg(i2o)%surface < itrois             .and. &
                  temp1 < lomai )
    endif

    if (condition) then

      if (kkdebug) then
          write(379,*)  "--------------------------------------------------------------------------"
          write(379,*) " Recombinaison type : i2o,i1o,i",i2o,i1o,i
          write(379,*) " We extend I2o until end of segment i (after displacement)"
          write(379,*) "--------------------------------------------------------------------------"
          call seginfo(i,"Recombinaison type : i2o,i1o,i")
      endif

      ! c'est beaucoup plus facile d'allonger i2o cote E que de faire l'inverse
      SEG(i2o)%NORME = temp1
      Ei(:) = Seg(i2o)%O(:) +  temp1 * bveclin(:,Li)
      Seg(I1o)%O(:) = Ei(:)
      Seg(i1o)%NORME = Izero
      Seg(I)%O(:) = Ei(:)
      Seg(i)%NORME = Izero
      Seg(I1e)%O(:) = Ei(:)
      if ( seg(i)%surface == IDEUX) then
        lonvead               = IZERO
        seg(i2o)%surface      = seg(i)%surface
        seg(i2o)%VarFreePlan  = seg(i)%VarFreePlan
        singulier(i2o)        = .true.
        singulier(i)          = .false.
        seg(i)%surface        = IZERO
        seg(i)%VarFreePlan    = IZERO
      endif

#ifdef MDC
      seg(i2o)%SIGFE(:)=seg(i1o)%SIGFE(:)
#endif
      ! VOISO voit sa longueur modifiee
      ! Ltemp:true if I1o was not considered as kneecap
      Ltemp = ((SEG(i1e)%NORME /= Izero) .or. seg(i1e)%jonc .or. seg(i1e)%gd > izero)
      SEG(i1e)%NORME = lonvead
      ! i1e etait rotule : il faut metre a jour ses voisins cote O,(cote E se fera plus tard)
      if(.not. Ltemp) call Voisinage(i1e,4)

      Ltemp = ((lonvead /= Izero) .or. seg(i1e)%jonc .or. seg(i1e)%gd > izero)
      ! temp1 is neighbor nn to be connected to i2e
      if(Ltemp) then
        temp1 = I1e
      else
        temp1 = seg(i1e)%vnne
        if (seg(i1e)%surface /= IZERO) then
          seg(i2o)%surface      = IDEUX
          seg(i2o)%VarFreePlan  = seg(i1e)%VarFreePlan
          seg(i1e)%surface      = IZERO
          seg(i1e)%VarFreePlan  = IZERO
        endif
      endif

      Ltemp = ((temp2 /= Izero) .or. seg(i2o)%jonc .or. seg(i2o)%gd > izero)
      if(.not. Ltemp)   call voisinage(i2o,24)

      GD2made = .false.
      call ConnecIJ(i2o,temp1,24)
      ! for future reactions
      i = i2o
      LONGSEGi = SEG(i)%NORME
      seg(i)%resdep = zero ; idep(i)= IZERO
      recomb=.TRUE.

    else

      !*******************************
      !*** Recombination i2e i1e I ***
      !*******************************

      temp2 = seg(i2e)%norme
      temp1 = lonad + temp2

      condition = .false.
      If (i2E /= NSEGMAX .and. domchange < IUN) then
        condition = ( lonvead == zero                        .and. &
                    Li == SEG(i2E)%VECLIN                  .and. &
                    i2o /= i2e                             .and. &
                    .not. SEG(i1e)%JONC                    .and. &
                    .not.(SEG(i2E)%jonc.and.(anseg  .or. seg(i)%surface > IZERO))  .and. &
                    seg(i2e)%gd < iun                      .and. &
                    .not. SEG(i)%JONC                      .and. &
                    seg(i2e)%vnne /= i2o                   .and. &
                    seg(i)%surface < itrois                .and. &
                    seg(i2e)%surface < itrois              .and. &
                    temp1 < lomai )
      endif

      if (condition) then

        if (kkdebug) then
          write(379,*) "--------------------------------------------------------------------------"
          write(379,*) "Recombinaison type i,i1e,i2e",i,i1e,i2e
          write(379,*) "  On recule l'origine de I2e "
          write(379,*) "--------------------------------------------------------------------------"
          call seginfo(i,"Recombinaison type i,i1e,i2e")
          write(379,*) "lonad",lonad,"lonvoad", lonvoad
        endif

        ! c'est beaucoup plus facile d'allonger i2e cote O que de faire l'inverse
        Oi(:) = Seg(i2E)%O(:) - lonad * bveclin(:,Li)
        Seg(I2E)%O(:) = Oi(:)
        SEG(i2E)%NORME = temp2 + lonad
        SEG(i2E)%resdep =  zero
        idep(i2E) =  Izero
        seg(i2e)%dom = seg(i)%dom
#ifdef MDC
        seg(i2e)%SIGFE(:)=seg(i1e)%SIGFE(:)
#endif
        SEG(i1E)%NORME = izero
        SEG(i)%NORME  = izero
        if (seg(i)%surface /= IZERO) then
            seg(i2e)%surface=seg(i)%surface
            seg(i2e)%VarFreePlan=seg(i)%VarFreePlan
            if (seg(i)%surface==IUN) seg(i2e)%dom=seg(i)%dom !transfers of domains if surface segment with origin
            singulier(i2e)=.true.
            singulier(i)=.false.
            seg(i)%surface = IZERO
            seg(i)%VarFreePlan=IZERO
        end if
        if ( seg(i2e)%surface == IUN) lonvoad=IZERO

        ! VOISO voit sa longueur modifiee
        ! Ltemp:true if I1o was not considered as kneecap
        Ltemp = ((SEG(i1o)%NORME /= Izero) .or. seg(i1o)%jonc .or. seg(i1o)%gd > izero)
        SEG(i1o)%NORME=lonvoad
        ! i1o etait rotule : il faut metre a jour ses voisins cote O,(cote E se fera plus tard)
        if(.not. Ltemp) call Voisinage(i1o,2)

        Ltemp = ((SEG(i1o)%NORME /= Izero) .or. seg(i1o)%jonc .or. seg(i1o)%gd > izero)
        ! temp1 is neighbor nn to be connected to i2e
        if(Ltemp) then
          temp1 = I1o
        else
          temp1 = seg(i1o)%vnno
          if (seg(i1o)%surface /= IZERO) then
            seg(i2e)%surface= IUN
            seg(i2e)%VarFreePlan=seg(i1o)%VarFreePlan
            seg(i2e)%dom=seg(i1o)%dom !transfers of domains if surface segment
            seg(i1o)%surface=IZERO
            seg(i1o)%VarFreePlan=IZERO
          endif
        endif

        Ltemp = ((temp2 /= Izero) .or. seg(i2e)%jonc .or. seg(i2e)%gd > izero)

        if(.not. Ltemp)   call voisinage(i2e,23)

        GD2made = .false.
        call ConnecIJ(temp1,i2e,28)
        ! for future reactions
        i = i2e
        LONGSEGi = SEG(i)%NORME
        idep(i)=0
        recomb=.TRUE.

      else

        !*************************************************************
        ! Simple displacement without recombination between segments *
        !*************************************************************

        if (kkdebug) then
          write(379,*) "--------------------------------------------------------------------------"
          write(379,*) " Simple displacement i1o,i,i1e ",i1o,i,i1e
          write(379,*) "--------------------------------------------------------------------------"
          call seginfo(i," Simple displacement i1o,i,i1e ")
        endif

        if (simpleDisp) then

          londiff         = lonvoad - seg(i1o)%norme
          seg(i1o)%norme  = lonvoad

          ! Be careful that coordinates of Oi must not be calculated from the shift of the end of i1o as
          ! this point may be the coordinate Oi modulo the simulation box. Problems with the obstacle
          ! distance appears if we don t pay attention to this!
          Oi(:)           = Oi(:) + BVECLIN(:,seg(i1o)%veclin) * londiff
          seg(I)%o(:)     = Oi(:)
          seg(I)%norme    = lonad
          Ei(:)           = Oi(:) + lonad * BVECLIN(:,Li)
          seg(i1e)%o(:)   = Ei(:)
          seg(i1e)%norme  = lonvead

          if(kkdebug) write(379,fmt='(" Oip = ",3I9," Eip = ",3I9," AbsDEp = ",I4,  &
                            & " lonvoad = ",I4," lonvead = ", I4," lonad = ",I4)')  &
                              Oi,Ei,AbsDep,lonvoad,lonvead,lonad

          ! After displacement the voiso and voise of i can be of zero length and we must reconstruct
          ! the connectivity info. Better redefining i1o and i1e since elimination of pivot segment
          ! loops may have modified such information
          voisi1o = seg(i)%voiso
          thecounter = IZERO
          do while (seg(voisi1o)%norme == izero .and. voisi1o /=I1o)
            thecounter=thecounter + 1
            if (thecounter > ISEPT) then
              write(379,*)  i1o , i ,i1e
              call seginfo(i,'the counter voiso')
              stop
            endif
            if (voisi1o /= i1o) seg(voisi1o)%O = seg(i)%O
            call voisinage(voisi1o,541)
            if (out(i)) cycle SurI   ! subroutine voisinage created a small loop to eliminate
                                     ! Nothing more to do
            voisi1o = seg(voisi1o)%voiso
          enddo
          ! subroutine voisinage created a small loop to eliminate
          ! Nothing more to do

          voisi1e = seg(i)%voise
          thecounter=IZERO
          do while (seg(voisi1e)%norme == izero .and. voisi1E /=i1e)
            thecounter=thecounter + 1
            if (thecounter > ISEPT) then
              write(379,*)  i1o , i ,i1e
              call seginfo(i,'the counter voiso')
              stop
            endif
            if (voisi1e /= i1e) seg(voisi1e)%O = seg(i1e)%O
            call voisinage(voisi1e,542)
            if (out(i)) cycle SurI  ! subroutine voisinage created a small loop to eliminate
                                    ! Nothing more to do
            voisi1e=seg(voisi1e)%voise
          enddo


          ! Configuration touching the surface must be corrected
          if (seg(i)%surface>ITROIS) then
            call surface_correct(i,i1o,i1e)
            surf_correct = .true.
          endif

        else

          !special displacement for segment touching the surfaces
          Oi(1:3) = SEG(I)%O(1:3)
          Ei(:)   = Oi(:) + lonad * BVECLIN(1:3,Li)
          if (seg(I)%VarFreePlan < IUN) seg(i)%surface = IZERO

        endif

        if (seg(i)%surface == ITROIS) then
          print *, 'surface segments of type 3 are not expected to move'
          stop
        elseif (seg(i)%surface > ITROIS) then
          print *, 'surface segments of type > 3 are not expected to exist here'
          stop
        endif

        ! If kneecaps exist between i et i1o, they must be moved
        itemp       = seg(i)%voiso
        itemp_loop  = itemp           ! mark the beginning of a possible loop
        itemp_mark  = .false.

        if (seg(i)%surface /= IUN .and. .not. OUT(I)) then

          Do while (seg(itemp)%norme == izero .and. itemp /= nsegmax)
            if (itemp_mark .and. itemp==itemp_loop) exit    !it s a closed loop made of kneecaps
            itemp_mark=.true.
            SEG(itemp)%O = Oi
            if (itemp == seg(itemp)%voiso) then
              print *, "STOP,KK=",kk, "contact: erreur de voisinage en O, KK = ",KK
              exit
            endif
            itemp = seg(itemp)%voiso
          enddo

        end if

        ! If kneecaps exist between i1e et i, they must be moved
        itemp = seg(i)%voise
        itemp_loop=itemp   ! mark the begining of a possible loop
        itemp_mark=.false.

        if (seg(i)%surface /= IDEUX .and. .not. OUT(I)) then
          Do while (seg(itemp)%norme == izero .and. itemp /= nsegmax)
            if (itemp_mark .and. itemp==itemp_loop) exit    !it s a closed loop made of kneecaps
            itemp_mark=.true.
            SEG(itemp)%O = Ei
            if (itemp == seg(itemp)%voise ) then
                print *, "STOP,KK=",kk, "contact: erreur de voisinage en E, KK ",KK
                exit
            endif
            itemp = seg(itemp)%voise
          enddo
        end if

        vnn2e=seg(i1e)%vnne
        vnn2o=seg(i1o)%vnno

        ! special case i1e = GD, seg(i1e)%vnne = screw segment with length > 0
        ! and no kneecap between i and i1e: superposition of two screw segments,
        ! annihilation of the superposed screw length
        condition = .false.
        If (seg(i1e)%gd > izero .and. seg(i1e)%surface < IUN) then
          condition = ( seg(i1e)%norme > izero                                  .and. &
                        TYSEG(seg(vnn2e)%veclin) == iun                         .and. &
                        deassoc(seg(vnn2e)%veclin) /= deassoc(seg(i1e)%veclin)  .and. &
                        (seg(i)%veclin == conec(seg(I1e)%veclin,1) .or.               &
                        seg(i)%veclin == conec(seg(I1e)%veclin,2))              .and. &
                        VDEi(1) == bveclin(1,seg(vnn2e)%veclin)                 .and. &
                        VDEi(2) == bveclin(2,seg(vnn2e)%veclin)                 .and. &
                        VDEi(3) == bveclin(3,seg(vnn2e)%veclin)                 .and. &
                        seg(vnn2e)%surface < IUN                                      &
                      )
        endif

        if (condition) then

          if (kkdebug) call seginfo(i,"before GD screw correction ")

          if (seg(i1e)%norme < seg(vnn2e)%norme) then

            seg(vnn2e)%O      = seg(i1e)%O
            seg(vnn2e)%norme  = seg(vnn2e)%norme-seg(i1e)%norme
            seg(vnn2e)%dom    = seg(i1e)%dom
            seg(i1e)%norme    = izero
            GD2made = .false.
            call connecIJ(I1e,vnn2e,178)

            if (kkdebug) call seginfo(i,"GD screw displacement I1E 1")

          elseif (seg(i1E)%norme == seg(vnn2e)%norme) then

            seg(vnn2e)%O      = seg(i1e)%O
            seg(vnn2e)%dom    = seg(i1e)%dom
            seg(vnn2e)%norme  = izero
            seg(i1e)%norme    = izero
            vnn2e             = seg(vnn2e)%vnne
            GD2made = .false.
            call connecIJ(I1e,vnn2e,179)
            call voisinage(vnn2e,1791)

            if (kkdebug) call seginfo(i,"GD screw displacement I1E 2")

          else

            seg(i1e)%norme    = seg(i1e)%norme-seg(vnn2e)%norme
            seg(vnn2e)%O      = seg(i1e)%O+seg(i1e)%norme*bveclin(:,seg(I1e)%veclin)
            seg(vnn2e)%norme  = izero
            seg(vnn2e)%dom    = seg(seg(vnn2e)%voise)%dom
            vnn2e             = seg(vnn2e)%vnne
            GD2made = .false.
            call connecIJ(I1e,vnn2e,180)
            call voisinage(vnn2e,1801)

            if (kkdebug) call seginfo(i,"GD screw displacement I1E 3")

          endif

          if (out(i)) cycle SurI    ! subroutine voisinage created a small loop to eliminate
                                    ! Nothing more to do

        endif

        ! special case i1o = GD, seg(i1o)%vnno = screw segment with length > 0
        ! and no rotule between i and i1o: superposition of a two screw segments,
        ! annihilation of the superposed screw length

        condition = .false.
        if (seg(i1o)%gd > izero .and. seg(i1o)%surface < IUN) then
          condition = ( seg(i1o)%norme > izero                                  .and. &
                        TYSEG(seg(vnn2o)%veclin) == iun                         .and. &
                        deassoc(seg(vnn2o)%veclin) /= deassoc(seg(i1o)%veclin)  .and. &
                        (seg(i)%veclin == conec(seg(I1o)%veclin,1) .or.               &
                        seg(i)%veclin == conec(seg(I1o)%veclin,2))              .and. &
                        VDOi(1) == -bveclin(1,seg(vnn2o)%veclin)                .and. &
                        VDOi(2) == -bveclin(2,seg(vnn2o)%veclin)                .and. &
                        VDOi(3) == -bveclin(3,seg(vnn2o)%veclin)                .and. &
                        seg(vnn2o)%surface < IUN                                      &
                      )
        endif

        if (condition) then

          if (kkdebug) call seginfo(i,"before GD screw correction ")

          if(seg(i1o)%norme < seg(vnn2o)%norme) then

            seg(seg(i1o)%vnno)%norme = seg(vnn2o)%norme-seg(i1o)%norme
            seg(i1o)%norme           = IZERO
            seg(I1o)%O               = seg(i)%O
            seg(I1o)%dom             = seg(i)%dom
            GD2made = .false.
            call connecIJ(vnn2o,I1o,181)

            if (kkdebug) call seginfo(i,"GD screw displacement I1O 1")

          elseif(seg(i1o)%norme == seg(vnn2o)%norme) then

            seg(vnn2o)%norme = IZERO
            seg(i1o)%norme   = IZERO
            seg(I1o)%O       = seg(i)%O
            seg(I1o)%dom     = seg(i)%dom
            vnn2o            = seg(vnn2o)%vnno
            GD2made = .false.
            call connecIJ(vnn2o,I1o,182)
            call voisinage(vnn2o,1821)

            if (kkdebug) call seginfo(i,"GD screw displacement I1O 2")

          else

            seg(i1o)%O       = seg(vnn2o)%O
            seg(i1o)%norme   = seg(i1o)%norme-seg(vnn2o)%norme
            seg(vnn2o)%norme = IZERO
            seg(i1o)%dom     = seg(vnn2o)%dom
            vnn2o            = seg(vnn2o)%vnno
            GD2made = .false.
            call connecIJ(vnn2o,i1o,183)
            call voisinage(vnn2o,1831)

            if (kkdebug) call seginfo(i,"GD screw displacement I1O 3")

          endif

          if (out(i)) cycle SurI    ! subroutine voisinage created a small loop to eliminate

        endif

        ! The variable I1Orotule(1) = .false. If i1o is not of zero length before moving
        ! The variable I1Orotule(2) = .false. if i1o is not of zero length after moving
        ! si il y a un changement , il faut refaire la connectivite
        I1Orotule(1) = (seg(i1o)%gd < iun .and. .not.seg(i1o)%jonc .and. DVNNO == 0 )
        I1Orotule(2) = (seg(i1o)%gd < iun .and. .not.seg(i1o)%jonc .and. lonvoad == 0 )

        if (kkdebug) write(379,*) "I1Orotule", I1Orotule

        ! The variable I1Erotule(1) = .false. if i1E is not of zero length before moving
        ! The variable I1Erotule(2) = .false. if i1E is not of zero length after moving
        ! si il y a un changement , il faut refaire la connectivite
        I1Erotule(1) = (seg(i1E)%gd < iun .and. .not.seg(i1E)%jonc .and. DVNNE == 0 )
        I1Erotule(2) = (seg(i1E)%gd < iun .and. .not.seg(i1E)%jonc .and. lonvEad == 0 )

        if (kkdebug) write(379,*) "I1Erotule", I1Erotule

        ! On intervient dans la connectivite si la condition suivante n'est pas verifie
        condition = ((I1Orotule(1) .eqv. I1Orotule(2))        .and. &
                      (I1Erotule(1) .eqv. I1Erotule(2))        .and. &
                      (lonad /= 0))

        if (kkdebug) write(379,*) "condition", condition

        ! si seg(i1o)%vnno == i1e  COLLAPSE
        if(.not. condition .and. (seg(i1o)%vnno /= i1e) .and. domchange < IUN ) then
          if (kkdebug) then
            write(379,*) 'OUT(i1o)=',OUT(i1o),'i1o',i1o
            write(379,*) 'OUT(i1e)=',OUT(i1e),'i1e',i1e
            write(379,*) 'OUT(i)=',OUT(i),'i',i
          endif

          if((I1orotule(1) .neqv.  I1orotule(2)) .and. .not. (OUT(i1o))) call voisinage(i1o,14)
          if(lonad == izero .and. .not. (OUT(I)))                        call voisinage(i,5)
          if((I1erotule(1) .neqv.  I1erotule(2)) .and. .not. (OUT(i1e))) call voisinage(i1e,4)

          if (kkdebug) call seginfo(i,"voisinage")

          temp1 = SEG(i)%voiso
          vnn1o = SEG(i)%vnno
          temp2 = SEG(i)%voise
          vnn1e = SEG(i)%vnne

          ! if I becomes a kneecap, one should redo the connection
          if (LONAD == IZERO .and. .not. OUT(i)) then
            if (seg(vnn1o)%vnno == vnn1e .and. seg(vnn1o)%vnno /= nsegmax) then
               call desbou(i)
               cycle
            endif
            if (seg(vnn1e)%surface == IZERO .and. seg(vnn1o)%surface == IZERO)  then
              GD2made = .false.
              call connecIJ(vnn1o,vnn1e,34)
              if(kkdebug) write(379,*) " I=",I,"  disparait apres deplacement, CYCLE SUR I"
            endif
            cycle
          else
            ! determination of the number of kneecaps
            Ivnn1o = IUN
            do while (temp1 /= vnn1o .and. .not. OUT(i))
              Ivnn1o = Ivnn1o + IUN
              if(temp1 == seg(temp1)%voiso) then
                print *, " fatal error : contact: infinite loop 11 at step kk = ",KK
                stop
              endif
              temp1 = seg(temp1)%voiso
            enddo

            Ivnn1e = IUN
            do while (temp2 /= vnn1e .and. .not. OUT(i))
              Ivnn1e = Ivnn1e + IUN
              if(temp2 == seg(temp2)%voise) then
                print *, " fatal error : contact: infinite loop 11 bis at step k= ",KK
                stop
              endif
              temp2 = seg(temp2)%voise
            enddo

            ! supression of accumulated pivot segments, if no annihilation
            !  The connectivity is tested on the O side of i
            if (vnn1o /= nsegmax) then
              temp1 = seg(vnn1o)%veclin
              temp2 = Li
              if (Ivnn1o - 1 > nbrot(temp1,temp2,1) .and. notannili2o .and. .not.OUT(i)) then
                GD2made = .false.
                call connecIJ(vnn1o,i,34)
              endif
            else
              ! if we touch a pinning point, we systematically reconstruct the connectivity to eliminate any accumulated pivot segment
              if(.not.OUT(i) .and. seg(i)%surface == IZERO .and. .not. surf_correct) then
                GD2made = .false.
                call connecIJ(vnn1o,i,3434)
              endif
            endif

            !  The connectivity is tested on the E side of i
            if (vnn1e /= nsegmax) then
              temp1 = seg(vnn1e)%veclin
              temp2 = Li
              if (Ivnn1e - 1 > nbrot(temp2,temp1,1) .and. notannili2o .and. .not.OUT(i)) then
                GD2made = .false.
                call connecIJ(i,vnn1e,35)
              endif
            else
              ! if we touch a pinning point, we systematically reconstruct the connectivity to eliminate any accumulated pivot segment
              if(.not.OUT(i) .and. seg(i)%surface == IZERO .and. .not. surf_correct) then
                GD2made = .false.
                call connecIJ(i,vnn1e,3535)
              endif
            endif

          endif

        elseif (seg(i1o)%vnno == i1e .and. i1e /= NSEGMAX) then

          ! Small loops annihilation
          if (kkdebug) call seginfo(I,"A small closed loop is eliminated starting from I")
          call desbou(i)
          jnct=.false.
          if (kkdebug) call seginfo(I," apres desbou          ")

        endif

      endif     ! End of the recombination test at the end

    endif      ! End of the recombination test at the origin

  endif       ! End of the recombination test

  if (recomb) then

    vnn2e=seg(i1e)%vnne
    vnn2o=seg(i1o)%vnno

    !special case i1e = GD, seg(i1e)%vnne = screw segment with length > 0
    !and no rotule between i and i1e: superposition of a two screw segments,
    !annihilation of the superposed screw length
    condition = .false.
    if (seg(i1e)%gd > izero .and. seg(i1O)%surface < IUN) then
      condition = ( seg(i1e)%norme > IZERO                                  .and. &
                    TYSEG(seg(vnn2e)%veclin)==IUN                           .and. &
                    deassoc(seg(vnn2e)%veclin) /= deassoc(seg(I1e)%veclin)  .and. &
                    (seg(i)%veclin == conec(seg(I1e)%veclin,1) .or.               &
                    seg(i)%veclin == conec(seg(I1e)%veclin,2))              .and. &
                    VDEi(1) == bveclin(1,seg(vnn2e)%veclin)                 .and. &
                    VDEi(2) == bveclin(2,seg(vnn2e)%veclin)                 .and. &
                    VDEi(3) == bveclin(3,seg(vnn2e)%veclin)                 .and. &
                    seg(vnn2o)%surface < IUN                                      &
                  )
    endif

    if (condition) then

      if (kkdebug) call seginfo(i,"before GD screw correction ")

      if (seg(i1E)%norme < seg(vnn2e)%norme) then

        seg(vnn2e)%O     = seg(i1e)%O
        seg(vnn2e)%norme = seg(vnn2e)%norme - seg(i1e)%norme
        seg(vnn2e)%dom   = seg(i1e)%dom
        seg(i1e)%norme   = IZERO
        GD2made = .false.
        call connecIJ(I1e,vnn2e,184)

        if (kkdebug) call seginfo(i,"GD screw displacement I1E 1")

      elseif (seg(i1E)%norme == seg(vnn2e)%norme) then

        seg(vnn2e)%O     = seg(i1e)%O
        seg(vnn2e)%norme = IZERO
        seg(vnn2e)%dom   = seg(i1e)%dom
        seg(i1e)%norme   = IZERO
        vnn2e            = seg(vnn2e)%vnne
        GD2made = .false.
        call connecIJ(I1e,vnn2e,185)
        call voisinage(vnn2e,1851)

        if (kkdebug) call seginfo(i,"GD screw displacement I1E 2")

      else

        seg(i1e)%norme   = seg(i1e)%norme - seg(vnn2e)%norme
        seg(vnn2e)%O     = seg(i1e)%O + seg(i1e)%norme*bveclin(:,seg(I1e)%veclin)
        seg(vnn2e)%norme = IZERO
        seg(vnn2e)%dom   = seg(seg(vnn2e)%voise)%dom
        vnn2e            = seg(vnn2e)%vnne
        GD2made = .false.
        call connecIJ(I1e,vnn2e,186)
        call voisinage(vnn2e,1861)

        if (kkdebug) call seginfo(i,"GD screw displacement I1E 3")

      endif

    endif

    !special case i1o = GD, seg(i1o)%vnno = screw segment with length > 0
    !and no rotule between i and i1o: superposition of a two screw segments,
    !annihilation of the superposed screw length
    condition = .false.
    if (seg(i1o)%gd > izero .and. seg(i1o)%surface < IUN ) then
      condition = ( seg(i1o)%norme > IZERO                                  .and. &
                    TYSEG(seg(vnn2o)%veclin)==IUN                           .and. &
                    deassoc(seg(vnn2o)%veclin) /= deassoc(seg(I1o)%veclin)  .and. &
                    (seg(i)%veclin == conec(seg(I1o)%veclin,1) .or.               &
                    seg(i)%veclin == conec(seg(I1o)%veclin,2))              .and. &
                    VDOi(1) == -bveclin(1,seg(vnn2o)%veclin)                .and. &
                    VDOi(2) == -bveclin(2,seg(vnn2o)%veclin)                .and. &
                    VDOi(3) == -bveclin(3,seg(vnn2o)%veclin)                .and. &
                    seg(vnn2o)%surface < IUN                                      &
                  )

    endif
    if (condition) then

      if (kkdebug) call seginfo(i,"before GD screw correction ")

      if (seg(i1o)%norme < seg(vnn2o)%norme) then

        seg(seg(i1o)%vnno)%norme = seg(vnn2o)%norme-seg(i1o)%norme
        seg(i1o)%norme           = IZERO
        seg(i1o)%dom             = seg(i)%dom
        seg(I1o)%O               = seg(i)%O
        GD2made = .false.
        call connecIJ(vnn2o,I1o,187)

        if (kkdebug) call seginfo(i,"GD screw displacement I1O 1")

      elseif (seg(i1o)%norme == seg(vnn2o)%norme) then

        seg(vnn2o)%norme = IZERO
        seg(i1o)%norme   = IZERO
        seg(i1o)%dom     = seg(i)%dom
        seg(I1o)%O       = seg(i)%O
        vnn2o            = seg(vnn2o)%vnno
        GD2made = .false.
        call connecIJ(vnn2o,I1o,188)
        call voisinage(vnn2o,1881)

        if (kkdebug) call seginfo(i,"GD screw displacement I1O 2")

      else

        seg(i1o)%O       = seg(vnn2o)%O
        seg(i1o)%norme   = seg(i1o)%norme-seg(vnn2o)%norme
        seg(vnn2o)%norme = IZERO
        seg(i1o)%dom     = seg(vnn2o)%dom
        vnn2o            = seg(vnn2o)%vnno
        GD2made = .false.
        call connecIJ(vnn2o,i1o,189)
        call voisinage(vnn2o,1891)

        if (kkdebug) call seginfo(i,"GD screw displacement I1O 3")

      endif

  endif

  endif

  if (kkdebug) then
      call seginfo(I," apres dep.   avant reaction                 ")
      if (notannili2e .and. notannili2o) call check_connec(" apres dep.   avant reaction       ")
  endif

  if (domchange > IZERO) cycle SurI

  if(obsta == 0 .or. seg(i)%surface > IZERO) then
    if (seg(i)%surface > IDEUX) stop "surface segment gt 2 is not allowed here!!! "
    cycle
  elseif(obsta > IZERO .and. seg(obsta)%surface > IZERO) then
    cycle
  endif
  !*******************************************************
  ! End of segment displacement When there is no obstacles
  !*******************************************************

  ! Definition of new characteristics of I after displacement
  Oi(1:3)   = SEG(i)%O(1:3)
  LongSegI  = seg(I)%norme
  Ei(1:3)   = Oi(1:3) + longsegI * VLI(1:3)
  i1o       = SEG(i)%voiso
  i1e       = SEG(i)%voise
  vnn1o     = SEG(i)%vnno
  vnn1e     = SEG(i)%vnne
  vnn2o     = SEG(vnn1o)%vnno
  vnn2e     = SEG(vnn1e)%vnne
  ovnn1o    = SEG(obsta)%vnno
  ovnn1e    = SEG(obsta)%vnne
  longsegJ  = seg(obsta)%norme
  Lj        = seg(obsta)%veclin

  !****************************************
  ! Beginning of the treatment of reactions
  !****************************************

  !*********Checking consistency of detection of reaction
  if (anseg .and. Jnct) then
      print *, " Fatal Error: junction and Annihilation ? "
      print *, 'kk:',kk,'i:',i,'obsTA:',obsTA
      STOP
  ENDIF

  !***********************
  !Treatment of junctions
  !***********************

  if (jnct) then
      !segments touching a junction or a gd segment.
      !in case the segment moves naturally reaching a junction (or a gd) we stop the segment at the junction,
      !skipping the treatment of junctions.
      !if the segment character is different from the junction (or a gd) character

      if (vnn1o /= nsegmax ) then

        if ((seg(vnn1o)%jonc .and. tyseg(seg(i)%veclin) /= tyseg(seg(vnn1o)%veclin) .and. seg(obsta)%jonc)  .or.  &
            ((seg(vnn1o)%gd * seg(obsta)%gd) > izero .and. tyseg(seg(i)%veclin) /= iun)                     .or.  &
            seg(i)%jonc ) cycle surI

      elseif (vnn1e /= nsegmax ) then

        if ((seg(vnn1e)%jonc .and. tyseg(seg(i)%veclin) /= tyseg(seg(vnn1e)%veclin) .and. seg(obsta)%jonc)  .or.  &
            ((seg(vnn1e)%gd * seg(obsta)%GD) > izero .and. tyseg(seg(i)%veclin) /= iun)                     .or.  &
            seg(i)%jonc ) cycle surI

      else

        stop "Something went wrong vnn1o and vnn1e == Nsegmax"

      endif

      !Test to track particular case of a collinear junction/annihilation attempt with sharing GD segment
      !Such reaction between two close segments on the same line is useless and therefore not permitted
      if ( seg(i1o)%norme == IZERO .and. seg(vnn1o)%gd > izero .and. vnn2o == obsta .or. &
         & seg(i1e)%norme == IZERO .and. seg(vnn1e)%gd > izero .and. vnn2e == obsta      &
         ) cycle surI

      ! if I and obsta are parallel with covering (i.e. ptobst(1) = -1), the treatment of the future is simple
      ! else, we should know if the intersection point is at the origin/end of I or obsta
      if (ptobst(1) /= -1) then
        ! distance of ptObst from Oi in VLi
        interi = (ptobst(InnI) - Oi(InnI)) / CInnI
        ! In some special cases ptobst is not on the segment I, then a modulo translation
        !is needed to calculate interi
        !print*,"i",ptobst(InnI),Oi(InnI),interi,abs(CInnI),VLi(InnI)
        if (interi < 0) then
          if (VLi(InnI) > 0) then
            interi = (ptobst(InnI) + modur(InnI) - Oi(InnI)) / CInnI
          else
            interi = (ptobst(InnI) - modur(InnI) - Oi(InnI)) / CInnI
          endif
        elseif (interi > longsegI) then
          if (VLi(InnI) > 0) then
            interi = (ptobst(InnI) - modur(InnI) - Oi(InnI)) / CInnI
          else
            interi = (ptobst(InnI) + modur(InnI) - Oi(InnI)) / CInnI
          endif
        endif
        if (interi < 0 .or. interi > longsegI) then
          write(*,*) "the modulo shift of ptobst is not correctly made for I at step",kk,interi,longsegI
          write(*,*) "Segment =",I,'  with direction= ',Li,'  and InnI=',InnI
          write(*,*) "ptobst =",ptobst(InnI),'  Oi = ',Oi(InnI),'  and CInnI=',CInnI
          write(*,*) "interi_ini =",(ptobst(InnI) - Oi(InnI)) / CInnI,"interi =",interi,'  longsegI = ',longsegI
          call seginfo(obsta, "debug junction problem")
          stop
        endif

        ! distance of ptObst from Oj in VLj
        InnJ = Inn(Lj)
        interj = (ptobst(InnJ) - seg(obsta)%o(InnJ)) / bveclin(InnJ,Lj)
        ! In some special cases ptobst is not on the segment obsta, then a modulo translation
        !is needed to calculate interi
        ! print*,"j",kk,obsta,ptobst(InnJ),seg(obsta)%o(InnJ),interj,abs(bveclin(InnJ,Lj)),Bveclin(InnJ,Lj)
        ! call seginfo(obsta, "debug obsta")
        if (interj < 0) then
          if (Bveclin(InnJ,Lj) > 0) then
            interj = (ptobst(InnJ) + modur(InnJ) - seg(obsta)%o(InnJ)) / bveclin(InnJ,Lj)
          else
            interj = (ptobst(InnJ) - modur(InnJ) - seg(obsta)%o(InnJ)) / bveclin(InnJ,Lj)
          endif
        elseif (interj > longsegJ) then
          if (Bveclin(InnJ,Lj) > 0) then
            interj = (ptobst(InnJ) - modur(InnJ) - seg(obsta)%o(InnJ)) / bveclin(InnJ,Lj)
          else
            interj = (ptobst(InnJ) + modur(InnJ) - seg(obsta)%o(InnJ)) / bveclin(InnJ,Lj)
          endif
        endif
        if (interj < 0 .or. interj > longsegJ) then
          write(*,*) "the modulo shift of ptobst is not correctly made for obsta at step",kk,interj,longsegJ
          stop
        endif

        ! cas ou l'intersection est sur les noeuds de i ou de j, on n'interdit plus lajonction
        ! que lorseque le nouveau caractere de la jonction est incompatible avec les voisin
        ! on appliquera plus tard le traitement suivant

        JoncOI = (interi == 0)  ! The obstacle cuts I in Oi
        JoncEI = (interi - longsegI == 0)  ! The obstacle cuts I in Ei
        JoncOJ = (interJ == 0)  ! I cuts the obstacle in OJ
        JoncEJ = (interj - longsegJ == 0) !  I cuts the obstacle in EJ
      else
        if(Lj /= axejoncj(Li,Lj) .or.  Li /= axejonci(Li,Lj)) then
            print *, " Bad diagnostic for junction between parallel segments kk,I",KK,I
            stop
        endif
        ! since ptobst(1) /= -1, we must have
        JoncOI = .false.
        JoncEI = .false.
        JoncOJ = .false.
        JoncEJ = .false.
      endif
      ! corresponds to the junction line between the systems of Li and Ljselec

      if (JoncOI) then
        if (vnn1o /= nsegmax) then
            temp1 = seg(vnn1o)%veclin
            if(nbrot(temp1,LJI,1) > 2) then
              IF (kkdebug) write(379,*) " annulation de jonction en Oi"
              CYCLE surI
            endif
        endif
      endif
      ! intersection sur origine de J
      if (JoncOJ) then
        if (ovnn1o /= nsegmax) then
            temp1 = seg(ovnn1o)%veclin
            if(nbrot(temp1,LJJ,1) > 2) then
              IF (KKDEBUG) write(379,*) " annulation de jonction en Oj"
              CYCLE surI
            endif
        endif
      endif
      ! if kneecap, then we know what exactly to do
      if( .not. kneecap) then
        ! intersection sur origine de I
        if (JoncEi) then
            if (vnn1e /= nsegmax) then
              temp1 = seg(vnn1e)%veclin
              if(nbrot(temp1,LJI,1) > 2) then
                  IF (KKDEBUG) write(379,*) " annulation de jonction en Ei"
                  CYCLE surI
              endif
            endif
        endif
        ! intersection sur origine de I
        if (JoncEj) then
            if (ovnn1e /= nsegmax) then
              temp1 = seg(ovnn1e)%veclin
              if(nbrot(temp1,LJJ,1) > 2) then
                  IF (KKDEBUG) write(379,*) " annulation de jonction en Ej"
                  CYCLE surI
              endif
            endif
        endif
      endif

      if (kkdebug) then
        write(379,*) "  Information before jonction "
        write(379,*) "codeobst =", codeobst
        write(379,*) "anseg  =", anseg, " jnct =", jnct
        write(379,*) "LJI =", LJi, "  LJJ =", LJJ,"  InnI =", InnI
        write(379,*) "ptobst =", ptobst(:)
        write(379,*) "(pt-Oi)Vd =",dot_product(ptobst(:)-OiP(:),Vdi)!
        write(379,*) "interi  =",interi ,",   interj  =", interj
        write(379,*) "JoncOI  =",JoncOI ,",   JoncEI  =", JoncEI
        write(379,*) "JoncOJ  =",JoncOJ ,",   JoncEJ  =", JoncEJ
      endif

      !*******************************************************
      !****a) traitment de I *********************************
      !*******************************************************
      !*i* si Li n'est pas la direction de jonction on l'introduit

      ! cas ou l'intersection est sur les extremites, on teste si les vl des voisins
      ! sont deja bien oriente pour la jonction
      if (JoncOi) then
        ! intersection sur l'origine de I
        PlusSeg = PlusSeg + 1
        I1 = NSEGM + PlusSeg
        jonc_cote_i = I1
        SEG(I1)%JONC = .TRUE.
        SEG(I1)%tJONC = 0
        SEG(I1)%O(1:3) = Oi(1:3)
        SEG(I1)%NORME = izero
        SEG(I1)%VECLIN = Lji
        SEG(I1)%GRAIN = SEG(i)%GRAIN
        seg(i1)%dom=seg(i)%dom

        Itemp = Iboite(I)
        IBoite(I1) = Itemp
        NsegBoite(itemp) =  NsegBoite(itemp) + IUN
        IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I1

        temp1 = seg(I)%voiso
        seg(i)%voiso = i1
        seg(I1)%voiso = temp1
        seg(I1)%vnno = seg(i)%vnno
        seg(I1)%voise = I
        seg(I1)%vnne = I
        if(temp1 /= nsegmax) seg(temp1)%voise = I1

        ! The correct connectivity is now reconstruct
        GD2made = .false.
        call connecij(vnn1o,i1,3541)
        call connecij(i1,I,3542)

      elseif (JoncEi) then
        ! intersection sur l'extremite de I
        PlusSeg = PlusSeg + 1
        I1 = NSEGM + PlusSeg
        jonc_cote_i = i1
        SEG(I1)%JONC = .TRUE.
        SEG(I1)%tJONC = 0
        SEG(I1)%O(1:3) = Oi(1:3) + LongSegI * bveclin(1:3,Li)
        SEG(I1)%NORME = izero
        SEG(I1)%VECLIN = Lji
        SEG(I1)%GRAIN = SEG(i)%GRAIN

        Itemp = Iboite(I)
        IBoite(I1) = Itemp
        NsegBoite(itemp) =  NsegBoite(itemp) + IUN
        IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I1

        temp1 = seg(I)%voise
        seg(i)%voise = i1
        seg(I1)%voiso = I
        seg(I1)%vnno = I
        seg(I1)%voise = temp1
        seg(I1)%vnne = seg(i)%vnne
        if(temp1 /= nsegmax) seg(temp1)%voiso = I1

        if (Nbcvxdom > IUN) then
          call assign_domain(i1,14001)
        else
           SEG(I1)%dom=IUN
        endif

        ! The correct connectivity is now reconstruct
        GD2made = .false.
        call connecij(I,i1,3543)
        call connecij(i1,vnn1e,3544)

      else
        ! in all other cases I will be the future junction
        jonc_cote_i = i
        SEG(i)%JONC = .TRUE.
        SEG(i)%tJONC = 0
        if (Li /= Lji) then
            ! ici on est sur que I deviendra jonc

            !* On decoupe et on introduit la direction de jonction
            !******************************************************
            PlusSeg = PlusSeg+2
            I2 = NSEGM+PlusSeg
            I1 = I2-1

            Loiseg(I1) = LoiSeg(I)
            Loiseg(I2) = LoiSeg(I)

            !* Le premier troncon :
            Oobs(1:3)        = SEG(I)%O(1:3)
            SEG(I1)%O(1:3)   = Oobs(1:3)
            SEG(I1)%NORME    = interi
            SEG(I1)%wait     = izero
            SEG(I1)%grain    = SEG(I)%grain
            seg(I1)%dom      = seg(I)%dom
            SEG(I1)%bloquer  = seg(I)%bloquer
            SEG(I1)%unload   = seg(I)%unload
            SEG(I1)%VECLIN   = Li
            if (seg(I)%surface == IUN) then !in case that the segment is touching the surface with the origin
              seg(I)%surface       = IZERO
              seg(I1)%surface      = IUN
              seg(I1)%Varfreeplan  = seg(i)%Varfreeplan
              seg(I1)%dom          = seg(i)%dom !transfers of domains if surface segment
              seg(I)%Varfreeplan   = IZERO
            endif
            !==================================================================================================
            ! remise du nouveau seg dans la boites correspondantes
            Oobs(:) = modulo( Oobs(:),modur(:))
            Eobs(:) = Oobs(:) + interi/ideux * bveclin(:,Li)
            ! le modulo permet de respecter les CPL pour les boite
            IO = modulo(int(Eobs(1)/tailleBoite(1),DPI),NboitesX) + IUN
            IE = modulo(int(Eobs(2)/tailleBoite(2),DPI),NboitesY) + IUN
            OO = modulo(int(Eobs(3)/tailleBoite(3),DPI),NboitesZ) + IUN
            itemp = B3D_1D(io,ie,oo)
            IBoite(I1) = itemp
            NsegBoite(itemp) =  NsegBoite(itemp) + IUN
            IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I1
            ! cette formule represente une bijonction entre Iboite et les trois indices des boites
            ! c'est utlie dans contact parce qu'il faut deplacer les segemnt par ordre predefini
            !==================================================================================================

            SEG(I1)%GRAIN    = SEG(i)%GRAIN
            IDEP(I1)         = IDEP(I)
            SEG(I1)%RESDEP   = SEG(I)%RESDEP
            seg(i1)%anglevis = seg(i)%anglevis
            !*** PLUS SAGE DE METTRE A ZERO
            TAUTOT(I1)       = TAUTOT(I)
            TAUAPP(I1)       = TAUAPP(I)

            !* Le deuxieme troncon :
            Oobs(:)          = seg(I)%O(1:3)+interi*bveclin(1:3,Li)
            SEG(I2)%O(1:3)   = Oobs(:)
            SEG(I2)%NORME    = seg(i)%norme-interi
            SEG(I2)%wait     = izero
            SEG(I2)%grain    = seg(I)%grain

            SEG(I2)%bloquer  = seg(I)%bloquer
            SEG(I2)%unload   = seg(I)%unload
            SEG(I2)%VECLIN   = Li                         !***
            if (seg(i)%surface == IDEUX) then !in case that the segment is touching the surface with the end
              seg(i)%surface       = IZERO
              seg(i2)%surface      = IDEUX
              seg(i2)%Varfreeplan  = seg(i)%Varfreeplan
              seg(i)%Varfreeplan   = IZERO
            endif
            !==================================================================================================
            ! remise du nouveau seg dans la boites correspondantes
            Oobs(:) = modulo( Oobs(:),modur(:))
            Eobs(:) = Oobs(:) + (seg(i)%norme-interi)/ideux * bveclin(:,Li)
            ! le modulo permet de respecter les CPL pour les boite
            IO = modulo(int(Eobs(1)/tailleBoite(1),DPI),NboitesX) + IUN
            IE = modulo(int(Eobs(2)/tailleBoite(2),DPI),NboitesY) + IUN
            OO = modulo(int(Eobs(3)/tailleBoite(3),DPI),NboitesZ) + IUN
            itemp = B3D_1D(io,ie,oo)
            IBoite(I2) = itemp
            NsegBoite(itemp) =  NsegBoite(itemp) + IUN
            IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I2
            !==================================================================================================

            SEG(I2)%GRAIN    = SEG(i)%GRAIN
            IDEP(I2)         = IDEP(I)
            SEG(I2)%RESDEP   = SEG(I)%RESDEP
            seg(i2)%anglevis = seg(i)%anglevis
            !seg(i2)%dom = seg(i)%dom
            !*** PLUS SAGE DE METTRE A ZERO
            TAUTOT(I2)= TAUTOT(I)
            TAUAPP(I2)= TAUAPP(I)

            !* Les segments rotules et le segment initial au milieu :

            SEG(I)%O(1:3) = SEG(I2)%O(1:3)
            SEG(I)%NORME = 0
            SEG(I)%VECLIN = Lji                         !***

            SEG(I1)%voiso = seg(i)%voiso
            SEG(I1)%vnno = seg(i)%vnno
            SEG(I1)%voise = i
            SEG(I1)%vnne = i
            SEG(I)%voiso = i1
            SEG(I)%vnno = i1
            SEG(I2)%voiso = i
            SEG(I2)%vnno = i
            SEG(I2)%voise = seg(i)%voise
            SEG(I2)%vnne = seg(i)%vnne
            SEG(I)%voise = i2
            SEG(I)%vnne = i2

            if (Nbcvxdom > IUN) then
              call assign_domain(i2,14002)
            else
              seg(i2)%dom=IUN
            endif
            seg(I)%dom=SEG(I2)%dom

            if (seg(I1)%voiso /= nsegmax) seg(SEG(I1)%voiso)%voise=i1
            if (seg(I2)%voise /= nsegmax) seg(SEG(I2)%voise)%voiso=i2

            ! The correct connectivity is now reconstruct
            GD2made = .false.
            call connecij(vnn1o,i1,45)
            call connecij(i1,i,46)

            GD2made = .false.
            call connecij(i,i2,48)
            call connecij(i2,vnn1e,49)

        endif
      ENDIF    ! FIN DE TRAITEMENT DE I

      !*******************************************************
      !****************************   fin de traitment de I **
      !*******************************************************
      !*b* Cas du segment j
      !*********************
      !*i* si Lj n'est pas la direction de jonction on l'introduit

      ! cas ou l'intersection est sur les extremites, on teste si les vl des voisins
      ! sont deja bien oriente pour la jonction
      IF (kneecap) then
        ! if a kneecap was already found in the good direction it is too simple
        jonc_cote_J = obsta
        SEG(obsta)%JONC = .TRUE.
        SEG(obsta)%tJONC = 0

        call voisinage(obsta,4746)

      elseif (JoncOJ) then
        ! intersection sur l'origine de I
        PlusSeg = PlusSeg + 1
        I1 = NSEGM + PlusSeg
        jonc_cote_J = i1
        SEG(I1)%JONC = .TRUE.
        SEG(I1)%tJONC = 0
        SEG(I1)%O(1:3) = seg(obsta)%O(1:3)
        SEG(I1)%NORME = izero
        SEG(I1)%VECLIN = Ljj                         !***
        SEG(I1)%GRAIN = SEG(obsta)%GRAIN
        seg(I1)%dom = SEG(obsta)%dom

        Itemp = Iboite(obsta)
        IBoite(I1) = Itemp
        NsegBoite(itemp) =  NsegBoite(itemp) + IUN
        IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I1

        temp1 = seg(obsta)%voiso
        seg(obsta)%voiso = I1
        seg(I1)%voiso = temp1
        seg(I1)%vnno = seg(obsta)%vnno
        seg(I1)%voise = obsta
        seg(I1)%vnne = obsta
        if(temp1 /= nsegmax) seg(temp1)%voise = I1

        ! The correct connectivity is now reconstruct
        GD2made = .false.
        call connecij(ovnn1o,i1,3741)
        call connecij(i1,obsta,42)

      elseif (JoncEj) then
        ! intersection sur l'extremite de I
        PlusSeg = PlusSeg + 1
        I1 = NSEGM + PlusSeg
        jonc_cote_J = i1
        SEG(I1)%JONC = .TRUE.
        SEG(I1)%tJONC = 0
        SEG(I1)%O(1:3) = seg(obsta)%O(1:3) + LongSegJ * bveclin(1:3,Lj)
        SEG(I1)%NORME = izero
        SEG(I1)%VECLIN = Ljj
        SEG(I1)%GRAIN = SEG(OBSTA)%GRAIN

        Itemp = Iboite(obsta)
        IBoite(I1) = Itemp
        NsegBoite(itemp) =  NsegBoite(itemp) + IUN
        IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I1

        temp1 = seg(obsta)%voise
        seg(obsta)%voise = i1
        seg(I1)%voiso = obsta
        seg(I1)%vnno = obsta
        seg(I1)%voise = temp1
        seg(I1)%vnne = seg(obsta)%vnne
        if(temp1 /= nsegmax) seg(temp1)%voiso = I1

        if (Nbcvxdom > IUN) then
          call assign_domain(i1,14003)
        else
          seg(i1)%dom=IUN
        endif

        ! The correct connectivity is now reconstruct
        GD2made = .false.
        call connecij(OBSTA,i1,354)
        call connecij(i1,Ovnn1e,3514)

      else
        ! in all other cases, Obsta will be the future junction
        jonc_cote_J = Obsta
        SEG(obsta)%JONC = .TRUE.
        SEG(obsta)%tJONC = 0
        if (Lj /= Ljj) then
            ! ici on est sur que J deviendra jonc
            PlusSeg  = PlusSeg + 2
            I2       = NSEGM + PlusSeg
            I1       = I2-1

            Loiseg(I1) = LoiSeg(I)
            Loiseg(I2) = LoiSeg(I)

            !* Le premier troncon :
            Oobs(:)          = SEG(obsta)%O(1:3)
            SEG(I1)%O(1:3)   = Oobs(:)
            SEG(I1)%NORME    = interj
            SEG(I1)%grain    = SEG(j)%grain
            seg(I1)%dom      = seg(obsta)%dom
            SEG(I1)%wait     = izero
            SEG(I1)%bloquer  = SEG(J)%bloquer
            SEG(I1)%unload   = SEG(J)%unload
            SEG(I1)%VECLIN   = Lj                         !***
            if (seg(J)%surface == IUN) then !in case that the segment is touching the surface with the origin
              seg(J)%surface       = IZERO
              seg(i1)%surface      = IUN
              seg(i1)%Varfreeplan  = seg(j)%Varfreeplan
              seg(i1)%dom          = seg(j)%dom !transfers of domains if surface segment
              seg(j)%Varfreeplan   = IZERO
            endif

            !==================================================================================================
            ! remise du nouveau seg dans la boites correspondantes
            Oobs(:) = modulo(Oobs(:),modur(:))
            Eobs(:) = Oobs(:) + (interJ)/ideux * bveclin(:,LJ)
            ! le modulo permet de respecter les CPL pour les boite
            IO = modulo(int(Eobs(1)/tailleBoite(1),DPI),NboitesX) + IUN
            IE = modulo(int(Eobs(2)/tailleBoite(2),DPI),NboitesY) + IUN
            OO = modulo(int(Eobs(3)/tailleBoite(3),DPI),NboitesZ) + IUN
            itemp = B3D_1D(io,ie,oo)
            IBoite(I1) = itemp
            NsegBoite(itemp) =  NsegBoite(itemp) + IUN
            IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I1
            !==================================================================================================
            IDEP(I1)         = IDEP(obsta)
            SEG(I1)%RESDEP   = SEG(obsta)%RESDEP
            seg(i1)%anglevis = seg(obsta)%anglevis
            !*** PLUS SAGE DE METTRE A ZERO
            TAUTOT(I1)       = TAUTOT(obsta)
            TAUAPP(I1)       = TAUAPP(obsta)

            !* Le deuxieme troncon :
            Oobs(:)            = SEG(obsta)%O(1:3)+interj*bveclin(1:3,Lj)
            SEG(I2)%O(1:3)     = Oobs(:)
            SEG(I2)%NORME      = seg(obsta)%norme-interj
            SEG(I2)%grain      = SEG(j)%grain

            SEG(I2)%wait       = izero
            SEG(I2)%bloquer    = SEG(J)%bloquer
            SEG(I2)%unload     = SEG(J)%unload
            SEG(I2)%VECLIN     = Lj                         !***
            if (seg(J)%surface == IDEUX) then !in case that the segment is touching the surface with the end
              seg(j)%surface       = IZERO
              seg(i2)%surface      = IDEUX
              seg(i2)%Varfreeplan  = seg(j)%Varfreeplan
              seg(j)%Varfreeplan   = IZERO
            endif
            !==================================================================================================
            ! remise du nouveau seg dans la boites correspondantes
            Oobs(:) = modulo( Oobs(:),modur(:))
            Eobs(:) = Oobs(:) + (seg(OBSTA)%norme-interJ)/ideux * bveclin(:,LJ)
            ! le modulo permet de respecter les CPL pour les boite
            IO = modulo(int(Eobs(1)/tailleBoite(1),DPI),NboitesX) + IUN
            IE = modulo(int(Eobs(2)/tailleBoite(2),DPI),NboitesY) + IUN
            OO = modulo(int(Eobs(3)/tailleBoite(3),DPI),NboitesZ) + IUN
            itemp = B3D_1D(io,ie,oo)
            IBoite(I2) = itemp
            NsegBoite(itemp) =  NsegBoite(itemp) + IUN
            IndexBoite(Itemp)%ListeSeg(NsegBoite(Itemp)) = I2
            !==================================================================================================
            IDEP(I2)         = IDEP(obsta)
            SEG(I2)%RESDEP   = SEG(obsta)%RESDEP
            seg(i2)%anglevis = seg(obsta)%anglevis
            !*** PLUS SAGE DE METTRE A ZERO
            TAUTOT(I2)       = TAUTOT(obsta)
            TAUAPP(I2)       = TAUAPP(obsta)

            !* Les segments rotules et le segment initial au milieu :

            SEG(obsta)%O(1:3) = SEG(I2)%O(1:3)
            SEG(obsta)%NORME  = 0
            SEG(obsta)%VECLIN = Ljj    !***

            SEG(I1)%voiso = seg(obsta)%voiso
            SEG(I1)%vnno = seg(obsta)%vnno
            SEG(I1)%voise = obsta
            SEG(I1)%vnne = obsta
            SEG(obsta)%voiso = i1
            SEG(obsta)%vnno = i1
            SEG(I2)%voiso = obsta
            SEG(I2)%vnno = obsta
            SEG(I2)%voise = seg(obsta)%voise
            SEG(I2)%vnne = seg(obsta)%vnne
            SEG(obsta)%voise = i2
            SEG(obsta)%vnne = i2

            if (seg(I1)%voiso /= nsegmax) seg(SEG(I1)%voiso)%voise = i1
            if (seg(I2)%voise /= nsegmax) seg(SEG(I2)%voise)%voiso = i2

            if (Nbcvxdom > IUN) then
              call assign_domain(i2,14004)
            else
              seg(i2)%dom=IUN
            endif
            seg(obsta)%dom    = SEG(I2)%dom

            ! The correct connectivity is now reconstruct
            GD2made = .false.
            call connecij(ovnn1o,i1,51)
            GD2made = .false.
            call connecij(i1,obsta,52)

            GD2made = .false.
            call connecij(obsta,i2,54)
            GD2made = .false.
            call connecij(i2,ovnn1e,55)

        endif
        if (.not. seg(jonc_cote_j)%jonc) then
            print* , " CONTACT: mauvaise memorisation de binom cote J, kk = ",kk
            stop
        endif

      endif !*** jonction        ! -------------------/3807/
      if(kkdebug) write(379,*) " binome de jonc: ", jonc_cote_j,jonc_cote_I
      seg(jonc_cote_I)%probadev = zero
      seg(jonc_cote_J)%probadev = zero
      seg(jonc_cote_I)%ijonc = jonc_cote_j
      seg(jonc_cote_J)%ijonc = jonc_cote_I
      IDEP(OBSTA)=IZERO
  endif


  !**************************
  !Treatement of annihilation
  !**************************
  if (anseg) then ! -------------/3810/
      !******************** *** Traitement pour les seg de vect de Burgers id
      !*   Annihilation   * **  et de vect ligne colineaires mais de signes oppose
      !******************** *   se recouvrant apres deplacement
      !*******************************************************
      !*******************************************************
      ! Le principe de l'annihilation est le suivant :
      !  On memorise les VNNO et VNNE des segments s'annihilant,! on determine la longeur de recouvrement et on restaure
      ! la connectivite
      !
      ! I -> VNNE de OBSTA et OBSTA -> VNNE de I
      !*** Les segments ont, lors d'un mvt conventionnel, ete amenes l'un sur l'autre
      !*** ON CONSERVE L ORIGINE DES SEGMENTS I ET OBSTA AINSI QUE LEUR VOISO ET VNNO

      !**************************************************************************************************************

      ! Initialization
      GD2made = .false.

      if (seg(i)%surface > IZERO) then
        call seginfo(i, "A segment touching a free surface is called for annihilation" )
        stop "A segment touching a free surface is called for annihilation - see debug file"
      endif

      if (seg(i)%jonc) then
        itemp = seg(I)%ijonc
        if (kkdebug) write(379,*) " annihilation and I=",i, " already junction with",itemp,": DESTROYED "
        SEG(I)%jonc = .FALSE.
        SEG(I)%ijonc = NSEGMAX
        SEG(I)%tjonc = 0
        SEG(itemp)%jonc = .FALSE.
        SEG(itemp)%ijonc = NSEGMAX
        SEG(itemp)%tjonc = 0
        if(seg(itemp)%norme == 0) call voisinage(itemp,342)
      endif
      if (seg(obsta)%jonc) then
        itemp = seg(obsta)%ijonc
        if (kkdebug) write(379,*) " annihilation and obsta=",obsta, " already junction with",itemp,": DESTROYED "
        SEG(obsta)%jonc = .FALSE.
        SEG(obsta)%ijonc = NSEGMAX
        SEG(obsta)%tjonc = 0
        SEG(itemp)%jonc = .FALSE.
        SEG(itemp)%ijonc = NSEGMAX
        SEG(itemp)%tjonc = 0
        if(seg(itemp)%norme == 0) call voisinage(itemp,343)
      endif
      !*** L'ANNIHILATION EST PRIORITAIRE SUR LES JONCTION ET LE GD
      !***************************************************************************************************************

      DELLONG = (seg(obsta)%o(InnI)-seg(i)%o(InnI))/CInnI

      ! Sometime PBCs induce errors in the DELLONG calculation
      if (abs(DELLONG*CInnI) > modur(InnI)*half) then
        if (seg(i)%o(InnI) < IUN .and. seg(obsta)%o(InnI) > IZERO) then
          DELLONG = (seg(obsta)%o(InnI)- (seg(i)%o(InnI)+ modur(InnI)))/CInnI
        elseif (seg(obsta)%o(InnI) < IUN .and. seg(i)%o(InnI) > IZERO) then
          DELLONG = ((seg(obsta)%o(InnI)+ modur(InnI))-seg(i)%o(InnI))/CInnI
        elseif (seg(i)%o(InnI) > (modur(InnI)+IUN) .and. seg(obsta)%o(InnI) < modur(InnI)) then
          DELLONG = (seg(obsta)%o(InnI)- (seg(i)%o(InnI)- modur(InnI)))/CInnI
        elseif (seg(obsta)%o(InnI) > (modur(InnI)+IUN) .and. seg(i)%o(InnI) < modur(InnI)) then
          DELLONG = ((seg(obsta)%o(InnI)- modur(InnI))-seg(i)%o(InnI))/CInnI
        elseif (seg(i)%o(InnI) >  seg(obsta)%o(InnI) ) then
          DELLONG = ((seg(obsta)%o(InnI)+ modur(InnI))-seg(i)%o(InnI))/CInnI
        elseif (seg(i)%o(InnI) <  seg(obsta)%o(InnI) ) then
          DELLONG = (seg(obsta)%o(InnI)- (seg(i)%o(InnI)+ modur(InnI)))/CInnI
          else
          write(*,*) 'there is something strange with the PBC in the Obstacle part'
          print*,'<---PBC--->',InnI,CInni,DELLONG
          print*,seg(i)%o(InnI)
          print*,seg(obsta)%o(InnI)
          stop
        endif
        if(kkdebug ) then
          write(379,*) '<---PBC---',InnI,CInni,DELLONG
          write(379,*) seg(i)%o(InnI)
          write(379,*) seg(obsta)%o(InnI)
        endif
      endif

      !*** sauvegarde des voisins...
      bufveclin = seg(i)%veclin
      bufvoise  = seg(i)%voise
      bufnorme = seg(i)%norme
      io = seg(i)%voiso
      ie = seg(i)%voise
      oo = seg(obsta)%voiso
      oe = seg(obsta)%voise
      gdHandIo = izero
      gdHandIe = izero
      gdHandOo = izero
      gdHandOe = izero

      !Determination des bras non nul de i et j :  hand pour bras d'annihilation
      HandIo = io
      do while (seg(HandIo)%norme == Izero)
          if (HandIo == seg(HandIo)%voiso) then
            print *, " UPDATE : boucle 62 : i = voiso. !"
            stop
          endif
          HandIoold = HandIo
          HandIo    = seg(HandIo)%voiso
          if (seg(HandIo)%gd > izero) gdHandIo = seg(HandIo)%gd  ! when passing a gd we remember its type
      enddo

      HandIe = ie
      do while (seg(HandIe)%norme == Izero)
          if (HandIe == seg(HandIe)%voise) then
            print *, " UPDATE : boucle 63 : i = voiso. !"
            stop
          endif
          HandIeold = HandIe
          HandIe    = seg(HandIe)%voise
          if (seg(HandIe)%gd > izero) gdHandIe = seg(HandIe)%gd  ! when passing a gd we remember its type
      enddo

      HandOo = Oo
      do while (seg(HandOo)%norme == Izero)
          if (HandOo == seg(HandOo)%voiso) then
            print *, " UPDATE : boucle 64 : i = voiso. !"
            stop
          endif
          if (seg(HandOo)%gd > izero) gdHandOo = seg(HandOo)%gd  ! when passing a gd we remember its type
          HandOoold = HandOo
          HandOo    = seg(HandOo)%voiso
          if (seg(HandOo)%gd > izero) gdHandOo = seg(HandOo)%gd  ! when passing a gd we remember its type
      enddo

      HandOe = Oe
      do while (seg(HandOe)%norme == Izero)
          if (HandOe == seg(HandOe)%voise) then
            print *, " UPDATE : boucle 64 : i = voiso. !"
            stop
          endif
          HandOeold = HandOe
          HandOe    = seg(HandOe)%voise
          if (seg(HandOe)%gd > izero) gdHandOe = seg(HandOe)%gd  ! when passing a gd we remember its type
      enddo

      ! l'annihilation est prioritaire : elle detruite la connectivite des segment
      ! c'est pour quoi il faut trouver les vnno et vnne de norme non nulle
      ! verification
      if (seg(HANDIo)%norme == 0 .or. seg(HANDIe)%norme == 0 .or. &
          seg(HANDOo)%norme == 0 .or. seg(HANDOe)%norme == 0 ) then
        print *," probleme d'annihilation un des hand est nul KK = ",KK
      endif

      !***************************************************************************************************************

      if(kkdebug ) then
        write(379,*) " long anni =", dellong, "  CInnI =", CInnI
        write(379,*) " Oj - Oi =", seg(obsta)%o(:) - seg(i)%o(:)
        write(379,*) " Oji-veclin =", (seg(obsta)%o(:)-seg(i)%o(:)) - dellong*bveclin(:, LI)
        write(379,*) " Oj - Oi * ni=", dot_product(seg(obsta)%o(:)-seg(i)%o(:),bvecnor(:,Li)), &
                      " Oi-Oj * di=", dot_product(seg(i)%o(:)-seg(obsta)%o(:),bvecdep(:,Li))
        write(379,*) " I   =",I , "    OBSTA  =",obsta
        call conf(i)
        call conf(obsta)
        write(379,*) " HIo =", HANDIo, " Io    =", Io, " Ie    =", Ie, " HIe    =", HANDIe
        write(379,*) " HOo =", HANDOo, " Oo    =", Oo, " Oe    =", Oe, " HOe    =", HANDOe
        write(379,*) "HANDIe  =",HANDIe , "HANDOe       =", HANDOe
        write(379,*) "NVNORME I=", dellong-seg(obsta)%norme," NVNORME OBSTA=", dellong-seg(I)%norme
        ! verification
        if (seg(HANDIo)%norme == 0 .or. seg(HANDIe)%norme == 0 .or. &
              seg(HANDOo)%norme == 0 .or. seg(HANDOe)%norme == 0 ) then
            write(379,*) " probleme d'annihilation un des  hand est nul"
            call conf(i)
            call conf(j)
            stop
        endif
      endif

      !*********************************************
      !********calcul des nouvelles normes   *******
      !*********************************************
      NVNORMEI = dellong - seg(obsta)%norme
      NVNORMEO = dellong - seg(I)%norme

      ! on modifie les normes de I et Obsta
      seg(i)%norme = int(abs(NVNORMEI),DPI)
      seg(obsta)%norme = int(abs(NVNORMEO),DPI)
      !*********************************************
      !*********************************************
      ! d'abord le traitement des cas speciaux
      ! cas 0 : boucle collapse
      if (obsta == handIe .and. I == handOe ) then
        if (kkdebug) then
          write(379,*) " ann type 0"
          call seginfo(i, " avant supression de boucle dans contact")
        endif
        itemp = ie
        do i1 = 0, 100
            if(ie == seg(itemp)%voise) exit
            ie = seg(itemp)%voise
            if (ie == itemp .or. itemp == nsegmax) stop " UPDATE : boucle 32 : interdiction d'effacer !"
            seg(itemp)%voiso = itemp
            seg(itemp)%voise = itemp
            out(itemp)       = .true.
            ! if the loop contains gd segments, stat on the gd eliminated is made
            if (seg(itemp)%gd == iun) oldntgd1 = oldntgd1 + 1
            if (seg(itemp)%gd == ideux) then
              oldntgd2 = oldntgd2 + 1
              if (kkdebug) write(379,*) 'In update - 1, oldntgd after icrement =', oldntgd1,oldntgd2
            endif
            seg(itemp)%gd    = izero
            if(seg(itemp)%jonc) then
              i2 = seg(itemp)%Ijonc
              seg(itemp)%jonc = .false.
              seg(itemp)%Ijonc = NSEGMAX
              seg(itemp)%Ijonc = 0
              seg(i2)%jonc = .false.
              seg(i2)%Ijonc = NSEGMAX
              seg(i2)%tjonc = 0
              if(seg(i2)%norme == izero) then
                temp1 = seg(i2)%vnno
                temp2 = seg(i2)%vnne
                GD2made = .false.
                call connecij(temp1,temp2,1306)
              endif
            endif
            itemp = ie
        enddo

      ! cas 1 : obsta = HandIe (as I and obsta overlap on some length, obsta will be eliminated)
      elseif (obsta == handIe) then
        ! de toute facon obsta sera suprime par connecij
        seg(obsta)%norme = izero
        ! cas : 1.1 : i ne change pas de sense : on connecte I a oe
        if (nvnormeI > izero) then
            if (kkdebug) write(379,*) " ann type 1.1, J annihilated, I doesnt change VL"
            GD2made = .false.
            call connecij(i,handOe,2000)
        ! cas : 1.2 : i change de sense : on change veclin de I et
        !  on refait la connection HandIO-I et I-HandOe
        elseif (nvnormeI < izero) then
            if (kkdebug) write(379,*) " ann type 1.2, J annihilated, VLi becomes VL obsta "
            seg(i)%veclin = seg(obsta)%veclin
          !  if (gdHandIe > iun .or. seg(seg(i)%vnne)%gd > iun .or.   &
          !      seg(obsta)%gd > iun .or. seg(seg(obsta)%vnne)%gd > iun) then ! At least 1 GD2 must be made
            if (gdHandIe > iun .or. seg(seg(i)%vnne)%gd > iun) then ! At least 1 GD2 must be made
              GD2made = .false.
              call connecij(HandIo,i,-28011)
              if (GD2made) then
                !GD2made = .false.
                call connecij(i,handOe,2802)
              else
                call connecij(i,handOe,-2802) ! The formation of a GD2 is needed here
              endif
            else
              call connecij(HandIo,i,25011)
              GD2made = .false.
              call connecij(i,handOe,2501)
            endif
        ! cas : 1.3 : i et J s'annule : on connecte handio a handoe
        else
            if (kkdebug) write(379,*) " ann type 1.3, i et J annihilate eatch other "
            GD2made = .false.
            call connecij(handIo,handOe,2002)
        endif

      ! cas 2 : I = handOe (as obsta and I overlap on some length, I will be eliminated)
      elseif (I == handOe) then
        ! de toute facon I sera suprime par connecij
        seg(i)%norme = izero
        ! cas : 2.1 : obsta ne change pas de sense : on connecte I a oe
        if (nvnormeO > izero) then
            if (kkdebug) write(379,*) " ann type 2.1, I annihilated, J doesnt change VL"
            GD2made = .false.
            call connecij(obsta,handIe,2003)
        ! cas : 2.2 : obsta change de sense : on change veclin de J
        !  on refait la connection HandOo-O et O-HandIe
        elseif (nvnormeO < izero) then
            if (kkdebug) write(379,*) " ann type 2.2, I annihilated, VLJ becomes VLI"
            seg(obsta)%veclin = seg(i)%veclin
          !  if (gdHandOe > iun .or. seg(seg(obsta)%vnne)%gd > iun &
          !      .or. seg(i)%gd > iun .or. seg(seg(i)%vnne)%gd > iun) then ! At least 1 GD2 must be made
            if (gdHandOe > iun .or. seg(seg(obsta)%vnne)%gd > iun) then ! At least 1 GD2 must be made
              GD2made = .false.
              call connecij(HandOo,obsta,-28044)
              if (GD2made) then
                ! GD2made = .false.
                call connecij(obsta,handIe,2805)
              else
                call connecij(obsta,handIe,-2805)   ! The formation of a GD2 is needed here
              endif
            else
              call connecij(HandOo,obsta,20044)
              GD2made = .false.
              call connecij(obsta,handIe,2004)
            endif
        ! cas : 2.3 : obsta s'annule : on connecte handOo a handIe
        else
            if (kkdebug) write(379,*) " ann type 2.3, i et J annihilate eatch other"
            GD2made = .false.
            call connecij(handOo,handIe,2005)
        endif

      ! cas general : I et Obsta appartiennent a deux bouts eloignes de dislocation
      ! C-a-d il faut faire le traitement des 4 bras apres l'annihilation
      else
        ! comme on change veclin de I et obsta, on memorise veclin de I
        itemp = seg(i)%veclin

        ! side I
        seg(i)%voise = oe
        if(oe /= nsegmax) seg(oe)%voiso = i
        ! cas 3.1 : nvnorme > 0 : on connect i  a HandOe
        if (nvnormeI > izero) then
            if (kkdebug) write(379,*) " ann en O type 3.1, I is shortened"
            GD2made = .false.
            if (seg(seg(obsta)%vnno)%gd > iun .or.  &
                seg(seg(obsta)%vnne)%gd > iun .or.  &
                seg(seg(i)%vnne)%gd > iun           ) then
              call connecij(i,handOe,-2006)    ! The formation of a GD2 is needed here
            else
              call connecij(i,handOe,2006)
            endif
        ! cas : 3.2 : i change de sense : on change veclin de I
        !  on refait la connection HIO-I et on connecte I a e
        elseif (nvnormeI < izero) then
            if (kkdebug) write(379,*) " ann en O type 3.2, VLi becomes VL obsta"
            seg(i)%veclin = seg(obsta)%veclin
            if (gdHandIe > iun .or.                 &
                seg(seg(i)%vnne)%gd > iun .or.      &
                seg(obsta)%gd > iun .or.            &
                seg(seg(obsta)%vnno)%gd > iun .or.  &
                seg(seg(obsta)%vnne)%gd > iun) then ! At least 1 GD2 must be made
              GD2made = .false.
              call connecij(HandIo,i,-20077)
              if (GD2made) then
                GD2made = .false.
                call connecij(i,handOe,2807)
              else
                GD2made = .false.
                call connecij(i,handOe,-2807) ! The formation of a GD2 is needed here
              endif
            else
              call connecij(HandIo,i,20077)
              GD2made = .false.
              call connecij(i,handOe,2007)
            endif
            if (seg(i)%surface /= IZERO) then
                seg(obsta)%surface=seg(i)%surface
                seg(obsta)%VarFreePlan= seg(i)%VarFreePlan
                if(seg(i)%surface==IUN) seg(obsta)%dom = seg(i)%dom !transfers of domains if surface segment with originn
                seg(i)%surface=IZERO
                seg(i)%VarFreePlan=IZERO
            endif
        ! cas : 3.3 : i s'annule : on connecte handio a handoe
        else
            if (kkdebug) write(379,*) " ann en O type 3.3, HandIo is connected to HandOe"
            if (handIo /= NSEGMAX .or. handOe /= NSEGMAX) then
              GD2made = .false.
              call connecij(handIo,handOe,2008)
            else
            !In case that two segment touching the surface (belonging to the same loop)
            !bump into each other we eliminate the two kneecaps left
              out(HandIoold)=.true.
              seg(HandIoold)%voiso = HandIoold
              seg(HandIoold)%voise = HandIoold
              out(HandIoold)       = .true.
              out(HandOeold)=.true.
              seg(HandOeold)%voiso=HandOeold
              seg(HandOeold)%voise=HandOeold
              out(HandOeold)       = .true.
              seg(i)%voiso=i
              seg(i)%voise=i
              out(i)=.true.
              seg(obsta)%surface=IZERO
              seg(i)%surface=IZERO
              seg(i)%VarFreePlan=IZERO
              seg(obsta)%VarFreePlan=IZERO
            endif
        endif

        ! side obsta
        seg(Obsta)%voise = Ie
        if(Ie /= nsegmax) seg(Ie)%voiso = obsta
        ! cas : 4.1 : obsta ne change pas de sense : on connecte I a oe
        if (nvnormeO > izero) then
            if (kkdebug) write(379,*) " ann en E type 4.1, Obsta is shortened"
            GD2made = .false.
            if (seg(seg(i)%vnne)%gd > iun .or.  &
                seg(seg(obsta)%vnne)%gd > iun   ) then
              call connecij(obsta,handIe,-2009)    ! The formation of a GD2 is needed here
            else
              call connecij(obsta,handIe,2009)
            endif
        ! cas : 4.2 : obsta change de sense : on change veclin de I
        !  on refait la connection HOo-O et on connecte O a Hie
        elseif (nvnormeO < izero) then
            if (kkdebug) write(379,*) " ann en E type 4.2, VLJ becomes VLI"
            seg(obsta)%veclin = itemp
            if (gdHandOe > iun .or. seg(seg(obsta)%vnne)%gd > iun &
                .or. seg(i)%gd > iun .or. seg(seg(i)%vnne)%gd > iun) then  ! At least 1 GD2 must be made
              GD2made = .false.
              call connecij(HandOo,obsta,-20101)
              if (GD2made) then
                GD2made = .false.
                call connecij(obsta,handIe,2710)
              else
                call connecij(obsta,handIe,-2710) ! The formation of a GD2 is needed here
              endif
            else
              GD2made = .false.
              call connecij(HandOo,obsta,20101)
              GD2made = .false.
              call connecij(obsta,handIe,2010)
            endif
            if (nvnormeI >= izero) then
              if (seg(obsta)%surface /= IZERO) then
                seg(i)%surface=seg(obsta)%surface
                seg(i)%VarFreePlan= seg(obsta)%VarFreePlan
                if(seg(i)%surface==IUN) seg(i)%dom = seg(obsta)%dom !transfers of domains if surface segment with origin
                seg(obsta)%surface=IZERO
                seg(obsta)%VarFreePlan=IZERO
              elseif (seg(i)%surface /= IZERO) then
                seg(obsta)%surface=seg(i)%surface
                seg(obsta)%VarFreePlan= seg(i)%VarFreePlan
                if(seg(i)%surface==IUN) seg(obsta)%dom = seg(i)%dom !transfers of domains if surface segment with origin
                seg(i)%surface=IZERO
                seg(i)%VarFreePlan=IZERO
              endif
            endif
        ! cas : 4.3 : obsta s'annule : on connecte handOo a handIe
        else
            if (kkdebug) write(379,*) " ann en E type 4.3, HandOo is connected to HandIe"
            if (handOo /= NSEGMAX .or. handIe /= NSEGMAX) then
              GD2made = .false.
              call connecij(handOo,handIe,2002)
            else
              !In case that two segment touching the surface (belonging to the same loop)
              !bump into each other we eliminate the two kneecaps left
              seg(HandOoold)%voiso=HandOoold
              seg(HandOoold)%voise=HandOoold
              out(HandOoold)=.true.
              seg(HandIeold)%voiso=HandIeold
              seg(HandIeold)%voise=HandIeold
              out(HandIeold)=.true.
              seg(obsta)%voiso=obsta
              seg(obsta)%voise=obsta
              out(obsta)=.true.
              seg(i)%surface=IZERO
              seg(obsta)%surface=IZERO
              seg(i)%VarFreePlan=IZERO
              seg(obsta)%VarFreePlan=IZERO
            endif
        endif
      endif
  endif

  !***************************************************************
  ! End of the reactions treatment, i.e. junction and annihilation
  !***************************************************************

  if (kkdebug) then
      call seginfo(I," apres REACTION ---------------------------     ")
      if (notannili2e .and. notannili2o) call check_connec(" apres            reaction       ")
  endif

enddo SURI

!******************************************************************
! End of the loop used to apply the displacement on every segment I
!******************************************************************

! Update of important quantities
dAIRE(1+NSEGM:NSEGM+PLUSSEG) = ZERO          ! aire balayee par le segment

Nsegm = Nsegm + PlusSeg
if(nsegm > nsegmax) then
  write(*,*) 'Problem in (8) as nsegm > nsegmax'
  write(*,*) nsegm, PlusSeg
  stop
endif
PlusSeg = 0

!****************************************************************************
! The segments are finally all replaced in the reference simulation cell (PBC)
!****************************************************************************
if(.not. decalageCPL)then
   seg(1:nsegm)%o(1)=modulo(seg(1:nsegm)%o(1),modur(1))
   seg(1:nsegm)%o(2)=modulo(seg(1:nsegm)%o(2),modur(2))
   seg(1:nsegm)%o(3)=modulo(seg(1:nsegm)%o(3),modur(3))
else
   write(*,*) "Boudary shift have not been used for a long time, they must be check before using it"
   stop
   do jtemp=1,nsegm
       NBDIM(1:3)=nint((seg(jtemp)%o(1:3)-modulo(int(seg(jtemp)%o(1:3),DPI),modur(1:3)))*InvModur(1:3),DPI)
       DECAL(1)=SUM(SHIFT(1,1:3)*NBDIM)
       DECAL(2)=SUM(SHIFT(2,1:3)*NBDIM)
       DECAL(3)=SUM(SHIFT(3,1:3)*NBDIM)
       seg(jtemp)%o(1:3)=modulo(seg(jtemp)%o(1:3)+DECAL(1:3),modur(1:3))
   enddo
endif

#ifdef OUTPUTSURF
  close(999)
#endif

#ifdef MDC

#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
#endif

    if (zcall) then
       sweptsurfdata(1)=nbswept

       if(nbswept > nsegm_INI*mdc_timestep_ratio) then
         write(*,*) 'problem with sweptsurfadata array'
         write(*,*) 'swept area number > size of the array' &
                ,nbswept,'>',nsegm_INI*mdc_timestep_ratio
         write(*,*) 'simulation step , nsegm', kk , nsegm
       endif

       call zebmpisend(sweptsurfdata, 1+16*nbswept)

       nbswept=0
       sweptsurfdata(:)=0
    endif

#if defined(PA) && defined(MDC)
!    ! Le proc zero a fini d ecrire, on libere tout le monde
!    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
#endif

#endif

!****************************************************
! Update of important quantities for the creep mode *
!****************************************************

!> \ingroup creep
!> \note The management of the two lists par_touch() and par_cut() in creep mode \{
 if (creep) then
  !--------------------------------------------------------------------
  !> - A new particle cutting event induced by climb must be proceeded when Timebreak=true \n
  !! 1) We first calculate the climbing time (time_climb) to reach the h_par_touch distance for all segments touching a particle \n
  !! 2) Form the list of touching segments, the segment with the lowest time_climb is selected to define a new par_cut dislocation section \n
  !! 3) The new segments h_par_touch value at the end of the timebreak is defined \n
  !! 4) The selected segment (touch2cut) is indexed in the par_cut(:,:) list \n
  !! 5) As the segment in plane touch2cut moved from par_touch to par_cut state, the corresponding index in the par_touch list must be eliminate \n
  !! 6) The timebreak info is saved in the file 'file_timebreak' to reconstruct the creep curve at the end of the calculation \n
  !--------------------------------------------------------------------
  if (timebreak) then

    allocate(time_climb(1:npar_touch))
    allocate(vitesse_climb(1:npar_touch))

    ! We first calculate the climbing time for the segments touching a particle to reach the h_par_touch distance
    do ipar_touch = 1, npar_touch

!       print*,"angle_par_touch(ipar_touch) =", angle_par_touch(ipar_touch),halfpii, halfpii - angle_par_touch(ipar_touch)
!       print*,"tau_par_touch(ipar_touch) =", tau_par_touch(ipar_touch)

      ! the climbing segment character
      theta = halfpii - angle_par_touch(ipar_touch)

      ! the climbing segment resolved stress in the climb direction
      tauclimb = tau_par_touch(ipar_touch)

      ! time to reach h_par_touch
      vitesse_climb(ipar_touch) = vclimb(tauclimb,theta)
      time_climb(ipar_touch) = h_par_touch(ipar_touch) / vitesse_climb(ipar_touch)

    enddo

    ! Form the list of touching line section, the dislocation with the lowest time_climb
    ! is selected to define a new par_cut dislocation section
    touch2cut = minloc(time_climb(1:npar_touch), dim = 1)

    ! The new segment h_par_touch values at the end of the timebreak are defined
    h_par_touch(1:npar_touch) = h_par_touch(1:npar_touch) - (vitesse_climb(1:npar_touch) * time_climb(touch2cut))

    ! the selected segment (touch2cut) is indexed in the par_cut(:,:) list
    npar_cut = npar_cut + iun
    par_cut(npar_cut,1) = par_touch(touch2cut,1)
    par_cut(npar_cut,2) = par_touch(touch2cut,2)
    par_cut(npar_cut,3) = par_touch(touch2cut,3)
    par_cut(npar_cut,4) = iun
    par_cut(npar_cut,5) = izero
    timebreak = .false.

    ! As the segment in plane touch2cut moved from par_touch to par_cut state,
    ! it must be removed from the par_touch list
    !     print*,'in Timebreak we selected par_touch(touch2cut,1)',par_touch(touch2cut,1)
    !     print*,"h=",h_par_touch(1:npar_touch)
    !     print*,"v=",vitesse_climb(1:npar_touch)
    !     print*,"t",time_climb(1:npar_touch)
    par_touch(touch2cut,4) = 0
    par_touch(touch2cut,5) = -101

#ifdef PA
    if (Mon_Rang_PairImpair ==  0) then
      ! Only the process zero of the communicator can write in output files, the others are waiting
#endif
      ! The timebreak info is saved in the file 'file_timebreak' to reconstruct the creep curve at the end of the calculation
      open(79,FILE=file_timebreak,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
      write(79,'(1X,I9,3X,2(E11.4,1X))') kk, time_climb(touch2cut), vitesse_climb(touch2cut)*avalue
      close(79)
#ifdef PA
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    else
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    endif
#endif

    deallocate(time_climb)
    deallocate(vitesse_climb)

  endif

  !-----------------------------------------------------------------------------
  !> - A  cleanup loop procedure to eliminate the empty touching planes in the par_touch list \n
  !! An element of the par_touch list is eliminate, if a segment blocked in this part_touch plane was moved in the par_cut list (see above) \n
  !! or when no segments was found in the corresponding part_touch plane during 100 steps (the dislocation line was able to shear the particle without climb) \n
  !-----------------------------------------------------------------------------
  icont2 = -1
  do ipar_touch = 1, npar_touch
    if ((ipar_touch + (icont2 + 1)) > npar_touch) exit      ! test to exit the cleanup procedure

    if (par_touch(ipar_touch,4) < 1) then                   ! test on the number of steps with no segment found in the touch plane
      par_touch(ipar_touch,5) = par_touch(ipar_touch,5) -1
    endif

    if (par_touch(ipar_touch,5) < -100) then                ! After 100 steps free of segments in the touching plane, the latter is removed

      ! This touching plane can be eliminated as all dislocations are now outside of the particle
      par_touch(ipar_touch,:) = izero
      h_par_touch(ipar_touch) = -un
      angle_par_touch(ipar_touch) = -100.
      tau_par_touch(ipar_touch) = zero
      icont2 = icont2 + 1                               ! number of touching plane eliminated

      if ((ipar_touch + (icont2 + 1)) > npar_touch) exit    ! test to exit the cleanup procedure

      iswitch = npar_touch - icont2
      par_touch(ipar_touch,:) = par_touch(iswitch,:)          ! we switch position in the par_touch tab
      h_par_touch(ipar_touch) = h_par_touch(iswitch)
      angle_par_touch(ipar_touch) = angle_par_touch(iswitch)
      tau_par_touch(ipar_touch) = tau_par_touch(iswitch)
    endif
  enddo

  ! the new number of touching planes
  npar_touch = npar_touch - (icont2 + 1)

  !--------------------------------------------------------------------------
  !> - A cleanup loop procedure to eliminate the empty cutting planes in the par_cut list \n
  !! If a par_cut plane was found empty of segments during 100 steps, then the corresponding index \n
  !! of the list can be removed as no segments are cutting the particle in this particular plane.
  !--------------------------------------------------------------------------
  icont2 = -1
  do ipar_cut = 1, npar_cut
    if ((ipar_cut + (icont2 + 1)) > npar_cut) exit      ! test to exit the cleanup procedure

    if (par_cut(ipar_cut,4) < 1) then                   ! test on the number of steps with no segment found in the cutting plane
      par_cut(ipar_cut,5) = par_cut(ipar_cut,5) -1
    endif

    if (par_cut(ipar_cut,5) < -100) then                ! After 100 steps free of segments in the cutting plane, the latter is removed

      ! This cutting plane can be eliminated as all dislocations are now outside of the particle
      par_cut(ipar_cut,:) = izero
      icont2 = icont2 + 1                               ! number of cutting plane eliminated

      if ((ipar_cut + (icont2 + 1)) > npar_cut) exit    ! test to exit the cleanup procedure

      iswitch = npar_cut - icont2
      par_cut(ipar_cut,:) = par_cut(iswitch,:)          ! we switch position in the par_cut tab
    endif
  enddo

  ! the new number of cutting planes
  npar_cut = npar_cut - (icont2 + 1)

endif
!> \}

deallocate (OBSTAB)
deallocate (MiniListeSeg)
#ifdef PA
if (TAILLE_PAIRIMPAIR /= 1) then
  deallocate(MiniNseg_Send)
  deallocate(MiniNseg_Recv)
endif
deallocate (MiniListe_Send)
deallocate (MiniListe_Recv)
#endif

if (NbDep /= 0) deallocate(DEP_REAC)

end subroutine UPDATE


!######################################################################
!# A function defining the probabity that an immobilized segment      #
!# glide in a Cr precipitate according to an Arrhenius                #
!# low acconting for thermally activated process.                     #
!######################################################################
function JUMP_Cr(taueff_Pa,longsegI_a,Ipar,radius_a)

implicit none

real(kind=DP), intent(in)       :: taueff_Pa, Radius_a,longsegI_a
real(kind=DP)                   :: deltaG,proba_th,tir, taueff,diameter,tau_max,longsegI
integer(kind=DPI), intent(in)   :: Ipar
Logical                         :: Jump_Cr, all_sizes

DeltaG = dix
taueff = taueff_Pa * 1.D-9  ! put tau_eff in GPa unit
LongsegI = longsegI_a * avalue * 1.D+9  ! put seg length in nm
!jump_rint *, taueff_Pa,taueff ,par(Ipar)%tau

if(Radius_a > Bdiva) then
   Diameter = Radius_a * DEUX * avalue * 1.D+9 ! diameter in nm
else
   Diameter = BDIVA * DEUX * avalue * 1.D+9 ! diameter in nm
endif

Jump_Cr = .false.
! There are philisophies for the computation of the activation energy DeltaG
! DeltaG must scale either with the segment length or with is particle size

! all_lengths = T means DeltaG is independent of the segment length but proportional to Diameter
! all_lengths = F means DeltaG is proportional to the segment length but independent of the Diameter
all_sizes = .false.
! to reproduce infinite activation volume region
tau_max =  par(Ipar)%tau  ! tau_max in GPa

! If par(Ipar)%tau < zero, this is the CONVENTION of Orowan particle.
if (tau_max > zero .and. taueff > zero) then
    if (taueff >= tau_max) then
      Jump_Cr = .true.
    else
        If(temperature > zero .and. taueff > 0.8) then
            ! all_sizes : means average of all MD data inluding the 1nm precipitate data
            if (all_sizes) then
                ! delta_G per nm
                deltaG = 1.85 - 2.26*taueff + 0.94 *taueff*taueff - 0.134*taueff*taueff*taueff
                ! delta_G  = delta_G per nm * length
                deltaG = deltaG * Diameter
                ! the frequency of jump is proportional to Debye and to (b/L)
                proba_th = Debye * half * Bdiva / radius_a * deltat * exp(-deltaG * UnkT)
                if(PROBA_th > 1.0D-5) then
                    TIR = taus88()
                    if (proba_th > tir) then
                      JUMP_Cr = .true.
                      !print *, "ACTIVATION :",deltaG,Bdiva/Lseg_a,taueff,PROBA_th, temperature
                      !read(*,*)
                      ! print *, " apriori thermal permission for seg ",I, " and par :",ktemp
                    endif
                endif
                ! delta_G0. per nm
            else
                ! The folowing deltaG function are not appropriate for the 1nm precipitate
                ! computation of Delta G per nm
                if(taueff <= 1.15) then
                    deltaG = 2.48 - 1.85 * taueff
                    !print *, " basse, taueff = ", taueff
                elseif(taueff > 1.15 .and. taueff <= 1.74) then
                    deltaG = 0.985 - 0.55 * taueff
                    !print *, " intermediate, taueff = ", taueff
                elseif(taueff > 1.74) then
                    deltaG = 0.13 - 0.06 * taueff
                    !print *, " haute, taueff = ", taueff
                endif

                if(deltaG < zero) DeltaG = zero
                ! computation of Delta G for the effective Cr precipitate siaze
                deltaG = deltaG * LongsegI
                !print *, " diameter (nm)", Diameter
                ! the frequency of jump is proportional to Debye and to (b/L)
                if(deltaG < 0.01) then
                    JUMP_Cr = .true.
                else
                    proba_th = Debye * half * Bdiva / radius_a * deltat * exp(-deltaG * UnkT)
                    if (PROBA_th > 1.0D-5) then
                        TIR = taus88()
                        if (proba_th > tir) then
                            JUMP_Cr = .true.
                            print *, " tau",taueff,"DG ",deltaG, "proba ", PROBA_th
                            !print *, "ACTIVATION :",deltaG,Bdiva/Lseg_a,taueff,PROBA_th, temperature
                            !read(*,*)
                            ! print *, " apriori thermal permission for seg ",I, " and par :",ktemp
                        endif
                    endif
                endif
            endif
        endif
    endif
endif
!if(jump_void) print *, "taueffred, proba", taueff_red,taueff/1000000000,proba_th

end function JUMP_Cr

!################################################################
!> \ingroup creep
!> This function provides the climb velocity at prescribed
!! temperature as a function of the climb force, dislocation
!! character and assumed remotely constant concentration.
!! \cite Bako:2011fk
!################################################################
function vclimb(tauclimb,theta)

implicit none

real(kind = DP)               :: vclimb     !< The climb velocity
real(kind = DP), intent(in)   :: tauclimb   !< climb component of PK force * b (computed in ELASTI)
real(kind = DP), intent(in)   :: theta      !< angle to the screw direction in radian
real(kind = DP)               :: prefac, Ceqnorm, arg

  ! Initialization
  if (tauclimb < 1.D-19) then
    vclimb = 1.D-20 ! Tauclimb is not defined (e.g. screw segments), hence the segment climb velocity is set to a tiny value
    return
  else
    vclimb = zero
  endif

  ! The default solution is edge dislocation segment
  prefac = un
  if (theta > 0.1) then     ! theta must always be > 0
    prefac = (un / dsin(theta))
  endif

  Ceqnorm = 1.0D13                      ! The saturation value we consider
  arg = UnkT_Joule * tauclimb * xmu * Omega * prefac
  if (arg < 30.D0) Ceqnorm = dexp(arg)  ! The real value

  ! The climb velocity in avalue/s
  ! As the direction of climb is not taken into account for the moment
  ! better to apply an abs function here !!!
  vclimb = abs(prefac2 * prefac * (Ceqnorm - satfac) / avalue)

  if (vclimb == zero) then
    write(*,*)  "problem in the vclimb function; vclimb = 0 !"
    print*,"tauclimb = ", tauclimb
    print*,"theta = ", theta
    print*,"prefac =", prefac
    print*,"prefac2 =", prefac2
    print*,"arg =", arg, UnkT, tauclimb, Omega
    print*,"Ceqnorm =", Ceqnorm
    print*,"satfac =", satfac
    print*,"xmu =",xmu
    stop
  endif

end function vclimb


!################################################################################################
!                                         DETECT_LOOP
!
! This subroutine is a loop on all possible obstacles segments
! the most important obstacles are listed in arrays OBSTAB
! the output of the DETECT is:
!  OBSTAB(nbobs)%numseg = indice of the obstacle segment
!  OBSTAB(nbobs)%ptobs  =  default : (0 0 0)
!          intersection point : planar and non planar obstacles with one point of intersection
!          case of parallel segments with are superposed: convention  (-1 -1 -1)
!  OBSTAB(nbobs)%plani = a key for coplanar obstacles
!  OBSTAB(nbobs)%reac = will be defined in the following except (sub-lattice(4) and junction(5)
!  DistIobs : Distance between I and obstacles  All obstacles are at the same distance
!#################################################################################################
subroutine Detect_loop(Nsegtester,nbobs,FirstTime)

implicit none

integer(kind=DPI) :: NsegTester  ! The number of obstacle to test
integer(kind=DPI) :: nbobs       ! The reduce number of obstacle found in the loop

integer(kind=DPI) ,dimension(3) :: Oj, Ej, Cj, Aj, Dj, Vlj, VLi1o
integer(kind=DPI) ,dimension(3) :: OiOj, OiEj, EiOj, EiEj, OjCj, OiCj, EiCj

integer(kind=DPI) :: DebList, EndList
integer(kind=DPI) :: j, itemp, ktemp
integer(kind=DPI) :: Lj, normeJ
integer(kind=DPI) :: OiCjscaVDi1o, EiCjscaVDi1e, OiAjscaVDi1o
integer(kind=DPI) :: EiAjscaVDi1e, OiDjscaVDi1o, EiDjscaVDi1e
integer(kind=DPI) :: OiOjscaN, OiEjscaN, OiOjscaVDi1o, OiEjscaVDi1o, EiOjscaVDi1e, EiEjscaVDi1e, OipCjscaVLi
integer(kind=DPI) :: normsqrOjCj, normsqrVLj, normsqrVLi
integer(kind=DPI) :: proj

real(kind=DP)     :: distij
real(kind=DP)     :: OiOjprojVDi, OiEjprojVDi, OiCjprojVDi, OiAjprojVDi, OiDjprojVDi
real(kind=DP)     :: temp

logical           :: pbssreso, jjonc, surrail, jparali, jplani
logical           :: FirstTime    ! Logical parameter noting first and second call of DETECT_LOOP
logical           :: TakeObst     ! Logical parameter noting real segment obstacle in the loop detection

#ifdef PA
logical           :: ParaComNeeded
#endif

!----------------------------------
! Initialization
nbobs           = 0
distij          = 99999999
MiniNsegTester  = 0
MiniListeSeg(1:dim_alloc) = 0

do ktemp = 1 ,dim_alloc
  OBSTAB(ktemp)%numseg = izero
  OBSTAB(ktemp)%react = izero
  OBSTAB(ktemp)%ptobs(:) = (/0,0,0/)
  OBSTAB(ktemp)%plani = .false.
enddo

#ifdef PA

! In parallel calculation the obstacles detection can be made with one or two pass in the subroutine
if (FirstTime .and. ((TAILLE_PAIRIMPAIR * 2) < NsegTester) .and. (TAILLE_PAIRIMPAIR /= 1)) then

  ! Two pass are useful in parallel calculation if NsegTester is large

  ! Comm buffers initialization
  MiniNseg_Send(1:TAILLE_PAIRIMPAIR) = 0
  MiniNseg_Recv(1:TAILLE_PAIRIMPAIR) = 0
  MiniListe_Send(1:dim_alloc)        = 0
  MiniListe_Recv(1:dim_alloc)        = 0

  ! The calculation dispatch between procs
  DebList = ((NsegTester/TAILLE_PAIRIMPAIR) * MON_RANG_PAIRIMPAIR) + 1
  EndList = DebList + (NsegTester/TAILLE_PAIRIMPAIR) - 1

  if ((MON_RANG_PAIRIMPAIR + 1) == TAILLE_PAIRIMPAIR)  EndList = NsegTester   ! we want to end correctly the loop

  ParaComNeeded = .true.

else

  ! The standard solution with one pass
  DebList = 1
  EndList = NsegTester

  ParaComNeeded   = .false.

endif

#else

  ! The standard solution with one pass
  DebList = 1
  EndList = NsegTester
  Firsttime = .true.

#endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! The detect loop of the subroutine
DETECT : do ktemp = DebList, EndList

  ! The tested segment (obstacle) index
  j = listeSegTester(ktemp)

  !================================================================================================
  ! les tests suivant visent a eliminer le plus vite possible les segments qui ne sont pas obstacles
  ! de j. Leur ordre et leur formulation ont ete optimises
  !================================================================================================
  ! the tested segments must not be a pivot segment or a segment to eliminate
  if (out(j)) cycle DETECT

  normeJ=seg(j)%norme
  if (normej == IZERO) cycle DETECT
  ! Segments in a different grain are not obstacle
  if ((GB == ITROIS) .and. (seg(i)%grain /= seg(j)%grain)) cycle DETECT
  ! Segments touching a free surface are for reason of simplicity not a possible obstacle
  if ((GB /= izero) .and. (seg(j)%surface /=  izero)) cycle DETECT

  ! A test to detect problems
  if (J > nsegm + plusseg .or. J == izero) then
    print *, "I =",I,"  indice boucle ",ktemp," segment SURJ",J
    stop
  endif

  ! Initialization
  EiCjscaVDi1e  = izero
  OiCjscaVDi1o  = izero
  OiCjprojVDi   = izero
  OjCj(1:3)     = 99999999
  OiCj(1:3)     = 99999999
  surrail       = .false.
  jparali       = .false.
  jplani        = .false.
  jjonc         = .false.
  TakeObst      = .true.
  Cj(:)         = (/0,0,0/)

  !================================================================================================
  ! Can segment J cut segment I swept area ?
  !================================================================================================
  ! If J cut I slip plane, Oj and Ej are over and under I slip plane
  Oj(1:3) = SEG(j)%O(1:3)
  Lj = SEG(j)%VECLIN
  Ej(1:3) = Oj(1:3)+normeJ*BVECLIN(1:3,Lj)
  OiOj(1:3) = Oj(1:3)-Oi(1:3)
  OiEj(1:3) = Ej(1:3)-Oi(1:3)
  OiOjscaN = OiOj(1) * Nvec(1)+ OiOj(2) * Nvec(2) + OiOj(3) * Nvec(3)
  OiEjscaN = OiEj(1) * Nvec(1)+ OiEj(2) * Nvec(2) + OiEj(3) * Nvec(3)

  temp =  real(OiOjscaN,DP) * real(OiEjscaN,DP)
  if (temp > zero) cycle DETECT

  ! Is J in front of I ?
  OiOjprojVDi = (OiOj(1)*VDi(1)+OiOj(2)*VDi(2)+OiOj(3)*VDi(3))/NormsqrVDi   !> Projection of OiOj on Vdi (VDi units)
  OiEjprojVDi = (OiEj(1)*VDi(1)+OiEj(2)*VDi(2)+OiEj(3)*VDi(3))/NormsqrVDi   !> Projection of OiEj on Vdi (VDi units)
  if (OiOjprojVDi <= ZERO .and. OiEjprojVDi <= ZERO) cycle DETECT           !> If both segment j extremities are in front of i
                                                                            !!(VDi direction) then j is not an obstacle - Cycle

  ! Is J behind previous obstacle ?
  if (OiOjprojVDi > distiobs  .and.  OiEjprojVDi > distiobs) cycle DETECT   !> distiobs to keep the closest J
  ! Segments whose extremities are behind distiobs are not obstacles

  ! Are J extremities in front of the rail in Oi ?
  OiOjscaVDi1o = OiOj(1)*VDi1o_or(1)+OiOj(2)*VDi1o_or(2)+OiOj(3)*VDi1o_or(3)
  OiEjscaVDi1o = OiEj(1)*VDi1o_or(1)+OiEj(2)*VDi1o_or(2)+OiEj(3)*VDi1o_or(3)
  if (OiOjscaVDi1o < IZERO  .and.  OiEjscaVDi1o < IZERO) cycle DETECT
  ! Segments whose extremities are behind rails are not obstacles

  ! IF J is VNN2 then J is not an obstacle if the dislocation is withdrawing
  if (j==vnn2o.and.((OiOjscaVDi1o<=IZERO.and.OiEjscaVDi1o<=IZERO).or.OiEjprojVDi>distiobs)) cycle DETECT

  ! Are J extremities in front of the rail in Ei ?
  EiOj(1:3) = Oj(1:3)-Ei(1:3)
  EiEj(1:3) = Ej(1:3)-Ei(1:3)
  EiOjscaVDi1e = EiOj(1)*VDi1e_or(1) + EiOj(2)*VDi1e_or(2) +EiOj(3)*VDi1e_or(3)
  EiEjscaVDi1e = EiEj(1)*VDi1e_or(1) + EiEj(2)*VDi1e_or(2) +EiEj(3)*VDi1e_or(3)
  if( EiOjscaVDi1e < IZERO .and. EiEjscaVDi1e < IZERO) cycle DETECT
  ! Segments whose extremities are behind rails are not obstacles

  ! IF J is VNN2 then J is not an obstacle if the dislocation is withdrawing
  if (j==vnn2e.and.((EiOjscaVDi1e<=IZERO.and.EiEjscaVDi1e<=IZERO).or.OiOjprojVDi>distiobs)) cycle DETECT

  ! Remaining segments are potentially obstacles closest than the previous obstacle
  ! The final intersection point is still needed to be defined
  if (kkdebug) then
    write(379,fmt='("Seg i : ",I10)') i
    write(379,fmt='("Seg j : ",I10)') j
    write(379,fmt='("  OiOjscaN : ",I10)') OiOjscaN
    write(379,fmt='("  OiEjscaN : ",I10)') OiEjscaN
    write(379,fmt='("  OiOjprojVDi : ",F10.1)') OiOjprojVDi
    write(379,fmt='("  OiEjprojVDi : ",F10.1)') OiEjprojVDi
    write(379,fmt='("  OiOjscaVDi1o : ",I10)')  OiOjscaVDi1o
    write(379,fmt='("  OiEjscaVDi1o : ",I10)')  OiEjscaVDi1o
    write(379,fmt='("  EiOjscaVDi1e : ",I10)')  EiOjscaVDi1e
    write(379,fmt='("  EiEjscaVDi1e : ",I10)')  EiEjscaVDi1e
  endif

  !=================================================================================================
  ! le segment j peut il etre obstacle de I
  !=================================================================================================
  ! on va eliminer un certain nombre de segments j de la recherche d'obstacles :
  ! on ne va rechercher comme obstacle seulement les segments de SG differents de
  ! celui de i qui peuvent conduire a une jonction (SG differents sauf les SG
  ! coplanaires) et les segments qui peuvent conduire a une annihilation

  ! les systemes coplanaires ne sont pas obstacle
  if (sysobst(Li,Lj) == izero) cycle DETECT

  ! i lui meme et ses rails de ne sont pas des obstacles
  if (j == i) cycle DETECT
  if (j == i1o) cycle DETECT
  if (j == i1e) cycle DETECT
  if (j == vnn1o) cycle DETECT
  if (j == vnn1e) cycle DETECT

  !si j est jonction avec un rail de i il n'est pas obstacle la jonction se zippe
  if (seg(i1o)%jonc .and. (j==seg(i1o)%ijonc)) cycle DETECT
  if (seg(i1e)%jonc .and. (j==seg(i1e)%ijonc)) cycle DETECT
  if (seg(vnn2o)%jonc.and. j == SEG(vnn2o)%ijonc) cycle DETECT
  if (seg(vnn2e)%jonc.and. j == SEG(vnn2e)%ijonc) cycle DETECT
  if (seg(vnn1o)%jonc.and. j == SEG(vnn1o)%ijonc) cycle DETECT
  if (seg(vnn1e)%jonc.and. j == SEG(vnn1e)%ijonc) cycle DETECT

  ! si vnn1 de i jonction avec un vnn1 de j, j n'est pas obstacle j est simplment un bras de la jonction
  if (seg(seg(j)%vnno)%jonc) then
    itemp = seg(j)%vnno
    if (itemp == seg(vnn1o)%Ijonc.or.itemp == seg(vnn1e)%Ijonc .or.itemp == seg(vnn2o)%Ijonc .or. &
        itemp == seg(vnn2e)%Ijonc) cycle DETECT
  endif

  if (seg(seg(j)%vnne)%jonc) then
    itemp = seg(j)%vnne
    if ((itemp == seg(vnn1o)%Ijonc).or.(itemp == seg(vnn1e)%Ijonc) .or.  itemp == seg(vnn2o)%Ijonc &
          .or. itemp == seg(vnn2e)%Ijonc) cycle DETECT
  endif


  VLj(:) = bveclin(:,Lj)

  !================================================================================================
  ! recherche du pt d'intersection entre i et J
  !================================================================================================
  ! ce pt ne depend que de la geometrie de i et j, la nature de la reaction
  ! n'est etudiee seult apres
  ! 2 cas sont possibles soit l'intersection entre i et j est une pt soit c'est un
  ! segment de droite, ce qui equivalent a savoir si j appartient ou non au PG de i
  if (OiOjscaN /= izero .or. OiEjscaN /= izero) then

    ! j n'appartient pas au PG de i, l'intersection entre j et PG de i est un pt
    jplani = .false.

    if (OiOjscaN /= izero .and. OiEjscaN /= izero) then

      ! ni Oj ni Ej n'appartiennent au PG de i il faut calculer le pt
      ! d'intersection nomme Cj
      Cj(1:3) = Oj(1:3)+abs(OiOjscaN)*SEG(j)%NORME*VLj(1:3)/(abs(OiOjscaN)+abs(OiEjscaN))

      ! on connait Cj, il faut verifier qu'il est bien dans l'aire balayee par i
      OiCj(:) = Cj(:)-Oi(:)
      EiCj(:) = Cj(:)-Ei(:)
      OiCjprojVDi = (OiCj(1)*VDi(1)+OiCj(2)*VDi(2)+OiCj(3)*VDi(3))/NormsqrVDi
      OiCjscaVDi1o = OiCj(1)*VDi1o_or(1)+OiCj(2)*VDi1o_or(2)+OiCj(3)*VDi1o_or(3)
      EiCjscaVDi1e = EiCj(1)*VDi1e_or(1)+EiCj(2)*VDi1e_or(2)+EiCj(3)*VDi1e_or(3)

      if (OiCjprojVDi <= zero .or. OiCjscaVDi1o < izero .or. EiCjscaVDi1e < izero) then

        cycle DETECT ! Cj est en dehors de l'aire de I

      else ! Cj est le pt obstacle

        distij = OiCjprojVDi
        OjCj(:)= Cj(:) - Oj(:)! necessaire pour les pb ss reseau

      endif

      if (kkdebug) write(379,*) 'PG j /= PG i, pt NbObsest Cj /= Oj ou Ej'

    else ! Oj ou Ej est pt d'intersection il faut trouver lequel

      if (OiOjscaN /= izero)then ! Ej est le pt d'intersection

        Cj(:) = Ej(:)
        OiCjprojVDi = OiEjprojVDi
        OiCjscaVDi1o = OiEjscaVDi1o
        EiCjscaVDi1e = EiEjscaVDi1e

      else ! Oj est le pt d'intersection

        Cj(:) = Oj(:)
        OiCjprojVDi = OiOjprojVDi
        OiCjscaVDi1o = OiOjscaVDi1o
        EiCjscaVDi1e = EiOjscaVDi1e

      endif ! on connait Cj, il faut verifier qu'il est bien dans l'aire balayee i

      if (OiCjprojVDi <= zero .or. OiCjscaVDi1o < izero .or. EiCjscaVDi1e < izero) then

        cycle DETECT ! Cj est en dehors de l'aire de I

      else ! Cj est le pt obstacle

        distij = OiCjprojVDi
        OiCj(:) = Cj(:) - Oi(:)
        OjCj(:) = Cj(:) - Oj(:)! necessaire pour les pb ss reseau

      endif

      if (kkdebug) write(379,*) 'PG j /= PG i, pt NbObsest Oi ou Ej'

    endif

    ! j appartient au PG de i
    !*****************************************************************************

  else

    ! 2 cas se presentent (j appartient a PG i) soit j est parallele a i soit
    ! il ne l'est pas et il faut rechercher le pt le plus proche
    jplani = .true.

    if (OiOjprojVDi /= OiEjprojVDi) then

      ! il existe un pt plus pres
      ! on recherche d'abord l'extremite de j la plus proche de i (suivant VDi)
      if (OiOjprojVDi < OiEjprojVDi) then ! oj est le pt le plus proche de i

        Aj(:) = Oj(:)! extemite de j la + proche de i
        OiAjprojVDi = OiOjprojVDi
        OiAjscaVDi1o = OiOjscaVDi1o
        EiAjscaVDi1e = EiOjscaVDi1e
        Dj(:) = Ej(:)! autre extremite
        OiDjprojVDi = OiEjprojVDi
        OiDjscaVDi1o = OiEjscaVDi1o
        EiDjscaVDi1e = EiEjscaVDi1e

      else ! ej est l'extremite la + proche de i

        Aj(:) = Ej(:)! extemite de j la + proche de i
        OiAjprojVDi = OiEjprojVDi
        OiAjscaVDi1o = OiEjscaVDi1o
        EiAjscaVDi1e = EiEjscaVDi1e
        Dj(:) = Oj(:)! autre extremite
        OiDjprojVDi = OiOjprojVDi
        OiDjscaVDi1o = OiOjscaVDi1o
        EiDjscaVDi1e = EiOjscaVDi1e
        VLj(:)=-VLj(:)

      endif

      ! l'extremite la +pres Aj est elle dans l'aire balayee de i?
      if (OiAjscaVDi1o >= izero .and. EiAjscaVDi1e >= izero) then! l'extremite la + pres est le pt d'intersection

        distij = OiAjprojVDi
        Cj(:) = Aj(:)
        OiCj(:) = Cj(:) - Oi(:)

      else

        ! Aj est en dehors de l'aire, le pt d'intersection est sur un rail
        surrail = .true.

        if (OiAjscaVDi1o < izero) then ! le pt d'intersection est sur le rail en Oi
            Cj(:) = Aj(:)+(abs(OiAjscaVDi1o)*SEG(j)%NORME*VLj(:))/(abs(OiAjscaVDi1o)+abs(OiDjscaVDi1o))
        else ! le pt d'intersection est sur le rail en Ei
            Cj(1:3) =Aj(1:3)+((abs(EiAjscaVDi1e)*SEG(j)%NORME*VLj(1:3))/(abs(EiAjscaVDi1e)+abs(EiDjscaVDi1e)))
        endif

      endif
      OiCj(:) = Cj(:) - Oi(:)
      distij = (OiCj(1)*VDi(1)+OiCj(2)*VDi(2)+OiCj(3)*VDi(3))/normsqrVDi
      OjCj(:) = Cj(:) - Oj(:)! necessaire pour les pb ss reseau

      if (kkdebug) write(379,*) 'j appartient au PG i, j secant a i'

    else ! j est colineaire a i
      distij = OiOjprojVDi

      ! y a t il recouvrement
      if ((OiOjscaVDi1o<=izero.and.OiEjscaVDi1o<=izero).or.(EiOjscaVDi1e<=izero.and.EiEjscaVDi1e<=izero)) then

        if (OiOjscaVDi1o==izero.or.EiOjscaVDi1e==izero) Cj(:) = Oj(:) ! j coupe l'aire balayee sur un rail
        if (OiEjscaVDi1o==izero.or.EiEjscaVDi1e==izero) Cj(:) = Ej(:)

        OiCj(:) = Cj(:) - Oi(:)
        OjCj(:) = Cj(:) - Oj(:)! necessaire pour les pb ss reseau

        if (kkdebug) write(379,*) 'j appartient au PG i'

      else

        jparali = .true.

        Cj(:) = (/-1,-1,-1/) ! signifie que y aura recouvrement
        if (kkdebug) write(379,*) 'j appartient au PG i, j // i'

      endif

    endif

    ! si j est un vnn2 de i, dans qlq cas particulier on le detect comme obstacle. j n est un obstacle
    ! seulment dans le cas ou la ligne de dislocation se replie sur elle meme
    if (j == vnn2o) then

      if (Cj(1)/=Ej(1).and.Cj(2)/=Ej(2).and.Cj(3)/=Ej(3)) then
        if (kkdebug) write(379,*) 'voisin est obstacles : ligne qui se replie sur elle meme'
      else
        cycle DETECT
      endif

    endif

    if (j == vnn2e) then
      if (Cj(1)/=Oj(1).and.Cj(2)/=Oj(2).and.Cj(3)/=Oj(3)) then
        if (kkdebug) write(379,*) 'voisin est obstacles : ligne qui se replie sur elle meme'
      else
        cycle DETECT
      endif

    endif

  endif

  if (kkdebug) then
    write(379,*) 'Cj(:) :',Cj(:)
    write(379,*) 'OiCj(:) :',OiCj(:)
    write(379,*) 'oiVDi:',OiCjprojVDi,'OiVDi1o :',OiCjscaVDi1o,'EiCVDi1e :',EiCjscaVDi1e
    write(379,*) 'i:',i,' j:',j,' a un pt d intersection'
  endif

  !==================================================================================================
  ! j est il obstacle?
  !==================================================================================================
  VLj(:) = bveclin(:,Lj) !redefinition necessaire car VLj peut avoir changer de signe lors calcul de Cj(:)

  if (distij > distiobs .or. distij <= zero) then

    cycle DETECT  ! j est + loin il n'est pas obstacle

  elseif (distij == distiobs) then

    NbObs = Nbobs + 1 ! j aussi proche de i que le precedent obstacle detecte
    ! qd il n'y a pas de pb de ss reseau on ne sais pas encore la nature de la reaction.
    if (nbobs > 100) then
        call seginfo(I,"  depassement d'obstacle pour le segment           ")
        stop
    endif

  elseif (distij < distiobs .and. distij > izero) then

    NbObs = 1 ! j est le nouvel obstacle le plus proche
    OBSTAB(1:50)%react = izero

  endif

  !==================================================================================================
  ! detection et gestion des pb de ss reseau
  !==================================================================================================
  pbssreso = .false.

  ! probleme de sous reseau en VDi
  if (.not. jparali) then !suivant si j est parallel a i le test pour les pb de ss reseau est different

    if (modulo((OiCj(1)*VDi(1)+OiCj(2)*VDi(2)+OiCj(3)*VDi(3)),(normsqrVDi*modep)) /= zero) pbssreso = .true. !pb ss reseau en VDi*modep

    if (kkdebug) write(379,*) 'modulo VDi - not jparali:',modulo(dot_product(OiCj,VDi),(normsqrVDi*modep))

  else

    if ((modulo((OiOj(1)*VDi(1)+OiOj(2)*VDi(2)+OiOj(3)*VDi(3)),(normsqrVDi*modep))) /= zero) pbssreso = .true.

    if (kkdebug) write(379,*) 'modulo VDi - jparali:',modulo(dot_product(OiOj,VDi),(normsqrVDi*modep))

  endif

  ! probleme de sous reseau en VLj
  if ((.not.pbssreso).and..not.jparali) then ! pas test ssreso VLi et VLj si

    normsqrOjCj = OjCj(1)*OjCj(1) + OjCj(2)*OjCj(2) + OjCj(3)*OjCj(3)
    normsqrVLj = VLj(1)*VLj(1) + VLj(2)*VLj(2) + VLj(3)*VLj(3)

    if (modulo(normsqrOjCj,normsqrVLj) /= zero) pbssreso = .true. !le modulo est calcule sur le carre de la
    ! norme car j'evite pb d arrondi du a la fct sqrt ou alors je peut faire le test du modulo sur chacun des 3
    ! composantes de ces 2 vecteurs (puisqu'il sont colineaires)

    if (kkdebug) write(379,*) 'modulo VLj',modulo(normsqrOjCj,normsqrVLj),normsqrOjCj,normsqrVLj

  endif

  ! probleme de sous reseau en VLi
  if ((.not.pbssreso).and..not.jparali.and.(.not.surrail)) then

    normsqrVLi = VLi(1)*VLi(1) + VLi(2)*VLi(2) + VLi(3)*VLi(3)
    VLi1o(1:3) = bveclin(1:3,Lio)

    if ((VLi1o(1)*VLi(1)+VLi1o(2)*VLi(2)+VLi1o(3)*VLi(3)) /= izero) then

        ! la projection de OiCj selon VLi il doit etre entier
        proj = ((OiCj(1)*VDi(1)+OiCj(2)*VDi(2)+OiCj(3)*VDi(3))/((VLi1o(1)*VDi(1)+VLi1o(2)*VDi(2)+VLi1o(3)*VDi(3))))

        ! puisque pas de pb de ss reseau jusque ici en VDi jusque ici
        OipCjscaVLi = (OiCj(1)*VLi(1)+OiCj(2)*VLi(2)+OiCj(3)*VLi(3))-proj*(VLi1o(1)*VLi(1)+VLi1o(2)*VLi(2)+VLi1o(3)*VLi(3))

        if (OipCjscaVLi/=zero) then
          if (modulo(OipCjscaVLi,normsqrVLi)/=zero) pbssreso=.true.
          if (kkdebug) write(379,*) 'modulo VLi :',modulo(OipCjscaVLi,normsqrVLi)
        endif

        !  oip(:) = proj * VLj(:) + Oi(:) ! calcul du pt oip
    else ! recherche ss reseau adapte qd VLi et VLi1o sont orthogonaux

        if (modulo((OiCj(1)*VLi(1)+OiCj(2)*VLi(2)+OiCj(3)*VLi(3)),normsqrVLi)/=zero) pbssreso=.true.
        if (kkdebug) write(379,*) 'modulo VLi :',modulo((OiCj(1)*VLi(1)+OiCj(2)*VLi(2)+OiCj(3)*VLi(3)),normsqrVLi)
    endif

  endif

  if (.not.pbssreso ) then
    if (seg(j)%jonc) jjonc = .true. ! j est deja une jonction
  endif

  ! ATTENTION on regarde les pb de ss reseau dans la detection
  ! si i petit et j obstacle + pb ss reso ou j jonc=>j n'est  plus obs !
  if (jjonc)    OBSTAB(nbobs)%react = icinq     !pb j deja jonc
  if (pbssreso) OBSTAB(nbobs)%react = iquatre   !pb de ss reseau

  ! final list of real obstacle is save in OBSTA
  if (nbobs > 0 .and. TakeObst) then
    OBSTAB(nbobs)%numseg = j
    OBSTAB(nbobs)%ptobs  = Cj(:)
    if (jplani) OBSTAB(nbobs)%plani = .true.
    distiobs = distij
  endif

  ! For calculations with two pass it is useful to build up a MiniListeSeg
  ! with a tiny ListSegTester to use on each procs in the second pass
  MiniNsegTester = MiniNsegTester + 1
  if (MiniNsegTester > dim_alloc) then
    write(379,*) 'Hoops in the detect loop, MiniNsegTester = ', MiniNsegTester, ' and dim_alloc = ', dim_alloc
    stop 'Problem with the dim_alloc and MiniNsegTester values inside the detect loop!'
  endif
  MiniListeSeg(MiniNsegTester) = int(J,4)

enddo DETECT
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if (kkdebug) then
  write(379,*) 'Firsttime =',Firsttime,'  Nbobs =', Nbobs
  write(379,*) 'MiniNsegTester =',MiniNsegTester
  write(379,*) 'MiniListeSeg =',MiniListeSeg(1:MiniNsegTester)
  write(379,*) ' '
endif

#ifdef PA

  ! We must gather and reduce the information we got on each procs before the end of the subroutine
  if (ParaComNeeded) then

    MiniNseg_Send(MON_RANG_PAIRIMPAIR + 1) = MiniNsegTester  ! The number of MiniNsegTester found in each proc

    ! We first reduce the number of MiniNsegTester
    CALL MPI_ALLREDUCE(MiniNseg_Send(1),MiniNseg_Recv(1),TAILLE_PAIRIMPAIR,MPI_INTEGER,MPI_SUM,COMM_PAIRIMPAIR,IERR)

    MiniNsegTester = sum(MiniNseg_Recv(1:TAILLE_PAIRIMPAIR))     ! The total number of MiniNsegTester found

    if (MiniNsegTester /= 0) then

      if (MiniNsegTester > dim_alloc) then
        write(379,*) 'Hoops, MiniNsegTester = ', MiniNsegTester, ' and dim_alloc = ', dim_alloc
        stop 'Problem with the dim_alloc and MiniNsegTester values !'
      endif

      ! The MiniListe_Send array is build and contains elements coming only from procs where obstacles were found
      if (MiniNseg_Send(MON_RANG_PAIRIMPAIR + 1) /= 0) then

        ! location of the information coming from each procs
        if (MON_RANG_PAIRIMPAIR /= 0) then
          array_deb = sum(MiniNseg_Recv(1:MON_RANG_PAIRIMPAIR)) + 1
        else
          array_deb = 1
        endif

        array_fin = array_deb + MiniNseg_Recv(MON_RANG_PAIRIMPAIR + 1) -1

        ! The corresponding array of MiniListeSeg
        MiniListe_Send(array_deb:array_fin) = MiniListeSeg(1:MiniNseg_Send(MON_RANG_PAIRIMPAIR + 1))

      endif

      ! Then we reduce the MiniListeSeg arrays to distribute the information
      CALL MPI_ALLREDUCE(MiniListe_Send(1),MiniListe_Recv(1),MiniNsegTester,MPI_INTEGER,MPI_SUM,COMM_PAIRIMPAIR,IERR)

      ! The final MiniListeSeg indentical on each procs
      MiniListeSeg(1:MiniNsegTester) = MiniListe_Recv(1:MiniNsegTester)

      if (kkdebug) then
        write(379,*) 'ParaComNeeded =',ParaComNeeded
        write(379,*) 'Firsttime =',Firsttime,'  Nbobs =', Nbobs
        write(379,*) 'MiniNsegTester =',MiniNsegTester
        write(379,*) 'MiniListeSeg =',MiniListeSeg(1:MiniNsegTester)
        write(379,*) ' '
      endif

    endif

  endif

#endif

end subroutine Detect_loop

end module CONTACT
!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
