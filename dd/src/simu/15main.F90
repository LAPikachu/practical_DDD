
!===================================================================================================
!========================    DEBUT  OF THE MAIN MODULE  ============================================
!===================================================================================================

include '00copyright.f90'

!> \brief This is the main program bloc.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English

#ifndef GX

program MICROMEGAS

#else

subroutine MICROMEGAS

#endif

use BRICAMAT
use VARGLOB
use DEBUG
use INIT
use microstructure
use ELASTI
use DYNAM
use TOPOLO
use CONTACT
use RESUL
use BIGSAVE
use nucleation
use initconftools

! Most subroutines needed by MDC
#ifdef MDC

use mdctools

#endif

#ifdef GX

use INTERGRA

#endif

!######################

implicit none

!######################

#ifdef GX
integer(DPI)    :: eqzer
#endif

#ifndef MDC
integer(DPI)        :: epspt
#endif

#ifdef PA
character (len=30) :: file_debug_loop     !< Name of the debug_loop file
character (len=2)  :: file_index          !< Index of the debug_loop file in parallel computation
#endif

integer(DPI)        :: KK0
integer             :: heure(8)
character (len =10) :: x(3)
real(DP)            :: rseconds,freq1,freq3

!! variable to manage arguments attached to the mM executable files
integer             ::narg
character(len=20)   ::name
character(len=256)  :: cmd
logical             ::file_exist



heure(1:8) = izero

Print *, " "
Print *, " "
Print *, "                                    MM    M  EEEE   GGGG    AAA    SSSS"
Print *, "                                    M M M M  EE    G       A   A  SS   "
Print *, "         MM MM  I  CCC  RRRR  OOOO  M  M  M  EEEE  G  GGG  AAAAA   SSS "
Print *, "         M M M  I  C    RRR   O  0  M     M  EE    GG  GG  A   A     SS"
Print *, "         M   M  I  CCC  R  R  0000  M     M  EEEE   GGGG   A   A  SSSS "
Print *, " "
Print *, " "

Print *, " ========================================================================================"
Print *, "                                CopyRight"
Print *, " "
Print *, " mM (for microMegas) is an open source program of DD (Dislocation Dynamics) simulation"
Print *, " originally developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE."
Print *, " Copyright (C) 1993, 1996, 2000, 2002, 2004, 2013   B. Devincre, L. Kubin, M. Condat,"
Print *, " C. Lemarchand, R. Madec, S. Groh, G. Monnet, J. Durinck, P. Carrez, C. de Sansal,"
Print *, " M-G. Tagorti, S. Queyreau., S. Naamane, A. Roos, A. Vattre, R. Gatti, O. Jamond,"
Print *, " F. Boioli, A. Roos, L. Korzeczek."
Print *, "  "
Print *, " This is a free software, and you are welcome to redistribute it under the"
Print *, " conditions of the GNU General Public License (see the 'Copyright' file in"
Print *, " the program distribution)."
Print *, " ========================================================================================="
Print *, "              "

!#########################################################################
!# Etapes de preparation de la simulation : Module INIT                  #
!#########################################################################

!!! looking for arguments attached to MicroMegas
narg=command_argument_count()

!! -------------------------------------------------------------------------
!! (DCM) In the new DCM algorihtm the time step of the simulation
!! is defined in Zebulon (see .inp file)
!! mM knows the time step but not the total number of step (Zebulon rules!!)
!! -------------------------------------------------------------------------

call LOAD               ! Lecture des donnees, tabulation...

if(narg>0)then

   if(narg > 1) stop "only one argument is permitted in Micromegas executables"

   call get_command_argument(1,name)

   select case(adjustl(name))
!First known args
    case("--domvtx") !to compute geometry and volume of domains defined in b_plan_conc
      domvtx=.TRUE.
      write(*,*) 'WARNING: Argument domvtx activated'
      call configuration
      call load_boundaries
      call vertex_calculations
      INQUIRE(FILE="../bin/volume.py", EXIST=file_exist)
      if (file_exist) then
        !write(cmd,'("../bin/volume.py")')
        cmd = "../bin/volume.py"
      else
        !write(cmd,'("volume.py")')
        cmd = "volume.py"
      endif
      call system(cmd)
      stop "domvtx: geometry computed and stored in ../out/simulationvolume.vtk!!"
    case("--initconf")
      initconf=.TRUE.
      write(*,*) '*** WARNING: Argument initconf activated ***'
      write(*,*) '*** New procedure to generate initial distribution of dislocation sources ***'
      call microconf
      stop
    case default
     write(*,*) ''
     write(*,*) 'STOP!!!! the option ', adjustl(name) ,'is not valid'
      write(*,*) ''
      stop
    end select
  endif

if (kkdebug) then
  write(379,*) "apres load "
  CALL FLUSH(379)
  !call segsinfo (" avant configuration de segments")
endif

call configuration

call initialisation_fichiers
! lecture and modification-treatment of the initail configuration
! of the simualtion. It belongs to the module : segments

!! -----------------------------------------------------------------------
!! (DCM) the area swept by initial dislocation loops is sent to Zebulon
!! to get the initial guess of FE stress field
!! (see subroutine init_plas_eigen called inside generation micorstructure)
!! ------------------------------------------------------------------------

call generation_microstructure

write(*, '("  Number of segments before net procedure:  ",I10)') Nsegm
call NET
write(*, '("  Number of segments after net procedure:  ",I10)') Nsegm

#ifdef MDC

!call time_step !Mm gets simulation time step from Zeb

call initialize_MDC_var

call allocation_initiale
! kk = 0

if (mode_deformation_key) then
call Init_Solli
endif

#else

call allocation_initiale

! Stat pour le calcul de la vitesse de deformation pour le controle

call Init_Solli

if (nsegm < 100) call disloinfo ("avant BOUCLE PRINCIPALE                         ")

! The nucleation procedure
if (key_nucleation) then
  call periodic_nucleation_input
  tot_loop=IZERO
  if (deltat_newloop < DELTAT) write(*,*)'Error  deltat_newloop is lower than DELTAT0'
endif
#endif
!=============================================================================
!KISAUVE=10    ! One may need to restart a simulation to debbug
!KSTAT=1       ! but with a different periodicity of writing
!KKIM=1        ! in the outputs of the simulation
!KPREDRAW=1
!=============================================================================

KK0 = KK

freq1 = deltat/normdep(1)
freq3 = deltat/normdep(3)
call date_and_time(x(1),x(2),x(3),heure)
rseconds = heure(6)*60 + heure(7) + real(heure(8))/1000.0 - 0.1

if (kk==IZERO) then
  print *, " "
  print *, ">>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<"
  print *, "  Begining of the simulation time loop  "
  print *, ">>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<"
  print *, " "
  write(*,'("      Date  : ",I2,"/",I2,"/",I4)') heure(3), heure(2), heure(1)
  write(*,'("      Heure : ", I2 ,"H" ,I2, " et " ,I2, "." ,I3, " seconds")') &
          heure(5), heure(6),heure(7),heure(8)
  print*, " "
endif

!**************************************
!**** Boucle generale sur le temps ****
!**************************************

#ifdef MDC

!! --------------------------------------------------------------------------------------
!! (DCM) during the first iteration we need to call Zebulon to get the stress field of
!! the initial configuration (kk=0) to compute forces at the step kk=1 (so zcall is set to
!! true before entering the do while loop).
!! Generally speaking, in the DCM algorithm, forces are evaluated at the step KK using
!! the FE stress field solution at the end of the step KK-1
!! --------------------------------------------------------------------------------------
zcall=.true.

BP1 : do while (endcalcul==0)

#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
#endif

  if(zcall) then
  write(*,*)
  write(*,*)'##################### New DCM time step ###########################'
  endif
#if defined(MDC) && defined(PA)
 endif
#endif

  kk = kk + 1

#else
BP1 : do KK = kk0 + 1, NSTEP

#endif


! steps of the time loop are written in a debug file when debug mode is active
#ifdef PA
  if(iterinfo > izero) then
    open(10,STATUS='SCRATCH')
      write  (10,'(I2)') Mon_Rang
      rewind (10)
      read   (10,*) file_index
    close(10)
    file_debug_loop= "../out/debug/test_boucle"//trim(file_index)//".txt"
    open(27,FILE=file_debug_loop,STATUS='unknown')
    write(27,*) "  Iteration numero " , kk
  endif
#else
  if(iterinfo > izero) then
    open(27,FILE='../out/debug/test_boucle.txt',STATUS='unknown')
    write(27,*) "  Iteration numero " , kk
  endif
#endif

  if (kk == IUN .and. relax_TL /= izero) then
        print*, " "
        print*, "Relaxation with only line tension, KK = ",kk
        print*, " "
  endif

  if (kk == Relax_TL+1  .and. relax_INT /= relax_tl) then
        print*, " "
        print *, "Relaxation with elastic interactions, KK = ",kk
        print*, " "
  endif

  if (kk == Relax_Int+1 .and. relax_reac /= relax_int) then
        print*, " "
        print *, "Relaxation with contact reactions, KK = ",kk
        print*, " "
  endif

  ! The debug activation
  kkdebug = (kk == iterinfo)

  if (kkdebug) then
    write(379,*) 'At the begining of the step nsegm, ntgd1, ntgd2, oldntgd1 and oldntgd2 = ' &
                  ,nsegm,ntgd1,ntgd2,oldntgd1,oldntgd2
    if(nsegm < 20000) call disloinfo ("---> Main loop beginning <----")
    call check_connec (" Main loop beginning ")
  endif

  ! on calcul les stat de carto sans (1) et avec (3) contrainte apres la relaxation
  ! en tension de ligne et a la fin de la simulation
  if ((kk == relax_TL .or. KK == NSTEP) .and. cartograph == ideux) call stat_carto
  if (kkdebug) then
    write(379,*) "  apres stat_carto "
    CALL FLUSH(379)
  endif


#ifdef MDC

! loading control procedures
if (mode_deformation_key) then
call Solli
endif

!! ---------------------------------------------------------------------
!! (DCM) loading conditions are defined in Zebulon control file (...)
!! ---------------------------------------------------------------------

  !*********************************
  !*** Conditions de chargement  ***
  !*********************************
#else

!   if (mode_deformation /= IHUIT) then

    call SOLLI

!   endif

  if (iterinfo > izero) write(27,*) "  after SOLLI "
  if(kkdebug) then
    if (nsegm < 20000) call disloinfo ("    after SOLLI                    ")
    CALL FLUSH(379)
  endif

#endif
  !****************************************************************
  !*** Si interfacage graphique : Appel au programme graphique  ***
  !****************************************************************

#ifdef GX

  eqzer=0

  if((MOD(kk-1,KPREDRAW).eq.IZERO) .or. kkdebug) then

    ! The graphical interface call to visualize the dislocation microstructure
    call MODEGRAF(eqzer)

    if (kkdebug) call mes(" after MODEGRAF ")

  endif

#endif

  !****************************************************************
  !*** DISCRETI mainly cut segments when needed by discretization *
  !****************************************************************
  if(kkdebug) then
    if (nsegm < 20000) call disloinfo ("    apres discreti  avant force")
    write(379,*) "  "
    CALL FLUSH(379)
  endif



  call DISCRETI

  if(iterinfo > izero) write(27,*) "  after  DISCRETI "


  if(kkdebug) then
    if (nsegm < 20000) call disloinfo ("    apres discreti  avant force                     ")
    write(379,*) "  "
    call check_connec ("   apres discreti  avant force                  ")
    CALL FLUSH(379)
  endif

  !*****************************************************************
  !*** Force : - Calcul de la force effective sur chaque segment ***
  !***         TAU.b.L (mu.a2)                                   ***
  !***         - Calcul de la cission resolue dans differentes   ***
  !***         directions de deplacement                         ***
  !*****************************************************************

  !!--------------------------------------------------------------------------------------
  !!(MDC) if Zcall = true in subroutine FORCE we get the information about FE stress field
  !!   at the center of segments with length /= 0
  !!--------------------------------------------------------------------------------------

  call FORCE

  if (kkdebug) then
    ! call disloinfo ("    apres  force  avant depredi                    ")
    call check_connec ("   apres force  avant depredi                  ")
    CALL FLUSH(379)
  endif

  !*************************************************************
  !*** Prediction des deplacements en fonction               ***
  !***               de la loi des vitesses                  ***
  !*************************************************************

  call DEPPREDIC

  if (iterinfo > izero) write(27,*) "  after DEPPREDIC "

  if (kkdebug) then
    !       call disloinfo ("    apres   depredi                    ")
    call check_connec ("   apres  depredi                  ")
    !       call seginfo(2730," probleme")
    !       call seginfo(3149," probleme")
    !       stop
  endif


#ifdef MDC
  !!--------------------------------------------------------------------------------------
  !!!! (DCM) Zebulon time step starts here !!!!
  !! the zcall switch is used to call Zebulon. !!!
  !!if Zcall = true the areas swept by dislocation are transferred to Zebulon !!
  !! The FE code solves the elastic problem, and then the computed stress field is used to
  !! calculates forces on dislocation segments at the step KK+1
  !! BE CAREFUL !!
  !! in the new version of the DCM Zebulon can be called with a given periodicity
  !!(MDC_timestep_ratio variable  in mM sources, **timestep_ratio in the DCM .inp input file)
  !! ----------------------------------------------------------------------------------------

   if (modulo(kk, mdc_timestep_ratio)==0) then
      zcall=.true.
   else
      zcall=.false.
   end if
#endif

  !***********************************************************
  !*** Recherche des obstacles, du plus proche obstacle    ***
  !*** traitement de reactions de contacts, deplacements.  ***
  !***********************************************************
  !! (DCM) if Zcall = True the area swept are sent to Zebulon.
  call UPDATE

  if(iterinfo > izero) write(27,*) "  after UPDATE "
  if(kkdebug) then
    if (nsegm < 20000) call disloinfo ("   after UPDATE                       ")
    !call seginfo (1,"   apres update                       ")
    !call check_connec ("  apres update                    ")
    CALL FLUSH(379)
  endif

  !********************************************************
  !*** Solve problems of connectivity on segs at O or E ***
  !*** 1. Check the junction status                     ***
  !*** 2. Eliminate not necessary GD segments           ***
  !*** ... plus many other things                       ***
  !********************************************************

  call corriger_config

  if(iterinfo > izero) write(27,*) "  apres corriger_config "
  if(kkdebug) then
    if (nsegm < 20000) call disloinfo ("   apres corriger_config  ")
    call check_connec ("  apres corriger_cofig    ")
    CALL FLUSH(379)
  endif

  !*****************************
  !* Procedure NET: cleanup  ***
  !*****************************

  call NET

  if(iterinfo > izero) write(27,*) "  apres NET "
  if (kkdebug) then
    if(nsegm < 20000) call disloinfo ("   - apres NET --                       ")
    call check_connec ("   - apres NET --                   ")
    CALL FLUSH(379)
  endif

  ! The nucleation procedure is periodically called if needed
  if (key_nucleation ) then
    accutime = KK * deltat
    if (kk == izero) last_newdislo_t = deltat_newloop

    if ((accutime - last_newdislo_t) >= deltat_newloop) then
      print *, 'call periodic nucleation ---> nsegm = ', nsegm
      call periodic_nucleation
      last_newdislo_t = accutime
      print *, 'exit periodic nucleation <--- nsegm = ', nsegm
    end if
  end if

!***********************************************************************
!*** Calcul des quantites moyenne (deformation, densite de disloc ...) *
!***********************************************************************

#ifdef MDC
  call STATISMDC
#else
  call STATIS
#endif

  if(iterinfo > izero) write(27,*) "  apres stats "
  !call checkit("*** VERIFICATION checkit ***")

  !********************************************************************
  !*** Sauvegarde de la configuration tous les "KKIM" iterations ***
  !***  pour faire un film de la trajectoire des dislocations    ***
  !********************************************************************

#ifndef MDC

  if (mod(KK,kkim) == izero .or. KK == iun) then

    call SAUVima

  endif

#endif

  !************************************************************
  !*** Calcul du nombre de jonction et de glissement devies ***
  !************************************************************

  call reacloc

  if(iterinfo > izero) write(27,*) "  apres reacloc "

  !****************************************************
  !*** Elementary information printed on the screen ***
  !****************************************************
#ifndef MDC
  if (modulo(kk,Kpredraw) == 0 .or. kk == iun) then
    !    if (modulo(kk,IUN) == 0) then
    epspt = nint(100.0*eps_pl_inst/deltat/epsilonpoint)
    ! D_brute(numero_loi(1)) : average of predicted displacement of screw (over all screw segments)
    ! D_Vraie(numero_loi(1)):  average of effective displacement of screw (over all screw segments)
    ! D_brute(numero_loi(3)):  average of predicted displacement of non-screw (over all non-screw segments)
    ! D_Vraie(numero_loi(3)):  average of effective displacement of non screw (over all non-screw segments)

    if(Numero_loi(1) /= Numero_loi(3)) then
      write(*,10) kk,nsegm,epspt,epso*100.0,sigma*xmu*1D-6*schmid, &
                    D_brute(numero_loi(1)),D_Vraie(numero_loi(1)), &
                    D_brute(numero_loi(3)),D_Vraie(numero_loi(3)), nbjonction,ntgd,ntfs
    else

      if (modulo(kk,20*Kpredraw) == 0.or. kk == IUN) Then
        write(*,*) ' Steps  Segms  Epsidot   Epsi(%)    Sig(MPa)   Tau(MPa)     DepBrut     DepTrue   N_Jonc  N_CS  N_Long  N_FS'
      endif

      if (mode_deformation /= 8) then
        write(*,11) kk,nsegm,epspt,epso*100.0,sigma*xmu*1D-6,sigma*xmu*1D-6*schmid, &
                    D_brute(numero_loi(2)), D_vraie(numero_loi(2)), nbjonction, ntgd, Nb_GroSeg, ntfs
      else
        TestSig(:,:) = stress_fact * SigappS(:,:)

        write(*,11) kk,nsegm,epspt,epso*100.0,maxval(TestSig)*xmu*1D-6,maxval(TestSig)*xmu*1D-6*schmid, &
                    D_brute(numero_loi(2)), D_vraie(numero_loi(2)), nbjonction, ntgd, Nb_GroSeg, ntfs
      endif

    endif

  endif

10  format(I6,"; ",I6,"s; ",I4,"%; ",ES10.2E2,"%;",ES10.2E2,"MPa;DBVV=",ES10.2E2,"/",ES10.2E2, &
         ";DBVC=",ES10.2E2,"/",ES10.2E2,"; ",I5," : ",I5," : ",I5)
11  format(I7," ",I6,"  ",I5,"  ",ES10.2E2,"  ",ES10.2E2," ",ES10.2E2,"   ",ES10.2E2,"  ",ES10.2E2,"  ",I5,"  ",I5,"  ",I5,"  ",I5)

  !******************************************
  !*** The segment configuration is saved ***
  !******************************************

  if(MOD(KK,KISAUVE).eq.IZERO) then

      call saveall
      call save_config

  endif
#endif

  !**************************************************************************
  !*** Validation of the segment configuration before we start a new step ***
  !**************************************************************************
  Plusseg=izero            ! Initialisation pour parano
  call check_connec ("      fin de la boucle principale                 ")
  if(iterinfo > izero) write(27,*) "  apres check_connec "

  !******************************************
  !*** The multipoles domaines allocation ***
  !******************************************
  if (allocation_dynamique_boites) call allocation_dynamique
  if(iterinfo > izero) write(27,*) "  apres allocation_dynamique "

  if(iterinfo > izero) close (27)
  if(kkdebug) then
    CALL FLUSH(379)
    stop
  endif

#ifdef MDC
  ! Ask Z what to do next
#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
#endif
  if (zcall) then
      write(*,*)'=======================    mM OUTPUT  =============================='
      write(*,*) ' Steps Segms DepBrut DepTrue  N_Jonc  N_CS  N_FS'
  endif
#if defined(MDC) && defined(PA)
 endif
#endif

  if (modulo(kk,kpredraw) == 0 .or. kk == iun) then
    write(*,12) kk,nsegm, &
         D_brute(numero_loi(2)), D_vraie(numero_loi(2)), nbjonction,ntgd,ntfs
12  format(I7," ",I7," ",F6.3,"  ",F6.3,"  ",I5,"  ",I5,"  ",I5)
  endif


  if(MOD(kk,KISAUVE).eq.IZERO) then
    call saveall
    call save_config
  endif

  if( mod(KK,kkim) == izero .or. KK == iun) then
    call SAUVima
  endif

  ! Z says continue (DCM time step end here --> linear problem solved)
  if (zcall) call whattodo ! Contains Z ->DD No. XX

#endif

  if (kkdebug) then
    if(nsegm < 20000) call disloinfo ("   - fin --                       ")
    CALL FLUSH(379)
  endif

enddo BP1

!***************************************
!**** Fin de la boucle sur le temps ****
!***************************************

call desallocation

call date_and_time(x(1),x(2),x(3),heure)

rseconds = heure(6)*60 + heure(7) + real(heure(8))/1000.0 - 0.1

#if defined(PA) && !defined(MDC)

!=================================
!===   Fin du parallelisme  ======
!=================================

print*,"cloture de MPI"

CALL MPI_FINALIZE(IERR)

if (IERR /= 0) stop "Problem with MPI ending"
#endif

print *, " " ; print *, " "
print *, " Fin des simulations ................................................."
print *, " " ; print *, " "
write(*,'(" Date  : ",I2,"/",I2,"/",I4)') heure(3), heure(2), heure(1)
write(*,'(" Heure : ", I2 ,"H" ,I2, " et " ,I2, "." ,I3, " seconds")') &
            heure(5), heure(6),heure(7),heure(8)
write(*,*) " En secondes : ", rseconds
print *, " "

#ifndef GX

#ifdef MDC
call MPI_Barrier(MPI_COMM_WORLD,ierr)
print *,'here'
call MPI_Finalize(zebmpierr)
if (zebmpierr /= 0) stop "Problem with MPI ending 2"
#endif

end program MICROMEGAS

#else

end subroutine MICROMEGAS

#endif

!##########################################################################
!# subroutine Infos                                                       #
!########################################################## 27/10/00 ######
subroutine infos

use debug

#ifdef GX
use INTERGRA
#endif
implicit none

integer(kind=DPI) :: I,toi(3),gogoo,gogoe
character :: KARAK

real(8) :: totoUaccess
!*** Pour eviter ce type de message etrange sur les DEC :
!    Unaligned access pid=8803 <Pm_Dec_GX>
!     va=0x14272e19c pc=0x3ff81a091b0 ra=0x3ff81a091a8 inst=0x9c090000

KARAK='h'

do while (KARAK.ne.'q'.and.KARAK.ne.'Q')

    write (*,*) 'Infos>';read (*,*) KARAK;

    if (KARAK.eq.'h'.or.KARAK.eq.'H') then
       write (*,*) "P : portion de boucle"
       write (*,*) "T : tout les segments"
       write (*,*) "S : infos generales sur un Segments"
       write (*,*) "D : infos sur la dynamique d'un segment"
    endif

    if (KARAK.eq.'p'.or.KARAK.eq.'P') then
       write (*,*) 'de I a J :'
       read (*,*) gogoo
       read (*,*) gogoe
       call confie(gogoo,gogoe)
    endif

    if (KARAK.eq.'t'.or.KARAK.eq.'T') then
       do i=1,nsegm,1
           write (*,*) '{',i,'}',&
           seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%voise, &
           seg(i)%jonc,seg(i)%ijonc,seg(i)%gd,seg(i)%vnno,seg(i)%vnne
       enddo
    endif

    if (KARAK.eq.'s'.or.KARAK.eq.'S') then
       write (*,*) 'Segment n >';read (*,*) I;
       write (*,*) '{',i,'}',&
       seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%voise, &
       seg(i)%jonc,seg(i)%ijonc,seg(i)%gd,seg(i)%vnno,seg(i)%vnne

!**** CHANGEMENT DE REPERE
       tOi(:) = seg(i)%o(:)+seg(i)%norme/2*bVecLin(:,seg(i)%veclin)-modur(:)/2
       seg(1:nsegm)%o(1)=modulo(seg(1:nsegm)%o(1)-toi(1),modur(1))
       seg(1:nsegm)%o(2)=modulo(seg(1:nsegm)%o(2)-toi(2),modur(2))
       seg(1:nsegm)%o(3)=modulo(seg(1:nsegm)%o(3)-toi(3),modur(3))
       write (*,*) '{',i,'}',&
       seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%voise, &
       seg(i)%jonc,seg(i)%ijonc,seg(i)%gd,seg(i)%vnno,seg(i)%vnne
#ifdef GX
!  call MODEGRAF(i)
#endif
!**** RETOUR AU REPERE INITIAL
       seg(1:nsegm)%o(1)=modulo(seg(1:nsegm)%o(1)+toi(1),modur(1))
       seg(1:nsegm)%o(2)=modulo(seg(1:nsegm)%o(2)+toi(2),modur(2))
       seg(1:nsegm)%o(3)=modulo(seg(1:nsegm)%o(3)+toi(3),modur(3))

    endif

    if (KARAK.eq.'d'.or.KARAK.eq.'D') then
       write (*,*) 'Deplacement de>';read (*,*) I;
       write (*,*) '{',i,'}'
       write (*,*) 'IDEP: ',idep(i)
       write (*,*) 'tautot:',tautot(i)
       write (*,*) 'TauTL:',TauTL(i)
       write (*,*) 'TAUap:',TauapP(i)
       write (*,*) 'Daire:',daire(i)
       totoUaccess = seg(i)%resdep
       write (*,*) '%RESDEP:',totoUaccess
    endif

enddo

end subroutine infos


!===================================================================================================
!========================    END OF THE MAIN MODULE  ===============================================
!===================================================================================================
