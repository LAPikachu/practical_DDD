
!===================================================================================================
!========================    DEBUT    MODULE     ===================================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to dislocation loop nucleation.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module nucleation

use constantes
use BRICAMAT
use VARGLOB
use VARBASE
use INIT
use CONNEC
use microstructure
use elasti
use initconftools

implicit none

real(kind=DP)                              :: deltat_newloop          !< New loop nucleation time test periodicity
real(kind=DP)                              :: last_newdislo_t
real(kind=DP)                              :: loopcritR               !< Initial critical loop radius (A unit)
real(kind=DP)                              :: refstress
real(kind=DP)                              :: Delta_loopcritR         !< Loop radius increment after a large number of unsuccessful nucleation events
real(kind=DP)                              :: Delta_refstress
real(kind=DP)                              :: echellesim
real(kind=DP)                              :: VecBurgerssim

integer(kind=DPI),allocatable              :: h_boxx(:),h_boxy(:),h_boxz(:) !< size of the region in X,Y,Z direction where nucleation events are tested
integer(kind=DPI),allocatable              :: box_position(:,:)       !< Origin of the box where nucleation events are tested
integer(kind=DPI)                          :: nloop                   !< Number of new loops to consider at each nucleation iteration
integer(kind=DPI)                          :: tot_loop
integer(kind=DPI)                          :: nbregion                !< Number of region to include/exclude for nucleation
integer(kind=DPI)                          :: kunsuccess_max          !< Maximum number of iterations for nucleation try at a given loop size
integer(kind=DPI),parameter                :: nsegloop = 8            !< number of segments forming the loop nucleus
integer(kind=4), allocatable               :: seedn(:)                !< a list of seeds
integer(kind=4)                          :: seedns                  !< the seed for the nucleation procedure

character(len=40)                          :: loop_nucleated

logical, allocatable                       :: inoroutside(:)          !< Logical for region definition (true = if inside, false = outside)

!******************************************
               contains
!******************************************

!#############################################################################
!#periodic nucleation  of loops : reading input parameters for dislocation
!#loop insertion
!#############################################################################
subroutine periodic_nucleation_input

implicit none

character(len=40)       :: file_Periodic_nucleation
integer(kind=DPI)       :: ii
integer(kind=4)         :: i

last_newdislo_t = zero

file_Periodic_nucleation = "../in/nucl_def"

open(2,file=file_Periodic_nucleation,STATUS='OLD')
read(2,*) deltat_newloop        !new loop nucleation test periodicity
read(2,*) loopcritR             !initial critical loop radius (A unit)
read(2,*) Delta_loopcritR       !loop radius increment after a large number of unsuccessful nucleation events
read(2,*) nbregion              !Number of region to include/exclude for nucleation
allocate (inoroutside(nbregion))
allocate (h_boxx(nbregion))
allocate (h_boxy(nbregion))
allocate (h_boxz(nbregion))
allocate (box_position(3,nbregion))
do ii=1,nbregion
  read(2,*) inoroutside(ii)           !logical for region definition (true = if inside, false = outside)
  read(2,*) h_boxx(ii)                !size of the box (X direction) where nucleation events are tested
  read(2,*) h_boxy(ii)                !size of the box (Y direction) where nucleation events are tested
  read(2,*) h_boxz(ii)                !size of the box (Z direction) where nucleation events are tested
  read(2,*) box_position(:,ii)           !origin of the box where nucleation events are tested
enddo
read(2,*) nloop                 !number of new loops to consider at each nucleation iteration
read(2,*) kunsuccess_max        !maximum number of iterations for nucleation at a given loop size
read(2,*) seedns            ! seed size to choose nucleation position

close(2)

print *, ' '
print *, 'Time interval for testing a new loop', deltat_newloop
print *, 'Initial loop radius (Angstrom)', loopcritR

! remetre les rayons d'inetraction en avalue
if (Microboucle >  (4 * loopcritR * 1.e-10 / avalue) ) then
  Microboucle = INT(4 * loopcritR * 1.e-10 / avalue,DPI)     ! microboucle en a
  write(*,*) ' '
  write(*,*) 'The microbloucle size is redefined to the value =',Microboucle, loopcritR, avalue
endif

! The file with the nucleated loop info is created
loop_nucleated = "../out/loop_nucleated"
open(22,file = loop_nucleated, STATUS='REPLACE')
close(22)

! Initialization of the random number generator
allocate(seedn(seedns))
seedn = seedns + 123 * (/(i, i = 1,seedns)/)
call random_seed(put=seedn)
!print*,'seedn',seedn
deallocate(seedn)

end subroutine periodic_nucleation_input


!#############################################################################
!#periodic nucleation of loops : dislocation loop insertion
!#############################################################################
subroutine periodic_nucleation

implicit none

integer(kind=DPI)     :: loopSEG(nsegloop,5), loop_edge, loop_slipsystem, nloopplus
integer(kind=DPI)     :: itemp,jj, ii, ikk, kunsuccess,loop_radius,dom
integer(kind=DPI)     :: Oitest(3),OiInside(3)
integer(kind=DPI)     :: count
integer(kind=DPI)     :: iloop

real(kind=DP)                 :: randnumb,signL
real(kind=DP)                 :: Tausegloop(nsegloop), Tautotloop, tautmp, tauTLtmp
real(kind=DP),dimension(3,3)  :: sigtmp     ! The stress tensor at r
real(kind=DP),dimension(3)    :: r          ! Coordinate of the  stress calculation

logical               :: positive_loop,inside,outside
logical               :: keep

! The file with the nucleated loop info is re-open
loop_nucleated = "../out/loop_nucleated"
open(22,file=loop_nucleated,STATUS='OLD',POSITION='APPEND')

!####################################
print *, 'Enter periodic nucleation at step KK=',KK

! dimension of the nucleation loop
loop_edge = nint(loopcritR*tan(PII/real(NBASERED))/(avalue*1.e10))*IDEUX
loop_radius = nint(loopcritR/(avalue*1.e10)) ! "inscribed" radius

! definition of the stating coordinates for each loop we want to introduce
do ikk = 1, nloop!*nsegloop, nsegloop

  keep = .false.
  kunsuccess = -iun

  do while (.not. keep)

    ! In the previous calculation too many iteration as been made it useless to continue
    ! with this value for  loopcritR
    if (kunsuccess > kunsuccess_max) then
      print*,'Any position for nucleation is found, the critical radius must be increased !'
      print*,'loopcritR = ',loopcritR
      exit
    endif

    jj=1

    ! The loop slip system is randomly defined
    call random_number(randnumb)
    loop_slipsystem = INT(randnumb * 12 + 1,DPI)

    inside=.false.

    do while (.not. inside)

      call random_number(randnumb)
      Oitest(1)=domboundsX(1)+nint(randnumb*(domboundsX(2)-domboundsX(1)),DPI)

      call random_number(randnumb)
      Oitest(2)=domboundsY(1)+nint(randnumb*(domboundsY(2)-domboundsY(1)),DPI)

      call random_number(randnumb)
      Oitest(3)=domboundsZ(1)+nint(randnumb*(domboundsZ(2)-domboundsZ(1)),DPI)


      ! Oitest is moved to the closer simulation lattice point
      Oitest(1:3) =noeud(Oitest(1:3),IUN)



      call check_all_boundaries(real(Oitest,DP),outside,dom)
      if (outside) cycle


      do ii=1,nbregion

        OiInside(1:3)=Oitest(1:3)-box_position(1:3,ii)
        if (inoroutside(ii)) Then !Definition of an inclusive region
          if ((OiInside(1) > IZERO .and. OiInside(1) < h_boxx(ii)) .and. &
              (OiInside(2) > IZERO .and. OiInside(2) < h_boxy(ii)) .and. &
              (OiInside(3) > IZERO .and. OiInside(3) < h_boxz(ii)) ) then


            inside = .true.
            !print *, 'INSIDE','ciao'
            write (15,*) Oitest(1:3)
            exit
          endif
        else !Definition of an exclusive region
          if ((OiInside(1) > IZERO .and. OiInside(1) < h_boxx(ii)) .and. &
              (OiInside(2) > IZERO .and. OiInside(2) < h_boxy(ii)) .and. &
              (OiInside(3) > IZERO .and. OiInside(3) < h_boxz(ii)) ) then
            inside = .false.
           else
            inside = .true.
            exit
          endif
         endif

       enddo

    enddo

    loopSEG(1,2:4)=Oitest(1:3)
    call random_number(randnumb)
    count=IZERO

    do while(.not. keep .and. count < IDEUX)
      ! The loop sign of rotation is defined
      if (nint(randnumb,DP) == 0) then
        randnumb=1.d0
        positive_loop = .true.
        ! The first screw segment of the selected slip system
        loopSEG(1,1) = IUN + (loop_slipsystem - IUN) * nbasered
        loopSEG(1,5) = nint(loop_edge / NORMLIN(loopSEG(1,1)))

        !Define the 7 complementary segments to close the loop
        do jj = 2,nsegloop
          loopSEG(jj,1) = loopSEG(jj-1,1) + IUN
          loopSEG(jj,5) = nint(loop_edge/ NORMLIN(loopSEG(jj,1)))
          do ii =1,3
            loopSEG(jj,ii+1) = loopSEG(jj-1,ii+1) + loopSEG(jj-1,5) * bveclin(ii,loopSEG(jj-1,1))
            if (loopSEG(jj,ii+1) < 0) loopSEG(jj,ii+1) = loopSEG(jj,ii+1) + modur(ii)
          end do
          Oitest(1:3)=loopSEG(jj,2:4)
          call check_all_boundaries(real(Oitest,DP),outside,dom)
          if (outside) exit

        end do
        if (outside) then
          count=count+1
          cycle
        endif

      else
        randnumb=0.d0
        positive_loop=.false.
        loopSEG(1,1) = nbasered + (loop_slipsystem - IUN) * nbasered
        loopSEG(1,5) = nint(loop_edge / NORMLIN(loopSEG(1,1)))

        do jj = 2,nsegloop
          loopSEG(jj,1) = loopSEG(jj-1,1) - IUN
          loopSEG(jj,5) = nint(loop_edge / NORMLIN(loopSEG(jj,1)))
          do ii =1,3
            loopSEG(jj,ii+1) = loopSEG(jj-1,ii+1) + loopSEG(jj-1,5) * bveclin(ii,loopSEG(jj-1,1))
            if (loopSEG(jj,ii+1) < 0) loopSEG(jj,ii+1) = loopSEG(jj,ii+1) + modur(ii)
          end do
          Oitest(1:3)=loopSEG(jj,2:4)
          call check_all_boundaries(real(Oitest,DP),outside,dom)
          if (outside) exit
        end do
        if (outside) then
          count=count+1
          cycle
        endif

      end if

      ! define the full information for the virtual new segments building the tested loop
      do jj=1,nsegloop

        nsegm = nsegm + 1
        call Init_seg(nsegm)

        seg(nsegm)%O(:) = loopSEG(jj,2:4)
        seg(nsegm)%veclin = loopSEG(jj,1)
        seg(nsegm)%norme = loopSEG(jj,5)

        if (jj == 1) then
          seg(nsegm)%voiso = nsegm + nsegloop - IUN
          seg(nsegm)%vnno = nsegm + nsegloop - IUN
        else
          seg(nsegm)%voiso = nsegm - IUN
          seg(nsegm)%vnno = nsegm - IUN
        endif
        if (jj == nsegloop) then
          seg(nsegm)%voise = nsegm - nsegloop + IUN
          seg(nsegm)%vnne = nsegm - nsegloop + IUN
        else
          seg(nsegm)%voise = nsegm + IUN
          seg(nsegm)%vnne = nsegm + IUN
        endif
        if (nbcvxdom > IUN) then
          call assign_domain(nsegm,1414)
        else
          seg(nsegm)%dom = IUN
        endif

      end do
      keep=.true.
    enddo

    !check that the stress is sufficiently high for the tested loop expansion
    if (keep) then

      kunsuccess = kunsuccess + IUN
      Tausegloop(:) = zero
      Tautotloop = zero
      nloopplus = izero

      ! The stress calculation is made for each segment part of the loop
      do iloop = 1,nsegloop

        itemp = nsegm - nsegloop + iloop

        r = (seg(itemp)%O(:) + seg(seg(itemp)%vnne)%O(:)) * half

        ! We first calculate the internal stress + applied stress
        call sigma_int(r,sigtmp,itemp,tautmp)
        Tausegloop(iloop) = Tausegloop(iloop) + tautmp

        ! We add the line tension contribution
        call loop_linetension(itemp,Tautltmp,positive_loop)
        Tausegloop(iloop) = Tausegloop(iloop) + tauTLtmp

        Tautotloop = Tautotloop + Tausegloop(iloop)
        if (Tausegloop(iloop)>zero) nloopplus = nloopplus + iun  ! The number of segment with a positive force

      end do

      ! The loop is possibly in expansion if the average effective force * applied force is positive
      ! and if at least 6 segments move in the good direction
      if ((positive_loop .and. Tautotloop > zero .and. nloopplus > 6) .or. &
          (.not. positive_loop .and. Tautotloop < zero .and. nloopplus < 2) )  then

        print *, 'positive segments',nloopplus
        print *, 'Tautotloop',Tautotloop
        print *, 'positive_loop',positive_loop
        !if (kkdebug) then
          write(22,*) 'positive segments',nloopplus
          write(22,*) 'Tautotloop',Tautotloop
          write(22,*) 'positive_loop',positive_loop
        !endif
        ! A new nucleation loop is defined
        keep=.true.
        tot_loop = tot_loop + IUN

        if (positive_loop) then
          signL = 1.d0
        else
          signL =-1.d0
        endif

        write(22,*) 'Number loops generated = ',tot_loop
        write(22,*)'Number of nucleation iteration = ', kunsuccess
        write(22,*) 'Area swept',AireSys(loop_slipsystem),AireSysInst(loop_slipsystem)

77      format(I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,f12.6)
        do jj = nsegm-nsegloop+1 , nsegm

            daire(jj)=  0.5 * float(seg(jj)%norme*loop_radius) * NORMLIN(seg(jj)%veclin) * signL

            AireSys(loop_slipsystem)      = AireSys(loop_slipsystem) + daire(jj)
            AireSysInst(loop_slipsystem)  = AireSysInst(loop_slipsystem) + daire(jj)



            ! Calcul of the mechanical work per step associate to the applied stress = tauApp*gamma
            ! attension: dans les simulation la force = tau*b*longueur(en a)/(mu*a)
            !temp                = mudivV* daire(i)
            !TrAppSysInst(itemp) = TrAppSysInst(itemp) + tauApP(i)*temp

            ! Calcul of the mechanical work per step associate to the internal stress
            !TrIntSysInst(itemp) = TrIntSysInst(itemp) + (tauInt(I)+tauTL(I))*temp

            if (jj == 1 .or. jj == 5) then
              airevisSys(loop_slipsystem) = airevisSys(loop_slipsystem) + daire(jj)*sign(1.D0,SchmidSys(loop_slipsystem))
            elseif (jj == 3 .or. jj ==  7) then
              airecoinSys(loop_slipsystem) = airecoinSys(loop_slipsystem) + daire(jj)*sign(1.D0,SchmidSys(loop_slipsystem))
            else
              AireMixtSys(loop_slipsystem) = AireMixtSys(loop_slipsystem) + daire(jj)*sign(1.D0,SchmidSys(loop_slipsystem))
            endif

             write(22,77)         &
            jj                ,&
            seg(jj)%O(:)      ,&
            seg(jj)%norme     ,&
            seg(jj)%veclin    ,&
            seg(jj)%voiso     ,&
            seg(jj)%vnno      ,&
            seg(jj)%voise     ,&
            seg(jj)%vnne      ,&
            daire(jj)


        enddo
        write(22,*) 'Area swept nucl',AireSys(loop_slipsystem),AireSysInst(loop_slipsystem)
        call flush(379)
        exit

      else

         ! This loop does not pass the force criteria and must be eliminated
         keep=.false.
         do jj = nsegm-nsegloop+1 , nsegm
            seg(jj)%O(:)=izero
            seg(jj)%veclin=izero
            seg(jj)%norme=IZERO
            seg(jj)%voiso=izero
            seg(jj)%voise=izero
            seg(jj)%vnno=izero
            seg(jj)%vnne=izero
            seg(jj)%dom=IZERO
         enddo
         nsegm=nsegm-nsegloop

      end if
    end if

  enddo !end while loop

end do !end cycle over number of loops to be inserted

print *,'Number loops generated = ',tot_loop
print *,'Number of nucleation iteration = ', kunsuccess

if (kunsuccess > kunsuccess_max) then
  loopcritR = loopcritR + Delta_loopcritR
!   print*,'loopcritR= ', loopcritR
endif

close(22)

end subroutine periodic_nucleation


!##########################################################################
!# Calculation of the internal stress associate to a dislocation          #
!# microstructure with nsegm segments at point r                          #
!# Itemp is the index of segment attached to point r                      #
!##########################################################################
subroutine sigma_int(r,sigtmp,Itemp,Tautmp)

implicit none

real(kind=dp),dimension(6)               :: sigint
real(kind=dp),dimension(3,nsegmax)       :: ra0, rb0
real(kind=dp),dimension(3)               :: ra,rb
real(kind=DP), dimension(3)              :: Tbis,Bbis,BPVTbis

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



integer(kind=dpi)                        :: i,k,carac,NX,NY,NZ,Nrep(1:3),IndexB,RX,RY,RZ
integer(kind=dpi),dimension(3)           :: Tr,ab,centrei,BoxCenter_here

real(kind=DP),dimension(3,nsegmax)       :: T      ,&      ! vecteur ligne normalise
                                            B      ,&      ! vecteur de Burgers normalise
                                            BPVT           ! vecteur normale a B et T

real(kind=dp),dimension(3,3)             :: sigtmp                    ! The stress tensor at r
real(kind=dp),dimension(3)               :: r,G_tmp                   ! Coordinate of the  stress calculation
integer(kind=dpi)                        :: j,IVA,VLI,Itemp,compt
real(kind=dpi)                           :: TAUtmp,tmp,fortmp




! Initialization
BoxCenter_here(1:3) = modur(1:3)/2    !<The center of the simulation cell
Nrep(1:3) =0

compt     = 0
do K = 1,nsegm

  compt = compt + 1

  if (seg(k)%norme == 0) cycle      ! pivot segments can be eliminated from the calculations

  T(1:3,K)  = bveclin(1:3,seg(k)%veclin)/normlin(seg(k)%veclin)    !< The segment vector direction
  indexB    = (int((seg(k)%veclin-1)/8))*8 + 1
  B(1:3,K)  = bveclin(1:3,indexB)/normlin(indexB)                  !< The Burgers vector of i

  bpvt(1,k) = b(2,k)*t(3,k)-b(3,k)*t(2,k)         !< B vector product T
  bpvt(2,k) = b(3,k)*t(1,k)-b(1,k)*t(3,k)
  bpvt(3,k) = b(1,k)*t(2,k)-b(2,k)*t(1,k)

enddo


!initialization of NSDSF parameters
m4pn=0.25d0/(1.0d0-dpoiss)
mn4pn=m4pn*dpoiss
a2m8p=bspread2*0.125d0


sigint(:)     = 0.0d0
sigtmp(:,:)   = 0.0d0
sigapps(:,:)  = 0.0d0

! The translation needed for periodic boundaries
! With this translation r is shifted at the center of the simulation box
Tr(1:3)   = BoxCenter_here(1:3)- nint(r(1:3))
do i = 1,nsegm
  ! Those segments don t need to be considered
  if (seg(i)%norme == 0) cycle

  ab(:) = seg(i)%norme*bveclin(:,seg(i)%veclin)
  centrei(:) = seg(i)%o(:) + ab(:)/2
  centrei(:) = modulo(centrei(:)+Tr(:), ModuR(:))
  ra0(:,i) = real(BoxCenter_here(:) - (centrei(:) - ab(:)/2),DP)
  rb0(:,i) = real(BoxCenter_here(:) - (centrei(:) + ab(:)/2),DP)
enddo

! Loop on the periodic replicas
do NX = -Nrep(1),Nrep(1)
  do NY = -Nrep(2),Nrep(2)
    do NZ = -Nrep(3),Nrep(3)

      ! The segments loop
      do i = 1, nsegm

        if (seg(i)%norme == 0) cycle          ! We eliminate pivot segments
        if (i == itemp .or. i == seg(itemp)%voiso .or. i == seg(itemp)%voise) cycle


        ! information relative to the segments
        ! The segment type (the tyseg type 1 <-> 4!)
        carac = tyseg(seg(i)%veclin)
        ! The vectors connecting r to O and E of segment i
        ra(1) = ra0(1,i) + NX*MODUR(1)
        ra(2) = ra0(2,i) + NY*MODUR(2)
        ra(3) = ra0(3,i) + NZ*MODUR(3)
        rb(1) = rb0(1,i) + NX*MODUR(1)
        rb(2) = rb0(2,i) + NY*MODUR(2)
        rb(3) = rb0(3,i) + NZ*MODUR(3)

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
 !                  SIGINTERN(:)=SIGINTERN(:)+ BdivPA *(SIGBV(:)*DIVY2B-SIGAV(:)*DIVY2A)
          SIGINT(:)=SIGINT(:) + I_03*s_03
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

        SIGINT(:) = SIGINT(:) + (I_03*s_03 + I_13*s_13 + I_05*s_05 + I_15*s_15 + I_25*s_25)

      enddo !i = 1,nsegm

    enddo !NVOLZ
  enddo !NVOLY
enddo !NVOLX

SIGINT(:) = SIGINT(:) * bdivpa

! Different solutions of loading
if (mode_deformation==IHUIT) then

  ! If the applied stress is provided by an external calculation, e.g.FEM
  ! We simply define the value of SigappS at Ro(i)

  RX = int(real((modulo(int(R(1)),modur(1))-FEgrid_or(1))/(real(FEgrid_size(1))/real(dimFE_Sig(1),DP)),DP))+IUN
  RY = int(real((modulo(int(R(2)),modur(2))-FEgrid_or(2))/(real(FEgrid_size(2))/real(dimFE_Sig(2),DP)),DP))+IUN
  RZ = int(real((modulo(int(R(3)),modur(3))-FEgrid_or(3))/(real(FEgrid_size(3))/real(dimFE_Sig(3),DP)),DP))+IUN

  SigappS(1,1) = FE_SigappS(RX,RY,RZ,1)
  SigappS(2,2) = FE_SigappS(RX,RY,RZ,2)
  SigappS(3,3) = FE_SigappS(RX,RY,RZ,3)
  SigappS(2,3) = FE_SigappS(RX,RY,RZ,4)
  SigappS(1,3) = FE_SigappS(RX,RY,RZ,5)
  SigappS(1,2) = FE_SigappS(RX,RY,RZ,6)
  SigappS(3,2) = SigappS(2,3)
  SigappS(3,1) = SigappS(1,3)
  SigappS(2,1) = SigappS(1,2)

  ! SigappS is finally multiplied by stress_fact
  SigappS(:,:) = stress_fact * SigappS(:,:)

elseif ( key_crack .and. ( mode_deformation == IQUATRE .or. mode_deformation == ITROIS ) ) then

  ! simulations with the crack field

  call calculate_crackfield(int(R,DPI))

  SigappS(1,1) = Sigapp(1,1) + Sigcrack(1)
  SigappS(2,2) = Sigapp(2,2) + Sigcrack(2)
  SigappS(3,3) = Sigapp(3,3) + Sigcrack(3)

  SigappS(1,2) = Sigapp(1,2) + Sigcrack(4)
  SigappS(2,3) = Sigapp(2,3) + Sigcrack(5)
  SigappS(1,3) = Sigapp(1,3) + Sigcrack(6)

  SigappS(3,2) = SigappS(2,3)
  SigappS(3,1) = SigappS(1,3)
  SigappS(2,1) = SigappS(1,2)

else

  ! All the other cases of loading

  SigappS(1,1) = Sigapp(1,1)
  SigappS(2,2) = Sigapp(2,2)
  SigappS(3,3) = Sigapp(3,3)

  SigappS(1,2) = Sigapp(1,2)
  SigappS(2,3) = Sigapp(2,3)
  SigappS(1,3) = Sigapp(1,3)

  SigappS(3,2) = SigappS(2,3)
  SigappS(3,1) = SigappS(1,3)
  SigappS(2,1) = SigappS(1,2)

endif

! The effective stress calculation
sigtmp(2,1) = (sigint(6)+SigappS(2,1))!*XMU/1.d6
sigtmp(3,2) = (sigint(4)+SigappS(3,2))!*XMU/1.d6
sigtmp(3,1) = (sigint(5)+SigappS(3,1))!*XMU/1.d6
sigtmp(1,1) = (sigint(1)+SigappS(1,1))!*XMU/1.d6
sigtmp(2,2) = (sigint(2)+SigappS(2,2))!*XMU/1.d6
sigtmp(3,3) = (sigint(3)+SigappS(3,3))!*XMU/1.d6
sigtmp(1,2) = (sigint(6)+SigappS(1,2))!*XMU/1.d6
sigtmp(2,3) = (sigint(4)+SigappS(2,3))!*XMU/1.d6
sigtmp(1,3) = (sigint(5)+SigappS(1,3))!*XMU/1.d6


! The PK-force calculation
G_tmp(:)  = zero
TAUtmp    = zero
fortmp    = zero
tmp       = zero
VLI       = seg(Itemp)%veclin

do J = 1,3
  G_tmp(1) = G_tmp(1)  + SIGtmp(1,J)*B(J,Itemp)
  G_tmp(2) = G_tmp(2)  + SIGtmp(2,J)*B(J,Itemp)
  G_tmp(3) = G_tmp(3)  + SIGtmp(3,J)*B(J,Itemp)
enddo

! The resolved shear stress calculation
do IVA=1,3
  J       = IPERM(IVA)
  K       = IPERM(IVA+1)
  FORTMP  = G_tmp(J) * T(K,Itemp) - G_tmp(K)* T(J,Itemp)
  tmp     = VecNorDep(IVA,VLI) ! vecteur dep dans le sys system devie numero J
  TAUtmp  = TAUtmp + FORTMP * tmp
enddo

end subroutine sigma_int

!##########################################################################
!# Calculation of the line tension on the nucleation loop                 #
!##########################################################################
subroutine loop_linetension(itemp,Tautltmp,positive_loop)

integer(kind=DPI) :: I,itemp  !*** DIFFERENTS
real(kind=DP)  :: TauTLtmp
real(kind=DP)  :: alpha,AdivB
real(kind=DP)  :: RAYON2,ANGLE,TFACT,TSIGN,LogTerm,rayon
real(kind=DP), dimension(3) :: VCOURB!,VUMONT

integer(kind=DPI) :: VLI
logical           :: positive_loop

! Initialization
TauTLtmp  = zero
TFACT     = un / (4.D0*Pii*(UN-dpoiss))
AdivB     = un / BdivA

! The inscribed loop radius
!rayon=real(nint(loopcritR*cos(PII/real(nbasered))/(avalue*1.e10))) ! "inscribed" radius
rayon=real(int(loopcritR)/(avalue*1.e10),DP) ! "inscribed" radius
rayon2=rayon*rayon

i=itemp

VLI = -100         ! VLI Should never be used without be redefined in the loop
if (seg(i)%norme/=0) then

  VLI = seg(i)%veclin
  ! Calcul du vecteur rayon de courbure
  if (positive_loop) then
      VCOURB(1:3)=-UN*bvecdep(1:3,VLI)
  else
      VCOURB(1:3)=UN*bvecdep(1:3,VLI)
  endif
  VCOURB(1:3)=normavec(VCOURB(1:3))*dsqrt(rayon2)

  !alpha is the angle between PM and PQ, where P is the half point of the segment i,
  ! M and Q are the half point of the non null neighbours of the segment i
  !and PM and PQ are the segments linking these points
  alpha=PII/2.d0

  ! Angle=angle de la tangente locale a la direction vis

  !scalar product between the screw direction and the curvature radius
  ANGLE=VCOURB(1)*vecnorlin(1,assoc(VLi,1))+   &
        VCOURB(2)*vecnorlin(2,assoc(VLi,1))+   &
        VCOURB(3)*vecnorlin(3,assoc(VLi,1))
  !cosi ne of the angle between the normal to I and b
  ANGLE=ANGLE/dsqrt(rayon2)

  if (dabs(angle) > UN) angle = dsign(UN,angle)
  !angle de la tangente locale a la direction vis
  ANGLE       = dble(PII)      * HALF - dacos(ANGLE)

  !TSIGN est la projection de VCOURB  dans la direction de deplacement
  TSIGN=(VECnorDEP(1,VLI)*VCOURB(1)+ &
         VECnorDEP(2,VLI)*VCOURB(2)+ &
         VECnorDEP(3,VLI)*VCOURB(3))

  TSIGN=SIGN(1.D0,TSIGN)

  RAYON2=dsqrt(rayon2)

  !!!Foreman line tension
  LogTerm = DLOG(alpha*RAYON2*AdivB)
        if (abs(LogTerm) < un) then
          LogTerm = 1.
        endif
  TauTLtmp= TSIGN * TFACT * BdivA /(RAYON2)*                        &
            ((UN + dpoiss - trois * dpoiss * DSIN(ANGLE)**2)*       &
            LogTerm - (dpoiss*COS(2.0D0*ANGLE)))

end if

end subroutine loop_linetension

end module nucleation
