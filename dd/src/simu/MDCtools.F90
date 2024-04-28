
!===================================================================================================
!========================    DEBUT    MODULE   "MDCTOOLS"  =========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the Discrete Continuous Model (DCM) simulation.
!> Mainly the procedure defining the connectivity between mM et the FE simulation code ZeBuLon.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module MDCTOOLS

use varglob
use bricamat
use constantes
!use DEBUG

implicit none

real (kind = DP)                          :: deltat_zeb
real (kind = DP),dimension(6,6)           :: kelastique
real (kind = DP)                          :: xmu_mpa
real (kind = DP)                          :: mdcvol
integer                                   :: endcalcul,restartstep

contains

!###############################################################################
!# Send to Zebulon some global data                                            #
!###############################################################################
subroutine send_globdata

  implicit none

  call zebmpipack(real(ModuR(1:3), DP), 3)
  call zebmpipacksend

  call zebmpipack(real(vecburgers, DP))
  call zebmpipack(real(halfthickness, DP))
  call zebmpipack(real(avalue, DP))
  call zebmpipacksend

end subroutine send_globdata

!###############################################################################
!# Get from Zebulon some global data                                           #
!###############################################################################
subroutine recv_globdata

  implicit none

  call zebmpirecv(mdcvol)

end subroutine recv_globdata

!################################################################
!# The restart switch:                                          #
!# restart_key = 0 --> new simulation                           #
!# restart_key = 1 --> restart and old simulation               #
!#                                                              #
!################################################################
subroutine restart

  implicit none
  integer :: restart_key

  !-MPI- Begin Z ->DD No. 4
  call zebmpirecv(restart_key)
  sideja = restart_key

  call zebmpisend(444)

end subroutine restart

!------------------------------------------------------------------------!
!  Micromegas receives its time step from Zebulon                        !
!------------------------------------------------------------------------!
subroutine time_step

  implicit none
  
  integer (kind = DPI) :: i,discr_param

   

deltat=zero

#if defined(MDC) && defined(PA)
  if (Mon_Rang == IZERO) then
#endif
  !-MPI- Begin Z ->DD No. 7
  !-MPI- Read from problem_static_mechanical_dislocation::make_increment
  call zebmpirecv(mdc_timestep_ratio)

  call zebmpirecv(deltat)

  call zebmpisend(777)
  !-MPI- End   Z ->DD No. 7

  if (kk==0 .or. kkdebug) write(379,*) "mM received deltat = ", deltat
  if (kk==0 .or. kkdebug) write(379,*) "Zebulon deltat = ", deltat*mdc_timestep_ratio

#if defined(MDC) && defined(PA)
  endif

CALL MPI_BCAST(deltat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(mdc_timestep_ratio,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

#endif

KRC = mdc_timestep_ratio

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
if (Mon_Rang ==  IZERO) then
#endif

nsegm_INI = nsegm

discr_param = IZERO
do i = 1,nsegm
discr_param = discr_param + seg(i)%norme
enddo

discr_param = discr_param/(2*loma(1))/nsegm



write(379,*) 'intial size of the sweptsurfdata array'
write(379,*) 'nsegm = ', nsegm
write(379,*) 'discretization parameter',discr_param  
write(379,*) 'DCM time step ratio', mdc_timestep_ratio
write(379,*) 'size of the sweptsurfdata array ', nsegm_INI*discr_param*IQUATRE*mdc_timestep_ratio*16+1
call flush()
! allocate sweptsurfdata array
allocate(sweptsurfdata(nsegm_INI*discr_param*IQUATRE*mdc_timestep_ratio*16+1))
nsegm_INI = nsegm_INI*IQUATRE*discr_param
sweptsurfdata(:)=0
#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
endif
#endif

end subroutine time_step

!----------------------------------------------------------------------!
! Zebulon tells to Micromegas what to do...
! endcalcul=0 --> new time step
! endcalcul=1 --> exit -- the simulation is finished
!----------------------------------------------------------------------!
subroutine whattodo

  implicit none

#if defined(MDC) && defined(PA)

  if (Mon_Rang == 0) then
        ! Only the process zero of the communicator can write in output files, the others are waiting
#endif

  !-MPI- Begin Z <-DD No. 9
  call zebmpirecv(endcalcul)

  !call zebmpirecv(restartstep)

  call zebmpisend(999)
  !-MPI- End   Z <-DD No. 9

#if defined(PA) && defined(MDC)
    ! Le proc zero a fini d ecrire, on libere tout le monde

  endif
  CALL MPI_BCAST(endcalcul,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

#endif

  if( endcalcul>0 ) then
    ! On va fermer les fichiers de stockage des differentes quantites

   close( 8)
   close(29)

    print *
    print *,'===============  ZEBULON TELLS MM TO EXIT  ==============='
    print *
  endif

end subroutine whattodo


subroutine dist2ps(p, c, t, l2,  buf, dn2, dt2, sdt2)
  ! IN  p: point considered
  ! IN  c: center of segment
  ! IN  t: segment vector ("end point" - "begin point")
  ! IN  l2: square of the norm ofthe segment
  ! IN  buf: buffer
  ! OUT dn2: square of the normal distance to the segment
  ! OUT dt2: square of the distance of the projection of the point on the segment line to c
  ! OUT sdt2: sign of dt2

  implicit none

  real (kind=DP), dimension(3)      :: p, c, t
  real (kind=DP), dimension(11)     :: buf
  real (kind=DP)                    :: l2, dn2, dt2, sdt2

  buf(1:3)  =p(1:3)-c(1:3) ! v
  buf(4)    =buf(1)*t(1)+buf(2)*t(2)+buf(3)*t(3) ! v.t
  buf(5)    =buf(4)/l2 ! (v.t)/(t.t)
  buf(6:8)  =buf(5)*t(1:3) ! vt=(v.t)/(t.t)*t
  buf(9:11) =buf(1:3)-buf(6:8) ! vn=v-vt

  dn2=buf(9)*buf(9)+buf(10)*buf(10)+buf(11)*buf(11)
  dt2=buf(6)*buf(6)+buf(7)*buf(7)+buf(8)*buf(8)
  sdt2=sign(1.D0, buf(4))

end subroutine dist2ps
! ENDTAG OJ



!#########################################################
!# The intial plastic deformation corresponding to the   #
!# Voltera loop process creation is calculated           #
!#########################################################
subroutine init_plas_eigen

  implicit none


  integer,dimension(3)        :: or1,or2,or3,or4
  integer,dimension(3)        :: disp,burger,loop_center
  integer                     :: idistance,eigensign
  integer(kind=DPI)                     :: i,k


  integer, parameter :: nbmaxinitseg=20000

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
if (Mon_Rang ==  IZERO) then
#endif

  allocate(sweptsurfdata(1+nbmaxinitseg*16))

  nbswept=0
  sweptsurfdata(1)=1



  do i = 1,nsegm

    ! Only the first segment of dipolar loop is needed
    if (MOD((i-1),iquatre) /= IZERO .or. seg(i)%voiso == seg(i)%voise) cycle



    ! The first edge segment of the dipolar loop is indexed as 'k'
    k = i
    if ((GB/=0 .and. seg(k)%vnno/= k+3) .or. (GB==0  .and. seg(k)%vnno/= -3)) then
      print *,"the DCM initial configuration supports only loops made by 4 segments"
      print *,seg(k)%vnno, k+3
    endif


    or1=int(seg(k)%O)
    or2=or1 + int(seg(k)%norme*bveclin(:,seg(k)%veclin))
    or3=or2 + int(seg(k+1)%norme*bveclin(:,seg(k+1)%veclin))
    or4=or3 + int(seg(k+2)%norme*bveclin(:,seg(k+2)%veclin))


    !sign of the plastic eigenstrain in the initial loop


    if (syseg(seg(k)%veclin) == syseg(seg(k+1)%veclin)) then

       loop_center = (or1 + or2 + or3 + or4)/4

       if (dot_product(bvecdep(1:3,seg(k)%veclin),(loop_center-or1)) < 0) then
          eigensign = 1
       else
         eigensign = -1
       endif

    else

      if (dot_product(bvecdep(1:3,seg(k)%veclin),&
         bvecdep(1:3,seg(k+1)%veclin)) > 0) then
          eigensign = 1
      else
         eigensign = -1
      endif


    endif

    ! The number of virtual displacement we must impose to edge 1 to build-up the Voltera loop
    ! Equal the norme of edge 2
    idistance = int(seg(k+1)%norme)



    ! The edge 1 displacement direction = line direction of edge 2
    disp(1:3) = int(bveclin(1:3,seg(k+1)%veclin))

    ! The dipolar loop Burgers vector
    burger(1:3) = int(bveclin(1:3,assoc(seg(k)%veclin,1)))

    sweptsurfdata(1+16*nbswept+1)=idistance*eigensign
    sweptsurfdata(1+16*nbswept+2:1+16*nbswept+4)=or4(1:3)
    sweptsurfdata(1+16*nbswept+5:1+16*nbswept+7)=or3(1:3)
    !sweptsurfdata(1+16*nbswept+8:1+16*nbswept+10)=-1*disp(1:3)
    !sweptsurfdata(1+16*nbswept+11:1+16*nbswept+13)=-1*disp(1:3)
    sweptsurfdata(1+16*nbswept+8:1+16*nbswept+10)=-1*disp(1:3)*idistance*eigensign
    sweptsurfdata(1+16*nbswept+11:1+16*nbswept+13)=-1*disp(1:3)*idistance*eigensign
    sweptsurfdata(1+16*nbswept+14:1+16*nbswept+16)=burger(1:3)


    nbswept=nbswept+1

  enddo           ! End the dislocation lines loop

  sweptsurfdata(1)=nbswept

  call zebmpisend(sweptsurfdata, 1+16*nbswept)

  nbswept = 0

  ! print *,'End of subroutine defplasconf'

  deallocate(sweptsurfdata)

#if defined(MDC) && defined(PA)
!Only Process 0 is in charge of transferring data to Zebulon
endif
#endif

end subroutine init_plas_eigen




!###########################################################
!# The tensor of elasticity is defined at each gauss point #
!# Z sends its components in MPa (as they are specified    #
!# in the Z material inputfile), and they are immediately  #
!# converted to the mM Pa)                                 #
!###########################################################
subroutine elasticity_matrix

integer      :: m
real(kind=8) :: mm,lambda

!-MPI- Begin Z ->DD No. 3
!-MPI- Talking to mcesd_stdi::matrix_of_elasticity_for_dislocation()
!-MPI- in mcesd_std_dislocation.c

do m = 1,6
  call zebmpipackrecv
  call zebmpiunpack(kelastique(m,1:6), 6)
  kelastique(m,1:6) = kelastique(m,1:6)*1.0d6
enddo

call zebmpisend(333)
!-MPI- End Z ->DD No. 3

lambda = kelastique(1,2)
mm     = kelastique(4,4)*0.5d0 !!k44 =2c44
dpoiss = 0.5d0*lambda/(lambda+mm)
write(66,*) "elastic constant [GPa]"
write(66,*)  kelastique(1,1:6)*1e-9
write(66,*)  kelastique(2,1:6)*1e-9
write(66,*)  kelastique(3,1:6)*1e-9
write(66,*)  kelastique(4,1:6)*1e-9
write(66,*)  kelastique(5,1:6)*1e-9
write(66,*)  kelastique(6,1:6)*1e-9

xmu = mm

write(66,*) 'Material parameter: '
write(66,*) 'shear modulus [GPa]: ',xmu * 1.0d-9
write(66,*) 'Poisson ration=',dpoiss

end subroutine elasticity_matrix


!########################################################################
!# The Reuss average is used to define the isotropic shear modulus      #
!# Called from subroutine parametrage (in 07init.f90, MDC only)         #
!########################################################################
subroutine FEshear_modulus

implicit none
double precision :: k55, k11, k12

k55 = kelastique(5,5)
k11 = kelastique(1,1)
k12 = kelastique(1,2)

! The average is defined as : gr=sqrt(c44(c11-c12)/2)
! For Zebulon, k55 = 2*c44, so gr=1/2*sqrt(k55(k11-k12))
xmu     = half*dsqrt( k55*(k11-k12) ) ! xmu en  Pa
xmu_mpa = xmu*1.0d-6                  ! xmu en MPa

write(66,*) 'Material parameter: '
write(66,*) 'Reuss average shear modulus [MPa]: ',xmu_mpa

end subroutine FEshear_modulus



end module mdctools
