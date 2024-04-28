!This program generates a list of vertices for each domain defined in b_plan_conc file.
!Once the vertices are computed a python program is used to compute the volume of
!the simulated geometrical object

module INITCONFTOOLS

use VARGLOB
use CONSTANTES
use VARBASE
use BRICAMAT

implicit none

integer,parameter  :: nsmax=20                    ! Maximum Number of systems
integer,parameter  :: nbmax=160                   ! Maximum Number of elementary character (nbasered x nsmax)
integer,parameter  :: ndismax=50000               ! Maximum number of segments

real,parameter     :: pisurdeux=1.570796327       ! Pi/2
real,parameter     :: pi=3.141592654              ! Pi
real,parameter     :: troispisurdeux=4.712388980  ! 3Pi/2

type(plans)  :: PlaneDef(100)

!new definition of planes!! derived type was slow better to define arrays...
integer(kind=DPI), allocatable :: PlaneDef_MillerI(:,:)
real(kind=DP), allocatable     :: PlaneDef_MillerR(:,:)
real(kind=DP), allocatable     :: PlaneDef_pos(:)
integer(kind=DPI), allocatable :: PlaneDef_knd(:)
integer(kind=DPI), allocatable :: PlaneDef_dom(:)

! These are unsigned integers in the C version
INTEGER, SAVE :: s1 = 123, s2 = 456, s3 = 789


contains


!###########################################################################################
!< This subroutine generates a list of vertices for each domain defined in b_plan_conc file.
!< Once the vertices are computed a python program is used to compute the volume of
!< the simulated geometrical object
!###########################################################################################
subroutine vertex_calculations

integer(kind=DPI)               :: i,ii,l,n,jj
integer(kind=DPI)               :: fixplaneid,freeplneid,ncomb,r
integer(kind=DPI)               :: plns(3),tmpvtx
integer(kind=DPI), allocatable  :: planes(:),list(:,:)

real(kind=DP), allocatable      :: vtx(:,:,:)
real(kind=DP)                   :: xminvtx(nbcvxdom),xmaxvtx(nbcvxdom)
real(kind=DP)                   :: yminvtx(nbcvxdom),ymaxvtx(nbcvxdom)
real(kind=DP)                   :: zminvtx(nbcvxdom),zmaxvtx(nbcvxdom)
real(kind=DP)                   :: coord(3)

logical                         :: solution,TestGood

open(38,file='../out/vertexlist',status='UNKNOWN')

!number of planes we need to intersect to compute a vertex (used to compute combinations without repetitions)
r = 3

do i=1,Nbplanmax+6
  write(*,fmt='(I8,I8,I8,F12.3)')     PlaneDef_MillerI(1:3,i),PlaneDef_pos(i)
  write(*,fmt='(F12.3,F12.3,F12.3)')  PlaneDef_MillerR(1:3,i)
enddo

fixplaneid=0
freeplneid=Nbplan

allocate(vtx(nbcvxdom,200,3))

write(38,*) Nbcvxdom, '  0  ','  0  '

do i=1,Nbcvxdom

  n=NbPlanCvxDom(1,i)+NbPlanCvxDom(2,i)+6

  !number of possible plane combinations
  ncomb=Factorial(n)/(Factorial(r)*Factorial(n-r))

  allocate(planes(n))
  allocate(list(3,ncomb))

  !print *,'ncomb',ncomb,n,r
  l=0
  do ii=fixplaneid+1,fixplaneid+NbPlanCvxDom(1,i)
     l=l+1
     planes(l)=ii
  enddo

  do ii=freeplneid+1,freeplneid+NbPlanCvxDom(2,i)
     l=l+1
     planes(l)=ii
  enddo

  do ii=nbplanmax+1,nbplanmax+6
      l=l+1
      planes(l)=ii
  enddo

  call combinations(ncomb,n,planes,list)

  tmpvtx=0

A1:do ii=1,ncomb

    plns=list(1:3,ii)

    call CramerRule(plns,solution,coord)

    if (solution) then

      Testgood=.true.
      call check_boundaries(i,coord,TestGood)

      if (.not. Testgood) then
        if (tmpvtx /=  0) then
          do jj=1,tmpvtx
            if(abs(vtx(i,jj,1)-coord(1)) < 1.d-5 .and. &
               abs(vtx(i,jj,2)-coord(2)) < 1.d-5 .and. &
               abs(vtx(i,jj,3)-coord(3)) < 1.d-5) cycle A1
          enddo

        endif
        tmpvtx=tmpvtx+1
        if (tmpvtx > 200) then
          write(*,*) 'Strange, Too many vertices!!(> 200)'
          stop
          endif
        vtx(i,tmpvtx,1:3)=coord

      endif
    endif
  enddo A1

  !now we can calculate the intersection using the Cramer formula..
  print *,'number of vertex domain ', i,' = ', tmpvtx
  write(38,*) tmpvtx,'  0  ','  0  '
  do ii=1,tmpvtx
     write(38,fmt='(F20.3,F20.3,F20.3)') vtx(i,ii,1:3)
  enddo

  fixplaneid=fixplaneid+NbPlanCvxDom(1,i)
  freeplneid=freeplneid+NbPlanCvxDom(2,i)

  xminvtx(i)=minval(vtx(i,1:tmpvtx,1))
  xmaxvtx(i)=maxval(vtx(i,1:tmpvtx,1))
  yminvtx(i)=minval(vtx(i,1:tmpvtx,2))
  ymaxvtx(i)=maxval(vtx(i,1:tmpvtx,2))
  zminvtx(i)=minval(vtx(i,1:tmpvtx,3))
  zmaxvtx(i)=maxval(vtx(i,1:tmpvtx,3))

  deallocate(planes)
  deallocate(list)

enddo  !Nbcvxdom

domboundsX(1)=nint(minval(xminvtx(:)),DPI)
domboundsX(2)=nint(maxval(xmaxvtx(:)),DPI)
domboundsY(1)=nint(minval(yminvtx(:)),DPI)
domboundsY(2)=nint(maxval(ymaxvtx(:)),DPI)
domboundsZ(1)=nint(minval(zminvtx(:)),DPI)
domboundsZ(2)=nint(maxval(zmaxvtx(:)),DPI)

close(38)

end subroutine vertex_calculations


!#############################################
!<
!#############################################
RECURSIVE function Factorial(n)  RESULT(Fact)

implicit none
Integer(kind=DPI) :: Fact
integer(kind=DPI), intent(IN) :: n

if (n == 0) then
   Fact = 1
else
   Fact = n * Factorial(n-1)
endif

end function Factorial


!############################################
!<
!############################################
subroutine combinations(ncomb,n,planes,list)

implicit none

integer(kind=DPI), intent(IN) :: ncomb,n,planes(n)
integer(kind=DPI), intent(OUT) :: list(3,ncomb)
integer(kind=DPI) :: i,j,k,ccomb
integer(kind=DPI) :: loc1,loc2,loc3

ccomb = 0
do i =1,n
  loc1 = planes(i)
  do j =1,n
    if (j .gt. i) then
      loc2 = planes(j)
      do k=1,n
        if (k .gt. j) then
          loc3 = planes(k)
          ccomb=ccomb+1
          list(1:3,ccomb)=(/loc1,loc2,loc3/)
        endif
      enddo !k
    endif
  enddo !j
enddo !i

end subroutine combinations


!##########################################
!<
!##########################################
subroutine check_boundaries(ii,O,TestGood)

implicit none

logical :: Testgood,bloquer

real(kind=DP)         :: O(3),TestO,normint(3)
integer(kind=DPI)     :: JJ,ii
integer(kind=DPI)     :: plane,planes

plane=0
PLANES=Nbplan

if (ii > 1) then
   do jj=1,ii-1
     plane= plane+NbPlanCvxDom(1,jj)
     planes= planes+NbPlanCvxDom(2,jj)
   enddo
endif

bloquer=.false.

do jj= plane+1,plane+NbPlanCvxDom(1,ii)

  ! Normale projection of OiP on the barriere
  !inormint(:)= PlaneDef(jj)%miller(:)
  !normint(:) = real(inormint(:),DP)
  normint(1:3) = PlaneDef_MillerR(1:3,jj)
  TestO=O(1)*normint(1)+O(2)*normint(2)+O(3)*normint(3)-PlaneDef_pos(jj)

  ! if (kkdebug) then
!       write (379,*) ""
!       write (379, fmt= '(">>>>> Detecting FIXED BARRIER ")')
!       write (379,*) ""
!       write (379,fmt='(5A7)') "VFP","kind", "","Miller",""
!       write (379,fmt='(5I7)') jj,PlaneDef(jj)%knd,PlaneDef(jj)%miller(:)
!       write (379,*) ""
!       write (379,fmt='("Coord",F20.3)'),O
!       write (379,fmt='("SProd (Vect1O,Miller)",F20.3)') TestO
!       write (379,*) ""
!     endif

  if (TestO > 1.e-6) then
    bloquer=.true.
    !write(379,*) "1 outside domain ******", ii,jj,testo
    exit
  endif

enddo

if (.not. bloquer)  then

  do jj =Planes+1,Planes+NbPlanCvxDom(2,ii)
     !inormint(:)= PlaneDef(jj)%miller(:)
     !normint(:) = real(inormint(:),DP)
     normint(1:3) = PlaneDef_MillerR(1:3,jj)
     TestO=O(1)*normint(1)+O(2)*normint(2)+O(3)*normint(3)-PlaneDef_pos(jj)
!        if (kkdebug) then
!          write (379,*) ""
!          write (379, fmt='(">>>>> Detecting  Free BARRIER ")')
!          write (379,*) ""
!          write (379,fmt='(5A7)') "VFP","kind", "","Miller",""
!          write (379,fmt='(5I7)') jj,PlaneDef(jj)%knd,PlaneDef(jj)%miller(:)
!          write (379,*) ""
!          write (379,fmt='("Coord",F20.3)'), O
!          write (379,fmt='("SProd (Vect1O,Miller)",F20.3)') TestO
!          write (379,*) ""
!       endif

    if (TestO > 1.e-6)  then
      bloquer=.true.
      !write(379,*) "2: outside domain *****", ii, jj ,testO
      exit
    endif

  enddo   ! End of loop jj

endif

if(bloquer) then
  Testgood = .true.
else
  TestGood = .false.
endif

end subroutine

!###############################
!<
!###############################
subroutine load_boundaries

implicit none

integer(kind=DPI)     :: I,J,JJ,II
integer(kind=DPI)     :: plane,planeFP
integer(kind=DPI)     :: itmp,jtmp

integer(kind=DPI),allocatable :: TMiller(:,:),Tknd(:)
real(kind=DP),allocatable     :: Tpos(:)

character(len=1)              :: carac

! The data file
open(44,file="../in/bond_domains_def",STATUS='OLD')

! Comments in the file header are removed
 carac = 'x'
do while (carac /= "#")
    read (44,*) carac
enddo

! information about the plane configuration
read(44,*) NbCvxDom     ! number of convex domains
read(44,*) key_Rot_Dom  ! Key to activate the crystal rotation in domains

NbPlanDom = IZERO

allocate(ListcvxDom(NbCvxDom))
allocate(ListcvxDomsegi(NbCvxDom))
allocate(DomConnec(NbCvxDom,NbCvxDom))

if (NbCvxDom /= IZERO) then
  allocate (NbPlanCvxDom(2,NbCvxDom))
  !print *,NBCVXDOM
  do I=1,NbCvxDom
    ListcvxDom(i) = i
    read(44,*) NbPlanCvxDom(1,i)
    read(44,*) NbPlanCVxDom(2,i)

    NbPlanDom=NbPlanDom+NbPlanCvxDom(1,i)+NbPlanCvxDom(2,i)
  enddo
else
  print *,"ERROR: at least one close domain must be defined!!!!"
  stop
endif

NbPlanMax = NbPlanDom
!print *,NbPlanMax
! The list of plane equations

allocate (TMiller(NbplanMax,3))
allocate (TPos(Nbplanmax))
allocate (TKnd(Nbplanmax))

allocate (Plane_MillerI(3,NbplanMax))
allocate (Plane_MillerR(3,NbplanMax))
allocate (Plane_pos(Nbplanmax))
allocate (Plane_knd(Nbplanmax))
allocate (Plane_dom(NbplanMax))

do I=1 ,NbPlanMax
  read(44,*) TMiller(i,:), Tpos(i), Tknd(i)
  !Print *, TMiller(i,:) , Tpos(i)
enddo

JJ=IZERO
NbPlan=IZERO

do I=1 ,NbCvxDom
  !Print *, NbCvxDom
  !print *, NbPlanCvxDom(1,i)
  do J=1,NbPlanCvxDom(1,i)
     NbPlan=NbPlan + IUN
     Plane_MillerI(1:3,Nbplan)=TMiller(JJ+NbPlan,:)
     Plane_MillerR(1:3,Nbplan)=real(TMiller(JJ+NbPlan,:),DP)
     Plane_pos(NbPlan)=Tpos(JJ+NbPlan)
     Plane_knd(NbPlan)=TKnd(JJ+NbPlan)
     Plane_dom(NbPlan)=i
  enddo
  JJ=JJ+NbPlanCvxDom(2,i)
enddo

NbFreeplan = IZERO
JJ= IZERO

do I=1 ,NbCvxDom
  JJ=JJ+NbPlanCvxDom(1,i)
  !Print *, NbCvxDom
  !print *, NbPlanCvxDom(2,i)
  do J=1,NbPlanCvxDom(2,i)
     NbFreePlan=NbFreePlan + IUN
     !Print *,JJ+NbFreePlan , NbFreeplan+NbPlan
     !read(*,*)
     Plane_MillerI(1:3,NbFreeplan+NbPlan)=TMiller(JJ+NbFreePlan,:)
     Plane_MillerR(1:3,NbFreeplan+NbPlan)=real(TMiller(JJ+NbFreePlan,:),DP)
     Plane_pos(NbFreeplan+NbPlan)=Tpos(JJ+NbFreePlan)
     Plane_knd(NbFreeplan+NbPlan)=TKnd(JJ+NbFreePlan)
     Plane_dom(NbFreeplan+NbPlan)=i
  enddo
  !JJ=JJ+NbPlanCvxDom(1,i)
  !print *, 'jj end',jj
enddo
!Print *,NbPlan
!Print *,NbFreePlan

!Connectivity matrix for domain interaction (Internal Stress only)
!True if domains should not see each other for internal stress calculation
DomConnec(:,:) = .True.

if (NbCvxDom > IUN) then
  do I=1 ,((NbCvxDom*NbCvxdom)-nbcvxdom)/2
    read(44,*) itmp,jtmp,DomConnec(itmp,jtmp)
  enddo
endif

do I=1 ,NbCvxDom
  do J=1 ,NbCvxDom
    if (I==J) DomConnec(I,I) = .False. ! Domain see itself
    if (.not. DomConnec(I,J)) DomConnec(J,I)=.False. !Symetrize
  enddo
enddo

write(*,*)
write(*,fmt='("!> Connectivity matrix for domain interaction (Internal Stress only)")')
write(*,fmt='("   X",50I4)') ( I, I=1,NbCvxDom)
do I=1 ,NbCvxDom
  write(*,fmt='(I4,50L4)') I,DomConnec(I,:)
enddo
write(*,*)

open(335,file='../out/debug/surflist.txt',STATUS='unknown')
  write(335,*) "------------------------------------"
  write(335,*) "Loaded boundaries"
  write(335,*) "------------------------------------"
  write (335,fmt='(A3,A5,3A20,A15,A7,A7,A7)') &
  "#","VFP","","Miller","","VPpos","Dom","Kind"
  write (335,fmt='("")')

    plane=0
    planeFP=0

  do ii=1,NbCvxDom
    write (335,fmt='("Domain : ", I4)') ii
    do jj= plane+1,plane+NbPlanCvxDom(1,ii) ! Loop over transparent and grain boundary surfaces
      write (335,fmt='(A3,I5,3I20,F15.3,I7,A6)') &
        "#" ,jj, Plane_MillerI(1:3,jj), Plane_pos(jj), Plane_knd(jj),""
    enddo

    if (NbPlanCvxDom(2,ii)/=IZERO) then ! Check if there is free surfaces for domain ii
      do jj= NbPlan+planeFP+1,NbPlan+planeFP+NbPlanCvxDom(2,ii) ! Loop over free surfaces
        write (335,fmt='(A3,I5,3I20,F15.3,I7,A6)') &
          "#" ,jj,Plane_MillerI(1:3,jj), Plane_pos(jj), Plane_knd(jj),"FS"
      enddo
    endif

    plane=plane+NbPlanCvxDom(1,ii)     !Update plane for domain ii+1
    planeFP=planeFP+NbPlanCvxDom(2,ii) !Update planeFP for domain ii+1
    write (335,fmt='("")')
  enddo
close(335)

close (44)

allocate (PlaneDef_MillerI(3,NbplanMax+ISIX))
allocate (PlaneDef_MillerR(3,NbplanMax+ISIX))
allocate (PlaneDef_pos(Nbplanmax+ISIX))
allocate (PlaneDef_knd(Nbplanmax+ISIX))
allocate (PlaneDef_dom(NbplanMax+ISIX))

PlaneDef_MillerI(1:3,1:NbplanMax)=Plane_MillerI(1:3,1:NbplanMax)
PlaneDef_MillerR(1:3,1:NbplanMax)=Plane_MillerR(1:3,1:NbplanMax)
PlaneDef_pos(1:Nbplanmax)=Plane_pos(1:NbplanMax)
PlaneDef_knd(1:Nbplanmax)=Plane_knd(1:NbplanMax)
PlaneDef_dom(1:Nbplanmax)=Plane_dom(1:NbplanMax)

do i=1,6
  select case (i)
    case(1,4)
      PlaneDef_MillerI(1:3,Nbplanmax+i)=(/1,0,0/)
      PlaneDef_MillerR(1:3,Nbplanmax+i)=real(PlaneDef_MillerI(1:3,Nbplanmax+i),DP)
    case(2,5)
      PlaneDef_MillerI(1:3,Nbplanmax+i)=(/0,1,0/)
      PlaneDef_MillerR(1:3,Nbplanmax+i)=real(PlaneDef_MillerI(1:3,Nbplanmax+i),DP)
    case(3,6)
      PlaneDef_MillerI(1:3,Nbplanmax+i)=(/0,0,1/)
      PlaneDef_MillerR(1:3,Nbplanmax+i)=real(PlaneDef_MillerI(1:3,Nbplanmax+i),DP)
  end select
  if (i<=3) then
    PlaneDef_pos(Nbplanmax+i)=real(modur(i),DP)
  else
    PlaneDef_pos(Nbplanmax+i)=zero
  endif

  PlaneDef_knd(Nbplanmax+i)=IZERO
  PlaneDef_dom(Nbplanmax+i)=IZERO

enddo

end subroutine



!################################################################
!< Cramer's rule for 3x3 matrices
!< to compute intersection between planes define in b_plan_conc
!<       2 x -3 y +4 z = 7
!<       5 x -6 y +8 z = 9
!<       11 x +13 y +14 z = 15
!#################################################################
subroutine CramerRule(plns,solution,coord)

implicit none

integer(kind = DPI), Dimension (3), intent(IN) :: plns
logical, intent(OUT) :: solution
Real(kind = DP), Dimension (3), intent(OUT) ::  coord
Real(kind = DP), Dimension (3) :: A, B, C ,D
Real(kind = DP) :: sysdet, xdet, ydet, zdet

solution =.false.
D=(/PlaneDef_pos(plns(1)),PlaneDef_pos(plns(2)),PlaneDef_pos(plns(3))/)

A = (/PlaneDef_millerR(1,plns(1)),PlaneDef_millerR(1,plns(2)),PlaneDef_millerR(1,plns(3))/)
B = (/PlaneDef_millerR(2,plns(1)),PlaneDef_millerR(2,plns(2)),PlaneDef_millerR(2,plns(3))/)
C = (/PlaneDef_millerR(3,plns(1)),PlaneDef_millerR(3,plns(2)),PlaneDef_millerR(3,plns(3))/)

sysdet = Determinant(A, B, C)

if (Abs(sysdet) < 10e-15) Then

    solution = .false.
    coord = (/-300000,-300000,-300000/)

else

  solution=.true.

  xdet = Determinant(D, B, C)
  ydet = Determinant(A, D, C)
  zdet = Determinant(A, B, D)

  coord(1)=xdet/sysdet
  coord(2)=ydet/sysdet
  coord(3)=zdet/sysdet

endif

Contains

Function Determinant(A, B, C)
  Real(kind=DP), Dimension (3), Intent (in) :: A, B, C
  Real(kind=DP) :: Determinant
  Determinant = - A(3)*B(2)*C(1) + A(2)*B(3)*C(1) + A(3)*B(1)*C(2) &
                - A(1)*B(3)*C(2) - A(2)*B(1)*C(3) + A(1)*B(2)*C(3)
End Function Determinant

End subroutine CramerRule


!###################################################
!< Calculation of the segment initial distribution #
!###################################################
subroutine microconf

implicit none

integer(kind=DPI) :: tnl,ix,iy,iz,acu,indb(3),oindb(3),nsys,NBinome
integer(kind=DPI) :: i,j,k,syst(nsmax),iv(3),it(3),nbud,Oi(3),Ei(3),OEi(3),segini(7,ndismax)
integer(kind=DPI) :: nsegm,nmax(3),longv(12),longc(12),pmil,pmilmun,nbtrue,incr,ivtemp(3)
integer(kind=DPI) :: nbtruels,stat(11),nbsfrparsys
integer(kind=DPI) :: UnCcoin,DeuxCcoin,UnCvis,DeuxCvis,DeuxBinome,poidsys
integer(kind=DPI) :: Lv,Lc,nbsfr,NbV,NbC,LsV,LsC,id,totosys,isys
integer(kind=DPI) :: inbo,nbo,alterno(1000,2),ccoin,cvis,cvisabs,premcara,deuxcara,premlong,deuxlong
integer(kind=DPI) :: id1,id2,id3,id4,pmil_1_3,pmil_2_4
integer(kind=DPI) :: sysdevnub,nbudx,nbudy,nbudz,dom,old_dom

real(kind=DP) :: long,dlong,xl,xnmax,modunbud(3),coefid(nbmax),disperl,xdisperl
real(kind=DP) :: dok,tot,ratiol,signaxe,tensapp(6),vec(3),vec2(3),dok_final
real(kind=DP) :: densfr,ldens,xlong,xdlong,para(3),volume_micron,rpara(3)
real(kind=DP) :: yadjamax,lpboite,nbpboite,Vpboite,xdim,ydim,zdim
real(kind=DP) :: veriflong,verifangle,verifangleold,angle,addangle

real(kind=DP), allocatable :: yopdja(:,:,:),yadja(:,:,:)

integer           :: ClefBinome,segchar,dokcompt,attempt
character (LEN=30) :: cmd

logical                 :: outside,boundaries
real(kind=DP)    ,parameter :: un  = 1.D0
integer(kind=DPI),parameter :: iun = 1

! Variables for the Fibonnaci random function
integer, parameter                      :: ipp=1279, iqq=1063
integer                                 :: fibost
real(kind=DP),dimension(1-ipp:ndismax)  :: FIBONNACI

! Variables for the flow stress evaluation
integer(DPI),dimension(12,12) ::  tabjonc   !< The interaction matrix
integer(DPI)                  ::  segsys    !< The segment system number
real(DP),dimension(0:5)       ::  aij       !< The interaction matrix coefficients
real(DP),dimension(12)        ::  rhosys    !< The effective density on each slip systems
real(DP)                      ::  seglength !< The segment length
real(DP),dimension(12)        ::  TauC      !< The critical flow stress on slip systems
real(DP)                      ::  Taylor    !< The Taylor coefficient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Warning
Print *, " "
Print *, " ==========================================================================="
Print *, "                                CopyRight"
Print *, " "
Print *, " microconf is an open source program, part of the ""mM"" project."
Print *, " It was originally developed at the 'Laboratoire d'Etude"
Print *, " des Microstructure), CNRS-ONERA, FRANCE."
Print *, " Copyright (C) 1993, 1996, 2000, 2002, 2004  Benoit Devincre, Ladislas Kubin"
Print *, " Marc Condat, Christophe Lemarchand, Ronan Madec, Sebastien Groh,"
Print *, " Ghiath Monnet."
Print *, "  "
Print *, " This is free software, and you are welcome to redistribute it under the"
Print *, " conditions of the GNU General Public License (see the 'Copyright' file in"
Print *, " the program distribution)."
Print *, " ==========================================================================="
Print *, " "
Print *, " "
Print *, " Results of the computation are writen in --->   ../in/New_conf "
Print *, " "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!
! Initializations !
!!!!!!!!!!!!!!!!!!!
nbtrue      = 0
nbinome     = 99
lsc         = 0
lc          = 0
angle       = 0.
dokcompt    = 0
nbudx = 0
nbudy = 0
nbudz = 0
dok = 0.

! Def of variables used for an evaluation of the flow stress
! FCC slip system interaction matrix
! ref: B. Devincre, L. Kubin / C. R. Physique 11 (2010) 274
aij(0) = 0.122  !< Self interaction
aij(1) = 0.070  !< Hirth interaction
aij(2) = 0.137  !< Glissile interaction
aij(3) = 0.127  !< Lomer interaction
aij(4) = 0.625  !< colinear annihilation
aij(5) = 0.122  !< coplanar interaction
tabjonc(1,:) = (/0,4,2,3,2,5,3,2,5,2,1,1/)
tabjonc(2,:) = (/4,0,5,2,3,2,2,5,2,3,1,1/)
tabjonc(3,:) = (/2,5,0,4,2,3,2,5,1,1,2,3/)
tabjonc(4,:) = (/3,2,4,0,5,2,3,2,1,1,5,2/)
tabjonc(5,:) = (/2,3,2,5,0,4,1,1,2,3,5,2/)
tabjonc(6,:) = (/5,2,3,2,4,0,1,1,5,2,2,3/)
tabjonc(7,:) = (/3,2,2,3,1,1,0,4,2,5,2,5/)
tabjonc(8,:) = (/2,5,5,2,1,1,4,0,3,2,3,2/)
tabjonc(9,:) = (/5,2,1,1,2,5,2,3,0,4,3,2/)
tabjonc(10,:)= (/2,3,1,1,3,2,5,2,4,0,2,5/)
tabjonc(11,:)= (/1,1,2,5,5,2,2,3,3,2,0,4/)
tabjonc(12,:)= (/1,1,3,2,2,3,5,2,2,5,4,0/)

! The input file
open(9,FILE="../in/confinit_def",STATUS='OLD')

! Commensurability between nbase and nbasered is checked
nsys = nbase/nbasered
if (nsys * nbasered /= nbase ) then
  print*," nbase is not proportionel to nbasered !"
  stop
endif

! Implementation of the tab containing the vectors characters
do j=1,nsys,1
    do k=1,nbasered,1
        i=(j-1)*nbasered+k
        !The norme of vectors
        coefid(i)=sqrt(real(bveclin(1,i)*bveclin(1,i)+bveclin(2,i)*bveclin(2,i)+bveclin(3,i)*bveclin(3,i),DP))
    enddo
enddo

!Initialisation du generateur aleatoire
call initialisation_fibonnaci(ipp,fibonnaci)
!Seed for the random number generator
read(9,*) fibost
fibost=modulo(fibost,ndismax)   ! fibost ne peut etre > a ndismax
call gene(ipp,iqq,fibost,fibonnaci)

! The elementary parameter of the simulation lattice
xl = avalue
xl = xl * 1.0D6
write(*,*)'Size of the simulation lattice parameter (micron):',xl


!Average size of the simulated volume V in micron
read(9,*) xnmax

! The shape of the simulated volume
read(9,*) rpara(1)
read(9,*) rpara(2)
read(9,*) rpara(3)

! Dimension de la boite de simulation en unite parametre du reseau de simulation
nmax(1) = (int(xnmax/xl*rpara(1))/(facteur_boite*2))*facteur_boite*2
nmax(2) = (int(xnmax/xl*rpara(2))/(facteur_boite*2))*facteur_boite*2
nmax(3) = (int(xnmax/xl*rpara(3))/(facteur_boite*2))*facteur_boite*2

modur(1:3)=nmax(1:3)

para(1) = nmax(1)*xl
para(2) = nmax(2)*xl
para(3) = nmax(3)*xl

! Calculation of the real simulated volume if GB=1
if (GB == 1) then

  boundaries=.true.
  call load_boundaries

  write(*,*) "------------------------------------"
  write(*,*) "         Loaded boundaries          "
  write(*,*) "------------------------------------"
  write (*,fmt='(A3,A5,3A20,A15,A7,A7,A7)') &
        "#","VFP","","Miller","","VPpos","Dom","Kind"

  do I=1 ,NbPlanMax
     write (*,fmt='(A3,I5,3I20,F15.3,2I7)') &
         "#" ,i, Plane_MillerI(1:3,i), Plane_pos(i), Plane_dom(i), Plane_knd(i)
  enddo

  call vertex_calculations

  !write(cmd,'("python ../bin/volume.py")')
  cmd = "python ../bin/volume.py"

  call system(cmd)

  open(411,file="../out/volumepy.txt",STATUS='OLD')
  read(411,*) VOLUME
  close(411)

  volume_micron=volume*xl*xl*xl

  ! The tmp file ../out/volumepy.txt can now be eliminated
  !write(cmd,'("rm ../out/volumepy.txt")')
  cmd = "rm ../out/volumepy.txt"

  call system(cmd)

  ! Warning to indicate a mistake in domain definition
  if (volume_micron > para(1)*para(2)*para(3)) then
    write(*,*)
    write(*,*) "!> WARNING - Boundary Domain volume greater than Simulation Box volume  BD =", &
              & volume_micron," SB=",para(1)*para(2)*para(3)
    write(*,*)
  endif

else

  boundaries      = .false.
  volume_micron   = (para(1)*para(2)*para(3))
  domboundsX(1:2) = (/IZERO,nmax(1)/)
  domboundsY(1:2) = (/IZERO,nmax(2)/)
  domboundsZ(1:2) = (/IZERO,nmax(3)/)

endif

!Number of slice used to optimized the dislocation density homogeneity
read(9,*) nbud

! if (boundaries .and. nbud > 1) then
!   write(*,*)'Inside a finite domain the number of regions'
!   write(*,*)'to optimize dislocation (nbud) has to be set to 1'
!   Write(*,*)'imposing nbud = 1'
!   nbud=1
!   read(*,*)
! endif
xdim = domboundsX(2)-domboundsX(1)
ydim = domboundsY(2)-domboundsY(1)
zdim = domboundsZ(2)-domboundsZ(1)

print *, domboundsX,domboundsY,domboundsZ

if (xdim >= ydim .and. xdim >= zdim) then
  nbudx = nbud
  nbudy = int(nbud*ydim/xdim)+1
  nbudz = int(nbud*zdim/xdim)+1
elseif(ydim >= xdim .and. ydim >= zdim) then
  nbudy = nbud
  nbudx = int(nbud*xdim/ydim)+1
  nbudz = int(nbud*zdim/ydim)+1
elseif(zdim >= ydim .and. zdim >= xdim) then
  nbudz = nbud
  nbudy = int(nbud*ydim/zdim)+1
  nbudx = int(nbud*xdim/zdim)+1
endif

allocate(yopdja(nbudx,nbudy,nbudz))   !< Temporary density tab
allocate(yadja(nbudx,nbudy,nbudz))    !< Definitive density tab

yopdja(:,:,:) = zero
yadja(:,:,:)  = zero
nbpboite      = nbudx * nbudy * nbudz                      ! nombre de boite total d homogenisation
modunbud(1)   = (xdim)/ real(nbudx,DP)
modunbud(2)   = (ydim)/ real(nbudy,DP)
modunbud(3)   = (zdim)/ real(nbudz,DP)
              ! dimension des boites d homogenisation en reel
Vpboite       = (xdim*ydim*zdim*xl*xl*xl)/(nbpboite)    ! volume of the optimization boite
! If the subdomains volume is larger than the simulated volume then the former must be reduced to the exact value
if (Vpboite > volume_micron) Vpboite = volume_micron

write(*,*)'Dimensions of the simulation box:', nmax(1:3)
write(*,*)'Dimensions of the simulated object:', xdim,ydim,zdim
write(*,*)'Total volume of the simulated object in micron**3 = ',volume_micron
write(*,*)'number of boxes for density optimization in x, y and z direction =',nbudx,nbudy,nbudz
write(*,*)'Dimensions of boxes for density optimization is simulation unit:',modunbud*xl
write(*,*)'Total number of subdomains for density optimization =',int(nbpboite)
write(*,*)'Volume of the optimization subdomains in micron**3 = ',Vpboite

!Expected dislocation density (10^12 per micron^2)
read(9,*) densfr
ldens=1./(sqrt(densfr))        ! La longueur caracteristique pour cette densite
write(*,*)'=> ',ldens,' reference length (rau-1/2)'

!Average length of the dislocation segments in the microstructure (microns)
read(9,*) long
xlong=long/xl     ! La longueur en unite reseau de simulation

!The segment length dispersion
read(9,*) disperl

!Discretization length use to build up different dislocation character (micron)
read(9,*) dlong
if(dlong > long) then
  print*,"PAUSE, The discretisation length cannot be larger than source"
  print*,"length, Automatically set to the source length"
  read(*,*)
  dlong=long-long*0.1  !dlong must be a little smaller than long
elseif(dlong < 1.d-12) then
  print*,"PAUSE, The discretisation length cannot be smaller than zero"
  print*,"length, Automatically set to the source length"
  dlong=long-long*0.1
endif

xdlong = dlong/xl

! reference length of segments per character for each slip systems
do j=1,nsys,1
  longv(j)=int(xlong/coefid((j-1)*8+1))
  longc(j)=int(xlong/coefid((j-1)*8+3))
enddo

nbsfr = int(densfr*volume_micron/long,DPI)  ! Approximate number of source to define

! La longueur de dislocation que l on doit retrouver en moyenne par boite d homogeneisations
lpboite = ((densfr*volume_micron)/nbpboite)
write (*,*) "Length of dislocation expected in each domains of size",lpboite," micron"

lpboite = lpboite/xl  ! on se remet en unite reseau

!Tolerate error in the microstructure distribution ( >>1 means any filtering )
read (9,*) dok_final

!Type of segment microstructure: 1->sources; 2->dipoles; 3->dipolar loops; 4-> Planar loop
read (9,*) ClefBinome

!The ratio between segment length in loops or dipolar loops (needed in  case(3:4))
read (9,*) ratiol

!Type of segments character in microstructure 1 and 2 ( 1 -> screw, 2 -> edge, 3 -> random)
read (9,*) segchar

!When micro=4, The loading tensor is provided to build up only non collapsing planar loops.
read (9,*) tensapp(1:6)        ! sig11, sig22, sig33, sig32, sig31, sig21

! Combien de fois faut il tourner la boucle des binomes
select case(clefBinome)
case(1)            ! Cas d une distribution de sources
  NBinome=1
case(2)            ! Cas d un distribution de paires de sources
  NBinome=2
case(3:4)          ! Loop or dipolar loop distribution
  NBinome=1
case(:0,5:)        ! Cas non definis
  print*,' STOP The select key for the type of source distribution is wrong'
  stop
end select

! In the case of FR sources with random orientation, we must check that the discretisation length is small enough
select case(clefBinome)
  case(1:2)
    select case(segchar)
      case(:0,3)
        if(dlong > long/4) stop 'The discretisation length must be much smaller than the source length'
        if (boundaries) stop 'FR sources can only be screew or edge with boundaries'
    end select
end select

! What slip systems must be considered
!write (*,*) "systemes"
totosys=0                  ! Compteur du nombre de systemes a prendre en compte
do isys=1,nsys,1
    read (9,*) syst(isys)
    if(syst(isys)> 1 .or. syst(isys)< 0) then
      print*,"PAUSE, One slip system has a non standart activation key number"
      read(*,*)
    endif
    totosys=totosys+syst(isys)
enddo

! Nombre de paire de source par systeme de glissement
select case(clefBinome)
case(1,2)
  if (xdlong > 500) then
    write(*,*) "WARNING, the discretization length is short !!!"
    write(*,*) "You may have problem with the angular dispersion"
    read(*,*)
  endif
  nbsfrparsys = int(nbsfr / (totosys*NBinome))
  nbsysdev = 1   ! Computation does not need to be repeated for each cross-slip planes
case(3)
  ! In dipolar loop case, since each loop of sys i generate 2 edge segments of i
  ! plus 2 edge of its cross slip system, we must consider 1/(2*(1+ratiol)) Binomes
  nbsfrparsys = int(nbsfr / (real(totosys*NBinome*nbsysdev)*2.*(1.+ratiol)))
case(4)
  ! In loop case, since each loop of sys i generate 4 segments of i
  ! with a differnce between screw and edge segments equal to ratio1,
  ! we must consider 1/(4*(1+ratiol)) Binomes
  nbsfrparsys = int(nbsfr / (real(totosys*NBinome*nbsysdev)*2.*(1.+ratiol)))
  nbsysdev = 1   ! Computation does not need to be repeated for each cross-slip planes
end select

if (nbsfrparsys < 1) then
    write(*,*) "The density is too small. STOP!!!"
    STOP
endif

! message de resume
write(*,*) " Number of slip systems considered =", totosys
select case(clefBinome)
case(1,2)
  write(*,*) " Number of sources per slip system = ", nbsfrparsys
case(3,4)
  write(*,*) " Number of loops per slip systems = ", nbsfrparsys
end select

close(9)

!************************************************************************************
!*** CALCUL...
!************************************************************************************

!*Initialisation
!****************
nsegm=0
addangle=0.
attempt = 0

!*Boucle sur les sources
!***********************
bsources:   do j = 1,nbsfrparsys,1

    !*Do loop on the number of cross-slip planes
    !***********************
    bsysdev:   do sysdevnub = 1,nbsysdev,1

        !*Boucle sur les systemes  ! La boucle sur les systemes est a l interieure
        !************************  ! pour mieux distribuer a densite sur chaque systemes
        bsys: do i = 1,nsys
            if (syst(i).eq.0) cycle bsys

            !*Boucle sur le poids du systeme    ! On reboucle plusieurs fois sur une
            !meme systeme si c est utile pour augmenter son poid
            !*******************************
            !print *, " i, syst(i)  = ",i, syst(I)
            bpoids : do poidsys = 1,syst(i),1

                write(*,'("Source :", I4,"/",I4," ;  System    :",I2,"/",I2)') j,nbsfrparsys,i,nsys

                ! il faut 6 tirage pour la suite de cette boucle
                call gene(ipp,iqq,fibost,fibonnaci) ! on re-melange pour les parano
                call gene(ipp,iqq,6,fibonnaci)

                ! La longueur elementaire vis et coin des segments formant les
                ! sources est defini aleatoirement, ainsi on defini un angle
                ! quelconque
                Lc = -999999
                Lv = -999999

                select case(clefBinome)

                case(1,2)  ! Cas des sources

                  select case(segchar)

                    case(1)   !We want only screw sources
                      Lv = longv(i)
                      Lc = 0

                    case(2)   !We want only edge sources
                      Lc = longc(i)
                      Lv = 0

                    case(3)   !We want sources of random character
                      angle = 2.*pi*fibonnaci(3)
                      Lc    = int(abs(sin(angle) * longc(i)))
                      Lv    = int(abs(cos(angle) * longv(i)))

                    case(:0) !Case if we impose the source angle character
                      angle = (real(abs(segchar)) * pi) / 180. !Conversion from degrees to radians
                      Lc    = int(abs(sin(angle) * longc(i)))
                      Lv    = int(abs(cos(angle) * longv(i)))

                    case default
                      write(*,*) "The value for the source character is not correct"
                      STOP

                  end select

                case(3)    ! Cas des boucles dipolaires coin
                  angle = 2.*pi*fibonnaci(3)
                  Lc    = longc(i)
                  Lv    = 0

                case(4)    ! Case of planar loops
                  angle = 0.  ! The sign of the loop must not be defined randomly
                  Lv    = longv(i)
                  Lc    = 0

                end select

                ! Sens des caracteres coin des deux sources a definir
                if (angle < pi) then
                  UnCcoin   = nbasered*(i-1)+(nbasered/4)+1
                  DeuxCcoin = nbasered*(i-1)+3*(nbasered/4)+1   ! la coin de sens contraire
                else
                  UnCcoin   = nbasered*(i-1)+3*(nbasered/4)+1
                  DeuxCcoin = nbasered*(i-1)+(nbasered/4)+1     ! la coin de sens contraire
                endif

                ! Sens des caracteres vis des deux sources a definir
                cvisabs = nbasered*(i-1)+1
                if (angle < pisurdeux .or. angle > troispisurdeux) then
                  UnCvis    = cvisabs
                  DeuxCvis  = nbasered*(i-1)+(nbasered/2)+1  ! La vis de sens contraire
                else
                  UnCvis    = nbasered*(i-1)+(nbasered/2)+1
                  DeuxCvis  = cvisabs                        ! La vis de sens contraire
                endif

                ! On calcul une dispersion
                if (fibonnaci(4) > 0.5d0) then
                  xdisperl = (1.0d0 + (fibonnaci(5)*disperl))
                else
                  xdisperl = (1.0d0 - (fibonnaci(6)*disperl))
                endif

                ! Nombre de segments Vis et coin pour former la source
                ! avant dispersion sur les longueur (critere angulaire seul)
                NbV = int((Lv*coefid(((i-1)*8)+1))/xdlong)
                NbC = int((Lc*coefid(((i-1)*8)+3))/xdlong)

                ! C est le caractere avec le plus grand segment qui fixe
                ! la forme de discretisation. En effet, c est la condition de
                ! discretisation la plus forte, pour l autre caractere on peut
                ! allonger la longueur des segments
                if (Nbv > Nbc) then
                  Nbv = Nbc + 1
                endif
                if (Nbc > Nbv) then
                  Nbc = Nbv + 1
                endif
                if (Nbc == Nbv) then
                  Nbv = Nbc + 1
                endif

                nbo=0

                select case(clefBinome)
                case (1,2)
                  Nbo = Nbv + Nbc
                case (3,4)
                  Nbo = 4           ! Loops are made of 4 segments
                end select

                ! Calcul de la longueur des segments compatible avec le nombre de
                ! segment reconstruisant la ligne (calcule ci dessus) et la disperstion
                ! Attention a nouveau la division entiere pour les questions de sous
                ! reseau va inevitablement introduire une dispersion sur la longueur
                ! total voulue (on minor toujours)
                if (NbV /= 0) then
!                   LsV=(int(Lv*xdisperl/NbV)/facteur)*facteur ! => multiple de 8 pour cfc
                  LsV = int(Lv*xdisperl/NbV)
                else
                  LsV = 0
                endif
                if (NbC /= 0) then
!                   LsC=(int(Lc*xdisperl/NbC)/facteur)*facteur ! => multiple de 8 pour cfc
                  LsC = int(Lc*xdisperl/NbC)
                else
                  LsC = 0
                endif

                !**** BOUCLE SUR LES BINOMES
                db : do DeuxBinome = 1,NBinome
                    !print *, "DeuxBinome   = ",DeuxBinome

                    if(DeuxBinome == 1) then
                      Ccoin = UnCcoin
                      Cvis  = UnCvis
                    else
                      Ccoin = DeuxCcoin
                      Cvis  = DeuxCvis
                    endif

                    dokcompt = 0
                    acu      = 0

                    !*** Boucle sur les tentatives de trouver une coordonnee Oi qui va bien
25                  acu = acu+1

                    if (acu > 100000 .and. boundaries) then
                      print*,acu,"attempt to introduce segments outside the simulation volume have been made"
                      print*,dokcompt,"attempt to introduce segments in a fulfilled subdomains have been made"
                      print*,"You should change some parameters in confinit_def!"
                      stop
                    endif

26                  if (dokcompt > 10000) then
                      if (attempt <= 3) then
                         dok = dok + 0.3333d0*dok_final
                         attempt=attempt + 1
                         dokcompt=0
                         print*,"incresing dok value by 0.33% to finalize microstructure generation, new dok value = ",dok
                         print*,"attempt ",attempt
                      else
                        print*,"The microstructure condition of homogeneity is too high,"
                        print*,"You should increase the dok value in confinit_def"
                        stop
                      endif
                    endif

                    ! We defined the coordinate of a segment origin compatible with the simulation lattice
                    outside = .true.
                    do while (outside)

                      ! Origin of a first segment is generated
                      call gene(ipp,iqq,6,fibonnaci)
                      Oi(1) = domboundsX(1) + int(fibonnaci(1)*(domboundsX(2)-domboundsX(1)),DPI)
                      Oi(2) = domboundsY(1) + int(fibonnaci(2)*(domboundsY(2)-domboundsY(1)),DPI)
                      Oi(3) = domboundsZ(1) + int(fibonnaci(3)*(domboundsZ(2)-domboundsZ(1)),DPI)

                      ! Oi is moved to the closer simulation lattice point
                      Oi(1:3) = noeud(Oi(1:3),IUN)

                      Oi= modulo(Oi(:),nmax(:))

                      ! The coordinate of the segment origine is tested
                      if (boundaries) then
                         ListCvxdomSegI(1:NbCvxdom) = ListCvxdom(1:NbCvxdom)
                         dom = 0

                         call check_all_boundaries(real(Oi,DP),outside,dom)

                         old_dom = dom

                      else

                        if     (Oi(1) < 0)        then
                        elseif (Oi(2) < 0)        then
                        elseif (Oi(3) < 0)        then
                        elseif (Oi(1) > nmax(1)) then
                        elseif (Oi(2) > nmax(2)) then
                        elseif (Oi(3) > nmax(3)) then
                        else
                            ! Oi is inside the volume
                            outside = .false.
                        endif

                      endif

                    enddo

                    ! At this point we know that the Origin is in the simulated volume

                    ! Variable utile a la reconstruction des lignes
                    if (Nbv > Nbc) then
                      premcara = cvis
                      premlong = LsV
                      deuxcara = ccoin
                      deuxlong = LsC
                    else
                      premcara = ccoin
                      premlong = LsC
                      deuxcara = cvis
                      deuxlong = LsV
                    endif

                    ! Longueur total de chaque caractere par source
                    Lv = LsV * NbV
                    Lc = LsC * NbC

                    ! Calcul de la position de l extremite
                    ei = oi + Lv*bveclin(1:3,cvis) + Lc*bveclin(1:3,ccoin)

                    if (boundaries) then
                      ListCvxdomSegI(1:NbCvxdom) = ListCvxdom(1:NbCvxdom)
                      ListCvxdomSegI(1) = old_dom
                      ListCvxdomSEGI(old_dom) = 1
                      call check_all_boundaries(real(ei(:),DP),outside,old_dom)
                          if (outside) then
!                             print*,"--> 3"
                            goto 25
                          endif
                    endif

                    ! Vecteur source (Origine - extremite)
                    oei = ei - oi

                    ! Calcul de la taille de la source
                    veriflong = sqrt(real((oei(1)**2+oei(2)**2+oei(3)**2),DP))

                    ! Calcul de l angle de la source
                    verifangle =  ((oei(1)*bveclin(1,cvisabs)+    &
                                    oei(2)*bveclin(2,cvisabs)+    &
                                    oei(3)*bveclin(3,cvisabs)) / (veriflong*coefid(((i-1)*8)+1)) )
                    verifangle = acos ( dsign(1.d0,verifangle)*min(1.D0,dabs(verifangle)) )

                    !<<<<<<<<<<<<<
                    ! Optimization of the dislocation density distribution + intersection of boundaries test

                    ! initialization
                    iv(:) = oi(:)                 ! le premier increment part de l origine du segment
                    yopdja(:,:,:) = yadja(:,:,:)  ! le tableau temporaire = le tableau permanent

                    ! print*,' oi1 = ',oi
                    id1 = 0
                    id2 = 0
                    id3 = 0
                    id4 = 0
                    pmil     = 0
                    pmil_2_4 = 0
                    pmil_1_3 = 0
                    nbtruels = 0

                    ! Boucle incremental sur la longueur des segments formant une ligne
                    do inbo = 1,nbo,1       ! case (1,2) FR source are made of nbo screw and edge segments
                                            ! case (3,4) loops are made of 4 segments
                        id = 0
                        pmil = 0
                        select case(clefBinome)
                        case (1,2)

                          ! On alterne le type de segment le long de la ligne
                          if (((inbo/2)*2) /= inbo) then
                            alterno(inbo,1) = premcara
                            alterno(inbo,2) = premlong
                          else
                            alterno(inbo,1) = deuxcara
                            alterno(inbo,2) = deuxlong
                          endif

                          ! Initialization
                          id        = alterno(inbo,1)           ! caractere du segment
                          pmil      = alterno(inbo,2)           ! longueur du segment
                          nbtruels  = 1                         ! le nombre de vecteur elementaire par boite
                          ivtemp(:) = modulo(iv(:),nmax(:))
                          oindb(1)  = int((ivtemp(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          oindb(2)  = int((ivtemp(2)-domboundsY(1))/modunbud(2))+1
                          oindb(3)  = int((ivtemp(3)-domboundsZ(1))/modunbud(3))+1

                        case (3)

                          ! Initialization
                          id        = premcara                  ! caractere du segment
                          pmil      = premlong                  ! longueur du segment
                          nbtruels  = 1                         ! le nombre de vecteur elementaire par boite
                          ivtemp(:) = modulo(iv(:),nmax(:))
                          oindb(1)  = int((ivtemp(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          oindb(2)  = int((ivtemp(2)-domboundsY(1))/modunbud(2))+1
                          oindb(3)  = int((ivtemp(3)-domboundsZ(1))/modunbud(3))+1

                          select case(inbo)
                          case (1)
                            ! We must select the signe of the loop = select id2, id3, id4
                            call gene(ipp,iqq,6,fibonnaci)
                            if (fibonnaci(1).gt.0.5d0) then
                              id2=nbasered*(sysdev(int(id/nbasered)+1,sysdevnub)-1)+(nbasered/4)+1    ! 1er edge CS
                              id3=DeuxCcoin                                ! Le segment coin binome
                              id4=nbasered*(sysdev(int(id/nbasered)+1,sysdevnub)-1)+3*(nbasered/4)+1  ! 2er edge CS
                            else
                              id2=nbasered*(sysdev(int(id/nbasered)+1,sysdevnub)-1)+3*(nbasered/4)+1  ! 1er edge CS
                              id3=DeuxCcoin                                ! Le segment coin binome
                              id4=nbasered*(sysdev(int(id/nbasered)+1,sysdevnub)-1)+(nbasered/4)+1    ! 2er edge CS
                            endif
                          case (2)
                            id = id2
                            pmil  = int(pmiL*ratiol)
                          case (3)
                            id = id3
                          case (4)
                            id = id4
                            pmil  = int(pmiL*ratiol)
                          end select

                        case (4)

                          ! Initialization
                          id        = premcara                  ! caractere du segment
                          pmil      = premlong                  ! longueur du segment
                          nbtruels  = 1                         ! le nombre de vecteur elementaire par boite
                          ivtemp(:) = modulo(iv(:),nmax(:))
                          oindb(1)  = int((ivtemp(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          oindb(2)  = int((ivtemp(2)-domboundsY(1))/modunbud(2))+1
                          oindb(3)  = int((ivtemp(3)-domboundsZ(1))/modunbud(3))+1

                          select case(inbo)
                          case (1)
                            ! We must select the signe of the loop
                            vec(1) = (tensapp(1)*bveclin(1,id) + tensapp(6)*bveclin(2,id) + tensapp(5)*bveclin(3,id))
                            vec(2) = (tensapp(6)*bveclin(1,id) + tensapp(2)*bveclin(2,id) + tensapp(4)*bveclin(3,id))
                            vec(3) = (tensapp(5)*bveclin(1,id) + tensapp(4)*bveclin(2,id) + tensapp(3)*bveclin(3,id))
                            vec2(1) = (vec(2)*bveclin(3,id)-vec(3)*bveclin(2,id))
                            vec2(2) = (vec(3)*bveclin(1,id)-vec(1)*bveclin(3,id))
                            vec2(3) = (vec(1)*bveclin(2,id)-vec(2)*bveclin(1,id))
                            signaxe= (vec2(1)*bveclin(1,unCcoin)+vec2(2)*bveclin(2,unCcoin)+vec2(3)*bveclin(3,unCcoin))

                            call gene(ipp,iqq,2,fibonnaci)
                            if (signaxe .gt. 0.d0) then
                              if (fibonnaci(1) .gt. 0.5d0) then
                                id1 = id         ! 1 segment to build the loop
                                id2 = DeuxCcoin  ! 2 segment to build the loop
                                id3 = DeuxCvis   ! 3 segment to build the loop
                                id4 = UnCcoin    ! 4 segment to build the loop
                                pmil_1_3 = pmiL
                                pmil_2_4 = int(pmiL*longc((id2/8)+1)/longv((id2/8)+1)*ratiol*0.5)*2

                              else
                                id1 = DeuxCvis   ! 1 segment to build the loop
                                id2 = UnCcoin    ! 2 segment to build the loop
                                id3 = id         ! 3 segment to build the loop
                                id4 = DeuxCcoin  ! 4 segment to build the loop
                                pmil_1_3 = pmiL
                                pmil_2_4 = int(pmiL*longc((id2/8)+1)/longv((id2/8)+1)*ratiol*0.5)*2
                              endif
                            else
                              if (fibonnaci(2).gt.0.5d0) then
                                id1 = id         ! 1 segment to build the loop
                                id2 = unCcoin    ! 2 segment to build the loop
                                id3 = DeuxCvis   ! 3 segment to build the loop
                                id4 = DeuxCcoin  ! 4 segment to build the loop
                                pmil_1_3 = pmiL
                                pmil_2_4 = int(pmiL*longc((id2/8)+1)/longv((id2/8)+1)*ratiol*0.5)*2
                              else
                                id1 = DeuxCvis   ! 1 segment to build the loop
                                id2 = DeuxCcoin  ! 2 segment to build the loop
                                id3 = id         ! 3 segment to build the loop
                                id4 = unCcoin    ! 4 segment to build the loop
                                pmil_1_3 = pmiL
                                pmil_2_4 = int(pmiL*longc((id2/8)+1)/longv((id2/8)+1)*ratiol*0.5)*2
                              endif
                            endif

                            id    = id1
                            pmil  = pmil_1_3
                          case (2)
                            id    = id2
                            pmil  = pmil_2_4
                          case (3)
                            id    = id3
                            pmil  = pmil_1_3
                          case (4)
                            id    = id4
                            pmil  = pmil_2_4
                          end select

                        end select

                        ! Boucle sur la longueur de chaque segment avec un pas
                        ! elementaire afin de localiser les troncons de ligne
                        do tnl = 1, pmil, 1

                            pmilmun = tnl-1
                            it(:) = iv(:) + pmilmun*bveclin(:,id)   ! coordonne de chaque increment du segment

                            ! if a section of the segment is found outside the boundaries we must stop
                            if (boundaries) then
                              ListCvxdomSegI(1:NbCvxdom) = ListCvxdom(1:NbCvxdom)
                              ListCvxdomSegI(1) = old_dom
                              ListCvxdomSegI(old_dom) = 1
                              call check_all_boundaries(real(it(:),DP),outside,old_dom)
                              if (outside) then
!                                 print*,"--> 1"

                                goto 25
                              endif
                            endif

                            it(:) = modulo(it(:),nmax(:))           ! les CLP

                            ! test de localisation (indice de boite)
                            indb(1)  = int((it(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            indb(2)  = int((it(2)-domboundsY(1))/modunbud(2))+1
                            indb(3)  = int((it(3)-domboundsZ(1))/modunbud(3))+1

                            ! y a t il un changement d indice = changement de boite, ou fin de segment
                            if (oindb(1) /= indb(1) .or.    &
                                oindb(2) /= indb(2) .or.    &
                                oindb(3) /= indb(3) .or.    &
                                (pmilmun == (pmil-1))) then
                              ! Si la densite dans une boite depasse la limite fixee, la ligne n est
                              ! plus acceptable et il faut refaire un tirage
                              if (((yopdja(oindb(1),oindb(2),oindb(3)) + (nbtruels*coefid(id)*xl)) / Vpboite) >    &
                                  (densfr * (1.+dok))) then
                                dokcompt = dokcompt + 1       ! Number of wrong attempt
                                goto 26
                              endif

                              ! on incremente la longueur dans la boite car on vient de
                              ! changer  de boite ou de terminer le segment
                              yopdja(oindb(1),oindb(2),oindb(3)) = yopdja(oindb(1),oindb(2),oindb(3)) +    &
                                                                   (nbtruels * coefid(id) * xl)
                              nbtruels = 0        ! on initialise si on change de boite ou si on est en fin de segment
                              oindb(:) = indb(:)  ! on memorise l indice de boite
                            endif

                            nbtruels = nbtruels + 1 ! conteur de tranche elementaire trouve dans oindb()

                        enddo

                        ! The end of the segment must also be tested with boundaries
                        if (boundaries) then
                          it(:) = iv(:) + pmil*bveclin(:,id)   ! coordonne de chaque increment du segment
                          ListCvxdomSegI(1:NbCvxdom) = ListCvxdom(1:NbCvxdom)
                          ListCvxdomSegI(1) = old_dom
                          ListCvxdomSEGI(old_dom) = 1
                          call check_all_boundaries(real(it(:),DP),outside,old_dom)
                          if (outside) then
!                             print*,"--> 3"

                            goto 25
                          endif
                        endif

                        ! Origine du nouveau segment
                        if (boundaries) then
                          iv(:) = iv(:)+pmil*bveclin(:,id)
                        else
                          iv(:) = modulo(iv(:)+pmil*bveclin(:,id),nmax(:))
                        endif

                    enddo

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! The dislocation made of nbo segments is OK, we can save the result !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    ! Iv is set again to the first segment origine
                    iv(:) = oi(:)

                    id = 0
                    oindb(:) = -1
                    ! Il faut boucler une deuxieme fois sur les nbo pour sauver la solution validee ci-dessus
                    do inbo = 1,nbo,1

                        select case(clefBinome)

                        !!!!!!!!!!!!
                        case(1,2)  ! Cas des sources

                          ! Initialization
                          id        = alterno(inbo,1)
                          pmil      = alterno(inbo,2)
                          nbtruels  = 1

                          ! The new iv value
                          if (inbo > 1) then
                             iv(:) = iv(:) + segini(4,nsegm) * bveclin(:,segini(5,nsegm))
                             iv(:) =  modulo(iv(:),nmax(:))
                          endif

                          oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                          oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1

                          nsegm = nsegm + 1
                          do k=1,3
                            segini(k,nsegm)=iv(k)
                          enddo

                          segini(4,nsegm) = pmiL
                          segini(5,nsegm) = id

                          if (inbo /= 1) then
                            segini(6,nsegm) = -1
                          else
                            segini(6,nsegm) = 0
                          endif

                          if (inbo /= nbo) then
                            segini(7,nsegm) = -1         ! version sans table de voisin
                          else
                            segini(7,nsegm) = 0
                          endif

                        !!!!!!!!!!!
                        case(3)   ! Dipolar loop, The asymetry is provided thanks to ratiol

                          ! Initialization
                          id        = premcara                  ! caractere du segment
                          pmil      = premlong                  ! longueur du segment
                          nbtruels  = 1                         ! le nombre de vecteur elementaire par boite
                          oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                          oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1

                          select case(inbo)
                          case (1)
                            !------------
                            ! 1st segment
                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = pmiL
                            segini(5,nsegm)   = id
                            segini(6,nsegm)   = -3
                            segini(7,nsegm)   = -1
                          case (2)
                            !-----------
                            ! 2d segment
                            ! notice that the length of edge segments must be taken pair for some crystal symmetry

                            iv(:) = modulo(iv(:) + pmil*bveclin(:,id),nmax(:))

                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = int(pmiL*ratiol)
                            segini(5,nsegm)   = id2
                            segini(6,nsegm)   = -1
                            segini(7,nsegm)   = -1

                            id    = id2
                            pmil  = int(pmiL*ratiol)
                            oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                            oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1
                          case (3)
                            !-----------
                            ! 3d segment

                            iv(:) = modulo((iv(:)+int(pmil*ratiol)*bveclin(:,id2)),nmax(:))

                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = pmiL
                            segini(5,nsegm)   = id3
                            segini(6,nsegm)   = -1
                            segini(7,nsegm)   = -1

                            id    = id3
                            oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                            oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1
                          case (4)
                            !-----------
                            ! 4d segment
                            ! notice that the length of edge segments must be taken pair for some crystal symetry

                            iv(:) = modulo(iv(:) + pmil*bveclin(:,id3),nmax(:))

                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = int(pmiL*ratiol)
                            segini(5,nsegm)   = id4
                            segini(6,nsegm)   = -1
                            segini(7,nsegm)   = -3

                            id    = id4
                            pmil  = int(pmiL*ratiol)
                            oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                            oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1
                          end select

                        !!!!!!!!!!!
                        case(4)   ! planar loop

                          ! Initialization
                          id        = id1                       ! caractere du segment
                          pmil      = pmiL_1_3                  ! longueur du segment
                          nbtruels  = 1                         ! le nombre de vecteur elementaire par boite
                          oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                          oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1

                          select case(inbo)
                          case (1)
                            !------------
                            ! 1st segment
                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = pmiL_1_3
                            segini(5,nsegm)   = id1
                            segini(6,nsegm)   = -3
                            segini(7,nsegm)   = -1
                          case (2)
                            !------------
                            ! 2d segment
                            ! notice that the length of edge segments must be taken pair for some crystal symmetry

                            iv(:)=modulo(iv(:)+pmiL_1_3*bveclin(:,id1),nmax(:))

                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = pmil_2_4
                            segini(5,nsegm)   = id2
                            segini(6,nsegm)   = -1
                            segini(7,nsegm)   = -1

                            id    = id2
                            pmil  = pmil_2_4
                            oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                            oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1
                          case (3)
                            !-----------
                            ! 3d segment

                            iv(:)=modulo(iv(:)+segini(4,nsegm)*bveclin(:,segini(5,nsegm)),nmax(:))

                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = pmiL_1_3
                            segini(5,nsegm)   = id3
                            segini(6,nsegm)   = -1
                            segini(7,nsegm)   = -1

                            id    = id3
                            pmil  = pmil_1_3
                            oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                            oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1
                          case (4)
                            !-----------
                            ! 4d segment
                            ! notice that the length of edge segments must be taken pair for some crystal symetry

                            iv(:)=modulo(iv(:)+segini(4,nsegm)*bveclin(:,segini(5,nsegm)),nmax(:))

                            nsegm = nsegm + 1
                            segini(1:3,nsegm) = iv(1:3)
                            segini(4,nsegm)   = pmil_2_4
                            segini(5,nsegm)   = id4
                            segini(6,nsegm)   = -1
                            segini(7,nsegm)   = -3

                            id    = id4
                            pmil  = pmil_2_4
                            oindb(1)  = int((iv(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                            oindb(2)  = int((iv(2)-domboundsY(1))/modunbud(2))+1
                            oindb(3)  = int((iv(3)-domboundsZ(1))/modunbud(3))+1
                          end select

                        end select


                        ! The definitive density distribution in subdomains is saved
                        do tnl = 1, pmil, 1

                          pmilmun  = tnl-1
                          it(:)    = iv(:) + pmilmun*bveclin(:,id)
                          it(:)    = modulo(it(:),nmax(:))

                          indb(1)  = int((it(1)-domboundsX(1))/modunbud(1))+1  ! numero boite contenant l origine du segment
                          indb(2)  = int((it(2)-domboundsY(1))/modunbud(2))+1
                          indb(3)  = int((it(3)-domboundsZ(1))/modunbud(3))+1

                          if (oindb(1) /= indb(1) .or.    &
                              oindb(2) /= indb(2) .or.    &
                              oindb(3) /= indb(3) .or.    &
                              (pmilmun == (pmil-1))) then

                            if(((yadja(oindb(1),oindb(2),oindb(3)) + (nbtruels*coefid(id)*xl)) / Vpboite) >    &
                                (densfr * (1+dok))) then
                              write (*,*) 'pb ici',(yadja(oindb(1),oindb(2),oindb(3))+(nbtruels*coefid(id)*xl)) /  &
                                  Vpboite, (densfr * (1+dok))
                              write (*,*) yadja(oindb(1),oindb(2),oindb(3)),(nbtruels*coefid(id)*xl)
                              stop
                            endif

                            ! The density in boite oindb is increased
                            yadja(oindb(1),oindb(2),oindb(3)) = yadja(oindb(1),oindb(2),oindb(3)) +   &
                                                                (nbtruels*coefid(id)*xl)

                            nbtruels = 0
                            oindb(:) = indb(:)

                          endif

                          nbtruels = nbtruels + 1
                        enddo

                    enddo

                    do ix=1,nbudx,1
                      do iy=1,nbudy,1
                        do iz=1,nbudz,1
                          !test calcul du taut de remplissage de la boite totale de simulation
                          if (yadja(ix,iy,iz) /= yopdja(ix,iy,iz)) then
                          print *, "problem with density calculation",yadja(ix,iy,iz),yopdja(ix,iy,iz)
                          stop

                          endif
                        enddo
                      enddo
                    enddo

                    if(DeuxBinome.eq.1) verifangleold = verifangle
                    addangle = addangle + verifangle

                enddo db
                !**** BOUCLE SUR LES BINOMES

                verifangle = verifangle + verifangleold
                !write (*,*) 'Somme des deux angles :',verifangle

            enddo bpoids
            !*Boucle sur le poids du systeme
            !*******************************

        enddo bsys
        !*Boucle sur les systemes
        !************************

    enddo bsysdev
    !*Do loop on the number of cross-slip planes
    !************************

enddo bsources
!*Boucle sur les sources
!***********************

!*** FORMAT
3 format(2x,10(i7,2x))

!*** ECRITURE DU FICHIERS DE SOURCES DE FRANCK-READ
open(20,file='../in/New_conf',status='unknown')

! write(20,'(30i2)') (syst(i),i=1,nsys)   ! Only the slip systems present in the microstructure are activated
write(20,'(30i2)') (iun,i=1,nsys)         ! All slip systems are activated
write(20,*) nsegm
write(20,*) nmax

if (crystal_structure=="HCP") then
  do i=1,nsegm
      if (segini(5,i) < 9)   write(2,3) (segini(j,i),j=1,7)
  enddo
  do i=1,nsegm
      if (segini(5,i) > 8)   write(2,3) (segini(j,i),j=1,7)
  enddo
else
  do i=1,nsegm
      write(20,61) i,segini(1,i),segini(2,i),segini(3,i),segini(4,i),segini(5,i),            &
          segini(6,i),segini(6,i),segini(7,i),segini(7,i),.false.,0,0
  enddo
61 format(1x,I5,2x,3I10,2x,I5,2x,I7,2x,4I7,2x,L3,2x,I7,2x,I3)
endif
!*** Rappel des parametres utilises en fin de fichier ***

write (20,*)"------------------------------------------------------------------"
write (20,*)"------------------------------------------------------------------"
write (20,*)'Latice simulation parameter = ',xl
write (20,*)'Size of the simulation box (microns) = ',xnmax
write (20,*)'Parallelepipede (Non : 0, Oui : PARA) : ',para
write (20,*)'Number of domains used to homogenize the density (NbUd) = ',nbud
write (20,*)'Density (*10^12) = ',densfr
write (20,*)'Length of sources (microns), dispertion), number of sources (total and per slip systems) = '&
            &,long,disperl,nbsfr,NBSFRPARSYS
write (20,*)'Tolerance : ',dok
write (20,*)'Angular dispersion of character ',addangle/((nbsfr/2)*2)

!*** Donnees statistiques ***
stat(:)=0
nbtrue=0
tot=0.

yadjamax=0.d0
do ix=1,nbudx,1
  do iy=1,nbudy,1
    do iz=1,nbudz,1
      ! calcul du taut de remplissage de la boite totale de simulation
      if (yadja(ix,iy,iz).ne.0.) nbtrue = nbtrue + 1
      !  on note la densite maxi vu dans une boite
      if (yadja(ix,iy,iz).gt.yadjamax) yadjamax = yadja(ix,iy,iz)
    enddo
  enddo
enddo

write (20,*)'Max : ',yadjamax/Vpboite

! calcul de l histograme de la densite min a max sur 10 tranches
do ix = 1,nbudx,1
  do iy = 1,nbudy,1
    do iz = 1,nbudz,1
        incr = int(yadja(ix,iy,iz)/yadjamax)*9 + 1
        if(incr > 10 .or. incr < 1) then
          write (*,*) "saturation",ix,iy,iz,incr
          stop
        endif
        stat(incr) = stat(incr) + 1
        tot        = tot + yadja(ix,iy,iz)
    enddo
  enddo
enddo

write (20,*) "Ratio of occupation in domains = ",float(nbtrue)/nbpboite
write (20,*) "Histogramme of density distribution in domains ",nbpboite
do ix = 1,10,1
  write (20,*) ix,'<', (real(ix)/9.*yadjamax)/Vpboite, float(stat(ix))/nbpboite*100,'%'
enddo
write (*,*) "Effective density =",tot/Volume_micron,' 10^12'
write (20,*) "Effective density =",tot/Volume_micron,' 10^12'

! Calculation of the effective density on each slip system
rhosys(1:12) = 0.
do i = 1,nsegm
  segsys = segini(5,i)/nbasered + 1
  seglength = segini(4,i) * coefid(segini(5,i))
  rhosys(segsys) = rhosys(segsys) + seglength
enddo
rhosys(1:12) = rhosys(1:12) * avalue / (Volume_micron * 1.d-18)

write (20,*) "Effective density on the slip systems (10^12) !"
do i=1,12
  write (20,*) rhosys(i) / 1.e12
enddo

! Calculation of the Taylor coefficient associated to the microstructure
do i = 1,12
  TauC(i) = 0.
  do j = 1,12
    TauC(i) = TauC(i) + aij(tabjonc(i,j)) * rhosys(j)
  enddo
  TauC(i) = sqrt(TauC(i))
enddo
! TauC must be calculated only on slip system with a non-zero density
do i = 1,12
  if (rhosys(i) == 0.) TauC(i) = 1.D20
enddo

print*,TauC(1:12)
Taylor = minval(TauC(1:12)) / sqrt(sum(rhosys(1:12)))
write (*,*) "The effective forest coefficient of the microstructure =", Taylor
write (20,*) "The effective forest coefficient of the microstructure =", Taylor

close(20)

end subroutine microconf


!###############################################
!<
!###############################################
subroutine check_all_boundaries(O,outside,dom)

implicit none

logical :: insideDO(NbCvxDom),insideDOfix(NbCvxDom),insideDOfree(NbCvxDom)
logical :: outside,outsharedO, insideO(NBplanDom)

real(kind=DP)         :: TestO,normint(3),O(3)
integer(kind=DPI)     :: dom,old_dom,JJ,ii,hh
integer(kind=DPI)     :: planes,freesurf

outside = .true.
insideO(:)=.false.
insideDO(:)=.false.
insideDOfix(:)=.false.
insideDOfree(:)=.false.

old_dom = dom

if (Nbplan > IZERO) then

  bc1: do hh = 1,Nbcvxdom

    ii = ListCvxdomsegi(hh)

    planes=IZERO

    if (ii > IUN) then
      do jj=1,ii-1
        planes   = planes  +  NbPlanCvxDom(1,jj)
      enddo
    endif

    if (kkdebug) then
      write(379,*) " "
      write (379,fmt='("####---|  TESTING DOMAIN ", I2)') ii
      write(379,*) ""
      write (379,fmt='("-- Checking FIXED boundaries")')
    endif

    outsharedO=.false.

    do jj= planes+1,planes+NbPlanCvxDom(1,ii)

      ! Normale projection of OiP on the barrier

      normint(:) = Plane_MillerR(1:3,jj)

      TestO = O(1)*normint(1)+O(2)*normint(2)+O(3)*normint(3)-Plane_pos(jj)

      if (kkdebug) then
        write (379,*) ""
        write (379,fmt='(A7,I7)') "VFP : ",JJ
        write (379,fmt='(3x," TestO",F30.3)') TestO
      endif

      if (Plane_knd(jj) == IZERO) then !If impenetrable barrier

        if (TestO >= -numtol_dotp) then ! if the segment is outside the barrier

          InsideDO(ii)  = .false.

          insideDOfix(ii) =.false.

          if (kkdebug) write (379,fmt='("!/ Segment cross fixed boundary ", I3 ," kind = 0 with both ends")') jj
          !cycle bc1

        elseif (TestO < -numtol_dotp) then
          insideO(jj) = .true.

        endif

      elseif (Plane_knd(jj) == IUN ) then !If shared barrier

        if ( TestO > zero  ) then

          InsideDO(ii)  = .false.

          insideDOfix(ii) =.false.

          outsharedO=.true.

          if (kkdebug) write (379,fmt='("!/ Segment cross shared boundary ", I3 ,"  kind = 1 with both ends")') jj

        elseif (TestO <= zero) then
          insideO(jj) = .true.

        endif

      endif

    enddo

    if (All(insideO(planes+1:planes+NbPlanCvxDom(1,ii))) &
        ) then
        insideDOfix(ii) = .true.
        if(NbfreePlan == IZERO) then
          insideDO(ii) = .true.
          dom = ii
          outside = .false.
          return
        endif
        if (kkdebug) write(379,*) "insideDOfix ii ==", ii ,"case 1"
    endif

    if ((.not. insideDOfix(ii) .and. .not. outsharedO .and. ii == Old_dom)) then
       outside = .true.
       if (kkdebug) then
          write(379,*) "BLOQUER: outside fixed barrier 1"
          write(379,*) .not. insideDOfix(ii) , outsharedO , ii, old_dom
       endif
    return
    endif
 enddo bc1

 if (.not. any(insideDOfix(:))) then
    outside = .true.
    if (kkdebug) write(379,*) "BLOQUER: outside fixed barrier 2"
    return
 endif

endif

bc2: do hh = 1,Nbcvxdom

  ii = ListCvxdomSegi(hh)

  if (NbfreePlan > IZERO) then

    freesurf=Nbplan

    if (ii > IUN) then
      do jj=1,ii-1
        freesurf  = freesurf +  NbPlanCvxDom(2,jj)
      enddo
    endif

    if (kkdebug) then
      write(379,*) " "
      write (379,fmt='("####---|  TESTING DOMAIN ", I2)') ii
      write(379,*) ""
      write (379,fmt='("-- Checking FREE boundaries")')
    endif

    do jj =Freesurf+1,Freesurf+NbPlanCvxDom(2,ii)

        if (kkdebug) then
          write(379,*) "freesurf",freesurf
          write(379,*) "NbPlanCvxDom(2,ii)",NbPlanCvxDom(2,ii)
        endif

       normint(:) = Plane_MillerR(1:3,jj)

       TestO = O(1)*normint(1)+O(2)*normint(2)+O(3)*normint(3)-Plane_pos(jj)

        if (kkdebug) then
          write (379,*) ""
          write (379,fmt='(A7,I5)') "VFP : ",JJ
          write (379,fmt='(3x," TestO",F30.3)') TestO
        endif

        if (TestO >= -numtol_dotp) then

        if (kkdebug) write (379,fmt='("!/ Segment cross free surface  ", I3 ," with both ends")') jj

          insideO(jj)=.false.

        elseif (TestO < -numtol_dotp) then

          insideO(jj) = .true.

        endif

      enddo

      if(kkdebug) then
        write(379,*) "freesurf",freesurf
        write(379,*) "NbPlanCvxDom(2,ii)",NbPlanCvxDom(2,ii)
        write(379,*) "Freesurf",All(insideO(Freesurf+1:Freesurf+NbPlanCvxDom(2,ii)))
      endif
      if (All(insideO(Freesurf+1:Freesurf+NbPlanCvxDom(2,ii)))) insideDOfree(ii)=.true.

    endif

    if ((insideDOfree(ii) .and. NbPlan == IZERO)     .or. &
        (insideDOfix(ii) .and.  insideDOfree(ii))) then
        dom = ii
        insideDO(ii) = .true.
        outside = .false.
        if (kkdebug) write(379,*) "insideDO ii ==", ii ,"case 4"
        return
    endif

enddo bc2

end subroutine


!################################################
!<
!################################################
subroutine initialisation_fibonnaci(ipp,harvest)
! attention les valeurs de ipp et iqq qui marche bien
! sont respectivement 1279 et 1063

integer ::ipp,i
real(kind=8) harvest(1-ipp:ndismax)

!----------Lecture du fichier aleatoire--------------------
open(10,file='../in/mM_Data_Files/random_dat',FORM='FORMATTED',STATUS='OLD')

read(10,1000)(harvest(i-ipp),i=1,ipp)

1000 format(10e14.7)
close(10)

call gene(1279,1063,ndismax,harvest)

end subroutine initialisation_fibonnaci


!####################################
!<
!####################################
subroutine gene(ipp,iqq,itet,random)
! attention les valeurs de ipp et iqq qui marche bien
! sont respectivement 1279 et 1063

integer :: ipp,iqq,itet,i
real(kind=8) random(1-ipp:itet)

do i=1,itet
  random(i)=mod((random(i-ipp)+random(i-iqq)),1.d0)
end do

do i=1,ipp
  random(i-ipp)=random(itet-ipp+i)
end do

end subroutine gene


!##########################################################################
!< subroutine noeud(o) : renvoie le point le plus proche de O situe sur
!< le reseau principal de la simulation
!##########################################################################
Function noeud(Oi,homo)

integer (kind=DPI) :: Oi(3),noeud(3)
integer (kind=DPI) :: coefX,coefY,coefZ,homo
real(DP)           :: Xnorme(3),Ynorme(3),Znorme(3)  ! les 3 vecteurs elementaire du reseau NORMES
real(DP)           :: n(3),x,y,z,XY(3),tmp(3)
real(DP)           :: rcompx(3), rcompy(3), rcompz(3)

!print*,x_reseau,y_reseau,z_reseau
!read(*,*)

Xnorme(:) = Real(X_reseau(:),DP) / norivect(x_reseau)
Ynorme(:) = Real(Y_reseau(:),DP) / norivect(Y_reseau)
Znorme(:) = Real(Z_reseau(:),DP) / norivect(Z_reseau)

n(:) = prodvec(Xnorme,Ynorme)
! calculation of component of (Oi) on Z_reseau
z =  dot_product(n,real(Oi,DP)) / dot_product(Znorme,n)  ! component on z_axis : unite = z_reseau normed
z = z / norivect(Z_reseau)              ! component on z_axis : unite = z_reseau
rcompZ(:) = z * z_reseau(:)

XY(:) = Oi(:) - rcompz(:)    ! calculation of component of (Oi) in the X-Y plane

tmp(:) = prodvec(n,Xnorme)             ! vector in XY perpendicular to I
if(abs(dot_product(XY,n)) > numtols) then
   print *," dot_product(xy,n) =",dot_product(xy,n)
   stop
endif

y = dot_product(tmp,XY) / dot_product(Ynorme,tmp)  ! component on y_axis : unite = y_reseau normed
y = y / norivect(y_reseau)              ! component on y_axis : unite = y_reseau
rcompY(:) = Y * Y_reseau(:)
rcompX(:) = XY(:) - rcompY(:)

if((abs(dot_product(rcompx,tmp)) > numtols) .or. &
(abs(dot_product(rcompx,n)) > numtols)) then
   print *," dot_product(rcompx,tmp) ",dot_product(rcompx,tmp)
   print *," dot_product(rcompx,n)",dot_product(rcompx,n)
   stop
endif

X = dot_product(rcompX,Real(X_reseau(:),DP)) / norivect(X_reseau)**2    ! component on x_axis : unite = y_reseau normed

! calculation of integer coeficient in the net elementary vectors
coefX = NINT(X,DPI)   ! integer component on x_axis : unite = z_reseau
coefY = NINT(Y,DPI)   ! integer component on y_axis : unite = z_reseau
coefZ = NINT(Z,DPI)   ! integer component on z_axis : unite = z_reseau

! Coefs have to be multiple of homo :
coefX =  (coefX/homo) * homo
coefY =  (coefY/homo) * homo
coefZ =  (coefZ/homo) * homo

! fournir noeud(:) : le point le plus proche de Oi
noeud(:) = coefX * x_reseau(:) + coefy * y_reseau(:) + coefz * z_reseau(:)
if(norivect(noeud(:)-Oi(:)) > trois * norivect(z_reseau)) then
   write(*,'(" Oi ini :", 3I15," HOMO = ",I4)') Oi(:) , HOMO
   write(*,'("x,rcompx : ",4F15.2)') x,rcompx(:)
   write(*,'("y,rcompy : ",4F15.2)') y,rcompy(:)
   write(*,'("z,rcompz : ",4F15.2)') z,rcompz(:)
   write(*,'("Oi reel : ",3F15.2)') rcompx(:) + rcompy(:) +rcompz(:)
   write(*,'("Node    : ",3I15)')noeud(:)
   print *, X, y ,z
   write(*,'("coef    : ",3I15)')coefX, coefy , coefz
   PRINT *, "noeud(:)-Oi(:) =",noeud(:)-Oi(:)
   stop
endif

end function noeud


end module INITCONFTOOLS
