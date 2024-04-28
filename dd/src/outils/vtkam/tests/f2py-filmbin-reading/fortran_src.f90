MODULE FILM

USE varglob

IMPLICIT NONE

CONTAINS

!-------------------------------------------
!---
!-------------------------------------------

subroutine readfilmheader(cristallo,avalue,bveclin,bvecnor)

implicit none

integer                           ,intent(out)  :: cristallo
real(kind=DP)                     ,intent(out)  :: avalue

integer(kind=DPI),dimension(3,96), intent(out)         :: bveclin
integer(kind=DPI),dimension(3,96), intent(out)         :: bvecnor

integer             :: nsegmaxf
integer(kind=DPI)   :: i
integer(kind=DPI)   :: kk
integer(kind=DPI)   :: nbase
integer(kind=DPI)   :: KODERR

integer(kind=DPI),dimension(3)            :: modur

character           :: fichier*50

kk=1

!*** Ouverture de l unite contenant le film
fichier = "film.bin"

!open(49,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
open(49,file='../../../../../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
read(49) cristallo    ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
read(49) nsegmaxf     ! The nsegmax value used in this simulation
read(49) avalue       ! The simulation lattice parameter (the simulation scaling)
read(49) nbase        ! Nb of vectors direction used for the material crystal symmetry
read(49) (bveclin(:,I),bvecnor(:,I), I=1,nbase)    ! table of vector direction and displacement for this simulation
read(49) modur(:)     ! Dimension of the 3D simulated volume in a value

if (nsegmax /= nsegmaxf) then
  write(*,*) "The value of nsegmax for film.bin is: ", nsegmaxf, " or nsegmax in vtkam/varglob is ", nsegmax
  write(*,*) "You must compile again 'vtkam' with the correct value of nsegmax : STOP !"
  stop
endif

end subroutine readfilmheader

!-------------------------------------------
!---
!-------------------------------------------

subroutine readfilmstep(kkextract, nsegm)

implicit none

integer(kind=DPI),intent(in)    :: kkextract
integer(kind=DPI),intent(out)   :: nsegm

integer(kind=DPI),allocatable,dimension(:)   :: junctt
integer(kind=DPI),allocatable,dimension(:)   :: surface

!integer(kind=DPI),dimension(nsegmax)   :: junctt
!integer(kind=DPI),dimension(nsegmax)   :: surface

integer(kind=DPI),allocatable,dimension(:,:) :: iseg

real(kind=DP)       :: xstress
real(kind=DP)       :: epso
real(kind=DP)       :: deltat

integer(kind=DPI)   :: kk
integer(kind=DPI)   :: j
integer(kind=DPI)   :: KODERR

integer(kind=DPI),allocatable,dimension(:,:)    :: TABVOIS

character           :: caractest

kk = kkextract - 1

Steps_loop: do while (kk.lt.kkextract)

  ! A new simulation step
  read(49,IOSTAT=KODERR) kk,deltat,NSEGM,XSTRESS,EPSO
  print *, kk
  if(KODERR < 0) then
    print*,"  End of the film.bin file"
    stop
  elseif(KODERR > 0) then
    print*,"Reading problem 1",KODERR,kk
    stop
  endif

  If (kk >= kkextract) then
    write(*,'(I7, ", Eps:",F8.4,"%, Sigma =", F8.2," MPa, Nseg = ",I5)') KK,EPSO,XSTRESS, nsegm
    exit Steps_loop
  endif

  allocate(junctt(nsegm))
  allocate(surface(nsegm))
  allocate(iseg(5,nsegm))
  allocate(TABVOIS(2,nsegm))

  ! The segments micro structure is loaded
  READ(49,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),surface(j), J=1,nsegm)

  deallocate(junctt)
  deallocate(surface)
  deallocate(iseg)
  deallocate(TABVOIS)

  if(KODERR < 0) then
    print*,"  End of the film.bin file"
    stop
  elseif(KODERR > 0) then
    print*,"Reading Error readfilmstep",KODERR,nsegm,j
    stop
  endif

  read(49,IOSTAT=KODERR) caractest

  if(KODERR < 0) then
    print*,"  End of the film.bin file"
    stop
  elseif(KODERR > 0) then
    print*,"Reading error 3",caractest
    stop
  endif

enddo Steps_loop

end subroutine readfilmstep

!-------------------------------------------
!---
!-------------------------------------------

subroutine readfilmstepdata(mode, bveclin,bvecnor, nsegm, avalue, ldis_act, slipsyscurv,curvact,  &
kneepoints, points, slipsystem,loop,ntotpoints,ntotlines,nptsctr,seg_center,curv_vect)

implicit none

integer(kind=DPI),intent(in)   :: mode
integer(kind=DPI),intent(in)   :: nsegm
integer(kind=DPI),dimension(3,96), intent(in)   :: bveclin
integer(kind=DPI),dimension(3,96), intent(in)   :: bvecnor
integer(kind=DPI)                , intent(in)   :: slipsyscurv
Logical                          , intent(in)   :: curvact
real(kind=DP)                    , intent(in)   :: Ldis_act,avalue

!integer(kind=DPI),dimension(nsegmax), intent(out)   :: O1,O2,O3
real(DP),dimension(3,2*nsegm), intent(out) :: points
real(DP),dimension(3,2*nsegm), intent(out) :: kneepoints
integer(DPI), dimension(nsegm+1,nsegm), intent(out) :: loop

integer(DPI), intent(out)                         :: nptsctr
real(DP)    , intent(out)  ,dimension(3,nsegm)   :: seg_center
real(DP)    , intent(out)  ,dimension(3,nsegm)   :: curv_vect

integer(DPI), intent(out)  :: ntotpoints
integer(DPI), intent(out)  :: ntotlines
integer(DPI), intent(out)  :: slipsystem(2*nsegm)

integer(kind=DPI),dimension(nsegm)   :: junctt
integer(kind=DPI),dimension(nsegm)   :: surface
integer(kind=DPI),dimension(5,nsegm) :: iseg

real(kind=DP)    ,dimension(3,nsegm)      :: Vcourb

integer(DPI)  :: ntotkneepoints
integer(DPI)  :: is
integer(DPI)  :: js
integer(DPI)  :: b1(3)
integer(DPI)  :: b2(3)
integer(DPI)  :: dislo

integer(kind=DPI)   :: j
integer(kind=DPI)   :: KODERR

integer(kind=DPI),dimension(2,nsegm)    :: TABVOIS

character           :: caractest

points(:,:)         = 0.
kneepoints(:,:)     = 0.
slipsystem(:)       = 0
ntotpoints          = 0
ntotkneepoints      = 0
ntotlines           = 0

VCourb(1:3,1:nsegm)   = 0.0            ! The curvature vector is initialized to zero

! The segments micro structure is loaded
READ(49,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),surface(j), J=1,nsegm)

if(KODERR < 0) then
  print*,"  End of the film.bin file"
  stop
elseif(KODERR > 0) then
  print*,"Reading Error readfilmstepdata",KODERR,nsegm,j
  stop
endif
read(49,IOSTAT=KODERR) caractest
if(KODERR < 0) then
  print*,"  End of the film.bin file"
  stop
elseif(KODERR > 0) then
  print*,"Reading error 3",caractest
  stop
endif

If (mode == 1) then

  do dislo = 1, nsegm
    if (Iseg(4,dislo) == 0) then
      ntotkneepoints = ntotkneepoints + 1
      kneepoints(:,ntotkneepoints + 1) = real(Iseg(1:3,dislo),DPI)
    else

      ntotlines = ntotlines + 1
      loop(ntotlines,1)= 2

      ntotpoints=ntotpoints+1
      points(:,ntotpoints) = real(iseg(1:3,dislo),DPI)
      loop(ntotlines,2)= ntotpoints-1

      ntotpoints=ntotpoints+1
      points(:,ntotpoints) = real(iseg(1:3,dislo),DPI) + real((iseg(4,dislo)-1)*bveclin(1:3,iseg(5,dislo)),DPI)
      loop(ntotlines,3)= ntotpoints-1

      !--- coloring
      if (junctt(dislo) == 0) then !we assign a color to each slip system
        slipsystem(ntotpoints) = int((iseg(5,dislo)-1)/8)+1
      else
        if (JUNCTT(dislo) == -1) then
          slipsystem(ntotpoints)= 20  !gd segment
        else
          is=(int((iseg(5,dislo)-1)/8)*8+1)
          js=(int((iseg(5,JUNCTT(dislo))-1)/8)*8+1)
          b1(1:3)=bveclin(1:3, is)
          b2(1:3)=bveclin(1:3, js)
          if (bvecnor(1,is)==bvecnor(1,js) .and. bvecnor(2,is)==bvecnor(2,js) &
            .and. bvecnor(3,is)==bvecnor(3,js)) then
            slipsystem(ntotpoints)=19 !coplanar interaction
          elseif (dot_product(b1,b2) == 0) then
            slipsystem(ntotpoints)=18 ! Hirth Lock
          elseif (dot_product((b1+b2),bvecnor(1:3,is)) == 0 .or. &
                  dot_product((b1+b2),bvecnor(1:3,js)) == 0) then
            if (dot_product((b1+b2),bveclin(1:3,iseg(5,dislo))) == 0) then
                slipsystem(ntotpoints)=16 !edge glissile
            else
                slipsystem(ntotpoints)=17 !60 glissile
            endif
          else
              if (dot_product((b1+b2),bveclin(1:3,iseg(5,dislo))) == 0) then
                slipsystem(ntotpoints)=14 !90 Lomer-Cotrell
              else
                slipsystem(ntotpoints)=15 !60 Lomer-Cotrell
              endif
          endif ! if (bvecnor(1,is)==bvecnor(1,js) .and. bvecnor(2,is)==bvecnor(2,js) &
        endif !if (junctt(i) == -1)
      endif ! if (junctt(i) == 0)
      !---
      slipsystem(ntotpoints-1) = slipsystem(ntotpoints)
    endif

  enddo
  kneepoints(1,1) = ntotkneepoints !Store number of kneecaps at (1,1) position

elseif (mode == 2) then
  if (curvact) call calcurv(Iseg,Junctt,tabvois,surface,nsegm,bveclin,avalue,Ldis_act,Vcourb)

  call reorder_microstructure(nsegm, slipsyscurv,bveclin,bvecnor,iseg,tabvois,slipsystem, junctt, &
points,loop,ntotpoints,ntotlines,nptsctr,seg_center,curv_vect,Vcourb,curvact)

endif

end subroutine readfilmstepdata

!-------------------------------------------
!---
!-------------------------------------------

subroutine reorder_microstructure(nsegm,slipsyscurv, bveclin,bvecnor,iseg,tabvois,slipsystem, &
  & junctt, points,loop, ntotpoints,ntotlines,nptsctr,seg_center,curv_vect, vcourb, curvact)

implicit none

integer(kind=DPI), intent(in)             :: nsegm    !< Number of segments (Replica or not)
integer(kind=DPI),  dimension(3,96)       :: bveclin
integer(kind=DPI),  dimension(3,96)       :: bvecnor
integer(kind=DPI),  dimension(5,nsegm)  :: iseg
integer(kind=DPI)                         :: slipsyscurv
real(kind=DP)    ,dimension(3,nsegm)      :: Vcourb

integer(kind=DPI),  dimension(2,nsegm)         :: tabvois
integer(kind=DPI),  dimension(nsegm)           :: junctt
integer(DPI)     ,  dimension(2*nsegm)         :: slipsystem
integer(DPI)     ,  dimension(nsegm+1,nsegm), intent(inout)  :: loop

integer(DPI)  :: i
integer(DPI)  :: Ei(3)
integer(DPI)  :: dislo

integer(DPI)  :: ntotpoints
integer(DPI)  :: ntotlines
integer(DPI)  :: linepoints
integer(DPI)  :: npointdislo
integer(DPI)  :: is
integer(DPI)  :: js
integer(DPI)  :: b1(3)
integer(DPI)  :: b2(3)
integer(DPI)  :: nptsctr

real(DP),dimension(3,2*nsegm)  :: points
real(DP),dimension(3)           :: EiR
real(DP),dimension(3,nsegm)    :: seg_center
real(DP),dimension(3,nsegm)    :: curv_vect
logical ,dimension(nsegm)      :: zipseg

Logical    :: curvact

! reduction of the maximum number of loops and segments to reduce the size memory allocation
print *, "fdeb5"
zipseg(1:nsegm) = .false.
slipsystem(1:2*nsegm) = 0
!loop(1:limloop,1:limlines) = 0
loop(:,:) = 0

!to write the .vtp file the connectivity is rebuilt

if (nsegm < nsegmax) then

  ! 1. We discard segments with zero lenght
  do dislo = 1, nsegm
    if(iseg(4,dislo) == 0 .and. junctt(dislo)==0)  zipseg(dislo) = .true. !kneecaps are not considered, they have zero length so we don't display it
    if(iseg(4,dislo) == 0 .and. junctt(dislo)/=0) then  !GD segments are treated as kneecaps, but in ths case we need to rebuild the connectivity
      zipseg(dislo)=.true.
      tabvois(1,tabvois(2,dislo)) = tabvois(1,dislo)
      tabvois(2,tabvois(1,dislo)) = tabvois(2,dislo)
    endif
  enddo

  ntotpoints = 0
  ntotlines = 0
  nptsctr = 0

  ! 2. connectivity of pinned lines (not closed loop)
  do dislo = 1, nsegm
    if (zipseg(dislo)) cycle
    if (iseg(4,dislo) /= 0 .and. tabvois(1,dislo) == nsegmax .and. .not. zipseg(dislo)) then
        linepoints=0
        i = dislo

        ntotlines=ntotlines+1
          Do while (i /= nsegmax)

            zipseg(i) = .true.


            ntotpoints=ntotpoints+1
            linepoints=linepoints+1
            points(1,ntotpoints)=iseg(1,i)
            points(2,ntotpoints)=iseg(2,i)
            points(3,ntotpoints)=iseg(3,i)

            print *,i,points(1:3,ntotpoints),iseg(4,i), nsegmax
            loop(ntotlines,linepoints+1)=ntotpoints-1

            if (junctt(i) == 0) then !we assign a color to each slip system
              slipsystem(ntotpoints) = int((iseg(5,i)-1)/8)+1

            else
              if (JUNCTT(i) == -1) then
                slipsystem(ntotpoints)= 20  !gd segment
              else
                is=(int((iseg(5,i)-1)/8)*8+1)
                js=(int((iseg(5,JUNCTT(i))-1)/8)*8+1)
                b1(1:3)=bveclin(1:3, is)
                b2(1:3)=bveclin(1:3, js)
                if (bvecnor(1,is)==bvecnor(1,js) .and. bvecnor(2,is)==bvecnor(2,js) &
                  .and. bvecnor(3,is)==bvecnor(3,js)) then
                  slipsystem(ntotpoints)=19 !coplanar interaction
                elseif (dot_product(b1,b2) == 0) then
                  slipsystem(ntotpoints)=18 ! Hirth Lock
                elseif (dot_product((b1+b2),bvecnor(1:3,is)) == 0 .or. &
                        dot_product((b1+b2),bvecnor(1:3,js)) == 0) then
                  if (dot_product((b1+b2),bveclin(1:3,iseg(5,i))) == 0) then
                      slipsystem(ntotpoints)=16 !edge glissile
                  else
                      slipsystem(ntotpoints)=17 !60 glissile
                  endif
                else
                    if (dot_product((b1+b2),bveclin(1:3,iseg(5,i))) == 0) then
                      slipsystem(ntotpoints)=14 !90 Lomer-Cotrell
                    else
                      slipsystem(ntotpoints)=15 !60 Lomer-Cotrell
                    endif
                endif ! if (bvecnor(1,is)==bvecnor(1,js) .and. bvecnor(2,is)==bvecnor(2,js) &
              endif !if (junctt(i) == -1)
            endif ! if (junctt(i) == 0)

            if (tabvois(2,i)/=nsegmax) then
              if ((int((iseg(5,i)-1)/8)+1 /= int((iseg(5,tabvois(2,i))-1)/8)+1 &
                .or. (junctt(i)==0 .and. junctt(tabvois(2,i)) > 0) .or. &
                  (junctt(i) > 0 .and. junctt(tabvois(2,i)) == 0)) .and. iseg(4,i)>1) then

                EiR(1:3) = real(iseg(1:3,i),DPI) + real((iseg(4,i)-1)*bveclin(1:3,iseg(5,i)),DPI)
                ntotpoints=ntotpoints+1
                linepoints=linepoints+1
                points(1:3,ntotpoints)=EiR(1:3)

                loop(ntotlines,linepoints+1)=ntotpoints-1

                slipsystem(ntotpoints) = slipsystem(ntotpoints-1)

              endif

              if ((slipsyscurv == int((iseg(5,i)-1)/8)+1 .or. slipsyscurv >12) .and. curvact) then
                nptsctr=nptsctr+1
                seg_center(1:3,nptsctr)=real((2.d0*iseg(1:3,i)+iseg(4,i)*bveclin(1:3,iseg(5,i)))*0.5d0,DPI)
                curv_vect(1:3,nptsctr) =VCOURB(1:3,i)
              endif

              Ei(1:3) = iseg(1:3,i) + iseg(4,i)*bveclin(1:3,iseg(5,i))
              i=tabvois(2,i)

              !to handle PBC, we cut the line, adding a point

              if (Ei(1) /= iseg(1,i) .or. Ei(2) /= iseg(2,i) .or. Ei(3) /= iseg(3,i) .or. i == nsegmax) then

                  ntotpoints=ntotpoints+1
                  linepoints=linepoints+1
                  points(1,ntotpoints)=Ei(1)
                  points(2,ntotpoints)=Ei(2)
                  points(3,ntotpoints)=Ei(3)
                  loop(ntotlines,linepoints+1)=ntotpoints-1
                  loop(ntotlines,1)=linepoints
                  linepoints=0
                  ntotlines=ntotlines+1
                  slipsystem(ntotpoints) = slipsystem(ntotpoints-1)


              endif
            else
              i=tabvois(2,i)
            endif


          enddo

        if (linepoints==0) then
          ntotlines=ntotlines-1
        else
          loop(ntotlines,1)=linepoints
        endif

      endif

  enddo

  !3. connectivity of closed loop

  do dislo = 1, nsegm

      if (zipseg(dislo)) cycle

      i = dislo
      linepoints=0
      npointdislo=ntotpoints+1  !to close the loop
      ntotlines=ntotlines+1

      Do while (.not. zipseg(i) .and. i /= nsegmax)
            zipseg(i) = .true.

            !origin

            ntotpoints=ntotpoints+1
            linepoints=linepoints+1
            points(1:3,ntotpoints)=real(iseg(1:3,i),DPI)

            loop(ntotlines,linepoints+1)=ntotpoints-1

            if (junctt(i) == 0) then !we assign a color to each slip system
                slipsystem(ntotpoints) = int((iseg(5,i)-1)/8)+1
              else
                if (JUNCTT(i) == -1) then
                slipsystem(ntotpoints)= 20  !gd segment
              else
                is=(int((iseg(5,i)-1)/8)*8+1)
                js=(int((iseg(5,JUNCTT(i))-1)/8)*8+1)
                b1(1:3)=bveclin(1:3, is)
                b2(1:3)=bveclin(1:3, js)
                if (bvecnor(1,is)==bvecnor(1,js) .and. bvecnor(2,is)==bvecnor(2,js) &
                  .and. bvecnor(3,is)==bvecnor(3,js)) then
                  slipsystem(ntotpoints)=19 !coplanar interaction
                elseif (dot_product(b1,b2) == 0) then
                  slipsystem(ntotpoints)=18 ! Hirth Lock
                elseif (dot_product((b1+b2),bvecnor(1:3,is)) == 0 .or. &
                        dot_product((b1+b2),bvecnor(1:3,js)) == 0) then
                  if (dot_product((b1+b2),bveclin(1:3,iseg(5,i))) == 0) then
                      slipsystem(ntotpoints)=16 !edge glissile
                  else
                      slipsystem(ntotpoints)=17 !60 glissile
                  endif
                else
                    if (dot_product((b1+b2),bveclin(1:3,iseg(5,i))) == 0) then
                      slipsystem(ntotpoints)=14 !90 Lomer-Cotrell
                    else
                      slipsystem(ntotpoints)=15 !60 Lomer-Cotrell
                    endif
                endif ! if (bvecnor(1,is)==bvecnor(1,js) .and. bvecnor(2,is)==bvecnor(2,js) &
              endif !if (junctt(i) == -1)
            endif

            if ((int((iseg(5,i)-1)/8)+1 /= int((iseg(5,tabvois(2,i))-1)/8)+1 &
              .or. (junctt(i)==0 .and. junctt(tabvois(2,i)) > 0) .or. &
                (junctt(i) > 0 .and. junctt(tabvois(2,i)) == 0)) .and. iseg(4,i) > 1) then

              EiR(1:3) = real(iseg(1:3,i),DPI) + real((iseg(4,i)-1)*bveclin(1:3,iseg(5,i)),DPI)
              ntotpoints=ntotpoints+1
              linepoints=linepoints+1
              points(1:3,ntotpoints)=EiR(1:3)

              loop(ntotlines,linepoints+1)=ntotpoints-1

              slipsystem(ntotpoints) = slipsystem(ntotpoints-1)

            endif


            if ((slipsyscurv == int((iseg(5,i)-1)/8)+1 .or. slipsyscurv >12) .and. curvact) then
              nptsctr=nptsctr+1
              seg_center(1:3,nptsctr)=real((2.d0*iseg(1:3,i)+iseg(4,i)*bveclin(1:3,iseg(5,i)))*0.5d0,DPI)
              curv_vect(1:3,nptsctr) =VCOURB(1:3,i)
            endif

            Ei(1:3) = iseg(1:3,i) + iseg(4,i)*bveclin(1:3,iseg(5,i))

            i=tabvois(2,i)

            !to handle PBC, we cut the line, adding a point

            if (Ei(1) /= iseg(1,i) .or. Ei(2) /= iseg(2,i) .or. Ei(3) /= iseg(3,i)) then !&

                ntotpoints=ntotpoints+1
                linepoints=linepoints+1
                points(1,ntotpoints)=Ei(1)
                points(2,ntotpoints)=Ei(2)
                points(3,ntotpoints)=Ei(3)
                loop(ntotlines,linepoints+1)=ntotpoints-1
                loop(ntotlines,1)=linepoints


                slipsystem(ntotpoints) = slipsystem(ntotpoints-1)


                linepoints=0
                ntotlines=ntotlines+1

            elseif (i==dislo) then
                ntotpoints=ntotpoints+1
                points(1:3,ntotpoints)=real(iseg(1:3,dislo),DPI)
                linepoints=linepoints+1
                loop(ntotlines,linepoints+1)=ntotpoints-1  !npointdislo-1 !to close the loop
                slipsystem(ntotpoints)=slipsystem(npointdislo)
            endif
      enddo

      if (linepoints==0) then
        ntotlines=ntotlines-1
      else
        loop(ntotlines,1)=linepoints
      endif

  enddo

  else
    print *, " Segments information are not listed since Nsegm > Nsegmax"
endif

end subroutine reorder_microstructure


subroutine calcurv(Iseg,Junctt,tabvois,surface,nsegm,bveclin,avalue,Ldis_act,Vcourb)

implicit none

integer(kind=DPI),intent(in)   :: nsegm

real(kind=DP),dimension(3,nsegm)       :: Vcourb

logical,dimension(nsegm)          :: tabGD    !< The GD marker
logical,dimension(nsegm)          :: tabjonc  !< The jonction marker

integer(kind=DPI),dimension(nsegm)   :: junctt
integer(kind=DPI),dimension(nsegm)   :: surface
integer(kind=DPI),dimension(5,nsegm) :: iseg
integer(kind=DPI),dimension(2,nsegm)    :: TABVOIS

! Variables used in the vector curvature calculation
integer(kind=DPI)   :: P1(3),P2(3),P3(3),Oi(3),Ei(3),VNN,compteur
integer(kind=DPI)   :: CompteurG, CompteurD   !< for radius of curvature statistics
integer(kind=DPI)   :: ivli, vli, long
integer(kind=DPI)   :: Rlocal, Rlocal2
integer(kind=DPI)   :: BUFTLRO(1:3,3), I1O, I1E, NBDIM(3), DECAL(3)
real(kind=DP)       :: Ldis_act, alpha
real(kind=DP)       :: norme        !< for radius of curvature statistics

integer (kind=DPI)                  :: XLOMAX           !< The discretization length

integer(kind=DPI),dimension(3,96)  :: bveclin

integer(kind=DPI)  :: i

real(kind=DP)       :: avalue

real, parameter     :: fac_rlocal = 0.5   !< The factor applied to the rlocal cste (see vector curvature calculation)

integer(kind=DPI),  dimension(3)    :: ROO, trAB, Tr
real(kind=DP),      dimension(3)    :: RA, RB

! Defintion of the regularization length used in the curvature definition for
! the local line tension. This length is scaled by the discretization length
! All the segments participating to the local line tension will be excludded
! from the interactions (self-interaction)
Ldis_act  = Ldis_act * 1D-6               ! the discretization length in meter
Xlomax    = NINT(Ldis_act/avalue,DPI)     ! the discretization length in the simulation units
rlocal    = INT(xlomax * fac_rlocal,DPI)  ! the reference length used for the curvature (line tension) calculation
rlocal2   = rlocal * rlocal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculation of the vector radius of curvature for the segments
  ! This calculation is essentially a copy and paste of the calculation made in mM
  ! to calculate the local line tension.
  ! The calculation is made only in the reference box without replicas

  ! The junction and GD information is recovered
  tabGD(1:nsegm)   = .false.
  tabjonc(1:nsegm) = .false.

  do i = 1, nsegm
    if (JUNCTT(i) < 0) tabGD(i)   = .true.   ! This segment is a GD segment
    if (JUNCTT(i) > 0) tabjonc(i) = .true.   ! this segment is part of a junction
  enddo

  do i = 1, nsegm

    long          = ISEG(4,I)
    VLI           = ISEG(5,I)
    i1o           = tabvois(1,I)
    i1e           = tabvois(2,I)

    if (long == 0) cycle     ! The curvature vector cannot be calculated on pivot segments

    ! DETERMINATION of the points needed for the calculation of the local radius of curvature:
    ! P2 is always the center of I
    P2(1:3) = iseg(1:3,i) + (long * bveclin(1:3,VLI)) / 2

    ! P1 is :
    ! - the origin of I, if I1o is a pinning point
    ! - the center of I1o, if I1o is a real segment (not a part of a long segment discretized in small pieces)
    ! - the center of the effective neighbor before discretization
    if (i1o  == nsegmax .or. tabjonc(i1o)) then

      P1(1:3) = iseg(1:3,i)

    else
      ! Oi : is the beginning of the effective I1o
      ! Ei : is the end of the effective I1o
      Ei(1:3) = iseg(1:3,i)

      ! to avoid PBC problem, all coordinates are calculated relatively to the origin of I
      Tr(1:3) = iseg(4,I1o) * bveclin(1:3,ISEG(5,i1o)) / 2
      Oi(1:3) = Ei(1:3) - Tr(1:3) * 2

      compteur = 1

      if (RLOCAL2 < 0.0) then
        ! P1 is the center of the vnn segment
        P1(1:3) = Oi(1:3) + Tr(1:3)
      else

        ! The vector connecting P2 and the center of i1o
        ROO(1:3) = (Oi(1:3) + Tr(1:3)) - P2(1:3)

        ! Do we need to look for more distant segment i10 to regularize the line tension curvature ?
        if (SUM(real(ROO,DP)*real(ROO,DP)) < Rlocal2) then

          !The half length of i1o segment
          TrAB(1:3) = Tr(1:3)

          ! Ra and Rb are vectors between P2 and the extremities of i1o
          RA(1:3)=real(ROO(1:3) + TrAB(1:3),DP)
          RB(1:3)=real(ROO(1:3) - TrAB(1:3),DP)

          vnn = i1o

          ! exact calculation of the line tension impose that the local radius of curvature is calculated
          ! with points P1 and P3 not very close from P2. A minimum distance of Rlocal is imposed!

          do while (tabvois(1,vnn) /= nsegmax  .and.          &
                    .not. tabjonc(vnn)       .and.          &
                    .not. tabgd(vnn)         .and.          &
                    SUM(RA*RA) < Rlocal2     .and.          &
                    SUM(RB*RB) < Rlocal2)

            ! The next i1o segment along the dislocation line is tested
            vnn = tabvois(1,vnn)
            compteur = compteur + 1

            ! The vector connecting P2 and the center of previous vnn
            Tr(1:3) = TrAB(1:3)

            !The half length of the new vnn segment
            IVLI = iseg(5,vnn)
            TrAB(1:3) = (iseg(4,vnn)*bveclin(1:3,IVLI)) / 2

            !The center of the new vnn segment
            ROO(1:3) = ROO(1:3) - (Tr(1:3) + TrAB(1:3))

            ! Ra and Rb are vectors between P2 and the extremities of new i1o
            RA(:)=real(ROO(:)+TrAB(:),DP)
            RB(:)=real(ROO(:)-TrAB(:),DP)

            ! We exit the loop if more than 10 segments were tested
            if(compteur > 10) exit

          enddo

          ! P1 is the center of the vnn segment
          P1(1:3) = P2(1:3) + ROO(1:3)
        endif

      endif

      ! Now we look for the problem of long straight lines cut in many segments,
      ! straight lines section must be treated as one long segment, whatever rlocal
      if (compteur == 1) then

        ! The initial vnn is redefined
        vnn = tabvois(1,i1o)

        do while (vnn /= nsegmax .and. iseg(5,vnn) == iseg(5,i1o))

          Oi(1:3) =  Oi(1:3) - iseg(4,vnn) * bveclin(1:3,iseg(5,vnn))
          compteur = compteur + 1
          vnn = tabvois(1,vnn)

          ! no more than 10 aligned segments are considered
          if(compteur > 10) exit

        enddo

        P1(1:3) = (Oi + Ei) / 2     ! The new P1

      endif

      compteurG = compteur  ! the compteur value is saved for statistics

    endif  ! End of the tests on the P1 side

    ! P3 is :
    ! - the end     of I if I1e is pinning point
    ! - the center  of I1e if I1e is a real segment (not a part of a long segment beeing descretised)
    ! - the center of the effctive neighbor before descretisation
    if (i1e  == nsegmax .or. tabjonc(i1e)) then

      P3(1:3) = iseg(1:3,i) + long * bveclin(1:3,VLI)  ! end of I

    else

      ! begining of effective I1e = end of I
      Oi(1:3) = iseg(1:3,i) + long * bveclin(1:3,VLI)

      ! starting end of I1e
      Tr(:) =   iseg(4,i1e) * bveclin(1:3,iseg(5,i1e)) / 2
      Ei(1:3) = Oi(1:3) + Tr(1:3) * 2

      compteur = 1

      if (Rlocal2 < 0.0) then

        ! P3 is the center of the vnn segment
        P3(1:3) = Oi(1:3) + Tr(:)

      else

        ! The vector connecting P2 and the center of i1e
        ROO(1:3) = (Oi(1:3) + Tr(1:3)) - P2(1:3)

        ! Do we need to look for more distant segment i1e to regularize the line tension curvature ?
        if (SUM(real(ROO,DP)*real(ROO,DP)) < Rlocal2) then

          !The half length of i1e segment
          TrAB(:) = Tr(1:3)

          ! Ra and Rb are vectors between P2 and the extremities of i1o
          RA(:) = real(ROO(:) + TrAB(:),DP)
          RB(:) = real(ROO(:) - TrAB(:),DP)

          vnn = i1e

          ! exact calculation of the line tension impose that the local radius of curvature is calculated
          ! with points P1 and P3 not very close from P2. A minimum distance of Rlocal is imposed!
          do while (tabvois(2,vnn) /= nsegmax       .and.     &
                    .not.tabjonc(vnn)               .and.     &
                    .not.tabgd(vnn)                 .and.     &
                    SUM(RA*RA) < Rlocal2            .and.     &
                    SUM(RB*RB) < Rlocal2)

            ! The next i1e segment along the dislocation line is tested
            vnn = tabvois(2,vnn)
            compteur = compteur + 1

            ! The vector connecting P2 and the center of previous vnn
            Tr(1:3) = TrAB(1:3)
            !The half length of the new vnn segment
            IVLI = iseg(5,vnn)
            TrAB(1:3) = (iseg(4,vnn)*bveclin(1:3,IVLI)) / 2

            ! The vector connecting P2 and the center of vnn
            ROO(1:3) = ROO(1:3) + (Tr(1:3) + TrAB(1:3))

            ! Ra et Rb distances entre le centre de i et les extremites de vnn
            RA(:)=real(ROO(:) - TrAB(:),DP)
            RB(:)=real(ROO(:) + TrAB(:),DP)

            if(compteur > 10) exit

          enddo

          ! P3 is the center of the vnn segment
          P3(1:3) = P2(1:3) + ROO(1:3)

        endif

      endif

      ! Now we look for the problem of long straight lines cut in many segments
      ! straight lines section must be treated as one long segment, whatever rlocal
      if (compteur == 1) then

        ! The initial vnn is redefined
        vnn = tabvois(2,i1e)

        do while(vnn /= nsegmax .and. iseg(5,vnn) == iseg(5,i1e))

          EI(1:3) = Ei(1:3) + iseg(4,vnn) * bveclin(1:3,iseg(5,vnn)) ! end of effective I1e
          compteur = compteur + 1

          vnn = tabvois(2,vnn)
          if(compteur > 10) exit

        enddo

        P3(1:3) = (Oi + Ei) / 2

      endif

      CompteurD = Compteur ! We save compteur value for statistics

    endif     ! End of the P3 tests

    if (surface(i) == 0) then    ! The curvature vector has no meaning for the segments touching a surface or interface

      ! test on the existence of a local curvature (VCourb /= infinity)
      Decal(1:3) = P1(1:3) - P2(1:3)
      NBdim(1:3) = P3(1:3) - P2(1:3)

      ! the following test is true only when P1,P2,P3 are alligned
      if ((Decal(2)*NBdim(3) /= Decal(3)*NBdim(2)) .or.          &
          (Decal(3)*NBdim(1) /= Decal(1)*NBdim(3)) .or.          &
          (Decal(1)*NBdim(2) /= Decal(2)*NBdim(1)))    then

        BUFTLRO(1:3,1) = P1(1:3)
        BUFTLRO(1:3,2) = P2(1:3)
        BUFTLRO(1:3,3) = P3(1:3)

        ! Calcul du vecteur rayon de courbure
        VCOURB(1:3,i) = VECCOURB(BUFTLRO(1:3,1),BUFTLRO(1:3,2),BUFTLRO(1:3,3), alpha)

        norme = sqrt(VCOURB(1,i)*VCOURB(1,i) + &
                     VCOURB(2,i)*VCOURB(2,i) + &
                     VCOURB(3,i)*VCOURB(3,i) ) ! avec la fonction intrinsèque dot_product

        !the following lines are commented, they are useful to visualize curvature
        !vector in paraview: in case of large curvature vectors you can impose a threshold (default=10 times the discretization length)
        !to reduce the maximum size of the vectors in the graphical interface, to improve the visualization.

        if (norme > Xlomax*20) then
           VCOURB(1:3,i)=VCOURB(1:3,i)*Xlomax*20/norme
        endif


      endif

    endif

    !213  format(2x,3(I4,2x),2x, E13.5,2x, L3, 9(I4,2x), 2(E13.5,2x))
    !write(28,213) i,              &           ! numéro du noeud
    !             compteurG,      &           ! nb voisins gauche
    !              compteurD,      &           ! nb voisins gauche
    !              norme                       ! norme of the curvature_radius vector
  enddo     ! End of the segment loop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine calcurv


FUNCTION VECCOURB(OM,OP,OQ,alpha)

implicit none

real(kind=DP), dimension(3)            :: VECCOURB
Integer(kind=DPI), dimension(3),intent(in) :: OM,OP,OQ
real(kind=DP)                          :: ALPHA1, ALPHA2
real(kind=DP),intent(out)              :: ALPHA
real(kind=DP), dimension(3)            :: PM,PQ
real(kind=DP), dimension(3)            :: A, B, VA,VB, PMvecPQ
real(kind=DP)                          :: m
real(kind=DP)                          :: InvnormPM, InvnormPQ, normVC
real(kind=DP), parameter               :: half = 0.5D0

PM(1:3) = real(OM(1:3)-OP(1:3),DP)
PQ(1:3) = real(OQ(1:3)-OP(1:3),DP)

  PMvecPQ(1) = (PM(2)*PQ(3) - PM(3)*PQ(2))
  PMvecPQ(2) = (PM(3)*PQ(1) - PM(1)*PQ(3))
  PMvecPQ(3) = (PM(1)*PQ(2) - PM(2)*PQ(1))
  VA(1) = (PMvecPQ(2)*PM(3) - PMvecPQ(3)*PM(2))
  VA(2) = (PMvecPQ(3)*PM(1) - PMvecPQ(1)*PM(3))
  VA(3) = (PMvecPQ(1)*PM(2) - PMvecPQ(2)*PM(1))
  VB(1) = (PMvecPQ(2)*PQ(3) - PMvecPQ(3)*PQ(2))
  VB(2) = (PMvecPQ(3)*PQ(1) - PMvecPQ(1)*PQ(3))
  VB(3) = (PMvecPQ(1)*PQ(2) - PMvecPQ(2)*PQ(1))
  A(1:3) = OP(1:3) + half * PM(1:3)
  B(1:3) = OP(1:3) + half * PQ(1:3)
  m = (VA(2)*(A(1)-B(1)) + VA(1)*(B(2)-A(2))) / ((VB(1)*VA(2))-(VB(2)*VA(1)))
  VECCOURB(1) = (B(1) + m * VB(1)) - OP(1)
  VECCOURB(2) = (B(2) + m * VB(2)) - OP(2)
  VECCOURB(3) = (B(3) + m * VB(3)) - OP(3)

  normVC = sqrt(VECCOURB(1)**2 + VECCOURB(2)**2 + VECCOURB(3)**2)


  InvNormPM  = 1.0 / sqrt(PM(1)**2 + PM(2)**2 + PM(3)**2)
  InvNormPQ  = 1.0 / sqrt(PQ(1)**2 + PQ(2)**2 + PQ(3)**2)
  PM(1:3) =   PM(1:3) * InvNormPM
  PQ(1:3) =   PQ(1:3) * InvNormPQ
  ALPHA1  = 2 * Asin(half/(InvNormPM*normVC))
  ALPHA2  = 2 * Asin(half/(InvNormPQ*normVC))
  ALPHA   = ALPHA1 + ALPHA2

END FUNCTION VECCOURB

END MODULE FILM
