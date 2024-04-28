!------------------------------
SUBROUTINE INFOS
write (*,*) 'This option is not a function of cam !'
END SUBROUTINE INFOS
!-------------------------------

SUBROUTINE MICROMEGAS

implicit none

integer,parameter   :: DPI = selected_int_kind(9)
integer,parameter   :: DP = selected_real_kind(p=14)
integer,parameter   :: nsegmax = 100000

real(kind=DP)       :: xstress,epso,avalue,deltat
integer             :: cristallo,nsegmaxf
integer(kind=DPI)   :: i,ijk,kk,j,nbsegdep,freqpause,debutaffich,freqaffich
integer(kind=DPI)   :: nsegm,iseg,bveclin(3,96),bvecnor(3,96),TABVOIS(2,nsegmax),  &
                       junctt,modur(3),modur_ref(3),nbase,surface
integer(kind=DPI)   :: GIMIX,GIMIY,GIMIZ,GISEGM,GIMODUR(213)
integer(kind=DPI)   :: affnboimagex,affnboimagey,affnboimagez
integer(kind=DPI)   :: KODERR
character           :: fichier*50,fichbase*50,caractest*1

dimension iseg(5,nsegmax),junctt(nsegmax),surface(nsegmax)

!Warning
Print *, " "
Print *, " "
Print *, "                                 MM    M  EEEE   GGGG    AAA    SSSS"
Print *, "                                 M M M M  EE    G       A   A  SS   "
Print *, "      MM MM  I  CCC  RRRR  OOOO  M  M  M  EEEE  G  GGG  AAAAA   SSS "
Print *, "      M M M  I  C    RRR   O  0  M     M  EE    GG  GG  A   A     SS"
Print *, "      M   M  I  CCC  R  R  0000  M     M  EEEE   GGGG   A   A  SSSS "
Print *, " "
Print *, " "
Print *, " ==========================================================================="
Print *, "                                CopyRight"
Print *, " "
Print *, " mM (for microMegas) is an open source program of DD (Dislocation Dynamics) simulation"
Print *, " originally developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE"
Print *, " Copyright (C) 1993, 1996, 2000, 2002, 2004   B. Devincre, L. Kubin, M. Condat,"
Print *, " C. Lemarchand, R. Madec, S. Groh, G. Monnet, J. Durinck, P. Carrez, C. de Sansal,"
Print *, " M-G. Tagorti, S. Queyreau. "
Print *, "  "
Print *, " This is free software, and you are welcome to redistribute it under the"
Print *, " conditions of the GNU General Public License (see the 'Copyright' file in"
Print *, " the program distribution)."
Print *, " ==========================================================================="
Print *, "              "

write (*,*)
write (*,*) "To visualize trajectory films computed with mM type the following infos!"
write (*,*)

!*** On commence par definir les conditions de visualisation de la boite de simulation
write (*,*) "Number of periodic replica used in direction: X  Y  Z  ?"
read (*,*) affnboimagex,affnboimagey,affnboimagez

write (*,*) "Visualisation start at mM step:  ?"
read (*,*) debutaffich

if (debutaffich <= 0)  debutaffich = 1

write (*,*) "Snapshot frequency (periodicity of mM steps visualized):  ?"
read (*,*) freqaffich
if (freqaffich <= 0)   freqaffich = 1

write (*,*) "Pause frequency (periodicity of mM steps we stop visualization):  ?"
read (*,*) freqpause
if (freqpause <= 0)    freqpause = 1

! initializations
kk        = 1

!*** Ouverture de l unite contenant le film
fichier = "film.bin"

open(49,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)

read(49) cristallo    ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
read(49) nsegmaxf     ! The nsegmax value used in this simulation
read(49) avalue       ! The simulation lattice parameter (the simulation scaling)
read(49) nbase        ! Nb of vectors direction used for the material crystal symmetry
read(49) (bveclin(:,I),bvecnor(:,I), I=1,nbase)    ! table of vector direction and displacement for this simulation
read(49) modur(:)     ! Dimension of the 3D simulated volume in a value

if (nsegmax /= nsegmaxf) then
  write(*,*) 'The standard value of nsegmax=100000 was modified in this film !'
  write(*,*) 'The nsegmax value used in this simulation = ',nsegmaxf
  write(*,*) 'You may have to compile again camera.f90 with the correct value !'
  write(*,*) 'Pause !!!'
  read(*,*)
endif

! en cas de probleme avec les redemarages par exemple
if (cristallo > 7 .or. cristallo < 1) then
  write(*,*) 'The file defining the crystal symmetry is unknown'
  write(*,*) 'Cristallo =',cristallo
  write(*,*) 'You can redefined the cristallo type: '
  write(*,*) 'CS=1, CC=2, CFC=3, HC=4, ORT=5, MGO=6, DC=7'
  write(*,*) 'Your choice ?'
  read (*,*) cristallo

  !attention il faut recommencer la lecture du debut
  !    close(50)
  !
  !    open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)
  ! !        open(50,file='../out/'//fichier,FORM='UNFORMATTED',convert="little_endian",STATUS='OLD',ERR=99,IOSTAT=KODERR)
  ! !        open(50,file='../out/'//fichier,FORM='UNFORMATTED',convert="big_endian",STATUS='OLD',ERR=99,IOSTAT=KODERR)
endif

if (cristallo == 1) fichbase = "CS"
if (cristallo == 2) fichbase = "BCC"
if (cristallo == 3) fichbase = "CFC"
if (cristallo == 4) fichbase = "HCP"
if (cristallo == 5) fichbase = "ORT"
if (cristallo == 6) fichbase = "MGO"
if (cristallo == 7) fichbase = "DC"

if (cristallo > 7 .or. cristallo < 1) then
  print*, "cristallo : Unknown value"          ! No way
  stop
endif

print *, "  The crystal symmetry is of type ",fichbase

close(49)

! Initilization are finiched we start the effective reading
open(50,file='../out/'//fichier,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=KODERR)

read(50) cristallo    ! The crystal symmetry (BVD file used) in the simulation at the origin of this film
read(50) nsegmaxf     ! The nsegmax value used in this simulation
read(50) avalue       ! The simulation lattice parameter (the simulation scaling)
read(50) nbase        ! Nb of vectors direction used for the material crystal symmetry
read(50) (bveclin(:,I),bvecnor(:,I), I=1,nbase)    ! table of vector direction and displacement for this simulation
read(50) modur(:)     ! Dimension of the 3D simulated volume in a value

ijk = 0

! We start the loop on simulation steps saved in film.bin
Steps_loop: do while (kk.gt.0)

  ijk = ijk + 1

  ! A new simulation step
  read(50,IOSTAT=KODERR) kk, deltat, NSEGM, XSTRESS, EPSO

  if(KODERR < 0) then
    print*,"  End of the film.bin file"
    stop
  elseif(KODERR > 0) then
    print*,"Reading problem 1",KODERR,ijk,kk, deltat, NSEGM, XSTRESS, EPSO
    cycle Steps_loop
  endif

  if ( (kk >= debutaffich) .and. (modulo(kk, freqaffich)== 0) ) then

      write(*,'(I8, ", Eps:",ES10.2E2,"%, Sigma =", ES10.2E2," MPa, Nseg = ",I7)') KK, EPSO, XSTRESS, nsegm

      GISEGM=NSEGM

      ! GIMODUR definitions
      ! This table is used to transfer data between the c and Fort codes

      ! Dimension of the simulation cell to visualise
      GIMODUR(1)=MODUR(1)*(2*affnboimagex+1)
      GIMODUR(2)=MODUR(2)*(2*affnboimagey+1)
      GIMODUR(3)=MODUR(3)*(2*affnboimagez+1)

      ! The shifting informations, for the moment set to zero
      ! POUR LE MOMENT LES DECALAGES SONT A ZERO, AU PASSAGE IL SERAIT
      ! PEUT ETRE BON D'UTILISER SIX PETITS SCALAIRES INITIALISES UNE BONNE
      ! FOIS POUR TOUTE AU LIEU DE FAIRE DES SUMS A CHAQUE OPERATION DE CLP...
      GIMODUR(4)=0 !shiftxy
      GIMODUR(5)=0 !shiftxz
      GIMODUR(6)=0 !shiftyx
      GIMODUR(7)=0 !shiftyz
      GIMODUR(8)=0 !shiftzx
      GIMODUR(9)=0 !shiftzy

      ! In the following a Base of 8 vectors per slip systems is assumed
      ! but the graphical interface is more generale than this and works
      ! with more or less sofisticated simulation base of discretization
      GIMODUR(10)=8                 ! nb of vecteurs per slip systems
      GIMODUR(11)=NBASE/8           ! nb of slip systems taken into account
      GIMODUR(12)=nbase             ! total number of vectors

      ! The space scalling factor
      GIMODUR(13)=int(0.000001D0/avalue)

      do I = 1,NBASE
        GIMODUR(13+i)=(i-1)/8+1     ! Index of slip systems (needed to the
      enddo                         ! graphical information function

      ! The segments micro structure is loaded
      READ(50,IOSTAT=KODERR) (Iseg(1:5,j),Junctt(j),tabvois(1:2,j),surface(j), J=1,nsegm)

      if(KODERR < 0) then
        print*,"  End of the film.bin file"
        stop
      elseif(KODERR > 0) then
        print*,"Reading Error 2",KODERR,nsegm,j
        cycle Steps_loop
      endif

      do J=1,nsegm
        ! The information in Junctt is reformatted for the C codes
        ! Junctt = 1 means a junction segment
        ! Junctt = 2 means a GD segment
        ! Junctt = 0 means a standard segment
        if (Junctt(j) > 0) junctt(J)= 1
        if (Junctt(j) < 0) junctt(J)= 2
      enddo

  !     Segs_loop: do j=1,NSEGM      !>>>>>>>>>>
  !         read(50,IOSTAT=KODERR) ISEG(1,J),ISEG(2,J),ISEG(3,J),ISEG(4,J),ISEG(5,J),JUNCTT(J)
  !         print*,'j=',j
  !         print*,ISEG(1,J),ISEG(2,J),ISEG(3,J),ISEG(4,J),ISEG(5,J),JUNCTT(J)
  !         if(KODERR /= 0) then
  !            print*,"probleme lecture 2",KODERR
  !            kk=1
  !            read(*,*)
  !            cycle Steps_loop
  !         endif
  ! !        if (JUNCTT(J).eq.0) junctt(J)=nsegmax
  !         if (JUNCTT(J).eq. J) junctt(J) = J -1
  !     enddo Segs_loop            !<<<<<<<<<<

      ! We save the exact shape (3) of the simulation box
      moduR_ref(3)=moduR(3)

      ! The steps we want to visualize
      if (affnboimagex.ne.0 .or. affnboimagey.ne.0 .or. affnboimagez.ne.0) then

        ! Boucles sur les images
        gimix=0; gimiy=0; gimiz=0
        do while(gimix.ne.(2*affnboimagex).or.gimiy.ne.(2*affnboimagey).or.gimiz.ne.(2*affnboimagez))
            if(gimix.ne.(2*affnboimagex)) then
              gimix=gimix+1
            else
              if(gimiy.ne.(2*affnboimagey)) then
                  gimix=0 ; gimiy=gimiy+1
              else
                  if(gimiz.ne.(2*affnboimagez)) then
                    gimix=0; gimiy=0;  gimiz=gimiz+1
                  endif
              endif
            endif
            JUNCTT(gISEGM+1:gISEGM+NSEGM) = JUNCTT(1:nsegm)
            ISEG(1,gISEGM+1:gISEGM+NSEGM) = ISEG(1,1:NSEGM)+gimix*modur(1)
            ISEG(2,gISEGM+1:gISEGM+NSEGM) = ISEG(2,1:NSEGM)+gimiy*modur(2)
            ISEG(3,gISEGM+1:gISEGM+NSEGM) = ISEG(3,1:NSEGM)+gimiz*modur(3)
            ISEG(4,gISEGM+1:gISEGM+NSEGM) = iSEG(4,1:nsegm)
            ISEG(5,gISEGM+1:gISEGM+NSEGM) = iSEG(5,1:nsegm)
            gISEGM=gISEGM+NSEGM
        enddo
      endif

      call PREDRAW(gisegm,iseg,TABVOIS,BVECLIN,kk,gimoduR,junctt,nbsegdep,nsegmax)

      !The next test is to take user input into account
      if(modulo(kk, freqpause) == 0 ) then

        print *," Pause : press any key to continue ..."
        read *
        call PREDRAW(gisegm,iseg,TABVOIS,BVECLIN,kk,gimoduR,junctt,nbsegdep,nsegmax)

      else

        call PREDRAW(gisegm,iseg,TABVOIS,BVECLIN,kk,gimoduR,junctt,nbsegdep,nsegmax)

      endif

      read(50,IOSTAT=KODERR) caractest

      if(KODERR < 0) then
        print*,"  End of the film.bin file"
        stop
      elseif(KODERR > 0) then
        print*,"Reading error 3",caractest
        cycle Steps_loop
      endif

  else

    READ(50,IOSTAT=KODERR)

    caractest = '@'
    do while (caractest /= '&')
      read(50,IOSTAT=KODERR) caractest

      if(KODERR < 0) then
        print*,"  End of the film.bin file"
        stop
      elseif(KODERR > 0) then
        print*,"Reading problem 666",KODERR
        stop
      endif
    enddo

  endif

enddo Steps_loop

close (50)

END SUBROUTINE MICROMEGAS
