!***************************************************************
!Calcul du produit scalaire de deux vecteurs normalise
function absprodscanorm(v1,v2)

implicit none
real*8 :: v1(3),v2(3),absprodscanorm

absprodscanorm=abs((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/&
(sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))*&
 sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))))

!write (*,*) absprodscanorm

end function absprodscanorm

!***************************************************************
!Calcul des facteur de Schmidt
function FS(z,b,n)

implicit none
real*8 :: z(3),b(3),n(3),FS,absprodscanorm


!write (*,*) z
!write (*,*) b
!write (*,*) n

FS=absprodscanorm(z,b)*absprodscanorm(z,n)

end function FS

!***************************************************************
!Produit scalaire de deux vecteur
function absprodsca(v1,v2)

implicit none
integer :: v1(3),v2(3),absprodsca

absprodsca=abs(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))

end function absprodsca

!***************************************************************
!Reconstitution des courbes a partir des coeficient d interactions
program MADHOC

implicit none

real*8 :: z(3),fs,tfs(12),coeffa(6),TAUpredSYS(12),AlphaSYSexp(12),&
          AlphaSYS(12),RauSYS(12),burg,b(3),n(3),hocus

real*8 :: accutime,RAUDIS,RAUDMO,RAUDJONC,VITMOY,EPSDOT,EPSO,TauFric,mu
real*8 :: XSTRESS,SQ_RAUDIS,EPSPLA(3,3),XSTRESS_SQ_RAUDIS  
real*8 :: RAUSYS_VOLUME_AVALUE_AVALUE(12)
real*8 :: GAMSYS(12)
real*8 :: TauCoeff(6)

integer :: i,j,matinter(12,12),absprodsca,tb(12,3),tn(12,3)

integer :: JK,toto

!***********
!Lecture des donnŽes de calcul
open(8,FILE='madhocSAR.in',STATUS='OLD',POSITION='rewind',action='read')

! axe de solicitation
read(8,*)z(1),z(2),z(3)
write (*,*) ' Z : ',z(1),z(2),z(3)

! info sur les systemes
do i=1,12
 read(8,*) tb(i,1),tb(i,2),tb(i,3),tn(i,1),tn(i,2),tn(i,3)
 b(1:3)=dfloat(tb(i,1:3))
 n(1:3)=dfloat(tn(i,1:3))
 hocus=fs(z,b,n)  ! calcul du facteur de schmidt
 tfs(i)=hocus
enddo

!lecture des coefs
do i=1,6
 read(8,*) coeffa(i)
 write (*,*) 'Coeff a',i-1,' -> ',coeffa(i)
enddo

!norme du vecteur de burgers
read(8,*) burg
write (*,*) 'Burgers : ',burg

!contrainte de friction de reseau
read(8,*) TauFric
write (*,*) 'TauFric : ',TauFric

!module de cisaillement
read(8,*) mu
write (*,*) 'mu : ',mu
close(8)

!creation de la matrice des a
matinter(:,:)=0

do i=1,12
 do j=1,12

  if (i.eq.j) then
    !  write(*,*) i,j,'identiques',matinter(i,j)
  else if (absprodsca(tn(i,:),tn(j,:)).eq.3) then
    matinter(i,j)=1
    ! write(*,*) i,j,'coplanaires',matinter(i,j)
  else if (absprodsca(tb(i,:),tb(j,:)).eq.0) then
    matinter(i,j)=3
    !  write(*,*) i,j,'hirth',matinter(i,j)
  else if (absprodsca(tb(i,:),tb(j,:)).eq.2) then
    matinter(i,j)=2
    !  write(*,*) i,j,'devie',matinter(i,j)
  else if ((absprodsca(tb(i,:)+tb(j,:),tn(i,:)).eq.0).or.&
          (absprodsca(tb(i,:)+tb(j,:),tn(j,:)).eq.0)) then
    matinter(i,j)=4
    !  write(*,*) i,j,'glissile',matinter(i,j)
  else !GLISSILE
    matinter(i,j)=5
    !  write(*,*) i,j,'lomer',matinter(i,j)

  endif
 
 enddo
enddo

!ecriture de la matrice
do i=1,12
 write (*,59) i,matinter(i,:)
enddo

!***********
open(18,FILE='out/sigeps',action='read')
open(28,FILE='out/rau',action='read')
open(38,FILE='out/gamma',action='read')

open(17,FILE='SigSys.dat',STATUS='REPLACE')
close(17)

toto=1

do while(toto.eq.1)

  !lecture de sigeps
  read(18,*) XSTRESS,EPSO,SQ_RAUDIS,EPSPLA(1,1),EPSPLA(2,2), &
              EPSPLA(3,3),EPSPLA(1,2),EPSPLA(1,3),EPSPLA(2,3),&
              XSTRESS_SQ_RAUDIS,accutime  

  !lecture de rau
  read(28,*) EPSO,(RAUSYS_VOLUME_AVALUE_AVALUE(JK),JK=1,12)
  
  !lecture de gamma
  read(38,*) EPSO,(GAMSYS(JK),JK=1,12)

  do i=1,12
    TAUpredSYS(i)=0.
    RauSYS(i)=0.
    if(tfs(i).ne.0) then
      do j=1,12
        TAUpredSYS(i)=TAUpredSYS(i)+coeffa(matinter(i,j)+1)*RAUSYS_VOLUME_AVALUE_AVALUE(j)
      enddo
      TAUpredSYS(i)=mu*burg/tfs(i)*(sqrt(TAUpredSYS(i))+TauFric)
    endif
  enddo

  open(17,FILE='SigSys.dat',STATUS='OLD',position='APPEND')
    write (17,55) accutime,epso,xstress,(TAUpredSYS(JK),JK=1,12)
  close(17)

enddo

26 format(2X,7(E12.5,1X))
27 format(2X,11(E12.5,1X))
58 format(2X,8(E12.5,1X))
37 format(2X,13(E12.5,1X))
55 format(2X,15(E12.5,1X))
57 format(2X,20(E12.5,1X))
56 format(100(I4,1X),E12.5)
59 format(20(I2,1X))

close( 8)
close(18)
close(19)
close(28)
close(38)
close(25)
close(48)
close(49)

end program MADHOC
