!***************************************************************
!produit scalaire de deux vecteur nomalise
function absprodscanorm(v1,v2)

implicit none
real(kind=8) :: v1(3),v2(3),absprodscanorm

absprodscanorm=abs((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/&
(sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))*&
 sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))))

!write (*,*) absprodscanorm

end function absprodscanorm

!***************************************************************
!calcul du facteur de Schmidt
function FS(z,b,n)

implicit none
real(kind=8) :: z(3),b(3),n(3),FS,absprodscanorm


!write (*,*) z
!write (*,*) b
!write (*,*) n

FS=absprodscanorm(z,b)*absprodscanorm(z,n)

end function FS

!***************************************************************
!Produit scalaire de deux vecteurs
function absprodsca(v1,v2)

implicit none
integer :: v1(3),v2(3),absprodsca

absprodsca=abs(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))

end function absprodsca

!***************************************************************
!determination des coeficient de la matrice d interaction
program depou_coef

implicit none

real(kind=8) :: accutime,EPSO,z(3),VecBurgers,tau0,densite_foret
real(kind=8) :: densite_ref,densite_foret_0,densite_ref_0
real(kind=8) :: XSTRESS,SQ_RAUDIS,EPSPLA(3,3),XSTRESS_SQ_RAUDIS
real(kind=8) :: rausys(12),sysfor(12),gammasys(12),raujonc(12)
real(kind=8) :: fs,tfs(12),b(3),n(3),hocus,mu,coef,tauresnorm
real(kind=8) :: rausysm1,gammasysm1,pente,coK,klpm

integer :: i,j,toto,absprodsca,tb(12,3),tn(12,3),sysmob,kk
integer,parameter  :: DP=8

open(8,FILE='in_depou')
open(18,file='../out/sigeps')
open(28,FILE='../out/rau',action='read')
open(38,FILE='../out/gamma',action='read')
open(48,FILE='../out/raujonc',action='read')
open(19,file='coef.dat')
open(20,file='rauga.dat')

!lecture de l axe de traction
read(8,*)z(1),z(2),z(3)
!Norme du vecteur de Burger
read(8,*)VecBurgers
!Friction de reseau
read(8,*)tau0
!module de cisaillement
read(8,*)mu
!identite du systeme mobile
read(8,*)sysmob
!composition de la foret (zero ou un suivant qu un systeme sois present ou absant)
read(8,*)(sysfor(j),J=1,12)

!infos sur les systemes
do i=1,12
 read(8,*) tb(i,1),tb(i,2),tb(i,3),tn(i,1),tn(i,2),tn(i,3)
 b(1:3)=real(tb(i,1:3),DP)
 n(1:3)=real(tn(i,1:3),DP)
 hocus=fs(z,b,n)
 tfs(i)=hocus   ! facteur de Schmidt
 write (*,*) i,hocus
enddo

!boucle de calculs
toto=1
kk=1
klpm=0

do while(toto.eq.1)

  read(18,*) XSTRESS,EPSO,SQ_RAUDIS,EPSPLA(1,1),EPSPLA(2,2), &
              EPSPLA(3,3),EPSPLA(1,2),EPSPLA(1,3),EPSPLA(2,3),&
              XSTRESS_SQ_RAUDIS,accutime
  read(28,*) epso,(rausys(j),J=1,12)
  read(38,*) epso,(gammasys(j),J=1,12) 
  read(48,*) epso,(raujonc(j),J=1,12) 

  densite_foret=0.
  densite_ref=0.

  !calcul de la densite de foret et totale
  ! rho_forest_sysmob = som_sys((rho_i-rho_jonc_i)*(1 si syst secant ou 0 si syst copla)
  do i=1,12
    densite_foret=densite_foret+(rausys(i)-raujonc(i))*sysfor(i)
    !densite_foret=densite_foret+rausys(i)*sysfor(i)
  enddo

  !Evolution of the reference density
  !densite_ref=densite_ref+(rausys(sysmob)-raujonc(sysmob))
  densite_ref=densite_ref+rausys(sysmob)

  if (kk == 50) then
    densite_foret_0=densite_foret
    densite_ref_0=densite_ref
  endif

  !calcul de la contrainte resolue
  tauresnorm=((1.E6*Xstress*tfs(sysmob)-tau0)/(mu*VecBurgers))

  !calcul du coeficient de la matrice d interaction
  !coef = ((tau-tau_0)/(mu*b))**2/rho_forest
  coef=tauresnorm*tauresnorm/densite_foret

  !calcul du lpm
  !K_lpm =
  !sqrt(coef)*gamma_i/(2*b*(sqrt(rho_forest_sysmob)-sqrt(rho_forest_sysmod_0)))
  if (kk /= 50) then 
          
    !Formule manipe en masse
!    klpm=sqrt(coef)*gammasys(sysmob)/(2.*VecBurgers*      &
!        (sqrt(densite_foret)-sqrt(densite_foret_0)))
    !Formule manipe modele
    klpm=sqrt(coef)*gammasys(sysmob)*sqrt(densite_foret)/  &
           (VecBurgers*(densite_ref-densite_ref_0))
  endif

  !sauvegarde
  write(19,*) accutime,sqrt(coef),klpm
  write(20,*) gammasys(sysmob),rausys(sysmob),densite_foret

  kk=kk+1

enddo

close(18)
close(19)
close(20)
close(8)
close(28)
close(38)

end
