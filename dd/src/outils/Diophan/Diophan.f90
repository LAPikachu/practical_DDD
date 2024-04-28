module bricamat

contains

!###########################################################################
!# Calcul du carre de la norme d'un vecteur entier                         #
!############################################################### 06/11/98 ##
function inorivect(i)

implicit none

integer (8):: i(3),inorivect

inorivect=i(1)*i(1)+i(2)*i(2)+i(3)*i(3)
return

end function inorivect

!###########################################################################
!# Calcul du produit vectoriel a^b  (a et b sont des vecteurs entier)      #
!############################################################### 06/11/98 ##
function prodivec(a,b)

implicit none

integer(8) :: prodivec(3),a(3),b(3)

prodivec(1) = a(2)*b(3)-a(3)*b(2)
prodivec(2) = a(3)*b(1)-a(1)*b(3)
prodivec(3) = a(1)*b(2)-a(2)*b(1)

return
end function prodivec

!###########################################################################
!# Produit scalaire entre deux vecteurs entiers a.b                        #
!############################################################### 06/11/98 ##
function iprodsca(a,b)

implicit none

integer(8) :: iprodsca
integer(8)  :: a(3),b(3)

iprodsca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

return
end function iprodsca

!###########################################################################
!# decrit la relation entre deux vecteurs x et y : si 
!    3 = egalite 
!    2 = parallels , meme sens (normes differentes)
!    1 = angle > 0 et < 90 degree
!    0 = perpendiculaires
!   -1 = angle > -90 et < 0
!   -2 = parallels et sense opposee (normes differentes)
!   -3 = oppose (somme = 0)
!############################################################### 01/02/01 ##

function etat(x,y)

implicit none

integer(kind=8),intent(in),dimension(3) :: x,y
integer(kind=8),dimension(3) :: somme,bufprodvec
integer(kind=8) :: etat,bufprodsca

etat = -10
somme(:) = x(:) - y(:)
if (inorivect(somme) == 0) then ! vecteurs identiques
   etat = 3
   return
else
   somme(:) = x(:) + y(:)
   if (inorivect(somme) == 0) then ! vecteurs opposes
      etat = -3
      return
   else
      bufprodsca=iprodsca(x,y)
      bufprodvec=prodivec(x,y)
      if (inorivect(bufprodvec) == 0) then  ! vecteurs parallels
         if (bufprodsca > 0) etat = 2   ! meme sens
         if (bufprodsca < 0) etat = -2  ! sens opposes
         return
      else
         if (bufprodsca == 0)  etat = 0 ! perpendiculaire
         if (bufprodsca > 0)   etat = 1 ! meme sens (angles obtus)
         if (bufprodsca < 0)   etat = -1! sens opposes (angles aigus)
         return
      endif
   endif
endif
!if (abs(etat) > 3) STOP ' erreur dans programme etat(x,y)'
return

end function etat

end module bricamat

!############################################################### 01/02/01 ##

program Diophan

use bricamat

implicit none
logical :: ok
integer(kind=4) :: ix,iy,iz,i
integer(kind=4) :: ix1,iy1,iz1
integer(kind=4) :: ix2,iy2,iz2
integer(kind=4) :: ix3,iy3,iz3
integer(kind=4) :: ix4,iy4,iz4
integer(kind=4) :: m,mm
integer(kind=8) :: h,k,l
integer(kind=8) :: Lx,Ly,Lz,dxy,dxz,dyx,dyz,dzx,dzy
integer(kind=8) :: bLx,bLy,bLz,bdxy,bdxz,bdyx,bdyz,bdzx,bdzy
integer(kind=8) :: iniLx,iniLy,iniLz,inidxy,inidxz,inidyx,inidyz,inidzx,inidzy
integer(kind=8) :: finLx,finLy,finLz,findxy,findxz,findyx,findyz,findzx,findzy
integer(kind=8) :: lprodsca,combi,combiv(100000,3),vect(3)
real(kind=8)    :: dd,d1,d2,d3,d4,dL,combid(100000)
real(kind=8)    :: dist,deltadist,tdeltadist,S,Ls

!*** Initialisation bricamat
!tabprem(:)=0;tabprem(1)=1;tabprem(2)=2;

combi=0
dist=0d0;deltadist=900000000000000000d0
bLx=0;bLy=0;bLz=0
bdxy=0;bdxz=0;bdyx=0;bdyz=0;bdzx=0;bdzy=0

!****** LECTURE *****

write(*,*)
write(*,*) 'Rappel des conditions de calculs !'
write(*,*)

open(1,file='Diophan.IN',STATUS='OLD') 
read (1,*) h,k,l
write(*,'(I4," h  ",I4," k  ",I4," l  ")') h,k,l
read (1,*) iniLx,iniLy,iniLz
write(*,'(I4," Lx_min   ",I4," Ly_min   ",I4," Lz_min   ")') iniLx,iniLy,iniLz
read (1,*) finLx,finLy,finLz
write(*,'(I4," Lxmax    ",I4," Ly_max   ",I4," Lz_max   ")') finLx,finLy,finLz
read (1,*) inidxy,inidxz
write(*,'(I4," dxy_min   ",I4," dxz_min   ")') inidxy,inidxz
read (1,*) findxy,findxz
write(*,'(I4," dxy_max   ",I4," dxz_max   ")') findxy,findxz
read (1,*) inidyx,inidyz
write(*,'(I4," dyx _min  ",I4," dyz _min  ")') inidyx,inidyz
read (1,*) findyx,findyz
write(*,'(I4," dyx_max   ",I4," dyz_max   ")') findyx,findyz
read (1,*) inidzx,inidzy
write(*,'(I4," dzx _min  ",I4," dzy _min  ")') inidzx,inidzy
read (1,*) findzx,findzy
write(*,'(I4," dzx_max   ",I4," dzy_max   ")') findzx,findzy
write(*,*)
read (1,*) mm
write(*,*) mm," nb iterations par entiers"
write(*,*)

write(*,*)
write(*,*) 'Debut des calculs de calculs !'
write(*,*)

!*** Evaluation de la longueur en isotrope

S=finlx*finly*finlz	! Pgcd supposé unitaire
Ls=sqrt(S)/(((finLx)**2+(finLy)**2+(finLz)**2)**(0.5d0))
! write (*,*) "L iso",Ls   ! Longueur des ref un peu mysterieuse
! write(*,*)

!*****************************************

do Lx=iniLx,finLx,1;do Ly=iniLy,finLy,1;do Lz=iniLz,finLz,1
do dxy=inidxy,findxy,1;do dxz=inidxz,findxz,1
do dyx=inidyx,findyx,1;do dyz=inidyz,findyz,1
do dzx=inidzx,findzx,1;do dzy=inidzy,findzy,1

if((dxy.ne.0.and.dyx.ne.0).or.(dxz.ne.0.and.dzx.ne.0).or.(dyz.ne.0.and.dzy.ne.0)) then
 write (*,*) "shifts incompatibles"
 cycle
endif

d1=900000000000000000d0;d2=d1;d3=d2;d4=d3
ix1=0;iy1=0;iz1=0
ix2=0;iy2=0;iz2=0
ix3=0;iy3=0;iz3=0
ix4=0;iy4=0;iz4=0
dL=((Lx)**2+(Ly)**2+(Lz)**2)**(0.5d0)
m=mm
ok=.false.

do ix=-m,m,1;do iy=-m,m,1;do iz=-m,m,1
 if(ix.eq.0.and.iy.eq.0.and.iz.eq.0)cycle
 lprodsca=h*(Lx*ix+iy*dxy+iz*dxz)+k*(Ly*iy+ix*dyx+iz*dyz)+l*(Lz*iz+ix*dzx+iy*dzy)
 if(lprodsca.eq.0) then
  dd=((Lx*ix+iy*dxy+iz*dxz*1.d0)**2+(Ly*iy+ix*dyx+iy*dyz*1.d0)**2+&
      (Lz*iz+ix*dzx+iy*dzy*1.d0)**2)**(0.5d0)*0.5D0/dL


vect(1)=(Lx*ix+iy*dxy+iz*dxz)
vect(2)=(Ly*iy+ix*dyx+iz*dyz)
vect(3)=(Lz*iz+ix*dzx+iy*dzy)

do i=1,combi,1
 if(abs(etat(vect(1:3),combiv(i,1:3))).ge.2)then
  if(dd.lt.combid(i))then
   combiv(i,1:3)=vect(1:3)
   combid(i)=dd
!write (*,*) dd
!write (*,*) vect(1:3)
!write (*,*) ix,iy,iz
  endif
  goto 36
 endif
enddo      

combi=combi+1
combiv(combi,1:3)=vect(1:3)
combid(combi)=dd
!write (*,*) dd
!write (*,*) vect(1:3)
!write (*,*) ix,iy,iz
      
36 if((ix+ix1).eq.0.and.(iy+iy1).eq.0.and.(iz+iz1).eq.0)cycle
 if((ix+ix2).eq.0.and.(iy+iy2).eq.0.and.(iz+iz2).eq.0)cycle
 if((ix+ix3).eq.0.and.(iy+iy3).eq.0.and.(iz+iz3).eq.0)cycle
 if((ix+ix4).eq.0.and.(iy+iy4).eq.0.and.(iz+iz4).eq.0)cycle
  if(dd.lt.d1) then
   d4=d3 ; ix4=ix3 ; iy4=iy3 ; iz4=iz3
   d3=d2 ; ix3=ix2 ; iy3=iy2 ; iz3=iz2
   d2=d1 ; ix2=ix1 ; iy2=iy1 ; iz2=iz1
   d1=dd ; ix1=ix  ; iy1=iy  ; iz1=iz
!   write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix,iy,iz,d1,d2,d3,d4
   ok=.true.
  else
   if(dd.lt.d2) then
    d4=d3 ; ix4=ix3 ; iy4=iy3 ; iz4=iz3
    d3=d2 ; ix3=ix2 ; iy3=iy2 ; iz3=iz2
    d2=dd ; ix2=ix  ; iy2=iy  ; iz2=iz
!    write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix,iy,iz,d1,d2,d3,d4
    ok=.true.
   else
    if(dd.lt.d3) then
     d4=d3 ; ix4=ix3 ; iy4=iy3 ; iz4=iz3
     d3=dd ; ix3=ix  ; iy3=iy  ; iz3=iz 
!     write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix,iy,iz,d1,d2,d3,d4
     ok=.true.
    else
     if(dd.lt.d4) then
      d4=dd ; ix4=ix ; iy4=iy ; iz4=iz
!      write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix,iy,iz,d1,d2,d3,d4
      ok=.true.
     endif
    endif
   endif
  endif
 endif
enddo;enddo;enddo

if(.not.ok)stop "Nombre d'itérations par entiers trop faible !!!"

!write (*,'(5(I5,1X)," => :")') Lx,Ly,lz,dyx,dyz
!write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix4,iy4,iz4,d4
!write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix3,iy3,iz3,d3
!write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix2,iy2,iz2,d2
!write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix1,iy1,iz1,d1

tdeltadist=(d2-d1)/d1

write(*,*) ' Progression des resultats :'
write (*,'(9(I4,1X),">",2(F10.2,1X))') Lx,Ly,Lz,dxy,dxz,dyx,dyz,dzx,dzy,d1,tdeltadist
write (*,'(9(I4,1X),">",2(F10.2,1X))') bLx,bLy,bLz,bdxy,bdxz,bdyx,bdyz,bdzx,bdzy,dist,deltadist
write(*,*) 

if((d1.gt.dist.and.tdeltadist.lt.deltadist).or.(d1.gt.dist.and.d1.lt.Ls)) then
 dist=d1;deltadist=tdeltadist
 bLx=Lx;bLy=Ly;bLz=Lz;bdxy=dxy;bdxz=dxz;bdyx=dyx;bdyz=dyz;bdzx=dzx;bdzy=dzy
 write (*,*) "Quatres plus petits libre parcourt",d1,d2,d3,d4
 write (*,*) "Plus petite distance et facteur d anisotropie ===>",dist,deltadist
 write(*,*) 
! pause
endif

!*** CC
enddo;enddo;enddo;enddo;enddo;enddo;enddo;enddo;enddo

write(*,*) '###################'
write(*,*) '# Resultat finale #'
write(*,*) '###################'
write(*,*)
write(*,'(I4," h  ",I4," k  ",I4," l  ")') h,k,l
write(*,*)
write(*,'(I4," Lx_min   ",I4," Ly_min   ",I4," Lz_min   ")') iniLx,iniLy,iniLz
write(*,'(I4," Lxmax    ",I4," Ly_max   ",I4," Lz_max   ")') finLx,finLy,finLz
write(*,'(I4," dxy_min   ",I4," dxz_min   ")') inidxy,inidxz
write(*,'(I4," dxy_max   ",I4," dxz_max   ")') findxy,findxz
write(*,'(I4," dyx _min  ",I4," dyz _min  ")') inidyx,inidyz
write(*,'(I4," dyx_max   ",I4," dyz_max   ")') findyx,findyz
write(*,'(I4," dzx _min  ",I4," dzy _min  ")') inidzx,inidzy
write(*,'(I4," dzx_max   ",I4," dzy_max   ")') findzx,findzy
write(*,*)
write(*,*) mm," nb iterations par entiers"
!write(*,*)
!write (*,*) "L iso",Ls
write(*,*)
write(*,*) 'Dimension de la boite optimum  ==> ',bLx,bLy,bLz
write(*,*) 'Solution de shift optimum      ==> ',bdxy,bdxz,bdyx,bdyz,bdzx,bdzy
write(*,*) 'Libre parcourt min et anisotropie  =====> ',dist,deltadist
write(*,*)

!write (*,'(3(I3,1X),"->",4(F8.3,1X))') ix,iy,iz
write(*,*) 'Solutions des equations Diophantines pour les quatres plus petites distances de libre parcourt'
write (*,'(3(I3,1X))') ix1,iy1,iz1
write (*,'(3(I3,1X))') ix2,iy2,iz2
write (*,'(3(I3,1X))') ix3,iy3,iz3
write (*,'(3(I3,1X))') ix4,iy4,iz4

! Le fichier ci dessous donne la liste des vecteurs entiers definissant une
! ligne boucle sur elle meme - Moebus line !!!
open(2,file='LignesInfinies.OUT',STATUS='REPLACE') 
write(2,*) combi
do i=1,combi,1
 write(2,'(f10.4,3(1x,I12))') combid(i),combiv(i,1:3)
enddo
write(2,*) 
write(2,'(I4," h",I4," k",I4," l")') h,k,l
write(2,'(I8," iLx",I8," iLy",I8," iLz")') iniLx,iniLy,iniLz
write(2,'(I8," idxy",I8," idxz")') inidxy,inidxz
write(2,'(I8," idyx",I8," idyz")') inidyx,inidyz
write(2,'(I8," idzx",I8," idzy")') inidzx,inidzy
write(2,*) mm," nb iterations par entiers"
close(2) 

write(*,*)

end program Diophan
