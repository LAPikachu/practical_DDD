! Programme pour tester les solutions de l'equation de Diaphantine
! qui controle les artefactes des CLP dans la simulation DD
program diaphantine

implicit none

integer,parameter  :: Lx=4560 , Ly=7296 , Lz=2736    ! dimension de la boite
integer,parameter  :: nx=1    , ny=1    , nz=1           ! direction du plan de glissement

real,dimension(1000,4)  :: solution
real,dimension(2,5)     :: projection,better

real     :: borne,norme,memo,dnorme
real     :: nnx,nny,nnz,fact,grandcote,cosphi
integer  :: kx,ky,kz,i,resu,maxx,j,sol,solbino

nnx=nx; nny=ny; nnz=nz
better(1,4)=10000000.0   ! On initialise la norme a une longueur giga

dnorme=sqrt(real(nnx**2+nny**2+nnz**2))
nnx=nnx/dnorme
nny=nny/dnorme
nnz=nnz/dnorme

! calcul de l'indice de plan n maximum dans la boite de simulation
maxx=abs(Lx*nx+Ly*ny+Lz*nz)

! longeur max des solutions recherchee
write(*,*) 'nombre d image maxi prise en compte'
read(*,*) fact
borne=fact*real(Lx*nnx+Ly*nny+Lz*nnz)
grandcote=max(Lx,Ly,Lz)

do i=0,maxx
print*,i,maxx

j=0 ! le nombre de solutions par plan est initialise a zero

! On cherche toutes les solutions proches
do kx=-50,50
  do ky=-50,50
    do kz=-50,50
       resu=kx*Lx*nx+ky*Ly*ny+kz*Lz*nz    ! Equation Diaphantine
       if (resu == i) then 
         norme=sqrt((real(kx*Lx))**2+(real(ky*Ly))**2+(real(kz*Lz))**2)  ! norme du vecteur mod
         if ( norme < borne ) then

           j=j+1
           if (j>1000) stop            !teste de debordement de tableau
           solution(j,1)=real(kx*Lx)   !on note toutes les solution OK
           solution(j,2)=real(ky*Ly)
           solution(j,3)=real(kz*Lz)
           solution(j,4)=norme

         endif
       endif
    enddo
  enddo
enddo

! Parmi les solutions on cherche le plus petit binome qui va bien
do sol=1,j-1
  do solbino=sol+1,j

    if (abs(solution(sol,4)-solution(solbino,4)) <  grandcote ) then  ! les deux solutions on presque meme norme

      cosphi=(solution(sol,1)*solution(solbino,1)+  &
              solution(sol,2)*solution(solbino,2)+  &
              solution(sol,3)*solution(solbino,3))/(solution(sol,4)*solution(solbino,4))

       if (cosphi < -0.95) then    ! les deux vecteur mod sont donc presque antiparallele

          if(solution(sol,4) < better(1,4)) then
            better(1,1:4)=solution(sol,1:4)
            better(2,1:4)=solution(solbino,1:4)
            better(1,5)=i
            better(2,5)=i
            print*,better(1,1:5)
            print*,better(2,1:5)
          endif

        endif

!      endif
    endif

  enddo
enddo

enddo

write(*,*) 'Le binome de vecteur mod donnant la plus petite solution'
write(*,*) better(1,1:5),better(1,4)/maxx
write(*,*) better(2,1:5)
write(*,*) better(1,1)/Lx,better(1,2)/Ly,better(1,3)/Lz
write(*,*) better(2,1)/Lx,better(2,2)/Ly,better(2,3)/Lz

end
