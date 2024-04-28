! programme permetant de tester simplement les solution de l'equation de Diaphantine
! qui controle une partie des artefactes des CLP de la simulation
program main

implicit none

integer  :: Nx , Ny , Nz   ! dimension de la boite
integer,parameter  :: nxx=2 , nyy=-1 , nzz=-1    ! direction du plan de glissement
integer  :: i,j,k,resu
integer  :: vecx,vecy
real     :: d1,d2




do Nx= 280, 320,2
    do Ny = 280, 320,2
        do Nz = 280, 320,2

            d1=1000000

            do i= -50, 50
                do j= -50, 50
                    do k= -50, 50
                        if(i == 0 .and. j == 0 .and. k == 0) cycle
                        resu = i * Nx * nxx + j * Ny * nyy + k * Nz * nzz
                        if (resu == 0) then
                           d2=sqrt(real(i*Nx)*real(i*Nx)+real(j*Ny)*real(j*Ny)+real(k*Nz)*real(k*Nz))*0.5
                           if (d2 < d1) then
                              d1=d2
                              if(abs(i-j) < i/2 .and. abs(i-k) < i/2 .and.abs(k-j) < i/2)then
                                 print *, "--------------------------------------------------------"
     write (*,'(" Nx, Ny, Nz = ",3I8," UVW = ",3I5," lambda = ", F10.6)') &
                                Nx*1584,Ny*1584,Nz*1584,i,j,k,d2/sqrt(real(Nx*Nx)+real(Ny*Ny)+real(Nz*Nz))
                              endif
                           endif
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo

END program main

! 
! 
!
!
!
!
!
!
!
!

