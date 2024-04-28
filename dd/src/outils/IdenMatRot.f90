program IdenMatRot

implicit none

integer,parameter  :: DP=selected_real_kind(p=14) !< Size for the variables with real type
integer,parameter  :: DPI=selected_int_kind(13)   !< Size for the variables with integer type

integer(kind=DPI) :: limit
integer(kind=DPI) :: x,y,z,xtmp,IUN,xtest
integer(kind=DPI) :: x2,y2,z2
integer(kind=DPI) :: nresul

real(kind=DP) :: error_vref,error_vref2,quality
real(kind=DP) ::  nv1, nv2, nv3

real(kind=DP), dimension(3) :: vref   ,nvref
real(kind=DP), dimension(3) :: vref2  ,nvref2
real(kind=DP), dimension(3) :: v1     ,v1n
real(kind=DP), dimension(3) :: v2     ,v2n
real(kind=DP), dimension(3) :: v3     ,v3n

!Initialization
limit = 50
error_vref = 0.999
error_vref2 = 0.5

IUN = 1


! The reference vector we want to approach as [100]:
vref =  (/-1, 2, 1 /)
! The reference vector we want to approach as [010]:
vref2 =  (/ -1, 0, -1 /)

xtest = int(vref(1))

nvref(:) = vref(:)/norvect(vref)
nvref2(:) = vref2(:)/norvect(vref2)


x=-limit

nresul=0

!Do  x = -limit,limit
!  if (sign(1.,float(x)) /= sign(1d0,vref(1)) ) cycle

Do  xtmp = 0,limit
  x = sign(IUN,xtest)*xtmp
  write(*,fmt='("Step ",I4,"/",I4)') xtmp,limit
  Do  y = -limit,limit

    Do  z = -limit,limit

      v1 =  (/ x, y, z /)
      nv1 = norvect(v1)

      if (nv1 <= 0) cycle
      v1n(:) = v1(:) /nv1
      !print *,x,y,z,dot_product(v1n,nvref), dot_product(v1n,nvref) > error, (nv1 == int(nv1))
      if ( (dot_product(v1n,nvref) >= error_vref) .and. (nv1 == int(nv1)) ) then

        !Loops to find a vector v2 compatible with v1

        Do  x2 = -limit,limit
          Do  y2 = -limit,limit
            Do  z2 = -limit,limit

                !the vector v2
                v2= (/ x2, y2, z2 /)
                nv2 = norvect(v2)
                if (nv2 <= 0) cycle
                v2n(:) = v2(:) /nv2
                !print *,int(v1),int(v2)
                if ( (dot_product(v2n,nvref2) >= error_vref2) .and. (nv1 == nv2) ) then

                  v3 = cross(v1,v2)/nv2
                  nv3 = norvect(v3)

                  if (nv3 <= 0) cycle
                  v3n(:) = v3(:) /nv3
                  if (nv1 == nv3) then
                      nresul=nresul+1
                      quality = (sqrt(dot_product(v1n,nvref)**2 + dot_product(v2n,nvref2)**2) -sqrt((error_vref+error_vref2)**2))*10000
                      write( *,fmt='(I4, ": [",3I4,"] ","[",3I4,"] "," [",3I4,"] "," Norm: ",I4," Quality: ", 2F11.6, F9.3)') &
nresul,int(v1(:)), int(v2(:)),int(v3(:)), int(nv3),dot_product(v1n,nvref),dot_product(v2n,nvref2),quality
                  endif
                endif

            enddo !x2
          enddo !y2
        enddo !z2


      endif
    enddo !z

  enddo !y

enddo !x

if (nresul ==0) then
  write(*,*) "!> Nothing was found. You could maybe increase limit or decrease accepted error value"
else
  write(*,*) "!> Results found : ", nresul
endif

contains

FUNCTION cross(a, b)
  REAL(kind=DP), DIMENSION(3) :: cross
  REAL(kind=DP), DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

function norvect(r)

implicit none

real(kind=DP) :: norvect,r(3)

norvect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

end function norvect

end program IdenMatRot

