
!===================================================================================================
!========================    DEBUT    MODULE  "BRICAMAT"   =========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module does not have a brief description yet
!! \todo  Change some comments into doxygen comments
!! \todo  Translate some comments into English
module BRICAMAT

use constantes
use varglob
#ifdef MDC
use varglob
#endif

implicit none


INTEGER :: seed1  ,&  !< Generateur uniforme de nombre pseudo-aleatoire
           seed2  ,&  !< Generateur uniforme de nombre pseudo-aleatoire
           seed3      !< Generateur uniforme de nombre pseudo-aleatoire

!===================================================================================================
!================================     INTERFACES     ===============================================
!===================================================================================================

#ifdef MDC

interface zebmpipack
  module procedure &
    zebmpipack_int, &
    zebmpipack_double, &
    zebmpipack_char, &
    zebmpipack_intarray, &
    zebmpipack_doublearray, &
    zebmpipack_chararray
end interface zebmpipack

interface zebmpiunpack
  module procedure &
    zebmpiunpack_int, &
    zebmpiunpack_double, &
    zebmpiunpack_char, &
    zebmpiunpack_intarray, &
    zebmpiunpack_doublearray, &
    zebmpiunpack_chararray
end interface zebmpiunpack

interface zebmpisend
  module procedure &
    zebmpisend_int, &
    zebmpisend_double, &
    zebmpisend_char, &
    zebmpisend_intarray, &
    zebmpisend_doublearray, &
    zebmpisend_chararray
end interface zebmpisend

interface zebmpirecv
  module procedure &
    zebmpirecv_int, &
    zebmpirecv_double, &
    zebmpirecv_char, &
    zebmpirecv_intarray, &
    zebmpirecv_doublearray, &
    zebmpirecv_chararray
end interface zebmpirecv

#endif

!===================================================================================================
!================================     FUNCTIONS      ===============================================
!===================================================================================================

contains

!###########################################################################
!> \brief   The modified combsort algo to order a list with non regular index.
!###########################################################################
!> \ingroup Vect_calc
subroutine combsort(id,val)

integer, intent(in out) :: id(:)
integer, intent(in out) :: val(:)
integer                 :: temp_id
integer                 :: temp_val
integer                 :: i
integer                 :: gap
logical                 :: swapped = .true.

gap = size(val)
do while (gap > 1 .or. swapped)
  gap = int(gap / 1.3)
  if (gap < 1) gap = 1
  swapped = .false.
  do i = 1, size(val)-gap
    if (val(i) > val(i+gap)) then
      temp_val = val(i)
      temp_id  = id(i)
      val(i) = val(i+gap)
      id(i)  = id(i+gap)
      val(i+gap) = temp_val
      id(i+gap) = temp_id
      swapped = .true.
    end if
  end do
end do

end subroutine combsort

!###########################################################################
!> \brief   Algorithme (vectoriel) de recherche de la distance entre l'origine
!!          d'un segment et le point d'intersection de ce segment avec une
!!          surface qu'il transperce.
!###########################################################################
!> \ingroup Vect_calc
!! \param [in] IDEB  Please describe parameter
!! \param [in] IFIN  Please describe parameter
!! \param [in] INCRE Please describe parameter
!! \param [in] NMAX  Please describe parameter
!! \param [in] IA    Please describe parameter
!! \param [in] VEC   Please describe parameter
function LINEA(IDEB,IFIN,INCRE,NMAX,IA,VEC)

implicit none

integer(DPI), parameter    :: IZERO = 0
integer(DPI), dimension(3) :: IA,VEC,IAN,ITEST
integer(DPI)               :: LINEA,IDEB,IFIN,INCRE,I,NMAX(3)

do I = IDEB,IFIN,INCRE
    IAN(1:3) = IA(1:3)+I*VEC(1:3)

!ITEST(1:3) = (IAN(1:3)-1)*(IAN(1:3)-NMAX) !*** Ancienne def du reseau
!if(ITEST(1).LE.IZERO.AND.ITEST(2).LE.IZERO.AND.ITEST(3).LE.IZERO) exit

!*** Il
    ITEST(1:3) = (IAN(1:3))*(IAN(1:3)-NMAX(1:3)) !*** Nouvelle def du reseau
    if(ITEST(1).LT.IZERO.AND.ITEST(2).LT.IZERO.AND.ITEST(3).LT.IZERO) exit

!*** On part d'un point a l'exterieur et on cherche un point a l'interieur

enddo

LINEA = I

end function linea

!###########################################################################
!> \brief         Pour extraire d'un tenseur, l'amplitude d'un champ dans une direction
!!                donnee. Et vis versa.
!###########################################################################
!> \ingroup Tens_calc
!> \todo          Has to be translated in english / description of parameters is needed
!! \date          end of 1999
!! \param AMPLI   Please describe parameter
!! \param Z       Please describe parameter
!! \param XTENS   Please describe parameter
!! \param SENS    SENS=.True.   /  Amplitude uniaxiale + axe -> tenseur.
!! \param SENS    SENS=.False.  /  Tenseur + axe             -> amplitude uniaxiale
!> \ingroup Tens_calc

subroutine ORIENTti(AMPLI,Z,XTENS,SENS)

implicit none

logical         :: SENS
integer (DPI)   :: I,J,K
real(kind=DP)   :: AMPLI,XTENS(3,3),XA(3,3),XB(3,3),XC(3,3)
real(kind=DP)   :: Z(3)

if(SENS) then

! XTENS(1:3,1:3) = 0.
! do I=1,3
!  XTENS(I,1) = Z(I)*Z(1)*AMPLI
!  XTENS(I,2) = Z(I)*Z(2)*AMPLI
!  XTENS(I,3) = Z(I)*Z(3)*AMPLI
! enddo

   XA(1:3,1:3)=0.
   XB(1:3,1:3)=0.
   XC(1:3,1:3)=0.
   XTENS(1:3,1:3)=0.
   XB(3,3)=AMPLI
   XA(3,1)=Z(1)
   XA(3,2)=Z(2)
   XA(3,3)=Z(3)
   DO I=1,3
       DO J=1,3
           DO K=1,3
               XC(I,J)=XC(I,J)+XA(K,I)*XB(K,J)
           enddo
       enddo
   enddo
   DO  I=1,3
       DO  J=1,3
           DO  K=1,3
               XTENS(I,J)=XTENS(I,J)+XC(I,K)*XA(K,J)
           enddo
       enddo
   enddo
else
   DO  I=1,3
       DO  J=1,3
           XA(I,J)=0.
           XB(I,J)=XTENS(I,J)
           XC(I,J)=0.
           XTENS(I,J)=0.
       enddo
   enddo
   XA(3,1)=Z(1)
   XA(3,2)=Z(2)
   XA(3,3)=Z(3)
   DO  I=1,3
       DO  J=1,3
           DO  K=1,3
               XC(I,J)=XC(I,J)+XA(I,K)*XB(K,J)
           enddo
       enddo
   enddo
   DO  I=1,3
       DO  J=1,3
           DO  K=1,3
               XTENS(I,J)=XTENS(I,J)+XC(I,K)*XA(J,K)
           enddo
       enddo
   enddo
   AMPLI=XTENS(3,3)

endif

end subroutine ORIENTti

!###########################################################################
!> \brief         Pour extraire d'un tenseur, l'amplitude d'un champ dans une direction
!!                donnee. Et vis versa.
!! \todo          Has to be translated in english / description of parameters is needed
!###########################################################################
!> \ingroup Tens_calc
!> \date          end of 1999
!! \param AMPLI   Please describe parameter
!! \param Z       Please describe parameter
!! \param XTENS   Please describe parameter
!! \param SENS    SENS=.True.   /  Amplitude uniaxiale + axe -> tenseur.
!! \param SENS    SENS=.False.  /  Tenseur + axe             -> amplitude uniaxiale

subroutine ORIENTDIEG(AMPLI,Z,XTENS,SENS)

implicit none

logical         :: SENS
integer(DPI)    :: I
real(kind=DP)   :: AMPLI,XTENS(3,3)
real(kind=DP)   :: Z(3)

if(SENS) then

   XTENS(1:3,1:3) = 0.
   do I=1,3
       XTENS(I,1) = Z(I)*Z(1)*AMPLI
       XTENS(I,2) = Z(I)*Z(2)*AMPLI
       XTENS(I,3) = Z(I)*Z(3)*AMPLI
   enddo

else

   AMPLI = 0.
   do I=1,3
       AMPLI=AMPLI+XTENS(I,1)*Z(I)*Z(1)
       AMPLI=AMPLI+XTENS(I,2)*Z(I)*Z(2)
       AMPLI=AMPLI+XTENS(I,3)*Z(I)*Z(3)
   enddo

endif

end subroutine ORIENTDIEG

!###########################################################################
!> \brief Based on a tensile axis: TensAxe, a new reference frame is defined where
!!        TensAxe is parallel to the Z axis. Any tensor TensIn can then be redefined
!!        either from the crystal frame into the loading axis frame or reversally.
!###########################################################################
!> \ingroup Tens_calc
!> \todo    Has to be translated in english / description of parameters is needed
!! \date    unknown
!! \param   TensAxe,TensIn,TensOut  Please describe parameters
!! \param   SENS                    Please describe parameter

Subroutine TurnTens(TensAxe,TensIn,TensOut,SENS)

implicit none

logical ,intent(in)                    :: SENS
integer                                :: I,J,K
real(kind=DP),dimension(3,3)           :: a,ap,TensIn,TensOut,TensTmp
real(kind=DP),dimension(3),intent(in)  :: TensAxe
real(kind=DP),dimension(3)             :: TensAxeN,Vec1,Vec2

TensOut(:,:) = zero
TensTmp(:,:) = zero

! print*,'TensIn'
! print*,TensIn(1,1),TensIn(1,2),TensIn(1,3)
! print*,TensIn(2,1),TensIn(2,2),TensIn(2,3)
! print*,TensIn(3,1),TensIn(3,2),TensIn(3,3)

! TensAxe must be a unitary vector
TensAxeN = normavec(TensAxe)

! We first define a reference frame where the Z axis is parallel to TensAxe
a(3,1) = TensAxeN(1)
a(3,2) = TensAxeN(2)
a(3,3) = TensAxeN(3)

! Depending on the value of TensAxe, 3 simple solutions are tested.
If (abs(TensAxeN(1)) > numtols) then
   ! We calculate the new reference frame with the X axis in the x0y crystal plane
   Vec1(3) = zero
   Vec1(2) = dsqrt((un/(un+(TensAxeN(2)/TensAxeN(1))**2)))
   Vec1(1) = - Vec1(2)*(TensAxeN(2)/TensAxeN(1))
   ! The fird vector is simply defined with a cross product
   Vec2 = prodvec(TensAxeN,Vec1)

   ! Tensor a is now completly defined
   a(1,1) = Vec1(1)
   a(1,2) = Vec1(2)
   a(1,3) = Vec1(3)
   a(2,1) = Vec2(1)
   a(2,2) = Vec2(2)
   a(2,3) = Vec2(3)
!     print*,'passe1',TensAxeN
elseif (abs(TensAxeN(2)) > numtols) then
   ! We calculate the new reference frame with the X axis in the y0z crystal plane
   Vec1(1) = zero
   Vec1(3) = dsqrt((un/(un+(TensAxeN(3)/TensAxeN(2))**2)))
   Vec1(2) = - Vec1(3)*(TensAxeN(3)/TensAxeN(2))
   ! The fird vector is simply defined with a cross product
   Vec2 = prodvec(TensAxeN,Vec1)

   ! Tensor a is now completly defined
   a(1,1) = Vec1(1)
   a(1,2) = Vec1(2)
   a(1,3) = Vec1(3)
   a(2,1) = Vec2(1)
   a(2,2) = Vec2(2)
   a(2,3) = Vec2(3)
!     print*,'passe2',TensAxeN
else
   ! Any transformation is needed since the ZPrime axis is parallel to the Z crystal axis
   a(1,1) = un
   a(2,1) = zero
   a(3,1) = zero
   a(1,2) = zero
   a(2,2) = un
   a(3,2) = zero
   a(1,3) = zero
   a(2,3) = zero
   a(3,3) = un
!     print*,'passe3',TensAxeN
endif

! print*,'a'
! print*,a(1,1),a(1,2),a(1,3)
! print*,a(2,1),a(2,2),a(2,3)
! print*,a(3,1),a(3,2),a(3,3)

! The conjugate tensor of a is define
ap(1,1) = a(1,1)
ap(2,1) = a(1,2)
ap(3,1) = a(1,3)
ap(1,2) = a(2,1)
ap(2,2) = a(2,2)
ap(3,2) = a(2,3)
ap(1,3) = a(3,1)
ap(2,3) = a(3,2)
ap(3,3) = a(3,3)

! According the the SENS value we do one frame transformation or the reverse
if (SENS) then
   do i = 1,3
      do j = 1,3
         do k = 1,3
            TensTmp(i,k) = TensTmp(i,k) + a(i,j)*TensIn(j,k)
         enddo
      enddo
   enddo
! print*,'tensTmp'
! print*,TensTmp(1,1),TensTmp(1,2),TensTmp(1,3)
! print*,TensTmp(2,1),TensTmp(2,2),TensTmp(2,3)
! print*,TensTmp(3,1),TensTmp(3,2),TensTmp(3,3)
   do i = 1,3
      do j = 1,3
         do k = 1,3
            TensOut(i,k) = TensOut(i,k) + TensTmp(i,j)*ap(j,k)
         enddo
      enddo
   enddo
!    print*,'true'
else
   do i = 1,3
      do j = 1,3
         do k = 1,3
            TensTmp(i,k) = TensTmp(i,k) + ap(i,j)*TensIn(j,k)
         enddo
      enddo
   enddo
   do i = 1,3
      do j = 1,3
         do k = 1,3
            TensOut(i,k) = TensOut(i,k) + TensTmp(i,j)*a(j,k)
         enddo
      enddo
   enddo
!    print*,'false'
endif

! print*,'tensout'
! print*,TensOut(1,1),TensOut(1,2),TensOut(1,3)
! print*,TensOut(2,1),TensOut(2,2),TensOut(2,3)
! print*,TensOut(3,1),TensOut(3,2),TensOut(3,3)

End Subroutine TurnTens

!###########################################################################
!> \brief Product between a reference tensor and a scalar defining a stress
!! amplitude (SENS = .true.)
!###########################################################################
!> \ingroup Tens_calc
!> \todo    Has to be translated in english / description of parameters is needed
!! \date    unknown
!! \param   AMPLI     Please describe parameter
!! \param   tensapp   Please describe parameter
!! \param   XTENS     Please describe parameter
!! \param   SENS      Please describe parameter

subroutine TensIncr(AMPLI,tensapp,XTENS,SENS)

implicit none

logical         :: SENS
real(kind=DP)   :: AMPLI,XTENS(3,3),tensapp(3,3)

if(SENS) then
   XTENS(:,:) = ampli * tensapp(:,:)
else
   stop " TensIncr with SENS = false, is not defined"
endif

end subroutine TensIncr

!################################################################
!> \brief Application d un matrice de rotation sur un tenseur ###
!!        TENS= ROT-1 * XTENS * ROT                           ###
!################################################################
!> \ingroup Tens_calc
!! \param   Xtens       The tensor to rotate
!! \param   Xrot        The rotation tensor
!! \param   Sigprim2    The rotated tensor
subroutine TensRotation (Xtens, Xrot,Sigprim2)

implicit none

real(kind=DP)       ::Xtens(3,3),Xrot(3,3),Sigprim(3,3),Sigprim2(3,3)


SigPrim(:,:)  = 0.d0
SigPrim2(:,:) = 0.d0

SigPrim(:,:)  = matmul(Xtens(:,:),Xrot(:,:))
Sigprim2(:,:) = matmul(transpose(Xrot(:,:)),SigPrim(:,:))

end subroutine TensRotation

!###################################################################
!> \brief Calculate polar coordinates from cartesian coordinates ###
!###################################################################
!> \ingroup Vect_calc
!! \param   X       Cartesian coordinate X
!! \param   Y       Cartesian coordinate Y
!! \param   R       Cartesian coordinate R
!! \param   THETA   Cartesian coordinate THETA

SUBROUTINE RPOLARS( X, Y, R, THETA)
! Subroutine transforming input  (X, Y) to output (R, THETA)
IMPLICIT NONE
REAL(kind=DP),INTENT(IN)  :: X, Y ! cartesiancoordinates (input)
REAL(kind=DP),INTENT(OUT) :: R, THETA ! polar coordinates (output)
!REAL*8,INTENT(OUT) :: R, THETA ! polar coordinates (output)
R = SQRT( X ** 2 + Y ** 2 ) ! radius
THETA = ATAN2( Y, X ) ! inverse tangent between -pi and pi
END SUBROUTINE RPOLARS


!###########################################################################
!> \brief
!!  Calcul du rayon de courbure pour un cercle passant par trois points
!!  M,P et Q
!!  Pour un rayon de courbure infinie le vecteur retourne par la
!!  fonction est nul.
!!  Alpha : angle entre M et Q sur le cercle
!!  T : vecteur tangant ie normal a R dans MPQ
!############################################################### 30/07/99 ##
!> \ingroup Vect_calc
!> \todo    Has to be translated in english / description of parameters is needed
!! \date    30/07/99
!! \param   OM        Distance between the circle center and one of the points of left vnn for the calculation of curvature radius
!! \param   OP        Distance between the circle center and the center of the considered segment for the calculation of curvature radius
!! \param   OQ        Distance between the circle center and one of the points of right vnn for the calculation of curvature radius
!! \param   alpha     Angle between OM and OQ

FUNCTION VECCOURBCRAMER(PMint,PQint,alpha) RESULT(VECCOURB)

implicit none

real(kind=DP), dimension(3)            :: VECCOURB
real(kind=DPI), dimension(3),intent(in) :: PMint,PQint
real(kind=DP)                          :: ALPHA1, ALPHA2,normVC2,invnormVC2
real(kind=DP)                          :: DET,DETx,DETy,DETz,arg1,arg2
real(kind=DP),intent(out)              :: ALPHA
real(kind=DP), dimension(4)            :: PM,PQ,PMvecPQ
real(kind=DP), dimension(3)            :: A,B,OM,OQ,OR
real(kind=DP), parameter               :: half = 0.5D0

PM(1:3) = PMint
A(1:3) = half * PM(1:3)
PM(4) = PM(1)*A(1) + PM(2)*A(2) + PM(3)*A(3)


PQ(1:3) = PQint
B(1:3) = half * PQ(1:3)
PQ(4) = PQ(1)*B(1) + PQ(2)*B(2) + PQ(3)*B(3)


PMvecPQ(1) = (PM(2)*PQ(3) - PM(3)*PQ(2))
PMvecPQ(2) = (PM(3)*PQ(1) - PM(1)*PQ(3))
PMvecPQ(3) = (PM(1)*PQ(2) - PM(2)*PQ(1))
PMvecPQ(4)= 0.d0

!Cramer's rule to determinate the center of the circumference defined by the three
!points P,M,Q (P1,P2,P3 in the subroutine FORCE)

DET = PM(1)*PQ(2)*PMvecPQ(3) - PM(1)*PQ(3)*PMvecPQ(2) &
    - PM(2)*PQ(1)*PMvecPQ(3) + PM(2)*PQ(3)*PMvecPQ(1) &
    + PM(3)*PQ(1)*PMvecPQ(2) - PM(3)*PQ(2)*PMvecPQ(1)



if (DET < 1.d-15) then
  !The points P,M,Q are aligned, no curvature radius
  if (DET<0.d0) print *, DET
  VECCOURB(1) = 0.d0
  VECCOURB(2) = 0.d0
  VECCOURB(3) = 0.d0
  alpha=0.d0
else

  DETx = PM(4)*PQ(2)*PMvecPQ(3) - PM(4)*PQ(3)*PMvecPQ(2) &
       - PM(2)*PQ(4)*PMvecPQ(3) + PM(2)*PQ(3)*PMvecPQ(4) &
       + PM(3)*PQ(4)*PMvecPQ(2) - PM(3)*PQ(2)*PMvecPQ(4)

  DETy = PM(1)*PQ(4)*PMvecPQ(3) - PM(1)*PQ(3)*PMvecPQ(4) &
       - PM(4)*PQ(1)*PMvecPQ(3) + PM(4)*PQ(3)*PMvecPQ(1) &
       + PM(3)*PQ(1)*PMvecPQ(4) - PM(3)*PQ(4)*PMvecPQ(1)

  DETz = PM(1)*PQ(2)*PMvecPQ(4) - PM(1)*PQ(4)*PMvecPQ(2) &
       - PM(2)*PQ(1)*PMvecPQ(4) + PM(2)*PQ(4)*PMvecPQ(1) &
       + PM(4)*PQ(1)*PMvecPQ(2) - PM(4)*PQ(2)*PMvecPQ(1)

  OR(1) = DETx/DET     !x coordinate of the center of the circumference
  OR(2) = DETy/DET     !y coordinate of the center of the circumference
  OR(3) = DETz/DET     !z coordinate of the center of the circumference
  VECCOURB(1:3)=OR(1:3)  !curvature radius corresponds to the center of the circumference,
                         ! indeed the point P is located at (0,0,0) and the curvature radius is OR-P

  normVC2 = VECCOURB(1)*VECCOURB(1)+VECCOURB(2)*VECCOURB(2)+VECCOURB(3)*VECCOURB(3)
  invnormVC2 = 1/normVC2

  OM(1:3)=OR(1:3)-PM(1:3)
  OQ(1:3)=OR(1:3)-PQ(1:3)
!   if (kkdebug) then
!     write(379,*) 'in CVECCOURBCRAMER'
!     write(379,*) OM(1)*OM(1)+OM(2)*OM(2)+OM(3)*OM(3)
!     write(379,*) OQ(1)*OQ(1)+OQ(2)*OQ(2)+OQ(3)*OQ(3)
!     write(379,*) OM
!     write(379,*) OQ
!     write(379,*) VECCOURB
!     write(379,*) normVC2
!     write(379,*) (OM(1)*VECCOURB(1)+OM(2)*VECCOURB(2)+OM(3)*VECCOURB(3))
!     write(379,*) (OQ(1)*VECCOURB(1)+OQ(2)*VECCOURB(2)+OQ(3)*VECCOURB(3))
!     write(379,*) (OM(1)*VECCOURB(1)+OM(2)*VECCOURB(2)+OM(3)*VECCOURB(3))/normVC2
!     write(379,*) (OQ(1)*VECCOURB(1)+OQ(2)*VECCOURB(2)+OQ(3)*VECCOURB(3))/normVC2
!     call FLUSH(379)
!   endif
  arg1=(OM(1)*VECCOURB(1)+OM(2)*VECCOURB(2)+OM(3)*VECCOURB(3))
  arg2=(OQ(1)*VECCOURB(1)+OQ(2)*VECCOURB(2)+OQ(3)*VECCOURB(3))

  if (abs(arg1+normVC2) < 1.d-6) then
     ALPHA1=Pii*dsqrt(normVC2)
  else
     ALPHA1=Acos(arg1*invnormVC2)
  endif

  if (abs(arg2+normVC2) < 1.d-6) then
     ALPHA2=Pii*dsqrt(normVC2)
  else
     ALPHA2=Acos(arg2*invnormVC2)
  endif

  ALPHA=ALPHA1+ALPHA2
endif

END FUNCTION VECCOURBCRAMER

!##################################################################################
!> \brief  The subroutine needed to define the seeds a of random number generator Tau88()
!##################################################################################
!> \ingroup Scal_calc
SUBROUTINE init_seeds

IMPLICIT NONE

! SEEDS PAR DEFAUT DU G.U. (OPTIMISEES)
! These are unsigned integers in the C version
seed1 = 123
seed2 = 456
seed3 = 789

IF (IAND(seed1,-2) == 0)  seed1 = 123 - 1023
IF (IAND(seed2,-8) == 0)  seed2 = 456 - 1023
IF (IAND(seed3,-16) == 0) seed3 = 789 - 1023

END SUBROUTINE init_seeds

!######################################
!> \brief A radom number generator
!######################################
!> \ingroup Scal_calc
!> \todo describe parameters
!> \details Generates a random number between 0 and 1.
!> \details Translated from C function in:
!! \cite l1996maximally
!> The cycle length is claimed to be about \f[2^{88}\f] or about \f[3E^{+26}.\f]
!! Actually - \f[(2^{31} - 1).(2^{29} - 1).(2^{28} - 1).\f]
function taus88() RESULT(random_numb)

IMPLICIT NONE

REAL (kind=4) :: random_numb

INTEGER(kind=4)       :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(seed1,13), seed1), -19)
seed1 = IEOR( ISHFT( IAND(seed1,-2), 12), b)
b  = ISHFT( IEOR( ISHFT(seed2,2), seed2), -25)
seed2 = IEOR( ISHFT( IAND(seed2,-8), 4), b)
b  = ISHFT( IEOR( ISHFT(seed3,3), seed3), -11)
seed3 = IEOR( ISHFT( IAND(seed3,-16), 17), b)
random_numb = IEOR( IEOR(seed1,seed2), seed3) * 2.3283064365E-10 + 0.5

end function taus88

!######################################################################
!> \brief  Pour tester l egaliter de deux reels avec une precision donnee
!######################################################################
!> \ingroup Scal_calc
function egalite (x,y)

implicit none

real (kind=DP) :: x,y,preci
logical        :: egalite

preci = (dabs(x) + dabs(y)) / 2 * numtols
if (preci < numtols) preci = numtols
!print *, x,y,preci
!print *, x*x,y*y,dabs(x*x - y*y)

if (dabs(x - y) <= preci) then
   egalite = .TRUE.
else
   egalite = .FALSE.
endif

end function egalite

!#################################################
!> \brief Pour tester si un nombre est entier
!#################################################
!> \ingroup Scal_calc
function entier(x)

implicit none

real (dp)     :: x,z,preci,xx
integer (DPI) :: n
logical       :: entier

z = dabs (x)
preci = z * numtols
n = int(z+PRECI,DPI)
xx = real(n)

if (egalite(z,xx)) then
   entier = .TRUE.
else
   entier = .FALSE.
endif

end function entier

!###########################################################################
!> \brief Calcul de la norme d'un vecteur reel
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function norvect(r)

implicit none

real(kind=DP) :: norvect,r(3)

norvect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

end function norvect

!###########################################################################
!> \brief Calcul de la norme d'un vecteur entier
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function norivect(i)

implicit none

integer(DPI)   :: i(3)
real(kind=DP)    :: r(3),norivect

r(1:3)=real((i(1:3)),DP)
norivect=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

end function norivect

!###########################################################################
!> \brief Calcul du carre de la norme d'un vecteur entier
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function inorivect(i)

implicit none

integer (DPI):: i(3),inorivect

inorivect=i(1)*i(1)+i(2)*i(2)+i(3)*i(3)

end function inorivect

!###########################################################################
!> \brief Calcul du vecteur normalise associe a un vecteur reel
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function normavec(r)

implicit none

real(kind=DP) :: normavec(3),r(3),deno

! deno=un/dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
! normavec(:)=r(:)*deno
deno=UN/(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))            !*** division plus lourde que multiplication
normavec(:)=dsign(un,r(:))*dsqrt(r(:)**2*deno)

end function normavec

!###########################################################################
!> \brief Calcul du vecteur normalise associe a un vecteur entier
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function normaivec(i)

implicit none

integer(DPI) :: i(3)
real(kind=DP) :: normaivec(3),r(3),deno

deno = UN/norivect(i)
r=real(i,DP)
normaivec(:)=r(:)*deno

end function normaivec

!###########################################################################
!> \brief Calcul du produit vectoriel a^b
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function prodvec(a,b)

implicit none

real(kind=DP) :: prodvec(3),a(3),b(3)

prodvec(1) = a(2)*b(3)-a(3)*b(2)
prodvec(2) = a(3)*b(1)-a(1)*b(3)
prodvec(3) = a(1)*b(2)-a(2)*b(1)

end function prodvec

!###########################################################################
!> \brief Cross product of two vectors [integer]
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function prodivec(a,b)

implicit none

integer(DPI) :: prodivec(3),a(3),b(3)

prodivec(1) = a(2)*b(3)-a(3)*b(2)
prodivec(2) = a(3)*b(1)-a(1)*b(3)
prodivec(3) = a(1)*b(2)-a(2)*b(1)

end function prodivec

!###########################################################################
!> \brief  Produit scalaire entre deux vecteurs reels a.b                          #
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function prodsca(a,b)

implicit none

real(kind=DP) :: prodsca,a(3),b(3)

prodsca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

end function prodsca

!###########################################################################
!> \brief  Produit scalaire entre deux vecteurs entiers a.b                        #
!############################################################### 06/11/98 ##s
!> \ingroup Vect_calc
function iprodsca(a,b)

implicit none

integer(DPI) :: iprodsca,a(3),b(3)

iprodsca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

end function iprodsca

!###########################################################################
!> \brief  Produit scalaire entre deux vecteurs entiers a.b -  le resultats est
!>  retourne dans un reel.
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function prodisca(a,b)

implicit none

integer(DPI) :: a(3),b(3)
real (DP) :: prodisca

prodisca=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

end function prodisca

!###########################################################################
!> \brief  Cosinus entre deux vecteurs
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function cose2v(a,b)   ! cos2v(a,b) serait sans doute mieux

implicit none

real(DP) :: cose2v,a(3),b(3),na,nb

na=norvect(a)
nb=norvect(b)
if (na.ne.(0.).and.nb.ne.(0.)) then
   cose2v = dot_product(a,b)/(na*nb)
   if (abs(cose2v) > UN) cose2v=UN*sign(UN,cose2v) !avoid numerical errors
else
   cose2v = 0.
endif

end function cose2v

!###########################################################################
!> \brief  Cosinus entre deux vecteurs entiers
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function icose2v(a1,b1)   ! cos2v(a,b) serait sans doute mieux

implicit none

real(DP) :: icose2v,xxx
integer(DPI) :: a(3),b(3),a1(3),b1(3)
a = reduire_vec(a1)
b = reduire_vec(b1)

xxx = dble(dot_product(a,b))/norivect(a)
xxx = xxx/norivect(b)
if (xxx < -1.0*UN) then
   icose2v = -1.0*UN
elseif (xxx > UN) then
   icose2v = UN
else
   icose2v = xxx
endif

end function icose2v

!###########################################################################
!> \brief Signe du cosinus entre deux vecteurs si non nul
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function sicose2v(a,b)   ! cos2v(a,b) serait sans doute mieux

implicit none

real(DP) :: ricose2v
integer (DPI):: a(3),b(3),sicose2v

ricose2v = dot_product(a,b)/(norivect(a)*norivect(b))
if (ricose2v.ne.0) then
   sicose2v = NINT(ricose2v/dabs(ricose2v),DPI)
else
   sicose2v=0
endif

end function sicose2v

!###########################################################################
!> \brief  Sinus entre deux vecteurs
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function sine2v(a,b)   ! idem sin2v(a,b)   ???

implicit none

real(kind=DP) :: sine2v,a(3),b(3)

sine2v = norvect(prodvec(a,b))/(norvect(a)*norvect(b))

end function sine2v

!###########################################################################
!> \brief  Sinus entre deux vecteurs entiers
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function isine2v(a,b)   ! idem sin2v(a,b)   ???

implicit none

real(kind=DP) :: isine2v
integer (DPI)       ::      a(3),b(3)

isine2v = norivect(prodivec(a,b))/(norivect(a)*norivect(b))

end function isine2v

!###########################################################################
!> \brief  Signe du sinus entre deux vecteurs
!############################################################### 06/11/98 ##
!> \ingroup Vect_calc
function sisine2v(a,b)   ! idem sin2v(a,b)   ???

implicit none

real(kind=DP) :: risine2v
integer(DPI) :: sisine2v,a(3),b(3)

!write(*,*) 'sin1'
risine2v = norivect(prodivec(a,b))/(norivect(a)*norivect(b))
!write(*,*) 'sin2'
if (risine2v.ne.0) then
   sisine2v = NINT(risine2v/dabs(risine2v),DPI)
else
   sisine2v=0
endif

end function sisine2v

!###########################################################################
!> \brief  Calcul du signe d'une variable entiere
!############################################################### 06/11/98 ##
!> \ingroup Scal_calc
function isigne(n)

implicit none

integer(kind=DPI) :: n,m,isigne

m = n / abs (n)
if (n == 0) then
   m = 0
else
   m = n / abs (n)
endif

if ((m > 1) .OR. (m < -1)) stop 'A problem of sign definition'

isigne = m

end function isigne

!###########################################################################
!> \brief  Calcul du signe d'une variable entiere
!############################################################### 06/11/98 ##
!> \ingroup Scal_calc
function signe(x)

implicit none

real (kind=DP)       :: x,xx
integer(kind=DPI)    :: n,signe,m

m = 0
xx= real (m)
if (egalite (x,xx)) then
   n = 0
else
   n = INT(x / dabs (x),DPI)
endif

signe = n

end function signe

!###########################################################################
!> \brief  Calcul du plus grand diviseur commun
!############################################################### 06/11/98 ##
!> \ingroup Scal_calc
function pgdc3(a,b,c)

implicit none

integer(kind=DPI)               :: entier,n,test,limite, pgdc3
integer(kind=DPI),intent (in)   :: a,b,c
integer(kind=DPI)               :: aa,bb,cc

n = -1
aa = abs (a)
bb = abs (b)
cc = abs (c)
!print *, "donnees :", aa,bb,cc
if (aa + bb + cc == 0) stop ' Pas de PGDC pour valeurs nulles '

limite = max (aa,bb,cc)

if (limite == 1) then
   n = 1

else

   do entier = 1, limite

       test = modulo(aa,entier)+ modulo(bb,entier)+ modulo(cc,entier)
       if (test == 0) then
          n = entier
       endif
   enddo

endif

test = modulo(aa,n)+modulo (bb,n)+ modulo (cc,n)

if (test /= 0) stop ' Probleme de calcul de PGDC '

pgdc3 = n

end function pgdc3

!###########################################################################
!> \brief  Calcul du plus petit cummun multiplicateur de trois nombre entiers
!########################################################Ghiath 15/04/2002#
!> \ingroup Scal_calc
function ppcm3(a,b,c)

implicit none

integer(kind=DPI)               :: i, ppcm3,aa,bb,cc
integer(kind=DPI),intent (in)   :: a,b,c

aa = abs (a)
bb = abs (b)
cc = abs (c)

ppcm3 = 0
if (aa * bb * cc == 0) stop ' Pas de PPCM3 avec valeurs nulles '

Do i = 1, aa * bb * cc
    if (modulo(i,aa) + modulo(i,bb) + modulo(i,cc) == 0 ) then
       ppcm3 = i
       return
    endif
enddo

if (ppcm3 == 0) stop " mauvaise procedure de calcul de ppcm3"

end function ppcm3

!###########################################################################
!> \brief  Calcul du plus petit entier multiplicateur d'un nombre rationnel
!########################################################Ghiath 15/04/2003#
!> \ingroup Scal_calc
function ppem(a)

implicit none

integer(DPI) :: i, ppem
Real(DP),intent (in) :: a

ppem = 0
Do i = 1, 1000000
    if (entier(a*i)) then
       ppem = nint(a*i)
       return
    endif
enddo

if (ppem == 0) print *, a,"  n'est pas un nombre rationnel"

end function ppem

!###########################################################################
!> \brief  renvoie le plus petit vecteur parallel au vecteur INput
!############################################################### 01/02/01 ##
!> \ingroup Vect_calc
function reduire_vec(vec)

implicit none

integer(DPI)                             :: n
integer(DPI),intent (in), dimension (3)  :: vec
integer(DPI), dimension (3)              :: reduire_vec

n = pgdc3(vec(1),vec(2),vec(3))
reduire_vec (1:3)= vec(1:3) / n

end function reduire_vec

!###########################################################################
!> \brief  decrit la relation entre deux vecteurs x et y : si
!>    3 = egalite
!>   2 = parallels , meme sens (normes differentes)
!>    1 = angle > 0 et < 90 degree
!>    0 = perpendiculaires
!>   -1 = angle > -90 et < 0
!>   -2 = parallels et sense opposee (normes differentes)
!>   -3 = oppose (somme = 0)
!############################################################### 01/02/01 ##
!> \ingroup Vect_calc
function etat(x,y)

implicit none

integer(kind=DPI), intent (in), dimension (3) :: x,y
integer(kind=DPI), dimension (3)              :: somme
integer(kind=DPI)                             :: etat

etat = 10
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
      if (inorivect(prodivec(x,y)) == 0) then  ! vecteurs parallels
         if (dot_product(x,y) > 0) etat = 2   ! meme sens
         if (dot_product(x,y) < 0) etat = -2  ! sens opposes
         return
      else
         if (dot_product(x,y) == 0)  etat = 0 ! perpendiculaire
         if (dot_product(x,y) > 0)   etat = 1 ! 0 < angle < 90
         if (dot_product(x,y) < 0)   etat = -1! -90 < angle < 0
         return
      endif
   endif
endif

if (abs(etat) > 3) stop ' erreur dans programme etat(x,y)'

end function etat

!##########################################
!> \brief To check that two points are identical
!##########################################
!> \ingroup Vect_calc
function meme_point(x,y)

implicit none

integer(kind=DPI),intent (in), dimension (3)    :: x,y
integer(kind=DPI)                               :: dif
logical :: meme_point

meme_point = .true.

dif = abs(x(1)-y(1)) + abs(x(2)-y(2)) + abs(x(3)-y(3))
if (dif /= 0) meme_point = .false.

end function meme_point

!##########################################################################
!> \brief Fonction qui donne les coordonnees de l intersection entre un plan et  #
!>  une droite (direction)                                                 #
!>      Miller = Miller indices of the plane                                #
!>     Pos = Index of the plane                                            #
!>     OrigSeg = One point of the line                                     #
!>     VecSeg = The line vector                                             #
!##########################################################################
!> \ingroup Vect_calc
function InterPlanSeg (Miller,Pos,OrigSeg,VecSeg)

implicit none

real(kind=DP),dimension(3),intent(in)     :: Miller !< Plane normal vector (a,b,c)
real(kind=DP),intent(in)                  :: pos    !< Plane position (d)
real(kind=DP),dimension(3),intent(in)     :: OrigSeg
Integer(kind=DPI),dimension(3),intent(in) :: VecSeg

real(kind=DP),dimension(3)                :: InterPlanSeg
real(kind=DP),dimension(3)                :: VecSegreal
real(kind=DP)                             :: Interm
real(kind=DP)                             :: Coeff
real(kind=DP)                             :: scalaire


VecSegReal(1:3) = real(VecSeg(1:3),DP)

scalaire =   Miller(1)*VecSegReal(1) &
           + Miller(2)*VecSegReal(2) &
           + Miller(3)*VecSegReal(3)


if (abs(scalaire) < numtold) then
  ! The line is parallel to the plane
  InterplanSeg(:)=(/0.,0.,0./)                  ! Il n y a pas d intersection, donc zero
else

  Interm = + OrigSeg(1)* Miller(1)   &
           + OrigSeg(2)* Miller(2)   &
           + OrigSeg(3)* Miller(3)   &
           - Pos

  Coeff = -Interm/scalaire

  ! Coordinates of the line and plane intersection point
  InterPlanSeg(:) = OrigSeg(1:3) + Coeff*VecSegReal(:)
endif

end function InterPlanSeg

!########################################################################################################
!> \brief Subroutine to find the glide component of the Peach-Koehler force (the resolved sheer stress)
!> on all the slip systems.
!>  Remarque : If StressTens is a unit tensor, we compute a generelized Schmid factor
!########################################################################################################
!> \ingroup Tens_calc
Subroutine GlidCompPK(StressTens, nbase, ntsg, BvecN, NvecN, GlidComp)

implicit none

integer(kind=DPI),intent(in)    :: nbase,ntsg

real(kind=DP),dimension(3,3),intent(in)       :: StressTens
real(kind=DP),dimension(NTSG),intent(out)     :: GlidComp
real(kind=DP),dimension(3,NBASE),intent(in)   :: BvecN,NvecN
real(kind=DP),dimension(3)                    :: Vec

integer(kind=DPI)                         :: k, l

! Initializations
GlidComp(1:NTSG) = 0.0

!The loop on the slip systems
do l = 1,NTSG
  k = (l-1)*8+1

  Vec(:) = MatMul(BvecN(:,k),StressTens)
  GlidComp(l) = Dot_product(vec(:),NvecN(:,k))
enddo

end Subroutine GlidCompPK


#ifdef MDC

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiinit
  implicit none
  if (.not.zebmpiinitialized) then
    !call MPI_Init(zebmpierr)
    call MPI_Comm_get_parent(zebmpicomm, zebmpierr)
    zebmpiinitialized=.TRUE.
  else
    stop 'FATAL ERROR: try to initialized already initialized mpi'
  endif

end subroutine zebmpiinit

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpireallocrecvbuf(newbufsize)
  implicit none
  integer :: newbufsize, ierr
  character, pointer :: tmpbuf(:)

  newbufsize=2*newbufsize

  if (associated(zebmpirecvbuf)) then
    allocate(tmpbuf(newbufsize), stat=ierr)
    if(ierr/=0) stop 'FATAL ERROR: problem while allocating zebmpi receive buffer'
    tmpbuf(1:zebmpirecvbufsize)=zebmpirecvbuf(1:zebmpirecvbufsize)
    deallocate(zebmpirecvbuf)
    zebmpirecvbuf=>tmpbuf
  else
    allocate(zebmpirecvbuf(newbufsize))
  endif
  zebmpirecvbufsize=newbufsize
end subroutine zebmpireallocrecvbuf

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpireallocsendbuf(newbufsize)
  implicit none
  integer :: newbufsize, ierr
  character, pointer :: tmpbuf(:)

  newbufsize=2*newbufsize

  if (associated(zebmpisendbuf)) then
    allocate(tmpbuf(newbufsize), stat=ierr)
    if(ierr/=0) stop 'FATAL ERROR: problem while allocating zebmpi receive buffer'
    tmpbuf(1:zebmpisendbufsize)=zebmpisendbuf(1:zebmpisendbufsize)
    deallocate(zebmpisendbuf)
    zebmpisendbuf=>tmpbuf
  else
    allocate(zebmpisendbuf(newbufsize))
  endif
  zebmpisendbufsize=newbufsize
end subroutine zebmpireallocsendbuf

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipacksend
  implicit none
  call MPI_Send(zebmpisendbuf, zebmpisendbufpos, MPI_PACKED, 0, 0, zebmpicomm, zebmpierr)
  zebmpisendbufpos=0
end subroutine zebmpipacksend

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipack_int(buf)
  implicit none
  integer :: buf
  integer :: packsizeinc
  call MPI_Pack_size(1, MPI_INTEGER, zebmpicomm, packsizeinc, zebmpierr)
  if (zebmpisendbufpos+packsizeinc>zebmpisendbufsize) then
    call zebmpireallocsendbuf(zebmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, 1, MPI_INTEGER, zebmpisendbuf, zebmpisendbufsize, zebmpisendbufpos, zebmpicomm, zebmpierr)
end subroutine zebmpipack_int

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipack_double(buf)
  implicit none
  real(kind=DP) :: buf
  integer :: packsizeinc
  call MPI_Pack_size(1, MPI_DOUBLE_PRECISION, zebmpicomm, packsizeinc, zebmpierr)
  if (zebmpisendbufpos+packsizeinc>zebmpisendbufsize) then
    call zebmpireallocsendbuf(zebmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, 1, MPI_DOUBLE_PRECISION, zebmpisendbuf, zebmpisendbufsize, zebmpisendbufpos, zebmpicomm, zebmpierr)
end subroutine zebmpipack_double

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipack_char(buf)
  implicit none
  character :: buf
  integer :: packsizeinc
  call MPI_Pack_size(1, MPI_CHARACTER, zebmpicomm, packsizeinc, zebmpierr)
  if (zebmpisendbufpos+packsizeinc>zebmpisendbufsize) then
    call zebmpireallocsendbuf(zebmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, 1, MPI_CHARACTER, zebmpisendbuf, zebmpisendbufsize, zebmpisendbufpos, zebmpicomm, zebmpierr)
end subroutine zebmpipack_char

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipack_intarray(buf,n)
  implicit none
  integer :: buf(:)
  integer :: n
  integer :: packsizeinc
  call MPI_Pack_size(n, MPI_INTEGER, zebmpicomm, packsizeinc, zebmpierr)
  if (zebmpisendbufpos+packsizeinc>zebmpisendbufsize) then
    call zebmpireallocsendbuf(zebmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, n, MPI_INTEGER, zebmpisendbuf, zebmpisendbufsize, zebmpisendbufpos, zebmpicomm, zebmpierr)
end subroutine zebmpipack_intarray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipack_doublearray(buf,n)
  implicit none
  real(kind=DP) :: buf(:)
  integer :: n
  integer :: packsizeinc
  call MPI_Pack_size(n, MPI_DOUBLE_PRECISION, zebmpicomm, packsizeinc, zebmpierr)
  if (zebmpisendbufpos+packsizeinc>zebmpisendbufsize) then
    call zebmpireallocsendbuf(zebmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, n, MPI_DOUBLE_PRECISION, zebmpisendbuf, zebmpisendbufsize, zebmpisendbufpos, zebmpicomm, zebmpierr)
end subroutine zebmpipack_doublearray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipack_chararray(buf,n)
  implicit none
  character(len=*) :: buf
  integer :: n
  integer :: packsizeinc
  call MPI_Pack_size(n, MPI_CHARACTER, zebmpicomm, packsizeinc, zebmpierr)
  if (zebmpisendbufpos+packsizeinc>zebmpisendbufsize) then
    call zebmpireallocsendbuf(zebmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, n, MPI_CHARACTER, zebmpisendbuf, zebmpisendbufsize, zebmpisendbufpos, zebmpicomm, zebmpierr)
end subroutine zebmpipack_chararray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpipackrecv
  implicit none
  integer :: packsize
  call MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
  call MPI_Get_count(zebmpistatus, MPI_PACKED, packsize, zebmpierr)
  if (packsize>zebmpirecvbufsize) then
    call zebmpireallocrecvbuf(packsize)
  endif
  call MPI_Recv(zebmpirecvbuf, packsize, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
  zebmpirecvbufpos=0
end subroutine zebmpipackrecv

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiunpack_int(buf)
  implicit none
  integer :: buf
  call MPI_Unpack(zebmpirecvbuf, zebmpirecvbufsize, zebmpirecvbufpos, buf, 1, MPI_INTEGER, zebmpicomm, zebmpierr)
end subroutine zebmpiunpack_int

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiunpack_double(buf)
  implicit none
  real(kind=DP) :: buf
  call MPI_Unpack(zebmpirecvbuf, zebmpirecvbufsize, zebmpirecvbufpos, buf, 1, MPI_DOUBLE_PRECISION, zebmpicomm, zebmpierr)
end subroutine zebmpiunpack_double

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiunpack_char(buf)
  implicit none
  character :: buf
  call MPI_Unpack(zebmpirecvbuf, zebmpirecvbufsize, zebmpirecvbufpos, buf, 1, MPI_CHARACTER, zebmpicomm, zebmpierr)
end subroutine zebmpiunpack_char

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiunpack_intarray(buf, n)
  implicit none
  integer :: buf(:)
  integer :: n
  call MPI_Unpack(zebmpirecvbuf, zebmpirecvbufsize, zebmpirecvbufpos, buf, n, MPI_INTEGER, zebmpicomm, zebmpierr)
end subroutine zebmpiunpack_intarray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiunpack_doublearray(buf, n)
  implicit none
  real(kind=DP) :: buf(:)
  integer       :: n
  call MPI_Unpack(zebmpirecvbuf, zebmpirecvbufsize, zebmpirecvbufpos, buf, n, MPI_DOUBLE_PRECISION, zebmpicomm, zebmpierr)
end subroutine zebmpiunpack_doublearray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpiunpack_chararray(buf, n)
  implicit none
  character(len=*) :: buf
  integer   :: n
  call MPI_Unpack(zebmpirecvbuf, zebmpirecvbufsize, zebmpirecvbufpos, buf, n, MPI_CHARACTER, zebmpicomm, zebmpierr)
end subroutine zebmpiunpack_chararray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpisend_int(buf)
  implicit none
  integer :: buf
  call MPI_Send(buf, 1, MPI_INTEGER, 0, 0, zebmpicomm, zebmpierr)
end subroutine zebmpisend_int

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpisend_double(buf)
  implicit none
  real(kind=DP) :: buf
  call MPI_Send(buf, 1, MPI_DOUBLE_PRECISION, 0, 0, zebmpicomm, zebmpierr)
end subroutine zebmpisend_double

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpisend_char(buf)
  implicit none
  character :: buf
  call MPI_Send(buf, 1, MPI_CHARACTER, 0, 0, zebmpicomm, zebmpierr)
end subroutine zebmpisend_char

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpisend_intarray(buf, n)
  implicit none
  integer :: buf(:)
  integer :: n
  call MPI_Send(buf, n, MPI_INTEGER, 0, 0, zebmpicomm, zebmpierr)
end subroutine zebmpisend_intarray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpisend_doublearray(buf, n)
  implicit none
  real(kind=DP) :: buf(:)
  integer       :: n
  call MPI_Send(buf, n, MPI_DOUBLE_PRECISION, 0, 0, zebmpicomm, zebmpierr)
end subroutine zebmpisend_doublearray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpisend_chararray(buf, n)
  implicit none
  character :: buf(:)
  integer   :: n
  call MPI_Send(buf, n, MPI_CHARACTER, 0, 0, zebmpicomm, zebmpierr)
end subroutine zebmpisend_chararray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpirecv_int(buf)
  implicit none
  integer :: buf
  call MPI_Recv(buf, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
end subroutine zebmpirecv_int

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpirecv_double(buf)
  implicit none
  real(kind=DP) :: buf
  call MPI_Recv(buf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
end subroutine zebmpirecv_double

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpirecv_char(buf)
  implicit none
  character :: buf
  call MPI_Recv(buf, 1, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
end subroutine zebmpirecv_char

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpirecv_intarray(buf, n)
  implicit none
  integer :: buf(:)
  integer :: n
  call MPI_Recv(buf, n, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
end subroutine zebmpirecv_intarray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpirecv_doublearray(buf, n)
  implicit none
  real(kind=DP) :: buf(:)
  integer       :: n
  call MPI_Recv(buf, n, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
end subroutine zebmpirecv_doublearray

!########################################################################################################
!> \brief “please describe”
!########################################################################################################
!> \ingroup DCM
subroutine zebmpirecv_chararray(buf, n)
  implicit none
  character :: buf(:)
  integer   :: n
  call MPI_Recv(buf, n, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, zebmpicomm, zebmpistatus, zebmpierr)
end subroutine zebmpirecv_chararray

#endif

subroutine progress (current, total)
 !Shows how to place a non-advancing status counter...
 IMPLICIT NONE
 INTEGER(kind=DPI)  :: current,total
 write(*,FMT= "(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &        !achar(13) = CR ’\r’ (carriage ret)
     & " Percent Complete: ", (real(current)/real(total))*100.0,"%"
END subroutine progress

end module bricamat

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
