
!===================================================================================================
!========================    DEBUT    MODULE   "MICROSTRUCTURE"  ===================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains procedures related to the interactions between dislocations and others
!> microstructure elements such as interfaces, precipitates, free surfaces, etc ...
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module microstructure

use VARGLOB
use BRICAMAT
use VARBASE
use DEBUG
use CONNEC
use INITCONFTOOLS

#ifdef MDC
use mdctools
#endif

implicit none

logical           :: CPL_par      !<
logical           :: SEG_save     !<
logical           :: dipolaire    !<
logical           :: rognseg      !<
Real(kind=DPI)    :: V_particules !<
integer(kind=DPI) :: size_par_tab !<

integer(kind=DPI), DIMENSION(3) :: NBDIM  !<
integer(kind=DPI), DIMENSION(3) :: DECAL  !<
integer(kind=DPI), allocatable :: DEP_REAC(:,:)


!******************************************
                  contains
!******************************************
!##################################################################################################
!# This procedure reads the file fichseg containing data of the initial microstructure.           #
!# If the number of segments read is zero, then the procedure reads the file ../in/jonction       #
!# which contains information about 1 or 2 FR sources to be put "by hand" in the simulation box.  #
!##################################################################################################
subroutine lire_segments

implicit none

real (DP):: taux1,taux2,boite(3),decal1, decal2,solli1,solli2,rnpas1,rnpas2

integer (kind=DPI) :: i,carac1,carac2,ivis1,ivis2,l(3)
integer (kind=DPI) :: npas1,npas2,transla(3),centre(3)
real    (kind=DP)  :: tmp

if (cartograph == ideux) file_fichseg = "../in/SEGS"

open(3,file=file_fichseg,STATUS='OLD')

! lecture de la table de sollicitation des differents systemes
! le principe : tauapp(definitive) = tauapp(calcule) * solli_sys(sys)
! solli_sys(i) = 0 on impose tau app = 0  pour le systeme I
! solli_sys(i) = 1 on n'intervient pas sur tau app

!read(3,*) Solli_Sys(1:NTSG) : already done in init
read(3,*) tmp

! Nombre de segments dans la configuration initiale
read(3,*) Nsegm

seg_save = .false.
dipolaire = .false.

if(Nsegm /= izero) then
   ! Dimention des cotes de la boite de simulation
   read(3,*) modur(1),modur(2),modur(3)
   !With this mode we need to initialize (load) the tab FE_SigappS
   if(domvtx) return
   !if an external stress field is applied we check that those dimension are compatible with the simulation
    if (mode_deformation == IHUIT) then
      if (FEgrid_size(1) < IZERO .or. FEgrid_size(2) < IZERO .or. FEgrid_size(3) < IZERO) then
         FEgrid_size(:) = modur(:)
         FEgrid_or(:)=(/0,0,0/)
      endif

      if (dimFE_Sig(1) > FEgrid_size(1)/facteur_cris .or.     &
          dimFE_Sig(2) > FEgrid_size(2)/facteur_cris .or.     &
          dimFE_Sig(3) > FEgrid_size(3)/facteur_cris ) then
         write(*,*) '1: The FE_Sig file is not correctly defined! the grid is too fine '
        stop
      endif
      if (FEgrid_size(1) + FEgrid_or(1) > modur(1) .or.     &
          FEgrid_size(2) + FEgrid_or(2) > modur(2) .or.     &
          FEgrid_size(3) + FEgrid_or(3) > modur(3)) then
         write(*,*) '2: The FE_Sig file is not correctly defined! the size of the grid is too large'
          write(*,*)' the grid is defined outside the simulation box '
        stop
      endif
    endif

   i = nsegmax - IUN
   carac1 = izero
   carac2 = nsegm
   Raudis = zero


   do npas1 = 1, NSEGM

       ! in order to modify the saved segment file, the number of line can be different from the
       ! number of the segments. Thus we read first the segment data and then we update the number
       ! of the segment
       read (3,*) npas2,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
            seg(i)%voise,seg(i)%vnne, seg(I)%JONC,seg(i)%Ijonc,seg(i)%gd

       Raudis = Raudis + SEG(i)%NORME * NormLin(SEG(i)%VECLIN)

       if(carac2 < npas2) carac2 = npas2
       if(npas2 == 1 .and. npas1 /= 1) carac1 = npas1 - 1

       ! updating the segment number
       if(npas2 > izero) then

          seg(npas2 + carac1)%O(1:3)  = seg(i)%O(1:3)
          seg(npas2 + carac1)%norme   = seg(i)%norme
          seg(npas2 + carac1)%veclin  = seg(i)%veclin
          seg(npas2 + carac1)%JONC    = seg(i)%JONC
          seg(npas2 + carac1)%gd      = seg(i)%gd

          if (seg(i)%voiso > izero) then
             seg(npas2 + carac1)%voiso =  seg(i)%voiso  + carac1
          else
             seg(npas2 + carac1)%voiso =  seg(i)%voiso
          endif

          if (seg(i)%vnno > izero) then
             seg(npas2 + carac1)%vnno =   seg(i)%vnno + carac1
          else
             seg(npas2 + carac1)%vnno =   seg(i)%vnno
          endif

          if (seg(i)%voise > izero) then
             seg(npas2 + carac1)%voise =  seg(i)%voise + carac1
          else
             seg(npas2 + carac1)%voise =  seg(i)%voise
          endif

          if (seg(i)%vnne > izero) then
             seg(npas2 + carac1)%vnne =   seg(i)%vnne + carac1
          else
             seg(npas2 + carac1)%vnne =   seg(i)%vnne
          endif

          if(seg(i)%Ijonc /= nsegmax .and. seg(i)%Ijonc /= izero) then
             seg(npas2 + carac1)%Ijonc =  seg(i)%Ijonc + carac1
          else
             seg(npas2 + carac1)%Ijonc =  seg(i)%Ijonc
          endif

          ! The default initial value for %ijonc is the simulation starting value for nsegmax
          if (seg(i)%Ijonc == izero) seg(i)%Ijonc = nsegmax

          if(seg(npas2 + carac1)%voiso > izero .or. seg(npas2 + carac1)%voise > izero) then
             seg_save = .true.
          endif

          if(seg(npas2+carac1)%voiso < -ideux .or. seg(npas2+carac1)%voise < -ideux) dipolaire  = .true.

          if (SEG(npas2 + carac1)%VECLIN > nbase) then
             print *, " Pb with segment i index :", i
             stop
          endif

       elseif (npas2 == 0) then

          carac1 = npas1 - 1

       endif

   enddo

   ! the initial total density
   Raudis = raudis/avalue/avalue/(DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3)))

   ! Segments are now listed from 1
   do i = 1,nsegm
     seg(i) = seg(carac2 - nsegm + i)
   enddo

   if(cartograph == 2) then
      read(3,*) axe1_carto,axe2_carto
      read(3,*) phi1_carto,phi2_carto
      read(3,*) fichier_carto
   endif

   close(3)

else

   close(3)
   !   solli_sys(1:NTSG) = zero: already read in init

   open (3,FILE='../in/junction_def',STATUS='OLD')

   read(3,*) Nsegm  ! nombre de segemnts initiaux
   read(3,*) boite(1) ; ! taille de la boite den micron
   read(3,*) boite(2) ; ! taille de la boite den micron
   read(3,*) boite(3) ; ! taille de la boite den micron
   read(3,*) carac1 ! number of the line vector of Source 1
   read(3,*) taux1  ! Size of Source 1 in micron
   cpl_par = (taux1 < zero)
   read(3,*) rnpas1
   npas1 = int(rnpas1 * 1.0D-6 / avalue / normdep(carac1))
   read(3,*) decal1
   read(3,*) solli1
   raudis = taux1/boite(1)/boite(2)/boite(3)
   !solli_sys(syseg(carac1)) = solli1   : already done in init
   if(nsegm == 2) then
      read(3,*) carac2 !  number of the line vector of Source 1
      read(3,*) taux2 !  Size of Source 1 in micron
      read(3,*) rnpas2
      npas2 = int(rnpas2 * 1.0D-6 / avalue / normdep(carac2))
      read(3,*) decal2
      read(3,*) solli2
      raudis = raudis + taux2/boite(1)/boite(2)/boite(3)
      !solli_sys(syseg(carac2)) = solli2: already read in init
   endif
   ! the initial dislocation density in m-2
   raudis = raudis * 1.D12
   close(3)

   l(:) = nint(((boite(:) * 1.D-6 ) / deux/avalue),DPI)
   centre(:) = (l(1:3) / facteur_boite) * facteur_boite
   Modur(:) = centre(:)*IDEUX
   InvmoduR(1:3)=UN/Modur(1:3)
   l(:) = centre-noeud(centre,IUN)
   if (abs(l(1))+abs(l(2))+abs(l(3)) /= izero) then
      print*, "centre  mal calcule"
      print*, "facteur_boite = ", facteur_boite
      print*, "centre-noeud = ", l(:)
      stop
   endif

   ivis1 = 1 + ((carac1-1)/NBASERED)*NBASERED ! indice du vecteur vis associe

   SEG(iun)%VECLIN = carac1

   ! remise de la lognueur en metre
   taux1 = taux1 * 1.0D-6
   !====================================================================================================
   ! pour l'usage des CLP, on peut imposer une longueur de telle sort que la dislo boucle sur elle-meme
   if(modulo(carac1,ideux) == izero) then
      tmp = SQRT(real(66,DP))
   elseif(modulo(carac1-1,iquatre) == izero) then
      tmp = SQRT(real(ideux,DP))
   elseif(modulo(carac1-3,iquatre) == izero) then
      tmp = SQRT(real(3,DP))
   endif
   if(taux1 < zero)then
      seg(iun)%norme = (((NINT((modur(1) * tmp) /normlin(carac1),DPI)-iDIX) / ideux)* ideux)
      taux1 = seg(1)%norme * avalue * normlin(carac1)
   else
      !====================================================================================================

      !print *, "taille =", taille, seg(1)%veclin
      seg(iun)%norme = ((nint(taux1/avalue/normlin(carac1),DPI) / 2) * 2)
   endif

   ! le vecteur transla relie l'origine souhaite du segemnt au centre de la boite
   if (modulo(carac1,IDEUX) == 0 .and. (carac1 - ivis1) < 4) then
      l (:) = bveclin(:,ivis1)
   elseif (modulo(carac1,IDEUX) == 0 .and. (carac1 - ivis1) > 4) then
      l (:) = bveclin(:,ivis1+4)
   else
      l(:) = bvecdep(1:3,carac1)
   endif
   l(:) = bvecdep(1:3,carac1)
   transla (1:3) = npas1*l(1:3) + (IUN+nint(decal1,DPI))*(seg(1)%norme/IQUATRE)*bveclin(1:3,carac1)*IDEUX
   transla(:)= noeud (transla,IDEUX)
   if (etat(bvecnor(:,carac1),transla(:)) /= izero) then
      print*, " mauvais transla 0"
   endif

   l(1:3) = centre(1:3) - transla(1:3)

   !print *, "decalage 1/O =",sens1*npas1*l(1:3)
   !print *, "transla =",transla (1:3)
   if (etat(transla,noeud(transla,IDEUX)) /= 3) stop  " mauvais transla1"

   !if (dot_product(transla(:),bvecnor(:,carac1)) /= 0 ) stop " mauvaise translation1 dans generer_segment"
   seg(iun)%O(:) = l(1:3) ! origine en a
   SEG(iun)%Voiso = izero
   SEG(iun)%Vnno = izero
   SEG(iun)%Voise = izero
   SEG(iun)%Vnne = izero
   SEG(iun)%jonc = .false.
   SEG(iun)%Ijonc = izero
   seg(iun)%surface = izero
   seg(iun)%varfreeplan = izero
   if (nbcvxdom > IUN) then
     call assign_domain(iun,3654)
   else
     seg(iun)%dom = IUN
   endif

   if(nsegm == 2 ) then
      ! calculer pour le segment 2 seg(1)%O
      SEG(ideux)%VECLIN = carac2
      ! norme du seg 2 en bvd
      ! remise de la lognueur en metre
      taux2 = taux2 * 1.0D-6
      !print *, "taille =", taille, seg(1)%veclin
      seg(ideux)%norme = ((nint(taux2/avalue/norivect(bveclin(1:3,carac2)),DPI) * 2) / 2)
      !   seg(2)%norme = seg(1)%norme
      ivis2 = 1 + ((carac2-1)/NBASERED)*NBASERED ! vecteur vis associe

      ! le vecteur transla relie l'origine souhaite du segemnt au centre de la boite
      if (modulo(carac2,IDEUX) == 0 .and. (carac2 - ivis2) < 4) then
         l (:) = bveclin(:,ivis2)
      elseif (modulo(carac2,IDEUX) == 0 .and. (carac2 - ivis2) > 4) then
         l (:) = bveclin(:,ivis2+4)
      else
         l(:) = bvecdep(1:3,carac2)
      endif

      l(:) = bvecdep(1:3,carac2)
      !if (modulo(carac2,Ideux) == 0) npas2 = (npas2/ 55) * 55
      transla (1:3) = npas2*l(1:3) + (IUN+nint(decal2,DPI))*(seg(ideux)%norme/IQUATRE)*bveclin(1:3,carac2)*IDEUX
      !print *, "decalage 1/O =",sens2*npas2*l(1:3)
      !print *, "transla =",transla (1:3)
      transla(:)= noeud (transla,IDEUX)
      if (etat(transla,bvecnor(:,carac2)) /= izero) stop  " l hors du plan 2"
      l(:) = centre(:) - transla(:)
      if (etat(transla,noeud(transla,IDEUX)) /= 3)  stop  " mauvais transla2"
      !if (dot_product(transla(:),bvecnor(:,carac2)) /= 0 ) stop " mauvaise translation2 dans generer_segment"
      !print *, "2 ",sens2*npas, ivis2
      seg(ideux)%O(:) = l(:)
      ! definition des voin : points d'encrage
      SEG(ideux)%Voiso = izero
      SEG(ideux)%Vnno = izero
      SEG(ideux)%Voise = izero
      SEG(ideux)%Vnne = izero
      SEG(ideux)%jonc = .false.
      SEG(ideux)%Ijonc = izero
      seg(ideux)%surface = izero
      seg(ideux)%varfreeplan = izero
      if (nbcvxdom > IUN) then
        call assign_domain(ideux,3655)
      else
        seg(ideux)%dom = IUN
      endif

   endif

endif

end subroutine lire_segments

!##################################################################
!# On stoppe les segments a un rayon donne dans un sphere         #
!# centree dans la boite il faut donc calculer la distance        #
!# au centre de leur origine et extremite predite et informer     #
!# si le deplacement est compatible avec la position du JG        #
!##################################################################
subroutine barriere_spherique(EiP,OiP,bloquer)

implicit none

real(kind=DP)                          :: TST,TT(3),OiC(3),EiC(3),GB2_Size_real
logical,intent(out)                    :: bloquer
integer(kind=DPI),intent(in)           :: EiP(3),OiP(3)

GB2_Size_real =GB2_size/(AVALUE*1D6)    ! passage du micron en unite metrique

! Distance au centre du point OiP
OiC(:)=real(OiP(:))-real(Modur(:))/2.
EiC(:)=real(EiP(:))-real(Modur(:))/2.
TT(1) = norvect(OiC)
! Distance au centre du point EiP
TT(2) = norvect(EiC)

TST = MAXVAL(TT(1:2)) ! Le point le plus distant
! Bloquer =.true. si trop grand
bloquer = (TST > GB2_Size_real)


!essai barrieres planes
!***********************
!TTT = MAX(EiP(2),OiP(2))
!bloquer = (TTT > 3000)
!***************************

end subroutine barriere_spherique


!subroutine updated in order to prevent any ambiguity...
subroutine domain_barriers(i,i1o,i1e,OiP,EiP,bloquer,NdomO,NdomE,lonad,lonvoad,lonvead,domchange,surfR,SPR)

implicit none

real(kind=DP),dimension(3)    :: IntersecS  !< Point of intersection between segment and the surface
real(kind=DP),dimension(3)    :: Vect1O     !< Vector segment origin / IntersecO
real(kind=DP),dimension(3)    :: Vect1E     !< Vector segment end  / IntersecE
real(kind=DP),dimension(3)    :: normint    !< Surface burger indices (Real)
real(kind=DP),dimension(3)    :: OipReal    !<
real(kind=DP),dimension(3)    :: EipReal    !<
real(kind=DP),dimension(3)    :: TestOR


real (kind=DP)                :: DotProdO    !< Scalar Product between the vector [Oi-Intersect0] and the surface normal

real(kind=DP)                 :: TestO  !< Scalar test to define segment Origin position with respect to a surface. If strickly positif, then outside, else inside
real(kind=DP)                 :: TestE  !< Scalar test to define segment End position with respect to a surface. If strickly positif, then outside, else inside

integer(kind=DPI)             :: Flag         !<
integer(kind=DPI)             :: OldFreePlan  !< Old free surface touched index
integer(kind=DPI)             :: Surftype     !< Old free surface index
integer(kind=DPI)             :: Planes       !< Store the number of transparent surface and grain boundary surface checked through loops
integer(kind=DPI)             :: Freesurf     !< Store the number of freesurfaces
integer(kind=DPI)             :: VLI          !< Vector number for segment i (VECLIN)
integer(kind=DPI)             :: hh           !< hh={1,2} to test first and second domain in testlist
integer(kind=DPI)             :: DomO         !< Segment domain at Origin
integer(kind=DPI)             :: DomE         !< Segment domain at End
integer(kind=DPI)             :: ll           !<
integer(kind=DPI)             :: k            !< Increment variable for testlist (i,i1o,i1e)
integer(kind=DPI)             :: vnni1e       !< 1ie first non-zero neighbor at segment end
integer(kind=DPI)             :: vnni1o       !< 1io first non-zero neighbor at segment origin
integer(kind=DPI)             :: ii           !< Incremental variable for domain index
integer(kind=DPI)             :: jj           !< Incremental variable for surfaces index
integer(kind=DPI)             :: n            !<

integer(kind=DPI)             :: Testlist(3,7)  !< Storage list use to determine is segment is crossing a corner

integer(kind=DPI),intent(in)  :: i1o      !< Segment Origin neighbor
integer(kind=DPI),intent(in)  :: i1e      !< Segment End Neighbor
integer(kind=DPI),intent(in)  :: I        !< Segment I
integer(kind=DPI),intent(in)  :: lonad    !< Segment length after displacemen
integer(kind=DPI),intent(in)  :: lonvoad  !< Origin Neighbor Segment length after displacement
integer(kind=DPI),intent(in)  :: lonvead  !< End Neighbor Segment length after displacement

integer(kind=DPI),intent(in)  :: EiP(3) !< Projected value for Segment I End after displacement
integer(kind=DPI),intent(in)  :: OiP(3) !< Projected value for Segment I Origin after displacement

integer(kind=DPI),intent(out) :: NdomO  !< New domain for segment Origin
integer(kind=DPI),intent(out) :: NdomE  !< New domain for segment End
integer(kind=DPI),intent(out) :: SurfR  !< New domain for segment Origin
integer(kind=DPI),intent(out) :: SPR  !< New domain for segment End

integer(kind=DPI)             :: domchange  !< if we change from domain in barriere_conc

logical,intent(out)           :: bloquer  !<

logical                       :: bloquerSii,outsharedO,OutsharedE !<
logical, dimension(NBplanDom) :: insideO    !< True if the segment Origin is inside the domain with respect to tested surface
logical, dimension(NBplanDom) :: insideE    !< True if the segment End is inside the domain with respect to tested surface
logical, dimension(Nbcvxdom)  :: InsideDO   !< True if the segment Origin belong to the domain after have tested all surfaces
logical, dimension(Nbcvxdom)  :: InsideDE   !< True if the segment End belong to the domain after have tested all surfaces
logical, dimension(Nbcvxdom)  :: InsideDOfix   !< True if the segment Origin belong to the domain after have tested all surfaces
logical, dimension(Nbcvxdom)  :: InsideDEfix   !< True if the segment End belong to the domain after have tested all surfaces
logical, dimension(Nbcvxdom)  :: InsideDOfree   !< True if the segment Origin belong to the domain after have tested all surfaces
logical, dimension(Nbcvxdom)  :: InsideDEfree   !< True if the segment End belong to the domain after have tested all surfaces
logical, dimension(Nbcvxdom)  :: Corner_block   !< True if segment extremity is inside the corner

logical                       :: Condition       !<

if (kkdebug) then
  write(379,*)
  write(379,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(379,*) ">>> Entering domain barrier "
  write(379,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(379,*)
endif

!----------------
! Initializations
!----------------
domchange = IZERO
bloquer  = .false.
BloquerSii = .false.
Corner_block(:)=.false.
InsideO(:)=.false.
InsideE(:)=.false.
InsideDO(:)=.false.
InsideDE(:)=.false.
InsideDOfix(:)=.false.
InsideDEfix(:)=.false.
InsideDOfree(:)=.false.
InsideDEfree(:)=.false.


vnni1e=seg(i1e)%vnne
vnni1o=seg(i1o)%vnno

VLI=seg(i)%veclin

OldFreePlan = seg(i)%VarFreePlan ! Store free surface index
Surftype = seg(i)%surface

SurfR = Surftype
SPR = OldFreePlan

OipReal(:)=real(OiP(:),DP)
EipReal(:)=real(EiP(:),DP)

!Domains before displacement
DomO=seg(i)%dom
DomE=seg(i1e)%dom

!New domain for O and E are initialize to -1
NDomO=-IUN
NdomE=-IUN

if (kkdebug) then
  write(379,fmt='(A9,x,A9,x,A9)') "i1o","i","i1e"
  write(379,fmt='(I9,x,I9,x,I9)') i1o,i,i1e
  write(379,fmt='("Dom",4x,I2,8x,I2,8x,I2 )') seg(i1o)%dom,DomO,DomE
  write(379,*) " "
  write(379,fmt='("Projected extremities of I after disp. : ")')
  write(379,fmt='(3x, "Oip  = ",3I9)') Oip
  write(379,fmt='(3x, "Eip  = ",3I9)') Eip
  write(379,fmt='(A9,x,A9,x,A9)') "Lonvoad","lonad","lonvead"
  write(379,fmt='(I9,x,I9,x,I9)') lonvoad,lonad,lonvead
  write(379,*) " "
  write(379,fmt='("All surfaces equations ")')
  write(379,*) "-- Fixed "

  Do JJ = 1, NbPlanMax
      if (JJ == NbPlanMax - NbFreePlan + 1 ) write(379,*) "-- Free "
      write (379,fmt='(I7,x,A6,x,I7,3I15,F15.0)') jj,"kind",Plane_knd(jj),Plane_MillerI(1:3,jj),Plane_pos(jj)
  enddo

  write(379,*) " "
  write(379,fmt='("Oldfreeplan : ",I3 )') OldFreePlan
  write(379,fmt='("Connection to surface : ",I3 )') surftype
  write(379,*) " "
  write(379,fmt='("Domain testing order = ",50I4)') ListCVXdomsegi(:)
endif

Flag  = IZERO

!First of all we test interface, fixed barriers and shared boundaries (in concave domain)
!ONLY if segments I is inside one of this domain, the reactions with free surfaces are allowed


if (Nbplan > IZERO) then

  bc1: do hh = 1,Nbcvxdom

    ii = ListCVXdomsegi(hh)

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
    outsharedE=.false.

    do jj= planes+1,planes+NbPlanCvxDom(1,ii)

      ! Normale projection of OiP on the barrier

      normint(:) = Plane_MillerR(1:3,jj)

      TestO = OipReal(1)*normint(1)+OipReal(2)*normint(2)+OipReal(3)*normint(3)-Plane_pos(jj)
      TestE = EipReal(1)*normint(1)+EipReal(2)*normint(2)+EipReal(3)*normint(3)-Plane_pos(jj)

      if (kkdebug) then
        write (379,*) ""
        write (379,fmt='(A7,I7)') "VFP : ",JJ
        write (379,fmt='(3x," TestO",F30.3)') TestO
        write (379,fmt='(3x," TestE",F30.3)') TestE
      endif

      if (Plane_knd(jj) == IZERO) then !If impenetrable barrier

        if (TestO > zero .and. TestE > zero) then ! if the segment is outside the barrier

          InsideDO(ii)  = .false.
          InsideDE(ii)  = .false.

          insideDOfix(ii) =.false.
          insideDEfix(ii) =.false.

          if (kkdebug) write (379,fmt='("!/ Segment cross fixed boundary ", I3 ," kind = 0 with both ends")') jj
          !cycle bc1


        elseif (TestO > zero) then
          insideE(jj) = .true.

        elseif (TestE > zero) then
          insideO(jj) = .true.

        elseif (TestO <= zero .and. TestE <=zero) then
          insideO(jj) = .true.
          insideE(jj) = .true.

        endif

      elseif (Plane_knd(jj) == IUN ) then !If shared barrier

        if ( TestO > zero .and. TestE > zero ) then

          InsideDO(ii)  = .false.
          InsideDE(ii)  = .false.

          insideDOfix(ii) =.false.
          insideDEfix(ii) =.false.

          outsharedO=.true.
          outsharedE=.true.

          if (kkdebug) write (379,fmt='("!/ Segment cross shared boundary ", I3 ,"  kind = 1 with both ends")') jj

        elseif (TestO > zero) then
          insideE(jj) = .true.
          outsharedO=.true.

        elseif (TestE > zero) then
          insideO(jj) = .true.
          outsharedE=.true.


        elseif (TestO <= zero .and. TestE <= zero) then
          insideO(jj) = .true.
          insideE(jj) = .true.

        endif

      endif

    enddo

    if (All(insideO(planes+1:planes+NbPlanCvxDom(1,ii))) &
        ) then
        insideDOfix(ii) = .true.
        if(NbfreePlan == IZERO) then
          insideDO(ii) = .true.
          NdomO = ii
        endif
        if (kkdebug) write(379,*) "insideDOfix ii ==", ii ,"case 1"
    endif
    if (All(insideE(planes+1:planes+NbPlanCvxDom(1,ii))) &
        ) then
        insideDEfix(ii) = .true.
        if(NbfreePlan == IZERO) then
          insideDE(ii) = .true.
          NdomE = ii
        endif
        if (kkdebug) write(379,*) "insideDEfix ii ==", ii, "case 2"
    endif

    if ( insideDO(ii) .and. insideDE(ii) .and. ii == DomO .and. ii == DomE) return

    if (((.not. insideDOfix(ii) .and. .not. outsharedO .and. ii == domO) .or. &
         (.not. insideDEfix(ii) .and. .not. outsharedE .and. ii == domE))) then
       bloquer = .true.
       if (kkdebug) then
          write(379,*) "BLOQUER: outside fixed barrier 1"
          write(379,*) .not. insideDOfix(ii) , outsharedO , ii, domO
          write(379,*) .not. insideDEfix(ii) , outsharedE , ii, domE
       endif
    return
    endif
 enddo bc1

 if (.not. any(insideDOfix(:)) .or. .not. any(insideDEfix(:))) then
    bloquer = .true.
    if (kkdebug) write(379,*) "BLOQUER: outside fixed barrier 2"
    return
 endif

endif



bc2: do hh = 1,Nbcvxdom

  if (NbfreePlan > IZERO) then

    ii = ListCVXdomsegi(hh)

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

       TestO = OipReal(1)*normint(1)+OipReal(2)*normint(2)+OipReal(3)*normint(3)-Plane_pos(jj)
       TestE = EipReal(1)*normint(1)+EipReal(2)*normint(2)+EipReal(3)*normint(3)-Plane_pos(jj)

        if (kkdebug) then
          write (379,*) ""
          write (379,fmt='(A7,I5)') "VFP : ",JJ
          write (379,fmt='(3x," TestO",F30.3)') TestO
          write (379,fmt='(3x," TestE",F30.3)') TestE
        endif


        if (TestO >= -numtol_dotp .and. TestE >= -numtol_dotp) then

        if (kkdebug) write (379,fmt='("!/ Segment cross free surface  ", I3 ," with both ends")') jj

          insideO(jj)=.false.
          insideE(jj)=.false.

        elseif (TestO >= -numtol_dotp) then

          insideE(jj) = .true.

        elseif (TestE >= -numtol_dotp) then

          insideO(jj) = .true.

        elseif (TestO < -numtol_dotp .and. TestE < -numtol_dotp) then

          insideO(jj) = .true.
          insideE(jj) = .true.

        endif

      enddo

      if(kkdebug) then
        write(379,*) "freesurf",freesurf
        write(379,*) "NbPlanCvxDom(2,ii)",NbPlanCvxDom(2,ii)
        write(379,*) "Freesurf",All(insideO(Freesurf+1:Freesurf+NbPlanCvxDom(2,ii)))
        write(379,*) "Freesurf",All(insideE(Freesurf+1:Freesurf+NbPlanCvxDom(2,ii)))
      endif

      if (All(insideO(Freesurf+1:Freesurf+NbPlanCvxDom(2,ii)))) insideDOfree(ii)=.true.
      if (All(insideE(Freesurf+1:Freesurf+NbPlanCvxDom(2,ii)))) insideDEfree(ii)=.true.

    endif

    if ((insideDOfree(ii) .and. NbPlan == IZERO)     .or. &
        (insideDOfix(ii) .and.  insideDOfree(ii))) then
        NdomO = ii
        insideDO(ii) = .true.
        if (kkdebug) write(379,*) "insideDO ii ==", ii ,"case 4"
    endif
    if ((insideDEfree(ii) .and. NbPlan == IZERO)     .or. &
        (insideDEfix(ii) .and.  insideDEfree(ii))) then
        NdomE = ii
        insideDE(ii) = .true.
        if (kkdebug) write(379,*) "insideDE ii ==", ii,"case 5"
    endif

  if ( insideDO(ii) .and. insideDE(ii) .and. NdomO == DomO .and. NdomE == DomE &
      .and. surftype == IZERO ) return


enddo bc2

if (kkdebug) then
    write (379,*) ""
    write (379,fmt='(" Plane equation        ", 50I3 )') (jj,jj=1,NbplanMax)
    write (379,fmt='(" origin (dot prod < 0) ", 50L3 )')  insideO
    write (379,fmt='(" end    (dot prod < 0) ", 50L3 )')  insideE
    write (379,*) ""
    write (379,fmt='(" Domain list             ", 50I2 )') (ii, ii=1,NBcvxdom)
    write (379,fmt='(" Origin Inside fixed bar ", 50L2 )') insideDOfix(:)
    write (379,fmt='(" End    Inside fixed bar ", 50L2 )') insideDEfix(:)
    write (379,fmt='(" Origin Inside free surf ", 50L2 )') insideDOfree(:)
    write (379,fmt='(" End    Inside free surf ", 50L2 )') insideDEfree(:)
    write (379,fmt='(" Origin Inside domain    ", 50L2 )') insideDO(:)
    write (379,fmt='(" End    Inside domain    ", 50L2 )') insideDE(:)
    write (379,*) ""
    write(379,fmt='("Segment i ",I7,": Free surface connection :", I2 )')  I, surftype
    write (379,*) ""
endif


if (Nbcvxdom > IUN .and. ((seg(i1o)%surface /= IZERO .and. lonvoad == IZERO) .or. &
  & (seg(i1e)%surface /= IZERO .and. lonvead == IZERO))) then


! Segment i1o was touching with the Origin but become of zero length, so I is going to touch the surface
! The Origin of I has to be outside oldfreesurface
! What about the end ? Is it outside or inside a domain
! If end outside any domain -> thereafter tests
! If end is inside a domain, two choices
!   1)  The end is inside the i1o oldfreesurface => domchange = ITROIS  (similar as domchange = IUN but with cleaning i1o side)
!   2)  The end is ouside the i1o oldfreesurface => domchange = IQUATRE (similar as domchange = IDEUX but with cleaning i1o side)
! If Origin is inside a domain (not the former domain obviously), then

  if (seg(i1o)%surface /= IZERO) then
    SurfR = seg(i1o)%surface
    SPR = seg(i1o)%Varfreeplan
  else
    SurfR = seg(i1e)%surface
    SPR = seg(i1e)%Varfreeplan
  endif

  select case (SurfR)

    case (IUN)

      if (any(insideDE(:))) then

        if (insideDEfree(Plane_dom(SPR))) then
          !May or may not get i1o freesurface, to avoid ambiguity, domchange is impose for fs search
          domchange = ITROIS
          return
        else
          domchange = IQUATRE
          return
        endif

      elseif (any(insideDO(:))) then


        SPR = -1
        surfR = ISIX
!         stop "Case we want to fix lonvoad == 0 with Origin inside a domain"

        return

      endif

    case (IDEUX)! The case should not be observed

      if (any(insideDO(:))) then

        if (insideDOfree(Plane_dom(SPR))) then
          domchange = ITROIS
          return
        else
          domchange = IQUATRE
          return
        endif

      elseif (any(insideDE(:))) then


        SPR = -1
        surfR = ISIX

        return
!         stop "Case we want to fix lonvead == 0 with End inside a domain"

      endif

    case DEFAULT

      print *,'surface segments with type < 3 must not appear here (domain_barriers)'
      print *,"I1o", I1o, "surface", seg(I1o)%surface
      print *,"I1e", I1e, "surface", seg(I1e)%surface
      stop

  end select



endif


!
! if (seg(i1o)%surface /= IZERO .and. lonvoad == IZERO .and. any(insideDO(:))) then
!
!   SurfR = seg(i1o)%surface
!   SPR = seg(i1o)%Varfreeplan
!
!   if (NdomE == domE .and. insideDO(Plane_dom(SPR)) ) then
!     domchange = ITROIS
!   else
!     domchange = IQUATRE
!   endif
!
!   return
!
! elseif (seg(i1e)%surface /= IZERO .and. lonvead == IZERO .and. any(insideDE(:)) ) then
!
!   SurfR = seg(i1e)%surface
!   SPR = seg(i1e)%Varfreeplan
!
!   if (NdomO == domO .and. insideDO(Plane_dom(SPR))) then
!     domchange = ITROIS
!   else
!     domchange = IQUATRE
!   endif
!
!   return
!
! endif

Freesurf = Nbplan

if (kkdebug) then
   write(379,*) "domO & dome",DomO,DomE
   write(379,*) "NdomO & Ndome",NDomO,NDomE
   write(379,*) ".not. any(insideDO(:))",.not. any(insideDO(:)),".not. any(insideDE(:))", .not. any(insideDE(:))
endif

Condition = ((DomO == DomE    .and.                                                         &
             ((.not. any(insideDO(:)) .and. .not. any(insideDE(:)) ) .or.                   &
              (.not. any(insideDO(:))  .and. NDomE == DomE .and.                            &
               (insideDEfix(domO) .or. NbPlanCvxDom(1,domO) == izero)) .or.                 &
              (.not. any(insideDE(:))  .and. NDomO == DomO .and.                            &
               (insideDEfix(domE) .or. NbPlanCvxDom(1,domE) == izero)))) .or.               &
              (surftype /= IZERO .and. NdomO == -iun .and. NdomE == -iun) )

if (condition) then

    if (surftype < IDEUX) then
      ii = DomO
    else
      ii = DomE
    endif

    if(.not. any(insideDO(:))) NdomO = ii
    if(.not. any(insideDE(:))) NdomE = ii

    if (kkdebug) then
      write(379,*) " "
      write (379,fmt='("-- Checking FREE surfaces - O and E both in domain ii" )')
    endif

    n=IZERO

    if (ii > IUN) then
      do jj=1,ii-1
        freesurf  = freesurf +  NbPlanCvxDom(2,jj)
      enddo
    endif

    do jj=freesurf+1,freesurf+NbPlanCvxDom(2,ii)
      if ( ( .not. insideE(jj) .and. .not. (surftype == IDEUX .and. jj == OldfreePlan ) ) .or. &
           ( .not. insideO(jj) .and. .not. (surftype == IUN .and. jj == OldfreePlan ) ) ) then
        flag = flag + 1
        n=jj
      endif
    enddo

    if ((flag == IZERO .or. n < IUN ) .and. surftype < IUN) then
      write (379,*) "flag=izero !!Problem!! kk=",kk
      stop "flag=izero !!Problem!!"
    endif

    if (n > IZERO) then
      if (Flag > iun) then
          surfR=-IUN
          bloquer = .true.
          if (seg(i)%surface == IZERO) seg(i)%diseg=.TRUE.
          if (kkdebug) then
            write (379,*) ""
            write (379,fmt='("!/ TWO SURFACES CROSSED ")')
            write (379,fmt='("Segment I is blocked = ",L2," because Flag = ", L2)') Bloquer,Flag
          endif



          return
      else

        !Test for the surface case 6
        if  (.not. insideO(n) .and. .not. insideE(n) .and. surftype < IUN ) then

          surfR = ISIX !the segment i is crossing the surface jj with both side
          SPR = -IUN

          if (kkdebug) write (379,*) "!/ Flag surface CASE 6 - crossing the surface with both extremities "

          return

        endif

        if (.not. insideO(n) .and. .not. insideE(n) .and. (lonad == IZERO .or. surftype /=IZERO)) then !surftype /= IZERO

          if (Surftype == IUN) then
            if(seg(i1e)%dom == Plane_dom(OldFreePlan) .or. lonad /= IZERO) then
              !print *, "SOLVING CASE IDOUZE"
              if (lonvead == IZERO) then
                out(i) = .True.
              else
                surfR = IDOUZE
                SPR = -IUN
              endif
              return
            else
              !print *, "Domchange ICINQ SURF IUN"
              SurfR = IDOUZE
              SPR = -IUN
!               SurfR = IUN
!               SPR = OldFreePlan
!               domchange = ICINQ
              return
            endif

          elseif(Surftype == IDEUX) then
            if (seg(i1o)%dom == Plane_dom(OldFreePlan) .or. lonad /= IZERO) then
              !print *, "SOLVING CASE ITREIZE"
              if (lonvoad == IZERO) then
                out(i) = .True.
              else
                surfR = ITREIZE
                SPR = -IUN
              endif
              return
            else
              !print *, "Domchange ICINQ SURF IDEUX"
              SurfR = ITREIZE
              SPR = -IUN
              !SurfR = IDEUX
              !SPR = OldFreePlan
              !domchange = ICINQ
              return
            endif

          endif

        endif

        if  ( .not. insideO(n)         .and. &
              surftype /= IUN) then

          !Test for the surface case 4
          if (surftype /= IDEUX) then

              TestE= dot_product(BVECLIN(:,seg(i1o)%veclin),Plane_MillerI(1:3,n))
              TestO= dot_product(BVECLIN(:,seg(i)%veclin),Plane_MillerI(1:3,n))
              if(TestO*TestE < IZERO) then

              surfR  = IQUATRE !the segment i is crossing the surface jj with the origin
              SPR = - IUN

              if (kkdebug) write (379,*) "!/ surface CASE 4 - Segment crossing the surface with the origin "

              return
            else

              surfR = IQUATRE
              SPR = -1
              if (kkdebug) then
                write (379,*) "!/ Flag surface CASE 4 or 8 - displacement reduced because of incompatible rails"
                write (379,*) "seg(i1o)%veclin",seg(i1o)%veclin, " Bveclin", BVECLIN(:,seg(i1o)%veclin)
                write (379,*) "seg(i)%veclin",VLI, " Bveclin", BVECLIN(:,VLI)
              endif
              return

            endif

            ! Test for the surface case 8
          elseif (seg(i1o)%surface < IUN) then !

            surfR = IHUIT
            SPR = n

            if (kkdebug) write (379,*) "!/ surface CASE 8 "
            return

          else
            !write(*,*) "case we want to fix 8, kk= ", kk ,"i = ", i

            bloquer = .true.
            surfR=-1
            if (kkdebug) write (379,*) "!/ surface CASE 8 (else): small loop touching free surface that will be eliminated in NET"
            return

          endif

        endif

        if  ( .not. insideE(n)          .and. &
             surftype < IDEUX ) then


          ! Test for the surface case 5
          if (surftype /= IUN) then

            TestE= dot_product(BVECLIN(:,seg(i)%veclin),Plane_MillerI(1:3,n))
            TestO= dot_product(BVECLIN(:,seg(i1e)%veclin),Plane_MillerI(1:3,n))
            if(TestO*TestE < IZERO) then

              surfR     =  ICINQ !the segment i is crossing the surface jj with the end
              SPR =  -IUN

              if (kkdebug) write (379,*) "!/ Flag surface CASE 5 - Crossing the surface with the end"

              return

            else

              surfR     =  ICINQ
              SPR =   -IUN
              if (kkdebug) then
                write (379,*) "!/ Flag surface CASE 5 or 7 - displacement reduced because of incompatible rails"
                write (379,*) "seg(i1e)%veclin",seg(i1e)%veclin, " Bveclin", BVECLIN(:,seg(i1e)%veclin)
                write (379,*) "seg(i)%veclin",VLI, " Bveclin", BVECLIN(:,VLI)
              endif
              !seg(i)%diseg=.true.
              !seg(i1e)%diseg=.true.

              return

            endif

          else

            ! Test for the surface case 7
            if (seg(i1e)%surface < IUN) then

              SurfR   = ISEPT
              SPR  = n

              if (kkdebug) write (379,*) "!/ Flag surface CASE 7 "

              return

            else

              bloquer = .true.
              surfR=-1
              if (kkdebug) write (379,*) &
                          "!/ surface CASE 7 (else): small loop touching free surface that will be eliminated in NET"
              return

            endif

          endif

          if (kkdebug) write (379,*) "!/ Flag surface CASE 5,7 final "

        endif

      endif
    endif

else

  if((.not. any(insideDO(:)) .or. .not. any(insideDE(:)))) then
    if (NBCVXdom == 1) then
      write(*,*) "kk = ",kk
      stop "not possible with only one domain"
    endif
    if (surftype == IZERO) then
      if (.not. any(insideDO(:)) .and. .not. any(insideDE(:))) then
        !write(*,*) "case we want to fix 1 - 6, kk= ", kk ,"i = ", i

        surfR= ISIX
        SPR = -IUN
        return
      elseif(.not. any(insideDO(:))) then

        surfR= IQUATRE
        SPR = -IUN
        return
      elseif(.not. any(insideDE(:))) then

        surfR= ICINQ
        SPR = -IUN

        return
      endif
    else
      if (.not. any(insideDO(:)) .and. surftype == IDEUX) then
         surfR= ITREIZE
         SPR = -IUN
         return
      endif
      if (.not. any(insideDE(:)) .and. surftype == IUN) then
         surfR= IDOUZE
         SPR = -IUN
         return
      endif
    endif
  endif

 !
!In case the segment i changes domain we check that each point of i1o,i,i1e is still
!inside the volume

  if  (  NDomO /= NDomE  .or.  &
         DomO  /= DomE   .or.  &
         DomO  /= NdomO  .or.  &
         DomE  /= NDomE        &
      ) then

    Testlist(1,1:7)=(/i1o,seg(i1o)%dom,NdomO,lonvoad,seg(i1o)%O(1),seg(i1o)%O(2),seg(i1o)%O(3)/)
    Testlist(2,1:7)=(/i,NDomO,NdomE,lonad,OiP(1),OiP(2),OiP(3)/)
    Testlist(3,1:7)=(/i1e,NDomE,seg(seg(i1e)%voise)%dom,lonvead,EiP(1),EiP(2),EiP(3)/)

    if (kkdebug) then
      write (379,fmt='(">> Testlist" )')
      write (379,fmt='(7A9)') "i1o","Domo","NdomO","lonvoad","O(1)","O(2)","O(3)"
      write (379,fmt='(7I9)') Testlist(1,1:7)
      write (379,fmt='(7A9)') "i","NDomO","NdomE","lonad","Oip(1)","Oip(2)","Oip(3)"
      write (379,fmt='(7I9)') Testlist(2,1:7)
      write (379,fmt='(7A9)') "i1e","NDomE","DomE","lonvead","EiP(1)","EiP(2)","EiP(3)"
      write (379,fmt='(7I9)') Testlist(3,1:7)
    endif

    do k = 1, 3
      if (kkdebug) then
         write(379,*) ''
         write(379,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         write(379,*) 'Enter k loop ( 1=i1o | 2=i | 3=i1e )',k

      endif
      if ((surftype == iun .and. (k==2 .or. k==1)) .or. &
        (surftype == ideux .and. (k==2 .or. k==3))) then
        if (kkdebug) write(379,*) " Cycle because over k because either &
              & domain are the same of surface is already defined"
        cycle
      endif

      VLI = seg(Testlist(k,1))%veclin
      if (kkdebug) write(379,*) "testing VECLIN,segment: ",VLI,Testlist(k,1)

      do ll= IUN,Testlist(k,4)-IUN !Increment length

        Corner_block(:)=.false. !Initialize Corner_block to false

        TestOR(:)=real(modulo(Testlist(k,5:7)+ll*bveclin(:,VLI),modur),DP)
        if (kkdebug) write(379,*) "testing point: ",TestOR


        do hh = 1,NBcvxdom !hh to check both domain in teslist
          !if (hh == 1) then
          !  ii=Testlist(k,2) !Take first Domain
          !else
          !  ii=Testlist(k,3) !Take second Domain
          !endif
          ii = hh
          if (kkdebug) write(379,*) ">> Testing the extremity for domain ",ii, "hh =", hh

          planes=IZERO   !Plane is the increment of fixed boundary index
          freesurf=Nbplan !Planes (with a 's') is the increment of free boundary index

          if (ii > IUN) then !if it is not the first domain, increment plane and planes variables
            do jj=1,ii-1
              planes= planes+NbPlanCvxDom(1,jj)
              freesurf= freesurf+NbPlanCvxDom(2,jj)
            enddo
          endif

          !Test for fixed boundaries
          do jj= planes+1,planes+NbPlanCvxDom(1,ii)
            normint(:) = PLane_MillerR(1:3,jj)
            TestO=TestOR(1)*normint(1)+TestOR(2)*normint(2)+TestOR(3)*normint(3)-Plane_pos(jj)

            if (kkdebug) write(379,*) "TestO for fixed surface",jj," is ",TestO
            if (testO > numtol_dotp) then
              if (kkdebug) write(379,*) "As TestO is positive, the extremity ",hh," is assume inside a corner"
              Corner_block(hh)=.true.
              exit
            endif

          enddo

          !Test for free surfaces
          if (.not. Corner_block(hh)) then
            do jj =freesurf+1,freesurf+NbPlanCvxDom(2,ii)

              normint(:) = PLane_MillerR(1:3,jj)
              TestO=TestOR(1)*normint(1)+TestOR(2)*normint(2)+TestOR(3)*normint(3)-Plane_pos(jj)
              if (kkdebug) write(379,*) "TestO for free surface",jj," is ",TestO

              if (testO >= -numtol_dotp) then
              if (kkdebug) write(379,*) "As TestO is positive, the extremity ",hh," is assume inside a corner"
                Corner_block(hh)=.true.
                if (k == 2 .and. surftype < IUN) then
                  surfR = INEUF
                elseif (k == 1 .and. surftype < IUN) then
                  surfR = IDIX
                elseif (k == 3 .and. surftype < IUN) then
                  surfR = IONZE
                !else
                !  surfR=-1
                endif
                exit
              endif

            enddo
          endif

        enddo  !!End of the hh loop

        if (All(Corner_block)) then
          exit
        elseif(surfR > IHUIT)  then
          surfR = surftype
        endif
      !if (All(Corner_block)) then
      ! exit
      !else

      enddo ! End of the LL loop

      if (All(Corner_block)) then
        exit
      elseif(surfR > IHUIT)  then
        surfR = surftype
      endif

    enddo ! End of the k loop

  endif

  if ( All(Corner_block) .and. (surfR < INEUF) ) then

    bloquer=.true.
    seg(i)%diseg=.true.
    if (kkdebug) then
      write(379,*) ""
      write(379,*) '!!!!!!!!!!!!!! OBSTACLE DETECTED DURING CORNER BLOCK!!!!!!!!!!!!!!!!!!!'
      write(379,*) ""
    endif
    print *,kk
    write(379,*) "!!!!!!!!!!!!!! OBSTACLE DETECTED DURING CORNER BLOCK!!!!!!!!!!!!!!!!!!!"
    return


  endif

endif

  ! Here, we calculate the position of the segment extremity in order to reconnect the freesurface ( Same Idea as freesurface hooking up)
  ! If a segment extremity is inside an other domain, then we want to reduce the displacement to the smallest displacement required to change domain

if (surftype > IZERO) then


  if (kkdebug) then
    write(379,*) ''
    write(379,*) 'TREATMENT OF SURFACE SEGMENT DOMAIN CHANGEMENT'
    write(379,*) ''
  endif

    ! Intersection point on the free surface
    VLI=seg(i)%veclin
    normint(:)    = Plane_MillerR(1:3,OldFreePlan)
    IntersecS(:)  = InterPlanSeg(normint(:),Plane_pos(OldFreePlan),OiPReal,Bveclin(:,VLI)) !intersection with free surface after displacement!!

  if (kkdebug) then
    write(379,*) 'intersectS', IntersecS(:)
    write(379,*) 'OipReal',OiPReal
    write(379,*) 'bveclin',Bveclin(:,seg(I)%VECLIN)
  endif


  Vect1O(:) = OipReal(:)-IntersecS(:)
  Vect1E(:) = EipReal(:)-IntersecS(:)

  if(dot_product(Vect1O,normint)>=-numtol_dotp .and. &
    (dot_product(Vect1E,normint)>=-numtol_dotp)) then

   domchange = IDEUX !segment crossed the free surface with both extremity


    if (kkdebug) then
      write (379,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write (379,*) "!/ Displacement reduction due to domain change of projected extremity"
      write (379,*) ""
      write (379,fmt='(" Plane equation        ", 50I3 )') (ii,ii=1,NbplanMax)
      write (379,fmt='(" origin (dot prod < 0) ", 50L3 )')  insideO
      write (379,fmt='(" end    (dot prod < 0) ", 50L3 )')  insideE
      write (379,*) ""
      write (379,fmt='(" Domain list             ", 50I2 )') (ii, ii=1,NBcvxdom)
      write (379,fmt='(" Origin Inside domain    ", 50L2 )') insideDO(:)
      write (379,fmt='(" End    Inside domain    ", 50L2 )') insideDE(:)
      write (379,*) ""
      write(379,fmt='("Segment ",I7,": Free surface connection :", I2 )')  I, surftype
      write (379,*) ""
      write(379,fmt='("************** change domain key *************** ",I2 )') domchange
      write (379,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    endif

  else

    Do jj= 1,NbPlanMax - NbFreePlan ! Loop over fixed surfaces

      !Test only surfaces in the old domain minus the connected free surface, cycle
      If ( Plane_dom(jj) /= Plane_dom(Oldfreeplan) .or. jj == Oldfreeplan .or. &
        Plane_knd(jj) == IZERO) cycle

      normint(:) = Plane_MillerR(1:3,jj)
      DotProdO = IntersecS(1)*normint(1)+IntersecS(2)*normint(2)+IntersecS(3)*normint(3)-Plane_pos(jj)
      if (kkdebug) write (379,fmt='(x,"VFP",x,I2,x,"dom",x,I2,x," DotProdO",F15.3)') &
                                     & jj,Plane_dom(jj),DotProdO

      if (DotProdO > numtol_dotp) then

        domchange = IUN !segment crossed the free surface with one extremity

        if (kkdebug) then
          write (379,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (379,*) "!/ Displacement reduction due to domain change of projected extremity"
          write (379,*) ""
          write (379,fmt='(" Plane equation        ", 50I3 )') (ii,ii=1,NbplanMax)
          write (379,fmt='(" origin (dot prod < 0) ", 50L3 )')  insideO
          write (379,fmt='(" end    (dot prod < 0) ", 50L3 )')  insideE
          write (379,*) ""
          write (379,fmt='(" Domain list             ", 50I2 )') (ii, ii=1,NBcvxdom)
          write (379,fmt='(" Origin Inside domain    ", 50L2 )') insideDO(:)
          write (379,fmt='(" End    Inside domain    ", 50L2 )') insideDE(:)
          write (379,*) ""
          write(379,fmt='("Segment ",I7,": Free surface connection :", I2 )')  I, surftype
          write (379,*) ""
          write(379,fmt='("************** change domain key *************** ",I2 )') domchange
          write (379,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        endif

        return
      endif

    Enddo

  endif

  if (domchange == IZERO) then


     if (surftype == IUN) NdomO = Plane_Dom(Oldfreeplan)
     if (surftype == IDEUX) NdomE = Plane_dom(Oldfreeplan)

    if (kkdebug) then
      write (379,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write (379,*) "!/ any domain changement detected"
      write (379,*) ""
      write (379,fmt='("Old DomO :", I3," New DomO ",I3 )') domO,NdomO
      write (379,fmt='("Old DomE :", I3," New DomE ",I3 )') domE,NdomE

      write(379,fmt='("************** change domain key *************** ",I2 )') domchange
      write (379,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      endif
  endif

endif

end subroutine domain_barriers

!################################################################################
!# dans le cas d'une barriere spherique de rayon R, il faut reprendre le        #
!# fichier des sources de dislo initial (contenu dans un cube de dimension R    #
!# et supprimer tous les segments qui sont a l exterieur du rayon.              #
!################################################################################
subroutine Selection_spherique

implicit none

Real(kind=DP)            :: TST(2),GB2_Size_real,alpha,CCC,DELTA,AAA,BBB,seg_real
Integer(kind=DPI)        :: Ei(3),modEi(3),modOi(3),NbSeg,O(3)
integer (kind=DPI)       :: i,TT1(3),TT2(3),VectSeg(3)
character(len=1)         :: carac
logical                  :: NOTOK

!** lecture du diametre de la sphere
open(444,file='../in/b_spher_def',status="old")

!*lecture des infos concernant les barrieres planes
carac = 'x'
do while (carac /= "#")
    read (444,*) carac
enddo

read(444,*) GB2_size

close (444)

GB2_Size_real =GB2_size/(AVALUE*1D6)    ! passage du micron en unite metrique
if (2*GB2_size_REAL >= modur(1).or.2*GB2_size_REAL >= modur(2).or.2*GB2_size_REAL >= modur(3)) &
   stop 'Stop: The spheric grain is larger than the simulation cell'

if (RognSeg) then ! for a restart there is no need to cut segments  RognSeg=.false.

  open(444,file='../in/fichseg_int',STATUS="REPLACE")

  Nbseg=NSEGM
  !pour chaque segment on calcule la distance au centre de son origine et extremite TT1(1) et TT1(2)
  do I = 1, NSEGM
!    if(seg(i)%voise == -IUN) seg(i)%voise = IUN
    modOi(:)=modulo(seg(I)%O(:),modur(:))
    TT1(:)=modOi(:)-(MODUR(:)/2)
    TST(1)=real(norivect(TT1))
    Ei(:)=Seg(I)%O(:)+((seg(I)%norme)*(BVECLIN(:,seg(I)%veclin)))
! (modulo)entre guillemets car c est la valeur de l extremite par rapport
! au modulo de l origine
    modEi(:)=modOi+((seg(I)%norme)*(BVECLIN(:,seg(I)%veclin)))
    VectSeg(:)=-seg(I)%O(:)+Ei(:)
    TT2(:)=modEi(:)-(MODUR(:)/2)
    TST(2)=real(norivect(TT2))
!print*,'de normes',TST,"barriere",GB2_size_real;pause

!si TT1 et TT2 sont derriere la barriere on ne met pas le segment ds le fichier intermediaire
    if ((TST(1) > GB2_size_real).AND.(TST(2) > GB2_size_real)) then
       Nbseg=Nbseg-1
!                  print*,NbSeg
!                  print*,'dehors',I;pause
       cycle  !on passe au segment suivant

!si TT2 (l extremite) est en dehors on raccourci la norme
    elseif (TST(2) > GB2_size_real) then
! resolution d une equation du second degre pour calculer l intersection
! entre la barriere spherique et le segment
       AAA=real(inorivect(VectSeg))
       BBB=real(2*dot_product(VectSeg,TT1))
       CCC=real(inorivect(TT1))-(GB2_size_real)**2
       DELTA=BBB**2-4*AAA*CCC
       alpha=(-BBB+sqrt(DELTA))/(2*AAA)
       seg_real=real(seg(I)%norme)
       Seg(I)%norme= int(alpha*seg_real)
!probleme si la norme est a zero il vaut mieux supprimer le segment
       if (seg(i)%norme == Izero) then
          NbSeg=NbSeg-1
!                 print*,'extr dehors mais annule';pause
          cycle  !on sort
       endif
! l extremite est consideree comme point d ancrage
       seg(I)%voise=IZERO
!print*,'extr dehors',I;pause

!si TT1 est en dehors on raccourci la norme et on place l origine sur la barriere
    elseif (TST(1) > GB2_size_real) then
       AAA=real(inorivect(-VectSeg))
       BBB=real(2*dot_product(-VectSeg,TT2))
       CCC=real(inorivect(TT2))-(GB2_size_real)**2
!         print*,'origine dehors calcul du determninant',AAA,BBB,CCC
       DELTA=BBB**2-4*AAA*CCC
       alpha=(-BBB+sqrt(DELTA))/(2*AAA)
       seg_real=real(seg(I)%norme)
       Seg(I)%norme=int(alpha*seg_real)
!probleme si la norme est a zero il vaut mieux supprimer le segment
       if (seg(i)%norme == Izero) then
          NbSeg=NbSeg-1
!                 print*,'origine dehors mais annule';pause
          cycle  !on sort
       endif
! attention aux modulos pour le calcul de la nouvelle origine
       NOTOK = .true.
       do while (NOTOK)
         Seg(I)%O(:)=modEi(:)-((seg(I)%norme)*(BVECLIN(:,seg(I)%veclin)))
         ! Pinned point must be on the simulation lattice
         O(:) = seg(i)%O(:)
         O(:) = noeud(O,IDEUX)
         print*,'notok',O,seg(i)%O(:)
         if (O(1) /=  seg(i)%O(1) .or. O(2) /=  seg(i)%O(2) .or. O(3) /=  seg(i)%O(3)) then
           seg(I)%norme = seg(I)%norme - 1
         else
           NOTOK = .false.  ! O is correctly on the lattice
         endif
       enddo
       if (seg(I)%norme <= 0) then
         NbSeg=NbSeg-1
         cycle  ! We stop since kneecaps are not allowed here
       endif
       seg(I)%voiso=IZERO
    endif

!on place le segment dans le fichier intermediaire
333   format(x,7(x,I7))
    write (444,333)seg(I)%O(1:3),seg(I)%norme,SEG(I)%VECLIN,seg(I)%voiso,seg(I)%voise

  enddo

  !nouvelle valeur de NSEGM
  NSEGM=Nbseg

  !on recopie les segments qui nous interressent
  rewind(444)
  do I = 1, NSEGM

    read(444,*)seg(I)%O(1:3),seg(I)%norme,SEG(I)%VECLIN,seg(I)%voiso,seg(I)%voise
    !print*,i,seg(I)%O(1:3),seg(I)%norme,SEG(I)%VECLIN,seg(I)%voiso,seg(I)%voise
    if (SEG(I)%voiso /= izero) then
       SEG(I)%Voiso=I-1
    else
       SEG(I)%Voiso=izero
    endif
    if (seg(I)%voise /= izero) then
       SEG(I)%Voise=I+1
    else
       SEG(I)%Voise=izero
    endif

  enddo

  !le fait d avoir elimine des segments de norme nulle amene a faire des
  !corrections au niveau des indices du voisinages
  do I = 2,NSEGM-1
    if(seg(i)%voiso /= izero .and. seg(I-1)%voise == izero) seg(i)%voiso =izero
    if(seg(i)%voise /= izero .and. seg(I+1)%voiso == izero) seg(i)%voise =izero
  enddo

  close(444,status="delete")

endif  !End of RognSeg

end subroutine Selection_spherique

!############################################################################
!# Dans le cas d un systeme de barrieres planes il faut d abord recuperer   #
!# les donnees relatives aux barrieres dans un fichier d entree puis        #
!# rogner tous les segments qui coupent ces barrieres                       #
!# corrections are added to allow dislocations loops in the initial         #
!# configuration                                                            #
!############################################################################

!**********************************************
!*                                            *
!**********************************************
subroutine domain_definition

implicit none

real(kind=DPI),DIMENSION(3)  :: modOi

integer(kind=DPI)               :: I
integer(kind=DPI)               :: J
integer(kind=DPI)               :: JJ
integer(kind=DPI)               :: iiold
integer(kind=DPI)               :: length
integer(kind=DPI)               :: itemp
integer(kind=DPI)               :: maxlength
integer(kind=DPI)               :: lsegm
integer(kind=DPI)               :: iref
!integer(kind=DPI)               :: Li
integer(kind=DPI),DIMENSION(3)  :: O

!integer(kind=DPI),DIMENSION(3)  :: int_gpFS
!integer(kind=DPI),DIMENSION(3)  :: normal_ext
!real(kind=DP)                   :: norm

integer(kind=DPI),DIMENSION(100)  :: loop

logical,allocatable           :: bloquer(:)

logical                       :: FRsource
logical                       :: FlagListType
logical,dimension(NSEGM)      :: skip
logical                       :: keep
logical                       :: outside
character(len=256)            :: cmd
logical                       :: file_exist


!Initialization
FlagListType = .false.

! Comments in the file header are removed
! carac = 'x'
!do while (carac /= "#")
!  read (44,*) carac
!enddo

! Boundaries definition is loaded
! call load_boundaries

!!!! subroutine to compute geometry and volume of geometry defined in domain_def
!!!! the geometry is stored in the out folder in the geometry.vtk file
#ifdef PA
  if (Mon_Rang==IZERO) then
    ! Only the process zero of the communicator can write in output files, the others are waiting

#endif

! The boundary vertex calculation
call vertex_calculations
INQUIRE(FILE="../bin/volume.py", EXIST=file_exist)
if (file_exist) then
  cmd = "../bin/volume.py"
else
  cmd = "volume.py"
endif

call system(cmd)

#ifdef PA
  endif
  ! Proc zero ended is task, we can now liberate everybody
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

#endif

! #ifdef MDC
! allocate(fsurfext(NTSG,nbFreePlan,3))
! allocate(fsurfext_norm(NTSG,nbFreePlan))
! !print *,"NTSG", NTSG , "NbFreeplan",NbFreeplan
!
!# Calculation of the vectors needed to extend area swept for segment touching the surface.
!  do i=1,NTSG
!    Li=1+8*(i-1)
!    print *, 'Li',Li
!    do j=nbplan+1,nbplan+nbFreePlan
!
!      int_gpFS (1) = bvecnor(2,Li)*VarPlan(j)%Miller(3)-bvecnor(3,Li)*VarPlan(j)%Miller(2)
!      int_gpFS (2) = bvecnor(3,Li)*VarPlan(j)%Miller(1)-bvecnor(1,Li)*VarPlan(j)%Miller(3)
!      int_gpFS (3) = bvecnor(1,Li)*VarPlan(j)%Miller(2)-bvecnor(2,Li)*VarPlan(j)%Miller(1)
!
!      normal_ext(1) = bvecnor(2,Li)*int_gpFS(3)-bvecnor(3,Li)*int_gpFS(2)
!      normal_ext(2) = bvecnor(3,Li)*int_gpFS(1)-bvecnor(1,Li)*int_gpFS(3)
!      normal_ext(3) = bvecnor(1,Li)*int_gpFS(2)-bvecnor(2,Li)*int_gpFS(1)
!
!      norm=dsqrt(real(normal_ext(1)**2+normal_ext(2)**2+normal_ext(3)**2,dp))
!      if (norm > numtolD) norm=1.d0/norm
!
!      if (dot_product(VarPlan(j)%Miller(1:3),normal_ext(1:3)) < IZERO)  normal_ext(1:3)=-normal_ext(1:3)
!
!      fsurfext(i,j,1:3) = normal_ext(1:3)
!      fsurfext_norm(i,j) = norm
!
!      !write(379,*) "glide plane ", i, " FS " ,  VarPlan(j)%Miller, " Fsurfext ", fsurfext(i,j,1:3), &
!      !            " Fsurfext norm", fsurfext_norm(i,j)
!
!    enddo
!  enddo
!
! #endif

allocate(bloquer(NbCvxDom))

! The segment microstructure is tested to check if segments are cutting barriers
! RognSeg=.false.   ! Debug
if (RognSeg) then ! for a restart there is no need to cut segments  RognSeg=.false.

  !pour chaque segment on va regarder ou il coupe les barrieres plane et on
  !va reduire sa norme si l origine et l extremite sont de part et d autre
  !d une barriere
  open(555,file='../out/debug/debug_freesurf',STATUS="REPLACE")
  ! Segments are tested one by one

  do I = 1, NSEGM

    ! This segment is useless
    if(seg(i)%voiso == i) CYCLE

    ! detection of the special format :
    ! when any of the segment neighbors is < 0, this means that the file was generated my microconf
    ! and that the references to neighbor should be updated
    if(seg(i)%voiso < iun .or. seg(i)%vnno < iun .or. &
       seg(i)%voise < iun .or. seg(i)%vnne < iun ) then

      if (seg(i)%voiso < izero .or. seg(i)%vnno < izero .or. &
          seg(i)%voise < izero .or. seg(i)%vnne < izero ) then

        ! The neighbors list is of implicit type
        FlagListType = .true.

        if(seg(i)%voiso>0 .or. seg(i)%vnno>0 .or. seg(i)%voise>0 .or. seg(i)%vnne>0) then
          print *, " ERROR 1 the neighbours of I =",I, " are incompatible"
          stop
        endif
      endif

      ! Pinning points are defined in the input file by setting neighbor to zero
      ! but in the simulation, pinning points are the nsegmax segment

      ! The O (origin) side of segment i
      if(seg(i)%vnno == izero) then
        seg(i)%vnno = NSEGMAX
        if(seg(i)%voiso == izero) then
          seg(i)%voiso = NSEGMAX
        endif
      elseif (seg(i)%voiso == -IUN) then
        seg(i)%voiso = I - IUN
        seg(i)%vnno = I - IUN
        ! In the case where the origin of the segment has been modified
        ! to be compatible with the sub-lattice, the neighbor origin must be modified correspondingly
        itemp = I - IUN
        seg(i)%O(:) = seg(itemp)%O(:) + seg(itemp)%norme * bveclin(:,seg(itemp)%veclin)
      elseif (seg(i)%voiso < -IUN) then
        seg(i)%voiso = I + abs(seg(i)%voiso)
        seg(i)%vnno  = I + abs(seg(i)%vnno)
      elseif (seg(i)%vnne /= izero .and. seg(i)%vnne /= nsegmax) then
        print *, " ERROR 2 the neighbours of I =",I, " are incompatible"
        stop
      endif

      ! The E (extremity) side of segment i
      if(seg(i)%vnne == izero) then
        seg(i)%vnne = NSEGMAX
        if(seg(i)%voise == izero) then
          seg(i)%voise = NSEGMAX
        endif
      elseif (seg(i)%voise == -IUN) then
        seg(i)%voise = I + IUN
        seg(i)%vnne = I + IUN
      elseif (seg(i)%voise < -IUN) then
        seg(i)%voise = I - abs(seg(i)%voise)
        seg(i)%vnne  = I - abs(seg(i)%vnne)
      elseif (seg(i)%vnno /= izero .and. seg(i)%vnno /= nsegmax) then
          print *, " ERROR 3 the neighbours of I =",I, " are incompatible"
          stop
      endif

    endif

  enddo

 !***********************************************************************
 !* PART 1 : In case of a close domain the segment are cut:             *
 !*          if ExclIncl = F  the segment inside the domain are kept    *
 !*          if ExclIncl = T  the segment outside the domain are kept   *
 !***********************************************************************
  J    = -IUN
  skip = .false.

  AA1bis:do I = 1, NSEGM

    if (skip(i)) cycle

    ! Initialization
    iref      = i
    O(1:3)    = modulo(seg(iref)%O(:),modur)
    FRsource  = .false.
    itemp     = NSEGMAX + IUN
    lsegm     = IZERO
    keep      = .true.
    loop(:)   = IZERO
    if (seg(i)%vnne == NSEGMAX .and. seg(i)%vnno == NSEGMAX) FRsource = .true.

    do while (itemp /= iref .or. FRsource)

      lsegm = lsegm + IUN

      if (lsegm == iun) itemp =i

      !print *, 'itemp', itemp
      skip(itemp) = .true.

!       print *, 'lsegm', lsegm
!       print *, 'itemp', itemp
!       print *, 'iref', iref
!       print *, 'keep',keep

      if (lsegm > 9) stop 'infinite loop 00'

      if (keep) then

        if(FRsource) then
           maxlength = 0
        else
           maxlength = seg(itemp)%norme
        endif

        outside = .false.
        iiold   = IZERO

        do length = 0, maxlength !!!we test the full length of the segment

          bloquer   = .false.

          modOi(:) = real(O(:)+(length)*BVECLIN(:,seg(itemp)%veclin),DP)
          ! print *,length,modOi

          if (FRsource)  modOi(:) = modOi(:) + real(seg(itemp)%norme*BVECLIN(:,seg(itemp)%veclin),DP)

          ListCvxdomSegI(1:NbCvxdom) = ListCvxdom(1:NbCvxdom)
          if (length /= IZERO) then
            ListCvxdomSegI(1) = iiold
            ListCvxdomSegI(iiold) = 1
          endif

          call check_all_boundaries(modOi(:),outside,iiold)


          if(outside) then

            keep = .false.
            !write(379,*) "keep false",keep
            !write (379,*) "outside the defined domains"
            !write (379,*) "keep",keep
            !write (379,*) "bloq",bloquer
            !write (379,*) "outside",outside
          !      Print *, 'keep', keep, itemp
          !      Print *, 'bloquer', bloquer

            exit
       !      Print *, 'keep', keep, itemp
       !      Print *, 'bloquer', bloquer
           else

            if (length == izero) then
               seg(itemp)%dom = iiold
               iiold = iiold
            endif

           endif

        enddo!length

        O(:) = O(:)+maxlength*BVECLIN(:,seg(itemp)%veclin)

      endif  !keep

      loop(lsegm) = itemp

      if (FRsource) then
        itemp = iref
        FRsource = .false.
      else
        if (lsegm == iun .and. seg(itemp)%vnne == NSEGMAX .and. seg(itemp)%vnno == NSEGMAX) then
          FRsource = .true.
        else
          itemp = seg(itemp)%voise
        endif
      endif

    enddo  !do while

    !print *, 'loop',Loop(1:lsegm)
     !write(379,*) keep
     !write(379,*) '######################'

    if (.not. keep) then
      do jj = 1 , lsegm
        itemp=loop(jj)
        seg(itemp)%norme=IZERO
        seg(itemp)%voiso=itemp
        seg(itemp)%voise=itemp
        seg(itemp)%vnno=itemp
        seg(itemp)%vnne=itemp
        Out(itemp)=.true.
        seg(itemp)%dom=IZERO
      enddo
    endif

  enddo AA1bis

endif   !End RognSeg

!if (.not. FlagListType) call check_connec ("   End of the selection_plane subroutine    ")

!write(cmd,'("./domvtx")')
!call system(cmd)
!write(cmd, "python ../out/volume.py")
!call system(cmd)
!stop

end subroutine domain_definition

!############################################################################
!# This procedure is needed for the APB calculation, but MUST BE TESTED
!############################################################################
subroutine Selection_plane_cisaillement

implicit none

integer(kind=DPI)    :: I,NbSeg
integer(kind=DPI)    :: segvoiso,segvoise,segnorme,segveclin
integer(kind=DPI)    :: Oi(3),Ei(3)
logical              :: eliminate,CFALSE,cub(6)
real(kind=DP),dimension(3)    :: normint,IntersecO,IntersecE,Vect1O,Vect1E

integer(kind=DPI)             :: jj,inormint(3)
character(len=1)     :: carac

! First we load all the information regarding the planes building the microstructure
! Such information will be used in this subroutine, but mostly in 'Barriere_plane'

file_b_plan = "../in/b_plan_def"
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_b_plan    =   file_b_plan(1:len_trim(file_b_plan))//"2"
endif

! The microstructure input data file
open(44,file=file_b_plan,STATUS='OLD')

! Comments in the file header are bypassed
 carac = 'x'
do while (carac /= "#")
    read (44,*) carac
enddo

! general information about the microstructure
read(44,*) longminseg    ! facteur determinant la longeur minimum des segment bloques qu il faut discretiser
read(44,*) NbplanDom     ! nb plan delimitant la zone close cf explication dans ContBarrieres
read(44,*) InclExcl      ! clef d inclusion ou d exclusion  cf ContBarrieres
read(44,*) NbPlanMax     ! nbre de barrieres totales

! The list of plane equations
do I=1 ,NbPlanMax
  read(44,*) VarPlan(I)%Miller(:),Varplan(I)%Pos
enddo

read(44,*) TAU_APB
read(44,*) SIG_APB


close (44)

! ################

! The segment microstructure will be tested to check if segments are cutting barriers
! as well as if they are on the right side of the microstructure
if (RognSeg) then ! for a restart there is no need to cut segments  RognSeg=.false.

  ! Initialization
  NbSeg =  0        ! Number of segments we accept
  cfalse = .false.
  eliminate = .false.

  ! Accepted segment are saved in a simulation file
  open(444,file='../in/fichseg_int',STATUS="UNKNOWN")

  ! Segments are tested one by one
AA: do I = 1, NSEGM
    ! Initial segments for the mdc simulations must be of edge type

    segveclin = seg(I)%veclin
    if(TYSEG(segveclin) /= 3) stop "Input segments for the MDC simulation must only be edge segments !!!"

    ! Initial segments for the mdc simulations must be sources
    segvoiso = seg(I)%voiso
    segvoise = seg(I)%voise
    if(segvoiso /= 0 .or. segvoise /= 0)    &
      stop  "Input segments for the MDC simulation must only be source segments !!!"

    ! In order to find a segment length in the good region we first look for an acceptable
    ! starting point fixing the begining of the segment, then we look for the endding point
    segnorme = seg(I)%norme

    if (segnorme < 400) cycle AA

    Oi(:)  = Seg(I)%O(:)
    Ei(:)  =Seg(I)%O(:)+((seg(I)%norme)*(BVECLIN(:,seg(I)%veclin)))
!    Oi(:)  = modulo(Oi(:),modur(:))
!    Ei(:)  = modulo(Ei(:),modur(:))

    cub(:) = .false.

    do jj= 1,nbplanDom
      inormint(:)= Varplan(jj)%miller(:)
      normint(:) = real(inormint(:),DP)

      IntersecO(:) = InterPlanSeg(normint(:),Varplan(jj)%pos,Real(Oi,DP),inormint(:))
      IntersecE(:) = InterPlanSeg(normint(:),Varplan(jj)%pos,Real(Ei,DP),inormint(:))
      ! Normale projection vector
      Vect1O(:) = (real(Oi(:),DP)-IntersecO(:))
      Vect1E(:) = (real(Ei(:),DP)-IntersecE(:))

      ! The segment is in front of this barrier
      if ((dot_product(Vect1O,normint) > zero) .or. (dot_product(Vect1E,normint) > zero)) cub(jj)= .true.

      enddo

    if (all(cub) .eqv. InclExcl) cycle AA

    ! The effective length of segment i is now defined we can save the result

    nbseg = nbseg + 1
    write(444,*) nbseg,Oi(:),segnorme,segveclin,izero,izero,izero,izero,cfalse,izero,cfalse

  enddo AA

  ! The new number of segments
  nsegm = Nbseg

  ! The modified segment microstructure is loaded
  rewind (444)
  do I=1,nsegm

    read (444,*) nbseg,seg(i)%O(1:3),seg(i)%norme,seg(i)%veclin,seg(i)%voiso,seg(i)%vnno, &
                 seg(i)%voise,seg(i)%vnne, seg(I)%JONC,seg(i)%Ijonc,seg(i)%gd
  enddo

  close (444,status="DELETE")

endif   !End RognSeg

end subroutine Selection_plane_cisaillement

!####################################################################################
!> Upload of the info defined in the file 'eulerangle_dat' for multi-grains simulations #
!####################################################################################
subroutine DesorientMat

implicit none

Integer (kind=DPI)          :: i,j,k,l,ll,kk
real(kind=DP),allocatable   :: teta(:,:)
real(kind=DP)               :: Zgrain(3)                   !< The grain tensile axis
real(kind=DP)               :: SigPrim2(3,3)               !< A tmp stress tensor
real(kind=DP)               :: SchmidSysPrim2(NTSG)        !< A tmp SchmidSys vector
real(kind=DP)               :: nrv,nrv1,rv(3),rv1(3),rvrot(3),rv1rot(3)
character(LEN=1)            :: emile

open(444,file='../in/mM_Data_Files/eulerangle_dat',STATUS='OLD')

! We look for the beginning of Data
emile = 'x'
do while (emile /= "#")
    read (444,*) emile
enddo

! The number of grains
read(444,*) nbgrains

! If nbgrains=0 no need to go further
if (nbgrains == izero) return

! Allocations
allocate (SigAppOriente(nbgrains,3,3))
allocate (teta(nbGrains,3))
allocate (MatRotGrain(nbGrains,3,3))
allocate (TensRotSysGrain(nbgrains,NTSG,3,3))
allocate (AIRESYSGRAIN(nbgrains,NTSG))
allocate (RAUGRAIN(Nbgrains))
allocate (rauGRAINsys(Nbgrains,NTSG))
allocate (RAUGRAINjonc(Nbgrains))
allocate (rausysjoncGRAIN(Nbgrains,NTSG))
allocate (TensRotGrain(nbgrains,3,3))
allocate (Gamsysgrain(Nbgrains,NTSG))
allocate (SchmidSysGrain(Nbgrains,NTSG))

! Initialization
DesorientGrain            = .true.  ! loading in the grains is rotated
SigAppOriente(:,:,:)      = zero
TensRotSysGrain(:,:,:,:)  = zero
RAUGRAIN(:)               = zero
RAUGRAINjonc(:)           = zero
RAUGRAINsys(:,:)          = zero
RAUsysjoncGRAIN(:,:)      = zero
TensRotGrain(:,:,:)       = zero
SchmidSysGrain(:,:)       = zero
if (sideja /= iun) AireSysGRAIN(1:nbgrains,1:NTSG) = zero ! not to be initialized with restart

! The Euler angles in the grains
do i=1,nbGrains
  read(444,*) teta(i,:)
  ! print*,i,'teta',teta(i,:)
enddo

close(444)

!*******************
! The rotation matrix definition
do i=1,nbgrains
  MatRotGrain(i,1,1)=(-sin(teta(i,3))*cos(teta(i,2))*sin(teta(i,1))+cos(teta(i,3))*cos(teta(i,1)))
  MatRotGrain(i,1,2)=(-sin(teta(i,3))*cos(teta(i,1))-sin(teta(i,1))*cos(teta(i,2))*cos(teta(i,3)))
  MatRotGrain(i,1,3)=(sin(teta(i,2))*sin(teta(i,1)))
  MatRotGrain(i,2,1)=(cos(teta(i,2))*sin(teta(i,3))*cos(teta(i,1))+cos(teta(i,3))*sin(teta(i,1)))
  MatRotGrain(i,2,2)=(cos(teta(i,1))*cos(teta(i,2))*cos(teta(i,3))-sin(teta(i,3))*sin(teta(i,1)))
  MatRotGrain(i,2,3)=(-sin(teta(i,2))*cos(teta(i,1)))
  MatRotGrain(i,3,1)=(sin(teta(i,3))*sin(teta(i,2)))
  MatRotGrain(i,3,2)=(cos(teta(i,3))*sin(teta(i,2)))
  MatRotGrain(i,3,3)=(cos(teta(i,2)))
  ! print*,i,'MatRotGrain',MatRotGrain(i,:,:)
enddo

!*******************
! Calculation of the elementary rotation tensor for each slip system and for each grain
!  TensRotSysGrain = MatRot-1 * TensRotsys * MatRot
do i=1,Nbgrains
  do l=1,NTSG

    rvrot(:)=0.D0
    rv1rot(:)=0.D0
    ll=(l-IUN)*nbasered+IUN    ! Index of the first vector of each slip system
    rv(:) =real(bveclin(:,ll),DP)
    rv1(:)=real(bvecnor(:,ll),DP)

    do k=1,3
      do kk=1,3
        rvrot(k)=rvrot(k)+rv(kk)*matrotgrain(i,k,kk)
        rv1rot(k)=rv1rot(k)+rv1(kk)*matrotgrain(i,k,kk)
      enddo
    enddo

    nrv =norvect(rvrot)
    nrv1=norvect(rv1rot)

    do j=1,3
      do k=1,3
        TensRotSysgrain(i,l,j,k) = half*((rvrot(j)*rv1rot(k))-(rvrot(k)*(rv1rot(j))))/(nrv*nrv1)
      enddo
    enddo

  enddo
enddo

! Calculation of the initial Schmid factor for the slip systems in each grains.
do i=1,Nbgrains

  ! Calculation of the grain tensile axis
  do j = 1,3
    ZGrain(j) = MatRotGrain(i,j,1)*Z(1) +   &
                MatRotGrain(i,j,2)*Z(2) +   &
                MatRotGrain(i,j,3)*Z(3)
  enddo

  ! Calculation of the stress tensor in the grain
  call  ORIENTDIEG(un,ZGrain,SigPrim2,.TRUE.)

  ! Calculation of the Schmid factor
  call GlidCompPK(SigPrim2, nbase, ntsg, VecNorLin, VecNorNor, SchmidSysPrim2)
  SchmidSysGrain(i,:) = SchmidSysPrim2(:)

  !  print*,'I=',I,'MatRotGrain=',MatRotGrain(i,:,:)
  !  print*,'I=',I,'Sigapp=',Sigapp(:,:)
  !  print*,'Z=',Z(:)
  !  print*,'ZGrain=',ZGrain(:)
  !  print*,'I=',I,'SigPrim2=',SigPrim2(:,:)
  !  print*,'I=',I,'SchmidSysGrain(i,:)=',SchmidSysGrain(i,:) ; read(*,*)

enddo


deallocate (teta)

end subroutine DesorientMat

!####################################################
!# cas du pavage:                                   #
!# subroutine de verification et de correction      #
!# chaque segment doit faire parti d un grain       #
!# si tel n est pas le cas ( suite a une mauvaise   #
!# initialisation ou autre..) la subroutine lui     #
!# donne pour grain celui de son voisin             #
!####################################################
subroutine CorrectifGrain (indice,ref)

Implicit none

Integer(kind=DPI), intent (in)  :: indice
Integer, intent (in)            :: ref

if (seg(indice)%grain /= IZERO) return      ! No problem we stop here

if (seg(indice)%dom /= IZERO) then

  seg(indice)%grain = seg(indice)%dom

else

  if (kk > 1) print*,'Attention the grain of segment',indice,'is not defined',ref

  if (seg(seg(indice)%Voiso)%grain/=IZERO) then
    seg(indice)%grain=seg(seg(indice)%Voiso)%grain
  else
    seg(indice)%grain=seg(seg(indice)%Voise)%grain
  endif

endif

! seg(indice)%bloquer=seg(seg(indice)%Vnno)%bloquer

end subroutine CorrectifGrain

!#############################################################################
!> \ingroup creep
!> The particles microstructure is defined. The latter can be:
!! - loaded from the pre-existing file '../in/particules'
!! - or generated with the microstructure parameters defined in the file ../in/defparticules \n
!! Important lists for the particles definition and properties are allocated
!! and initialized. Also important lists needed in the creep loading mode are allocated and initialized.
!#############################################################################

subroutine generer_particules

implicit none

integer, parameter  :: Npmax = 10
real (kind=DP)      :: rand,density(Npmax),tauc(Npmax),perco(Npmax),Diameter(Npmax),tmp
integer (kind=DPI)  :: isys,i,nsystem_par,numero,seed(Npmax),npar_sys(npmax)
integer(kind=DPI)   :: centre(3),distmin, npar_ini,indiceP,parI,Oi(3),Tr(3),tire(3),max_par_radius
logical             :: testpar(npmax)

V_particules = zero
max_par_radius = izero
print *, ""
print *, "=============   particles information   ============="
print *, ""

file_defparticules = '../in/particules_def'

! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 2) then
  file_defparticules    =   file_defparticules(1:len_trim(file_defparticules))//"2"
   print *,file_defparticules
endif

file_particules = '../in/particules'

! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 2) then
  file_particules = file_particules(1:len_trim(file_particules))//"2"
endif

! We first check the computation status in the defparticules file
open (3,FILE=file_defparticules,STATUS='OLD')

read(3,*) Nsystem_par  ! nombre de systemes de particules

! The simulation take a pre-existing particules file
if (Nsystem_par == izero) then

  print *, " ======================================================================"
  print *, "          data on particles are not generated in the code"
  print *, "             They are read from file : particules "
  print *, " ======================================================================"

  open (303,File=file_particules,status="old")

  read (303,*) npar

  allocate (par(Npar))
  allocate (IboiteP(Npar))
  if (creep) then
    size_par_tab = 4*Npar
    if (size_par_tab < 100) size_par_tab = 100
    allocate (par_cut(size_par_tab,5))   ! The maximum number of cutting plane per particles is guessed to 4
    allocate (par_touch(size_par_tab,5)) ! The maximum number of touching plane per particles is guessed to 4
    allocate (h_par_touch(size_par_tab))
    allocate (angle_par_touch(size_par_tab))
    allocate (tau_par_touch(size_par_tab))
  endif

  do parI = 1 , npar
    read(303,*) Par(parI)%C(1:3),par(parI)%tau,par(parI)%R
       if(max_par_radius < par(parI)%R) max_par_radius = par(parI)%R
  enddo
  close (303)

  write(*,'(" Nombre of particles read in the file : = ",I10 )') npar
  indiceP = npar

else
  ! A particules file must be generated
  write(*,'(" Number of particle families = ", I5)') Nsystem_par

  do isys = 1 , Nsystem_par
       read(3,*) numero
    if (isys /= numero) call mes("ERREUR 1 du fichier biphase")
       read(3,*) TestPar(isys)  ! Key for testing. If True only one particle is generated in the center
       read(3,*) seed(isys) ! initiation of random number generation
       read(3,*) density(isys) ! number density in m-3
       read(3,*) Diameter(isys) ! diameter en nm
       read(3,*) TauC(isys)  ! maximum resistance in GPa
       read(3,*) perco(isys)  ! maximum resistance in GPa
    write(*,'(" Particle family number : ",I5)') isys
       if (TestPar(isys)) then
          print *, " TEST Mode : Only one particle is generated in the box center"
          write(*,'(" Diameter = ", F10.2, " nm")') diameter(isys)
          write(*,'(" Critical shear resistance stress = ", F6.2, " GPa")')  TauC(isys)
          npar_sys(isys) = IUN
       else
          write(*,'(" Diameter = ", F10.2, " nm")') diameter(isys)
          write(*,'(" Density = ",ES10.2E2, " m-3")')  Density(isys)
          write(*,'(" Percolation distance = ", F10.2, " nm")')  perco(isys)
          tmp =  un/DSQRT(density(isys)*1.0D-27* Diameter(isys))
          write(*,'(" Average planner spacing = ", F10.2, " nm")') tmp
          write(*,'(" Number of intersected particles = ", I10)') &
               INT(modur(2)*modur(3)*avalue*avalue/(tmp*tmp*1E-18))
          write(*,'(" Critical shear resistance stress = ", F6.2, " GPa")')  TauC(isys)
    density(isys) = density(isys) * (avalue**3)
          perco(isys) = perco(isys) * 1.E-9 / avalue
          npar_sys(isys) = INT(density(isys) * volume,DPI)
          write(*,'(" Number of particles to be generated = ", I8)') npar_sys(isys)
    endif
       Diameter(isys) =    Diameter(isys) * 1.0D-9 / avalue
       TauC(isys) = TauC(isys) * 1.0D9
       print *, "------------------------------------------------------------------------------"
  enddo

  close(3)

  !=================================
  ! the file particules is generated

  centre(:) = modur(:)/2
   ! The total number of particles
   Npar  = sum(npar_sys(1:nsystem_par))

   write(*,'(" TOTAL Nombre of particles to be generated = ", I10)') npar
   if (npar > npar_max) then
      print *, "depassement du nombre Max de particles : ",npar_max
      stop
      endif
  ! dynamical allocation of the particles array
  allocate (par(Npar))
  allocate (IboiteP(Npar))
  if (creep) then
    size_par_tab = 4*Npar
    if (size_par_tab < 100) size_par_tab = 100
    allocate (par_cut(size_par_tab,5)) ! The maximum number of cutting plane per particles is guessed to 4
    allocate (par_touch(size_par_tab,5)) ! The maximum number of touching plane per particles is guessed to 4
    allocate (h_par_touch(size_par_tab))
    allocate (angle_par_touch(size_par_tab))
    allocate (tau_par_touch(size_par_tab))
  endif

   tmp = (volume / npar)**(tiers)
   write(*,'(" Average 3D distance between particules : ",F10.2," nm")') tmp * avalue * 1E+9

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !============  definition of percolation distance ==================
   ! we compute the free spacing BETWEEN particles
   tmp = tmp - maxval(diameter(1:nsystem_par))
   write(*,'(" Average FREE 3D spacing between particules : ",F10.2," nm")') tmp * avalue * 1E+9

  IndiceP = izero
  do isys = 1 , Nsystem_par
    ! remetre les rayons d'inetraction en avalue
    npar_ini = IndiceP
    if(testpar(isys)) then
      indiceP = indiceP + IUN
      par(IndiceP)%tau  = tauc(Isys)
      !centre de la boite
          par(IndiceP)%R =   nint(diameter(isys)*half,DPI)
          if(max_par_radius < par(indicep)%R) max_par_radius = par(indicep)%R
      V_particules = V_particules + 4.0/trois * PII * (par(IndiceP)%R)**3
      par(IndiceP)%C = Centre
          !       par(IndiceP)%C(:) = Centre(:)

      ! if we want to put the par shifted / slip plane then ...
      OI(:) = BVECNOR(:,1)
      tmp = real(par(IndiceP)%R,DP)/DSQRT(real(Oi(1)**2,DP) + real(Oi(2)**2,DP) + real(Oi(3)**2,DP))
      ! tmp is the size of the precipitate in the unite of length of bvecnor
      par(IndiceP)%C(:) =  NINT(Centre(:) + Oi(:) * tmp * 0.0 )
      !       print *, tmp, par(IndiceP)%R,   NINT(BVECNOR(:,1) * (par(IndiceP)%R / tmp) / 5.0 )
      OI(:) = par(IndiceP)%C(:) - Centre(:)
      !       print *, DSQRT(real(Oi(1)**2,DP) + real(Oi(2)**2,DP) + real(Oi(3)**2,DP))

    else

      ! number of particules of system I
      if(npar_sys(Isys) == izero) then
         print *, " For particle family ", isys, " the particle number to be generated is ZERO"
         cycle
      endif
      ! loop for the generation of particules of the system Isys
      numero = izero
      distrib : do while (IndiceP < npar_ini + npar_sys(isys))
          numero = numero + iun
          if (Numero > 1000000) then
             print * , "The number of generated particles is =", indiceP
             print * , "Cannot find distribution respecting the percolation distance "
             print * , " Try to decrease the percolation distance. Sorry "
             stop
          endif

          ! random pic of the center of the particles
          do npar = 1, seed(isys)
            call random_number(rand)
          enddo
          call random_number(rand)
          rand = real(rand,DP) * real(modur(1),dp)
          Tire(1) = nint(rand)
          call random_number(rand)
          rand = real(rand,DP) * real(modur(2),dp)
          tire(2) = nint(rand)
          call random_number(rand)
          rand = real(rand,DP) * real(modur(3),dp)
          tire(3) = nint(rand)

          ! definition
          Tr(:) = Centre(:) - tire(:)
          ! problem of percolation: checkout overlap
          overlap : do parI = 1 , IndiceP
              distmin =  int(perco(isys))
              Oi(1:3) = modulo(Tr(1:3) + par(parI)%C(1:3),modur(1:3)) - Tr(1:3)
              if(norivect(tire-Oi) < distmin)  cycle distrib
          enddo overlap

          ! no overlap: the obstacle is selected
          indiceP = indiceP + IUN
          par(indiceP)%C(:) = tire(:)
          par(IndiceP)%R    = nint(diameter(isys)*half,DPI)
          if(max_par_radius < par(indiceP)%R) max_par_radius = par(indiceP)%R
          ! Determination of the obstacle size
          par(IndiceP)%tau = tauc(Isys)
          V_particules = V_particules + 4.0/trois * PII * (par(IndiceP)%R)**3
                ! write(*,'(3(2x,I8),2X,E9.1,2x,I7)') Par(indicep)%C(1:3), par(indicep)%tau , par(indicep)%R

      enddo distrib
    endif

    write(*,'("Cumulated volume fraction of particle : ",F6.3," %")') V_particules / Volume * 100E0
    !dist_moyenne_par = V_particules / npar / quatre * trois / PII
    !dist_moyenne_par = (dist_moyenne_par)**(un/trois)
    !dist_moyenne_par = UN/DSQRT(dist_moyenne_par * npar / volume)
    !write(*,'("Average planar distance between particules : ",F6.1," nm")')  dist_moyenne_par * avalue *1E+9

  enddo

endif

print *, ".................................................................."
npar = indiceP
if (npar /= IUN) write(*,'(" TOTAL NUMBER of generated particle= ",I10 )') npar
if(npar == izero) then
   stop " Erreur: microstructure: no obstacle detected"
elseif (npar> npar_max) then
   print *, "depassement du nombre Max de particles : ",npar_max
   stop
endif

write(*,'("Maximum particule radius : ",F8.1," nm")') max_par_radius * deux * avalue * 1E9

if (Microboucle >  max_par_radius / 8 ) then
   Microboucle = max_par_radius / 8
   write (*,'("*** Changing the Size of colapsing loops in GENERATION_Particules : ", I6," a = ",F5.2," nm")') &
        Microboucle, Microboucle * avalue * 1D9
endif

! The particles microtructure generated is saved in file_particules
open (3,File=file_particules,status="UNKNOWN")

write(3,*) "  ",npar
do i = 1, npar
    write(3,'(3(2x,I10),2X,E15.5,2x,I7)') par(i)%C(1:3), par(i)%tau , par(i)%R
enddo

close(3)

if (creep) then
  ! any particles can be cut in the initial state
  npar_cut = izero
  npar_touch = izero
  par_cut(:,:) = -iun
  par_touch(:,:) = -iun
  h_par_touch(:) = -un
  angle_par_touch(:) = 100.
  tau_par_touch(:) = zero
endif

end subroutine generer_particules

!###########################################################################
!# subroutine for the generation of prismatic dislocation loops            #
!# formed of one slip system with one of its collinear systems.            #
!# data has to be provided in file ../in/loops.txt                         #
!###########################################################################
subroutine generation_loops

implicit none

integer, parameter  :: Npmax    = 10
integer, parameter  :: NDL_max  = 100000
real (kind=DP)      :: rand,density(Npmax),sizedl(Npmax),temp,plansys1(3),volume_layer,proj
integer (kind=DPI)  :: sys1(Npmax), sys2(Npmax),Nloops_ini(Npmax),Lseg1, Lseg2, thickness
integer (kind=DPI)  :: isys,nsystem_DL,i1,i2,i3,i4, facteur_loop
integer (kind=DPI)  :: centre(3),percolation,Oi(3),Ei(3),count,seed
integer (kind=DPI)  :: centreDL(NDL_max,1:3),itemp,iloop,nloops,Vcirculation(3), Vburgers(3)
logical             :: testDL(npmax),repeat,layer

! the layer key indicates to the code to distribute loops only within a small layer
! passing by the center and paralell to the slip plane of the first system defined in loop
!layer = .true.
layer = .false.
centreDL(1:NDL_max,1) = izero
centreDL(1:NDL_max,2) = izero
centreDL(1:NDL_max,3) = izero

print *, ""
print *, "======================== loops information in mode Sideja = 3  ======================"
print *, ""

! We first check the computation status in the loops file
open (3,FILE="../in/loops_def",STATUS='OLD')

read(3,*) Nsystem_DL  ! nombre de systemes de loops

! The simulation take a pre-existing loops file
if (Nsystem_DL == izero) then
  stop " Erreur : number of loop systems = 0"
else
  ! A loops file must be generated
  write(*,'(" Number of loop families = ", I5)') Nsystem_DL
endif

read(3,*) seed  ! initialisation of random function

do isys = 1 ,Nsystem_DL
  read(3,*) TestDL(isys)    ! A key to generate only one particle at the center of the simulation
  read(3,*) sys1(isys)      ! Density of particles (m-3)
  read(3,*) sys2(isys)      ! Density of particles (m-3)
  read(3,*) sizedl(isys)    ! Radius of the particle(s) (nm)
  read(3,*) density(isys)   ! Density of particles (m-3)
enddo

close(3)

! initialization of the random function
do isys = 1 , seed
  call random_number(rand)
enddo

! percolation is the minimum distance between random loops
Thickness   = izero
percolation = izero

! facteur_loop is the ratio between the layer thinkness and the biggest loop size
Facteur_loop = 4

do isys = 1 , Nsystem_DL
  write(*,'(" Loop family number : ",I5)') isys
  write(*,'(" Density = ", E9.2, " m-3")')  Density(isys)
  proj = UN/dsqrt(Density(isys)*sizeDL(isys)*1.0D-9)*1.D9
  write(*,'(" Average and FREE spacing of loops = ",F7.2,"  and  ",F7.2," nm")') proj, proj - sizeDL(isys)
  density(isys) = density(isys) * (avalue**3)
  write(*,'(" Loop size = ", F7.2, " nm")') sizeDL(isys)
  sizeDL(isys) =  SizeDL(isys) * 1.0D-9 / avalue

  if(layer) then
    if(facteur_loop * sizeDL(isys) > thickness)  thickness = int(facteur_loop * sizeDL(isys))
    write(*,'("Loops are distributed randomly within a layer surrouning the box center")')
    write(*,'("Thickness of the layer for loop distribution : ",I5," nm")') nint(thickness * avalue * 1E9)
  else
    write(*,'("Loops are distributed randomly within the simulation box")')
  endif

  if(percolation > int(deux * sizeDL(isys))) percolation = int(sizeDL(isys)) * ideux
  ! Since loops are of the same nature as loops removed from the code in subroutine NET,
  !we have to reduce microboucle to a size lower than the length of prismatic loops.
  !otherwize, all generated loops will be removed.

  if (Microboucle >  sizeDL(isys) * itrois ) then
    Microboucle = int(sizeDL(isys)) * itrois
    write (*,'("  Changing Dimension of colapsing loops in GENERATION_LOOPS : ", I6," a = ",F5.2," nm")') &
        Microboucle, Microboucle * avalue * 1D9
  endif

  Vburgers(:) = reduire_vec(bveclin(:,assoc(sys1(isys),1)))
  Vcirculation(:)=  reduire_vec(prodivec(bveclin(:,sys1(isys)),bveclin(:,sys2(isys))))
  write(*,'("Burgers Vector of the loop    : [",3I3," ]")') Vburgers(1:3)
  write(*,'("Circulation vector of the loop: [",3I3," ]")') Vcirculation(1:3)
  print *, "========================================"

  if(etat(Vburgers,Vcirculation) == 3) then  ! vectors are equal
    print *, " Detected Type of the loop : VACANCY loop"
  elseif(etat(Vburgers,Vcirculation) == -3) then ! vectors are opposite
    print *, " Detected Type of the loop : INTERSTITIAL loop"
  else
    stop " Impossible detection of the Type of the loop. SORRY"
  endif

  print *, "========================================"

enddo

! nsegp is the current nomber of segments introduced beyond NSEGM
centre(:) = modur(:)/ideux
PlusSeg = 0
Nloops = 0
! normalized Miller indices of slip plan for the first slip sys
plansys1(:) = vecnornor(:,sys1(1))

if(layer) then

  if(modur(1) /= modur(2) .or. modur(1) /= modur(3)) then
    print *, " Procedure of loops generation within thin layer does not work with non cubic simulation box"
    print *, " Use a cubic simulation box. Sorry"
    stop
  endif

  if (BCC) then
    !surface of a 110 plan passing by the center is  1*sqrt(2)
    volume_layer = dsqrt(deux) * modur(1) * modur(1) * thickness
  elseif(CFC) then
    !surface of a 111 plan passing by the center is  3sqrt(2)/2 * sqrt(6)/4
    volume_layer = 3*dsqrt(deux)*half * dsqrt(6.0D0)/quatre * modur(1) * modur(1) * thickness
  else
    print *, "Undefined volume of loop layer for this crystallographic structure"
  endif

  write(*,'("Volume of the loop layer  : ",F10.1," nm3")') volume_layer * avalue**3 * 1E27

else

  volume_layer = volume

endif

do isys = 1 , Nsystem_DL

  nloops_ini(isys) =  nint(density(isys) * volume_layer,DPI)

  if (testDL(isys))then
    print *, " Loop TEST : only ONE loop is introduced in the box center "
  else
    if(nloops_ini(isys) /= izero) then
      write(*,'(" Expected number of loops for this system : ", I8)') nloops_ini(isys)
    else
      stop " Loop density too low"
    endif
  endif

enddo

do isys = 1 , Nsystem_DL
  print *, " Number of segments befor isys :", isys," is :",nsegm
  Lseg1 = int(sizedl(isys)/normlin(sys1(isys)))
  Lseg2 = int(sizedl(isys)/normlin(sys2(isys)))
  if (testDL(isys))then
    Nloops = Nloops + 1
    ! fr the position of the test loop, two choices. Comment/uncomment the specific line as needed
    !center of the loop in the center of the box
    !Oi(:) = Centre(:) - (Lseg1 * bveclin(:,sys1(isys)) + Lseg2 * bveclin(:,sys2(isys)))/IDEUX
    Oi(:) = Centre(:)
    Oi = noeud(Oi,IUN)
    centreDL(Nloops,:) =Oi(:)
    ! segment 1 is in the center of th box. Oi is the origin of the first segment of the loop
    ! for segments will be introduced
    i1 = nsegm + 1
    i2 = i1 + 1
    i3 = i2 + 1
    i4 = i3 + 1
    nsegm = nsegm + 4
    ! data for the first segment
    seg(i1)%O(:) = Oi(:)
    seg(i1)%veclin = sys1(isys)
    seg(i1)%norme= Lseg1
    seg(i1)%voiso = i4
    seg(i1)%voise = i2
    seg(i1)%vnno  = i4
    seg(i1)%vnne  = i2
    seg(i1)%surface = izero
    seg(i1)%varfreeplan = izero
    if (nbcvxdom > IUN) then
     call assign_domain(i1,3351)
    else
    seg(i1)%dom = IUN
    endif
    Ei(:) = seg(i1)%O(:) + Lseg1 * Bveclin(:,sys1(Isys))
    ! data for the second segment
    seg(i2)%O(:) = Ei(:)
    seg(i2)%veclin = sys2(isys)
    seg(i2)%norme= Lseg2
    seg(i2)%voiso = i1
    seg(i2)%voise = i3
    seg(i2)%vnno  = i1
    seg(i2)%vnne  = i3
    seg(i2)%surface = izero
    seg(i2)%varfreeplan = izero
    if (nbcvxdom > IUN) then
     call assign_domain(i2,3352)
    else
    seg(i2)%dom = IUN
    endif
    Ei(:) = Ei(:) + Lseg2 * Bveclin(:,sys2(Isys))
    ! data for the 3d segment
    seg(i3)%O(:) = Ei(:)
    if (modulo(sys1(isys),nbasered) < 5) then
      seg(i3)%veclin = sys1(isys) + 4
    else
      seg(i3)%veclin = sys1(isys) - 4
    endif
    seg(i3)%norme= Lseg1
    seg(i3)%voiso = i2
    seg(i3)%voise = i4
    seg(i3)%vnno  = i2
    seg(i3)%vnne  = i4
    seg(i3)%surface = izero
    seg(i3)%varfreeplan = izero
    if (nbcvxdom > IUN) then
     call assign_domain(i3,3353)
    else
    seg(i3)%dom = IUN
    endif
    Ei(:) = Ei(:) - Lseg1 * Bveclin(:,sys1(Isys))
    ! data for the 4d segment
    seg(i4)%O(:) = Ei(:)
    if (modulo(sys2(isys),nbasered) < 5) then
      seg(i4)%veclin = sys2(isys) + 4
    else
      seg(i4)%veclin = sys2(isys) - 4
    endif
    seg(i4)%norme= Lseg2
    seg(i4)%voiso = i3
    seg(i4)%voise = i1
    seg(i4)%vnno  = i3
    seg(i4)%vnne  = i1
    seg(i4)%surface = izero
    seg(i4)%varfreeplan = izero
    if (nbcvxdom > IUN) then
     call assign_domain(i4,3354)
    else
    seg(i4)%dom = IUN
    endif

    call connecIJ(i1,i2, -11)
    call connecIJ(i2,i3, -12)
    call connecIJ(i3,i4, -13)
    call connecIJ(i4,i1, -14)

    Nsegm = nsegm + Plusseg
    print *, " Number of segments generated during loop test :", nsegm
    PlusSeg = 0

  else

    do iloop = 1,Nloops_ini(isys)
      repeat = .true.
      count = 0
      do while(repeat)
        ! we contunue to select random position of the loop while it overlaps other loops
        call random_number(rand)
        rand = real(rand,DP) * real(modur(1),dp)
        Oi(1) = nint(rand)
        call random_number(rand)
        rand = real(rand,DP) * real(modur(2),dp)
        Oi(2) = nint(rand)
        call random_number(rand)
        rand = real(rand,DP) * real(modur(3),dp)
        Oi(3) = nint(rand)
        Oi = Noeud(Oi,iun)
        if (layer) then
          proj = (Oi(1)-centre(1))*plansys1(1) + (Oi(2)-centre(2))*plansys1(2)+(Oi(3)-centre(3))*plansys1(3)
          ! if the loop is out of the layer wr do not consider it
          if(2 * abs(proj) > thickness)  cycle
        endif
        count = count +1
        ! problem of percolation: checkout overlap
        repeat = .false.
        overlap : do itemp = 1 , Nloops
           temp = norivect(Oi(:)-centreDL(itemp,:))
           if(temp < percolation) then
              repeat = .true.
              exit
           endif
        enddo overlap
        if (repeat) cycle
        if (count > 1000) stop " no solution to avoid loop overlap"
        ! we know now that the loop is not overlapping others and that it is within the layer
        Nloops = Nloops + 1
        centreDL(Nloops,:) =Oi(:)
        ! segment 1 is in the center of th box. Oi is the origin of the first segment of the loop
        ! for segments will be introduced
        i1 = nsegm + 1
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 +1
        nsegm = nsegm + 4
        ! data for the first segment
        Oi(:) = Oi(:) - (Lseg1 * bveclin(:,sys1(isys)) + Lseg2 * bveclin(:,sys2(isys)))/IDEUX
        seg(i1)%O(:) = Oi(:)
        seg(i1)%veclin = sys1(isys)
        seg(i1)%norme= Lseg1
        seg(i1)%voiso = i4
        seg(i1)%voise = i2
        seg(i1)%vnno = i4
        seg(i1)%vnne = i2
        seg(i1)%surface = izero
        seg(i1)%varfreeplan = izero
        if (nbcvxdom > IUN) then
         call assign_domain(i1,3354)
        else
        seg(i1)%dom = IUN
        endif
        Ei(:) = seg(i1)%O(:) + Lseg1 * Bveclin(:,sys1(Isys))
        ! data for the second segment
        seg(i2)%O(:) = Ei(:)
        seg(i2)%veclin = sys2(isys)
        seg(i2)%norme= Lseg2
        seg(i2)%voiso = i1
        seg(i2)%voise = i3
        seg(i2)%vnno = i1
        seg(i2)%vnne = i3
        seg(i2)%surface = izero
        seg(i2)%varfreeplan = izero
        if (nbcvxdom > IUN) then
         call assign_domain(i2,3354)
        else
        seg(i2)%dom = IUN
        endif
        Ei(:) = Ei(:) + Lseg2 * Bveclin(:,sys2(Isys))
        ! data for the 3d segment
        seg(i3)%O(:) = Ei(:)
        if (modulo(sys1(isys),nbasered) < 5) then
          seg(i3)%veclin = sys1(isys) + 4
        else
          seg(i3)%veclin = sys1(isys) - 4
        endif
        seg(i3)%norme= Lseg1
        seg(i3)%voiso = i2
        seg(i3)%voise = i4
        seg(i3)%vnno = i2
        seg(i3)%vnne = i4
        seg(i3)%surface = izero
        seg(i3)%varfreeplan = izero
        if (nbcvxdom > IUN) then
         call assign_domain(i3,3354)
        else
        seg(i3)%dom = IUN
        endif
        Ei(:) = Ei(:) - Lseg1 * Bveclin(:,sys1(Isys))
        ! data for the 4d segment
        seg(i4)%O(:) = Ei(:)
        if (modulo(sys2(isys),nbasered) < 5) then
          seg(i4)%veclin = sys2(isys) + 4
        else
          seg(i4)%veclin = sys2(isys) - 4
        endif
        seg(i4)%norme= Lseg2
        seg(i4)%voiso = i3
        seg(i4)%voise = i1
        seg(i4)%vnno = i3
        seg(i4)%vnne = i1
        seg(i4)%surface = izero
        seg(i4)%varfreeplan = izero
        if (nbcvxdom > IUN) then
          call assign_domain(i4,3354)
        else
          seg(i4)%dom = IUN
        endif

        call connecIJ(i1,i2, -21)
        call connecIJ(i2,i3, -22)
        call connecIJ(i3,i4, -23)
        call connecIJ(i4,i1, -24)

        Nsegm = nsegm + Plusseg
        PlusSeg = 0

      enddo
    enddo
    print *, " Number of segments After LOOPS :", nsegm
  endif
enddo

print *, ""
print *, "====================================================  end of loop information   ============="
print *, ""

end subroutine generation_loops

!####################################################################
!# procedure to test the initial configuration, before we start     #
!####################################################################
subroutine verification_configuration

implicit none

integer (kind = DPI) :: i,ii
integer (kind = DPI) :: O(3),Li(3)

write(*, '("  Number of segments before connection :  ",I10)') Nsegm

!=====================================================================================
!============ The simulation box dimensions are controlled    ========================
!=====================================================================================

!The box dimension must be compatible with the Facteur_boite
O(:) = (modur(:) / Facteur_boite) * Facteur_boite - modur(:)

if (abs(O(1)) + abs(O(2)) + abs(O(3)) /= izero) then
   if (CFC) then
      print *, " The simulation box dimensions : ",modur
      print *, " Such box dimensions are not proportional to facteur_boite : ", facteur_boite
      stop
   else
      Modur(:) = (modur(:) / Facteur_boite) * Facteur_boite
      print *, " =================================================================================="
      write (*,'("   Attention : the simulation box size was modified : ",3I8)') modur(3)
      print *, " =================================================================================="
   endif
endif

if ((GB == itrois) .and. modur(1) /= Dimsyst) then
   print*,'modur a ete modifie'
   print*,'pour etre coherent il faut changer la variable Dimsyst', &
   'dans le fichier infonoeud et relancer le prgm pavage'
   stop
endif

! remettre les dimensions de la boite sur les noeuds du reseau
! les calcul montrent que les coordonnees de la boite doivent etre des multiples
! du facteur  Facteur_boit , defini dans le module de constantes
!call segsinfo(" avnt modulo                   ")
! si nsegm = 0, alors on appelle generer_segements pour creer une configuration
! pour l'etude des jonction

!==========================================================================================
!============= verification des dimensions de conformite==================================
! si le programme n'est pas arrete, c'est que les poits suivants :
! [ modur(1) ; 0 ; 0 ],  [ 0 ; modur(2) ; 0 ] et [ 0 ; 0 ; modur(3) ]
! sont sur le reseau de la simu
!==========================================================================================
O(1:3) = izero
O(1) = modur(1)
if (etat(O,noeud(O,IUN)) /= 3) then
   print*, " O = ",O(:)
   print*, "  noeud(O,IDEUX= ",noeud(O,IDEUX)
   print*, " load:mauvais facteur de multiplication selon X"
   stop
endif
O(1:3) = izero
O(2) = Modur(2)
if (etat(O,noeud(O,IUN)) /= 3) stop " load:mauvais facteur de multiplication selon Y"

O(1:3) = izero
O(3) = Modur(3)
if (etat(O,noeud(O,IUN)) /= 3) stop " load:mauvais facteur de multiplication selon Z"

! teste ultime
if (etat(modur,noeud(modur,IUN)) /= 3) stop " load:le point modur(1:3) n'est pas sur le reseau"

!=====================================================================================
!============ The segments distribution is controlled    ========================
!=====================================================================================

! We compt the number of sources in the initial configuration
Nombre_FR = IZERO
do i= 1,nsegm
  if (seg(i)%voiso == IZERO) nombre_FR = Nombre_FR + 1
  if (seg(i)%voise == IZERO) nombre_FR = Nombre_FR + 1
  if (seg(i)%voiso == NSEGMAX) nombre_FR = Nombre_FR + 1
  if (seg(i)%voiso == NSEGMAX) nombre_FR = Nombre_FR + 1
enddo

nombre_FR = nombre_FR / ideux  ! each source must have an initial and a ending pinning point

print*, " "
print*, "Initial configuration specificity"
write(*, '("  Number of Frank_Read sources: ",I10)') Nombre_FR

! Origin of segments is moved at correct coordinates on the simulation lattice
if(seg(1)%voiso > -4) then ! When we start from complex configuration extracted from SEG-save
                           ! files, we cannot check the connectivity and we must trust the simulation

  ! this is an initial segment configuration
  if(.not. seg_save) then

    do i = 1,nsegm

        ! The segments we must not modify
        if (seg(i)%voiso == i .and. seg(i)%voise == i) cycle    ! Segments to eliminate
        if (seg(i)%surface /= 0 ) cycle                         ! Segments touching a surface

        ! The segments we test and reconstruct if needed
        if (seg(i)%voiso == 0 .or. seg(i)%voiso == nsegmax .or. seg(i)%voiso == -3) then

          ! The segment coordinates are controlled
          O(:)        = seg(i)%O(:)
          seg(i)%O(:) = noeud(O,IUN)

          Li(:) = O(:)-seg(i)%O(:)  ! The shift vector imposed by the new segment origin location
          ii=IZERO

          ! Shift may be repeated many time to find a correct location
          do while (Li(1) /= izero .or. Li(2) /= izero .or. Li(3) /= izero)

            O(:) = O(:) + BVECLIN(:,seg(I)%veclin)
            seg(i)%O(:) =  noeud(O,IUN)

            Li(:) = O(:)-seg(i)%O(:)
            ii=ii+IUN

            if (ii > 100) then
              write(*,*) " "
              write(*,*) "The segment", i, " origin is found out of the material simulation sub-lattice !!!"
              write(*,*) "Correct location in a closed area could not be fixed..."
              Stop
            endif

          enddo

          seg(i)%norme = seg(i)%norme - ii

        endif

        ! The result of the shift is finally plotted
        O(:)  = seg(i)%O(:)
        O(:)  = noeud(O,IUN)
        LI(:) = O(:) -  seg(i)%O(:)

        if (kkdebug) then
          if(Li(1) /= izero .or. Li(2) /= izero .or. Li(3) /= izero) then
            write(379,*) " Origine du seg :",I, " n est pas sur le reseau "
            write(379,*) " seg(i)%O(:) ",seg(i)%O(:)
            write(379,*) " noeud ", O(:)
            write(379,*) " seg(i-1)%Ei(:) ",seg(i-1)%O(:) + seg(i-1)%norme * bveclin(1:3,seg(i-1)%veclin)
          endif
        endif

    enddo

  else    ! this is a restart segment configuration

    ! Pinning points are redefined with the NSEGMAX value
    ! This simple modification is needed when we wants to restart a simulation with
    ! a larger value of NSEGMAX
    do i = 1,nsegm
        if(seg(i)%voiso == izero) seg(i)%voiso = NSEGMAX
        if(seg(i)%vnno  == izero) seg(i)%vnno  = NSEGMAX
        if(seg(i)%voise == izero) seg(i)%voise = NSEGMAX
        if(seg(i)%vnne  == izero) seg(i)%vnne  = NSEGMAX
    enddo

  endif

endif

end subroutine  verification_configuration

!############################################################################
!# This procedure is called to make the appropriate connection between      #
!# segments generated in the pre-treatment or post treatment.              #
!# Microconf generate neighbors with negative or zero values               #
!# while the microstructure saved by the code is already with the accurate  #
!# neighbours.                                                              #
!############################################################################
subroutine connecxion_segments

implicit none

integer (kind=DPI)  :: i,itemp,iNew
logical             :: rough  ! key : True for segments generated by Microconf, false for saved segments

! Initialization
rough = .false.

#ifdef MDC
! For the MDC simulation the segment connectivity must always be checked
! since the data in the segment input file are modified to build-up the Voltera loop
rough = .true.
#endif

! First we modify neighbors from negative values to the real values with the
! following significations:
! 0   : pinning points (neighbor = nsegmax)
! -1  : the segment is connected to the previous / next segment in the file)
! -3  : the segment is connected to the 3d previous / next segment in the file)
itemp = IZERO


if(GB < IDEUX .or. GB == ICINQ) then !before it was if (GB == IZERO)

  do I = 1, nsegm

      ! detection of the special format :
      ! when any of the segment neighbors is < 0, this means that the file was generated my microconf
      ! and that the references to neighbors should be updated

      if(seg(i)%voiso == i) cycle     ! This segment is to eliminate, nothing to do

      ! The junction segments special case
      if(seg(i)%Ijonc > Izero)  then
        if( .not. seg(i)%jonc) then
          print *, "Error1: connection_segemnt: bad junction infortation for I=", I
          stop
        endif
      else
        if(seg(i)%jonc) then
          print *, "Error2: connection_segemnt: bad junction infortation for I =", I
          stop
        else
          seg(i)%Ijonc = nsegmax
        endif
      endif

      ! The segments we must test and reconstruct if needed
      if (seg(i)%voiso <= 0 .or. seg(i)%vnno <= 0 .or.    &
          seg(i)%voise <= 0 .or. seg(i)%vnne <= 0 )           then

        rough = .true.

        ! a test to eliminate not anticipated reconstructions
        if (seg(i)%voiso < -3 .or. seg(i)%voiso < -3 .or.   &
            seg(i)%voiso < -3 .or. seg(i)%voiso < -3)           then
          print *, " This type of initial segment configuration is not possible"
          stop
        endif

        !!!!!!!!!! The segments at the origin side

        ! Pinning point are noted by 0 or nsegmax as neighbors in the initial configuration
        ! but must be noted with nsegmax during the simulation
        if(seg(i)%voiso == izero) then

          seg(i)%voiso  = NSEGMAX
          seg(i)%vnno   = NSEGMAX

        ! The first neighbor reconstruction at the origin of segment i
        elseif (seg(i)%voiso == -1 .and. seg(i)%vnno == -1) then

          seg(i)%voiso = I + seg(i)%voiso
          seg(i)%vnno  = I + seg(i)%vnno

          ! In the case where the origin of the FR source has been modified
          ! to be compatible with the sub-lattice, the neighbor origin should be modified
          itemp       = I - IUN
          seg(i)%O(:) = seg(itemp)%O(:) + seg(itemp)%norme * bveclin(:,seg(itemp)%veclin)

        ! The first neighbor reconstruction at the origin of segment i
        elseif (seg(i)%voiso == -3) then

          seg(i)%voiso = I - seg(i)%voiso
          seg(i)%vnno  = I - seg(i)%vnno

        elseif (seg(i)%voiso == -1 .and. seg(i)%vnno /= -1) then

          seg(i)%voiso = I + seg(i)%voiso
          seg(i)%vnno  = I + seg(i)%vnno

        endif


        !!!!!!!!!! The segments at the end side

        ! Pinning point are noted by 0 or nsegmax as neighbors in the initial configuration
        ! but must be noted with nsegmax during the simulation
        if(seg(i)%voise == izero) then

          seg(i)%voise = NSEGMAX
          seg(i)%vnne  = NSEGMAX

        elseif (seg(i)%voise == -1 .and. seg(i)%vnne == -1) then

          seg(i)%voise = I - seg(i)%voise
          seg(i)%vnne  = I - seg(i)%vnne

        elseif (seg(i)%voise == -3) then

          seg(i)%voise = I + seg(i)%voise
          seg(i)%vnne  = I + seg(i)%vnne

        elseif (seg(i)%voise == -1 .and. seg(i)%vnne /= -1) then

          seg(i)%voise = I - seg(i)%voise
          seg(i)%vnne  = I - seg(i)%vnne

        endif

      endif

  enddo

else

  ! if the GB /= 0, neighbors have been modified according to procedures selection
  do I = 1, NSEGM

    if(seg(i)%voiso == izero ) seg(i)%voiso = NSEGMAX
    if(seg(i)%voise == izero ) seg(i)%voise = NSEGMAX

    seg(i)%vnno = seg(i)%voiso
    seg(i)%vnne = seg(i)%voise

    ! in the case the origin of a segment is properly set, it is connected to its neighbor
    itemp = seg(i)%voiso
    if (itemp < IUN) then
      print *, "ERROR of indexing neihbours in GB"
    endif

    if(itemp /= NSEGMAX) then
      seg(i)%O(:) = seg(itemp)%O(:) + &
      seg(itemp)%norme * bveclin(:,seg(itemp)%veclin)
    endif

  enddo

endif

! do I=1,NSEGM
!   write(*,fmt='(13I8)') i,seg(i)%O,seg(i)%O+seg(I)%norme*BVECLIN(:,seg(i)%veclin),seg(i)%veclin,      &
!                         seg(i)%norme,seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne
! enddo
! read(*,*)

! The connectivity is controlled when fichseg is not saved from simulation
if ((rough .or. (GB /= izero))) then

  PlusSeg = 0

  if (seg_save) then !Kneecap must be added to segment touching free surface in case of restart from saved configuration

    do I= 1,NSEGM

      if (seg(i)%surface==IUN .and. seg(i)%voiso==NSEGMAX) then

        PlusSeg = PlusSeg + IUN
        iNew    = Nsegm + PlusSeg

        call init_seg(Nsegm+PlusSeg)

        seg(iNew)%O       = seg(i)%O
        seg(iNew)%voise   = i
        seg(iNew)%vnne    = i
        seg(iNew)%gd      = izero
        seg(iNew)%veclin  = nbrot(seg(i)%veclin,seg(i)%veclin,2)
        seg(i)%voiso      = iNew
        seg(iNew)%dom     = seg(i)%dom

      endif

      if (seg(i)%surface==IDEUX .and. seg(i)%voise==NSEGMAX) then

        PlusSeg = PlusSeg + IUN
        iNew    = Nsegm + PlusSeg

        call init_seg(Nsegm+PlusSeg)

        seg(iNew)%O       = seg(i)%O + seg(i)%norme*BVECLIN(:,seg(i)%VECLIN)
        seg(iNew)%voiso   = i
        seg(iNew)%vnno    = i
        seg(iNew)%gd      = izero
        seg(iNew)%veclin  = nbrot(seg(i)%veclin,seg(i)%veclin,2)
        seg(i)%voise      = iNew
        seg(iNew)%dom     = seg(i)%dom

      endif

    enddo

  else

    do I=1,NSEGM

      if (seg(i)%voiso == i .and. seg(i)%voise == i) cycle

      if(seg(i)%norme == izero) then

        !Pivot segments may be present in initial configuration, but it is exeptionnal.
        print *, "Warning : found kneecap number : ",i," in the microconf generate file"

#ifdef MDC
        print *,'**********************************************************************************'
        print *,'Check the dispersion parameter in confinit_def   (maybe you have to decrease it..)'
        print *,'**********************************************************************************'
        stop
#endif

      endif

      itemp = seg(i)%vnno
      if(itemp == nsegmax) then
        call connecIJ(itemp,I, 1477)
      endif

      ! Connectivity between segment I and its %VoisE is systematically tested
      itemp = seg(i)%vnne
      if(seg(i)%norme /= 0) then
        call connecIJ(I,itemp, 1547)
      endif

    enddo

  endif

  nsegm   = nsegm + PlusSeg
  PlusSeg = 0

endif

write(*, '("  Number of segments after connection :  ",I10)') Nsegm

! do I=1,NSEGM
!   write(*,fmt='(13I8)') i,seg(i)%O,seg(i)%O+seg(I)%norme*BVECLIN(:,seg(i)%veclin),seg(i)%veclin,      &
!                         seg(i)%norme,seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne
! enddo
! read(*,*)

end subroutine connecxion_segments

!####################################################################
!# The subroutine used to read the Bigsave file. The one needed to  #
!# restart a simulation as it is.                                   #
!####################################################################
subroutine loadall

implicit none

integer(kind = DPI) :: Nstatold,i,nsegmax_old
real(kind = DP)     :: rap,deltatini
real(kind = DP)     :: EPSOLDold(2000),GAMMAOLDSYSold(NTSG_MAX,2000),SEPSOLDold,SGAMMAOLDSYSold(NTSG_MAX)
character           :: carac

!***************************************
! Some GB special cases
!***************************************
select case (GB)

case(itrois)

  open(46,file="../in/b_poly_def",STATUS='old')
  carac='x'
  do while (carac /= "#")
      read (46,*) carac
  enddo

  read(46,*) longminseg
  read(46,*) DesorientGrain
  read(46,*) DimSyst
  read(46,*) NbGrains

  close(46)

end select
!***************************************

file_bigsave  =   "../out/bigsave"
! A second set of input files is needed when Nb_phase = 2
if (Nb_phase == 2 .and. Ma_Couleur == 1) then
  file_bigsave  =   file_bigsave(1:len_trim(file_bigsave))//"2"
endif
file_bigsave  =   file_bigsave(1:len_trim(file_bigsave))//".bin"


open(99,FILE=file_bigsave,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')

Read(99) KK,nsegm,nsegmax_old,Npar,Nstatold,Sigma,EPSO,EPSDOT,DELSIG,accutime
if(Nstatold > 2000) stop " Ancien NSTATCONTROL > 2000 , problem "
Read(99) SIGAPP,SIGPLUS,DELTATini,SONDMS,VITMOY,ModuR(1:3),FSIG
Read(99) Solli_Sys(1:NTSG)
Read(99) SEG(1:nsegm)
Read(99) CranSys(1:NTSG)
Read(99) RAUSYS(1:NTSG),AireSYS(1:NTSG),GAMMADOTSYS(1:NTSG),TrAppSys(1:NTSG),TrIntSys(1:NTSG)
Read(99) RAUDMO_I,RAUDIS,RAUDMO,AireVisSys(1:NTSG),AireCoinSys(1:NTSG),AireMixtSys(1:NTSG)
Read(99) EPSOLDold(1:Nstatold),GAMMAOLDSYSold(1:NTSG,1:Nstatold)
Read(99) NTGD,NTGD1,NTGD2,OLDNTGD1,OLDNTGD2,airevis,airecoin,airemixte,Stress_fact

if (desorientgrain) Read(99) AireSysGrain(1:nbgrains,1:NTSG)

if (npar > izero) then
  allocate(par(npar))
  read(99) par(1:npar)
  particules = .true.
  allocate (IboiteP(Npar))
  if (creep) then
    size_par_tab = 4*Npar
    if (size_par_tab < 100) size_par_tab = 100
    allocate (par_cut(size_par_tab,5))    ! The maximum number of cutting plane per particles is guessed to 4
    allocate (par_touch(size_par_tab,5))  ! The maximum number of touching plane per particles is guessed to 4
    allocate (h_par_touch(size_par_tab))
    allocate (angle_par_touch(size_par_tab))
    allocate (tau_par_touch(size_par_tab))
  endif
endif

if ((.not. allocation_dynamique_boites) .and. calculate_gammabox) then
  read(99) NBoites
  allocate (Gammabox(Nboites,NTSG))   ! In order to calculate the gamma on each box
  read(99) Gammabox(1:NBoites,1:NTSG)
endif

close(99)

! The simulation is restarted with a new value of Nsegmax
! This modification must be applied to the seg array values = nsegmax
if(nsegmax_old /= nsegmax) then
  do i=1,nsegm
    if (seg(i)%voiso == nsegmax_old) seg(i)%voiso = nsegmax
    if (seg(i)%voise == nsegmax_old) seg(i)%voise = nsegmax
    if (seg(i)%vnno  == nsegmax_old) seg(i)%vnno  = nsegmax
    if (seg(i)%vnne  == nsegmax_old) seg(i)%vnne  = nsegmax
    if (seg(i)%ijonc == nsegmax_old) seg(i)%ijonc = nsegmax
  enddo
endif

#ifdef MDC

#if defined(MDC) && defined(PA)
 if (Mon_Rang ==  IZERO) then
#endif

   call send_globdata

   call recv_globdata

#if defined(MDC) && defined(PA)
 endif

  CALL MPI_BCAST(mdcvol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  ! Only the process zero of the communicator can write in output files, the others are waiting
#endif

call time_step
if (sideja==1) deltat0=deltat

#endif

! If we changed the Deltat0 or NstatControl in the control file
! values of epsold and gammaoldsys must be changed to account these modifications
rap = DELTATini/Deltat0

! print *," rap,Nstatcontrol",rap,Nstatcontrol
! deallocate(EPSOLD)
! deallocate(GAMMAOLDSYS)
! allocate (EPSOLD(NstatControl))
! allocate (GAMMAOLDSYS(NTSG,NstatControl))

SEPSOLDold = ((EPSOLDold(Nstatold)*DELTATini - EPSOLDold(1)*DELTATini)/rap)/Deltat0

do i = 1,NTSG
  SGAMMAOLDSYSold(i) = (GAMMAOLDSYSold(i,Nstatold)-GAMMAOLDSYSold(i,1))/rap
enddo

do i = 1, Nstatcontrol
    epsold(i) = (EPSOLDold(Nstatold)*rap-SEPSOLDold)+ i*SEPSOLDold/Nstatcontrol
    GAMMAOLDSYS(1:NTSG,i)= (GAMMAOLDSYSold(1:NTSG,i)-SGAMMAOLDSYSold(1:NTSG))    &
                          + i*SGAMMAOLDSYSold(1:NTSG)/Nstatcontrol
enddo

! Variables needed for the PBC
InvmoduR(1:3)   = un / Modur(1:3)
HalfModur(1:3)  = Modur(1:3) / 2
ModurMax        = maxval(Modur(1:3))

! Variables needed for mode_deformation = IHUIT
if (mode_deformation == IHUIT) then
  if (FEgrid_size(1) < IZERO .or. FEgrid_size(2) < IZERO .or. FEgrid_size(3) < IZERO) then
    FEgrid_size(:) = modur(:)
    FEgrid_or(:)=(/0,0,0/)
  endif
endif

end subroutine loadall

!############################################
!
!############################################
subroutine double_config

Integer(kind=DPI) :: dislo,vl, O(3),s,vec(3),i

PLUSSEG = 0

Do dislo = 1, nsegm
  if (seg(dislo)%voise == 0) then
    vl = seg(dislo)%veclin
    O(:) = seg(dislo)%O(:)+ seg(dislo)%norme*bveclin(:,vl)

    ! composante de decalage parallele a vis1
    Vec(:) = Bvecdep(:,assoc(vl,3))*(Nint(Ldedou/avalue/normdep(assoc(vl,3)))/4)*4
    ! on ajoute une composante de decalage parallele a coin1 ou coin2

    O(:) = O(:) + Vec(:)
    s = dislo
    do while (S /= 0)
      vl = seg(s)%veclin
      PLUSSEG = PLUSSEG + 1
      i = NSEGM + PLUSSEG
      !           print *,"s =", s,"      i=", i
      seg(i)%O = O(:)
      seg(i)%norme = seg(s)%norme
      if(seg(s)%voise == 0) then
        seg(i)%voiso = 0
      else
        seg(i)%voiso = i-1
      endif
      if(seg(s)%voiso == 0) then
        seg(i)%voise = 0
      else
        seg(i)%voise = i+1
      endif
      if(vl+4 <= nbase .and. syseg(vl) == syseg(vl+4)) then
        seg(i)%veclin = vl + 4
      else
        seg(i)%veclin = vl - 4
      endif
      s = seg(s)%voiso
      O(:) = seg(i)%O(:)+ seg(i)%norme*bveclin(:,seg(i)%veclin)
    enddo
  endif
enddo

print *, " segments initiaux =", nsegm
print *, " segments ajoute =", PLUSSEG
NSEGM = NSEGM + PLUSSEG

PLUSSEG = 0

end subroutine double_config

!##############################################
!
!##############################################
subroutine generation_microstructure

implicit none

integer,parameter   :: DPI_S=selected_int_kind(9)   ! The integer format needed for the film.bin file
integer             :: IOS                          ! Read iostat error message (0 = OK, > 0 means error)
integer(kind=dpi)   :: i
integer(kind=dpi_s) :: kompt, kompt2, komptlast
real(kind=dp)       :: dist_moyenne_par
character           :: caractest
character(len=256)  :: cmd

! initialization
SEG(1:NSEGMAX)%VECLIN   = izero
SEG(1:NSEGMAX)%O(1)     = izero
SEG(1:NSEGMAX)%O(2)     = izero
SEG(1:NSEGMAX)%O(3)     = izero
SEG(1:NSEGMAX)%NORME    = izero
SEG(NSEGMAX)%NORME      = NSEGMAX
SEG(NSEGMAX)%VOISO      = NSEGMAX
SEG(NSEGMAX)%VOISE      = NSEGMAX
SEG(NSEGMAX)%VNNO       = NSEGMAX
SEG(NSEGMAX)%VNNE       = NSEGMAX
do i = 1 ,nsegmax -1
  SEG(I)%VOISO          = I
  SEG(I)%VOISE          = I
  SEG(I)%VNNO           = I
  SEG(I)%VNNE           = I
enddo
SEG(1:NSEGMAX)%IJONC    = nsegmax
SEG(1:NSEGMAX)%RESDEP   = zero
SEG(1:NSEGMAX)%JONC     = .false.
SEG(1:NSEGMAX)%GD       = izero
SEG(1:NSEGMAX)%DISEG    = .false.
SEG(1:NSEGMAX)%bloquer  = .false.
SEG(1:NSEGMAX)%unload   = .false.
SEG(1:NSEGMAX)%WAIT     = IZERO
SEG(1:NSEGMAX)%GRAIN    = IZERO
SEG(1:NSEGMAX)%SURFACE  = IZERO
SEG(1:NSEGMAX)%VARFREEPLAN = IZERO
SEG(1:NSEGMAX)%DEPINST  = ZERO
SEG(1:NSEGMAX)%probadev = ZERO
SEG(1:NSEGMAX)%taudev   = ZERO
SEG(1:NSEGMAX)%Nsysdev  = IZERO
SEG(1:NSEGMAX)%dom      = IZERO

! Initialization
dom_effective_deb(1) = zero
dom_effective_deb(2) = zero
dom_effective_deb(3) = zero
dom_effective_end(1) = real(Modur(1),DP)
dom_effective_end(2) = real(Modur(2),DP)
dom_effective_end(3) = real(Modur(3),DP)

rau_ini = -un
dist_moyenne_par = zero

if(Sideja < izero .or. sideja > itrois) then
  print *, " ERROR : Sideja = ",sideja
  print *, "STOP: This Option SIDEJA is not supported !"
  stop
endif

! A restart or not ?
select case (SIDEJA)

case (izero,ideux,itrois)

  write (*,*) " ======================================================================"
  write (*,*) " ================         New  simulations          ==================="
  call lire_segments
  if(seg_save) write (*,*) " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  if(seg_save) write (*,*) " ===      Loaded configuration resulted from previous simulation  ====="
  write (*,*) " ======================================================================"

#ifdef PA
  if ((Mon_Rang_PairImpair + IUN) >  IUN) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  else
#endif
    !====================================================================================
    !==============          Pour le fichier film     ===================================
    !====================================================================================
    write (93) avalue                  ! unite de simulation: relation espace reel-valeurs entieres (real)
    write (93) int(nbase,DPI_S)        ! nombre de vecteurs de discretisation (entier)
    ! les 3 composantes du vecteur ligne et plan (3x2 entier))
    write (93) (int(bveclin(:,I),DPI_S),int(bvecnor(:,I),DPI_S), I=1,nbase)
    write (93) int(modur(:),DPI_S)     ! dimensions de la boite (3 entier)
    !====================================================================================
    !====================================================================================
#ifdef PA
    ! The writing process ended their task, we can now liberate everybody
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
#endif

  if (cartograph == 5) then
    call double_config
    write(*,*) "Dedoublement de la configuration initial"
  endif

  if(GB /= izero) then
    RognSeg = .true.  ! first start, we need to cut segments that cross boundaries

    print *,  " ======================================================================"
    print *,  " ==========================   ATTENTION    ============================"

    select case (GB)

      case (1)
        print *," ===========  You are using domain option  (GB=1)  ============"
        call load_boundaries
        call domain_definition
        if (key_Rot_Dom) call DesorientMat

        if (nbgrains /= 0) then
          print *," ===========    with the grains disorientation option !    ============"
          print *," ===========    Nb of grains = ",nbgrains,"      ============"
          do i = 1,nsegm
            call CorrectifGrain(i,222)
          enddo
        endif

      case (2)
        print *," ========  You are using spherical boundaries option  (GB=2)  ========="

        call Selection_spherique

    end select

    print *,  " ======================================================================"

  else

    ! The default solution for 1 domain
    ! All segments are in domain 1
    seg(1:nsegm)%dom = iun

    ! All segments are interacting one with each other
    allocate(DomConnec(1,1))
    DomConnec(1,1) = .false.

  endif

  ! Random generator is initialized
  call init_seeds

  ! The input configuration is checked
  call verification_configuration

#ifdef MDC

#if defined(MDC) && defined(PA)
  mdcvol = ZERO
  if (Mon_Rang ==  IZERO) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
#endif

    call send_globdata

    call recv_globdata

    ! Voltera loop process is calculated
    call init_plas_eigen


#if defined(PA) && defined(MDC)
  endif

  CALL MPI_BCAST(mdcvol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

#endif

  call time_step

  !!!ifdef MDC close
#endif

  !==========Initiation de la methode des boites    =============================
  ! il y a trois modes de fonctionnement de la methode des boites:
  ! methode_boite = .false.  : methode non appliquee : NBOITES = zero
  ! allocation_dynamique_boites = .true. l'allocation est dynamque
  ! methode_boite = .true.  et allocation_dynamique_boites = .false. methode de boites fixe
  !

#ifdef MDC
  !in case of MDC calculations the simulated volume is defined from the FE mesh
  VOLUME=mdcvol
#else

#ifdef PA
  if (Mon_Rang ==  IZERO) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
#endif

    if (GB /= 1) then
      VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
    else
      !with GB = 1 it is possible to compute the volume of finite domain launching mm with --domvtx option
      open(411,file="../out/volumepy.txt",STATUS='OLD',iostat=IOS)

      if (IOS > 0) then
        write(*,*) 'GB = 1 and a dislocated subdomain is not defined !'
        write(*,*) 'Simulation will proceed by assuming dislocations in the total simulation volume  !'
        write(*,*) 'pause !'
        read *
        VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
      else
        read(411,*) VOLUME
      endif

      close(411)

      if (VOLUME < IZERO) then
        write(*,*) 'GB = 1 and a dislocated subdomain was not correctly defined !'
        write(*,*) 'Simulation will proceed by assuming dislocations in the total simulation volume  !'
!         write(*,*) 'pause !'
!         read *
        VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
      endif

      if (IOS == 0) then
        ! The tmp file ../out/volumepy.txt can now be eliminated
        write(cmd,'("rm ../out/volumepy.txt")')
        call system(cmd)
      endif

    endif

#ifdef PA
  endif

  CALL MPI_BCAST(volume,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

#endif

#endif

  ! The effective volume boundaries are defined for the Greengard Algo
  if (GB == 1) call effective_boundaries

  ! Variables needed for the PBC
  InvmoduR(1:3)   = un / Modur(1:3)
  HalfModur(1:3)  = Modur(1:3) / 2
  ModurMax        = maxval(Modur(1:3))

  ! transformation de la taille des boites en avalue
  L_Boite = L_boite * 1.0D-6 / avalue

  ! si L_boite est >>> , alors pas de methode de boite , c'est-a-dire NBoites =IUN
  methode_boites = (deux * L_Boite <= real(minval(modur(:))))

  ! par convention si L_boite est negatif, les boites sont allouees dynamiquement.
  allocation_dynamique_boites = (L_boite < ZERO )

  ! verification de la compatibilite entre le
  if (cartograph == 2 .and. nombre_fr /= 2) stop " impossible de faire la carto pour N source FR /= 2"

  write(*,*)
  write(*, '(" Dimensions of the volume in AVALUE (",3(I15,1x),")")') MODUR(1:3)
  write(*, '(" Coordinates of the volume center   (",3(I15,1x),")")') MODUR(1:3)/2
  write(*, '(" Dimensions of the volume in micron (",3(F15.5,1x),")")') MODUR(1:3)*avalue*1D6
  write(*,*)
  write(*, '(" The total simulation volume in micron^3 =",1(F15.5,1x),"")') &
      modur(1)*avalue*1.D6*modur(2)*avalue*1.D6*modur(3)*avalue*1.D6
  write(*, '(" The dislocated simulation volume in micron^3 =",1(F15.5,1x),"")') &
      volume*(avalue*1.D6)**3
  write(*,*)

  if (particules)  then
    call generer_particules
  write(*,'("Fraction volumique des Particules : ",F6.3," %")')V_particules / Volume * 100E0
    dist_moyenne_par = V_particules / npar / quatre * trois / PII
    dist_moyenne_par = (dist_moyenne_par)**(un/trois)
    dist_moyenne_par = UN/DSQRT(dist_moyenne_par * npar / volume)
    write(*,'("Average planar distance between particules : ",F6.1," nm")')  dist_moyenne_par * avalue *1E+9
  endif

  print *,'Number of segments before connection', nsegm

  call connecxion_segments

  !when the loop key is detected, the subroutine below generates loops and the corresponding
  ! segments are added to the NSEGM
  if(loops) call generation_loops

case (IUN) ! The sideja case = 1 means restart
  if (GB == 1) then
     call load_boundaries
      if (key_Rot_Dom) then
          call DesorientMat
      endif
  endif

  call loadall

  !#ifdef MDC
  !
  !#if defined(MDC) && defined(PA)
  ! if (Mon_Rang ==  IZERO) then
  !#endif
  !
  !   call send_globdata
  !
  !   call recv_globdata
  !
  !#if defined(MDC) && defined(PA)
  ! endif
  !
  !  CALL MPI_BCAST(mdcvol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !  ! Only the process zero of the communicator can write in output files, the others are waiting
  !#endif
  !
  !call time_step
  !if (sideja==1) deltat0=deltat
  !
  !#endif

  write(*,*)
  write(*,*) "======================================================================"
  write(*,*) "============      restart of a simulation      ======================="
  write(*,10) kk,epso*100.0,sigma*xmu*1D-6,nsegm
10 format("   Iteration : ", I6," Eps =",F7.4,"%; Sigma = ",F6.2," MPa; Nseg = ",I6)
  write(*,*) "======================================================================"
  write(*,*) "======================================================================"

  ! Relaxation steps are set to zero since it is a restart
  write(*,*) "Relaxation steps are set to zero since it is a restart ...."
  Relax_TL   = Izero
  Relax_INT  = Izero
  Relax_Reac = Izero
  Nstep = Nstep + KK

  ! Re activer le calcul des matrice de desorientation dans le cas du pavage
  !   if (desorientgrain.and.GB==ITROIS)  call DesorientMat

#ifdef PA
  if ((Mon_Rang_PairImpair + IUN) >  IUN) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  else
#endif

    write(*,*)
    write(*,*) '   !!!   Please wait the trajectory file  film.bin  is read again   !!!'
    write(*,*)

    ! The header of the film.bin file is passed
    read(93,iostat=koderr)            ! read cristallo
    read(93,iostat=koderr)            ! read nsegmaxf
    read(93,iostat=koderr)            ! read avalue
    read(93,iostat=koderr)            ! read nbase
    read(93,iostat=koderr)            ! read (bveclin(:,I),bvecnor(:,I), I=1,nbase)
    read(93,iostat=koderr)            ! read modur(:)

    ! The file film.bin is stop at the restart step KK (the latter is not saved at restart)
    kompt     = -1
    komptlast = -1
    do while (kompt <= kk)
      caractest = '!'
      read(93,iostat=koderr) kompt
      read(93,iostat=koderr)
      read(93,iostat=koderr) caractest

      if (caractest /= "&") then
        ! The previous stop of the simulation corrupted the film.bin file
        ! We must clean up the file by removing the last corrupted step
        write(*,*)
        write(*,*) '   The film.bin file was corrupted, the last saved configuration is at step  =', komptlast
        write(*,*)
        rewind(93)
        ! The header of the film.bin file is passed
        read(93,iostat=koderr)            ! read cristallo
        read(93,iostat=koderr)            ! read nsegmaxf
        read(93,iostat=koderr)            ! read avalue
        read(93,iostat=koderr)            ! read nbase
        read(93,iostat=koderr)            ! read (bveclin(:,I),bvecnor(:,I), I=1,nbase)
        read(93,iostat=koderr)            ! read modur(:)
        kompt2  = -1
        do while (kompt2 /= komptlast)
          read(93,iostat=koderr) kompt2
          read(93,iostat=koderr)
          read(93,iostat=koderr)
        enddo
        exit !We went at the last correctly written step, exit do-while loop
      else
        komptlast = kompt    ! The last step correctly written
      endif

    enddo

#ifdef PA
    ! The writing process ended their task, we can now liberate everybody
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  endif
#endif


  if(GB /= izero) then
    RognSeg = .false.   ! no need to cut segments for restart simulation

    print *,  " ======================================================================"
    print *,  " ==========================   ATTENTION    ============================"

    select case (GB)

    case (1)
      print *," ===========  You are using domain option  (GB=1)  ============"
      call domain_definition

      if (nbgrains /= 0) then
        if (key_Rot_Dom) print *," ===========    with the grains disorientation option !    ============"
        print *," ===========    Nb of grains = ",nbgrains,"      ============"
      endif

    case (2)
      print *," ========  You are using spherical boundaries option  (GB=2)  ========="
      call Selection_spherique

    end select

    print *,  " ======================================================================"

  else

    ! If 1 domain, every segments are interacting
    allocate(DomConnec(1,1))
    DomConnec(1,1) = .false.

  endif

!#ifdef MDC

!#if defined(MDC) && defined(PA)
! if (Mon_Rang ==  IZERO) then
!#endif
!
!   call send_globdata
!
!   call recv_globdata
!
!#if defined(MDC) && defined(PA)
! endif
!
!  CALL MPI_BCAST(mdcvol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  ! Only the process zero of the communicator can write in output files, the others are waiting
!#endif
!
!
!#endif

  InvmoduR(1:3)=UN/Modur(1:3)
  !==========Initiation de la methode des boites    =============================
  ! il y a trois modes de fonctionnement de la methode des boites:
  ! methode_boite = .false.  : methode non appliquee : NBOITES = zero
  ! allocation_dynamique_boites = .true. l'allocation est dynamque
  ! methode_boite = .true.  et allocation_dynamique_boites = .false. methode de boites fixe
  !

  ! transformation de la taille des boites en avalue
  L_Boite = L_boite * 1.0D-6 / avalue

  ! si L_boite est >>> , alors pas de methode de boite , c'est-a-dire NBoites =IUN
  methode_boites = (deux * L_Boite <= real(minval(modur(:))))

  ! par convention si L_boite est negatif, les boites sont allouees dynamiquement.
  allocation_dynamique_boites = (L_boite < ZERO )
  !==========================================================================================

#ifdef MDC
  !in case of MDC calculations the simulated volume is defined from the FE mesh
  VOLUME=mdcvol
#else

#ifdef PA
  if (Mon_Rang ==  IZERO) then
    ! Only the process zero of the communicator can write in output files, the others are waiting
#endif

  if (GB /= 1) then
    VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
  else
    !with GB = 1 it is possible to compute the volume of finite domain launching mm with --domvtx option
    open(411,file="../out/volumepy.txt",STATUS='OLD',iostat=IOS)

    if (IOS > 0) then
      write(*,*) 'GB = 1 and a dislocated subdomain is not defined !'
      write(*,*) 'Simulation will proceed by assuming dislocations in the total simulation volume  !'
      write(*,*) 'pause !'
      read *
      VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
    else
      read(411,*) VOLUME
    endif

    close(411)

    if (VOLUME < IZERO) then
      write(*,*) 'GB = 1 and a dislocated subdomain was not correctly defined !'
      write(*,*) 'Simulation will proceed by assuming dislocations in the total simulation volume  !'
      write(*,*) 'pause !'
      read *
       VOLUME = DFLOAT(ModuR(1))*DFLOAT(Modur(2))*DFLOAT(ModuR(3))
    endif

    if (IOS == 0) then
      ! The tmp file ../out/volumepy.txt can now be eliminated
      write(cmd,'("rm ../out/volumepy.txt")')
      call system(cmd)
    endif

  endif

#ifdef PA
  endif

  CALL MPI_BCAST(volume,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

#endif

#endif

  ! The effective volume boundaries are defined for the Greengard Algo
  if (GB == 1) call effective_boundaries

  Nombre_FR = IZERO
  do i = 1 , nsegm
    if( seg(i)%voiso == nsegmax ) nombre_FR = Nombre_FR + 1
  enddo
  print*, " "
  print*, "Initial configuration specificity"
  write(*, '("  Number of Frank_Read sources: ",I10)') Nombre_FR

  print *, " "
  write(*, '(" Dimensions of the volume in AVALUE (",3(I10,1x),")")') MODUR(1:3)
  write(*, '(" Coordinates of the volume center   (",3(I10,1x),")")') MODUR(1:3)/2
  write(*, '(" Dimensions of the volume in micron (",3(F10.3,1x),")")') MODUR(1:3)*avalue*1D6

end select

if (methode_boites) then
  ! The number of replica to apply in the Greengard method if not defined by user
  if (Min(PBCIxDIM,PBCIyDIM,PBCIzDIM) < izero) then
    MAXMODUR = Max(Modur(1), Modur(2), Modur(3))
    PBCIxDIM = INT((MAXMODUR-1)/(Modur(1)*Ideux)) ! The number of times the simulated volume must
    PBCIyDIM = INT((MAXMODUR-1)/(Modur(2)*Ideux)) ! be repeated in the 3D space to calculate the LR
    PBCIzDIM = INT((MAXMODUR-1)/(Modur(3)*Ideux)) ! stress in a volume more or less isotropic
  endif
  print *, " "
  write(*, '(" Number of replicas applied in the Greengard method (",3(I3,1x),")")') PBCIxDIM,PBCIyDIM,PBCIzDIM
  print *, " "
endif

end subroutine generation_microstructure

!#####################################################################
! The effective volume boundaries are defined for the Greengard Algo #
!#####################################################################
subroutine effective_boundaries

implicit none

integer             :: IOS                          ! Read iostat error message (0 = OK, > 0 means error)
integer(kind=dpi)   :: i
integer             :: Nbvertex
real                :: vertex(3)

! Initialization
dom_effective_end(1) = 0
dom_effective_end(2) = 0
dom_effective_end(3) = 0
dom_effective_deb(1) = Modur(1)
dom_effective_deb(2) = Modur(2)
dom_effective_deb(3) = Modur(3)

open(412,file="../out/vertexlist",STATUS='OLD',iostat=IOS)

if (IOS > 0) then

  write(*,*) 'GB = 1 and the file ../out/vertexlist is not defined!'
  write(*,*) 'Simulation will proceed by considering in the Multipole Algo the total periodic simulation volume !'
  write(*,*) 'pause !'
  read *

  dom_effective_deb(1) = zero
  dom_effective_deb(2) = zero
  dom_effective_deb(3) = zero
  dom_effective_end(1) = real(Modur(1),DP)
  dom_effective_end(2) = real(Modur(2),DP)
  dom_effective_end(3) = real(Modur(3),DP)

else

  koderr = 0
  read(412,*,iostat=koderr)   ! This line contains the number of vertex domains, but is useless here

  do while(koderr == 0)

    ! Vertex are organized in paquet of Nbvertex elements
    read(412,*,iostat=koderr) Nbvertex

    if (koderr == 0) then

      ! The nbvertex are read and tested to see if they are boundaries of the effective simulation volume
      do i = 1, Nbvertex

        read(412,*,iostat=koderr) vertex(1:3)

        if (vertex(1) < real(dom_effective_deb(1))) dom_effective_deb(1) = int(vertex(1))
        if (vertex(2) < real(dom_effective_deb(2))) dom_effective_deb(2) = int(vertex(2))
        if (vertex(3) < real(dom_effective_deb(3))) dom_effective_deb(3) = int(vertex(3))
        if (vertex(1) > real(dom_effective_end(1))) dom_effective_end(1) = int(vertex(1))
        if (vertex(2) > real(dom_effective_end(2))) dom_effective_end(2) = int(vertex(2))
        if (vertex(3) > real(dom_effective_end(3))) dom_effective_end(3) = int(vertex(3))

      enddo

    endif

  enddo

endif

close(412)

! results of the identification
write(*,*)
write(*,*) 'The boundaries of the simulation effective volume are :'
write(*,*) 'X direction --> deb, end   ', dom_effective_deb(1), dom_effective_end(1)
write(*,*) 'Y direction --> deb, end   ', dom_effective_deb(2), dom_effective_end(2)
write(*,*) 'Z direction --> deb, end   ', dom_effective_deb(3), dom_effective_end(3)

end subroutine effective_boundaries

!############################################################
!> Managing segment touching a free surfaces
!! seg(i)%surface=1 the origin of i is touching the surface
!! seg(i)%surface=2 the end of i is touching the surface
!! At the end, the segment O or E must be on a surface or
!! the first position outside the interface.
!! If SigDep is set to NSEGMAX, no correction of swept area is made
!############################################################
subroutine resolve_freeSurface_hooking_up(i,i1o,i1e,lonad,lonvoad,lonvead,SIGdep,absdep,debugref)

implicit none

integer(kind=DPI)  :: LONAD   !< Segment length after displacement
integer(kind=DPI)  :: LONVOAD !< Origin Neighbor Segment length after displacement
integer(kind=DPI)  :: LONVEAD !< End Neighbor Segment length after displacement
integer(kind=DPI)  :: SIGdep  !< Sign of displacement (-1,1, or set to NSEGMAX)
integer(kind=DPI)  :: Li      !< Segment line direction
integer(kind=DPI)  :: itemp   !< Temporary variable
integer(kind=DPI)  :: I       !< Segment indice
integer(kind=DPI)  :: J       !< Free Surface indice
integer(kind=DPI)  :: absdep  !< Absolute displacement of segment I
integer(kind=DPI)  :: ii      !< Counter for while loops
integer(kind=DPI)  :: i1o     !< First Neighbor at Origin
integer(kind=DPI)  :: i1e     !< First Neighbor at End
integer(kind=DPI)  :: jj      !< Loop counter over free surfaces
integer(kind=DPI)  :: gg, irot    !< counter

integer  :: debugref  !< Debugging Reference

integer(kind=DPI),dimension(3)  :: OldOi    !< Store segment Origin when entering subroutine
integer(kind=DPI),dimension(3)  :: OldEi    !< Store segment End when entering subroutine
integer(kind=DPI),dimension(3)  :: Oi       !< Segment Origin
integer(kind=DPI),dimension(3)  :: Ei       !< Segment End
integer(kind=DPI),dimension(3)  :: tempOi   !< Temporary variable storing segment origin
integer(kind=DPI),dimension(3)  :: tempEi   !< Temporary variable storing segment end
integer(kind=DPI),dimension(3)  :: Trans    !< A translation vector use in area swept correction
integer(kind=DPI),dimension(3)  :: VecI     !< Translation vector of Oi
integer(kind=DPI),dimension(3)  :: OldVecI  !< Translation vector of OldOi
integer(kind=DPI),dimension(3)  :: inormint !< Surface burger indices (Integer)

real(kind=DP),dimension(3)      :: Oireal   !< Segment Origin
real(kind=DP),dimension(3)      :: Eireal   !< Segment End
real(kind=DP),dimension(3)      :: Intersec !< Point of intersection between segment and the surface
real(kind=DP),dimension(3)      :: VecO     !< Vector segment origin / surface intersection point
real(kind=DP),dimension(3)      :: VecE     !< Vector segment end / surface intersection point
real(kind=DP),dimension(3)    :: IntersecO  !< Point of intersection between segment and the surface when entering with Origin
real(kind=DP),dimension(3)    :: IntersecE  !< Point of intersection between segment and the surface when entering with End
real(kind=DP),dimension(3)    :: Vect1O     !< Vector segment origin / IntersecO
real(kind=DP),dimension(3)    :: Vect1E     !< Vector segment end  / IntersecE
real(kind=DP),dimension(3)    :: normint    !< Surface burger indices (Real)

real(kind=DP)      :: TriArea !< Triangle area swept
real(kind=DP)      :: a,b,c,s !< The Heron 's' formula variables
real(kind=DP)      :: temp    !< Temporary variable for force calculation

logical                       :: KeepEndCase  !< True if we have to modify segment length

if(kkdebug) then
  write(379,fmt='(" ")')
  write(379,fmt='("*** resolve_freeSurface_hooking_up, ref: ",I7," ****")') debugref
  write(379,fmt='("--> In resolve_freeSurface_hooking_up we modify i = ",I7,":")') i
  write(379,fmt='(" ")')
  call seginfo(i,'Begining resolve_freeSurface_hooking_up subroutine')
endif

! Test on segment i which as nothing to do here
if (seg(i)%VarFreePlan == izero) then
  seg(i)%surface = izero
  return
endif

! Information on segment i before displacement
Li    = seg(i)%veclin
OldOi = seg(i)%O
OldEi = OldOi(1:3)+seg(i)%norme*BVECLIN(1:3,Li)

! Information on segment i after displacement
!Oi(1:3)       = SEG(i1o)%O(1:3) + BVECLIN(1:3,SEG(i1o)%VECLIN)*lonvoad
Oi(1:3)       = oldOi(:) + BVECLIN(1:3,SEG(i1o)%VECLIN)*(lonvoad-seg(i1o)%norme)
SEG(I)%O(1:3) = Oi(:)
! SEG(I)%NORME  = lonad
Ei(1:3)       = Oi(1:3) + lonad * BVECLIN(1:3,Li)
OiReal = real(Oi(:),DP)
Eireal = real(Ei(:),DP)

! Point of intersection between segment and the surface
J=seg(i)%VarFreePlan
Intersec(1:3) = InterPlanSeg(Plane_MillerR(1:3,j),Plane_pos(j),Oireal,BVECLIN(:,Li))
! Distance to the interface
VecO(1:3) = Oireal(:)-Intersec(:)
VecE(1:3) = Eireal(:)-Intersec(:)

! Initialisations
tempOi(:) = -9999
tempEi(:) = -9999
keepEndCase = .true.
ii = -IUN

! Length correction for the segment touching the surface

if (kkdebug) then
  write (379,*) "Lengths information"
  write (379,fmt='(4A8)') "lonad","lonvoad", "lonvead", "SIGdep"
  write (379,fmt='(4I8)') lonad,lonvoad,lonvead,SIGdep
  write (379,*) ""
  write (379,fmt='(5A8)') "seg","vo","ve", "surf", "VFP"
  write (379,fmt='(5I8)') i, i1o, i1e, seg(i)%surface, seg(i)%VarFreePlan
  write (379,*) ""
  write (379,fmt='("Information on segment ",I7,": BEFORE displacement")') i
  write (379,fmt='(9A8)') "seg","veclin",'',"O(1:3)",'',"norme",'',"E(1:3)",''
  write (379,fmt='(9I8)') i, Li, OldOi, seg(i)%norme, OldEi
  write (379,*) ""
  write (379,fmt='("Information on segment ",I7,": AFTER displacement")') i
  write (379,fmt='(9A8)') "seg","veclin",'',"O(1:3)",'',"norme",'',"E(1:3)",''
  write (379,fmt='(9I8)') i, seg(i)%veclin, seg(i)%o, lonad, Ei(1:3)
endif

select case (seg(i)%surface)

  !*********
  case (IUN) !the origin of i is touching the surface

    a=lonvoad
    lonvoad=IZERO !the kneecap at the origin moves with the segment i

    if (dot_product(VecO(1:3),Plane_MillerR(1:3,j)) >= ZERO) then
      ! Oi is outside the volume and we look for the first coordinate outside
      ! The segment length is decreased
      do while (dot_product(VecO(1:3),Plane_MillerR(1:3,j)) >= ZERO)
        ii = ii + IUN
        if (ii > 1000) then
         write(*,fmt='("pb of infinit do while loop 1 in resolve_freeSurface_hooking_up", F20.3, "kk=",I7)') &
          dot_product(VecO(1:3),Plane_MillerR(1:3,j)),kk
          call seginfo(i,"Connection to surface attempt FS_hooking_up 1")
          write(*,*) "Oireal(:)-Intersec(:)",Oireal(:)-Intersec(:), "dotprod", &
          & dot_product(Oireal(:)-Intersec(:),Plane_MillerR(1:3,j))
          write(*,*) "Failed to connect surface ",J,"Plane_MillerR(1:3,j)",Plane_MillerR(1:3,j)," debugref: ", debugref
          call flush(379)
          stop
        endif
        tempOi(1:3) = Oi(:)
        Oi(1:3) = Oi(:)+BVECLIN(1:3,Li)
        VecO(1:3) = real(Oi(:)-Intersec(:))
      enddo

      if(kkdebug) write(379,*) 'case 1: the segment length after displacement is decreased by', ii
      lonad = lonad - ii
      Oi(:) = TempOi(:)
      seg(i)%norme = lonad
      seg(i)%O(:) = Oi(:)
      seg(i1o)%O(:) = Oi(:)
      Ei(:) = Oi(:) + lonad * BVECLIN(:,Li)
      SEG(i1e)%O(:) = Ei(:)
      seg(i1e)%norme = lonvead

        if (seg(i)%voise /=i1e) then
            gg=seg(i)%voise
            do while (gg /= i1e)
              irot=seg(gg)%voise
              seg(gg)%O =Ei

              gg=irot
            enddo
          endif

    else
      ! Oi is outside and we look for the closest coordinate to the surface
      ! the segment length is increased
      do while (dot_product(VecO(1:3),Plane_MillerR(1:3,j)) < ZERO)
        ii = ii + IUN
        if (ii > 1000) then
          write(*,fmt='("pb of infinit do while loop 2 in resolve_freeSurface_hooking_up", F20.3, "kk=",I7)') &
          dot_product(VecO(1:3),Plane_MillerR(1:3,j)),kk
          call seginfo(i,"Connection to surface attempt FS_hooking_up 2")
          write(*,*) "Oireal(:)-Intersec(:)",Oireal(:)-Intersec(:), "dotprod", &
          & dot_product(Oireal(:)-Intersec(:),Plane_MillerR(1:3,j))
          write(*,*) "Failed to connect surface ",J,"Plane_MillerR(1:3,j)",Plane_MillerR(1:3,j)," debugref: ", debugref
          call flush(379)
          stop
        endif
        tempOi(1:3) = Oi(:)
        Oi(1:3) = Oi(:) - BVECLIN(:,Li)
        VecO(1:3) = real(Oi(:)-Intersec(:),DP)
      enddo

      ! We test the possibility that the new coordinate is outside a fix boundary
      Oireal=Oi

      do jj= 1,nbplanMax-NbFreePlan

        if (Plane_dom(jj) /=  Plane_dom(j) .or. Plane_knd(jj) == IUN) cycle

        inormint(:)= Plane_MillerI(1:3,jj)
        normint(:) = Plane_MillerR(1:3,jj)
        IntersecO(:) = InterPlanSeg(normint(:),Plane_pos(jj),Oireal,inormint(:))

        ! Normale projection vector
        Vect1O(:) = Oireal- IntersecO(:)

        ! The segment displacement is reduced to zero since we touch a fix boundary
        if (dot_product(Vect1O,normint) >= zero) then
          lonad   = seg(i)%norme
          lonvoad = seg(i1o)%norme
          lonvead = seg(i1e)%norme
          SEG(I)%O(1:3) = OldOi
          TriArea = zero
          daire(i) = zero
          KeepEndCase = .false.
          absdep = IZERO
          if (kkdebug) write(379,*) "segment at a free surface is stopped as it touch a fix boundary",jj," at O"
        endif
      enddo

      if (KeepEndCase) then
        ii = ii + iun
        if(kkdebug) write(379,*) 'case 1bis: the segment length after displacement is increased by', ii
        lonad = lonad + ii
        seg(i)%norme = lonad
        seg(i)%O(1:3) = Oi(1:3)
        seg(i1o)%O(1:3) = Oi(1:3)
        Ei(:) = Oi(:) + (lonad) * BVECLIN(1:3,Li)
        SEG(i1e)%O(1:3) = Ei(1:3)
        seg(i1e)%norme = lonvead

          if (seg(i)%voise /=i1e) then
            gg=seg(i)%voise
            do while (gg /= i1e)
              irot=seg(gg)%voise
              seg(gg)%O =Ei

              gg=irot
            enddo
          endif
      endif

    endif

    if (KeepEndCase .and. Sigdep /= nsegmax) then
      if (ii /= izero) then
        ! we first calculate the length of the three sides of the triangle we must substract to daire
        a =  a * NORMLIN(seg(i1o)%veclin)
        ! We translate coordinates to the center of the simulated volume to
        ! avoid modulo problems when computing (Veci - OldVeci)
        Trans(:) = int(Modur(:) * half - OldOi(:),DPI)
        OldVecI(:)=modulo(OldOi(:)+Trans(:),modur(:))
        VecI(:)=modulo(Oi(:)+Trans(:),modur(:))
        b = sqrt(real((Veci(1)-OldVeci(1))*(Veci(1)-OldVeci(1))+  &
                      (Veci(2)-OldVeci(2))*(Veci(2)-OldVeci(2))+  &
                      (Veci(3)-OldVeci(3))*(Veci(3)-OldVeci(3))))
        c = ii * NORMLIN(Li)
        ! The Heron s formula
        s  = 0.5 * (a + b + c)

        if ( abs(s * (s-a) * (s-b) * (s-c)) < numtold ) then !to avoid negative values close to zero
          TriArea=zero
        else
          TriArea = sqrt(s * (s-a) * (s-b) * (s-c)) * SIGdep * sign(1.,real(ii))
        endif

      else
        TriArea = zero
      endif

    endif

  !***********
  case (IDEUX) !the end of i is touching the surface

    a=lonvead
    lonvead=IZERO !the kneecap at the end moves with the segment i
    if (kkdebug) write(379,*) j, Li, BVECLIN(:,Li)
    if (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) >= ZERO) then
      ! Ei is outside and we look for the closest coordinate to the surface
      ! the segment length is decreased
      do while (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) >= ZERO)
        ii = ii + IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 3 in resolve_freeSurface_hooking_up",  &
                dot_product(VecE(1:3),Plane_MillerR(1:3,j)), "kk=",kk
          stop
        endif
        tempEi(1:3) = Ei(:)
        Ei(1:3) = Ei(:)-BVECLIN(:,Li)
        if (kkdebug) write(379,*) Ei(:)
        VecE(1:3) = real(Ei(:)-Intersec(:),DP)
      enddo

      if(kkdebug) write(379,*) 'case 2: the segment length after displacement is decreased by', ii
      lonad = lonad - ii
      Ei(:) = tempEi(:)
      seg(i)%norme = lonad
      seg(i1e)%O(:) = Ei(:)
      seg(i1o)%norme = lonvoad

      if (seg(i)%voiso /=i1o) then
            gg=seg(i)%voiso
            do while (gg /= i1o)
              irot=seg(gg)%voiso
              seg(gg)%O = seg(i)%O(:)

              gg=irot
            enddo
          endif

    else
      ! Ei is in the volume and we look for the first coordinate outside
      ! The segment length is increased
      do while (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) < ZERO)
        ii = ii + IUN
        if (ii > 10000) then
          write(379,*) "pb of infinit do while loop 4 in resolve_freeSurface_hooking_up",  &
                dot_product(VecE(1:3),Plane_MillerR(1:3,j)) ,"kk=",kk
          write(379,*) "j",i,j,Plane_MillerR(1:3,j)
          write(379,*) "Intersec", Intersec
          write(379,*) "Ei initial", Ei(:)
          write(379,*) "BVECLIN(:,Li)", BVECLIN(:,Li)
          stop "See ../out/debug files !!!"
        endif
        tempEi(1:3) = Ei(:)
        Ei(1:3) = Ei(:) + BVECLIN(:,Li)
        VecE(1:3) = real(Ei(:)-Intersec(:))
      enddo

      ! We test the possibility that the new coordinate is outside a fix boundary
      Eireal=Ei

      do jj= 1,nbplanMax-NbFreePlan

        if (Plane_dom(jj) /=  Plane_dom(j) .or. Plane_knd(jj) == IUN) cycle

        inormint(:)= Plane_MillerI(1:3,jj)
        normint(:) = Plane_MillerR(1:3,jj)
        IntersecE(:) = InterPlanSeg(normint(:),Plane_pos(jj),Eireal,inormint(:))

        ! Normale projection vector
        Vect1E(:) = Eireal(:)-IntersecE(:)

        ! The segment displacement is reduced to zero since we touch a fix boundary
        if (dot_product(Vect1E,normint) >= zero) then
          lonad   = seg(i)%norme
          lonvoad = seg(i1o)%norme
          lonvead = seg(i1e)%norme
          SEG(I)%O(1:3) = OldOi
          TriArea = zero
          daire(i) = zero
          KeepEndCase = .false.
          absdep = IZERO
          if (kkdebug) write(379,*) "segment at a free surface is stopped as it touch a fix boundary",jj," at E"
        endif
      enddo

      if (keepEndCase) then
        ii = ii + iun
        if(kkdebug) write(379,*) 'case 2bis: the segment length after displacement is increased by', ii
        lonad = lonad + ii
        seg(i)%norme = lonad
        seg(i1e)%O(1:3) = Ei(1:3)
        seg(i1o)%norme = lonvoad

          if (seg(i)%voiso /=i1o) then
            gg=seg(i)%voiso
            do while (gg /= i1o)
              irot=seg(gg)%voiso
              seg(gg)%O = seg(i)%O(:)

              gg=irot
            enddo
          endif

      endif

    endif

    if (kkdebug) write(379,fmt='(A15,L7)') "KeepEndCase",KeepEndCase

    if (KeepEndCase .and. Sigdep /= nsegmax) then
      if (ii /= izero) then
        ! we first calculate the length of the three sides of the triangle we must substract to daire
        a =  a * NORMLIN(seg(i1e)%veclin)
        ! We translate coordinates to the center of the simulated volume to
        ! avoid modulo problems when computing (Veci - OldVeci)
        Trans(:) = int(Modur(:) * half - OldEi(:),DPI)
        OldVecI(:)=modulo(OldEi(:)+Trans(:),modur(:))
        VecI(:)=modulo(Ei(:)+Trans(:),modur(:))
        b = sqrt(real((Veci(1)-OldVeci(1))*(Veci(1)-OldVeci(1))+  &
                      (Veci(2)-OldVeci(2))*(Veci(2)-OldVeci(2))+  &
                      (Veci(3)-OldVeci(3))*(Veci(3)-OldVeci(3))))
        c = ii * NORMLIN(Li)
        ! The Heron s formula
        s  = 0.5 * (a + b + c)
        if ( abs(s * (s-a) * (s-b) * (s-c)) < numtold ) then !to avoid negative values close to zero
          TriArea=zero
        else
          TriArea = sqrt(s * (s-a) * (s-b) * (s-c)) * SIGdep * sign(1.,real(ii))
        endif
      else
        TriArea = zero
      endif
    endif

  !***********
  case DEFAULT
    print *, 'segments with surface type >2 must not appear in this part (resolve_freeSurface_hooking_up)'
    stop

end select

if (kkdebug) then
  write(379,*) 'Oldoi=',oldOi,'OldEi=',oldEi
  write(379,*) 'Oi=',Oi,'Ei',Ei
  write(379,*) 'seg(i)%o', seg(i)%o
  write(379,fmt='(A11,I7,A11,I7)') "SIGdep = ", SIGdep, "NSEGMAX = ", NSEGMAX
endif

if (SIGdep /= NSEGMAX) then

  if (kkdebug)  write(379,*) "AREA SWEPT CALCULATION IN FSSURFACE_HOOKINGUP"
  ! Correction of the area swept by segment touching the surface in unit a2
  daire(i) = daire(i) + TriArea

  ! This quantity is attribute to the slip systems activity
  itemp = syseg(Li)
  AireSys(itemp) = AireSys(itemp) + daire(i)
  AireSysInst(itemp) = AireSysInst(itemp) + daire(i)

  ! Calculation for the polycrystal case
  if (desorientgrain) then
    AireSysGRAIN(seg(I)%grain,itemp) = AireSysGRAIN(seg(i)%grain,itemp) + daire(i)
  endif

  ! mobile dislocation density calculation
  RAUDMO_I = RAUDMO_I+lonad*NORMLIN(Li)
  densitemob=RAUDMO_I        ! utile pour le control en metallofute2

  ! Calcul of the mechanical work per step associated to the applied stress = tauApp*gamma
  ! Attention: in the simulation, force is calculated as tau*b*longueur(en a)/(mu*a)
  temp = xmu/VOLUME*daire(i)
  TrAppSysInst(itemp) = TrAppSysInst(itemp)+tauApP(I)*temp

  ! Calcul of the mechanical work per step associate to the internal stress
  TrIntSysInst(itemp)=TrIntSysInst(itemp)+(tauInt(I)+tauTL(I))*temp

  ! A paranoid test to check that anglevis is always defined
  if (seg(i)%anglevis < zero) then
    print*,'A problem is found in the calculation of anglevis: microstructure', kk, i
    stop
  endif

  if (DesorientGrain) stop 'DesorientGrain grain with free surface is not yet considered'

  ! calculation of the area swep by screw and non-screw dislocations
  if (abs(seg(i)%anglevis) < Critic_angle) then
    AireVisSys(itemp) = AireVisSys(itemp) + daire(i)*sign(1.D0,SchmidSys(itemp))
  elseif (abs(halfpii-seg(i)%anglevis) < Critic_angle) then
    AireCoinSys(itemp) = AireCoinSys(itemp) + daire(i)*sign(1.D0,SchmidSys(itemp))
  else
    AireMixtSys(itemp) = AireMixtSys(itemp) + daire(i)*sign(1.D0,SchmidSys(itemp))
  endif

endif


if (kkdebug) then
  write(379,*) "Modifications has been made in resolve_freeSurface_hooking_up"
  write(379,fmt='(7A8)') "i","i1o","i1e","lonad","lonvoad","lonvead","SIGdep"
  write(379,fmt='(7I8)') i,i1o,i1e,lonad,lonvoad,lonvead,SIGdep
  call seginfo(i,'Ending resolve_freeSurface_hooking_up')
endif

end subroutine resolve_freeSurface_hooking_up

!########################################################################
!> Managing segment going out of a free surfaces (close loops)          #
!> seg(i)%surface=4 the segment i is going out with the origin          #
!> seg(i)%surface=5 as in  case 4 but now i is going out with the end   #
!> seg(i)%surface=6 the segment i is going out with both extremity,     #
!>        i will be eliminated, i1o e i1e cut to lie at the surface.    #
!> seg(i)%surface=7 segment touching the surface with the origin        #
!>                              crossing surface with the end.          #
!> seg(i)%surface=8 segment touching the surface with the end           #
!>                              crossing surface with the origin        #
!########################################################################
subroutine surface_correct(i,i1o,i1e)

integer(kind=DPI)  :: I     !< Segment I index
integer(kind=DPI)  :: Ivnne !< Segment I Vnne index
integer(kind=DPI)  :: Ivnno !< Segment I Vnne index
integer(kind=DPI)  :: J     !< Segment I VarFreePlan
integer(kind=DPI)  :: ii    !<
integer(kind=DPI)  :: jj    !<
integer(kind=DPI)  :: i1o   !<
integer(kind=DPI)  :: i1e   !<
integer(kind=DPI)  :: iNew, iNew2, iNew3  !< Newly created segment
integer(kind=DPI)  :: irot  !<
integer(kind=DPI)  :: idestru  !<
integer(kind=DPI)  :: itemp  !<
integer(kind=DPI)  :: i1otemp  !<
integer(kind=DPI)  :: i1etemp  !<
integer(kind=DPI)  :: VLI   !< Segment line direction
integer(kind=DPI)  :: DomE  !<
integer(kind=DPI)  :: ve !<
integer(kind=DPI)  :: SPR
integer(kind=DPI)  :: binome !< Junction binome
integer(kind=DPI)  :: SURFTYPE    !< seg%surface
integer(kind=DPI)  :: counter     !< do-while loop counter

integer(kind=DPI),DIMENSION(3)  :: Oi     !<
integer(kind=DPI),DIMENSION(3)  :: Ei     !<
integer(kind=DPI),DIMENSION(3)  :: Oi1o   !<
integer(kind=DPI),DIMENSION(3)  :: Ei1o   !<
integer(kind=DPI),DIMENSION(3)  :: Oi1e   !<
integer(kind=DPI),DIMENSION(3)  :: tempOi !<
integer(kind=DPI),DIMENSION(3)  :: tempEi !<

real(kind=DP),DIMENSION(3)      :: Intersec   !<
real(kind=DP),DIMENSION(3)      :: VecO       !<
real(kind=DP),DIMENSION(3)      :: VecE       !<

real(kind=DP)  :: norm

!Case INEUF
real(kind=DP)      :: DistanceS   !<  FS Distance
real(kind=DP)      :: DistanceP   !<  Point Distance
real(kind=DP)      :: DistDiff    !< Store distanceS-distanceP value
real(kind=DP)      :: tmpDistDiff !< Store temporary distanceS-distanceP value

integer(kind=DPI)  :: Normi  !< Storage of initial I norm

real(kind=DP),dimension(3)      :: Oireal   !< Segment Origin
real(kind=DP),dimension(3)      :: Eireal   !< Segment End

logical      :: Condition       !<
logical      :: outside         !< True = outside domains
j=seg(i)%VarFreePlan
SPR=j

if(kkdebug) then
  write(379,*) ' '
  write(379,*) 'In subroutine surface_correct we modify i=',i
  write(379,*) 'seg',i,'surf reaction type',seg(i)%surface,'Fsurface', j ,'kk=', kk
  write(379,*) 'segment i1o =',i1o,'Oi1o =',seg(i1o)%O,'Ei1o =',seg(i1o)%O+seg(i1o)%norme*BVECLIN(:,seg(i1o)%veclin)
  write(379,*) 'segment i1e =',i1e,'Oi1e =',seg(i1e)%O,'Ei1o =',seg(i1e)%O+seg(i1e)%norme*BVECLIN(:,seg(i1e)%veclin)
  call seginfo(i,'debut surf correct')
endif

! Close segment displacement are stop to avoid complex configuration
IDEP(i1o)=IZERO
IDEP(i1e)=IZERO
IDEP(seg(i)%vnno)=IZERO
IDEP(seg(i)%vnne)=IZERO

! The following tmp variables are intialized with stupid values for the benefit of debuging
tempOi(:) = -9999
tempEi(:) = -9999
itemp=-9999
i1otemp=-9999
i1etemp=-9999

SURFTYPE = seg(i)%surface
!##########################################
! The different cases of emerging segments #
!##########################################
select case (SURFTYPE)

  !*************
  case (IQUATRE,IDOUZE) !the segment i is going out with the origin

    if (seg(I)%surface == IQUATRE) then

      itemp = i
      i1otemp = i1o
      i1etemp = i1e

    elseif (seg(I)%surface == IDOUZE) then

      itemp = i1e
      i1otemp = i
      i1etemp = seg(i1e)%voise

      seg(i1otemp)%surface      = IZERO
      seg(i1otemp)%VarFreePlan  = IZERO
      seg(i1otemp)%voiso        = NSEGMAX
      seg(i1otemp)%vnno         = NSEGMAX

      !Cleaning i1o side
      ii=I1o
      counter = IZERO
      do while (ii /= NSEGMAX)

          if (counter > 3) then
            call seginfo(ii,"Problem do-while loop surf correct cleaning ii SURTYPE IDOUZE")
            call seginfo(i1o,"Problem do-while loop surf correct cleaning i1o SURTYPE IDOUZE")
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voiso ",seg(ii)%voiso,":"
            stop "Problem do-while loop surf correct cleaning i1o SURTYPE IDOUZE"
          endif
          counter = counter + IUN

          idestru=seg(ii)%voiso

          out(ii) = .true.
          seg(ii)%voiso        = ii
          seg(ii)%voise        = ii

          ii=idestru
      enddo

!       out(i1o) = .true.
!       seg(i1o)%voiso        = i1o
!       seg(i1o)%voise        = i1o

      seg(itemp)%surface      = IUN
      seg(itemp)%VarFreePlan  = j
      seg(itemp)%vnno         = NSEGMAX

    endif

    ! Segments on the O side of i can be eliminated
    if (seg(itemp)%voiso /= i1otemp) then
      ii=seg(itemp)%voiso
      counter = IZERO
      do while (ii /= i1otemp)
        if (counter > 10) stop "Problem do-while loop surf correct itemp-t1otemp kneecaps SURTYPE IQUATRE/IDOUZE"
        counter = counter + IUN
        irot=seg(ii)%voiso
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        !seg(ii)%surface=IZERO  !We want to keep surface information in case of surface correct junction segment
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd          = izero
        out(ii)             = .true.
        ii = irot
      enddo
    endif

    ! Segment i is modified to account for the surface
    VLI = seg(itemp)%veclin
    seg(itemp)%surface = IUN
    Oi = seg(itemp)%O

    if (SPR < IZERO) then
      tmpdistdiff=huge(tmpDistDiff)
      ! New free surface for the i1e segment
      do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

          if(dot_product(BVECLIN(:,VLI),Plane_MillerI(1:3,jj)) >= IZERO) cycle

          Intersec(1:3) =InterPlanSeg(Plane_MillerR(1:3,jj),Plane_pos(jj),real(Oi,DP),BVECLIN(:,VLI))
          VecO(1:3)=real(Intersec(:)-Oi,DP)

          if (dot_product(VecO,BVECLIN(:,VLI)) < ZERO) cycle

          norm = norvect(VecO)


          if (norm > seg(itemp)%norme*norivect(bveclin(:,VLI))) cycle

          !print *, norm, seg(i)%norme*norivect(bveclin(:,VLI))

          call check_boundaries(Plane_dom(jj),Intersec,outside)
          !print *
          !print *,Intersec,outside,jj

          if (outside) cycle

          if(norm < tmpdistdiff) then
            tmpdistdiff=norm
            j = jj
          endif

        enddo
        seg(itemp)%dom=Plane_dom(j)
    endif ! j<IZERO

    seg(itemp)%varfreeplan=j

    !intersection with the free surface
    Intersec(1:3) = InterPlanSeg(Plane_MillerR(1:3,j),Plane_pos(j),real(Oi,DP),BVECLIN(:,VLI)) !intersection with the free surface
    VecO(1:3) = real(Oi(:)-Intersec(:),DP)
    !We look for the new length of segment i
    ii = -IUN
    tempOi(1:3) = Oi(:)
    do while (dot_product(VecO(1:3),Plane_MillerR(1:3,j)) >= ZERO)
      ii = ii + IUN
      if (ii > 10000) then
        print*,"pb of infinit do while loop 5 in surface_correct",  &
              dot_product(VecO(1:3),Plane_MillerR(1:3,j)), "kk=",kk, 'segment=',itemp
        stop
      endif
      tempOi(1:3) = Oi(:)
      Oi(1:3) = Oi(:)+BVECLIN(:,VLI)
      VecO(1:3) = real(Oi(:)-Intersec(:),DP)
    enddo
    ! update of the length and the origin: now the segment i is at the free surface
    seg(itemp)%norme = seg(itemp)%norme-ii
    seg(itemp)%O(1:3)= tempOi(1:3)
    ! The vnno and vnne of i may be modified to account for the displacement of i
    if (seg(seg(itemp)%voiso)%norme /= IZERO) seg(itemp)%vnno = seg(itemp)%voiso
    if (seg(seg(itemp)%voise)%norme /= IZERO) seg(itemp)%vnne = seg(itemp)%voise

    ! i1e case IDOUZE GD? YEs : eliminate GD
    if (SURFTYPE == IDOUZE .and. seg(itemp)%gd > izero ) then
      ! Change itemp system to voise system
      seg(itemp)%veclin = SegDevIJ(VLI,seg(seg(itemp)%voise)%Veclin)
      !VLI has changed !
      VLI = seg(itemp)%veclin
      ! Remove itemp GD
      seg(itemp)%gd = izero
      !Minimize kneecap number
      call connecij(itemp, seg(itemp)%vnne,3571)
    endif

    ! A kneecap is added at the origin of the segment i/i1O at the surface
    plusseg = plusseg + IUN
    iNew = NSEGM + plusseg
    seg(itemp)%voiso = iNew
    seg(itemp)%vnno = nsegmax
    call init_seg(iNew)
    seg(iNew)%O = tempOi(1:3)
    seg(iNew)%voise = itemp
    seg(iNew)%vnne = itemp
    seg(iNew)%veclin = nbrot(VLI,VLI,2)
    seg(iNew)%dom = seg(itemp)%dom

    condition = .False.

    if (seg(i1otemp)%vnno == NSEGMAX ) then
       condition = .True.
    endif

    if (condition) then
      ii=i1otemp
      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Condition do-while loop surf correct SURTYPE IQUATRE/IDOUZE"
        counter = counter + IUN
        irot = seg(ii)%voiso
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        !seg(ii)%surface=IZERO  !We want to keep surface information in case of surface correct junction segment
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        ii = irot
      enddo

    else
      VLI=seg(i1otemp)%veclin !the segment i belongs to a loop, so the vnno is going out of the free surface too, but with the end
      Oi1o=seg(i1otemp)%O
      Ei1o=Oi1o(1:3)+seg(i1otemp)%norme*BVECLIN(:,VLI)
      seg(i1otemp)%surface=IDEUX

      if (SPR < IZERO) then
        tmpdistdiff=huge(tmpDistDiff)
        ! New free surface for the i1e segment
        do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

          if(dot_product(BVECLIN(:,VLI),Plane_MillerI(1:3,jj)) <= IZERO) cycle

          normsurf(1:3)=Plane_MillerR(1:3,jj)
          Intersec(1:3) = InterPlanSeg(normsurf,Plane_pos(jj),real(Ei1o,DP),BVECLIN(:,VLI))
          VecO(1:3)=real(Ei1o(:)-Intersec(:))

          if (dot_product(VecO,BVECLIN(:,VLI)) < IZERO) cycle

          norm = norvect(VecO)

          if (norm > seg(i1otemp)%norme*norivect(bveclin(:,VLI))) cycle

          !print *, norm, seg(i1o)%norme*norivect(bveclin(:,VLI))

          call check_boundaries(Plane_dom(jj),Intersec,outside)
          !print *
          !print *,Intersec,outside,jj
          if (outside) cycle

          if(norm < tmpdistdiff) then
            tmpdistdiff=norm
            j = jj
          endif

        enddo

      endif ! j<IZERO


      seg(i1otemp)%VarFreePlan = j !surface touched
      Intersec(1:3) =InterPlanSeg(Plane_MillerR(1:3,j),plane_pos(j),real(Oi1o,DP),BVECLIN(:,VLI))
      VecE(1:3) = real(Ei1o(:)-Intersec(:),DP)

      ii=-IUN
      tempEi(1:3)=Ei1o(:)
      do while (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) >= ZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 6 in surface_correct",  &
                dot_product(VecE(1:3),Plane_MillerR(1:3,j)), "kk=",kk, 'segment=', i1otemp
          stop
        endif
        tempEi(1:3)=Ei1o(:)
        Ei1o(1:3)=Ei1o(:)-BVECLIN(:,VLI)
        VecE(1:3)=real(Ei1o(:)-Intersec(:),DP)
      enddo

      ! i1o case IQUATRE GD? YEs : eliminate GD
      if (SURFTYPE == IQUATRE .and. seg(i1otemp)%gd > izero ) then
        ! Change itemp system to voiso system
        seg(i1otemp)%veclin = SegDevIJ(VLI,seg(seg(i1otemp)%voiso)%Veclin)
        !VLI has changed !
        VLI = seg(i1otemp)%veclin
        ! Remove itemp GD
        seg(i1otemp)%gd = izero
        !Minimize kneecap number
        call connecij(seg(i1otemp)%vnno,i1otemp,3572)
      endif

      seg(i1otemp)%norme=seg(i1otemp)%norme-ii
      plusseg=plusseg+ IUN
      iNew=NSEGM + plusseg
      seg(i1otemp)%voise = iNew! kneecap added
      seg(i1otemp)%vnne = nsegmax
      call init_seg(iNew)
      seg(iNew)%O = tempEi(1:3) !updte of the list of neighbour
      seg(iNew)%voiso = i1otemp
      seg(iNew)%vnno = i1otemp
      seg(iNew)%dom = Plane_dom(j)
      seg(iNew)%veclin = nbrot(VLI,VLI,2)
    endif

    if (seg(itemp)%voise /= i1etemp .and. .not. out(i1etemp)) then
      ii=seg(itemp)%voise
      counter = IZERO
      do while (ii /= i1etemp)
        if (counter > 10) stop "Problem loop surf correct at end SURTYPE IQUATRE/IDOUZE"
        counter = counter + IUN
        irot=seg(ii)%voise
        seg(ii)%O=seg(i1etemp)%O
        ii=irot
      enddo
    endif

  !***********
  case (ICINQ,ITREIZE) !!as in  case 4 but now i is going out with the end

    if (seg(I)%surface == ICINQ) then

      itemp = i
      i1otemp = i1o
      i1etemp = i1e

    elseif (seg(I)%surface == ITREIZE) then

      itemp = i1o
      i1etemp = i
      i1otemp = seg(i1o)%voiso

      seg(i1etemp)%surface      = IZERO
      seg(i1etemp)%VarFreePlan  = IZERO
      seg(i1etemp)%voise        = NSEGMAX
      seg(i1etemp)%vnne         = NSEGMAX

      !Cleaning i1e side
      ii=I1e
      counter = IZERO

      do while (ii /= NSEGMAX)
          if (counter > 3) then
            call seginfo(ii,"Infinite do-while loop surf correct cleaning ii SURTYPE ITREIZE")
            call seginfo(i1e,"Infinite do-while loop surf correct cleaning i1e SURTYPE ITREIZE")
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
            stop "Infinite do-while loop surf correct cleaning i1e SURTYPE ITREIZE"
          endif
          counter = counter + IUN

          idestru=seg(ii)%voise

          out(ii) = .true.
          seg(ii)%voiso        = ii
          seg(ii)%voise        = ii

          ii=idestru
      enddo
!
!       out(i1e) = .true.
!       seg(i1e)%voiso        = i1e
!       seg(i1e)%voise        = i1e

      seg(itemp)%surface      = IUN
      seg(itemp)%VarFreePlan  = j
      seg(itemp)%vnne         = NSEGMAX

    endif

!We erase zero neighboring segments
    if (seg(itemp)%voise /= i1etemp) then
      ii=seg(itemp)%voise
      counter = IZERO
      do while (ii /= i1etemp)
        if (counter > 10) stop "Infinite do-while loop surf correct cleaning i1etemp-i1e SURTYPE ICINQ/ITREIZE"
        counter = counter + IUN
        irot = seg(ii)%voise
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        seg(ii)%surface     = izero
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd    = izero
        seg(ii)%dom   = izero
        out(ii) = .true.
        ii = irot
      enddo
    endif

    DomE = seg(i1etemp)%dom

    condition = .False.

    if (seg(i1etemp)%vnne == NSEGMAX .or. seg(seg(i1etemp)%vnne)%surface == IDEUX) then
      condition = .True.
    endif

    if (condition) then

      ii=i1etemp

      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Condition do-while loop surf correct SURTYPE ICINQ/ ITREIZE"
        counter = counter + IUN
        irot = seg(ii)%voise
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        !seg(ii)%surface=IZERO    !We want to keep surface information in case of surface correct junction segment
        seg(ii)%VarFreePlan = izero
        seg(ii)%dom         = izero
        seg(ii)%gd          = izero
        out(ii) = .true.
        ii = irot
      enddo

    else
      VLI=seg(i1etemp)%veclin
      seg(i1etemp)%surface=IUN
      Oi1e=seg(i1etemp)%O

      if (SPR < IZERO) then
        tmpdistdiff=huge(tmpDistDiff)
        ! New free surface for the i1e segment
        do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

          if(dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,jj)) >= IZERO) cycle

          normsurf = Plane_millerR(1:3,jj)
          Intersec(1:3) = InterPlanSeg(normsurf,Plane_pos(jj),real(Oi1e,DP),BVECLIN(:,VLI))
          VecO(1:3)=real(Intersec(:)-Oi1e(:),DP)

          if (dot_product(VecO,BVECLIN(:,VLI)) < ZERO) cycle

          norm = norvect(VecO)

          if (norm > seg(i1etemp)%norme*norivect(bveclin(:,VLI))) cycle

          !print *, norm, seg(i1e)%norme*norivect(bveclin(:,VLI))

          call check_boundaries(Plane_dom(jj),Intersec,outside)
          !print *,Intersec,outside,jj

          if (outside) cycle

          if(norm < tmpdistdiff) then
            tmpdistdiff=norm
            j = jj
          endif
          seg(i1etemp)%dom=Plane_dom(j)

        enddo

      endif ! j<IZERO




      seg(i1etemp)%VarFreePlan=j
      Intersec(1:3) =InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),real(Oi1e,DP),BVECLIN(:,VLI))
      VecO(1:3) = real(Oi1e(:)-Intersec(:),DP)
      ii=-IUN
      tempOi(1:3)=Oi1e(:)
      do while (dot_product(VecO(1:3),Plane_millerR(1:3,j)) >= ZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 7 in surface_correct - DotProd =", &
                dot_product(VecO(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=',i1e
          stop
        endif
        tempOi(1:3)=Oi1e(:)
        Oi1e(1:3)=Oi1e(:)+BVECLIN(:,VLI)
        VecO(1:3)=real(Oi1e(:)-Intersec(:),DP)
      enddo

      ! i case ICINQ, i1e GD? YEs : eliminate GD
      if (SURFTYPE == ICINQ .and. seg(i1etemp)%gd > izero) then
        ! Change itemp system to voise system
        seg(i1etemp)%veclin = SegDevIJ(VLI,seg(seg(i1etemp)%voise)%Veclin)
        !VLI has changed !
        VLI = seg(i1etemp)%veclin
        ! Remove itemp GD
        seg(i1etemp)%gd = izero
        !Minimize kneecap number
        call connecij(i1etemp, seg(i1etemp)%vnne,3573)
      endif

      seg(i1etemp)%norme=seg(i1etemp)%norme-ii
      seg(i1etemp)%O(1:3)=tempOi(1:3)
      plusseg=plusseg+ IUN
      iNew=NSEGM + plusseg
      seg(i1etemp)%voiso = iNew
      seg(i1etemp)%vnno = nsegmax
      call init_seg(iNew)
      seg(iNew)%O = tempOi(1:3)
      seg(iNew)%voise = i1etemp
      seg(iNew)%vnne = i1etemp
      seg(iNew)%dom = seg(i1etemp)%dom
      seg(iNew)%veclin = nbrot(VLI,VLI,2)
    endif

    VLI=seg(itemp)%veclin
    seg(itemp)%surface=IDEUX
    Oi=seg(itemp)%O
    Ei=Oi(1:3)+seg(itemp)%norme*BVECLIN(:,VLI)

    if (SPR < IZERO) then
      tmpdistdiff=huge(tmpDistDiff)
      ! New free surface for the i1e segment
      do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

        if(dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,jj)) <= IZERO) cycle


        normsurf(1:3)=Plane_millerR(1:3,jj)
        Intersec(1:3) =InterPlanSeg(normsurf,Plane_pos(jj),real(Ei,DP),BVECLIN(:,VLI))
        VecO(1:3)=real(Ei(:)-Intersec(:),DP)

        if (dot_product(VecO,BVECLIN(:,VLI)) < ZERO) cycle

        norm = norvect(VecO)

        if (norm > seg(itemp)%norme*norivect(bveclin(:,VLI))) cycle

        !print *, norm, seg(i)%norme*norivect(bveclin(:,VLI))

        call check_boundaries(Plane_dom(jj),Intersec,outside)
        !print *,Intersec,outside,jj
        if (outside) cycle

        if(norm < tmpdistdiff) then
          tmpdistdiff=norm
          j = jj
        endif

      enddo

      endif ! j<IZERO

    seg(itemp)%varfreeplan=j
    Intersec(1:3) =InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),Real(Ei,DP),BVECLIN(:,VLI))
    VecE(1:3) = real(Ei(:)-Intersec(:),DP)
    ii=-IUN
    tempEi(1:3)=Ei(:)
    do while (dot_product(VecE(1:3),Plane_millerR(1:3,j)) >= ZERO)
      ii=ii+IUN
      if (ii > 10000) then
        print*,"pb of infinit do while loop 8 in surface_correct",  &
              dot_product(VecE(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=', itemp
        stop
      endif
      tempEi(1:3)=Ei(:)
      Ei(1:3)=Ei(:)-BVECLIN(:,VLI)
      VecE(1:3)=real(Ei(:)-Intersec(:),DP)
    enddo

    ! i case ITREIZE, i1o GD? YEs : eliminate GD
    if (SURFTYPE == ITREIZE .and. seg(itemp)%gd > izero) then
      ! Change itemp system to voise system
      seg(itemp)%veclin = SegDevIJ(VLI,seg(seg(itemp)%voiso)%Veclin)
      !VLI has changed !
      VLI = seg(itemp)%veclin
      ! Remove itemp GD
      seg(itemp)%gd = izero
      !Minimize kneecap number
      call connecij(seg(itemp)%vnno, itemp,3574)
    endif

    seg(itemp)%norme=seg(itemp)%norme-ii
    ! The vnno of i must sometime be modified to account for the displacement of i
    if (seg(seg(itemp)%voiso)%norme /= IZERO) seg(itemp)%vnno = seg(itemp)%voiso
    plusseg=plusseg+ IUN
    iNew=NSEGM + plusseg
    seg(itemp)%voise = iNew
    seg(itemp)%vnne = nsegmax
    call init_seg(iNew)
    seg(iNew)%O =Oi(1:3)+seg(itemp)%norme*BVECLIN(:,seg(itemp)%veclin)
    seg(iNew)%voiso = itemp
    seg(iNew)%vnno = itemp
    seg(iNew)%dom = Plane_dom(j)
    seg(iNew)%veclin = nbrot(VLI,VLI,2)

   if (seg(i)%voiso /= i1otemp .and. .not. out(i1otemp)) then
      ii=seg(itemp)%voiso
      counter = IZERO
      do while (ii /= i1otemp)
        if (counter > 10) stop "Infinite do-while loop surf correct at end SURTYPE ICINQ/ITREIZE"
        counter = counter + IUN
        irot=seg(ii)%voiso
        seg(ii)%O=seg(itemp)%O
        ii=irot
      enddo
    endif


  !**********
  case (ISIX) !the segment i is going out with both extremity, i will be eliminated, i1o e i1e cut to touch the surface.

    if (seg(i)%voise /= i1e) then
      ii=seg(i)%voise
      counter = IZERO
      do while (ii /= i1e)
        if (counter > 10) stop "Problem do-while loop surf correct ISIX n1"
        counter = counter + IUN
        irot = seg(ii)%voise
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        seg(ii)%surface     = izero
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        ii = irot
      enddo

    endif

    if (seg(i)%voiso /=i1o) then
      ii=seg(i)%voiso
      counter = IZERO
      do while (ii /= i1o)
        if (counter > 10) stop "Problem do-while loop surf correct ISIX n2"
        counter = counter + IUN
        irot = seg(ii)%voiso
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        seg(ii)%surface     = izero
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        ii = irot
      enddo
    endif

    seg(i)%voiso  = i
    seg(i)%voise  = i
    seg(i)%vnno   = i
    seg(i)%vnne   = i
    seg(i)%surface      = izero
    seg(i)%gd           = izero
    seg(i)%VarFreePlan  = izero !this segment will be eliminated
    out(i)  = .true.

    if (seg(i1e)%vnne == NSEGMAX) then

      ii=i1e
      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Problem do-while loop surf correct ISIX n3"
        counter = counter + IUN
        irot = seg(ii)%voise
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        seg(ii)%surface     = izero
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        ii = irot
      enddo

    else
      VLI=seg(i1e)%veclin
      seg(i1e)%surface=IUN
      Oi1e=seg(i1e)%O

      if (SPR < IZERO) then
        tmpdistdiff=huge(tmpDistDiff)
        ! New free surface for the i1e segment
        do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

          if(dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,jj)) >= IZERO) cycle

          normsurf = Plane_millerR(1:3,jj)
          Intersec(1:3) =InterPlanSeg(normsurf,Plane_pos(jj),real(Oi1e,DP),BVECLIN(:,VLI))
          VecO(1:3)=real(Intersec(:)-Oi1e(:),DP)

          if (dot_product(VecO,BVECLIN(:,VLI)) < IZERO) cycle

          norm = norvect(VecO)


          if (norm > seg(i1e)%norme*norivect(bveclin(:,VLI))) cycle

          !print *, norm, seg(i1e)%norme*norivect(bveclin(:,VLI))

          call check_boundaries(Plane_dom(jj),Intersec,outside)
          !print *,Intersec,outside,jj

          if (outside) cycle

          if(norm < tmpdistdiff) then
            tmpdistdiff=norm
            j = jj
          endif
          seg(i1e)%dom=Plane_dom(j)

        enddo

      endif ! j<IZERO

      seg(i1e)%varfreeplan=j
      Intersec(1:3) = InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),real(Oi1e,DP),BVECLIN(:,VLI))
      VecO(1:3) = real(Oi1e(:)-Intersec(:),DP)
      ii=-IUN
      do while (dot_product(VecO(1:3),Plane_millerR(1:3,j)) >= ZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 9 in surface_correct",  &
               dot_product(VecO(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=',i1e
          stop
        endif
        tempOi(1:3)=Oi1e(:)
        Oi1e(1:3)=Oi1e(:)+BVECLIN(:,VLI)
        VecO(1:3)=real(Oi1e(:)-Intersec(:),DP)
      enddo
      if (seg(seg(i1e)%voise)%vnno/=i1e) then
        seg(seg(i1e)%voise)%vnno=i1e
      endif

      ! i case ISIX,  i1e GD? YEs : eliminate GD
      if (seg(i1e)%gd > izero) then
        ! Change itemp system to voise system
        seg(i1e)%veclin = SegDevIJ(VLI,seg(seg(i1e)%voise)%Veclin)
        !VLI has changed !
        VLI = seg(i1e)%veclin
        ! Remove itemp GD
        seg(i1e)%gd = izero
        !Minimize kneecap number
        call connecij(i1e, seg(i1e)%vnne,3575)
      endif

      seg(i1e)%norme=seg(i1e)%norme-ii
      seg(i1e)%O(1:3)=tempOi(1:3)
      plusseg=plusseg+ IUN
      iNew=NSEGM + plusseg
      seg(i1e)%voiso = iNew
      seg(i1e)%vnno = nsegmax
      call init_seg(iNew)
      seg(iNew)%O = tempOi(1:3)
      seg(iNew)%voise = i1e
      seg(iNew)%vnne = i1e
      seg(iNew)%dom = seg(i1e)%dom
      seg(iNew)%veclin = nbrot(VLI,VLI,2)
    endif

    if (seg(i1o)%vnno == NSEGMAX) then

      ii=i1o
      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Problem do-while loop surf correct ISIX n4"
        counter = counter + IUN
        irot = seg(ii)%voiso
        seg(ii)%voiso = ii
        seg(ii)%vnno  = ii
        seg(ii)%voise = ii
        seg(ii)%vnne  = ii
        seg(ii)%surface     = izero
        seg(ii)%VarFreePlan = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        ii = irot
      enddo

    else
      VLI=seg(i1o)%veclin
      seg(i1o)%surface=IDEUX
      Oi1o=seg(i1o)%O
      Ei1o=Oi1o(1:3)+seg(i1o)%norme*BVECLIN(:,VLI)

      if (SPR < IZERO) then
        tmpdistdiff=huge(tmpDistDiff)
        ! New free surface for the i1e segment
        do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

          if(dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,jj)) <= IZERO) cycle


          normsurf(1:3)=Plane_millerR(1:3,jj)
          Intersec(1:3) =InterPlanSeg(normsurf,Plane_pos(jj),real(Ei1o,DP),BVECLIN(:,VLI))
          VecO(1:3)=real(Ei1o(:)-Intersec(:),DP)

          if (dot_product(VecO,BVECLIN(:,VLI)) < IZERO) cycle

          norm = norvect(VecO)

          if (norm > seg(i1o)%norme*norivect(bveclin(:,VLI))) cycle


           ! print *, norm, seg(i1o)%norme*norivect(bveclin(:,VLI))

          call check_boundaries(Plane_dom(jj),Intersec,outside)
          if (outside) cycle

          if(norm < tmpdistdiff) then
            tmpdistdiff=norm
            j = jj
          endif

        enddo

      endif ! j<IZERO

      seg(i1o)%VarFreePlan = j

      Intersec(1:3) = InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),real(Oi1o,DP),BVECLIN(:,VLI))
      VecE(1:3) = real(Ei1o(:)-Intersec(:),DP)
      ii=-IUN
      do while (dot_product(VecE(1:3),Plane_millerR(1:3,j)) >= IZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 10 in surface_correct",  &
                dot_product(VecE(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=', i1o
          stop
        endif
        tempEi(1:3)=Ei1o(:)
        Ei1o(1:3)=Ei1o(:)-BVECLIN(:,VLI)
        VecE(1:3)=real(Ei1o(:)-Intersec(:),DP)
      enddo
      if (seg(seg(i1o)%voiso)%vnne/=i1o) then
        seg(seg(i1o)%voiso)%vnne=i1o
      endif

      ! i case ISIX,  i1o GD? YEs : eliminate GD
      if (seg(i1o)%gd > izero) then
        ! Change itemp system to voise system
        seg(i1o)%veclin = SegDevIJ(VLI,seg(seg(i1o)%voiso)%Veclin)
        !VLI has changed !
        VLI = seg(i1o)%veclin
        ! Remove itemp GD
        seg(i1o)%gd = izero
        !Minimize kneecap number
        call connecij( seg(i1o)%vnno, i1o,3576)
      endif

      seg(i1o)%norme=seg(i1o)%norme-ii
      plusseg=plusseg+ IUN
      iNew=NSEGM + plusseg
      seg(i1o)%voise = iNew
      seg(i1o)%vnne = nsegmax
      call init_seg(iNew)
      seg(iNew)%O = tempEi(1:3)
      seg(iNew)%voiso = i1o
      seg(iNew)%vnno = i1o
      seg(iNew)%dom = Plane_dom(j)
      seg(iNew)%veclin = nbrot(VLI,VLI,2)
    endif

    !call seginfo(i1o,'check new procedure case 6 surface_correct i1o')
    !call seginfo(i1e,'check new procedure case 6 surface_correct i1e')
    !stop

  !***********
  case (ISEPT) !segment touching the surface with the origin and crossing surface with the end

    if (seg(i1e)%surface > IZERO) then
      ! Destruction of a micro-loop
      ii=i1o
      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Problem do-while loop surf correct ISEPT n1"
        counter = counter + IUN
        i1etemp = seg(ii)%voise
        seg(ii)%voiso = ii
        seg(ii)%voise = ii
        seg(ii)%VarFreeplan = izero
        seg(ii)%surface     = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        if (kkdebug) write (379,fmt='(x,"Delete microloop ISEPT segment ",x,I9,":" )') ii
        ii = i1etemp
      enddo

    elseif (seg(i1e)%norme /= IZERO) then
      VLI=seg(i1e)%veclin

      ! i case ISEPT,  i1e GD? YEs : eliminate GD
      if (seg(i1e)%gd > izero) then
        ! Change itemp system to voise system
        seg(i1e)%veclin = SegDevIJ(VLI,seg(seg(i1e)%voise)%Veclin)
        seg(i)%veclin = conec(seg(i1e)%veclin,1)
        !VLI has changed !
        VLI = seg(i1e)%veclin
        ! Remove itemp GD
        seg(i1e)%gd = izero
        !Minimize kneecap number
        call connecij( i1e, seg(i1e)%vnne ,3578)
      endif

      ! i1e is the new segment touching the surface
      seg(i1e)%surface = IUN
      seg(i1e)%VarFreePlan = seg(i)%VarFreePlan
      singulier(i1e) = .True.
      seg(i1e)%vnno = NSEGMAX

      Oi1e=seg(i1e)%O

      Intersec(1:3) =InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),Real(Oi1e,DP),BVECLIN(:,VLI))
      VecO(1:3) = real(Oi1e(:)-Intersec(:),DP)
      ii=-IUN
      do while (dot_product(VecO(1:3),Plane_millerR(1:3,j)) >= IZERO)
        ii=ii+IUN
        if (ii > 10000) then
           print*,"pb of infinit do while loop 11 in surface_correct",  &
                dot_product(VecO(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=',i1e
           stop
        endif
        tempOi(1:3)=Oi1e(:)
        Oi1e(1:3)=Oi1e(:)+BVECLIN(:,VLI)
        VecO(1:3)=real(Oi1e(:)-Intersec(:),DP)
      enddo
      seg(i1e)%norme=seg(i1e)%norme-ii
      seg(i1e)%O(1:3)=tempOi(1:3)

      ! i becomes kneecap of i1e
      seg(i)%VarFreeplan  = izero
      seg(i)%surface      = izero
      seg(i)%norme        = izero
      seg(i)%O            = seg(i1e)%O
      seg(i)%voiso        = NSEGMAX
      seg(i)%vnno         = NSEGMAX
      seg(i)%vnne         = i1e
      singulier(i)        = .False.
      seg(seg(i1e)%vnne)%vnno = i1e
      ! i1o can be elimanted
      seg(i1o)%voiso      = i1o
      seg(i1o)%voise      = i1o
      seg(i1o)%surface    = izero
      seg(i1o)%VarFreePlan = izero
      seg(i1o)%gd         = izero
      out(i1o) = .true.
    else
      VLI = seg(seg(i1e)%vnne)%veclin
      ! i1e is a kneecap, we must consider the vnne segments instead
      seg(seg(i1e)%vnne)%surface      =IUN
      seg(seg(i1e)%vnne)%VarFreePlan  = seg(i)%VarFreePlan
      seg(seg(i1e)%vnne)%vnno         = NSEGMAX
      singulier(seg(i1e)%vnne)        = .True.

      Oi1e = seg(seg(i1e)%vnne)%O

      Intersec(1:3) =InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),real(Oi1e,DP),BVECLIN(:,VLI))
      VecO(1:3) = real(Oi1e(:)-Intersec(:),DP)
      ii=-IUN
      do while (dot_product(VecO(1:3),Plane_millerR(1:3,j)) >= IZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 12 in surface_correct",  &
                dot_product(VecO(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=', i1o
          stop
        endif
        tempOi(1:3)=Oi1e(:)
        Oi1e(1:3)=Oi1e(:)-BVECLIN(:,VLI)
        VecO(1:3)=real(Oi1e(:)-Intersec(:),DP)
      enddo
      seg(seg(i1e)%vnne)%norme=seg(seg(i1e)%vnne)%norme-ii
      seg(seg(i1e)%vnne)%O(1:3)=tempOi(1:3)

      ! A surface kneecap is redefined
      i1e=seg(seg(i1e)%vnne)%voiso
      ii=seg(i1e)%voiso
      seg(i1e)%O=tempOi(1:3)
      seg(i1e)%VarFreeplan = IZERO
      seg(i1e)%surface = IZERO
      seg(i1e)%norme = IZERO
      seg(i1e)%voiso = NSEGMAX
      seg(i1e)%vnno = NSEGMAX
      ! All the other segments must be eliminated
      counter = IZERO
      do while (ii /= i)
        if (counter > 10) stop "Problem do-while loop surf correct ISEPT n2"
        counter = counter + IUN
        i1e = seg(ii)%voiso
        seg(ii)%voiso = ii
        seg(ii)%voise = ii
        seg(ii)%surface     = izero
        seg(ii)%gd          = izero
        seg(ii)%VarFreeplan = izero
        out(ii) = .true.
        ii = i1e
      enddo
      seg(i)%norme          = izero
      seg(i)%voiso          = i
      seg(i)%voise          = i
      seg(i)%VarFreeplan    = izero
      seg(i)%surface        = izero
      seg(i)%gd             = izero
      out(i) = .true.
      seg(i1o)%voiso        = i1o
      seg(i1o)%voise        = i1o
      seg(i1o)%VarFreeplan  = izero
      seg(i1o)%surface      = izero
      seg(i1o)%gd           = izero
      out(i1o) = .true.
    endif


  !***********
  case (IHUIT) !segment touching the surface with the end and crossing surface with the origin

    if (seg(i1o)%surface > IZERO) then
      ! Destruction of a micro-loop
      ii=i1e
      Counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Problem do-while loop surf correct IHUIT n1"
        counter = counter + IUN
        i1otemp = seg(ii)%voiso
        seg(ii)%voiso       = ii
        seg(ii)%voise       = ii
        seg(ii)%VarFreeplan = izero
        seg(ii)%surface     = izero
        seg(ii)%dom         = izero
        seg(ii)%gd          = izero
        out(ii) = .true.
        if (kkdebug) write (379,fmt='(x,"Delete microloop IHUIT segment ",x,I9,":" )') ii
        ii = i1otemp
      enddo

    elseif (seg(i1o)%norme /=IZERO) then
      VLI=seg(i1o)%veclin

      ! i case IHUIT,  i1o GD? YEs : eliminate GD
      if (seg(i1o)%gd > izero) then
        ! Change itemp system to voiso system
        seg(i1o)%veclin = SegDevIJ(VLI,seg(seg(i1o)%voiso)%Veclin)
        seg(i)%veclin=conec(seg(i1o)%veclin,1)
        !VLI has changed !
        VLI = seg(i1o)%veclin
        ! Remove itemp GD
        seg(i1o)%gd = izero
        !Minimize kneecap number
        call connecij(seg(i1o)%vnno, i1o, 3579)
      endif

      ! i1o is the new segment touching the surface
      seg(i1o)%surface = IDEUX
      seg(i1o)%VarFreePlan = seg(i)%VarFreePlan
      seg(i1o)%vnne = NSEGMAX
      ! i becomes kneecap of i1o
      seg(i)%VarFreeplan = IZERO
      seg(i)%surface = IZERO
      seg(i)%norme = IZERO
      seg(i)%vnno = i1o
      seg(i)%voise = NSEGMAX

      Oi1o=seg(i1o)%O
      Ei1o=Oi1o(1:3)+seg(i1o)%norme*BVECLIN(:,VLI)

      Intersec(1:3) =InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),Real(Ei1o,DP),BVECLIN(:,VLI))
      VecE(1:3) = real(Ei1o(:)-Intersec(:),DP)
      ii=-IUN
      do while (dot_product(VecE(1:3),Plane_millerR(1:3,j)) >= ZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 13 in surface_correct",  &
                dot_product(VecE(1:3),Plane_millerR(1:3,j)), "kk=",kk,'segment=', i1o
          stop
        endif
        tempEi(1:3)=Ei1o(:)
        Ei1o(1:3)=Ei1o(:)-BVECLIN(:,VLI)
        VecE(1:3)=real(Ei1o(:)-Intersec(:),DP)
      enddo
      seg(i1o)%norme=seg(i1o)%norme-ii
      seg(i)%O=Oi1o(1:3)+seg(i1o)%norme*BVECLIN(:,VLI)
      seg(i)%vnne = NSEGMAX
      seg(seg(i1o)%vnno)%vnne=i1o

      ! i1e can be elimanted
      seg(i1e)%voiso        = i1e
      seg(i1e)%voise        = i1e
      seg(i1e)%surface      = izero
      seg(i1e)%VarFreePlan  = izero
      seg(i1e)%gd           = izero
      out(i1e) = .true.

    else
      VLI=seg(seg(i1o)%vnno)%veclin
      ! i1o is a kneecap, we must consider the vnno segments instead
      seg(seg(i1o)%vnno)%surface      = IDEUX
      seg(seg(i1o)%vnno)%VarFreePlan  = seg(i)%VarFreePlan
      seg(seg(i1o)%vnno)%vnne         = NSEGMAX

        Oi1o=seg(seg(i1o)%vnno)%O
        Ei1o=Oi1o(1:3)+seg(seg(i1o)%vnno)%norme*BVECLIN(:,VLI)

      if (kkdebug) then
        write(379,*) ''
        write(379,*) 'BVECLIN(:,VLI) ', BVECLIN(:,VLI)
        write(379,*) 'Plane_millerI', Plane_millerI(1:3,j)
        write(379,*) 'dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,j)) ', dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,j))
      endif

      Intersec(1:3) =InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),real(Ei1o,DP),BVECLIN(:,VLI))
      VecE(1:3) = real(Ei1o(:)-Intersec(:),DP)
      ii=-IUN
      do while (dot_product(VecE(1:3),Plane_millerR(1:3,j)) >= ZERO)
        ii=ii+IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop 14 in surface_correct",  &
                dot_product(VecE(1:3),Plane_millerR(1:3,j)), "kk=",kk, 'segment=', i1o
          stop
        endif
        tempEi(1:3)=Ei1o(:)
        Ei1o(1:3)=Ei1o(:)-BVECLIN(:,VLI)
        VecE(1:3)=real(Ei1o(:)-Intersec(:),DP)
      enddo
      seg(seg(i1o)%vnno)%norme=seg(seg(i1o)%vnno)%norme-ii

      ! A surface kneecap is redefined
      i1o=seg(seg(i1o)%vnno)%voise
      ii=seg(i1o)%voise
      seg(i1o)%O=Oi1o(1:3)+seg(seg(i1o)%vnno)%norme*BVECLIN(:,VLI)
      seg(i1o)%VarFreeplan = IZERO
      seg(i1o)%surface = IZERO
      seg(i1o)%norme = IZERO
      seg(i1o)%voise = NSEGMAX
      seg(i1o)%vnne = NSEGMAX
      ! All the other segments must be eliminated
      counter = IZERO
      do while (ii /= i)
        if (counter > 10) stop "Problem do-while loop surf correct IHUIT n2"
        counter = counter + IUN
        i1o = seg(ii)%voise
        seg(ii)%voiso       = ii
        seg(ii)%voise       = ii
        seg(ii)%surface     = izero
        seg(ii)%gd          = izero
        seg(ii)%VarFreeplan = izero
        out(ii) = .true.
        ii = i1o
      enddo
      seg(i)%norme        = izero
      seg(i)%voiso        = i
      seg(i)%voise        = i
      seg(i)%VarFreeplan  = izero
      seg(i)%surface      = izero
      seg(i)%gd           = izero
      out(i) = .true.
      seg(i1e)%voiso        = i1e
      seg(i1e)%voise        = i1e
      seg(i1e)%VarFreeplan  = izero
      seg(i1e)%surface      = izero
      seg(i1e)%gd           = izero
      out(i1e) = .true.
    endif

  !***********
  case (INEUF,IDIX,IONZE) !Segment I crossing 2 fs qnd extremities in two different domains

    if (seg(I)%surface == INEUF) then
      itemp = i
    elseif (seg(I)%surface == IDIX) then
      itemp = i1o
      seg(i)%surface = IZERO
      seg(i)%VarFreePlan = IZERO

      ! i case IDIX,  i1o GD? YEs : eliminate GD
      if (seg(itemp)%gd > izero) then
        if(kkdebug) write(379,*) "i case IDIX,  i1o GD? YEs : eliminate"
        ! Change itemp system to voiso system
        seg(itemp)%veclin = SegDevIJ(seg(itemp)%veclin,seg(seg(itemp)%voiso)%Veclin)
        ! Remove itemp GD
        seg(itemp)%gd = izero
        !Minimize kneecap number
        call connecij(seg(itemp)%vnno, itemp, 3580)
      endif

    elseif   (seg(I)%surface == IONZE) then
      itemp = i1e
      seg(i)%surface = IZERO
      seg(i)%VarFreePlan = IZERO

      ! i case IONZE,  i1e GD? YEs : eliminate GD
      if (seg(itemp)%gd > izero) then
        if(kkdebug) write(379,*) " i case IONZE,  i1e GD? YEs : eliminate GD"
        ! Change itemp system to voise system
        seg(itemp)%veclin = SegDevIJ(seg(itemp)%veclin,seg(seg(itemp)%voise)%Veclin)
        ! Remove itemp GD
        seg(itemp)%gd = izero
        !Minimize kneecap number
        call connecij(itemp, seg(itemp)%vnne, 3581)
      endif

    endif

    Ivnne     = seg(itemp)%VNNE
    Normi     = seg(itemp)%norme
    VLI       = seg(itemp)%veclin
    ve        = seg(itemp)%voise
    Oi(:)     = SEG(itemp)%O(1:3)
    Ei(:)     = Oi(:) + Normi * BVECLIN(:,VLI)
    OiReal    = real(Oi(:),DP)
    Eireal    = real(Ei(:),DP)

    if (kkdebug) then
      write (379,fmt='(" ")')
      write (379,fmt='(">- Surface_correct case ", I4)') SURFTYPE
      write (379,fmt='(">- Segment coordinates before cutting procedure ",6I9)') Oi(:),Ei(:)
      write (379,fmt='(">- Cutting O side -< ")')
    endif

    DistanceS = huge(distanceS)
    DistanceP = huge(distanceS)

    DistDiff  = -huge(distanceS)
    j=-iun

    !--------------------------------------------------------------------
    !-- Seeking the closest freesurface
    tmpdistdiff=huge(tmpDistDiff)
    ! New free surface for the i1e segment
    do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

      if(dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,jj)) <= IZERO) cycle

      normsurf(1:3)=Plane_millerR(1:3,jj)
      Intersec(1:3) =InterPlanSeg(normsurf,Plane_pos(jj),Oireal,BVECLIN(:,VLI))
      VecO(1:3)=Oireal(:)-Intersec(:)

      if (dot_product(VecO,BVECLIN(:,VLI)) > IZERO) cycle

      norm = norvect(VecO)

      if (norm > normi*norivect(bveclin(:,VLI))) cycle


       ! print *, norm, seg(i1o)%norme*norivect(bveclin(:,VLI))

      call check_boundaries(Plane_dom(jj),Intersec,outside)
      if (outside) cycle

      if(norm < tmpdistdiff) then
        tmpdistdiff=norm
        j = jj
      endif

    enddo

    Intersec(1:3) = InterPlanSeg(Plane_millerR(1:3,j),Plane_pos(j),Oireal,BVECLIN(:,VLI))
    ! Distance to the interface
    VecE(1:3) = Eireal(:)-Intersec(:)

    if (dot_product(VecE(1:3),Plane_millerR(1:3,j)) >= ZERO) then
      ! Ei is outside and we look for the closest coordinate to the surface
      ! the segment length is decreased
      ii = IZERO
      do while (dot_product(VecE(1:3),Plane_millerR(1:3,j)) >= IZERO)
        ii = ii + IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop case INEUF 1 in surface_correct",  &
                dot_product(VecE(1:3),Plane_millerR(1:3,j)), "kk=",kk
          stop
        endif
        tempEi(1:3) = Ei(:)
        Ei(1:3) = Ei(:)-BVECLIN(:,VLI)
        if (kkdebug) write(379,*) Ei(:)
        VecE(1:3) = real(Ei(:)-Intersec(:),DP)
      enddo

      if(kkdebug) write(379,*) 'Case ',SURFTYPE,' Surface_correct step 1 :&
        & the segment length after displacement is decreased by', ii
      Ei(:) = tempEi(:)
      seg(itemp)%norme = seg(itemp)%norme - ii

    endif

    if (kkdebug) write (379,fmt='(">- New segment coordinates on O side ",6I9)') Oi(:),Ei(:)

    !Kneecap creation at end
    plusseg=plusseg+IUN
    INew=Nsegm+plusseg
    call init_seg(iNew)

    !Update I
    seg(itemp)%norme        = seg(itemp)%norme + IUN
    seg(itemp)%VOISE        = INew
    seg(itemp)%VNNE         = NSEGMAX
    seg(itemp)%Surface      = IDEUX
    seg(itemp)%VarfreePlan  = J
    singulier(itemp) = .True.

    !Update iNew
    seg(iNew)%O(:)          = Oi(:)+seg(itemp)%norme * BVECLIN(:,VLI)
    seg(iNew)%veclin        = CONEC(VLI,1)
    seg(iNew)%VOISO         = Itemp
    seg(iNew)%VNNO          = Itemp
    seg(iNew)%VOISE         = NSEGMAX
    seg(iNew)%VNNE          = NSEGMAX
    seg(iNew)%dom           = Plane_dom(J)

    !call connecij(I,NSEGMAX,987)
    !call seginfo(I,'Case INEUF,IDIX,IONZE  - New connection to surface - O side')


    !--------------------------------------------------------------------
    !-- Seeking the closest freesurface
    if (kkdebug) then
      write (379,fmt='(" ")')
      write (379,fmt='(">- Surface_correct case ",I4," - E side -< ")') SURFTYPE
    endif


    DistanceS = huge(distanceS)
    DistanceP = huge(distanceS)

    DistDiff  = -huge(distanceS)
    j = -IUN

    tmpdistdiff=huge(tmpDistDiff)

    ! New free surface for the i1e segment
    do jj=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

      if(dot_product(BVECLIN(:,VLI),Plane_millerI(1:3,jj)) >= IZERO) cycle


      normsurf(1:3)=Plane_millerR(1:3,jj)
      Intersec(1:3) =InterPlanSeg(normsurf,Plane_pos(jj),Eireal,BVECLIN(:,VLI))
      VecE(1:3)=Eireal(:)-Intersec(:)

      if (dot_product(VecE,BVECLIN(:,VLI)) < ZERO) cycle

      norm = norvect(VecE)

      if (norm > normi*norivect(bveclin(:,VLI))) cycle

      !print *, norm, seg(i)%norme*norivect(bveclin(:,VLI))

      call check_boundaries(Plane_dom(jj),Intersec,outside)
      !print *,Intersec,outside,jj
      if (outside) cycle

      if(norm < tmpdistdiff) then
        tmpdistdiff=norm
        j = jj
      endif

    enddo

    if (kkdebug) write (379,*) "surface ",j," selected"

    !New segment creation
    plusseg=plusseg+IUN
    INew2=Nsegm+plusseg
    call init_seg(iNew2)

    Oi(:) = Oi(:) + normi * BVEClIN(:,VLI)

    Intersec(1:3) = InterPlanSeg(Plane_MillerR(1:3,j),Plane_pos(j),Eireal,BVECLIN(:,VLI))
    ! Distance to the interface
    VecE(1:3) = Eireal(:)-Intersec(:)

    if (kkdebug) then
      write(379,*)"Eireal",Eireal
      write(379,*)"Intersec", Intersec
      write(379,*)"dot_product(VecE(1:3),Plane_MillerR(1:3,j))",dot_product(VecE(1:3),Plane_MillerR(1:3,j))
    endif

    if (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) < ZERO) then
      ! Oi is outside and we look for the closest coordinate to the surface
      ! the segment length is decreased
      ii = IZERO
      do while (dot_product(VecE(1:3),Plane_MillerR(1:3,j)) < ZERO)
        ii = ii + IUN
        if (ii > 10000) then
          print*,"pb of infinit do while loop case ",SURFTYPE," step 2 in surface_correct",  &
                dot_product(VecE(1:3),Plane_MillerR(1:3,j)), "kk=",kk
          stop
        endif
        tempOi(1:3) = Oi(:)
        Oi(1:3) = Oi(:) - BVECLIN(:,VLI)
        if (kkdebug) write(379,*) Oi(:)
        VecE(1:3) = real(Oi(:)-Intersec(:),DP)
      enddo

      if(kkdebug) write(379,*) 'Case ',SURFTYPE,' Surface_correct step 2 : the segment length after displacement is increase by', ii
      !Oi(:) = tempOi(:)
      seg(iNew2)%norme = ii

    endif

    if (kkdebug) write (379,fmt='(">- New segment coordinates on E side ",6I9)') &
                Oi(:),Oi(:) + seg(iNew2)%norme * BVECLIN(:,VLI)

    !Kneecap creation at Orgin
    plusseg=plusseg+IUN
    INew3=Nsegm+plusseg
    call init_seg(iNew3)

    !Update iNew2
    seg(iNew2)%VECLIN       = VLI
    seg(iNew2)%O(:)         = Oi(:)
    seg(iNew2)%VOISO        = iNew3
    seg(iNew2)%VNNO         = NSEGMAX
    seg(iNew2)%VOISE        = ve
    seg(iNew2)%VNNE         = Ivnne
    seg(iNew2)%Surface      = IUN
    seg(iNew2)%VarfreePlan  = J
    seg(iNew2)%dom          = Plane_dom(J)
    singulier(iNew2) = .True.

    !Update iNew3
    seg(iNew3)%O(:)          = Oi(:)
    seg(iNew3)%veclin        = CONEC(VLI,1)
    seg(iNew3)%VOISO         = NSEGMAX
    seg(iNew3)%VNNO          = NSEGMAX
    seg(iNew3)%VOISE         = inew2
    seg(iNew3)%VNNE          = inew2
    seg(iNew3)%dom           = Plane_dom(J)

    !Update voise
    seg(ve)%Vnno          = Inew2
    seg(ve)%voiso         = Inew2

    seg(Ivnne)%Vnno       = Inew2


    if (ve /= IVNNe .and. seg(ve)%voise /= IVNNE) then
      itemp = seg(ve)%voise

      do while(seg(itemp)%vnno /= seg(IVNNE)%vnno)
         call voisinage(itemp,372)
         itemp = seg(itemp)%voise
      enddo

    endif

    !Update vnne if needed
    !if (seg(i1e)%norme == IZERO) seg(Ivnne)%Vnno = Inew2

    !call seginfo(iNew2,'Case INEUF,IDIX,IONZE - New connection to surface - E side')
    !call connecij(NSEGMAX,iNew2,987)
    !call seginfo(iNew2,'Case INEUF - New connection to surface - E side')

  !***********
  case DEFAULT

    print *,'surface segments with type < 4 must not appear here (surface_correct)'
    print *,i,seg(i)%surface
    stop

end select

if (kkdebug) then
  write(379,*) 'Before gd, jonc surface correct'
  write(379,*) 'segment i =',i,'surface type ',seg(i)%surface,'jonc ',seg(i)%jonc,'gd ',seg(i)%gd,'out',out(i)
  write(379,*) 'segment i1o =',i1o,'surface type',seg(i1o)%surface,'jonc ',seg(i1o)%jonc,'gd ',seg(i1o)%gd,'out',out(i1o)
  write(379,*) 'segment i1e =',i1e,'surface type',seg(i1e)%surface,'jonc ',seg(i1e)%jonc,'gd ',seg(i1e)%gd,'out',out(i1e)
endif

! The case of GD segment is now updated
if (seg(i)%gd > izero) then
  if (kkdebug) write(379,*) 'GD segment surface correct treatment I',i
  call seginfo(I,"  stop GD surface correct I")
  stop "GD surface correct I"
endif

if (seg(i1o)%gd > izero .and. seg(i1o)%surface /= IZERO) then
  if (kkdebug) write(379,*) 'GD segment surface correct treatment I1O',i1o
  call seginfo(I1o,"  stop GD surface correct I")
  stop "GD surface correct I1o"
endif

if (seg(i1e)%gd > izero .and. seg(i1e)%surface /= IZERO) then
  if (kkdebug) write(379,*) 'GD segment surface correct treatment I1e',i1e
  call seginfo(I1e,"  stop GD surface correct I1e")
  stop "GD surface correct I1e"
endif

! The case of Junction segments is now updated
if (seg(i)%jonc .and. (seg(i)%surface > Izero .or. OUT(i))) then

  binome = seg(i)%Ijonc

  if (kkdebug) write(379,*) 'Junction segment surface correct treatment I ',i,&
    & ": norms : I= ", seg(i)%norme, " binome= ",seg(binome)%norme

  !We liberate the junction segment i
  if (seg(i)%norme /= 0) then
    seg(i)%JONC = .false.       !simple case
    seg(i)%IJONC = nsegmax
    seg(i)%tJONC = 0
  else
    call stopjonc(i,111)
  endif

  !We liberate the complementary junction segment k
  if (seg(binome)%norme /= 0) then
    seg(binome)%JONC = .false.       !simple case
    seg(binome)%IJONC = nsegmax
    seg(binome)%tJONC = 0
  else
    call stopjonc(binome,112)
  endif

endif

! The case of Junction segments is now updated
if (seg(i1o)%jonc .and. (seg(i1o)%surface > Izero .or. out(i1o))) then

  binome = seg(i1o)%Ijonc

  if (kkdebug) write(379,*) 'Junction segment surface correct treatment I1o ',i1o,&
    & ": norms : I= ", seg(i1o)%norme, " binome= ",seg(binome)%norme

  !We liberate the junction segment i
  if (seg(i1o)%norme /= 0) then
    seg(i1o)%JONC = .false.       !simple case
    seg(i1o)%IJONC = nsegmax
    seg(i1o)%tJONC = 0
  else
    call stopjonc(i1o,221)
  endif

  !We liberate the complementary junction segment k
  if (seg(binome)%norme /= 0) then
    seg(binome)%JONC = .false.       !simple case
    seg(binome)%IJONC = nsegmax
    seg(binome)%tJONC = 0
  else
    call stopjonc(binome,222)
  endif

endif


! The case of Junction segments is now updated
if (seg(i1e)%jonc .and. (seg(i1e)%surface > Izero .or. OUT(i1e))) then

  binome = seg(i1e)%Ijonc

  if (kkdebug) write(379,*) 'Junction segment surface correct treatment I1e ',i1e,&
    & ": norms : I= ", seg(i1e)%norme, " binome= ",seg(binome)%norme

  !We liberate the junction segment i
  if (seg(i1e)%norme /= 0) then
    seg(i1e)%JONC = .false.       !simple case
    seg(i1e)%IJONC = nsegmax
    seg(i1e)%tJONC = 0
  else
    call stopjonc(i1e,331)
  endif

  !We liberate the complementary junction segment k
  if (seg(binome)%norme /= 0) then
    seg(binome)%JONC = .false.       !simple case
    seg(binome)%IJONC = nsegmax
    seg(binome)%tJONC = 0
  else
    call stopjonc(binome,332)
  endif

endif

!---------------
!Kneecap junction and gd cleaning at surface segment
!---------------

!Segment I
if (seg(i)%surface > IZERO) then

  SURFTYPE = seg(i)%surface
  Condition = .false.
  if (kkdebug) write (379,fmt='(x,"Cleaning surface_correct I-side ",x,I9,":" )') I

  select case (SURFTYPE)

    case (IUN)

      ! E-side cleaning
      ii=i
      Ivnne = seg(i)%vnne

      counter = IZERO
      do while (ii /= Ivnne)

          if (counter > 1000) then
            call seginfo(I,"Infinite do-while loop surf correct cleaning I case IUN")
            write(379,*) "Check that seg(I)%vnne is correct : ", Ivnne, ":"
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
            stop "Infinite do-while loop surf correct cleaning I case IUN"
          endif

          counter = counter + IUN

          i1etemp=seg(ii)%voise

          if (seg(ii)%gd > izero) then
            if (kkdebug) write (379,fmt='(x,"Delete segment GD status",x,I9,":" )') ii
            seg(ii)%gd  = izero
            condition = .true.
          endif

          if (seg(ii)%jonc) then
            binome  = seg(ii)%Ijonc
            if (kkdebug) write (379,fmt='(x,"Delete segment and binome JONC status",x,I9,": binome",x,I9,":" )') ii,binome
            call stopjonc(ii,311)
            call stopjonc(binome,312)
            condition = .true.
          endif
          ii=i1etemp
      enddo

      if (condition) call connecij(I,Ivnne,313)

    case (IDEUX)

      ! E-side cleaning
      ii=i
      Ivnno = seg(i)%vnno

      counter = IZERO
      do while (ii /= Ivnno)

          if (counter > 1000) then
            call seginfo(I,"Infinite do-while loop surf correct cleaning I case IDEUX")
            write(379,*) "Check that seg(I)%vnno is correct : ", Ivnno, ":"
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voiso ",seg(ii)%voiso,":"
            stop "Infinite do-while loop surf correct cleaning i case IDEUX"
          endif

          counter = counter + IUN


          i1otemp=seg(ii)%voiso

          if (seg(ii)%gd > izero) then
            if (kkdebug) write (379,fmt='(x,"Delete segment GD status",x,I9,":" )') ii
            seg(ii)%gd  = izero
            condition = .true.
          endif

          if (seg(ii)%jonc) then
            binome  = seg(ii)%Ijonc
            if (kkdebug) write (379,fmt='(x,"Delete segment and binome JONC status",x,I9,": binome",x,I9,":" )') ii,binome
            call stopjonc(ii,321)
            call stopjonc(binome,322)
            condition = .true.
          endif
          ii=i1otemp
      enddo

      if (condition) call connecij(Ivnno,I,323)

    case DEFAULT

      print *,'surface segments with type < 4 must not appear here (surface_correct)'
      print *,i,seg(i)%surface
      stop

  end select

endif

!Segment I1o
if (seg(I1o)%surface > IZERO .and. .not. out(I1o) ) then

  SURFTYPE = seg(I1o)%surface
  Condition = .false.
  if (kkdebug) write (379,fmt='(x,"Cleaning surface_correct I1o-side ",x,I9,":" )') I1o

  select case (SURFTYPE)

    case (IUN) ! The case should not be observed
      ! E-side cleaning
      ii=I1o
      Ivnne = seg(I1o)%vnne

      counter = IZERO
      do while (ii /= Ivnne)

          if (counter > 1000) then
            call seginfo(i1o,"Infinite do-while loop surf correct cleaning i1o case IUN")
            write(379,*) "Check that seg(I1o)%vnne is correct : ", Ivnne, ":"
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
            stop "Infinite do-while loop surf correct cleaning i1o case IUN"
          endif

          counter = counter + IUN

          i1etemp=seg(ii)%voise

          if (seg(ii)%gd > izero) then
            if (kkdebug) write (379,fmt='(x,"Delete segment GD status",x,I9,":" )') ii
            seg(ii)%gd  = izero
            condition = .true.
          endif

          if (seg(ii)%jonc) then
            binome  = seg(ii)%Ijonc
            if (kkdebug) write (379,fmt='(x,"Delete segment and binome JONC status",x,I9,": binome",x,I9,":" )') ii,binome
            call stopjonc(ii,411)
            call stopjonc(binome,412)
            condition = .true.
          endif
          ii=i1etemp
      enddo

      if (condition) call connecij(I1o,Ivnne,413)

    case (IDEUX)

      ! E-side cleaning
      ii=I1o
      Ivnno = seg(I1o)%vnno

      counter = IZERO
      do while (ii /= Ivnno)

          if (counter > 1000) then
            call seginfo(i1o,"Infinite do-while loop surf correct cleaning i1o case IDEUX")
            write(379,*) "Check that seg(I1o)%vnno is correct : ", Ivnno, ":"
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voiso ",seg(ii)%voiso,":"
            stop "Infinite do-while loop surf correct cleaning i1o case IDEUX"
          endif

          counter = counter + IUN

          i1otemp=seg(ii)%voiso

          if (seg(ii)%gd > izero) then
            if (kkdebug) write (379,fmt='(x,"Delete segment GD status",x,I9,":" )') ii
            seg(ii)%gd  = izero
            condition = .true.
          endif

          if (seg(ii)%jonc) then
            binome  = seg(ii)%Ijonc
            if (kkdebug) write (379,fmt='(x,"Delete segment and binome JONC status",x,I9,": binome",x,I9,":" )') ii,binome
            call stopjonc(ii,421)
            call stopjonc(binome,422)
            condition = .true.
          endif
          ii=i1otemp
      enddo

      if (condition) call connecij(Ivnno,I1o,423)

    case DEFAULT

      print *,'surface segments with type < 4 must not appear here (surface_correct)'
      print *,I1o,seg(I1o)%surface
      stop

  end select

endif

!Segment I1e
if (seg(I1e)%surface > IZERO .and. .not. out(I1e)) then

  SURFTYPE = seg(I1e)%surface
  Condition = .false.
  if (kkdebug) write (379,fmt='(x,"Cleaning surface_correct I1e-side ",x,I9,":" )') I1e

  select case (SURFTYPE)

    case (IUN)

      ! E-side cleaning
      ii=I1e
      Ivnne = seg(I1e)%vnne

      counter = IZERO
      do while (ii /= Ivnne)

          if (counter > 1000) then
            call seginfo(i1e,"Infinite do-while loop surf correct cleaning i1e case IUN")
            write(379,*) "Check that seg(I1e)%vnne is correct : ", Ivnne, ":"
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
            stop "Infinite do-while loop surf correct cleaning i1e case IUN"
          endif

          counter = counter + IUN

          i1etemp=seg(ii)%voise

          if (seg(ii)%gd > izero) then
            if (kkdebug) write (379,fmt='(x,"Delete segment GD status",x,I9,":" )') ii
            seg(ii)%gd  = izero
            condition = .true.
          endif

          if (seg(ii)%jonc) then
            binome  = seg(ii)%Ijonc
            if (kkdebug) write (379,fmt='(x,"Delete segment and binome JONC status",x,I9,": binome",x,I9,":" )') ii,binome
            call stopjonc(ii,511)
            call stopjonc(binome,512)
            condition = .true.
          endif
          ii=i1etemp
      enddo

      if (condition) call connecij(I1e,Ivnne,513)

    case (IDEUX)! The case should not be observed
      ! E-side cleaning
      ii=I1e
      Ivnno = seg(I1e)%vnno

      counter = IZERO
      do while (ii /= Ivnno)

          if (counter > 1000) then
            call seginfo(i1e,"Infinite do-while loop surf correct cleaning i1e case IDEUX")
            write(379,*) "Check that seg(I1e)%vnno is correct : ", Ivnno, ":"
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voiso ",seg(ii)%voiso,":"
            stop "Infinite do-while loop surf correct cleaning i1e case IDEUX"
          endif

          counter = counter + IUN

          i1otemp=seg(ii)%voiso

          if (seg(ii)%gd > izero) then
            if (kkdebug) write (379,fmt='(x,"Delete segment GD status",x,I9,":" )') ii
            seg(ii)%gd  = izero
            condition = .true.
          endif

          if (seg(ii)%jonc) then
            binome  = seg(ii)%Ijonc
            if (kkdebug) write (379,fmt='(x,"Delete segment and binome JONC status",x,I9,": binome",x,I9,":" )') ii,binome
            call stopjonc(ii,521)
            call stopjonc(binome,522)
            condition = .true.
          endif
          ii=i1otemp
      enddo

      if (condition) call connecij(Ivnno,I1e,523)

    case DEFAULT

      print *,'surface segments with type < 4 must not appear here (surface_correct)'
      print *,I1e,seg(I1e)%surface
      stop

  end select

endif

if (kkdebug) then
  write(379,*) 'New segments are going to touch free surface',j
  write(379,*) 'segment i =',i,'surface type',seg(i)%surface
  write(379,*) 'segment i1o =',i1o,'surface type',seg(i1o)%surface
  write(379,*) 'segment i1e =',i1e,'surface type',seg(i1e)%surface
  call seginfo(i,'end surface correct segment i')
  call seginfo(i1o,'end surface correct segment i1o')
  call seginfo(i1e,'end surface correct i1e')
  !call check_connec('end surface_correct')
endif

end subroutine surface_correct

!> This function aims to correct the freesurfaces attribution
!> and to handle free surface connection in concave domains!!
subroutine correct_freeSurface_connection(i,i1o,i1e,lonad,lonvoad,lonvead,domchange,absdep)

implicit none

integer(kind=DPI)  :: I               !< Segment to Reconnect
integer(kind=DPI)  :: absdep          !< Absolute displacement of segment I
integer(kind=DPI)  :: VLI             !< Segment to Reconnect Line Vector
integer(kind=DPI)  :: VLIvois         !< Neighbor Segment Line Vector for segendO
integer(kind=DPI)  :: j               !< FS Counter
integer (kind=DPI) :: jj              !< Plane Counter
integer (kind=DPI) :: tmpvo           !< Temporary storage for neighbourg at Origin
integer (kind=DPI) :: tmpve           !< Temporary storage for neighbourg at End
integer (kind=DPI) :: ii              !< Temporary storage for a segment
integer(kind=DPI)  :: SurfType        !< Type of connection between the segment and the FS
integer(kind=DPI)  :: OldFreePlan     !< Storage of former free surface connected
integer(kind=DPI)  :: tmpVarFreePlan  !< Storage of temporary free surface selected
integer(kind=DPI)  :: TestTouch       !< TestTouch aims to know if veclin is going inwards the surface

integer            :: ref             !< Reference
integer(kind=DPI)  :: counter         !< Counter

integer(kind=DPI)  :: LONAD       !< Segment length after displacement
integer(kind=DPI)  :: LONVOAD     !< Origin Neighbor Segment length after displacement
integer(kind=DPI)  :: LONVEAD     !< End Neighbor Segment length after displacement
integer(kind=DPI)  :: i1o         !< First Neighbor at Origin rail
integer(kind=DPI)  :: i1e         !< First Neighbor at End rail
integer(kind=DPI)  :: voisi1o         !< First Neighbor at Origin
integer(kind=DPI)  :: voisi1e         !< First Neighbor at End
integer(kind=DPI)  :: i2o         !< Second Neighbor at Origin
integer(kind=DPI)  :: i2e         !< Second Neighbor at End
integer(kind=DPI)  :: modurx      !< Periodic boundary on X
integer(kind=DPI)  :: modury      !< Periodic boundary on X
integer(kind=DPI)  :: modurz      !< Periodic boundary on X
integer(kind=DPI)  :: domchange   !< change of domain during displacement (0=no change, 1=direct connection with new free surface,2=check each segment produced during displacement and then connect   )
integer(kind=DPI)  :: domcount    !< Number of domain segendO belongs


integer (kind=DPI), dimension(3) :: inormsurf !< Surface normal (Miller indices , integers)

real (kind=DP)     :: DotProdO    !<  Scalar Product between the vector [Oi-Intersect0] and the surface normal
real(kind=DP)      :: DistanceS   !<  FS Distance
real(kind=DP)      :: DistanceP   !<  Point Distance
real(kind=DP)      :: DistDiff    !< Store distanceS-distanceP value
real(kind=DP)      :: tmpDistDiff !< Store temporary distanceS-distanceP value
real(kind=DP)      :: planpos     !< Plane 'd' constant (ax+by+cz=d)

real (kind=DP), dimension(3)     :: normsurf   !< Surface normal (Miller indices , reals)
real (kind=DP), dimension(3)     :: IntersecO  !< Intersection coordinates with the surface if we extend the segment
real (kind=DP), dimension(3)     :: VectO      !< Vector [Oi-Intersect0]
real (kind=DP), dimension(3)     :: Oireal     !< real DP origin
real (kind=DP), dimension(3)     :: Oirealvois !< Neighbor real DP origin
real (kind=DP), dimension(3)     :: Eireal     !< real DP origin
real (kind=DP), dimension(3)     :: SegendO    !< real(true..) intersection with free surface OR point we want to consider to check domain

integer(kind=DP), dimension(3)  :: tempOi1e   !< temporary Origin i1e
integer(kind=DP), dimension(3)  :: tempEi1o   !< temporary End i1o
integer(kind=DPI),dimension(3)  :: Oi     !< Segment extremity coord
integer(kind=DPI),dimension(3)  :: Oivois !< Neighbor Segment extremity coord
integer(kind=DPI),dimension(3)  :: Ei     !< Segment extremity coord


integer(kind=DPI),dimension(2)  :: VLIchoice !< array / len, VLI / for shortest reconnection
integer(kind=DPI)               :: VLItest   !< VLI to test
integer(kind=DPI)               :: i1olen   !< VLI to test
integer(kind=DPI)               :: i1elen   !< VLI to test

logical :: Connection     !< True if segment touching extremity can connect the surface

logical,dimension(NbCvxDom) :: InsideDomains  !< True if segment touching extremity is inside the domain i

if (kkdebug) then
  write (379,fmt='("Function correct_freeSurface_connection ")')
  call seginfo(i,'Begining correct_freeSurface_connection')
  write (379,fmt='("i = ",I7)') i
endif

!Modulo operation to avoid problems
!------
modurx  = modur(1)
modury  = modur(2)
modurz  = modur(3)
seg(I)%o(1)=modulo(seg(I)%o(1),modurx)
seg(I)%o(2)=modulo(seg(I)%o(2),modury)
seg(I)%o(3)=modulo(seg(I)%o(3),modurz)
seg(I1o)%o(1)=modulo(seg(I1o)%o(1),modurx)
seg(I1o)%o(2)=modulo(seg(I1o)%o(2),modury)
seg(I1o)%o(3)=modulo(seg(I1o)%o(3),modurz)
seg(I1e)%o(1)=modulo(seg(I1e)%o(1),modurx)
seg(I1e)%o(2)=modulo(seg(I1e)%o(2),modury)
seg(I1e)%o(3)=modulo(seg(I1e)%o(3),modurz)
i1olen = seg(i1o)%norme
i1elen = seg(i1e)%norme
!------


Surftype        = seg(i)%surface
OldFreePlan     = seg(i)%Varfreeplan
tmpvarfreePlan  = OldFreePlan
VLI             = seg(i)%veclin
VLIvois         = HUGE(VLI)
VLIchoice(1)    = HUGE(VLI)
VLIchoice(2)    = IUN

!Domchange == ITROIS or IQUATRE : Segment touching the surface whose neigbor is also touching the surface and become of zero length

if (domchange <= IUN .or. domchange == ITROIS .or. domchange == ICINQ) then

  if (domchange <= IUN) then

    Oi          = seg(i)%O
    OiReal      = real(Oi,DP)
    normsurf(:) = Plane_MillerR(1:3,OldFreePlan)
    SegendO(:)  = InterPlanSeg(normsurf(:),Plane_pos(OldFreePlan),OiReal,Bveclin(:,VLI))
    if (domchange == IUN) then
      if (surftype == IUN   .and. (lonvead == IZERO .or. i1elen == IZERO)) call voisinage(i1e,1021)
      if (surftype == IDEUX .and. (lonvoad == IZERO .or. i1olen == IZERO)) call voisinage(i1o,1021)
    endif

  elseif (domchange == ITROIS) then

    call resolve_freeSurface_hooking_up(i,i1o,i1e,lonad,lonvoad,lonvead,NSEGMAX,absdep,3131)

    select case (SURFTYPE)

      case (IUN)

        if (lonvead == IZERO .or. i1elen == IZERO) call voisinage(i1e, 3001)
        i2o = seg(i1o)%voiso

        singulier(i1o)        = .false.
        seg(i1o)%Varfreeplan  = IZERO
        seg(i1o)%surface      = IZERO
        seg(i1o)%norme        = lonvoad !lonvoad should be zero here
        seg(i1o)%voiso        = NSEGMAX
        seg(i1o)%vnno         = NSEGMAX
        seg(i1o)%dom          = Plane_dom(OldFreePlan)

        !Elimination of i2o
        out(i2o) = .True.
        seg(i2o)%VOISO        = i2o
        seg(i2o)%VOISE        = i2o

        seg(i)%vnno = NSEGMAX
        seg(i)%dom  = Plane_dom(OldFreePlan)

        if (seg(i1o)%voise /= i ) then

          ii = seg(i1o)%voise

          counter = IZERO
          do while (ii /= i)

              if (counter > 10) then
                call seginfo(i1o,"Infinite do-while loop correct freesurface connection domachange ITROIS case IUN")
                write(379) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
                stop "Infinite do-while loop correct freesurface connection domachange ITROIS case IUN"
              endif

              counter = counter + IUN

              tmpve = seg(ii)%voise

              seg(ii)%voiso = ii
              seg(ii)%voise = ii
              out(ii) = .true.

              ii = tmpve

          enddo


        seg(i1o)%voise = I
        seg(i)%voiso   = I1o

        endif

      case (IDEUX)

        if (lonvoad == IZERO .or. i1olen == IZERO) call voisinage(i1o,3002)
        i2e = seg(i1e)%voise

        !I1e become a surface segment kneecap
        singulier(i1e)        = .false.
        seg(i1e)%Varfreeplan  = IZERO
        seg(i1e)%surface      = IZERO
        seg(i1e)%norme        = lonvead !lonvead should be zero here
        seg(i1e)%voise        = NSEGMAX
        seg(i1e)%vnne         = NSEGMAX
        seg(i1e)%dom          = Plane_dom(OldFreePlan)

        !Elimination of i2e
        out(i2e)              = .True.
        seg(i2e)%VOISO        = i2e
        seg(i2e)%VOISE        = i2e

        seg(i)%vnne = NSEGMAX

        if (seg(i1e)%voiso /= i ) then

          ii = seg(i1e)%voiso

          counter = IZERO
          do while (ii /= i)

              if (counter > 10) then
                call seginfo(i1e,"Infinite do-while loop correct freesurface connection domachange ITROIS case IDEUX")
                write(379) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voiso,":"
                stop "Infinite do-while loop correct freesurface connection domachange ITROIS case IUN"
              endif

              counter = counter + IUN

              tmpvo = seg(ii)%voiso

              seg(ii)%voiso = ii
              seg(ii)%voise = ii
              out(ii)       = .true.

              ii = tmpvo

          enddo

        seg(i1e)%voiso = I
        seg(i)%voise   = I1e

        endif

      case DEFAULT

        print *,'surface segments with type < 4 must not appear here (surface_correct)'
        print *,I1e,seg(I1e)%surface
        stop

    end select

    if (kkdebug) then
      call seginfo(i,'Domchange ITROIS - Cleaning loop connection to the freesurface')
      write (379,fmt='("!>  ")')
      !write (379,fmt='("!> Domchange move from ITROIS to IUN ")')
      !write (379,fmt='("!>  ")')
    endif

    Oi          = seg(i)%O
    OiReal      = real(Oi,DP)
    normsurf(:) = Plane_MillerR(1:3,OldFreePlan)
    SegendO(:)  = InterPlanSeg(normsurf(:),Plane_pos(OldFreePlan),OiReal,Bveclin(:,VLI))

  elseif (domchange == ICINQ) then

    if (surftype > IUN) then

      VLIvois       = seg(i1o)%veclin
      Oivois        = seg(i1o)%O
      OiRealvois    = real(Oivois,DP)
      normsurf(:)   = Plane_MillerR(1:3,OldFreePlan)
      SegendO(:)    = InterPlanSeg(normsurf(:),Plane_pos(OldFreePlan),OiRealvois,Bveclin(:,VLIvois))

      seg(i1o)%surface      = surftype
      seg(i1o)%varfreeplan  = OldFreePlan
      seg(i1o)%vnne         = NSEGMAX

      seg(i)%surface        = IZERO
      seg(i)%varfreeplan    = IZERO
      seg(i)%dom            = Plane_dom(OldFreeplan)

      ! I coordinates move to i1o end
      seg(i)%O        = SEG(i1o)%O(1:3) + BVECLIN(1:3,VLIvois)*lonvoad

      seg(i1o)%norme  = lonvoad
      seg(i)%norme    = lonad
      seg(i1e)%norme  = lonvead

      !if (lonad == izero) call voisinage(i,5021)
      !if (lonvead == IZERO .or. i1elen == IZERO) call voisinage(i1e, 5022)

      !lonvoad = IZERO
      lonad   = IZERO
      lonvead = IZERO

      ii = IZERO

      !if we are not already attached to the freesurface, move i1e Origin
      DotprodO = dot_product( real(SEG(i1o)%O(1:3)+seg(i1o)%norme*BVECLIN(:,VLIvois),DP), normsurf(:) )

      tempEi1o(1:3) = SEG(i1o)%O(1:3)+seg(i1o)%norme*BVECLIN(:,VLIvois)

      do while (DotprodO - Plane_pos(OldFreeplan) < ZERO)
        ii=ii+IUN
        if (ii > 100) then
          print*,"pb of infinit do while loop 5251 in surface_correct",  &
                DotprodO - Plane_pos(OldFreeplan), "kk=",kk, 'segment=', i1o
          stop
        endif
        tempEi1o(1:3) = tempEi1o(1:3) + BVECLIN(:,VLIvois)
        DotprodO      = dot_product( real(tempEi1o(1:3),DP), normsurf(:) )
      enddo

      if (ii > IZERO) SEG(i1o)%norme = ii

      seg(i)%norme  = IZERO
      lonvoad       = seg(i1o)%norme
      seg(i)%voise  = NSEGMAX
      seg(i)%vnne   = NSEGMAX

      !Elimination of i1e and eventually kneecaps
      counter = IZERO
      ii = i1e
      do while (ii /= Nsegmax)

          if (counter > 10) then
            call seginfo(i,"Infinite do-while loop correct freesurface connection domchange ICINQ case IDEUX")
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
            stop "Infinite do-while loop correct freesurface connection domachange ICINQ case IDEUX"
          endif

          counter = counter + IUN

          tmpve = seg(ii)%voise

          seg(ii)%voiso = ii
          seg(ii)%voise = ii
          out(ii)       = .true.

          ii = tmpve

      enddo

      if (i1olen == IZERO .or. lonvoad == IZERO) call voisinage(i1o, 5022)

    else

      VLIvois      = seg(i1e)%veclin
      Oivois       = seg(i1e)%O
      OiRealvois   = real(Oivois,DP)
      normsurf(:)  = Plane_MillerR(1:3,OldFreePlan)
      SegendO(:)   = InterPlanSeg(normsurf(:),Plane_pos(OldFreePlan),OiRealvois,Bveclin(:,VLIvois))

      seg(i1e)%surface      = surftype
      seg(i1e)%varfreeplan  = OldFreePlan
      seg(i1e)%vnno         = NSEGMAX

      seg(i)%surface      = IZERO
      seg(i)%varfreeplan  = IZERO
      seg(i)%dom          = Plane_dom(OldFreeplan)

      ! I coordinates move to i1o end
      seg(i)%O              = SEG(i1o)%O(1:3) + BVECLIN(1:3,SEG(i1o)%VECLIN)*lonvoad

      seg(i1o)%norme        = lonvoad
      seg(i)%norme          = lonad
      seg(i1e)%norme        = lonvead

      seg(i1e)%O(:)   = seg(i)%O
      !if (lonad == izero) call voisinage(i,5023)
      !if (lonvoad == IZERO .or. i1olen == IZERO) call voisinage(i1o, 5024)

      lonvoad = IZERO
      lonad   = IZERO
      !lonvead = IZERO

      !if we are not already attached to the freesurface, move i1e Origin
      DotprodO = dot_product( real(SEG(i1e)%O(1:3),DP), normsurf(:) )

      ii = IZERO
      VLItest = SEG(i1e)%VECLIN

      tempOi1e(1:3) = SEG(i1e)%O(1:3)

      do while (DotprodO - Plane_pos(OldFreeplan) < ZERO)
        ii=ii+IUN
        if (ii > 100) then
          print*,"pb of infinit do while loop 5151 in surface_correct",  &
                DotprodO - Plane_pos(OldFreeplan), "kk=",kk, 'segment=', i1e
          stop
        endif
        tempOi1e(1:3) = tempOi1e(1:3) - BVECLIN(:,VLIvois)
        DotprodO      = dot_product( real(tempOi1e(1:3),DP), normsurf(:) )
      enddo

      if (ii > IZERO) then
        SEG(i1e)%O(1:3) = tempOi1e(1:3) + BVECLIN(:,VLIvois)
        SEG(i1e)%norme  = lonvead + ii
      endif
      ! resolve_freeSurface_hooking_up : Correct i1e length to oldfreesurface. I Origin should be at i1e Origin
      !call resolve_freeSurface_hooking_up(i1e,i,seg(i1e)%voise,lonvead,lonad,seg(seg(i1e)%voise)%norme,NSEGMAX,absdep,5151)

      seg(i)%O      = SEG(i1e)%O(1:3)

      seg(i)%norme  = IZERO
      lonvead       = seg(i1e)%norme
      seg(i)%voiso  = NSEGMAX
      seg(i)%vnno   = NSEGMAX


      !Elimination of i1o and eventually kneecaps
      counter = IZERO
      ii = i1o
      do while (ii /= Nsegmax)

          if (counter > 10) then
            call seginfo(i,"Infinite do-while loop correct freesurface connection domchange ICINQ case IDEUX")
            write(379,*) "Info current state : We were at ii= : ", ii, ": with voiso ",seg(ii)%voiso,":"
            stop "Infinite do-while loop correct freesurface connection domachange ICINQ case IDEUX"
          endif

          counter = counter + IUN

          tmpvo=seg(ii)%voiso

          seg(ii)%voiso = ii
          seg(ii)%voise = ii
          out(ii)       = .true.

          ii = tmpvo

      enddo

      if (i1elen == IZERO .or. lonvead == IZERO) call voisinage(i1e, 5024)

    endif

    Oi     = seg(i)%O
    OiReal = real(Oi,DP)

    if (kkdebug) then
      call seginfo(i,'Domchange ICINQ -   Acting on i1o or i1e')
      write (379,fmt='("!>  ")')
    endif

  endif

else

  Oi(1:3)       = SEG(i1o)%O(1:3) + BVECLIN(1:3,SEG(i1o)%VECLIN)*lonvoad
  Ei(1:3)       = Oi(1:3) + lonad * BVECLIN(1:3,VLI)


  OiReal = real(Oi(:),DP)
  Eireal = real(Ei(:),DP)

  if (surftype == IUN) then

    segendO(:)      = Eireal(:)
    SEG(i1e)%O(:)   = Ei(:)
    seg(i1e)%norme  = lonvead
    seg(i)%O        = Oi
    seg(i)%norme    = lonad
    lonvoad         = IZERO
    if (lonvead == IZERO .or. i1elen == IZERO) call voisinage(i1e,4201)

    if (domchange == IQUATRE) then

      if (seg(i)%voiso /= i1o) then
        !Elimination of kneecaps
        counter = IZERO
        ii = seg(i)%voiso
        do while (ii /= i1o)

            if (counter > 10) then
              call seginfo(i,"Infinite do-while loop correct freesurface connection domchange IQUATRE case IDEUX")
              write(379,*) "Info current state : We were at ii= : ", ii, ": with voiso ",seg(ii)%voiso,":"
              stop "Infinite do-while loop correct freesurface connection domachange IQUATRE case IDEUX"
            endif

            counter = counter + IUN

            tmpvo = seg(ii)%voiso

            seg(ii)%voiso = ii
            seg(ii)%voise = ii
            out(ii)       = .true.

            ii = tmpvo

        enddo
        seg(i)%voiso   = i1o
        seg(i1o)%voise = i
      endif

      i2o = seg(i1o)%voiso

      singulier(i1o)        = .false.
      seg(i1o)%Varfreeplan  = IZERO
      seg(i1o)%surface      = IZERO
      seg(i1o)%norme        = lonvoad !lonvoad should be zero here
      seg(i1o)%voiso        = NSEGMAX
      seg(i1o)%vnno         = NSEGMAX
      seg(i1o)%dom          = Plane_dom(OldFreePlan)

      !Elimination of i2e
      out(i2o)              = .True.
      seg(i2o)%VOISO        = i2o
      seg(i2o)%VOISE        = i2o

      seg(i)%vnno = NSEGMAX

      if (kkdebug) then
        call seginfo(i,'Domchange IQUATRE SURF IUN - Cleaning loop connection to the freesurface')
        write (379,fmt='("!>  ")')
      endif

    endif

  else

    segendO(:)      = Oireal(:)
    seg(i1e)%O(1:3) = Ei(1:3)
    seg(i)%O        = Oi
    seg(i)%norme    = lonad
    seg(i1o)%norme  = lonvoad
    lonvead         = IZERO
    if (lonvoad == IZERO .or. i1olen == IZERO) call voisinage(i1o,4203)

    if (domchange == IQUATRE) then

      if (seg(i)%voise /= i1e) then
        !Elimination of kneecaps
        counter = IZERO
        ii = seg(i)%voise
        do while (ii /= i1e)

            if (counter > 10) then
              call seginfo(i,"Infinite do-while loop correct freesurface connection domchange IQUATRE case IDEUX")
              write(379,*) "Info current state : We were at ii= : ", ii, ": with voise ",seg(ii)%voise,":"
              stop "Infinite do-while loop correct freesurface connection domachange IQUATRE case IDEUX"
            endif

            counter = counter + IUN

            tmpve = seg(ii)%voise

            seg(ii)%voiso = ii
            seg(ii)%voise = ii
            out(ii)       = .true.

            ii = tmpve

        enddo
        seg(i)%voise = i1e
        seg(i1e)%voiso = i
      endif

      i2e = seg(i1e)%voise

      singulier(i1e)        = .false.
      seg(i1e)%Varfreeplan  = IZERO
      seg(i1e)%surface      = IZERO
      seg(i1e)%norme        = lonvead !lonvead should be zero here
      seg(i1e)%voise        = NSEGMAX
      seg(i1e)%vnne         = NSEGMAX
      seg(i1e)%dom          = Plane_dom(OldFreePlan)

      !Elimination of i2e
      out(i2e)              = .True.
      seg(i2e)%VOISO        = i2e
      seg(i2e)%VOISE        = i2e

      seg(i)%vnne = NSEGMAX

      if (kkdebug) then
        call seginfo(i,'Domchange IQUATRE SURF IDEUX - Cleaning loop connection to the fresurface')
        write (379,fmt='("!>  ")')
      endif

    endif

  endif
endif


! After displacement the voiso and voise of i can be of zero length and we must reconstruct
! the connectivity info. Better redefining i1o and i1e since elimination of pivot segment
! loops may have modified such information
if (surftype == IDEUX) then
  voisi1o = seg(i)%voiso
  counter = IZERO
  do while (seg(voisi1o)%norme == izero .and. voisi1o /= I1o)
    counter=counter + 1
    if (counter > ISEPT) then
      write(379,*)  i1o , i ,i1e
      call seginfo(i,'the counter voiso')
      stop
    endif
    if (voisi1o /= i1o) seg(voisi1o)%O = seg(i)%O
    call voisinage(voisi1o,541)
    voisi1o = seg(voisi1o)%voiso
  enddo
endif

if (surftype == IUN) then
  voisi1e = seg(i)%voise
  counter=IZERO
  do while (seg(voisi1e)%norme == izero .and. voisi1E /=i1e)
    counter=counter + 1
    if (counter > ISEPT) then
      write(379,*)  i1o , i ,i1e
      call seginfo(i,'the counter voiso')
      stop
    endif
    if (voisi1e /= i1e) seg(voisi1e)%O = seg(i1e)%O
    call voisinage(voisi1e,542)
    voisi1e=seg(voisi1e)%voise
  enddo
endif





if (domchange > IZERO) then
   InsideDomains(:)=.True. !Segment extremity is assumed by default inside the volume
else
   InsideDomains(:) = .false.
   InsideDomains(Plane_dom(OldFreePlan))=.true.
endif

DistanceS = huge(distanceS)
DistanceP = huge(distanceS)

if (domchange > IZERO) then
  DistDiff  = -huge(distanceS)
  tmpDistDiff= IUN
else
  DistDiff  = -huge(distanceS)
  tmpDistDiff= -IUN
endif

!--------------------------------------------------------------------
!-- Is the touching extremity inside the Volume ?
if (kkdebug) then
    write (379,fmt='(">- Is the touching extremity inside the Volume ? -< ")')
    write (379,fmt='("SegendO : " , 3F20.3)') SegendO
endif

if (domchange > IZERO) then

  Do jj= 1,NbPlanMax ! Loop over surfaces
    If ( .Not. InsideDomains(Plane_dom(jj)) ) cycle !If we already know that segment extremity is outside the domain of the plane, cycle

    normsurf(:) = Plane_MillerR(1:3,jj)
    DotProdO=SegendO(1)*normsurf(1)+SegendO(2)*normsurf(2)+SegendO(3)*normsurf(3)-Plane_pos(jj)

    if (kkdebug) write (379,fmt='("VFP",x,I2,x,"dom",x,I2,x," DotProdO",F15.3)') &
                                  & jj,Plane_dom(jj), DotProdO

    if (DotProdO > numtol_dotp) then
      InsideDomains(Plane_dom(jj)) =.False.
    endif

  Enddo

endif !domchange > IZERO

if (kkdebug) write(379,fmt='("# Inside Domains  ", 50L3)') InsideDomains(:)

If ( any(InsideDomains)) then

  domcount = count(InsideDomains(:))
  !--------------------------------------------------------------------
  !-- Seeking the closest freesurface
  if (kkdebug) write (379,fmt='(">- Segment is Inside - Seeking the closest free surface -< ")')

  if (kkdebug) then
    write (379,fmt='("#-")')
    write (379,fmt='("Old Free surface : ", I4," Equation : ", 3I20,F20.3)') OldFreePlan,Plane_MillerI(1:3,OldFreePlan),&
    Plane_pos(oldFreePlan)
    write (379,fmt='("Domain : ", I4," kind : ", I4 )') Plane_dom(OldFreeplan), Plane_knd(OldFreeplan)
  endif

  do J=NbPlanMax-NbFreePlan+1,NbPlanMax !Loop over freesurfaces

      if  (.not. InsideDomains(Plane_dom(J)) .or. &
          (j == Oldfreeplan .and. domchange /= ITROIS .and. domchange /= ICINQ )) cycle !We also want to test oldFreesurf for domchange == 3   as cleaning may replace the segment in the old domain

      normsurf(1:3)=normavec(Plane_MillerR(1:3,j))

      !Projected distance to an axe parallele to the surface going through the origin of the reference coordinate system
      if (normsurf(1) /= zero) then
        distanceS=Plane_pos(j)/Plane_MillerR(1,j)*normsurf(1)
      elseif (normsurf(2) /= zero) then
        distanceS=Plane_pos(j)/Plane_MillerR(2,j)*normsurf(2)
      elseif (normsurf(3) /= zero) then
        distanceS=Plane_pos(j)/Plane_MillerR(3,j)*normsurf(3)
      endif

      !Projected distance to an axe parallele to the surface going through the origin of the reference coordinate system
      distanceP = dot_product(normsurf,SegendO)

      tmpDistDiff = DistanceP - DistanceS

      if (kkdebug) then
        write (379,fmt='("#-")')
        write (379,fmt='("Free surface tested : ", I4," Equation : ", 3I20,F20.3)') J,Plane_MillerI(1:3,j),Plane_pos(j)
        write (379,fmt='("Domain : ", I4," kind : ", I4 )') Plane_dom(j), Plane_knd(j)
        write (379,fmt='("distanceS : ",F20.3," distanceP : ", F20.3," abs(DistanceS - DistanceP) : ", F20.3)') &
        & distanceS,distanceP,tmpDistDiff
      endif

      if (domchange > IZERO) then

        if ( tmpDistDiff < numtol_dotp .and. tmpDistDiff > DistDiff) then

          if (domcount > IUN) then !We are at a corner, we want to check segment direction
            !Select correct VLI for testing
            if (domchange /= ICINQ) then
              VLItest = VLI
            else
              VLItest = VLIvois
            endif
            !Test direction
            TestTouch   = Plane_MillerI(1,j)*BVECLIN(1,VLItest) &
                    & + Plane_MillerI(2,j)*BVECLIN(2,VLItest) &
                    & + Plane_MillerI(3,j)*BVECLIN(3,VLItest)

            if (kkdebug) write (379,fmt='("!> SegendO at corner, we check the segment direction = ", I4)') VlItest

            select case (SURFTYPE)

              case (IUN)

                if (TestTouch < IZERO) then
                  DistDiff = tmpDistDiff
                  tmpVarFreePlan = J
                  if (kkdebug) write (379,fmt='("===> LAST SELECTED SURFACE = ", I4)') tmpVarFreePlan
                endif

              case (IDEUX)

                if (TestTouch > IZERO) then
                  DistDiff = tmpDistDiff
                  tmpVarFreePlan = J
                  if (kkdebug) write (379,fmt='("===> LAST SELECTED SURFACE = ", I4)') tmpVarFreePlan
                endif

              case DEFAULT

                print *,'surface segments with type < 4 must not appear here (Surface selection) domchange /= IZERO'
                print *,'Surftype', SURFTYPE
                stop

            end select

          else
            DistDiff = tmpDistDiff
            tmpVarFreePlan = J
            if (kkdebug) write (379,fmt='("===> LAST SELECTED SURFACE = ", I4)') tmpVarFreePlan
          endif

        endif

      else

        if ( tmpDistDiff > -numtol_dotp .and. tmpDistDiff > DistDiff) then

          DistDiff = tmpDistDiff
          tmpVarFreePlan = J
          if (kkdebug) write (379,fmt='("===> LAST SELECTED SURFACE = ", I4)') tmpVarFreePlan

        endif

      endif

  enddo

  !--------------------------------------------------------------------
  !-- Attempt to reconnect free surface selected

  if (InsideDomains(Plane_dom(OldFreeplan)) .and. plane_dom(OldFreeplan) == plane_dom(tmpVarfreeplan)) then
  !simplest case: the domain is the same, looking for the new intersection with the closest free surface if tmpVarFreeplan /= OldFreeplan

    if (tmpVarFreeplan /= OldFreeplan) then
      if (kkdebug) then
        write (379,fmt='(">- Segment i is going to touch a new surface inside domain ", I4 )') &
              Plane_dom(OldFreeplan)
      endif
      seg(i)%Varfreeplan = tmpVarFreeplan
      call resolve_freeSurface_hooking_up(i,i1o,i1e,lonad,lonvoad,lonvead,NSEGMAX,absdep,6620)
    endif

    return

  else
    !# Attempt to extend segment I
    if (kkdebug) then
      write (379,fmt='(">- Attempt to reconnect free surface selected by extending segment I  -< ", I7 )') I
      Call seginfo(i,"Attempt to reconnect free surface selected by extending segment I")
    endif

    !VLItest=(/CONEC(VLI,1),CONEC(VLI,2), CONEC(CONEC(VLI,1),1), CONEC(CONEC(VLI,1),2),CONEC(CONEC(CONEC(VLI,1),1),1),CONEC(CONEC(CONEC(VLI,1),1),2)/)

    Connection = .FALSE.

    !Normal of the newly selected surface
    inormsurf=Plane_MillerI(1:3,tmpVarFreePlan)


    !TestTouch aims to know if veclin is going inwards the surface
    !In case of domchange ICINQ, We try to extend I which is the kneecap at fressurface
    TestTouch=inormsurf(1)*BVECLIN(1,VLI)+inormsurf(2)*BVECLIN(2,VLI)+inormsurf(3)*BVECLIN(3,VLI)

      if (kkdebug) then
        write (379,fmt='("#-")')
        write (379,fmt='("Free surface tested : ", I4)') tmpVarFreePlan
        write (379,fmt='("Equation : ", 3I20,F20.3)')inormsurf,Plane_pos(tmpVarFreePlan)
        write (379,fmt='("bveclin(:,VLI) : ",4I20 )') VLI, bveclin(:,VLI)
        write (379,fmt='("TestTouch : ",I20 )') TestTouch
        write (379,fmt='("Surftype : ",I3 )') Surftype
      endif

    If (domchange > IZERO) then

      If (domchange == ICINQ)   seg(I)%dom = Plane_dom(tmpVarFreePlan)

      If ( TestTouch < ZERO .and. Surftype==IUN) then

        normsurf= Plane_MillerR(1:3,tmpVarFreePlan)
        planpos = Plane_pos(tmpVarFreePlan)

        !Calculation of the intersection between Segment I and the Surface
        IntersecO(:) = InterPlanSeg(normsurf(:),planpos,Oireal,Bveclin(:,VLI))
        VectO(:) = Oireal(:)-IntersecO(:)

        if (dot_product(VectO,normsurf) > -numtol_dotp) then !Segment is already connected to the selected surface.
          !len = IZERO

          if (domchange /= ICINQ) then
            !Update segment
            seg(I)%VarFreePlan  = tmpVarFreePlan
            seg(I)%dom          = Plane_dom(tmpVarFreePlan)

            !Update I1o
            seg(I1o)%dom        = Plane_dom(tmpVarFreePlan)
            seg(I1o)%O(:)       = seg(I)%O(:)
          else
            !Update segment
            seg(I1e)%VarFreePlan  = tmpVarFreePlan
            seg(I1e)%dom          = Plane_dom(tmpVarFreePlan)
            seg(I)%dom            = Plane_dom(tmpVarFreePlan)

          endif

          if ( norvect(VectO(:)) >= (norivect(Bveclin(:,VLI))-numtol_dotp) ) then

            if (domchange == IDEUX .or. domchange == IQUATRE .or. domchange == ICINQ ) then
              if (norvect(VectO(:)) > norivect(lonad*Bveclin(:,VLI))) then
                if (domchange /= ICINQ) then
                  write (379,fmt='("VectO : ",3F9.3)') VectO
                  write (379,fmt='("norvect(VectO(:)) : ",F9.3 )') norvect(VectO(:))
                  write (379,fmt='("lonad : ",I3 )') lonad
                  Call seginfo(i,"Segment Origin and End are ouside freesurface - should not happend here - check debug")
                  stop "Segment Origin and End are ouside freesurface - should not happend here - check debug"
                endif
              else
                if (domchange == ICINQ) then
                  !Transfers of segment information from I1e to I
                  seg(I)%VarFreePlan  = tmpVarFreePlan
                  seg(I)%surface      = Surftype

                  seg(I1e)%VarFreePlan  = IZERO
                  seg(I1e)%surface      = IZERO
                  seg(I1e)%vnno         = I

                endif
                call resolve_freeSurface_hooking_up(i,i1o,i1e,lonad,lonvoad,lonvead,NSEGMAX,absdep,6630)
              endif
            else
              !Segment should be correctly connected to the surface, no more than one bveclin
              Call seginfo(i,"Segment is far from surface in correct_fs_connection with its Origin - check debug")
              stop "Segment is far from surface in correct_fs_connection with its Origin - check debug"
            endif
          endif

          !Connection to freesurface is accepted
          Connection = .TRUE.
          if (domchange == ICINQ) Call seginfo(i,"Final job domchange ICINQ SURF IUN")
        !else > go through connect_kneecap_fressurface to reconnect the segment by extending either the segment or one of its kneecap

        endif

      Elseif (TestTouch > ZERO .and. Surftype==IDEUX ) then

        normsurf= Plane_MillerR(1:3,tmpVarFreePlan)
        planpos = Plane_pos(tmpVarFreePlan)

        !Calculation of the intersection between Segment I and the Surface
        IntersecO(:) = InterPlanSeg(normsurf(:),planpos,Oireal,Bveclin(:,VLI))
        VectO(:) = Oireal+real(seg(i)%norme*Bveclin(:,VLI),DP)-IntersecO(:)

        if (dot_product(VectO,normsurf) > -numtol_dotp) then
          !len = IZERO

          !Update i1e
          if (domchange /= ICINQ) then
            !Update segment
            seg(I)%VarFreePlan  = tmpVarFreePlan

            seg(I1e)%dom       = Plane_dom(tmpVarFreePlan)
            seg(I1e)%O(:)      = Oi +  seg(I)%norme*Bveclin(:,VLI) ! update in case resolve freesurface hookinup not needed for domchange==2
          else
            !Update segment
            seg(I1o)%VarFreePlan  = tmpVarFreePlan
          endif

          if ( norvect(VectO(:)) >= (norivect(Bveclin(:,VLI))-numtol_dotp) ) then
            if (domchange == IDEUX .or. domchange == IQUATRE .or. domchange == ICINQ) then
              if (norvect(VectO(:)) > norivect(lonad*Bveclin(:,VLI))) then
                if (domchange /= ICINQ) then
                  write (379,fmt='("VectO : ",3F9.3 )') VectO
                  write (379,fmt='("norvect(VectO(:)) : ",F9.3 )') norvect(VectO(:))
                  write (379,fmt='("IntersecO(:) : ",3F15.3 )') IntersecO(:)
                  write (379,fmt='("Oireal+real(seg(i)%norme*Bveclin(:,VLI),DP) ",3F15.3 )') &
  Oireal+real(seg(i)%norme*Bveclin(:,VLI),DP)
                  write (379,fmt='("lonad : ",I3 )') lonad
                  Call seginfo(i,"Segment Origin and End are ouside freesurface - should not happend here - check debug")
                  stop "Segment Origin and End are ouside freesurface - should not happend here - check debug"
                endif
              else
                if (domchange == ICINQ) then
                  !Transfers of segment information from I1e to I
                  seg(I)%VarFreePlan  = tmpVarFreePlan
                  seg(I)%surface      = Surftype

                  seg(I1o)%VarFreePlan  = IZERO
                  seg(I1o)%surface      = IZERO
                  seg(I1o)%vnne         = I
                endif
                call resolve_freeSurface_hooking_up(i,i1o,i1e,lonad,lonvoad,lonvead,NSEGMAX,absdep,6640)
              endif
            else
              !Segment should be correctly connected to the surface, no more than one bveclin
              Call seginfo(i,"Segment is far from surface in correct_fs_connection with its End - check debug")
              stop "Segment is far from surface in correct_fs_connection with its End - check debug"
            endif
          endif

          !Connection to freesurface is accepted
          Connection = .TRUE.
          if (domchange == ICINQ) Call seginfo(i,"Final job domchange ICINQ SURF IDEUX")

        !else > go through connect_kneecap_fressurface to reconnect the segment by extending either the segment or one of its kneecap

        endif

      !Else : Connection is still False
      Endif

    endif !domchange

  endif !(InsideDomains(Plane_dom(OldFreeplan)%dom))


Else !If segment extremity is ouside all domains, then we have a problem

 write(379,*) "No domain found in correct_free_surface_connection!!!! kk =", kk
 write(379,*) " segment =", i
 write(379,*) " surftype =", surftype
 write(379,*) " Free surface touched = ", OldFreeplan
 write(379,*) " intersection =", SegendO

 stop "No domain found in correct_free_surface_connection!!!! "

endif !( any(InsideDomains) )


If ( .not. Connection ) then

  !# Deletion of all kneecaps at the extremity
  if (kkdebug) write (379,fmt='(">- Deletion of all kneecaps at the extremity  -< " )')

  if (Surftype == IUN) then

    if (domchange == ICINQ) then
      !Transfers of segment information from I1e to I
      seg(I)%VarFreePlan  = tmpVarFreePlan
      seg(I)%surface      = Surftype

      seg(I1e)%VarFreePlan  = IZERO
      seg(I1e)%surface      = IZERO
      seg(I1e)%vnno         = I

    else

      ii=seg(I)%voiso
      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Problem do-while loop correct_freeSurface_connection n1"
        counter = counter + IUN
        tmpvo = seg(ii)%voiso
        seg(ii)%voiso       = ii
        seg(ii)%voise       = ii
        seg(ii)%VarFreeplan = izero
        seg(ii)%surface     = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        if (kkdebug) write (379,fmt='(x,"Deleting kneecap",x,I9,":" )') ii
        ii = tmpvo
      enddo
      !Update neighborhood in order to call seginfo

    endif

    seg(I)%voiso = NSEGMAX
    seg(I)%vnno  = NSEGMAX

  else

    if (domchange == ICINQ) then
      !Transfers of segment information from I1e to I
      seg(I)%VarFreePlan  = tmpVarFreePlan
      seg(I)%surface      = Surftype

      seg(I1o)%VarFreePlan  = IZERO
      seg(I1o)%surface      = IZERO
      seg(I1o)%vnne         = I

    else

      ii=seg(I)%voise
      counter = IZERO
      do while (ii /= NSEGMAX)
        if (counter > 10) stop "Problem do-while loop correct_freeSurface_connection n2"
        counter = counter + IUN
        tmpve = seg(ii)%voise
        seg(ii)%voiso       = ii
        seg(ii)%voise       = ii
        seg(ii)%VarFreeplan = izero
        seg(ii)%surface     = izero
        seg(ii)%gd          = izero
        seg(ii)%dom         = izero
        out(ii) = .true.
        if (kkdebug) write (379,fmt='(x,"Deleting kneecap",x,I9,":" )') ii
        ii = tmpve
      enddo

    endif

    !Update neighborhood in order to call seginfo
    seg(I)%voise = NSEGMAX
    seg(I)%vnne  = NSEGMAX
  endif

  !-- Extension of segment - Test to reconnect the surface
  counter = IZERO
  do while (.not. Connection)
    !If (.not. Connection) then
    if (counter > 100) stop "Problem do-while loop correct_freeSurface_connection test to reconnect to the surface"

    !new procedure to find connection
    if(surftype == IUN) then
      Oi = Oi - counter * bveclin(:,VLI)
      seg(i)%O = Oi
    endif
    seg(i)%norme = seg(i)%norme + counter


    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to extend the segment attempt ", I3 , " -< ")') counter
    endif
    ref= 123
    VLItest = VLI
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    !-- First-First Kneecap - Test to reconnect the surface
    !If (.not. Connection) then
    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to connect the First-First Kneecap  -< ")')
    endif
    ref= 11
    VLItest = CONEC(VLI,1)
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    !-- First-Second Kneecap - Test to reconnect the surface
    !If (.not. Connection) then
    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to connect the First-Second Kneecap  -< ")')
    endif
    ref= 12
    VLItest = CONEC(VLI,2)
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    !-- Second-First Kneecap - Test to reconnect the surface
    !If (.not. Connection) then
    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to connect the Second-First Kneecap  -< ")')
    endif
    ref= 21
    VLItest = CONEC(CONEC(VLI,1),1)
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    !--Second-Second Kneecap - Test to reconnect the surface
    !If (.not. Connection) then
    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to connect the Second-Second Kneecap  -< ")')
    endif
    ref= 22
    VLItest = CONEC(CONEC(VLI,2),2)
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    !-- Third-First Kneecap - Test to reconnect the surface
    !If (.not. Connection) then
    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to connect the Third-First Kneecap  -< ")')
    endif
    ref= 31
    VLItest = CONEC(CONEC(CONEC(VLI,1),1),1)
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    !-- Third-Second Kneecap - Test to reconnect the surface
    !If (.not. Connection) then
    if (kkdebug) then
      write (379,*)
      write (379,fmt='(">- Trying to connect the Third-Second Kneecap  -< ")')
    endif
    ref= 32
    VLItest = CONEC(CONEC(CONEC(VLI,2),2),2)
    call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
    !Endif

    If (VLIchoice(1) /= huge(VLIchoice(1))) then
      Connection = .True.
      if (kkdebug) then
        write (379,*)
        write (379,fmt='(">- Final Connection with length ",I4," - Using VLI ",I4," at kk ",I9,"-< ")') VLIchoice(:),kk
        if (domchange == ICINQ) call seginfo(i,"begin domchange CINQ extended")
        write (379,*) I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref
      endif
      ref= 44
      VLItest = VLIchoice(2)
      call connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,Oi,Connection,VLIchoice,ref)
      lonad = seg(i)%norme
      if (surftype == 1) then
         call voisinage(i1e,9000)
      else
         call voisinage(i1o,9001)
      endif
      if (VLIchoice(1) > IDIX) &
         write (*,fmt='(">- Warning : Final Connection with length ",I4," - Using VLI ",I4," at kk ",I9,"-< ")') VLIchoice(:),kk
    endif
    counter= counter + IUN
  enddo

  if (counter > IUN) lonad = seg(i)%norme !Lonad is updated if its length had been modified to find a connection to the FS

  If (.not. Connection) then
    call seginfo(I,"! Could not connect the freesurface !")
    print *, "at kk = ", kk
    stop "! Could not connect the freesurface !"
  endif

endif !( .not. Connection )


if (kkdebug) call seginfo(I,"End of correct_freesurface_connection")

end subroutine correct_freeSurface_connection

subroutine connect_kneecap_freesurface(I,I1o,I1e,lonVOad,lonVEad,VLItest,SurfType,tmpVarFreePlan,O,Connection,VLIchoice,ref)

implicit none

integer(kind=DPI)  :: I               !< Segment to Reconnect
integer(kind=DPI)  :: VLItest         !< Line Vector tested
integer(kind=DPI)  :: SurfType        !< Type of connection between the segment and the FS
integer(kind=DPI)  :: tmpVarFreePlan  !< Storage of temporary free surface selected
integer(kind=DPI)  :: oldVarFreePlan  !< Storage of old free surface touched
integer            :: ref             !< Reference

integer(kind=DPI)  :: LONVOAD !< Origin Neighbor Segment length after displacement
integer(kind=DPI)  :: LONVEAD !< End Neighbor Segment length after displacement
integer(kind=DPI)  :: i1o     !< First Neighbor at Origin
integer(kind=DPI)  :: i1e     !< First Neighbor at End

integer(kind=DPI),dimension(3)  :: O !< Segment extremity coord
integer(kind=DPI),dimension(3)  :: E !< Segment extremity coord
integer(kind=DPI),dimension(3)  :: Otest !< Segment extremity coord
integer(kind=DPI),dimension(3)  :: Etest !< Segment extremity coord

integer(kind=DPI)  :: TestTouch       !< TestTouch aims to know if veclin is going inwards the surface
integer(kind=DPI)  :: TestDir         !< Testdir aims to know if veclin is going inwards the  oldsurface
integer(kind=DPI)  :: iNew,iNew2      !< New segment
integer(kind=DPI)  :: len             !< Multiplication factor of Veclin in order to reconnect the selected surface


!integer(kind=DPI)  :: lenold          !< Multiplication factor of Veclin in order to reconnect the selected surface

integer (kind=DPI), dimension(3) :: inormsurf !< Surface normal (Miller indices , integers)
integer (kind=DPI), dimension(3) :: inormoldsurf !< Surface normal (Miller indices , integers)

real(kind=DP)      :: NormO       !< Norm of VectO
real(kind=DP)      :: NormVeclin  !< Norm of considered Veclin
real(kind=DP)      :: planpos     !< Plane 'd' constant (ax+by+cz=d)
real(kind=DP)      :: DotProd     !< dot product result

real(kind=DP),dimension(3)      :: Oreal !< real DP Segment extremity coord
real(kind=DP),dimension(3)      :: Ereal !< real DP Segment extremity coord
!real(kind=DP)                   :: lenoldreal !< Length to the old Freesurface in the tested vectlin direction
real(kind=DP), dimension(3)     :: IntersecO  !< Intersection coordinates with the surface if we extend the segment
real(kind=DP), dimension(3)     :: VectO      !< Vector [Oi-Intersect0]

logical :: Connection     !< True if segment touching extremity can connect the surface

integer(kind=DPI),dimension(2)  :: VLIchoice !< array / len, VLI / for shortest reconnection

oldVarFreePlan = seg(i)%VarFreePlan

inormsurf=Plane_MillerI(1:3,tmpVarFreePlan)
inormoldsurf=Plane_MillerI(1:3,oldVarFreePlan)
!TestTouch aims to know if veclin is going inwards the surface
TestTouch=inormsurf(1)*BVECLIN(1,VLItest)+inormsurf(2)*BVECLIN(2,VLItest)+inormsurf(3)*BVECLIN(3,VLItest)
Testdir=inormoldsurf(1)*BVECLIN(1,VLItest)+inormoldsurf(2)*BVECLIN(2,VLItest)+inormoldsurf(3)*BVECLIN(3,VLItest)


if (kkdebug) then
  write (379,fmt='("Subroutine connect_kneecap_freesurface ", I7 )') ref
  write (379,fmt='(x,"VLItest ", I20 )') VLItest
  write (379,fmt='(x,"BVECLIN(:,VLItest) ", 3I20 )') BVECLIN(:,VLItest)
  write (379,fmt='(x,"TestTouch ", I20 )') TestTouch
  write (379,fmt='(x,"Testdir ", I20 )') Testdir
  write (379,fmt='(x,"SurfType ", I2 )') SurfType
endif

If ( TestTouch < IZERO .and. Surftype==IUN) then

  OReal(:)  = real(O(:),DP)

  planpos = Plane_pos(tmpVarFreePlan)

  !Calculation of the intersection between Segment I and the Surface
  IntersecO(:) = InterPlanSeg(Plane_MillerR(1:3,tmpVarFreePlan),planpos,Oreal(:),Bveclin(:,VLItest))
  VectO(:) = Oreal(:)-IntersecO(:)
  NormO = norvect(VectO)
  NormVeclin = norivect(Bveclin(:,VLItest))

  if (kkdebug) then
    write (379,fmt='(x,"IntersecO ", 3F20.3 )') IntersecO
    write (379,fmt='(x,"O ", 3I20 )') O
    write (379,fmt='(x,"VectO ", 3F20.3 )') VectO
    write (379,fmt='(x,"NormO ", F20.3 )') NormO
    write (379,fmt='(x,"NormVeclin ", F20.3 )') NormVeclin

  endif


  len = CEILING( NormO / NormVeclin - numtol_dotp) !length to add from the origin to the surface

  if (Testdir > IZERO) then
    planpos = Plane_pos(oldVarFreePlan)
    Otest(:) = O(:) - len * bveclin(:,VLItest)
    !Calculation of the intersection between Segment I and the old Surface
    Dotprod = Plane_MillerR(1,oldVarFreePlan)*Otest(1) + Plane_MillerR(2,oldVarFreePlan)*Otest(2) &
              + Plane_MillerR(3,oldVarFreePlan)*Otest(3) - planpos

    if (kkdebug) then
      write (379,fmt='(x,"O ", 3I20 )') Otest
      write (379,fmt='(x,"Dotprod ", F20.3 )') DotProd
    endif

  endif

  if (kkdebug) write (379,fmt='(x,"Length to add ", I20 )') len

  if ((Testdir <= IZERO) .or. &
      (Testdir > IZERO .and. Dotprod > numtol_dotp)) then


    if(real(len,DP) *  NormVeclin < real(VLIchoice(1),DP) * norivect(Bveclin(:,VLIchoice(2))))  then
      VLIchoice(1) = len
      VLIchoice(2) = VLItest
      !write (379,fmt='("At kk = ",I8,"!> Found a choice for reconnection : Length to add ", I5, " - VLItest ", I5 )') kk,len,VLItest
      !if (kkdebug) write (379,fmt='("!> Found a choice for reconnection : Length to add ", I5, " - VLItest ", I5 )') len,VLItest
    endif

    If (Connection) then
      !Segment creation and extension
      plusseg=plusseg+IUN
      INew=Nsegm+plusseg
      call init_seg(iNew)

      !Update I1o
      I1o = Inew

      if (VLItest == seg(i)%veclin) then

         seg(i)%norme = seg(i)%norme + len
         seg(i)%O(:) = O(:) - len * Bveclin(:,VLItest)
         seg(iNew)%O(:) = seg(i)%O(:)
         seg(i)%voiso = iNew
         seg(i)%vnno = NSEGMAX
         seg(i)%zerotl = .true.
         seg(i)%VarFreeplan = tmpVarFreePlan
         seg(I)%dom  =  Plane_dom(tmpVarFreePlan)
         seg(iNew)%VOISO = NSEGMAX
         seg(iNew)%VNNO = NSEGMAX
         seg(iNew)%VOISE = I
         seg(iNew)%VNNE = I
         seg(iNew)%veclin = Conec(VLItest,1)
         seg(INew)%dom  =  Plane_dom(tmpVarFreePlan)
         !Update LonVOad
         lonvoad = IZERO


      else

        !Kneecap creation
        plusseg=plusseg+IUN
        INew2=Nsegm+plusseg
        call init_seg(iNew2)

        !Update I
        seg(I)%VOISO     = iNew
        seg(I)%VNNO      = iNew
        seg(I)%Surface   = IZERO
        seg(I)%VarfreePlan   = IZERO
        seg(I)%dom   = Plane_dom(tmpVarFreePlan)
        singulier(I) = .false.
        if (seg(i1e)%norme == IZERO) seg(i)%vnne = seg(i1e)%vnne

        !Update iNew
        seg(iNew)%O(:)          = O(:) - len * Bveclin(:,VLItest)
        seg(iNew)%norme         = len
        seg(iNew)%veclin        = VLItest
        seg(iNew)%VOISO         = iNew2
        seg(iNew)%VNNO          = NSEGMAX
        seg(iNew)%VOISE         = I
        seg(iNew)%VNNE          = I
        seg(iNew)%dom           = Plane_dom(tmpVarFreePlan)
        Seg(i)%dom              = seg(iNew)%dom
        seg(iNew)%VarFreePlan   = tmpVarFreePlan
        seg(iNew)%surface       = surftype
        singulier(iNew)         = .true.
        seg(iNew)%zerotl        = .true.

        !Update LonVOad
        lonvoad = seg(i1o)%norme

        !Update iNew2
        seg(iNew2)%O(:)      = seg(iNew)%O
        seg(iNew2)%veclin    = Conec(VLItest,1)
        seg(iNew2)%VOISE     = iNew
        seg(iNew2)%VNNE      = iNew
        seg(iNew2)%dom       = Plane_dom(tmpVarFreePlan)

        !Check connectivity and add required kneecaps
        call connecij(iNew,I,ref)
        call voisinage(i,7001)
      endif
    endif

    !Else : Connection is still False
  endif
Elseif (TestTouch > IZERO .and. Surftype==IDEUX ) then

  E = O(:) + seg(i)%norme * Bveclin(:,seg(i)%veclin)
  Ereal = real(E,DP)


  planpos = Plane_pos(tmpVarFreePlan)

  !Calculation of the intersection between Segment I and the Surface
  IntersecO(:) = InterPlanSeg(Plane_MillerR(1:3,tmpVarFreePlan),planpos,EReal(:),Bveclin(:,VLItest))
  VectO(:) = Ereal(:)-IntersecO(:)
  NormO = norvect(VectO)
  NormVeclin = norivect(Bveclin(:,VLItest))

  if (kkdebug) then
    write (379,fmt='(x,"IntersecO ", 3F20.3 )') IntersecO
    write (379,fmt='(x,"E ", 3I20 )') E
    write (379,fmt='(x,"VectO ", 3F20.3 )') VectO
    write (379,fmt='(x,"NormO ", F20.3 )') NormO
    write (379,fmt='(x,"NormVeclin ", F20.3 )') NormVeclin
  endif

  len = CEILING( NormO / NormVeclin - numtold) !length to add from the origin to the surface

  if (Testdir < IZERO) then
    planpos = Plane_pos(oldVarFreePlan)
    Etest(:) = E(:) + len * Bveclin(:,VLItest)
    Dotprod = Plane_MillerR(1,oldVarFreePlan)*Etest(1) + Plane_MillerR(2,oldVarFreePlan)*Etest(2) &
              + Plane_MillerR(3,oldVarFreePlan)*Etest(3) - planpos
    !Calculation of the intersection between Segment I and the old Surface
    if (kkdebug) then
      write (379,fmt='(x,"Etest ", 3I20 )') Etest
      write (379,fmt='(x,"Dotprod ", F20.3 )') Dotprod
    endif

  endif
  if (kkdebug) write (379,fmt='(x,"Length to add ", I20 )') len
  if ((Testdir >= IZERO) .or. &
      (Testdir < IZERO .and. Dotprod > numtol_dotp)) then

    if (real(len,DP) *  NormVeclin < real(VLIchoice(1),DP) * norivect(Bveclin(:,VLIchoice(2))) ) then
      VLIchoice(1) = len
      VLIchoice(2) = VLItest
      !write (379,fmt='("At kk = ",I8,"!> Found a choice for reconnection : Length to add ", I5, " - VLItest ", I5 )') kk,len,VLItest
      !if (kkdebug) write (379,fmt='("!> Found a choice for reconnection : Length to add ", I5, " - VLItest ", I5 )') len,VLItest
    endif


    If (Connection) then
      !Segment creation and extension
      plusseg=plusseg+IUN
      INew=Nsegm+plusseg
      call init_seg(iNew) !Creation

      !Update I1o
      I1e = INew

      if (Vlitest == seg(i)%veclin) then
         seg(i)%norme = seg(i)%norme + len
         seg(I)%VOISE = iNew
         seg(I)%VNNE = NSEGMAX
         seg(I)%zerotl  = .true.
         seg(i)%VarFreeplan = tmpVarFreePlan
         seg(iNew)%VOISO  = I
         seg(iNew)%VNNO  = I
         seg(iNew)%veclin  = Conec(VLItest,1)
         seg(iNew)%O(:) = E(:) + len * Bveclin(:,VLItest)
         seg(INew)%dom  =  Plane_dom(tmpVarFreePlan)
         lonvead = IZERO

      else

        !Kneecap creation
        plusseg=plusseg+IUN
        INew2=Nsegm+plusseg
        call init_seg(iNew2) !Creation

        !Update I
        seg(I)%VOISE     = iNew
        seg(I)%VNNE      = iNew
        seg(I)%Surface   = IZERO
        seg(I)%VarfreePlan   = IZERO
        singulier(I) = .false.
        if (seg(i1o)%norme == IZERO) seg(i)%vnno = seg(i1o)%vnno

        !Update iNew
        seg(iNew)%O(:)          = E(:)
        seg(iNew)%norme         = len
        seg(iNew)%veclin        = VLItest
        seg(iNew)%VOISO         = I
        seg(iNew)%VNNO          = I
        seg(iNew)%VOISE         = iNew2
        seg(iNew)%VNNE          = NSEGMAX
        seg(iNew)%dom           = Plane_dom(tmpVarFreePlan)
        seg(iNew)%VarFreePlan   = tmpVarFreePlan
        seg(iNew)%surface       = surftype
        singulier(iNew)         = .true.
        seg(iNew)%zerotl        = .true.

        !Update LonVEad
        lonvead = seg(i1e)%norme

        !Update iNew2
        seg(iNew2)%O(:)      = E(:) + len * Bveclin(:,VLItest)
        seg(iNew2)%veclin    = Conec(VLItest,1)
        seg(iNew2)%VOISO     = iNew
        seg(iNew2)%VNNO      = iNew
        seg(iNew2)%dom       = Plane_dom(tmpVarFreePlan)

        !Check connectivity and add required kneecaps
        call connecij(I,iNew,ref)
        call voisinage(i,7002)
      endif
    endif
  endif
Endif

end subroutine connect_kneecap_freesurface

end MODULE microstructure
