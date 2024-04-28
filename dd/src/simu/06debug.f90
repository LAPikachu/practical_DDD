!===================================================================================================
!========================    DEBUT    MODULE  "DEBUG"   ============================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module contains many procedures helpful to debug mM.
!! \todo  Change some comments into doxygen comments
!! \todo  Some comments in this module have to be more specific and/or translated in English
module DEBUG

use constantes
use VARbase
use VARGLOB
use BRICAMAT

implicit none

contains

!#########################################################
!# subroutine Conf: Display the properties of segment i  #
!#########################################################
subroutine Conf(i)

implicit none

integer(kind=DPI) :: i      !< The segment i to display
integer(kind=DPI) :: ibr    !< The veclin index of segment i

ibr = seg(i)%veclin

61 format("========= conf de:",I7," ; KK=",I7," ; Long=",I7," ; Jonc=",L2,"; Ij = ",I7,"; Tj = ",I2)
62 format(" VD : ",I2," n[",3(I3,","),"], b[",3(I3,","),"], L[",3(I4,","),"]")
63 format("O(",I10,",",I10,",",I10,"); Vo=",I7,"(",I2,") ;VNNo=",I7,"(",I2,")")
64 format("E(",I10,",",I10,",",I10,"); Ve=",I7,"(",I2,") ;VNNe=",I7,"(",I2,")")
65 format("GD :", I1, "  OUT :", L2, "  RESDEP =",F6.3,"   NSEGM = ",I6," =======================================")

write(379,61) i,kk,seg(i)%norme,seg(i)%jonc,seg(i)%ijonc,seg(I)%tjonc
write(379,62) ibr, bvecnor(1:3,ibr),bveclin(1:3,assoc(ibr,1)),bveclin(1:3,ibr)
write(379,63) seg(i)%O(1:3),seg(i)%voiso,seg(seg(i)%voiso)%veclin,seg(i)%vnno,    &
              seg(seg(i)%vnno)%veclin
write(379,64) seg(i)%O(1:3)+seg(i)%norme*bveclin(1:3,ibr),seg(i)%voise,           &
              seg(seg(i)%voise)%veclin,seg(i)%vnne,seg(seg(i)%vnne)%veclin
write(379,65) seg(I)%gd, out(i), seg(i)%resdep, nsegm

end subroutine Conf

!######################################################################
!# subroutine segsinfo : Liste basic information on the segments tab  #
!######################################################################
subroutine segsinfo(message)

implicit none

integer(DPI)      :: i,m,Ei(3)
character (len=*) :: message


print *, " "
print *,"DEBUT <<<<<<  Infos <<<<<<<<<<<<<<<< : ", message

61 format(I5,":",I2,"(",I7,I7,I7,")",I6,"(",I7,I7,I7,")",I4,":",1x,                   &
          I5,"/",I5,":",I5,"/",I5,":",L1,":",I1,":",                                  &
          I5,1x,I3,1x,I5,1x,F6.3,1x,I4,I4,I3,I3,L1)

write (*,'(//," SegsInfo : ",I4," ; There is:",I4," segments to move",/)') kk,Nsegm
print*, message

if (nsegm < 20000) then

  write (*,*)"seg veclin O(1:3) norme E(1:3) idep vo/vnno ve/vnne jonc gd ijonc wait iboite rsesdep dom grain surf VFP singulier"

   Do i = 1,nsegm + PlusSeg

    m = seg(i)%veclin
    Ei(1:3) = seg(i)%O(1:3) + seg(i)%norme*bveclin(1:3,m)

    if(sysinfo < 0 .or. syseg(m) == sysinfo) then
      write(*,61) i,m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),INT(seg(i)%resdep+idep(i)),       &
            seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%gd,       &
            seg(i)%ijonc,seg(i)%wait,IBoite(i),seg(i)%resdep,seg(i)%dom,seg(i)%grain,      &
            seg(i)%surface,seg(i)%VarFreePlan,singulier(i)
    endif

       if (i > 100) exit

   enddo

else

   print *, " Segment info not printed as NSEGM > 20000 !!"

endif

print *, "########################################################################################################### "

end subroutine segsinfo

!#############################################################################
!> subroutine seginfo : List the segments connected to a given segment segini#
!#############################################################################
subroutine seginfo(segini,message)

implicit none

integer(DPI)  :: m        !< Segment index storage
integer(DPI)  :: i        !< Segment index storage
integer(DPI)  :: dislo    !< Segment index storage
integer(DPI)  :: segini   !< Segment of interest
integer(DPI)  :: compteur !< Segment counter
integer(DPI)  :: lim      !< Maximum number of segments to display

integer(DPI),DIMENSION(3)  :: Ei    !< Segment Coordinates at End

character (len=*) :: message  !< Message to display when entering subroutine

logical :: boucle             !< True if segment belongs to a loop

write(379,*) " "
write(379,*) "<<<<< SegInfo   : ", message
write(379,fmt='(" <<<<< kk = ",I7,"     Nsegm =",I7,"     Seg = ",I7,":")') kk,Nsegm+plusseg,segini
write(379,*) "--------------------------------------------------------------------------------------------------------------"
61 format(1I7,A1,8I7,I5,4I7,L5,I3,I7,I5,I7,F7.3,2I5,I4,L6,F7.3,L3)

! Maximum number of segment to show
lim = 25
if (nsegm <= 100) lim = 100

! First loop : We check if segini is a pinned segment or a section of closed loop
dislo = seg(segini)%voiso
boucle = .false.

! Search of the last connected segment on the O side
compteur = 0
if (dislo == segini .and. dislo == seg(segini)%voise) then
  write(379,*) " This segment must be eliminate :", segini
  RETURN

elseif (dislo == nsegmax) then

   dislo = segini
   compteur = 1

else

  do while(seg(dislo)%voiso /= nsegmax .and. dislo /= segini)

    m = seg(dislo)%voiso

    if(seg(m)%voise /= dislo) then
      write(379,*) " In seginfo of i=", segini,' at compteur= ',compteur
      write(379,*) " Problem of connectivity between",m," and",dislo
      write(379,*) " VoisO and VoisE of: ",m,seg(m)%voiso,seg(m)%voise
      write(379,*) " VoisO and VoisE of: ",dislo,seg(dislo)%voiso,seg(dislo)%voise
      stop "see ../out/debug file !!!"
    endif

    dislo     = seg(dislo)%voiso
    compteur  = compteur + 1
    if (compteur > lim) exit

  enddo

  if (dislo == segini) boucle = .true.

endif

write (379,fmt='(A8,8A7,A5,4A7,A5,A3,A7,A5,2A7,A5,A5,A4,A6,A7,A4)') &
      "segment","veclin",'',"O(1:3)",'',"norme",'',"E(1:3)",'',"idep", &
      "vo","vnno","ve","vnne","jonc","gd","ijonc","wait","iboite","resdep", &
      'dom',"surf", "VFP", "sing.", "angvis", "out"

i = dislo
compteur = 1

if (boucle) then            ! segini is part of a closed loop

  Do while (i /= nsegmax)

    m = seg(i)%veclin
    Ei(1:3) = seg(i)%O(1:3) + seg(i)%norme*bveclin(1:3,m)

    write(379,61) i,":",m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),idep(i),                    &
        seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%gd,         &
        seg(i)%ijonc,seg(i)%wait,IBoite(i),seg(i)%resdep,seg(i)%dom,seg(i)%surface,      &
        seg(i)%VarFreePlan,singulier(i),seg(i)%anglevis,out(i)

    if (i == seg(i)%voise) then

      write(379,*) "-"
      write(379,*) "DisloInfo :  i = seg(i)%voise (1), has a problem"
      exit

    else

      i = seg(i)%voise
      compteur = compteur + 1

      if (i == segini) exit
      if (compteur > 3*lim) then
        write(379,*) "-"
        write(379,*) "Exit inf loop seginfo boucle",compteur,i,seg(i)%voise,seg(seg(i)%voise)%voise
        exit
      endif

    endif

  enddo

else              ! segini is part of a pinned segment or a very long loop

  Do while (i /= nsegmax)

    m = seg(i)%veclin
    Ei(1:3) = seg(i)%O(1:3) + seg(i)%norme*bveclin(1:3,m)

    write(379,61) i,":",m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),idep(i),                    &
          seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%gd,       &
          seg(i)%ijonc,seg(i)%wait,IBoite(i),seg(i)%resdep,seg(i)%dom,seg(i)%surface,    &
          seg(i)%VarFreePlan,singulier(i),seg(i)%anglevis,out(i)

    if (i == seg(i)%voise) then

      write(379,*) "-"
      write(379,*) "DisloInfo :  i = seg(i)%voise (1), has a problem"
      exit

    else

      i = seg(i)%voise
      compteur = compteur + 1

      if (compteur > 3*lim) then
        write(379,*) "-"
        write(379,*) "exit inf loop bis seginfo boucle",compteur,i,seg(i)%voise,seg(seg(i)%voise)%voise
        exit
      endif

    endif

  enddo

  write(379,*) "--------------------------------------------------------------------------------------------------------------"

endif

write(379,*) "######################################################### Fin SEGinfo>"

end subroutine seginfo


!##########################################################################
!# subroutine disloinfo : decrire dislocation par dislocation tout les    #
!# segs constituant la micro dans l'ordre de connectivite                 #
!##########################################################################
subroutine disloinfo(message)

implicit none

integer(kind=DPI)   :: m
integer(kind=DPI)   :: i
integer(kind=DPI)   :: Ei(3)
integer(kind=DPI)   :: dislo
integer(kind=DPI)   :: nj
integer(kind=DPI)   :: j
integer(kind=DPI)   :: ldislo
character (len=*)   :: message
logical,allocatable :: zipseg(:)
logical             :: affich

allocate(zipseg(nsegm + Plusseg))
zipseg(1:nsegm+PlusSeg) = .false.

write(379,*) " "
write(379,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<Start Disloc info<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "
write(379,*) " <<<<<<<<<<<<<<<< : ", message

61 format(I6,":",I2,"(",I7,I7,I7,")",I6,"(",I7,I7,I7,")",I4,":",1x,                   &
          I6,"/",I6,":",I6,"/",I6,":",L1,":",I1,":",                                  &
          I5,1x,I3,1x,I5,1x,F6.3,1x,I3,I3,I3,L2)

write(379,*) " "
write(379,*) "Iteration = ",KK, "   NSEGM + PlusSeg =",Nsegm,"+",PlusSeg,"=",NSEGM+PlusSeg
write(379,*) " "
! First loop : The dislocation is ended by a pinned segment
do dislo = 1, nsegm + PlusSeg

     ldislo = seg(dislo)%veclin

     if(seg(dislo)%norme == 0 .and. seg(dislo)%voiso == dislo .and. seg(dislo)%voise == dislo) then

        zipseg(dislo) = .true.

     elseif (seg(dislo)%voiso == nsegmax) then

        affich = .false.
        if (sysinfo <= 0 .or. syseg(ldislo) == sysinfo) affich = .true.
    if (affich) write(379,*)"seg veclin O(1:3) norme E(1:3) idep vo/vnno ve/vnne jonc gd ijonc&
                          & wait iboite rsesdep dom surf VFP singulier"
        i = dislo

        Do while (i /= nsegmax)

            m = seg(i)%veclin
            Ei(1:3) = seg(i)%O(1:3) + seg(i)%norme*bveclin(1:3,m)
            if (i > nsegm + PlusSeg) write(379,*) " erreur 0 d'indice I dans disloinfo."

            if (zipseg(i)) then
               write(379,*) " DisloInfo : double voisinage 1 "
               exit
            else
               zipseg(i) = .true.
        if (affich) write(379,61) i,m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),idep(i),                       &
                          seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%gd,      &
                          seg(i)%ijonc,seg(i)%wait,IBoite(i),seg(i)%resdep,seg(i)%dom,seg(i)%surface,   &
                          seg(i)%VarFreePlan,singulier(i)
               if(i == seg(i)%voise) then
          write(379,*) "DisloInfo :  i = seg(i)%voise (1), has a problem"
                  exit
               else
                  i = seg(i)%voise
               endif
            endif

        enddo

     endif

 enddo

! Second loop : Segments are part of a closed loop
write(379,*) " "
write(379,*)" Segments part of a closed loop      :::::::::::::::::::::::::::::"
write(379,*)"seg veclin O(1:3) norme E(1:3) idep vo/vnno ve/vnne jonc gd ijonc&
                          & wait iboite rsesdep dom surf VFP singulier"

do dislo = 1, nsegm + PlusSeg

  if (zipseg(dislo)) cycle
  i = dislo

  Do while (.not. zipseg(i) .and. i /= nsegmax)
    m = seg(i)%veclin
    Ei(1:3) = seg(i)%O(1:3) + seg(i)%norme*bveclin(1:3,m)

    write(379,61) i,m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),INT(seg(i)%resdep+idep(i)), &
      seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%gd,     &
      seg(i)%ijonc,seg(i)%wait,IBoite(i),seg(i)%resdep,seg(i)%dom,seg(i)%surface,  &
      seg(i)%VarFreePlan,singulier(i)

    zipseg(i) = .true.

    if(i == seg(i)%voise) then
      write(379,*) "DisloInfo :  i = seg(i)%voise (2), pas normal"
      call conf(i)
      cycle
    else
      if (seg(i)%voise > nsegm + PlusSeg) then
        write(379,*) " erreur 1 d'indice I dans disloinfo.", I
        write(379,61) i,m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),INT(seg(i)%resdep+idep(i)), &
        seg(i)%voiso,seg(i)%vnno,seg(i)%voise,seg(i)%vnne,seg(I)%JONC,seg(i)%gd,       &
        seg(i)%ijonc,seg(i)%wait,IBoite(i),seg(i)%resdep,seg(i)%dom,seg(i)%surface,    &
        seg(i)%VarFreePlan,singulier(i)
        stop "see ../out/debug file !!!"
      endif

      i = seg(i)%voise
      if (zipseg(i)) write(379,*) "  .............................................."

    endif

  enddo

enddo

nj = 0

do i = 1,nsegm + PlusSeg

  if (seg(i)%jonc .and. seg(i)%voiso /= i ) then

    j = seg(i)%Ijonc

    if(seg(j)%jonc) then
      nj = nj + 1
    else
      write(379,*) message
      write(379,*) kk,"kk, disloinfo assymetric junction between i=",i," and j :", j
      call conf(i)
      call conf(j)
    endif

  endif

enddo

write(379,*) "Information about segment nsegmax =",nsegmax
i = nsegmax
m = seg(i)%veclin
Ei(1:3) = seg(i)%O(1:3)
write (379,61) i,m,seg(i)%O(1:3),seg(i)%norme,Ei(1:3),idep(I),seg(i)%voiso,seg(i)%vnno, &
seg(i)%voise, seg(i)%vnne, seg(I)%jonc, seg(i)%gd, seg(I)%ijonc

write(379,*) " KK = ", KK, "    ", message,"     --  fin ] >>>>>>"
write(379,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Fin disloinfo>"

deallocate(zipseg)

end subroutine disloinfo

!##########################################################################
!# subroutine Confie: Affiche une portion de boucle                       #
!###########################################################14/04/99#######
subroutine Confie(i,e)

integer (DPI) :: I,E,B

b=i

do while(b.ne.e)
    call conf(b)
    b=seg(b)%voise
    if (b.eq.nsegmax.or.b.eq.i) then
       write (379,*) "kk,ConfIE boucle : sortie forcee,i,e"
       write (379,*) kk,"ConfIE boucle : sortie forcee",i,e
       exit
    endif
enddo
call conf(b)

end subroutine Confie


!###################################################################################
!# subroutine check_connec: Verification de la connectivite de tous les segs       #
!# Techniquement : il s'agit de tester la connectivite entre i et son voisin en e  #
!###################################################################################
subroutine check_connec(message)

implicit none

character(len=*), OPTIONAL      :: message
logical                         :: cond
integer(kind=DPI)               :: i,ve,vo,VLe,vnne,J,vnno,VLI, Ti,Te,nsegmtot
integer(kind=DPI),dimension(3)  :: Oj,Ei
integer(kind=DPI)               :: JS                !< Temporary variable to store segment varfreeplan
real(kind=DP)       :: IntersectionJ(3)  !< Intersection point coordinates of a normal projection to the free surface
real(kind=DP)       :: DotProdTest       !< tested dot product

nsegmtot=nsegm + plusseg

do i = 1 , nsegmtot
    ve = seg(i)%voise
    vo = seg(i)%voiso
    ! These segments does not need to be tested as they are going to be eliminated
    if(i == ve .and. i == vo) cycle

    vnne = seg(i)%vnne
    vnno = seg(i)%vnno
    VLi = seg(i)%veclin
    VLe = seg(ve)%veclin
    Ti = seg(i)%norme
    Te = seg(ve)%norme
    if(seg(I)%probadev > zero .and. tyseg(Vli) /= iun)then
      write(379,*) message
      write(379,*) " check_connec_erreur de CS non vis, I=",I
      call seginfo(I, " avant erreur          ")
      stop  "see ../out/debug file !!!"
    endif

    if (abs(seg(I)%resdep) > 2*modep_max) then
      write(379,*) "error :debug check_connec: depassement en resdep kk,i=", kk,i,seg(I)%resdep
      stop  "see ../out/debug file !!!"
    endif

    ! We check the dimensions of tables
    cond = ((ve > nsegmtot .and. ve /= nsegmax) .or.        &
            (vo > nsegmtot .and. vo /= nsegmax) .or.        &
            (vnne > nsegmtot .and. vnne /= nsegmax) .or.    &
            (vnno > nsegmtot .and. vnno /= nsegmax))
    if (cond) then
      write(379,*) message, " nsegmtot: ", nsegmtot
      call seginfo(i," le seg i a des voisins d indice > nsegm+plusseg")
      stop  "see ../out/debug file !!!"
    endif

    ! Integrity of cross-slip segments is tested
    if (seg(i)%gd > izero) then
       if (tyseg(VLi) /= IUN) then
          write(379,*) message
          call seginfo(i," check-connec : GD et non VIS")
          stop  "see ../out/debug file !!!"
       endif
       IF(SYSEG(SEG(VE)%veclin) == SYSEG(SEG(Vo)%veclin)) then
         if (seg(i)%surface /= izero) then
           ! Segment touching a surface and gd should be relaxed from gd status
           seg(i)%gd = izero
         else
          write(379,*) "Segment ", I,": is declare Cross-sliping (GD) but its neigbours belong to the same slip system"
          stop  "see ../out/debug file !!!"
       endif
    endif

       ! Pivot segment at surface and gd should be relaxed from gd status
       if (seg(seg(i)%voiso)%surface /= IZERO .and. seg(i)%voise == NSEGMAX) then
         if (kkdebug) then
           write(379,*) 'Check_connec : a segment pivot at surface is found as gd - 1'
           write(379,*) 'kk = ',kk, 'i = ',i
           call seginfo(i," check-connec : GD and surface pivot e")
         endif
         seg(i)%gd = izero
       elseif (seg(seg(i)%voise)%surface /= IZERO .and. seg(i)%voiso == NSEGMAX) then
         if (kkdebug) then
           write(379,*) 'Check_connec : a segment pivot at surface is found as gd - 1'
           write(379,*) 'kk = ',kk, 'i = ',i
           call seginfo(i," check-connec : GD and surface pivot o")
         endif
         seg(i)%gd = izero
       endif
    endif

    ! integrity of junctions is tested
    if (seg(i)%jonc) then
       J = seg(i)%Ijonc
       if(.not. seg(J)%jonc .or. seg(J)%ijonc /= i) then
          write(379,*) message
          write(379,*) " Check-connec : kk =",kk," Jonction assymetrique entre i=",i," et j :", J
          if (kkdebug) stop "see ../out/debug file !!!"
       endif
    endif

    ! Connectivity of segment touching surface is tested
    if (seg(i)%surface == 1 .and. seg(i)%vnno /= nsegmax) then
      write(379,*) message
      write(379,*) " Check-connec : kk =",kk,"surface segment i=",i,"losts its pivot at the o side"
      if (kkdebug) then
        call seginfo(i, " dans check_connec: lost pivot o")
        stop "see ../out/debug file !!!"
      endif
    endif
    if (seg(i)%surface == 2 .and. seg(i)%vnne /= nsegmax) then
      write(379,*) message
      write(379,*) " Check-connec : kk =",kk,"surface segment i=",i,"losts its pivot at the e side"
      if (kkdebug) then
        call seginfo(i, " dans check_connec: lost pivot e")
        stop "see ../out/debug file !!!"
      endif
    endif

    ! The segments length is tested
    if (Ti < Izero) then
       write(379,*) message
       call seginfo(i, " dans check_connec: norme < 0 ???")
       stop "see ../out/debug file !!!"
    endif

    if(ve == nsegmax) cycle                    ! Seg ancre

    ! Segments connected at both ends to a different slip system are tested
    if (syseg(VLi) /= syseg(VLe) .and. (seg(i)%gd + seg(ve)%gd) < iun) then
       write(379,*) message
       write(379,*) " check-connec:sys dif. et pas de gd entre ", i," et  ",ve, " kk = ",kk
       call seginfo(i," check-connec: sys dif. et pas de gd entre")
       stop "see ../out/debug file !!!"
    endif

! Verification du sens de la circulation de la ligne
! on ne fait pas de teste si sys /= , ou i ou j est rotule GD
    if(modulo(int(abs(VLe -VLi)) - 1,nbasered) /= izero .and. &
       modulo(int(abs(VLi -VLe)) - 7,nbasered) /= izero) then
       write(379,*) message
       write(379,*) " check-connec:mauvaise circulation de ligne de i =", i, " kk = ",kk
       call seginfo(i," check-connec:mauvaise circulation de la ligne")
       stop "see ../out/debug file !!!"
    endif

    ! Nature of vnne (Voisin Non Nul a l Extremite) is tested, they must be non zero
    if (vnne/=nsegmax .and. seg(vnne)%norme == izero .and.  &
        .not.seg(vnne)%jonc .and. seg(vnne)%GD < iun) then
       write(379,*) message
       write(379,*) " check_connec : le vnne de i = ",i, "est nul .... KK=", KK
       call seginfo(i, " check_connec: le vnne de i est nul            ")
       stop "see ../out/debug file !!!"
       if (kkdebug) stop "see ../out/debug file !!!"
    endif

    ! Segments to be eliminated, but with a problem
    if (vo == I) then
       if (ve /= I) then        ! vo et vi ne sont pas idem a i
          write(379,*) message
          write(379,*) "Info sur VO" ; call conf(vo)
          write(379,*) "Info sur I" ; call conf(I)
          call seginfo(i, "check_connec: problem elimination 1")
          stop "see ../out/debug file !!!"
       elseif (seg(I)%norme /= 0) then    ! La norme de i n est pas nule
          write(379,*) message
          write(379,*) "Info sur VO" ; call conf(vo)
          write(379,*) "Info sur I" ; call conf(I)
          call seginfo(i, "check_connec: problem elimination 2")
          stop "see ../out/debug file !!!"
       else
          cycle
       endif

    elseif (ve == I) then
       ! cas normal on verifie la qualite des VNNO et VNNE
       ! se sont des voisin non nul sauf si GD ou Jonc
       write(379,*) message
       write(379,*) "Info sur VO" ; call conf(vo)
       write(379,*) "Info sur I" ; call conf(I)
       call seginfo(i, "check_connec: problem elimination 3")
       stop "see ../out/debug file !!!"
    endif

! 1.1) Verification que les indices des voisins pour i et ie sont compatibles
    if (seg(ve)%voiso /= i) then
       write(379,*) message
       write(379,*) " ERREUR 1.1 : Check_connec : voisin en E de i /= voisin en O de ie "
       write(379,*) " i=",i,"  ------- et -------- ve=",ve
       call seginfo(i, "voisin en E de i /= voisin en O de ie" )
       stop "see ../out/debug file !!!"
    endif

! 2) verification que le point Ei est egal au point Oie
    Oj(1:3) = seg(ve)%O(1:3)
    if (VLi < 1) then
       write(379,*) message
       call seginfo(i, " dans debug                 ")
       stop " check_connec: veclin =0      ??????"
    endif
    if (VLi > nbase) then
       write(379,*) message
       call seginfo(i, " dans debug                 ")
       stop " check_connec: veclin > nbase      ??????"
    endif
    Ei(1:3) = SEG(i)%O(1:3)+seg(i)%norme*bveclin(1:3,seg(i)%veclin)
    J = int(abs(modulo(Oj(1)-Ei(1),modur(1))),DPI) +    &
        int(abs(modulo(Oj(2)-Ei(2),modur(2))),DPI) +    &
        int(abs(modulo(Oj(3)-Ei(3),modur(3))),DPI)

    if (J /= 0) then
       write(379,*) message
       write(379,*) " ERREUR 2 : Check_connec : debug : Ei(x,y,z) /= Oie(x,y,z)    "
       write(379,*) " i=",i,"  ------- et -------- ie=",ve
       call seginfo(i, " dans debug                 ")
       stop "see ../out/debug file !!!"
    endif

!------------------------------------------------
!Check if segment is attached to the free surface
!------------------------------------------------

  if (seg(i)%VARFREEPLAN/=0) then !We test only segments touching a free surface

    if (seg(i)%surface/=UN) then !segment touching with end

      JS=seg(i)%varfreeplan
      !normal projection to the free surface
      IntersectionJ(:) = InterPlanSeg(Plane_MillerR(1:3,js),Plane_pos(js),&
         real(modulo( (seg(i)%O + seg(i)%norme*BVECLIN(:,seg(i)%VECLIN) ),modur(:) ),DP),BVECLIN(:,seg(i)%VECLIN))
      DotProdTest=dot_product(Plane_MillerR(1:3,js),&
         real(modulo(seg(i)%O + seg(i)%norme*BVECLIN(:,seg(i)%VECLIN),modur(:)),DP) - IntersectionJ(:))

        if (DotProdTest < -numtol_dotp) then
          call seginfo(i,"Segment is attached inside the volume with its end")
          write(379,*) "Segment is attached inside the volume with its end"
          write(379,*) "seg(i)%varfreeplan",seg(i)%varfreeplan
          write(379,*) "Plane_MillerI",Plane_MillerI(1:3,js)
          write(379,*) "IntersectionJ(:)",IntersectionJ(:)
          write(379,*) "DotProdTest",DotProdTest
          stop  "see ../out/debug file !!!"
        endif

    else !Segment touching with origin

      Js=seg(i)%varfreeplan
      !normal projection to the free surface
      IntersectionJ(:) = InterPlanSeg(Plane_MillerR(1:3,js),Plane_pos(js), &
             real(modulo(seg(i)%O,modur(:)),DP),BVECLIN(:,seg(i)%VECLIN))
      DotProdTest=dot_product(Plane_MillerR(1:3,js),real(modulo(seg(i)%O,modur(:)),DP) - IntersectionJ(:))

        if (DotProdTest < -numtol_dotp) then
          call seginfo(i,"Segment is attached inside the volume with its origin")
          write(379,*) "Segment is attached inside the volume with its origin"
          write(379,*) "seg(i)%varfreeplan",seg(i)%varfreeplan
          write(379,*) "Plane_MillerI",Plane_MillerI(1:3,js)
          write(379,*) "IntersectionJ(:)",IntersectionJ(:)
          write(379,*) "seg(i)%O - IntersectionJ(:)",modulo(seg(i)%O,modur(:)) - IntersectionJ(:)
          write(379,*) "DotProdTest",DotProdTest
          stop "see ../out/debug file !!!"
        endif

    endif

  endif

enddo

! Tests on the nsegmax segment
if (seg(nsegmax)%vnno.ne.nsegmax.or.seg(nsegmax)%vnne.ne.nsegmax.or.&
seg(nsegmax)%voiso.ne.nsegmax.or.seg(nsegmax)%voise.ne.nsegmax.or.&
seg(nsegmax)%o(1).ne.0.or.seg(nsegmax)%o(2).ne. 0 .or. seg(nsegmax)%veclin /= izero .or. &
seg(nsegmax)%o(3).ne.0 .or. seg(nsegmax)%norme /= nsegmax) then
   write(379,*) " KK ", kk
   write(379,*) " erreur 6 : check_connec : debug : caracteristique du segment NSEGMAX ont change"
   write(379,*) message
   call conf(nsegmax)
   stop "see ../out/debug file !!!"
endif

end subroutine check_connec

!##############################################################################
!# subroutine checkit : Selection de criteres juges oportuns a un temps donne #
!################################################################### 28/01/00 #
subroutine mes(message)

implicit none

character(len=*), OPTIONAL :: message

write(379,*) " ---------------- pause - debbug ------------------------------"
write(379,*) " "
write(379,*) message
write(379,*) " "
write(379,*) " ----Appuyer sur une touche pour continuer----------------------"
write(379,*) " ---------------------------------------------------------------"

write(*,*) " ---------------- pause - debbug ------------------------------"
write(*,*) " "
write(*,*) message
write(*,*) " "
write(*,*) " ----Appuyer sur une touche pour continuer----------------------"
write(*,*) " ---------------------------------------------------------------"

read (*,*)

end subroutine mes

!#########################################################################
!> This subroutine prints the segments information of a segment type object
!> Numbering starts from 1 to n (where n has to be less or equal than the segment type object dimension)
!#########################################################################
subroutine segdata(SEGTMP, n)

implicit none

integer(DPI)  :: n                           !< First n segments
integer(DPI)  :: ii                          !< Counter
type(SEGMENT), dimension(n) :: SEGtmp                      !< Segment type object

99 format(1I7,A1,8I7,4I7,L5,I3,I7,I5,F7.3,2I5,I4)

write (379,fmt='(A8,8A7,4A7,A5,A3,A7,A5,A7,A5,A5,A4)') &
      "segment","veclin",'',"O(1:3)",'',"norme",'',"E(1:3)",'', &
      "vo","vnno","ve","vnne","jonc","gd","ijonc","wait","resdep", &
      'dom',"surf", "VFP"
do ii = 1,n
  write(379,99) ii,":",segtmp(ii)%veclin,segtmp(ii)%O(1:3),segtmp(ii)%norme, &
          segtmp(ii)%O(1:3) + segtmp(ii)%norme*bveclin(1:3,segtmp(ii)%veclin), &
          segtmp(ii)%voiso,segtmp(ii)%vnno,segtmp(ii)%voise,segtmp(ii)%vnne,segtmp(ii)%JONC,segtmp(ii)%gd,      &
          segtmp(ii)%ijonc,segtmp(ii)%wait,segtmp(ii)%resdep,segtmp(ii)%dom,segtmp(ii)%surface,   &
          segtmp(ii)%VarFreePlan
enddo

end subroutine segdata

!#########################################################################
!> This subroutine prints no regular Field values for visualization in vtk format
!> tabDim = Dimension of the field tab
!> XCooR,YCoor,ZCoor = Coordinates of the field points in real format
!> FieldAmp = The field amplitude in real format
!#########################################################################
subroutine PrintFieldVTK(tabDim,XCooR,YCoor,ZCoor,FieldAmp,file2print)

implicit none

integer(kind=DPI), intent(in)           :: TabDim                 !< Dimension of the field tab
integer(kind=DPI)                       :: NbPoints               !< Dimension of the field tab
real(kind=DP), dimension(TabDim)        :: XCooR,YCoor,ZCoor      !< Coordinates of the field points
real(kind=DP), dimension(TabDim)        :: FieldAmp,FieldAmpTmp   !< The field amplitude
character(len=*)                        :: file2print             !< The name and place of the file to print

integer(kind=DPI)                       :: i,iplus
real(kind=DPI)                          :: NormFact, FieldMin, FieldMax

open (unit=10, status='unknown',file=file2print)

write(10,'(A)') '# vtk DataFile Version 2.0'
write(10,'(A)') 'Unstructured Grid Field'
write(10,'(A)') 'ASCII'
write(10,'(A)') 'DATASET UNSTRUCTURED_GRID'

! We eliminate from the field data point with amplitude == zero
iplus = 0
do i = 1, TabDim
  if (FieldAmp(i) /= zero) then
    iplus = iplus + 1
    XCooR(iplus) = XCooR(i)
    YCoor(iplus) = YCoor(i)
    ZCoor(iplus) = ZCoor(i)
    FieldAmpTmp(iplus) = FieldAmp(i)
  endif
enddo
NbPoints = iplus  ! The new Tabdim with useless field points removed

! iff all data points was removed, we go back to the initial field
if (iplus == 0) then
  NbPoints = TabDim
else
  FieldAmp(1:NbPoints) = FieldAmpTmp(1:NbPoints)
endif

! The field coordinates
write(10,'(A,I10,A)') 'POINTS', NbPoints, '  float'
do i = 1, NbPoints
  write(10,*) XCooR(i),YCoor(i),ZCoor(i)
enddo

! The field amplitude in real
write(10,'(A,I10)') 'POINT_DATA', NbPoints
write(10,'(A)') 'SCALARS FieldData float'
write(10,'(A)') 'LOOKUP_TABLE default'
do i = 1, NbPoints
  write(10,*) FieldAmp(i)
enddo

! The field amplitude in renormalized in int format
FieldMin = minval(FieldAmp(1:NbPoints))
FieldMax = maxval(FieldAmp(1:NbPoints))
if ((FieldMax - FieldMin) /= zero) then
  NormFact = 100.0/(FieldMax - FieldMin)
else
  NormFact = un
endif


! The field amplitude in real
write(10,'(A)') 'FIELD FieldData 1'
write(10,'(A,I10,A)') 'FieldRange  1', NbPoints, '  int'
do i = 1, NbPoints
  write(10,*) int((FieldAmp(i) - FieldMin) * NormFact)   ! an int value between 0 and 100
enddo

close(10)

end subroutine PrintFieldVTK

end module DEBUG

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
