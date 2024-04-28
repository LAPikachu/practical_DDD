!*******************************************************
!*** Construction d un pavage d octaedre tronque *******
!*** a partir d un reseau periodique BCC         *******
!*** (1*1)ou(2*2) en utilisant la methode        *******
!***  Voronoi
!*** fichier de sortie: un tableau de dimension3 ******
!*** donnant pour chaque grain 14 plans et pour  ******
!*** chaque plan ses indices de Miller et sa     ******
!*** position
!*******************************************************


program main

implicit none

!*Declaration des constantes

Integer, parameter  ::DP=8, DPI=8
Real(kind=DP),Dimension(1:2,1:3), parameter :: &
      POSun=reshape((/0.D0,0.5D0,0.D0,  &
                   0.5D0,0.D0,0.5D0/),(/2,3/))
Real(kind=DP),Dimension(1:16,1:3),parameter ::&
      POSdeux=reshape((/0.D0,0.D0,1.D0,1.D0,0.D0,1.D0 &
                       ,0.D0,1.D0,0.5D0,0.5D0,1.5D0  &
                       ,0.5D0,1.5D0,1.5D0,1.5D0,0.5D0 &
                       ,0.D0,1.D0,0.D0,0.D0,1.D0,1.D0,0.D0 &
                       ,1.D0,0.5D0,0.5D0,0.5D0,1.5D0,1.5D0 &
                       ,1.5D0,0.5D0,1.5D0,1.D0,0.D0,0.D0,1.D0, &
                       1.D0,0.D0,0.D0,1.D0,0.5D0,1.5D0,0.5D0, &
                       1.5D0,1.5D0,0.5D0,1.5D0,0.5D0/),(/16,3/))


Integer(kind=DPI), dimension(1:14,1:3),parameter::&
            Miller=reshape((/1,1,0,0,0,0, &
                             1,1,1,1,-1,-1,&
                             -1,-1,0,0,1,1,&
                             0,0,1,1,-1,-1,&
                             1,1,-1,-1,0,0,&
                             0,0,1,1,1,1,-1,-1,&
                             -1,-1,1,1/),shape=(/14,3/))


!* Declaration des variables


Integer(kind=DPI)    ::DimSyst, Nbgrains,i,j,k,IPtintersec(3),Fact,numb
real(kind=DP)        ::Ptintersec(3),Coefplan(3),distplan(14),indiceplan

Character(len=5)              ::POS
Character(len=1)              ::Carac

Real(kind=DP),allocatable     :: CoordGrain(:,:)
Integer(kind=DPI),allocatable :: TabFinal(:,:,:)
real(kind=DP),allocatable     :: TabFinal2(:,:)
logical                       ::KEY
                       
!*lecture des donnes d entree
open(51,file="../in/b_poly1",status='OLD')

carac = 'x'
do while (carac /= "#")
    read (51,*) carac
enddo
!*Dimension du systeme celui ci sera forcement cubique (X=Y=Z)
read(51,*) numb ! ne sert pas ici
read(51,*) KEY
read(51,*) DimSyst

!* Choix entre deux systemes
!*  deux grains
!*  seize grains
read(51,*) Nbgrains

if(NbGrains==2) then
        Fact=1
else
        Fact=2
endif

!*deux tableaux de maille primitive possible
!* POSun ou POSdeux
read(51,*) POS

close(51)

allocate (CoordGrain(Nbgrains,3))
allocate (TabFinal(Nbgrains,14,5))
allocate (TabFinal2(Nbgrains,14))

if(POS=='POSun')then
!print*,'syst1'
!print*,posun(1,:)
!print*,posun(2,:)
CoordGrain(:,:)=POSun(:,:)
!print*,'coordgrain vaut '
!print*,Coordgrain(1,:)
!print*,Coordgrain(2,:);pause
else
CoordGrain(:,:)=POSdeux(:,:)
endif
!*distance au centre de chacun des 14 plans
DistPlan=(/0.5D0,-0.5D0, &
                 0.5D0,-0.5D0, &
                 0.5D0,-0.5D0, &
                 0.25D0,-0.25D0, &
                 0.25D0,-0.25D0, &
                 0.25D0,-0.25D0, &
                 0.25D0,-0.25D0/)
              
!* pour chaque centre de grains, on trace la normale a chacun des 14 plans le 
!* constituant
!* les coordonnees du point d intersection plan/normale permettront d identifier
!* chaque plan
do i=1,Nbgrains
   do j=1,14
       Ptintersec(:)=(CoordGrain(I,:)+Miller(J,:)*DistPlan(J))
!       print*,'grain',i,'plan',j
!       print*,'grain de coordonnees',CoordGrain(i,:)
!       print*,'indice du plan',miller(j,:)       
!       print*,Ptintersec 
       do k=1,3
          IPtintersec(k)=Int(Ptintersec(k)*real(DimSyst)/real(Fact))
!          print*,IPtintersec
          if(IPtintersec(k)/=DimSyst) then
!                  print*,'avant mod',IPtintersec
                  IPtintersec(k)=modulo(IPtintersec(k),DimSyst)
!                  print*,'apres mod',IPtintersec ;pause
          endif
       enddo
!          print*,IPtintersec
       do k=1,3
          if (Miller(J,k)==0) then
                  CoefPlan(k)= 0
          else
                  Coefplan(k)= 1/Miller(J,K)
          endif
       enddo
!       print*,Coefplan 
       IndicePlan=(CoefPlan(1)*IPtintersec(1)+CoefPlan(2)*IPtintersec(2)+CoefPlan(3)*IPtintersec(3))
!       write(*,*) 'indice du plan',j,Indiceplan ;pause

       TabFinal(i,j,3:5)=Miller(j,1:3)
       TabFinal2(i,j)=IndicePlan
       write(*,*) 'grain',i,'  plan',j,'d indice:',miller(j,1:3),'de position',IndicePlan
    enddo
enddo

open(52,file="../in/b_poly2",status='Replace')
open(53,file="../in/b_poly3",status='Replace')

115 format(x,3(x,I2),x,f18.6)
 write(52,*)'fichier d entree de la simulation pour le pavage'
 write(52,*)
 write(52,*)'  Miller    Position'
 write(52,*)'  h  k  l'
 write(52,*)'-----------------------'
do i=1,Nbgrains
  do j=1,14
    write(52,115) Tabfinal(i,j,3:5),Tabfinal2(i,j)
  enddo
enddo

114 format(x,3(x,f18.6))
write(53,*)'fichier d entree de la simulation pour l initialisation du pavage'
write(53,*)'Coordonnees du centre de chaque grain'
write(53,*)'     X                      Y                     Z'
 write(53,*)'--------------------------------------------------------------'
do i=1,Nbgrains
  write(53,114) CoordGrain(i,:)*DimSyst/Fact
enddo
close(53)
close(52)
deallocate(TabFinal)
deallocate(Coordgrain)
deallocate(TabFinal2)

end


                       

