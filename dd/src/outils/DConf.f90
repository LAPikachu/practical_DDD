!########################################################
!#                                                      #
!# CONFIGURATION INITIALE AVEC DISTRIBUTION ANGULAIRE   #
!# DES POINTS D'ANCRAGE ALEATOIRE ET DISCRETISATION   	#
!#                                                      #
!#                           Ronan MADEC le 07//03/2001 #
!#                                                      #
!########################################################

program MyConfInit

use Ecuyer_random

!*** DECLARATIONS HISTORIQUE (MyConfInit)
implicit none  	
integer ::  id2,ids,idbs1,idbs2,tnl,ix,iy,iz,acu
integer ::  i,j,k,syst(12),inito,iv(3),nb(3),it(3),nbud,indb(3),oindb(3)
integer ::  iorig(3),ie(3),kk,nsegm,nmax(3),longv,longc,longm,pmil,longs,nbtrue,incr
real(kind=8) long,dlong,x,xl,volume,xnmax,rnmax(3),modunbud(3)
real(kind=8) random3,coefpmmixte,coefpmvis,coefpmcoin
real(kind=8) densfr,nbsfr,ldens,xlong,xdlong,nbboitparsfr,para(3)
real(kind=8) yadja(10,10,10),lpboite,nbpboite,aufinal
integer ::  nbtruels,stat(11),nbsfrparsys
real(kind=8) dok,tot
integer ::  cok,caracsys(2000),nbase

!*** NOUVELLES VARIABLES
integer :: Oi(3),Ei(3),OEi(3),iseg(7,50000)
integer :: Lv,Lc,NbV,NbC,LsV,LsC
integer :: inbo,nbo,alterno(100,2),ccoin,cvis,cvisabs,premcara,deuxcara,premlong,deuxlong
integer :: nbasered,bveclin(3,96),id,totosys,isys
real(kind=8) :: veriflong,verifangle,coefid(96)

!*** FORMAT
3 format(2x,10i7)

!*** INITIALISATIONS
nbasered=8
yadja(:,:,:)=0.
nbtrue=0

!*** LECTURE DE LA BASE DE VECTEURS
open(1,FILE='../../../out/BaseDeVect',STATUS='OLD')
read(1,*) nbase
do i=1,nbase,1
    read(1,*) bveclin(1:3,i)
    coefid(i)=(bveclin(1,i)**2+bveclin(2,i)**2+bveclin(3,i)**2)**(0.5)
enddo
close(1)

coefpmvis=coefid(1)
coefpmcoin=coefid(((nbasered/4)+1))
coefpmmixte=coefid(2)

!*** INITIALISATION DU GU de Nb aleatoire ***
write(*,*)'Initialisation du generateur aleatoire'
write(*,*)'Entrez un nombre entier positif :'
read(*,*) inito
do i=1,inito
    x=taus88()
enddo


!*** LECTURE DU FICHIER DE CONFIGURATION
write(*,*)'Entrez le parametre de Reseau en micron (25.98d-4) :'
read(*,*) xl
write(*,*)'Entrez la taille du cube en microns (15) :'
read(*,*) xnmax 
write(*,*)'Parallelepipede (1. 1. 1. => cube) ?'
read(*,*) para(1)
read(*,*) para(2)
read(*,*) para(3) 
write(*,*)'Nombre de boite par dimension pour uniformiser la densite (NbUd) ?'
read(*,*) nbud
nmax(1) = (int(xnmax/xl*para(1))/8)*8	!De maniere generale pour eviter tout pb de sous reseau
nmax(2) = (int(xnmax/xl*para(2))/8)*8  !on 'assure que l'on manipule toujours des multiples de
nmax(3) = (int(xnmax/xl*para(3))/8)*8  !huit
!*volume en microns cube
volume=(xnmax**3)*(para(1)*para(2)*para(3))
nbpboite=nbud**3
rnmax(:)=dble(nmax(:))
modunbud(:)=rnmax(:)/dble(nbud)  
write(*,*)'=> a b c :',nmax 
write(*,*)'=> ',volume,' volume en mic**3' 
write(*,*)'=> ',nbpboite,' boites pour optimiser la densite locale'
write(*,*)'Entrez la densité (*10^12) :'
read(*,*) densfr
ldens=1./(sqrt(densfr))
write(*,*)'=> ',ldens,' distance carac' 
write(*,*)'Entrez la longueur des sources (microns) :'
read(*,*) long
xlong=long/xl
write(*,*)'Entrez la longueur de discretisation (microns) :'
read(*,*) dlong
xdlong=dlong/xl
longv=(int(xlong/coefpmvis)/4)*4 !220 => *4	multiple de 8 toujours
longc=(int(xlong/coefpmcoin)/4)*4 !422 => *4
longm=(int(xlong/coefpmmixte)/2)*2 !440 => *2
nbsfr=densfr*volume/long

lpboite=((densfr*volume)/nbpboite)
write (*,*) "=> longueur de source par boite",lpboite," micron"
lpboite=lpboite/xl
write (*,*) "Optimisation de la densite des boites"
write (*,*) "(0.15 par ex pour 15%, 100 = pas d'optimisation) :"
read (*,*) dok

write (*,*) "systemes"
totosys=0
do isys=1,12,1
    read (*,*) syst(isys)
    totosys=totosys+syst(isys)
enddo
nbsfrparsys=nbsfr/totosys

write(*,*)'=> ',nbsfr,' sources ',nbsfrparsys,' sources par sys'

!************************************************************************************
!*** CALCUL...
!************************************************************************************

!*Initialisation
!****************
nsegm=0

!*Boucle sur les systemes
!*************************
be:	do i=1,12	

    if (syst(i).eq.0) CYCLE be

    write(*,*) 'Systeme :',i

!*Tentatives
!************
    bi:   DO j=1,nbsfrparsys*syst(i),1	

        acu=0

25      acu=acu+1

!*** origine
        do k=1,3
            oi(k)=(int(taus88()*rnmax(k))/8)*8 !modulo 8 => decomposable ds les 6 base v c c' du cfc
        enddo

!*** Orientation de la ligne de la ligne
        if (taus88().gt.0.5d0) then	
           ccoin=nbasered*(i-1)+(nbasered/4)+1
        else
           ccoin=nbasered*(i-1)+3*(nbasered/4)+1
        endif
        cvisabs=nbasered*(i-1)+1
        if (taus88().gt.0.5d0) then
           cvis=cvisabs
        else
           cvis=nbasered*(i-1)+(nbasered/2)+1
        endif

!*** Longueur et discretisation
        Lc=int(taus88()*longc)
        Lv=int((sqrt(xlong**2-(Lc*coefpmcoin)**2))/coefpmvis)

        NbV=(Lv*coefpmvis)/xdlong
        NbC=(Lc*coefpmcoin)/xdlong

        if (Nbv.gt.Nbc) then
           Nbv=Nbc+1
        endif
        if (Nbc.gt.Nbv) then
           Nbc=Nbv+1
        endif
        if (Nbc.eq.Nbv) then
           Nbv=Nbc+1
        endif
        nbo=nbv+nbc	 

        if (NbV.ne.0) then
           LsV=(int(Lv/NbV)/4)*4 ! => multiple de 8
        else
           LsV=0
        endif
        if (NbC.ne.0) then
           LsC=(int(Lc/NbC)/4)*4 ! => multiple de 8
        else
           LsC=0
        endif

        if (Nbv.gt.Nbc) then
           premcara=cvis
           premlong=LsV
           deuxcara=ccoin
           deuxlong=LsC
        else
           premcara=ccoin
           premlong=LsC
           deuxcara=cvis
           deuxlong=LsV		
        endif

        Lv=LsV*NbV	 
        Lc=LsC*NbC

        ei=oi+Lv*bveclin(1:3,cvis)+Lc*bveclin(1:3,ccoin)
        oei=ei-oi	 
        veriflong=(oei(1)**2+oei(2)**2+oei(3)**2)**0.5

        verifangle =  ( ( oei(1)*bveclin(1,cvisabs)+&
        oei(2)*bveclin(2,cvisabs)+&
        oei(3)*bveclin(3,cvisabs)    ) / (veriflong*coefpmvis) )

        verifangle = acos ( dsign(1.d0,verifangle)*min(1.D0,dabs(verifangle)) )

        oi(:)=modulo(oi(:),nmax(:))
        ei(:)=modulo(ei(:),nmax(:))
        iv=oi
        do inbo=1,nbo,1
            if (((inbo/2)*2).ne.inbo) then
               alterno(inbo,1) = premcara
               alterno(inbo,2) = premlong
            else
               alterno(inbo,1) = deuxcara
               alterno(inbo,2) = deuxlong
            endif
!******************** Optimisation de la densite *******************************
            id=alterno(inbo,1)
            pmil=alterno(inbo,2)
            nbtruels=0
            it(:)=modulo(iv(:),nmax(:))		!coordonnee avec CLP
            oindb(:)=int(it(:)/modunbud(:))+1	!numero boite
            do tnl=0,pmil-1,1
                it(:)=iv(:)+tnl*bveclin(:,id)	 
                it(:)=modulo(it(:),nmax(:))
                indb(:)=int(it(:)/modunbud(:))+1
                if (oindb(1).ne.indb(1).or.oindb(2).ne.indb(2).or.oindb(3).ne.indb(3).or.(tnl.eq.(pmil-1))) then
                   if((yadja(indb(1),indb(2),indb(3))+(nbtruels*coefid(id))).gt.(lpboite*(1+dok))) goto 25
                   nbtruels=0
                   oindb(:)=indb(:)
                endif
                nbtruels=nbtruels+1
            enddo
            iv(:)=modulo(iv(:)+pmil*bveclin(:,id),nmax(:))
!*******************************************************************************		
        enddo

        iv=oi	 
        do inbo=1,nbo,1

!******************** Optimisation de la densite *******************************
            id=alterno(inbo,1)
            pmil=alterno(inbo,2)
            nbtruels=0
            it(:)=modulo(iv(:),nmax(:))
            oindb(:)=int(it(:)/modunbud(:))+1
            do tnl=0,pmil-1,1
                it(:)=iv(:)+tnl*bveclin(:,id)	 
                it(:)=modulo(it(:),nmax(:))
                indb(:)=int(it(:)/modunbud(:))+1
                if (oindb(1).ne.indb(1).or.oindb(2).ne.indb(2).or.oindb(3).ne.indb(3).or.(tnl.eq.(pmil-1))) then
                   yadja(indb(1),indb(2),indb(3))=yadja(indb(1),indb(2),indb(3))+(nbtruels*coefid(id))
                   nbtruels=0
                   oindb(:)=indb(:)
                endif
                nbtruels=nbtruels+1
            enddo

            nsegm=nsegm+1
            do k=1,3
                iseg(k,nsegm)=iv(k)
            enddo
            iseg(4,nsegm)=pmiL
            iseg(5,nsegm)=id

            if (inbo.ne.1) then
               iseg(6,nsegm)=nsegm-1	
            else
               iseg(6,nsegm)=0	 
            endif

            iv(:)=modulo(iv(:)+pmil*bveclin(:,id),nmax(:))

            if (inbo.ne.nbo) then
               iseg(7,nsegm)=nsegm+1	
            else
               iseg(7,nsegm)=0	 
            endif
            nbtrue=0
            do ix=1,nbud,1
                do iy=1,nbud,1
                    do iz=1,nbud,1
                        if(yadja(ix,iy,iz).gt.0.) nbtrue=nbtrue+1
                    enddo
                enddo
            enddo
            write (*,*) nsegm+1,' : ',iseg(1:7,nsegm)
!*******************************************************************************		
        enddo

        write (*,*) i,j,acu,veriflong/xlong,verifangle,dble(nbtrue)/nbpboite*100,'%'

    enddo bi
enddo be

!*** ECRITURE DU FICHIERS DE SOURCES DE FRANCK-READ
open(2,file='FRKOUT_pm',status='unknown')
write(2,*) nsegm,nmax 
do i=1,nsegm
    write(2,3) (iseg(j,i),j=1,7)
enddo
!*** Rappel des parametres utilises en fin de fichier ***
write (2,*)'Parametre de reseau :',xl
write (2,*)'Taille du cube en microns :',xnmax
write (2,*)'Parallelepipede (Non : 0, Oui : PARA) : ',para
write (2,*)'Nombre de boite par dimension pour uniformiser la densite (NbUd) : ',nbud
write (2,*)'Densité (*10^12) : ',densfr
write (2,*)'Longueur des sources (microns) et nombre de sources: ',long,nbsfr,NBSFRPARSYS
write (2,*)'Tolerance : ',dok
!*** Donnees statistiques ***
stat(:)=0
nbtrue=0
tot=0.
do ix=1,nbud,1
    do iy=1,nbud,1
        do iz=1,nbud,1
            if(yadja(ix,iy,iz).ne.0.) nbtrue=nbtrue+1
            incr=(yadja(ix,iy,iz)/(lpboite*(1+dok)))*9+1
            if(incr.gt.11) then
               incr=11
               write (*,*) "saturation",ix,iy,iz
            endif
            stat(incr)=stat(incr)+1
            tot=tot+yadja(ix,iy,iz)
        enddo
    enddo
enddo
write (2,*) "taux de remplissage des boites : ",float(nbtrue)/nbpboite
write (2,*) "Histogramme de la repartition des ",nbpboite," boites en densite"
do ix=1,11,1
    write (2,*) ix,(ix-1)/9.*densfr*(1+dok),&
    float(stat(ix))/nbpboite*100,'%'	     
enddo
write (2,*) "Densite obtenue :",tot/volume*xl,' 10^-12'
close(2)

end program MyConfInit

