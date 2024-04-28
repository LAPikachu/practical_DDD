!
!                                       #   #    #   #
!                                       ## ##    ## ##
!                                       # # #    # # #
!                                       #   #    #   #
!                                       #   #    #   #
!                                       #   #ICRO#   #EGAS
!__________________________________________________________________________________________
!    
! MM (in the place of MicroMegas) is program of DDD (Discrete Dislocation Dynamics)
! initially developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
! Copyright (C) 1993, 1998, 2000, 2002, 2003, 2004 Benoit Devincre, Ladislas Kubin, Marc Condat,
! Ronan Madec, Ghiath Monnet
!
! This file is part of MM.
!
! MM is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! any later version.
!
! MM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MM; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
!
! For information contact by Email ronan.madec@cea.fr 
!__________________________________________________________________________________________

!
! Ronan MADEC fin Janvier 2003 : Determine la frontière entre le plan
! et ses images par pavage avec les polygones de la boite d'un seul plan
!__________________________________________________________________________
!
! Signe des décalages : vérifier la cohérence avec micromegas si temps
! mais doit être sans concéquences

#include"constantes_precision.F90"
#include"bricamat.F90"

!****** VARIABLE *****
module variables
use bricamat
 integer(kind=dpiv),parameter :: dpiv_2=2,dpiv_3=3,dpiv_4=4
 logical :: fillit
 integer(kind=dpiv) :: h,k,l,ni,Lx,Ly,Lz,dxy,dxz,dyx,dyz,dzx,dzy,llx,lly,llz,bufachkl,hklmax,hklmin,hklminmax,n0,nn0,pgtrans
 integer(kind=dpiv) :: bni,poly_blx,poly_bly,poly_blz,bdxy,bdxz,bdyx,bdyz,bdzx,bdzy
 integer(kind=dpiv) :: nit(25000),nnit(25000),ptit(3,25000),it
 integer(kind=dpiv) :: nit2(25000),nnit2(25000),ptit2(3,25000),it2,iit2,couple(25000),couple2(25000)
 integer(kind=dpiv) :: i,j,ik,x,y,z,dx,dy,dz,pg1,pg2,nbit,enx,eny,enmax,planexp(-2500000:2500000,7,6)
 integer(kind=dpiv) :: poly_front(75000,2,5),buf_i,jk
 integer(kind=dpiv) :: tab_point(3,10),ttab_point(2,10),nb,ox,oy,ex,ey,sx,sy,nbpoly,xmin,ymin,xmax,ymax,oox,ooy
 real(kind=dpr) :: ax(3),ay(3),az(3),nax,nay,naz
 real(kind=spr) :: fx,fy
 real(kind=dpr) :: buf_pt_front(3),buf_dist,buf_buf_dist,la_plus_faible_dist,la_plus_grande_dist,dist_carac,dcd10
 real(kind=dpr) :: buf_pt_front2(3),buf_dist2,buf_buf_dist2,la_plus_faible_dist2,la_plus_grande_dist2,dist_carac2
 integer(kind=dpiv) :: pt_dlpf(5),pt_dlpg(5),buf_ik,pt_dlpf2(5),pt_dlpg2(5),nbitsat,demi
 integer(kind=dpi) :: v1(3),v2(3),v3(3),v4(3),vv1(3),vv2(3),vv3(3),vv4(3),Helvind(2),modpol,fachkl,fachkl_cris
 real(kind=dpr) :: nor1,nor2,nor3,nor4,nor5,bary(5),bary2(5)
end module variables

!****** MANIPULATION DE POLYGONES *****
!__________________________________________________________________________

module polygones

use bricamat

use variables

contains

!__________________________________________________________________________

function plan(xx,yy,zz,tdx,tdy,tdz) ! en entree un pt*dechkl en sortie l'effet des clp + numero plan

 integer(kind=dpiv) :: plan,xx,yy,zz,tdx,tdy,tdz
 integer(kind=dpi) :: chg_boite,xxx,yyy,zzz,pplan

 chg_boite=0
 xxx=xx
 yyy=yy
 zzz=zz
 tdx=0
 tdy=0
 tdz=0
 
! première vague de décalages
 
 if (xx.eq.llx) then
  chg_boite=chg_boite+1
  tdy=tdy-fachkl*dyx
  tdz=tdz-fachkl*dzx
 else
  if (xx.eq.0) then
   chg_boite=chg_boite+1
   tdy=tdy+fachkl*dyx
   tdz=tdz+fachkl*dzx
  endif
 endif

 if (yy.eq.lly) then
   chg_boite=chg_boite+1
   tdx=tdx-fachkl*dxy
   tdz=tdz-fachkl*dzy
 else
  if (yy.eq.0) then
   chg_boite=chg_boite+1
   tdx=tdx+fachkl*dxy
   tdz=tdz+fachkl*dzy
  endif
 endif

 if (zz.eq.llz) then
   chg_boite=chg_boite+1
   tdx=tdx-fachkl*dxz
   tdy=tdy-fachkl*dyz
 else
  if (zz.eq.0) then
   chg_boite=chg_boite+1
   tdx=tdx+fachkl*dxz
   tdy=tdy+fachkl*dyz
  endif
 endif

! 2ème

 xx=xx+tdx
 if(xx.ge.llx) then
   xx=xx-llx
 else
  if(xx.le.0) then
    xx=xx+llx
  endif
 endif

 yy=yy+tdy 
 if(yy.ge.lly) then
   yy=yy-lly
 else
  if(yy.le.0) then
    yy=yy+lly
  endif
 endif

 zz=zz+tdz 
 if(zz.ge.llz) then
   zz=zz-llz
 else
  if(zz.le.0) then
    zz=zz+llz
  endif
 endif

 if (xx.lt.0.or.xx.gt.llx.or.yy.lt.0.or.yy.gt.lly.or.zz.lt.0.or.zz.gt.llz.or.chg_boite.ne.1) then
  write (*,*) "P ",xxx,yyy,zzz
  write (*,*) "S ",tdx,tdy,tdz
  write (*,*) "-> ",xx,yy,zz
  write (*,*) "L ",llx,lly,llz
  write (*,*) "D ",xxx-xx,yyy-yy,zzz-zz
 endif  
 if (xx.lt.0.or.xx.gt.llx.or.yy.lt.0.or.yy.gt.lly.or.zz.lt.0.or.zz.gt.llz) then
  write (*,*) "Pb de CLP",chg_boite
  stop
 endif 
 if (chg_boite.eq.0) then
  write (*,*) "Pb de changement de boite",chg_boite
  stop
 endif
 if (chg_boite.eq.2) then
  write (*,*) "Changement de boite sur une arete",chg_boite
!  stop
!  pause
 endif
 if (chg_boite.eq.3) then
  write (*,*) "Changement de boite dans un coin ie ponctuel",chg_boite
  stop
 endif

 pplan=h*(xx)+k*(yy)+l*(zz)
 
 if(mod(pplan,fachkl).ne.0)then
  write (*,*) h,xx,k,yy,l,zz,"->",pplan,"/",fachkl
  write (*,*) "P ",xxx,yyy,zzz
  write (*,*) "S ",tdx,tdy,tdz
  write (*,*) "-> ",xx,yy,zz
  write (*,*) "L ",llx,lly,llz
  write (*,*) "D ",xxx-xx,yyy-yy,zzz-zz
  stop"Problème de calcul d'indice de plan"
 endif
 
 plan=pplan/fachkl

! VERIFICATION DES PLANS
! write (*,*) h,xx,k,yy,l,zz,"->",pplan,"/",fachkl
 
 tdx=(xxx-xx)/fachkl
 tdy=(yyy-yy)/fachkl
 tdz=(zzz-zz)/fachkl
 
 return

end function plan

!__________________________________________________________________________

subroutine polygone(n,nb_point) ! plan (entree), nb_point(sortie)
                                ! tab_point(:,ti) coordonnée 3D
integer(kind=dpiv) :: n,nb_point,tshift(3)

integer(kind=dpi) :: xx,yy,zz,perm,buf(3),tti,ti

! n=hx+ky+lz
! x=y=0 => z=n/l ...

nb_point=0

if(l.ne.0) then
 zz = n/l
 if (zz.ge.0.and.zz.le.Lz) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=0
  tab_point(2,nb_point)=0
  tab_point(3,nb_point)=zz
 endif
 zz = (n-h*Lx-k*Ly)/l
 if (zz.ge.0.and.zz.le.Lz) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=Lx
  tab_point(2,nb_point)=Ly
  tab_point(3,nb_point)=zz
 endif
 zz = (n-h*Lx)/l
 if (zz.ge.0.and.zz.le.Lz) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=Lx
  tab_point(2,nb_point)=0
  tab_point(3,nb_point)=zz
 endif
 zz = (n-k*Ly)/l
 if (zz.ge.0.and.zz.le.Lz) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=0
  tab_point(2,nb_point)=Ly
  tab_point(3,nb_point)=zz
 endif
endif

if (k.ne.0) then
 yy = n/k
 if (yy.ge.0.and.yy.le.Ly) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=0
  tab_point(2,nb_point)=yy
  tab_point(3,nb_point)=0
 endif
 yy = (n-h*Lx-l*Lz)/k
 if (yy.ge.0.and.yy.le.Ly) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=Lx
  tab_point(2,nb_point)=yy
  tab_point(3,nb_point)=Lz
 endif
 yy = (n-h*Lx)/k
 if (yy.ge.0.and.yy.le.Ly) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=Lx
  tab_point(2,nb_point)=yy
  tab_point(3,nb_point)=0
 endif
 yy = (n-l*Lz)/k
 if (yy.ge.0.and.yy.le.Ly) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=0
  tab_point(2,nb_point)=yy
  tab_point(3,nb_point)=Lz
 endif
endif

if(h.ne.0) then
 xx = n/h
 if (xx.ge.0.and.xx.le.Lx) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=xx
  tab_point(2,nb_point)=0
  tab_point(3,nb_point)=0
 endif
 xx = (n-k*Ly-l*Lz)/h
 if (xx.ge.0.and.xx.le.Lx) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=xx
  tab_point(2,nb_point)=Ly
  tab_point(3,nb_point)=Lz
 endif
 xx = (n-k*Ly)/h
 if (xx.ge.0.and.xx.le.Lx) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=xx
  tab_point(2,nb_point)=Ly
  tab_point(3,nb_point)=0
 endif
 xx = (n-l*Lz)/h
 if (xx.ge.0.and.xx.le.Lx) then
  nb_point=nb_point+1
  tab_point(1,nb_point)=xx
  tab_point(2,nb_point)=0
  tab_point(3,nb_point)=Lz
 endif
endif

tti=1
ti=1
do while (tti.le.nb_point)
do while (ti.le.nb_point)
  if(tti.ne.ti) then
   if(tab_point(1,ti).eq.tab_point(1,tti).and.&
     tab_point(2,ti).eq.tab_point(2,tti).and.&
     tab_point(3,ti).eq.tab_point(3,tti)) then
    tab_point(:,ti:(nb_point-1))=tab_point(:,(ti+1):nb_point)
    nb_point=nb_point-1
   else
    ti=ti+1
   endif
  else
    ti=ti+1  
  endif
 enddo
 ti=1
 tti=tti+1
enddo

perm=1
do while (perm.ne.0)
 perm=0
 do ti=1,nb_point,1
  tti=ti+1
  if(ti.eq.nb_point) tti=1
  
  if(h.ne.0.or.l.ne.0) then ! TOUS LES CAS SAUF 010
  if(.not.(((tab_point(1,ti).eq.tab_point(1,tti)).and.&
            (tab_point(1,ti).eq.0.or.tab_point(1,ti).eq.lx)).or.&
           ((tab_point(2,ti).eq.tab_point(2,tti)).and.&
            (tab_point(2,ti).eq.0.or.tab_point(2,ti).eq.ly)).or.&
           ((tab_point(3,ti).eq.tab_point(3,tti)).and.&
            (tab_point(3,ti).eq.0.or.tab_point(3,ti).eq.lz)))) then
   buf(:)=tab_point(:,nb_point)
   tab_point(:,(tti+1):nb_point)=tab_point(:,tti:(nb_point-1))
   tab_point(:,tti)=buf(:)
   perm=1
   exit
  endif
  else ! CAS 010 : y=cste pour tous les points d'un plan
  if(.not.(((tab_point(1,ti).eq.tab_point(1,tti)).and.&
            (tab_point(1,ti).eq.0.or.tab_point(1,ti).eq.lx)).or.&
           ((tab_point(3,ti).eq.tab_point(3,tti)).and.&
            (tab_point(3,ti).eq.0.or.tab_point(3,ti).eq.lz)))) then
   buf(:)=tab_point(:,nb_point)
   tab_point(:,(tti+1):nb_point)=tab_point(:,tti:(nb_point-1))
   tab_point(:,tti)=buf(:)
   perm=1
   exit
  endif  
  endif
 enddo
enddo

!*** VERIFICATION DES SEGMENTS INTERSECTION PLAN BOITE
!write (*,*) "polygone-nb_point->",nb_point
!do ti=1,nb_point,1
! write (*,*) "->",tab_point(:,ti)
!enddo

end subroutine polygone

!__________________________________________________________________________

end module polygones

!
!********************************************************************************************
!

program PolyFront

use polygones

!*** Initialisation bricamat
tabprem(:)=0;tabprem(1)=1;tabprem(2)=2;

!****** LECTURE *****
open(1,file='PolyFront.IN',STATUS='OLD') 
read (1,*) modpol
write(*,*) "mode (0 ok, 1 bug !) ->",modpol
read (1,*) nbit,enx,eny
write(*,*) "nbit enx eny ->",nbit,enx,eny
read (1,*) h,k,l
write(*,*) "h k l ->",h,k,l

!*** Determination des vecteurs pour projeter sur plan hkl
az(1)=h
az(2)=k
az(3)=l
naz=sqrt(az(1)*az(1)+az(2)*az(2)+az(3)*az(3))
az=az/naz
if(h*k.ne.0) then
 ax(1)=-1.d0/h
 ax(2)=1.d0/k
 ax(3)=0.d0
else
 if(h*l.ne.0) then
  ax(1)=-1.d0/h
  ax(2)=0.d0
  ax(3)=1.d0/l
 else
  if(k*l.ne.0) then
   ax(1)=0.d0
   ax(2)=-1.d0/k
   ax(3)=1.d0/l
  else
   if(h.eq.0.and.k.eq.0) then
    ax(1)=1.d0
    ax(2)=1.d0
    ax(3)=0.d0
   else
    if(h.eq.0.and.l.eq.0) then
     ax(1)=1.d0
     ax(2)=0.d0
     ax(3)=1.d0
    else
     if(k.eq.0.and.l.eq.0) then
      ax(1)=0.d0
      ax(2)=1.d0
      ax(3)=1.d0
     else
      stop "PB DE NORMALE -> INDICES TOUS NULS"
     endif
    endif
   endif
  endif
 endif
endif
nax=sqrt(ax(1)*ax(1)+ax(2)*ax(2)+ax(3)*ax(3))
ax=ax/nax
ay=prodvec(az,ax)

read (1,*) Lx,Ly,Lz
write(*,*) 'Lx Ly Lz ->',Lx,Ly,Lz
!write(*,*) 'modulo(0,Lx Ly Lz) ->',modulo(0,Lx),modulo(0,Ly),modulo(0,Lz)
!write(*,*) 'modulo(Lx Ly Lz,Lx Ly Lz) ->',modulo(Lx,Lx),modulo(Ly,Ly),modulo(Lz,Lz)
read (1,*) dxy,dxz
write(*,*) 'dxy dxz ->',dxy,dxz
read (1,*) dyx,dyz
write(*,*) 'dyx dyz ->',dyx,dyz
read (1,*) dzx,dzy
write(*,*) 'dzx dzy ->',dzx,dzy

if((dxy.ne.0.and.dyx.ne.0).or.(dxz.ne.0.and.dzx.ne.0).or.(dyz.ne.0.and.dzy.ne.0)) then
 write (*,*) "shifts incompatibles"
 stop
endif

read (1,*) ni
write(*,*) 'ni ->',ni
close (1)
write(*,*)
write(*,*)
write(*,*)

!*** RECHERCHE DE GEOMETRIE EQUIVALENTE
if(h*k.ne.0) then
 bufachkl=ppmc_dp_m(h,k)
else
 bufachkl=h+k
endif
if(bufachkl*l.ne.0) then
 fachkl_cris=ppmc_dp_m(bufachkl,l)
else
 fachkl_cris=(bufachkl+l)
endif
write (*,*) "Facteur correctif h k l -> ",fachkl_cris

poly_blx=lx*h
poly_bly=ly*k
poly_blz=lz*l
bdxy=dxy*h
bdxz=dxz*h
bdyx=dyx*k
bdyx=dyx*k
bdzx=dzx*l
bdzy=dzy*l
bni=ni

pg1=Epgdc_dp_m(poly_blx,poly_bly)
pg2=Epgdc_dp_m(pg1,poly_blz)
pg1=Epgdc_dp_m(pg2,bdxy)
pg2=Epgdc_dp_m(pg1,bdxz)
pg1=Epgdc_dp_m(pg2,bdyx)
pg2=Epgdc_dp_m(pg1,bdyz)
pg1=Epgdc_dp_m(pg2,bdzx)
pg2=Epgdc_dp_m(pg1,bdzy)
pg1=Epgdc_dp_m(pg2,ni)
pgtrans=pg1

poly_blx=lx
poly_bly=ly
poly_blz=lz
bdxy=dxy
bdxz=dxz
bdyx=dyx
bdyz=dyz
bdzx=dzx
bdzy=dzy
bni=ni

write (*,*) "pgdc des translations : ",pgtrans

! la simplif marche plus
lx=(lx*fachkl_cris)/pgtrans
ly=(ly*fachkl_cris)/pgtrans
lz=(lz*fachkl_cris)/pgtrans
dxy=(dxy*fachkl_cris)/pgtrans
dxz=(dxz*fachkl_cris)/pgtrans
dyx=(dyx*fachkl_cris)/pgtrans
dyz=(dyz*fachkl_cris)/pgtrans
dzx=(dzx*fachkl_cris)/pgtrans
dzy=(dzy*fachkl_cris)/pgtrans
ni=(ni*fachkl_cris)/pgtrans

if(lx.lt.100.or.ly.lt.100.or.lz.lt.100) then
write (*,*) "Augmentation de la taille pour avoir un trace correct"
lx=lx*100
ly=ly*100
lz=lz*100
dxy=dxy*100
dxz=dxz*100
dyx=dyx*100
dyz=dyz*100
dzx=dzx*100
dzy=dzy*100
ni=ni*100
endif

write (*,*) "lx  ",poly_blx,lx
write (*,*) "ly  ",poly_bly,ly
write (*,*) "lz  ",poly_blz,lz
write (*,*) "dxy ",bdxy,dxy
write (*,*) "dxz ",bdxz,dxz
write (*,*) "dyx ",bdyx,dyx
write (*,*) "dyz ",bdyz,dyz
write (*,*) "dzx ",bdzx,dzx
write (*,*) "dzy ",bdzy,dzy
write (*,*) "ni ",bni,ni

hklmax=0
hklmin=0
if(h.gt.0) then
hklmax=hklmax+h*lx
else
hklmin=hklmin+h*lx
endif
if(k.gt.0) then
hklmax=hklmax+k*ly
else
hklmin=hklmin+k*ly
endif
if(l.gt.0) then
hklmax=hklmax+l*lz
else
hklmin=hklmin+l*lz
endif
hklminmax=hklmax-hklmin
write (*,*) "indice entre ->",hklmin," et ",hklmax,"(",hklminmax,")"
if(abs(hklmin).gt.2500000.or.abs(hklmax).gt.2500000)then
 write (*,*) "tableaux sous dimensionnés : 2500000 <",abs(hklmin),"ou",abs(hklmax)
 stop
endif
if(ni.lt.hklmin.or.ni.gt.hklmax) then
 ni=modulo(ni-hklmin,hklminmax)+hklmin
 write (*,*) "Correction de l'indice du plan initial : ",ni
endif

!Si h*k*l<>0 il y a un plan point
!  1 1 1 0  0 0
! -1 1 1 Lx 0 0
!  1-1 1 Lx 0 Lz
!  1 1-1 0  0 Lz
n0=hklmax+1
nn0=hklmax+1
if(h.gt.0.and.k.gt.0.and.l.gt.0)n0=0
if(h.lt.0.and.k.lt.0.and.l.lt.0)n0=0
if(h.lt.0.and.k.gt.0.and.l.gt.0)n0=h*lx
if(h.gt.0.and.k.lt.0.and.l.lt.0)n0=h*lx
if(h.gt.0.and.k.lt.0.and.l.gt.0)n0=h*lx+l*lz
if(h.lt.0.and.k.gt.0.and.l.lt.0)n0=h*lx+l*lz
if(h.gt.0.and.k.gt.0.and.l.lt.0)n0=l*lz
if(h.lt.0.and.k.lt.0.and.l.gt.0)n0=l*lz
nn0=(h*lx+k*ly+l*lz)-n0
write (*,*) "Ne pas utiliser des plans ponctuels ->",n0,nn0
if(ni.eq.n0.or.ni.eq.nn0) then
 ni=modulo(ni+pgtrans-hklmin,hklminmax)+hklmin
 write (*,*) "Correction de l'indice du plan initial : ",ni
endif

fachkl=2
llx=fachkl*lx
lly=fachkl*ly
llz=fachkl*lz

dist_carac=sqrt(1.d0*lx*lx+1.d0*ly*ly+1.d0*lz*lz)
dcd10=dist_carac
write(*,*) 'dcarac ->',dist_carac
la_plus_faible_dist=dist_carac*nbit*nbit
la_plus_grande_dist=0
xmin=la_plus_faible_dist
ymin=la_plus_faible_dist
xmax=-la_plus_faible_dist
ymax=-la_plus_faible_dist

!****** OUVERTURE FICHIER DE SORTIES ET ENTETE POSTSCRIPT *****
open (1,FILE='PolyFront.epsf',STATUS='UNKNOWN') 
open (2,FILE='PolyFront.eps',STATUS='UNKNOWN')
sx=(abs(ax(1)*lx)+abs(ax(2)*ly)+abs(ax(3)*lz))*(enx)
sy=(abs(ay(1)*lx)+abs(ay(2)*ly)+abs(ay(3)*lz))*(eny)
ox=0
oy=0
if(enx.ge.eny)then
 enmax=enx
else
 enmax=eny
endif
ex=nint(600*enx*1.d0/enmax)
ey=nint(600*eny*1.d0/enmax)
!ex=nint(600*enx*1.d0/enmax)
!ey=nint(600*eny*1.d0/enmax)
if((ex*1.d0/sx).lt.(ey*1.d0/sy))then
fx=ex*1.d0/sx
fy=ex*1.d0/abs(sx)*(sy/abs(sy))
ey=ey/(ey*1.d0/sy)*(ex*1.d0/sx)
else
fx=ey*1.d0/abs(sy)*(sx/abs(sx))
fy=ey*1.d0/sy
ex=ex/(ex*1.d0/sx)*(ey*1.d0/sy)
endif
sx=sx/2
sy=sy/2
write (2,'("%!PS-Adobe-3.0 EPSF-3.0")')
write (2,'("%%BoundingBox:",4(1X,I5))') ox,oy,ex,ey
write (2,'("%%BeginDocument: tests.eps")')
write (2,'("%%Title: graph_micmeg.epsf")')
write (2,'("%%DocumentFonts: Helvetica")')
write (2,'("%%EndComments")')
write (2,*) fx,fy," scale"
write (2,'("/Helvetica findfont ",I5," scalefont setfont")') nint(2/fy)

!*** INITIALISATION
planexp(:,:,:)=0 ! n,1,1 -> 0 pas exp, -1 exp, -2 images | n,1,2 nbpt | n,2:7,1:3 points
nit(:)=0 ! plan
nit2(:)=0 ! it suivante
nnit(:)=hklmax+1 ! it
nnit2(:)=0 ! d'av
ptit(:,:)=0 ! shift
ptit2(:,:)=0 ! it suivante
it=1 !nb de polygones
it2=0 !nb de segments des polygones
couple(:)=0 ! indice n -> nsuivant
couple2(:)=0

!on demarre
nit(1)=ni
nbpoly=0

do i=1,nbit !*** Iteration du pavage
 write (*,'(I3,1X,"iterations ++++++++++++++++++++++++++++++++++++ nb polygones : ",I6,1X)') i,nbpoly

bj: do j=1,it !*** Nombre de polygone (ie de segments du dernier polygone)
!write (*,'(I4,1X,"polygones -------------------------------------")') j


!write (*,'("Plan A Explore | Shift ->",1X,I5,"|",3(1X,I5))') nit(j),ptit(:,nit(j))
 call polygone(nit(j),nb)
  
  if(i.eq.1.and.j.eq.1) then
   bary(:)=0
   do ik=1,nb,1
    bary(1)=bary(1)+tab_point(1,ik)
    bary(2)=bary(2)+tab_point(2,ik)
    bary(3)=bary(3)+tab_point(3,ik)
   enddo
    bary(1)=bary(1)/nb
    bary(2)=bary(2)/nb
    bary(3)=bary(3)/nb
    bary(4)=ax(1)*bary(1)+ax(2)*bary(2)+ax(3)*bary(3)
    bary(5)=ay(1)*bary(1)+ay(2)*bary(2)+ay(3)*bary(3)
  endif
  
  fillit=.false.! n,1,1 -> 0 pas exp, -1 exp, -2 images, n,2:7,1:3 points
  
  if(planexp(nit(j),1,1).eq.(-2)) then
!   write (*,*) "Plan Deja Explore, translation differentes bis et cette fois ?"
   fillit=.false.
   do ik=1,nb,1
    if((tab_point(1,ik)+ptit(1,j)).ne.planexp(nit(j),1+ik,1)) fillit=.true.
    if((tab_point(2,ik)+ptit(2,j)).ne.planexp(nit(j),1+ik,2)) fillit=.true.
    if((tab_point(3,ik)+ptit(3,j)).ne.planexp(nit(j),1+ik,3)) fillit=.true.
   enddo
   if(fillit) then
!    write (*,*)"Plan Deja Explore, translation differentes"
   else
!    write (*,*) "Plan Deja Explore mais meme translation"
    cycle bj   
   endif
  endif
  
  if(planexp(nit(j),1,1).eq.(-1)) then
   fillit=.false.
   do ik=1,nb,1
    if((tab_point(1,ik)+ptit(1,j)).ne.planexp(nit(j),1+ik,1)) fillit=.true.
    if((tab_point(2,ik)+ptit(2,j)).ne.planexp(nit(j),1+ik,2)) fillit=.true.
    if((tab_point(3,ik)+ptit(3,j)).ne.planexp(nit(j),1+ik,3)) fillit=.true.
   enddo
   if(fillit) then
    planexp(nit(j),1,1)=-2
!    write (*,*)"Plan Deja Explore, translation differentes"
   else
!    write (*,*) "Plan Deja Explore mais meme translation"
    cycle bj   
   endif
  endif
  if(planexp(nit(j),1,1).eq.0) then
!    write (*,*) "NV Plan"
    fillit=.false.
    planexp(nit(j),1,1)=-1
    planexp(nit(j),1,2)=nb
   planexp(nit(j),1,3)=i
   nbpoly=nbpoly+1
  endif

!*** projections et trace
  do ik=1,nb,1

!write (*,*) "Decalage",j,ptit(1,j),ptit(2,j),ptit(3,j)

   ttab_point(1,ik)=ax(1)*(tab_point(1,ik)+ptit(1,j))+&
                       ax(2)*(tab_point(2,ik)+ptit(2,j))+&
                       ax(3)*(tab_point(3,ik)+ptit(3,j))
   ttab_point(2,ik)=ay(1)*(tab_point(1,ik)+ptit(1,j))+&
                       ay(2)*(tab_point(2,ik)+ptit(2,j))+&
                       ay(3)*(tab_point(3,ik)+ptit(3,j)) 
   ttab_point(1,ik)=sx+ttab_point(1,ik)
   ttab_point(2,ik)=sy+ttab_point(2,ik)
  enddo
  write (2,'(2(I8,1X),"moveto")') ttab_point(1:2,1)
  write (2,'(3(f8.6,1x),"setrgbcolor")') mod(i,dpiv_2)*1.0d0,mod(i,dpiv_3)*0.5d0,mod(i,dpiv_4)*0.33d0/nbit
  write (2,'(I4," setlinewidth")') min(9999,max(nint((nbit-i)*(nbit-i)/abs(fx*nbit*nbit)),1))
  Helvind(:)=0
  do ik=1,nb,1 ! ses voisins
   write (2,'(2(I8,1X),"lineto")') ttab_point(1:2,ik)
   Helvind(:)=Helvind(:)+ttab_point(:,ik)
  enddo  
  Helvind(1:2)=Helvind(1:2)/nb
  Helvind(1)=Helvind(1)-3*nint(2/fy)
  write (2,'(2(I8,1X),"lineto")') ttab_point(1:2,1)  
  if(fillit.and.modpol.eq.0) then
!   write (2,'("0.75  0.75  0.75  setrgbcolor")')
   write (2,'("fill")')
   write (2,'("0.0  0.0  0.0  setrgbcolor")')
  endif
  write (2,'("stroke")')
  write (2,'(2(I8,1X),"moveto")') Helvind(1:2)
  write (2,'("(",(I8,1X),") show")') nit(j)
  write (2,'("stroke")')
  
  if(planexp(nit(j),1,1).eq.(-1).or.(modpol.eq.1.and.planexp(nit(j),1,1).eq.(-2))) then
    
   do ik=1,nb,1 ! segment -> voisins
    planexp(nit(j),1+ik,1)=tab_point(1,ik)+ptit(1,j) !sauvegarde pour comparer coordonnees
    planexp(nit(j),1+ik,2)=tab_point(2,ik)+ptit(2,j)
    planexp(nit(j),1+ik,3)=tab_point(3,ik)+ptit(3,j)

!write (*,*) ik,"|",planexp(nit(j),1+ik,1),planexp(nit(j),1+ik,2),planexp(nit(j),1+ik,3)

    planexp(nit(j),1+ik,4)=ax(1)*(tab_point(1,ik)+ptit(1,j))+&
                           ax(2)*(tab_point(2,ik)+ptit(2,j))+&
                           ax(3)*(tab_point(3,ik)+ptit(3,j))
    planexp(nit(j),1+ik,5)=ay(1)*(tab_point(1,ik)+ptit(1,j))+&
                           ay(2)*(tab_point(2,ik)+ptit(2,j))+&
                           ay(3)*(tab_point(3,ik)+ptit(3,j)) 

    if(xmin.gt.planexp(nit(j),1+ik,4)) xmin=planexp(nit(j),1+ik,4)
    if(ymin.gt.planexp(nit(j),1+ik,5)) ymin=planexp(nit(j),1+ik,5)
    if(xmax.lt.planexp(nit(j),1+ik,4)) xmax=planexp(nit(j),1+ik,4)
    if(ymax.lt.planexp(nit(j),1+ik,5)) ymax=planexp(nit(j),1+ik,5)

    ikk=ik+1
    if(ik.eq.nb) ikk=1
    dx=(tab_point(1,ikk)-tab_point(1,ik))
    dy=(tab_point(2,ikk)-tab_point(2,ik))
    dz=(tab_point(3,ikk)-tab_point(3,ik))

    if (dx.ne.0) then
     if (dy.ne.0) then
      pg1=Epgdc_dp_m(dx,dy)
     else
      pg1=dx
     endif
    else
     if (dy.ne.0) then
      pg1=dy
     else
      pg1=0
     endif
    endif
    if (pg1.ne.0) then
     if (dz.ne.0) then
      pg2=Epgdc_dp_m(pg1,dz)
     else
      pg2=abs(pg1)
     endif
    else
     if (dz.ne.0) then
      pg2=abs(dz)
     else
      pg2=2
      write (*,*) "Pb de pg2 : plan ponctuel"
      write (*,*) "éviter les plans ",n0*pgtrans," et ",nn0*pgtrans
      stop
     endif
    endif
    
    if(pg2.gt.1) then
      demi=fachkl
    else
     if(pg2.eq.1) then
      write (*,*) tab_point(1,ikk)*h+tab_point(2,ikk)*k+tab_point(3,ikk)*l,"Segment unitaire -:(->",dx,dy,dz,"|",pg2
      demi=1    
     else   
      write (*,*) tab_point(1,ikk)*h+tab_point(2,ikk)*k+tab_point(3,ikk)*l,"Plan ponctuel -:(->",dx,dy,dz,"|",pg2
      write (*,*) "éviter les plans ",n0*pgtrans," et ",nn0*pgtrans
      stop
     endif
    endif
    
    dx=dx/pg2*demi
    dy=dy/pg2*demi
    dz=dz/pg2*demi

    x=tab_point(1,ik)*fachkl+dx
    y=tab_point(2,ik)*fachkl+dy
    z=tab_point(3,ik)*fachkl+dz
    if (x.lt.0.or.x.gt.llx.or.y.lt.0.or.y.gt.lly.or.z.lt.0.or.z.gt.llz) then
     x=tab_point(1,ik)*fachkl-dx
     y=tab_point(2,ik)*fachkl-dy
     z=tab_point(3,ik)*fachkl-dz
    endif

    if (x.lt.0.or.x.gt.llx.or.y.lt.0.or.y.gt.lly.or.z.lt.0.or.z.gt.llz) then     
     write (*,*) "segment -> nv polygone",pg1,pg2,demi
     write (*,*) tab_point(1,ikk),tab_point(1,ik),dx,x
     write (*,*) tab_point(2,ikk),tab_point(2,ik),dy,y
     write (*,*) tab_point(3,ikk),tab_point(3,ik),dz,z
     stop "pb de choix du point sur le segment"
    endif

!if(nit(j).ne.0.or.(h*k*l).eq.0) then
if(nit(j).ne.n0.and.nit(j).ne.nn0) then
 n = plan(x,y,z,dx,dy,dz)

!write (*,*) "PLAN",n,"(",nnit(j),")"

    if (nnit(j).ne.n) then !pour ne pas revenir sur nos pas
      it2=it2+1
      nit2(it2)=n
      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
      couple2(it2)=ik
      ptit2(1,it2)=dx+ptit(1,j)
      ptit2(2,it2)=dy+ptit(2,j)
      ptit2(3,it2)=dz+ptit(3,j)
    else
     planexp(n,couple(j)+1,6)=1
     planexp(nit(j),ik+1,6)=1
    endif
else

write (*,*) "PLAN PONCTUEL",j,ptit(1,j),ptit(2,j),ptit(3,j)

      it2=it2+1
dx=-Lx
dy=-dxy
if(dy.gt.0)dy=-Ly+dxy
dz=-dxz
if(dz.gt.0)dz=-Lz+dxz
      nit2(it2)=-(h*dx+k*dy+l*dz)
      if(nit2(it2).lt.hklmin)nit2(it2)=nit2(it2)+hklminmax
      if(nit2(it2).gt.hklmax)nit2(it2)=nit2(it2)-hklminmax
      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
      couple2(it2)=ik
      ptit2(1,it2)=dx+ptit(1,j)
      ptit2(2,it2)=dy+ptit(2,j)
      ptit2(3,it2)=dz+ptit(3,j)

write (*,*) "Plan PONCTUEL dim x -h*Lx->",nit2(it2)
stop
!*** Suite : pose encore quelques pbs et sert uniquement si le plan initial est no ou nno
!    à régler si du temps à perdre...

      it2=it2+1
dy=-Ly
dx=-dyx
if(dx.gt.0)dx=-Lx+dyx
dz=-dyz
if(dz.gt.0)dz=-Lz+dyz
      nit2(it2)=-(h*dx+k*dy+l*dz)
      if(nit2(it2).lt.hklmin)nit2(it2)=nit2(it2)+hklminmax
      if(nit2(it2).gt.hklmax)nit2(it2)=nit2(it2)-hklminmax
      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
      couple2(it2)=ik
      ptit2(1,it2)=dx+ptit(1,j)
      ptit2(2,it2)=dy+ptit(2,j)
      ptit2(3,it2)=dz+ptit(3,j)

write (*,*) "Plan PONCTUEL dim y -k*Ly->",nit2(it2)

      it2=it2+1
dz=-Lz
dx=-dzx
if(dx.gt.0)dx=-Lx+dzx
dy=-dzy
if(dy.gt.0)dy=-Ly+dzy
      nit2(it2)=-(h*dx+k*dy+l*dz)
      if(nit2(it2).lt.hklmin)nit2(it2)=nit2(it2)+hklminmax
      if(nit2(it2).gt.hklmax)nit2(it2)=nit2(it2)-hklminmax
      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
      couple2(it2)=ik
      ptit2(1,it2)=dx+ptit(1,j)
      ptit2(2,it2)=dy+ptit(2,j)
      ptit2(3,it2)=dz+ptit(3,j)

write (*,*) "Plan PONCTUEL dim z -l*Lz->",nit2(it2)

!      it2=it2+1
!      nit2(it2)=h*Lx+k*dxy+l*dxz+l*Lz+k*dzy+h*dzx
!      if(nit2(it2).lt.hklmin)nit2(it2)=nit2(it2)+hklminmax
!      if(nit2(it2).gt.hklmax)nit2(it2)=nit2(it2)-hklminmax
!dx=-Lx-dzx
!dy=-dxy-dzy
!dz=-Lz-dxz
!      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
!      couple2(it2)=ik
!      ptit2(1,it2)=dx+ptit(1,j)
!      ptit2(2,it2)=dy+ptit(2,j)
!      ptit2(3,it2)=dz+ptit(3,j)
!
!write (*,*) "Plan POINT dim x h*Lx->",nit2(it2)

!      it2=it2+1
!      nit2(it2)=k*Ly+h*dyx+l*dyz+h*Lx+k*dxy+l*dxz
!      if(nit2(it2).lt.hklmin)nit2(it2)=nit2(it2)+hklminmax
!      if(nit2(it2).gt.hklmax)nit2(it2)=nit2(it2)-hklminmax
!dx=-Lx-dyx
!dy=-Ly-dxy
!dz=-dyz-dxz
!      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
!      couple2(it2)=ik
!      ptit2(1,it2)=dx+ptit(1,j)
!      ptit2(2,it2)=dy+ptit(2,j)
!      ptit2(3,it2)=dz+ptit(3,j)
!
!write (*,*) "Plan POINT dim y k*Ly->",nit2(it2)

!      it2=it2+1
!      nit2(it2)=l*Lz+h*dzx+k*dzy+k*Ly+h*dyx+l*dyz
!      if(nit2(it2).lt.hklmin)nit2(it2)=nit2(it2)+hklminmax
!      if(nit2(it2).gt.hklmax)nit2(it2)=nit2(it2)-hklminmax
!dx=-dzx-dyx
!dy=-dzy-Ly
!dz=-Lz-dyz
!      nnit2(it2)=nit(j) !pour ne pas revenir sur nos pas a la prochaine it
!      couple2(it2)=ik
!      ptit2(1,it2)=dx+ptit(1,j)
!      ptit2(2,it2)=dy+ptit(2,j)
!      ptit2(3,it2)=dz+ptit(3,j)
!
!write (*,*) "Plan POINT dim z l*Lz->",nit2(it2)

endif
    
   enddo
  endif
 enddo bj

 write (*,'("=> ",I4," polygones pour la prochaine iterations")') it2
! do ik=1,it2,1 ! ecriture des voisins
!  write (*,*) nit2(ik)
! enddo

! passage au nouveau polygone
  do ik=1,it2,1 ! ecriture des voisins
   nit(ik)   = nit2(ik)
   nnit(ik)  = nnit2(ik)   
   ptit(:,ik)  = ptit2(:,ik)   
   couple(ik) = couple2(ik)
   couple2(ik)= 0
   nit2(ik)  = 0
   nnit2(ik) = 0
   ptit2(:,ik) = 0   
  enddo
  it=it2
  it2=0
enddo

if(it.ne.0)then
 write (*,'("PAVAGE INACHEVE")')
 write (*,'("=> ",I4," polygones pour la prochaine iterations")') it2
 write (*,'("AUGMENTER SA TAILLE ET LE NB D ITERATIONS")')
endif

!*** Fin postscript de travail
write (2,'("showpage")')
write (2,'("%EOF")')
write (2,'("%EndDocument")')
close (2)

write (*,*) "nbpoly -> ",nbpoly
!****** POLYGONE FRONTIERE
poly_front(:,:,:)=0
buf_i=0
i=0

!* mémorisation

write (*,*) "Construction du Polygone Frontière",hklmin,hklmax
do j=hklmin,hklmax,1
 if(planexp(j,1,1).ne.0) then
  do ik=2,planexp(j,1,2)+1,1
    i=i+1
    poly_front(i,1,1:5)=planexp(j,ik,1:5)
    buf_ik=ik+1
    if(buf_ik.eq.planexp(j,1,2)+2)buf_ik=2
    poly_front(i,2,1:5)=planexp(j,buf_ik,1:5)
  enddo
 endif
enddo

!* tri
  do jk=1,i,1
   do ik=jk+1,i,1
     v1(1:3)=poly_front(jk,2,1:3)-poly_front(jk,1,1:3)
     v2(1:3)=poly_front(ik,2,1:3)-poly_front(ik,1,1:3)

     v3(1:3)=poly_front(jk,1,1:3)-poly_front(ik,1,1:3)
!     vv1(1:3)=poly_front(jk,1,1:3)-poly_front(ik,1,1:3)

     nor1=norivect(prodivec(v1,v2))
     nor2=norivect(prodivec(v1,v3))
     if(nor1.eq.0.and.nor2.eq.0) then
     vv2(1:3)=poly_front(jk,2,1:3)-poly_front(ik,2,1:3)

     v4(1:3) =poly_front(jk,2,1:3)-poly_front(ik,1,1:3)
!     vv3(1:3)=poly_front(jk,2,1:3)-poly_front(ik,1,1:3)

     vv4(1:3)=poly_front(jk,1,1:3)-poly_front(ik,2,1:3) 
      nor1=norivect(v3)
      nor2=norivect(vv2)
      nor3=norivect(v4)
      nor4=norivect(vv4)
      if((nor1.eq.0.and.nor2.eq.0).or.(nor3.eq.0.and.nor4.eq.0))then
       poly_front(ik,:,4:5)=0
       poly_front(jk,:,4:5)=0
      endif     
      nor3=prodisca(v1,v2)
      nor4=prodisca(v1,v3)
      nor5=prodisca(v1,v4)
      if(((nor4*nor5).lt.0).or.(nor4.eq.0.and.(nor3*nor5).gt.0).or.(nor5.eq.0.and.(nor3*nor4).gt.0).or.&
         ((nor3*nor4).gt.0.and.(abs(nor4).lt.abs(nor3))).or.&
         ((nor3*nor5).gt.0.and.(abs(nor5).lt.abs(nor3))))then
       poly_front(ik,:,4:5)=0
       poly_front(jk,:,4:5)=0   
      endif
     endif
   enddo  
  enddo

  buf_i=i
  jk=1
  do while (jk.ne.buf_i+1)
   if (jk.eq.(buf_i+1)) exit
   if(poly_front(jk,1,4).eq.0.and.poly_front(jk,1,5).eq.0.and.&
      poly_front(jk,2,4).eq.0.and.poly_front(jk,2,5).eq.0) then
    if(jk.ne.buf_i) poly_front(jk:(buf_i-1),1:2,1:5)=poly_front((jk+1):buf_i,1:2,1:5)
    buf_i=buf_i-1
   else
    jk=jk+1
   endif
  enddo
  i=buf_i

bary2(:)=0
do ik=1,i,1
 bary2(1)=bary2(1)+poly_front(ik,1,1)+poly_front(ik,2,1)
 bary2(2)=bary2(2)+poly_front(ik,1,2)+poly_front(ik,2,2)
 bary2(3)=bary2(3)+poly_front(ik,1,3)+poly_front(ik,2,3)
enddo
bary2(1)=bary2(1)/(2*i)
bary2(2)=bary2(2)/(2*i)
bary2(3)=bary2(3)/(2*i)
bary2(4)=ax(1)*bary2(1)+ax(2)*bary2(2)+ax(3)*bary2(3)
bary2(5)=ay(1)*bary2(1)+ay(2)*bary2(2)+ay(3)*bary2(3)

la_plus_faible_dist2=dist_carac*nbit*nbit
la_plus_grande_dist2=0
i=buf_i 
do ik=1,i,1 
 do jk=1,2,1 
  buf_pt_front(1)=poly_front(ik,jk,1)-bary2(1)
  buf_pt_front(2)=poly_front(ik,jk,2)-bary2(2)
  buf_pt_front(3)=poly_front(ik,jk,3)-bary2(3)
  buf_dist=sqrt(buf_pt_front(1)*buf_pt_front(1)+buf_pt_front(2)*buf_pt_front(2)+buf_pt_front(3)*buf_pt_front(3))
  buf_pt_front2(1)=poly_front(ik,jk,1)-bary(1)
  buf_pt_front2(2)=poly_front(ik,jk,2)-bary(2)
  buf_pt_front2(3)=poly_front(ik,jk,3)-bary(3)
  buf_dist2=sqrt(buf_pt_front2(1)*buf_pt_front2(1)+buf_pt_front2(2)*buf_pt_front2(2)+buf_pt_front2(3)*buf_pt_front2(3))
  if(la_plus_faible_dist2.gt.buf_dist2) then
   la_plus_faible_dist2 = buf_dist2
   pt_dlpf2(1)=poly_front(ik,jk,1)
   pt_dlpf2(2)=poly_front(ik,jk,2)
   pt_dlpf2(3)=poly_front(ik,jk,3)
   pt_dlpf2(4)=ax(1)*(poly_front(ik,jk,1))+&
	       ax(2)*(poly_front(ik,jk,2))+&
	       ax(3)*(poly_front(ik,jk,3))
   pt_dlpf2(5)=ay(1)*(poly_front(ik,jk,1))+&
	       ay(2)*(poly_front(ik,jk,2))+&
	       ay(3)*(poly_front(ik,jk,3))
  endif
  if(la_plus_faible_dist.gt.buf_dist) then
   la_plus_faible_dist = buf_dist
   pt_dlpf(1)=poly_front(ik,jk,1)
   pt_dlpf(2)=poly_front(ik,jk,2)
   pt_dlpf(3)=poly_front(ik,jk,3)
   pt_dlpf(4)=ax(1)*(poly_front(ik,jk,1))+&
 	      ax(2)*(poly_front(ik,jk,2))+&
 	      ax(3)*(poly_front(ik,jk,3))
   pt_dlpf(5)=ay(1)*(poly_front(ik,jk,1))+&
 	      ay(2)*(poly_front(ik,jk,2))+&
 	      ay(3)*(poly_front(ik,jk,3))
  endif
  if(la_plus_grande_dist2.lt.buf_dist2) then
   la_plus_grande_dist2 = buf_dist2
   pt_dlpg2(1)=poly_front(ik,jk,1)
   pt_dlpg2(2)=poly_front(ik,jk,2)
   pt_dlpg2(3)=poly_front(ik,jk,3)
   pt_dlpg2(4)=ax(1)*(poly_front(ik,jk,1))+&
 	      ax(2)*(poly_front(ik,jk,2))+&
 	      ax(3)*(poly_front(ik,jk,3))
   pt_dlpg2(5)=ay(1)*(poly_front(ik,jk,1))+&
 	       ay(2)*(poly_front(ik,jk,2))+&
 	       ay(3)*(poly_front(ik,jk,3))
  endif   
  if(la_plus_grande_dist.lt.buf_dist) then
   la_plus_grande_dist = buf_dist
   pt_dlpg(1)=poly_front(ik,jk,1)
   pt_dlpg(2)=poly_front(ik,jk,2)
   pt_dlpg(3)=poly_front(ik,jk,3)
   pt_dlpg(4)=ax(1)*(poly_front(ik,jk,1))+&
	      ax(2)*(poly_front(ik,jk,2))+&
	      ax(3)*(poly_front(ik,jk,3))
   pt_dlpg(5)=ay(1)*(poly_front(ik,jk,1))+&
	      ay(2)*(poly_front(ik,jk,2))+&
	      ay(3)*(poly_front(ik,jk,3))
  endif   
 enddo  
enddo

write (*,*) "BARYCENTRE DU POLYGONE DE DEPART*************"

write (*,*) "PLUS FAIBLE DISTANCE ->",la_plus_faible_dist2,"(",la_plus_faible_dist2/dist_carac,")"
write (*,*) "sa direction ->",pt_dlpf2(1)-bary(1),pt_dlpf2(2)-bary(2),pt_dlpf2(3)-bary(3)
write (*,*) "PLUS GRANDE DISTANCE ->",la_plus_grande_dist2,"(",la_plus_grande_dist2/dist_carac,")"
write (*,*) "sa direction ->",pt_dlpg2(1)-bary(1),pt_dlpg2(2)-bary(2),pt_dlpg2(3)-bary(3)

write (*,*) "BARYCENTRE DU GRAND POLYGONE*****************"

write (*,*) "PLUS FAIBLE DISTANCE ->",la_plus_faible_dist,"(",la_plus_faible_dist/dist_carac,")"
write (*,*) "sa direction ->",pt_dlpf(1)-bary2(1),pt_dlpf(2)-bary2(2),pt_dlpf(3)-bary2(3)
write (*,*) "PLUS GRANDE DISTANCE ->",la_plus_grande_dist,"(",la_plus_grande_dist/dist_carac,")"
write (*,*) "sa direction ->",pt_dlpg(1)-bary2(1),pt_dlpg(2)-bary2(2),pt_dlpg(3)-bary2(3)

!****** ECRITURE POSTSCRIPT PROPRE *****_______________________________________________________________________________

planexp(:,2:7,4)=planexp(:,2:7,4)-xmin
planexp(:,2:7,5)=planexp(:,2:7,5)-ymin
poly_front(:,:,4)=poly_front(:,:,4)-xmin
poly_front(:,:,5)=poly_front(:,:,5)-ymin
pt_dlpf(4)=pt_dlpf(4)-xmin
pt_dlpf(5)=pt_dlpf(5)-ymin
pt_dlpg(4)=pt_dlpg(4)-xmin
pt_dlpg(5)=pt_dlpg(5)-ymin
pt_dlpf2(4)=pt_dlpf2(4)-xmin
pt_dlpf2(5)=pt_dlpf2(5)-ymin
pt_dlpg2(4)=pt_dlpg2(4)-xmin
pt_dlpg2(5)=pt_dlpg2(5)-ymin
xmax=xmax-xmin
ymax=ymax-ymin
ox=bary(4)-xmin
oy=bary(5)-ymin
oox=bary2(4)-xmin
ooy=bary2(5)-ymin
xmin = 0
ymin = 0
if(xmax.ge.ymax)then
 fx=550.d0/xmax
else
 if(ymax*1.d0/xmax.gt.1.33) then
  fx=650.d0/ymax
 else
  fx=450.d0/xmax 
 endif
endif
fy=fx*5*dcd10/20.d0
ex=xmax*fx+100
ey=ymax*fx+100
planexp(:,2:7,4)=planexp(:,2:7,4)+(50/fx)
planexp(:,2:7,5)=planexp(:,2:7,5)+(50/fx)
poly_front(:,:,4)=poly_front(:,:,4)+(50/fx)
poly_front(:,:,5)=poly_front(:,:,5)+(50/fx)
pt_dlpf(4)=pt_dlpf(4)+(50/fx)
pt_dlpf(5)=pt_dlpf(5)+(50/fx)
pt_dlpg(4)=pt_dlpg(4)+(50/fx)
pt_dlpg(5)=pt_dlpg(5)+(50/fx)
pt_dlpf2(4)=pt_dlpf2(4)+(50/fx)
pt_dlpf2(5)=pt_dlpf2(5)+(50/fx)
pt_dlpg2(4)=pt_dlpg2(4)+(50/fx)
pt_dlpg2(5)=pt_dlpg2(5)+(50/fx)
ox=ox+(50/fx)
oy=oy+(50/fx)
oox=oox+(50/fx)
ooy=ooy+(50/fx)

! ENTETE

write (1,'("%!PS-Adobe-3.0 EPSF-3.0")')
write (1,'("%%BoundingBox:",4(1X,I5))') 0,0,ex,ey
write (1,'("%%BeginDocument: PolyFront.eps")')
write (1,'("%%Title: PolyFront.epsf")')
write (1,'("%%DocumentFonts: Helvetica")')
write (1,'("%%EndComments")')
write (1,*) fx,fx," scale"
write (1,'("/Helvetica findfont ",I8," scalefont setfont")') nint(10*dcd10/fy/(0.5d0*nbpoly)**(1.d0/3.d0))
write (1,'(3(f8.6,1x),"setrgbcolor")') 0.d0,0.d0,0.d0

! PAVAGE

do j=1,nbit,1
do i=hklmin,hklmax,1
 if(i.ne.0.and.i.ne.planexp(i,1,1).and.planexp(i,1,3).eq.j) then
!  write (1,'(3(f8.6,1x),"setrgbcolor")') mod(j,dpiv_2)*0.6d0+0.2d0,mod(j,dpiv_3)*0.3d0+0.2d0,mod(j,dpiv_4)*0.2d0+0.2d0
  write (1,'(I6," setlinewidth")') nint(0.25d0/fx)
  write (1,'(2(I8,1X),"moveto")') planexp(i,2,4),planexp(i,2,5)
  Helvind(:)=planexp(i,2,4:5)
  do ik=3,planexp(i,1,2)+1,1 ! ses voisins
   write (1,'(2(I8,1X),"lineto")') planexp(i,ik,4),planexp(i,ik,5)
   Helvind(:)=Helvind(:)+planexp(i,ik,4:5)
  enddo  
  Helvind(1:2)=Helvind(1:2)/planexp(i,1,2)
  Helvind(1)=Helvind(1)-2*nint(dcd10/fy)
  write (1,'(2(I8,1X),"lineto")') planexp(i,2,4),planexp(i,2,5)
  write (1,'("stroke")')
! NUMERO DU PLAN  
  write (1,'(2(I8,1X),"moveto")') Helvind(1:2)
  write (1,'("(",(I7),") show")') i
  write (1,'("stroke")')
 endif
enddo
enddo

!*** LPM

write (1,'(I6," setlinewidth")') nint(0.75d0/fx)
write (1,'(3(f8.6,1x),"setrgbcolor")') 1.d0,0.d0,0.d0
write (1,'(2(I8,1X),"moveto")') ox,oy
write (1,'(2(I8,1X),"lineto")') pt_dlpf2(4),pt_dlpf2(5)
write (1,'("stroke")')
write (1,'(3(f8.6,1x),"setrgbcolor")') 0.d0,1.d0,0.d0
write (1,'(2(I8,1X),"moveto")') ox,oy
write (1,'(2(I8,1X),"lineto")') pt_dlpg2(4),pt_dlpg2(5)
write (1,'("stroke")')

write (1,'(I6," setlinewidth")') nint(0.5d0/fx)
write (1,'(3(f8.6,1x),"setrgbcolor")') 1.d0,0.d0,0.d0
write (1,'(2(I8,1X),"moveto")') oox,ooy
write (1,'(2(I8,1X),"lineto")') pt_dlpf(4),pt_dlpf(5)
write (1,'("stroke")')
write (1,'(3(f8.6,1x),"setrgbcolor")') 0.d0,1.d0,0.d0
write (1,'(2(I8,1X),"moveto")') oox,ooy
write (1,'(2(I8,1X),"lineto")') pt_dlpg(4),pt_dlpg(5)
write (1,'("stroke")')

! CONTOUR

write (1,'(3(f8.6,1x),"setrgbcolor")') 0.d0,0.d0,0.d0
write (1,'(I6," setlinewidth")') nint(2.d0/fx)

do jk=1,buf_i,1
 write (1,'(2(I8,1X),"moveto")') poly_front(jk,1,4),poly_front(jk,1,5)
 write (1,'(2(I8,1X),"lineto")') poly_front(jk,2,4),poly_front(jk,2,5)
 write (1,'("stroke")')
enddo  

! LEGENDE

write (1,'(I6," setlinewidth")') nint(0.5d0/fx)
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(19*dcd10/fy)
write (1,'(2(I8,1X),"lineto")') nint(7*dcd10/fy),nint(19*dcd10/fy)
write (1,'("stroke")')
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(19*dcd10/fy)
write (1,'(2(I8,1X),"lineto")') nint(5*dcd10/fy),nint(21*dcd10/fy)
write (1,'("stroke")')
write (1,'(2(I8,1X),"moveto")') nint(10*dcd10/fy),nint(19*dcd10/fy)
write (1,'(2(I8,1X),"lineto")') nint(10*dcd10/fy+dist_carac),nint(19*dcd10/fy)
write (1,'("stroke")')

write (1,'("/Helvetica findfont ",I8," scalefont setfont")') nint(2*dcd10/fy)

write (1,'(2(I8,1X),"moveto")') nint(8*dcd10/fy),nint(19*dcd10/fy)
write (1,'("(x) show")')
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(22*dcd10/fy)
write (1,'("(y) show")')
write (1,'(2(I8,1X),"moveto")') nint((10*dcd10+dcd10)/fy+dist_carac),nint(19*dcd10/fy)
write (1,'("(lxyz) show")')

write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(16*dcd10/fy)
write (1,'("(x ",3(f8.4,1X),") show")') ax(1),ax(2),ax(3)
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(14*dcd10/fy)
write (1,'("(y ",3(f8.4,1X),") show")') ay(1),ay(2),ay(3)
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(12*dcd10/fy)
write (1,'("(z (hkl) ",3(I8,1X),") show")') h,k,l
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(10*dcd10/fy)
write (1,'("(lx/y/z ",3(I8,1X),") show")') lx,ly,lz
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(8*dcd10/fy)
write (1,'("(dxy/xz dyx/yz dzx/zy",6(I8,1X),") show")') dxy,dxz,dyx,dyz,dzx,dzy
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(5*dcd10/fy)
write (1,'("(lmin ",(f8.4,1X),") show")') la_plus_faible_dist2/dist_carac
write (1,'(2(I8,1X),"moveto")') nint((5*dcd10+8*2*dcd10)/fy),nint(3*dcd10/fy)
write (1,'("(nbpolygone ",(I8,1X),") show")') nbpoly
write (1,'(2(I8,1X),"moveto")') nint(5*dcd10/fy),nint(3*dcd10/fy)
write (1,'("(lmax ",(f8.4,1X),") show")') la_plus_grande_dist2/dist_carac
write (1,'("stroke")')

write (1,'("showpage")')
write (1,'("%EOF")')
write (1,'("%EndDocument")')
close (1)


end program PolyFront
