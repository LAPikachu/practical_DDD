! MobiDiC (in the place of Mobile Dislocation Colony) is program of DD (Dislocation Dynamics)
! developed at the 'CEA, DAM, DIF' CEA, FRANCE
! Copyright (C) 2003, 2006, 2007, 2008 Ronan MADEC, Hervé TROADEC, Vincent LIARD,
! Aymeric CANTON, Anthony DYAN
!
! MobiDiC is initially based on mM (in the place of microMEGAS)
! developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
! Copyright (C) 1993, 1998, 2000 Benoit DEVINCRE, Ladislas KUBIN, Marc CONDAT,
! Ronan MADEC, Ghiath MONNET
!
! This file is part of MobiDiC
!
! MobiDiC is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! any later version.
!
! MobiDiC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MobiDiC; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
!
! For information contact by Email Ronan MADEC : ronan.madec@cea.fr 

!#####################################
!# Decomposition en nombres premiers #
!#####################################

function dec_prem_m(a,tabdec)

implicit none

integer(DPIM),intent(in) :: a
integer(DPIM):: aa,aaa,n,i,k,j,dec_prem_m,nb
integer(DPIM),intent(out) :: tabdec(dimPREM)

dec_prem_m=-1

!tabprem(:)=0;tabprem(1)=1;tabprem(2)=2;
!
!SONGER A INITIALISER TABPREM

aa = abs (a)
aaa = aa
if(aaa.eq.0)then
write (*,*) "dec_prem : entiers > 0"
stop
endif
n = 2
i = 1
tabdec(1:dimPREM)=0

do while(aaa.ne.1)
nb=0
if(modulo(aaa,tabprem(n)).eq.0) then
tabdec(i)=tabprem(n)
aaa=aaa/tabdec(i)
i=i+1
nb=nb+1
else
nb=0
endif
if(nb.eq.0) n=n+1
if(tabprem(n).eq.0) then
k=tabprem(n-1)
j=0
do while(j.ne.(n-1))
k=k+1
j=2
do while(modulo(k,tabprem(j)).ne.0.and.j.ne.(n-1))
j=j+1
enddo
enddo
tabprem(n)=k

if(n.eq.dimPrem) then
write (*,*) "prem(",n,")=",k," pour ",a,"->",tabdec(1:i-1),aaa
stop "bricamat_dec_prem.inc : tableau trop petit"
endif

endif
enddo

dec_prem_m=i-1

return

end function dec_prem_m

!###########################################################################

function dec_prem(a,tabdec)

implicit none

integer(DPI),intent(in) :: a
integer(DPI) :: aa,n,i,j,dec_prem,nb
integer(DPI),intent(out) :: tabdec(dimPREM)
integer(DPIM) :: k,aaa

dec_prem=-1

!tabprem(:)=0;tabprem(1)=1;tabprem(2)=2;
!
!SONGER A INITIALISER TABPREM

aa = abs (a)
aaa = aa
if(aaa.eq.0)then
write (*,*) "dec_prem : entiers > 0"
stop
endif
n = 2
i = 1
tabdec(1:dimPREM)=0

do while(aaa.ne.1)
nb=0
if(modulo(aaa,tabprem(n)).eq.0) then
tabdec(i)=tabprem(n)
aaa=aaa/tabdec(i)
i=i+1
nb=nb+1
else
nb=0
endif
if(nb.eq.0) n=n+1
if(tabprem(n).eq.0) then
k=tabprem(n-1)
j=0
do while(j.ne.(n-1))
k=k+1
j=2
do while(modulo(k,tabprem(j)).ne.0.and.j.ne.(n-1))
j=j+1
enddo
enddo
tabprem(n)=k

if(n.eq.dimPrem) then
write (*,*) "prem(",n,")=",k," pour ",a,"->",tabdec(1:i-1),aaa
stop "bricamat_dec_prem.inc : tableau trop petit"
endif

endif
enddo

dec_prem=i-1

return

end function dec_prem

!***********************************
!* PGCD via l'algorithme d'Euclide *
!***********************************
!
! Divise par trois le temps de la phase d'initialisation
! (cas test 12 systemes de 12 vecteurs sur Xeon)
! NB Ppgcd etait deja bcp plus efficace que la version de dd
!

function Epgdc_dp(aa,bb)

implicit none

integer(DPI),intent (in)  :: aa,bb
integer(DPI)              :: a,b,c,Epgdc_dp

if(bb.gt.aa)then
a=bb
b=aa
else
a=aa
b=bb
endif
Epgdc_dp=1

if(a*b.eq.0) then
Epgdc_dp=max(abs(a),abs(b))
if(Epgdc_dp.eq.0)then
write (*,*)"Uniquement des zero comme argument de pgdc_dp"
stop
endif
else

do while (b.ne.0)
c=mod(a,b)
a=b
b=c
enddo
Epgdc_dp=abs(a)

endif

return

end function Epgdc_dp

!###########################################################################

function Epgdc_dp_m(aaa,bbb)

implicit none

integer(DPIM),intent (in)  :: aaa,bbb
integer(DPIM)              :: a,b,c,Epgdc_dp_m,aa,bb

!write (*,*) "Epgdc_dp_m(",aa,bb,")"

aa=abs(aaa)
bb=abs(bbb)

if(bb.gt.aa)then
a=bb
b=aa
else
a=aa
b=bb
endif
Epgdc_dp_m=1

if(a*b.eq.0) then
Epgdc_dp_m=max(abs(a),abs(b))
if(Epgdc_dp_m.eq.0)then
write (*,*) "Uniquement des zero comme argument de pgdc_dp"
stop
endif
else

do while (b.ne.0)
c=mod(a,b)
a=b
b=c
enddo
Epgdc_dp_m=abs(a)

endif

return

end function Epgdc_dp_m

!**********************************************************************
!* Verifie si un vecteur deplacement est compatible avec la cristallo *
!**********************************************************************
!
! v=a*va+b*vb
!

subroutine decvec2d(a,b,iva,ivb,iv)

implicit none

integer(kind=DPI) :: iva(3),ivb(3),iv(3)
real(kind=DPR) :: va(3),vb(3),v(3)
real(kind=DPR) :: a,b

if(iv(1).eq.0.and.iv(2).eq.0.and.iv(3).eq.0)write (*,*) "decvec2d : iv = 0 !!!"

va=iva
vb=ivb
v=iv

if(iva(1)*ivb(2).eq.iva(2)*ivb(1)) then
va(1)=iva(3)
vb(1)=ivb(3)
v(1)=iv(3)
endif

if(va(1).ne.0.and.vb(2)*va(1).ne.vb(1)*va(2)) then
b=(v(2)-v(1)*va(2)/va(1))/(vb(2)-vb(1)/va(1)*va(2))
a=(v(1)-b*vb(1))/va(1)
else
if(vb(1).ne.0.and.va(2)*vb(1).ne.va(1)*vb(2)) then
a=(v(2)-v(1)*vb(2)/vb(1))/(va(2)-va(1)/vb(1)*vb(2))
b=(v(1)-a*va(1))/vb(1)
else
if(vb(1).eq.0.and.va(1).eq.0)then
stop "bricamat_decvec2d.inc : Les vecteurs doivent former une base."
endif
endif
endif

if(nint(a*va(1)+b*vb(1)).ne.v(1).or.nint(a*va(2)+b*vb(2)).ne.v(2).or.nint(a*va(2)+b*vb(2)).ne.v(2)) then

write (*,*) a,"|",iva
write (*,*) b,"|",ivb
write (*,*) a*iva+b*ivb,"|",iv

stop "bricamat_decvec2d.inc : pb de resolution."

endif

end subroutine decvec2d

!##############################################################################
!#                                                                            #
!# trouve un multiple entier pour transformer un reel en entier s'il existe   #
!#                                                                            #
!# si 0 -> ca a pas marche !!!                                                #
!#                                                                            #
!##############################################################################

function intfrac(sr)

implicit none

integer(DPI)          :: intfrac
real(DPR),intent (in) :: sr
real(DPR)             :: r,rb,rbb,rbbb

r=abs(sr) ! quelque soit le signe le pb est le meme

if (abs(PartEnt(r)-r).le.PRECISION) then
intfrac = 1
else
rb = max(1/abs(r-PartEnt(r)),1/abs(1-r+PartEnt(r)))
if (abs(PartEnt(rb)-rb).le.PRECISION) then
intfrac=PartEnt(rb)
else
rbb = rb*PartEnt(rb)
rbb = PartEnt(rbb)
if (abs(PartEnt(rbb*r)-(r*rbb)).le.PRECISION) then
intfrac=rbb
else
rbbb = rbb*PartEnt(rbb)
rbbb = PartEnt(rbbb)
if (abs(PartEnt(rbbb*r)-(r*rbbb)).le.PRECISION) then
intfrac=rbbb
else
intfrac = 0
write (*,*) 'IntFrac : Pas de solution pour ',r,rbbb,PartEnt(rbbb*r),(r*rbbb)
stop
endif
endif
endif
endif
return

end function intfrac

!##############################################################################

function int_red_frac(a,b)

implicit none

integer(DPI)          :: int_red_frac,a,b,aa,bb,aaa,bbb,divpar

aa=abs(a) ! quelque soit le signe le pb est le meme
bb=abs(b)
divpar=Epgdc_dp(aa,bb)
aaa=aa/divpar
bbb=bb/divpar
int_red_frac=bbb

!if(int_red_frac.lt.0)then
! write (*,*) "int_red_frac",int_red_frac
! write (*,*) "divpar",divpar
! stop
!endif

return

end function int_red_frac


!####################################################
!# Calcul du plus grand diviseur commun a 3 entrees #
!####################################################

function pgdc3(a,b,c)

implicit none

integer(kind=DPI) :: pgdc3
integer(kind=DPI),intent(in) :: a,b,c

if(a.eq.0)then
pgdc3=Epgdc_dp(b,c)
else
if(b.eq.0)then
pgdc3=Epgdc_dp(a,c)
else
if(c.eq.0)then
pgdc3=Epgdc_dp(a,b)
else
pgdc3=Epgdc_dp(a,Epgdc_dp(b,c))
endif
endif
endif

return

end function pgdc3

!#############################################################################
!# Calcul du plus grand diviseur commun par decomposition de facteur premier #
!#############################################################################

function Ppgdc_dp(a,b)

implicit none

integer(DPI),intent (in)  :: a,b
integer(DPI) :: tab1(dimPREM),tab2(dimPREM),i,k,ii,kk,Ppgdc_dp

Ppgdc_dp=1

if(a*b.eq.0) then
Ppgdc_dp=max(abs(a),abs(b))
if(Ppgdc_dp.eq.0)then
write (*,*) "Uniquement des zero comme argument de pgdc_dp"
stop
endif
else
ii=dec_prem(a,tab1)
kk=dec_prem(b,tab2)
i=1
k=1

do while (i.le.ii.and.k.le.kk)
if(tab1(i).lt.tab2(k)) i=i+1
if(tab1(i).gt.tab2(k)) k=k+1
if(tab1(i).eq.tab2(k)) then
Ppgdc_dp=Ppgdc_dp*tab1(i)
i=i+1
k=k+1
endif
enddo
endif

return

end function Ppgdc_dp

!###########################################################################

function Ppgdc_dp_m(a,b)

implicit none

integer(DPIM),intent (in)  :: a,b
integer(DPIM) :: tab1(dimPREM),tab2(dimPREM),i,k,ii,kk,Ppgdc_dp_m

Ppgdc_dp_m=1

if(a*b.eq.0) then
Ppgdc_dp_m=max(abs(a),abs(b))
if(Ppgdc_dp_m.eq.0)then
write (*,*) "Uniquement des zero comme argument de pgdc_dp"
stop
endif
else
ii=dec_prem_m(a,tab1)
kk=dec_prem_m(b,tab2)
i=1
k=1

do while (i.le.ii.and.k.le.kk)
if(tab1(i).lt.tab2(k)) i=i+1
if(tab1(i).gt.tab2(k)) k=k+1
if(tab1(i).eq.tab2(k)) then
Ppgdc_dp_m=Ppgdc_dp_m*tab1(i)
i=i+1
k=k+1
endif
enddo
endif

return

end function Ppgdc_dp_m

!########################################
!# Calcul du plus petit multiple commun #
!########################################

function ppmc_dp_m(a,b)

implicit none

integer(DPIM),intent (in)  :: a,b
integer(DPIM) :: ppmc_dp_m,aa,bb

aa=abs(a)
bb=abs(b)

if(aa*bb.eq.0) stop "ppmc_dp : un zero en argument"
ppmc_dp_m=(aa/Epgdc_dp_m(aa,bb))*bb

return

end function ppmc_dp_m

!#################################################

function ppmc_dp(a,b)

implicit none

integer(DPI),intent (in)  :: a,b
integer(DPI) :: ppmc_dp,aa,bb

aa=abs(a)
bb=abs(b)

if(aa*bb.eq.0) stop "ppmc_dp : un zero en argument"
ppmc_dp=aa*bb/Epgdc_dp(aa,bb)

return

end function ppmc_dp

!###########################################################
!# renvoie le plus petit vecteur parallel au vecteur INput #
!###########################################################

function reduire_vec2D(vec)

implicit none

integer(kind=DPI) :: n
integer(kind=DPI),intent(in),dimension(2) :: vec
integer(kind=DPI),dimension(2) :: reduire_vec2D

n = Epgdc_dp(vec(1),vec(2))
reduire_vec2D (1:2)= vec(1:2) / n

return
end function reduire_vec2D

!###########################################################
!# renvoie le plus petit vecteur parallel au vecteur INput #
!###########################################################

function reduire_vec(vec)

implicit none

integer(kind=DPI) :: n
integer(kind=DPI),intent(in),dimension(3) :: vec
integer(kind=DPI),dimension(3) :: reduire_vec

if(vec(1).eq.0.and.vec(2).eq.0.and.vec(3).eq.0) then
n=1
else
n = pgdc3(vec(1),vec(2),vec(3))
endif
!if(n.eq.0)then
! reduire_vec (1:3)= vec(1:3)
! write (*,*) "reduire_vec PB pgcd",vec(1),vec(2),vec(3)
!else
 reduire_vec (1:3)= vec(1:3) / n
!endif

return
end function reduire_vec

!##########################################################################################
!# renvoie le type de plan avec des indices toujours positifs et dans le sens decroissant #
!##########################################################################################

function type_vec(vec)

implicit none

integer(kind=DPI) :: n,permuty
integer(kind=DPI),intent(in),dimension(3) :: vec
integer(kind=DPI),dimension(3) :: type_vec

if(vec(1).eq.0.and.vec(2).eq.0.and.vec(3).eq.0) then
n=1
else
n = pgdc3(vec(1),vec(2),vec(3))
endif

type_vec (1:3)= abs(vec(1:3)) / n

if(type_vec(3).gt.type_vec(2))then
permuty=type_vec(3)
type_vec(3)=type_vec(2)
type_vec(2)=permuty
endif
if(type_vec(2).gt.type_vec(1))then
permuty=type_vec(2)
type_vec(2)=type_vec(1)
type_vec(1)=permuty
endif
if(type_vec(3).gt.type_vec(2))then
permuty=type_vec(3)
type_vec(3)=type_vec(2)
type_vec(2)=permuty
endif

return
end function type_vec

!########################################
!# Calcul de la norme d'un vecteur reel #
!########################################

function entier(x)

implicit none

real(kind=DPR) :: x,z,preci,xx
integer(kind=DPI):: n
LOGICAL :: entier

z = dabs (x)
preci = z * Precision
n = int(z+PRECI)
xx = real(n)
if (z.eq.xx) then
entier = .TRUE.
else
entier = .FALSE.
endif
return
end function entier

!****************************************************************
!* Calcul le facteur multiplicatif pour etre sur le sous reseau *
!****************************************************************
!
! v=a*va+b*vb
!

function FacSr(iva,ivr)

implicit none

integer(kind=DPI) :: iva(3),ivr(3),FacSr

! real(kind=DPR) :: va(3),vr(3)
! real(kind=DPR)    :: fa
integer(kind=DPI) :: ia

!write (*,*) "FacSr(",iva,",",ivr,")"

!va=iva
!vr=ivr

! Pour se positionner sur le sous reseau

!fa=va(1)/vr(1)
!write (*,*) "fa",fa
!ia=intfrac(fa)
!write (*,*) "pgcd",Epgdc_dp(abs(iva(1)),abs(ivr(1)))
ia=abs(ivr(1))/Epgdc_dp(abs(iva(1)),abs(ivr(1)))
!write (*,*) "ia",ia
!fa=va(2)/vr(2)
!write (*,*) "fa",fa
!ia=ppmc_dp(ia,intfrac(fa))
ia=ppmc_dp(ia,abs(ivr(2))/Epgdc_dp(abs(iva(2)),abs(ivr(2))))
!write (*,*) "ia",ia
!fa=va(3)/vr(3)
!write (*,*) "fa",fa
!FacSr=ppmc_dp(ia,intfrac(fa))
FacSr=ppmc_dp(ia,abs(ivr(3))/Epgdc_dp(abs(iva(3)),abs(ivr(3))))
!write (*,*) "ia",FacSr

end function FacSr
