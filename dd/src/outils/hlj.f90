! Pour Reformater Histogramme des longueurs de jonctions
!===================== R. MADEC, le 03/04/01 : =========
program hlj

integer*8:: ios,InHist(1000,100),nbhisto,&
           reghis,regtra,OutHist(100,50),it,ih,buf
real :: tin(1000),tout(1000)
!=======================================================

open(49,FILE='in',STATUS='OLD',POSITION='REWIND')
open(50,FILE='out',STATUS='REPLACE')

ios=0
nbhisto=1

do while(ios.eq.0)

 read(49,55,iostat = ios) (InHist(nbhisto,JK),JK=1,100),tin(nbhisto)
 
 if (ios.eq.0) then
!  write (*,*) nbhisto,tin(nbhisto),(InHist(nbhisto,JK),JK=1,100)  
  nbhisto=nbhisto+1
 else
  nbhisto=nbhisto-1
 endif
 
enddo

write (*,*) nbhisto,' histogrammes de 100 tranches'
write (*,*) 'regroupement des histogrammes par : '
read (*,*) reghis
write (*,*) 'regroupement des tranches par : '
read (*,*) regtra

OutHist(:,:)=0
do ih=1,nbhisto,1
 do it=1,100,1
!  write (*,*) ih,it,((ih-1)/reghis)+1,((it-1)/regtra)+1,InHist(ih,it)
  buf=OutHist(((ih-1)/reghis)+1,((it-1)/regtra)+1)
  OutHist(((ih-1)/reghis)+1,((it-1)/regtra)+1)=&
   OutHist(((ih-1)/reghis)+1,((it-1)/regtra)+1)+InHist(ih,it)
   if(buf.gt.OutHist(((ih-1)/reghis)+1,((it-1)/regtra)+1)) then
    write (*,*) 'saturation interger*8'
    stop
   endif
 enddo
  tout((ih/reghis)+1)=tin(ih)
  write (*,*) ((ih-1)/reghis)+1,tout(((ih-1)/reghis)+1)
enddo

write (*,*) ((nbhisto-1)/reghis)+1,' histogrammes ',(99/regtra)+1,' tranches'

write(50,56) (JK,JK=0,100)

do it=1,(99/regtra)+1,1

 write(50,56) it,(OutHist(JK,it),JK=1,100)
 
enddo

write(50,*) 'X temps/histo -> dt=',tout(1) 
write(50,*) 'Y longeurs en Ljonc/Ldiagomax'

55 format(100(I4,1X),E12.5)
56 format(101(I8,1X))

close(49)
close(50)

!========================================================
end program hlj
