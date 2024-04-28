import film

cristallo, avalue , bveclin, bvecnor = film.film.readfilmheader()
print "cristallo", cristallo
print "avalue", avalue

slipsyscurv = 100
ldis_act = 0.5
curvact = 'F'
mode = 2
kk = 80

#for i in range(1000):

nsegm = film.film.readfilmstep(kk)

kneepoints, looppoints, slipsystem, loop,ntotpoints,ntotlines,nptsctr,seg_center,curv_vect = film.film.readfilmstepdata(mode, bveclin,bvecnor, nsegm,avalue,ldis_act, slipsyscurv,curvact)

print nptsctr
#raw_input("")

