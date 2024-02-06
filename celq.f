       program celq
c
c      pair potential md and/or cell shape md
c      written by renata m. wentzcovitch (iloop *?#!$)
c      current version: 4/26/91
c
c      This version uses the Beeman algorithim of integration
c
       implicit double precision (a-h,o-z)
c
       parameter (pi=3.141592653589793d0)
       parameter (mxdtyp = 1, mxdatm =512)
       parameter (mxdnat = mxdtyp * mxdatm)
       parameter (mxdnst =4000)
       parameter (mxdnr = 400)
c
       character*1 ic,iio
       character*2 calc
       character*2 nameat(mxdtyp)
c
       dimension avec(3,3),avecd(3,3),avec2d(3,3),avec2di(3,3)
       dimension avec0(3,3),sig0(3,3)
       dimension g(3,3),gm1(3,3),gd(3,3),sigma(3,3),pim(3,3)
       dimension rat(3,mxdatm,mxdtyp)
       dimension rataux(3,mxdatm,mxdtyp),nataux(mxdtyp)
       dimension ratd(3,mxdatm,mxdtyp)
       dimension rat2d(3,mxdatm,mxdtyp)
       dimension natom(mxdtyp),atmass(mxdtyp)
       dimension ncell(3)
       dimension car(3,mxdatm,mxdtyp)
       dimension v(3,mxdatm,mxdtyp),a(3,mxdatm,mxdtyp)
       dimension vint(3,mxdatm,mxdtyp),aint(3,mxdatm,mxdtyp)
       dimension rati(3,mxdatm,mxdtyp),rat2di(3,mxdatm,mxdtyp)
       dimension rand(3,mxdatm,mxdtyp)
       dimension carn(3,mxdnat),card(3,mxdnat)
       dimension inda(mxdnat),indt(mxdnat)
       dimension f(3,mxdnat),fint(3,mxdnat),frr(3,3)
       dimension vx2(mxdtyp),vy2(mxdtyp),vz2(mxdtyp)
       dimension rms(mxdtyp),vmean(mxdtyp),ekin(mxdtyp)
       dimension fs(3,mxdnat),fsint(3,mxdnat)
       dimension ran(3,mxdnat), rad(3,mxdnat)
       dimension theta(3,3),avmod(3),scp(3)
       dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
       dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
c
c      the following matrices are used for i/o purposes
c
       dimension utam(mxdnst),ekam(mxdnst),etam(mxdnst)
       dimension utlm(mxdnst),eklm(mxdnst),etlm(mxdnst)
       dimension utm(mxdnst),ekintm(mxdnst),etotm(mxdnst)
       dimension pm(mxdnst),pvm(mxdnst),tnewm(mxdnst),vcellm(mxdnst)
       dimension avkm(mxdnst),avum(mxdnst)
       dimension avpm(mxdnst),avpvm(mxdnst)
       dimension thetam(3,mxdnst),bmodm(3,mxdnst)
       dimension tt(3),il(mxdnst)
c
      data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
c
      factem = 1.57889d5
      boltz = um / factem
      ptmass = .5d0 * 1837.36
c
c      open files
c
        open(unit=5,file='inp',status='old',form='formatted')
        open(unit=6,file='out',status='new',form='formatted')
        open(unit=7,file='io')
        open(unit=20,file='car',status='new',form='formatted')
        open(unit=21,file='e',status='new',form='formatted')
        open(unit=22,file='eal',status='new',form='formatted')
        open(unit=23,file='ave',status='new',form='formatted')
        open(unit=24,file='p',status='new',form='formatted')
        open(unit=25,file='avec',status='new',form='formatted')
        open(unit=26,file='tv',status='new',form='formatted')
        open(unit=31,file='sip',status='old',form='formatted')
c       open(unit=32,file='ve',status='new',form='formatted')
c
       call tpage(calc)
c
       call crstl(alatt,avec,avecd,avec2di,vcell,g,gm1,calc,
     1 ic,iio,sigma,ntype,natom,nataux,nameat,atmass,
     2 rat,rataux,ratd,rat2di,car,v,cmass,avec0,sig0,v0,
     3 rcut,ncell,nstep,ntcheck,ntimes,dt,temp,ttol,press,
     4 mxdtyp,mxdatm,mxdnst,mxdnr,
     5 vpp,fpp,r0,dr,nr)
c
c     zero accumulators
c
      acp = zero
      acu = zero
      ack = zero
      acpv = zero
      avp = zero
      avu = zero
      avk = zero
      avpv = zero
c
      call init(alatt,vcell,
     1          ntype,natom,natot,atmass,cmass,press,temp,
     2          mxdtyp,mxdatm,mxdnat,mxdnr,
     3          rcut,ncell,nstep,calc,ic,
     4          rat,ratd,rat2d,rat2di,car,v,a,avec0,sig0,v0,
     5          avec,avecd,avec2d,avec2di,g,gd,gm1,sigma,pim,
     6          rand,vmean,rms,vx2,vy2,vz2,ekin,
     7          card,carn,rad,ran,vpp,fpp,r0,dr,nr,
     8          indt,inda,f,fint,fs,fsint,frr,vint,aint,
     9          uta,eka,eta,utl,ekl,etl,ut,ekint,etot)
c
c      start tapes 20 -- 25
c
       write(20,102) nstep
       write(21,102) nstep
       write(22,102) nstep
       write(23,102) nstep
       write(24,102) nstep
       write(25,102) nstep
       write(26,102) nstep
c      write(32,102) nstep
c
c      start md loop
c
      nzero = 0
c
      do 100 nst = 1,nstep
       call move(alatt,vcell,
     1           ntype,natom,natot,atmass,cmass,press,
     2           mxdtyp,mxdatm,mxdnat,mxdnr,
     3           rcut,ncell,nst,nstep,calc,dt,
     4           avec,avecd,avec2d,avec2di,avmod,theta,scp,
     5           avec0,sig0,v0,g,gd,gm1,sigma,pim,
     6           rat,ratd,rat2d,rati,rat2di,car,v,a,
     7           card,carn,ran,rad,vpp,fpp,r0,dr,nr,
     8           indt,inda,f,fint,fs,fsint,frr,vint,aint,
     9           uta,eka,eta,utl,ekl,etl,ut,p,ekint,etot)
c
c       update accumulators and averages
c
        nzero = nzero + 1
        acu = acu + ut
        ack = ack + ekint
        acp = acp + p
        avu = acu / dfloat(nzero)
        avk = ack / dfloat(nzero)
        avp = acp / dfloat(nzero)
c
c       calculate pv
c
        pv = p * vcell
        acpv = acpv + pv
        avpv = acpv / dfloat(nzero)
c
c       update output matrices
c
        utm(nst) = ut
        ekintm(nst) = ekint
        etotm(nst) = etot
        utam(nst) = uta
        ekam(nst) = eka
        etam(nst) = eta
        utlm(nst) = utl
        eklm(nst) = ekl
        etlm(nst) = etl
        pm(nst) = p
        pvm(nst) = pv
        avpm(nst) = avp
        avkm(nst) = avk
        avum(nst) = avu
        avpm(nst) = avp
        avpvm(nst) = avpv
c
c
        vcellm(nst) = vcell
        if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
        do 75 k=1,3
          bmodm(k,nst) = avmod(k)
 75     continue
        thetam(1,nst) = theta(1,2)
        thetam(2,nst) = theta(2,3)
        thetam(3,nst) = theta(3,1)
        endif
c
c       choose # of degrees of freedom and calculate tnew
c
        if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
        tnew = dois / tres / dfloat(natot + 1) * avk / boltz
        else
        tnew = dois / tres / dfloat(natot - 1) * avk / boltz
        endif
c
c       careful with zero temperature
c
        if (temp .lt. 1d-14) then
        ts = zero
        ttol = zero
        else
        ts =  dabs(tnew / temp - um)
        endif
c
        tnewm(nst) = tnew
c
c       rescale velocities
c
        if (mod(nst,ntcheck) .eq. zero) then
        if ((ts .gt. ttol) .and. (ntimes .gt. 0)) then
        acu = zero
        ack = zero
        acp = zero
        acpv = zero
        avu = zero
        avk = zero
        avp = zero
        avpv = zero
c
c
        if (tnew .le. 0.1d-12) then
        alpha = zero
        else
        alpha = dsqrt(temp / tnew)
        endif
        do 90 nt = 1,ntype
        do 90 na = 1,natom(nt)
        do 90 k=1,3
          v(k,na,nt) = alpha * v(k,na,nt)
          ratd(k,na,nt) = alpha * ratd(k,na,nt)
  90    continue
        if ((calc .eq. 'cd') .or. (calc. eq. 'nd')
     .                      .or. (calc. eq. 'sd')) then
        do 95 k = 1,3
        do 95 l = 1,3
          avecd(l,k) = alpha * avecd(l,k)
 95     continue
        endif
        ntimes = ntimes - 1
        nzero = 0
        endif
        endif
c
        if ((((calc .eq. 'mm') .or. (calc .eq. 'cm')) .or.
     .        (calc .eq. 'nm') .or. (calc .eq. 'sm'))) then
c       write(6,109) alpha,nst
        do 97 nt = 1,ntype
        do 97 na = 1,natom(nt)
        do 97 k=1,3
          xx = rat2di(k,na,nt)*rat2d(k,na,nt)
          if ( xx .lt. zero ) then
            alpha = zero
            v(k,na,nt) = alpha * v(k,na,nt)
            ratd(k,na,nt) = alpha * ratd(k,na,nt)
          endif
  97    continue
        if ((calc .eq. 'cm') .or. (calc .eq. 'nm') .or.
     .     (calc .eq. 'sm')) then
        do 98 k = 1,3
        do 98 l = 1,3
          xx = avec2d(l,k)*avec2di(l,k)
          if (xx .lt. zero) then
            alpha = zero
            avecd(l,k) = alpha * avecd(l,k)
          endif
 98     continue
        endif
        endif
c
c       optional output
c 
       write(21,101) utm(nst),ekintm(nst),
     ,                etotm(nst),pvm(nst),nst
       write(22,103) utam(nst),ekam(nst),etam(nst),
     ,                utlm(nst),eklm(nst),etlm(nst),nst
       write(23,104)  avum(nst),avkm(nst),nst
       write(24,105)  press,pm(nst),avpm(nst),nst
       if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
       write(25,103) (bmodm(k,nst),k=1,3),
     ,                (thetam(k,nst),k=1,3),nst
       endif
       write(26,104) vcellm(nst),tnewm(nst),nst
c      write(32,114) (avec(j,1),j=1,3),(avec(j,2),j=1,3),
c    ,               (avec(j,3),j=1,3)
 100   continue
c
c      end of md loop
c
c
c      i/o for next run
c
        rewind 7
        if (iio .eq. 'r') then
c        
c      rotate lattice vectors
c 
        write(7,113) avmod(1)
        tt(1) = um
        tt(2) = zero
        tt(3) = zero
        write(7,106) (tt(i),i=1,3)
        tt(1) = avmod(2) / avmod(1) * scp(1)
        tt(2) = avmod(2) / avmod(1) * dsqrt(um - scp(1)**2)
        tt(3) = zero
        write(7,106) (tt(i),i=1,3)
        tt(1) = scp(3)
        tt(2) = (scp(2) - scp(1)*scp(3)) / dsqrt(um - scp(1)**2) 
        tt(3) = dsqrt( um - tt(1) * tt(1) - tt(2) * tt(2))
        tt(1) = avmod(3) / avmod(1) * tt(1)
        tt(2) = avmod(3) / avmod(1) * tt(2)
        tt(3) = avmod(3) / avmod(1) * tt(3)
        write(7,106) (tt(i),i=1,3)
        write(7,106) zero,zero,zero
        write(7,106) zero,zero,zero
        write(7,106) zero,zero,zero
        write(7,106) zero,zero,zero
        write(7,106) zero,zero,zero
        write(7,106) zero,zero,zero
        else
c
c      do not rotate lattice vectors
c
        write(7,113) um
        write(7,106) (avec(i,1),i=1,3)
        write(7,106) (avec(i,2),i=1,3)
        write(7,106) (avec(i,3),i=1,3)
        write(7,106) (avecd(i,1),i=1,3)
        write(7,106) (avecd(i,1),i=1,3)
        write(7,106) (avecd(i,2),i=1,3)
        write(7,106) (avec2di(i,3),i=1,3)
        write(7,106) (avec2di(i,2),i=1,3)
        write(7,106) (avec2di(i,3),i=1,3)
        endif
c
        write(7,102) ntype
        do 110 nt = 1,ntype
        write(7,107) natom(nt),nameat(nt),atmass(nt)/ptmass
        do 110 na = 1,natom(nt)
        ra1 = rat(1,na,nt) - rat(1,1,1)
        ra2 = rat(2,na,nt) - rat(2,1,1)
        ra3 = rat(3,na,nt) - rat(3,1,1)
        write(7,106) ra1,ra2,ra3
        write(7,106) (ratd(i,na,nt),i=1,3)
        write(7,106) (rat2di(i,na,nt),i=1,3)
 110    continue
c
c       output for later analysis
c
       do 120 i=1,nstep
       il(i) = i
120    continue 

c      write(21,101) (utm(nst),ekintm(nst),
c    ,                etotm(nst),pvm(nst),il(nst),nst=1,nstep)
c      write(22,103) (utam(nst),ekam(nst),etam(nst),
c    ,                utlm(nst),eklm(nst),etlm(nst),il(nst),nst=1,nstep)
c      write(23,104) (avum(nst),avkm(nst),il(nst),nst=1,nstep)
c      write(24,105) (press,pm(nst),avpm(nst),il(nst),nst=1,nstep)
c      if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
c      write(25,103) ((bmodm(k,nst),k=1,3),
c    ,                (thetam(k,nst),k=1,3),il(nst),nst=1,nstep)
c      endif
c      write(26,104) (vcellm(nst),tnewm(nst),il(nst),nst=1,nstep)
c
1000    stop
 101   format(1x,4d12.5,i6)
 102   format(3i5)
 103   format(1x,6d12.5,i6)
 104   format(1x,2d12.5,i6)
 105   format(1x,3d12.5,i6)
 106   format(3(3x,d20.12))
 107   format(i5,3x,a2,f12.6)
 108   format(1x,'alpha = ',f7.4,' nstep = ',i4,
     , ' tnew = ',f11.4,'temp = ',f11.4,/)
 109   format(1x,'at quench alpha = ',f7.4,' nstep = ',i4,/)
 111   format(1x,'time used for cpu,io,sys,mem = ',4d8.4,/)
 113   format(f12.6)
 114   format(1x,9f8.4)
      end
c*
c*
       SUBROUTINE SPCEL(NSC,
     1                  NTYPE,NATOM,RAT,
     2                  RATAUX,NATAUX,
     3                  MXDTYP,MXDATM)
C
C      GENERATES THE ATOMIC BASIS OF A SUPERCELL GIVEN
C      THE ORIGINAL ATOMIC BASIS AND THE SUPERCELL DIMENSIONS
C            RMW 8/2/90
C
C      INPUT:
C      NSC(I)      SUPERCELL DIMENSION ALONG THE LATTICE VECTOR I
C      RAT         READ IN ATOMIC BASIS IN LATTICE COORDINATES
C      NTYPE       # OF ATOMIC SPECIES
C      NATOM(NT)   # OF ATOMS FOR THE ATOMIC SPECIES NT
C      MXDTYP      ARRAY DIMENSION FOR TYPE OF ATOMS
C      MXDATM      ARRAY DIMENSION FOR ATOMS OF A GIVEN TYPE
C
C      OUTPUT:
C      NATOM(I)    NEW NUMBER OF ATOMS OF TYPE I
C      RAT(K,N,I)  NEW K-TH COMPONENT (IN LATTICE COORDINATES) OF THE
C                  POSITION OF THE N-TH ATOM OF TYPE I
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
       PARAMETER ( ZERO = 0.0D0, UM = 1.0D0, TRES = 3.0D0)
C
       DIMENSION NATOM(MXDTYP),RAT(3,MXDATM,MXDTYP)
       DIMENSION NSC(3),ISC(3)
C
       DIMENSION NATAUX(MXDTYP),RATAUX(3,MXDATM,MXDTYP)
C
C      STORE TEMPORARY RAT'S AND NATOM'S
C
       DO 100 NT=1,NTYPE
       NATAUX(NT) = NATOM(NT)
       NATOM(NT) = 0
       DO 100 NA=1,NATAUX(NT)
       DO 100 K=1,3
       RATAUX(K,NA,NT) = RAT(K,NA,NT)
100    CONTINUE
C
C      GENERATES SUPERCELL ATOMIC BASIS
C
       ISC(1) = NSC(1) - 1
       ISC(2) = NSC(2) - 1
       ISC(3) = NSC(3) - 1
       DO 200 I1SC = 0,ISC(1)
       DO 200 I2SC = 0,ISC(2)
       DO 200 I3SC = 0,ISC(3)
       DO 200 NT=1,NTYPE
       DO 150 NA = 1,NATAUX(NT)
       NATOM(NT) = NATOM(NT) + 1
       NAN = NATOM(NT)
       RAT(1,NAN,NT) = RATAUX(1,NA,NT) / DFLOAT(NSC(1)) +
     +                 DFLOAT(I1SC) / DFLOAT(NSC(1))
       RAT(2,NAN,NT) = RATAUX(2,NA,NT) / DFLOAT(NSC(2)) +
     +                 DFLOAT(I2SC) / DFLOAT(NSC(2))
       RAT(3,NAN,NT) = RATAUX(3,NA,NT) / DFLOAT(NSC(3)) +
     +                 DFLOAT(I3SC) / DFLOAT(NSC(3))
 150   CONTINUE
C
       IF (NATOM(NT) .GT. MXDATM .OR. NATOM(NT) .LE. 0)
     1                CALL FATAL(52,XDUM,NATOM(NT))
C
 200   CONTINUE
       RETURN
       END
C*
C*
       subroutine tpage(calc)
c
c      finds version, time and date
c      written 2/5/90. rmw
c
       implicit double precision (a-h,o-z)
c
       character*2 calc
       character*4 vers
       character*40 title
c
       read(5,101) title
       write(6,102) title
       vers='celq'
       write(6,103) vers
c
c      types of calculations:
c
c      calc = md (molecular dynamics)
c      calc = cd (cell dynamics)
c      calc = nd (new cell dynamics)
c      calc = sd (strain dynamics)
c      calc = mm (minimization of energy wrt molecular configuration)
c      calc = cm (minimization of energy wrt cell configuration)
c      calc = nm (min. of energy wrt cell conf. by new dynamics)
c      calc = sm (min. of energy wrt cell conf. by strain dynamics)
c
       read(5,104) calc
       write(6,105) calc
       return
 101   format(a20)
 102   format(5x,a40//)
 103   format(/,5x,'program version = ',a4)
 104   format(a2)
 105   format(/,5x,'calculation type = ',a2)
       end
c*
c*
       subroutine crstl(alatt,avec,avecd,avec2di,vcell,g,gm1,calc,
     1 ic,iio,sigma,ntype,natom,nataux,nameat,atmass,
     2 rat,rataux,ratd,rat2di,car,v,cmass,avec0,sig0,v0,
     3 rcut,ncell,nstep,ntcheck,ntimes,dt,temp,ttol,press,
     4 mxdtyp,mxdatm,mxdnst,mxdnr,
     5 vpp,fpp,r0,dr,nr)
c
c      subroutine reads and computes various crystal parameters.
c      adapted from sverre froyen plane wave program
c      modified 2/5/90. rmw
c
c      input:
c      mxdtyp      array dimension for type of atoms
c      mxdatm      array dimension for atoms of a given type
c      mxdnst      array dimension for max. # of steps
c      mxdnr       array dimension for max. # of points in pair potential mesh
c
c      output:
c      alatt       lattice constant
c      avec(i,j)   i-th component of j-th direct lattice vector
c      avec0(i,j)  same as above but for 1st timestep only
c      avecd(i,j)  time derivative of lattice vector (read from unit 7)
c      g(i,j)      avec' x avec
c                   _1
c      gm1(i,j)    g
c      sigma(i,j)  {avec(2)*avec(3),avec(3)*avec(1),avec(1)*avec(2)}/vcell
c      sig0(i,j)    same as above but for 1st timestep only
c      vcell       cell volume
c      v0          cell volume (1st timestep)
c      ntype       number of types of atoms
c      natom(i)    number of atoms of type i
c      nameat(i)   chemical symbol for the type i
c      rat(k,n,i)  k-th component (in lattice coordinates) of the
c                  position of the n-th atom of type i
c      ratd(k,n,i)  time derivative of rat (read from unit 7)
c
c      ic    continuation id:
c               s = start new run
c               c = continues previous run (input from unit 7)
c                   but randonizes initial atomic velocities
c               f = follow the previous trajectory initally (read from io)
c
c      iio   mode to print file io
c                      r = rotate lattice vectors and zero lattice velocities
c                      n = print latt. vectors and latt. velocities "as are"
c      simulation parameters:
c                            rcut = pair potential cut-off radii
c                            ncell(3) = # of cells in each direction
c                            atmass(nt) = atomic masses in proton mass unit
c                            cmass = cell mass in "strange" units
c                            nstep = # of steps for the simulation
c                            ntcheck = check temperature after ntcheck steps
c                            ntimes =   "      "      ntimes times
c                            dt = time stepin a.u.
c                            temp = initial (mawellian) temperature in k
c                            ttol = relative temp. deviation tolerance
c                            press = pressure in a. u. units
c
c      *** pair potentials and forces are read and checked in rdpp
c          vpp = pair potential
c          fpp = force derived from pair potential
c          r0 = first point in the radial mesh
c          dr = distance between points in the radial grid
c          nr = # of points in the radial grid
c
       implicit double precision (a-h,o-z)
       parameter ( pi = 3.141592653589793d0, twopi = 2.0d0*pi)
       parameter ( zero = 0.0d0, um = 1.0d0, tres = 3.0d0)
       parameter ( eps = 1.0d-10)
c
       character*1 ic,iio
       character*1 iop(3,3),jop(3)
       character*2 nameat(mxdtyp),calc
       dimension avec(3,3),avecd(3,3),avec2di(3,3)
       dimension nsc(3),avec0(3,3),sig0(3,3)
       dimension g(3,3),gm1(3,3),sigma(3,3)
       dimension natom(mxdtyp),atmass(mxdtyp)
       dimension ncell(3),v(3,mxdatm,mxdtyp)
       dimension rat(3,mxdatm,mxdtyp),car(3,mxdatm,mxdtyp)
       dimension ratd(3,mxdatm,mxdtyp),test(3,3)
       dimension rat2di(3,mxdatm,mxdtyp)
       dimension rataux(3,mxdatm,mxdtyp),nataux(mxdtyp)
       dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
       dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
c
      ptmass = .5d0 * 1837.36
      a0 = .529177d0
      ry = 13.6058d0
      ev = 1.602d-12
c
c     pcv converts pressure from Mbar to Ry/(au)**3
c
      pcv = a0*a0*a0/ry/1.602d0
c
c      read continuation id
c
       read(5,106) ic,iio
       if (ic .eq. 's') then
c
c      read scale factor
c
       read(5,100) alatt
       write(6,200) alatt
c
c      read supercell dimensions
c
       read(5,104) (nsc(i),i=1,3)
       write(6,104) (nsc(i),i=1,3)
c
c      read the basis vectors for the lattice (cartesian coordinates and
c      atomic units divided by ascale). modify according to iop...
c      multiply by alatt.
c
       read(5,101) (iop(i,1),avec(i,1),i=1,3)
       read(5,101) (iop(i,2),avec(i,2),i=1,3)
       read(5,101) (iop(i,3),avec(i,3),i=1,3)
       do 9 j=1,3
       do 9 i=1,3
         sgn = um
         if (avec(i,j) .lt. zero) sgn = -sgn
         avec(i,j) = dabs(avec(i,j))
         if (iop(i,j) .eq. 's') avec(i,j) = sqrt(avec(i,j))
         if (iop(i,j) .eq. 'c') avec(i,j) = avec(i,j)**(um/tres)
         if (iop(i,j) .eq. 't') avec(i,j) = avec(i,j)/tres
         if (iop(i,j) .eq. 'h') avec(i,j) = avec(i,j)/sqrt(tres)
         avec(i,j) = sgn*alatt*avec(i,j)*nsc(j)
  9    continue
       else
c
c      read avec and avecd from unit 7
c
       read(5,112) 
       read(5,104) 
       read(7,105) alatt
       read(7,112) (avec(i,1), i = 1,3)
       read(7,112) (avec(i,2), i = 1,3)
       read(7,112) (avec(i,3), i = 1,3)
       write(6,112) (avec(i,1), i = 1,3)
       write(6,112) (avec(i,2), i = 1,3)
       write(6,112) (avec(i,3), i = 1,3)
       read(5,112) 
       read(5,112) 
       read(5,112) 
       read(7,112) (avecd(i,1), i = 1,3)
       read(7,112) (avecd(i,2), i = 1,3)
       read(7,112) (avecd(i,3), i = 1,3)
       write(6,112) (avecd(i,1), i = 1,3)
       write(6,112) (avecd(i,2), i = 1,3)
       write(6,112) (avecd(i,3), i = 1,3)
       read(7,112) (avec2di(i,1), i = 1,3)
       read(7,112) (avec2di(i,2), i = 1,3)
       read(7,112) (avec2di(i,3), i = 1,3)
       write(6,112) (avec2di(i,1), i = 1,3)
       write(6,112) (avec2di(i,2), i = 1,3)
       write(6,112) (avec2di(i,3), i = 1,3)
c
       nsc(1) = 1
       nsc(2) = 1
       nsc(3) = 1
c
c      rescale avec
c 
       do 10 j = 1,3
       do 10 i = 1,3
         avec(i,j) = avec(i,j) * alatt
 10    continue
       endif
c
       write(6,201)
       do 11 j=1,3
         write(6,202) j,(avec(i,j),i=1,3),(avec(i,j)/alatt,i=1,3)
 11    continue
c
c      sigma = lattice wave-vectors/twopi * cell volume
c
       sigma(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
       sigma(2,1) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
       sigma(3,1) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
       sigma(1,2) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
       sigma(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
       sigma(3,2) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
       sigma(1,3) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
       sigma(2,3) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
       sigma(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)
c
c      cell volume
c
       vcell = sigma(1,1)*avec(1,1) + sigma(2,1)*avec(2,1) +
     +         sigma(3,1)*avec(3,1)
       if(vcell .lt. eps) call fatal(50,vcell,idum)
       vcell = dabs(vcell)
c
c      calculate g and gm1 matrices
c
       do 22 j = 1,3
       do 22 i = 1,3
         g(i,j) = zero
         gm1(i,j) = zero
 22    continue
       do 26 j = 1,3
       do 26 i = 1,3
         do 25 m = 1,3
           g(i,j) = g(i,j) + avec(m,i) * avec(m,j)
           gm1(i,j) = gm1(i,j) + sigma(m,i) * sigma(m,j)
 25      continue
         gm1(i,j) = gm1(i,j) / vcell / vcell
 26    continue
c
c
c      printout
c
c
       write(6,205) vcell
c
c
c      read the number of atoms of each type and their names.
c        read atomic positions in unit cell (in units of
c        basis vectors). modify according to iop.
c
       read(5,105) cmass,press
       press = press * pcv
       write(6,212) cmass,press
       cmass = cmass * ptmass
c
       if (ic .eq. 's') then
c
c      read from regular input
c
       read(5,102) ntype
       if(ntype .gt. mxdtyp .or. ntype .le. 0)
     1         call fatal(51,xdum,ntype)
       do 40 nt=1,ntype
         read(5,103) natom(nt),nameat(nt),atmass(nt)
         atmass(nt) = atmass(nt) * ptmass
         if(natom(nt) .gt. mxdatm .or. natom(nt) .le. 0)
     1          call fatal(52,xdum,natom(nt))
         jmax = natom(nt)
         do 41 ja=1,jmax
           read(5,101) (jop(i),rat(i,ja,nt),i=1,3)
           do 42 i=1,3
             sgn = um
             if (rat(i,ja,nt) .lt. zero) sgn = -sgn
             rat(i,ja,nt) = dabs(rat(i,ja,nt))
             if (jop(i) .eq. 's') rat(i,ja,nt) = sqrt(rat(i,ja,nt))
             if (jop(i) .eq. 'c') rat(i,ja,nt) = rat(i,ja,nt)**(um/tres)
             if (jop(i) .eq. 't') rat(i,ja,nt) = rat(i,ja,nt)/tres
             rat(i,ja,nt) = sgn*rat(i,ja,nt)
 42        continue
 41      continue
 40    continue
c
       else
c
c      read output from last run
c
       read(5,102) 
       read(7,102) ntype
       write(6,102) ntype
       if(ntype .gt. mxdtyp .or. ntype .le. 0)
     1         call fatal(51,xdum,ntype)
       do 44 nt = 1,ntype
         read(5,103) 
         read(7,103) natom(nt),nameat(nt),atmass(nt)
         write(6,103) natom(nt),nameat(nt),atmass(nt)
         atmass(nt) = atmass(nt) * ptmass
         if(natom(nt) .gt. mxdatm .or. natom(nt) .le. 0)
     1          call fatal(52,xdum,natom(nt))
         do 44 na = 1,natom(nt)
           read(5,107) 
           read(7,107) (rat(i,na,nt),i=1,3)
           read(7,107) (ratd(i,na,nt),i=1,3)
           read(7,107) (rat2di(i,na,nt),i=1,3)
           write(6,107) (rat(i,na,nt),i=1,3)
           write(6,107) (ratd(i,na,nt),i=1,3)
           write(6,107) (rat2di(i,na,nt),i=1,3)
 44    continue
       endif
c
c      generates supercell
c
       call spcel(nsc,
     1             ntype,natom,rat,
     2             rataux,nataux,
     3             mxdtyp,mxdatm)
c
c      printout
c
       write(6,206)
       ntt = 0
       do 45 nt=1,ntype
         jmax = natom(nt)
         do 46 ja=1,jmax
           ntt = ntt + 1
           do 461 k=1,3
             car(k,ja,nt) = avec(k,1)*rat(1,ja,nt) +
     +                      avec(k,2)*rat(2,ja,nt) +
     +                      avec(k,3)*rat(3,ja,nt)
             v(k,ja,nt) = avec(k,1)*ratd(1,ja,nt) +
     +                    avec(k,2)*ratd(2,ja,nt) +
     +                    avec(k,3)*ratd(3,ja,nt)
  461      continue
           write(6,207) ntt,nameat(nt),atmass(nt)/ptmass,
     1                 (rat(i,ja,nt),i=1,3),
     2                 (car(i,ja,nt),i=1,3)
 46      continue
 45    continue
c
c     define quantities related to strain diynamics
c
      do 47 j = 1,3
      do 47 i = 1,3
        avec0(i,j) = avec(i,j)
        sig0(i,j) = sigma(i,j)
 47   continue
      v0 = vcell
c
c     simulation parameters
c
       read(5,100) rcut
       write(6,208) rcut
       read(5,104) (ncell(i),i=1,3)
       write(6,209) ((2 * ncell(i) + 1),i=1,3)
       read(5,102) nstep,ntcheck,ntimes
       write(6,210) nstep,ntcheck,ntimes
       read(5,105) temp,ttol,dt
       write(6,211) temp,ttol,2*dt
c
c
       if (nstep .gt. mxdnst) call fatal(53,xdum,nstep)
c
c      read pair potential
c
       ilj = 1
       if (ilj .eq. 0) then
       call rdpp(mxdtyp,mxdnr,ntype,vpp,fpp,r0,dr,nr)
       endif
c
c
       return
c
 100   format(f12.6)
 101   format(3(2x,a1,f12.6))
 102   format(3i5)
 103   format(i5,3x,a2,f12.6)
 104   format(3i5)
 105   format(1x,3f12.6)
 106   format(1x,a1,3x,a1)
 107   format(3(3x,d20.12))
 112   format(1x,3(3x,d20.12))
 200   format(//,'  crystal structure:',//,'  lattice constant ',
     1        f12.5,' (a.u.)',/)
 201   format(/,'  primitive translation vectors',/,25x,'in a.u.',
     1        20x,'in units of the lattice constant')
 202   format('  a',i1,'=',3(2x,e12.6),5x,3(2x,f7.3))
 204   format('  b',i1,'=',3(2x,e12.6),5x,3(2x,f7.3))
 205   format(/,'  cell volume = ',f10.2)
 206   format(/,'  no.   type     mass     position(lat. coord.)  ',
     1          '    position(cartesian coordinates)'/)
 207   format(3x,i3,3x,a2,3x,f8.2,3(3x,f7.3),5x,3(3x,e12.6))
 208   format(/,2x,'pair potential cut off = ',f10.5)
 209   format(/,2x,'no. of cells in each direction =  ',3i5)
 210   format(/,2x,'no. of steps = ', i5,'   , check t each',i5,
     ,        ' steps for',i5,' times')
 211   format(/,2x,'temperature = ',f9.2,' k','       ttol = ',
     ,        f5.3,'      time step = ',f10.4)
 212   format(/,2x,'cell mass = ',f15.10,'  pressure = ',f9.5,' a.u.')
       end
c*
c*
      subroutine init(alatt,vcell,
     1          ntype,natom,natot,atmass,cmass,press,temp,
     2          mxdtyp,mxdatm,mxdnat,mxdnr,
     3          rcut,ncell,nstep,calc,ic,
     4          rat,ratd,rat2d,rat2di,car,v,a,avec0,sig0,v0,
     5          avec,avecd,avec2d,avec2di,g,gd,gm1,sigma,pim,
     6          rand,vmean,rms,vx2,vy2,vz2,ekin,
     7          card,carn,rad,ran,vpp,fpp,r0,dr,nr,
     8          indt,inda,f,fint,fs,fsint,frr,vint,aint,
     9          uta,eka,eta,utl,ekl,etl,ut,ekint,etot)
c
c     subroutine written by rmw on 2/5/90   .               ..
c     initialize the dynamical variables q, q, and computes q
c     they are all given in cartesian coordinates
c
c     input:
c     mxdtyp = array dimension for type of atoms
c     mxdatm = array dimension for atoms of a given type
c     mxdnat = array dimension for total number of atoms
c     calc = calculation type
c     ic = continuation id
c     ntype = number of types of atoms
c     natom(nt) = number of atoms of type nt
c     atmass(nt) = atomic masses for atoms of type nt (in proton masses)
c     cmass = cell mass
c     avec(3,3) = lattice vectors
c     avec0(3,3) = initial lattice vectors
c                  t
c     g(3,3) = avec * avec
c                  _1
c     gm1(3,3) = g
c     sigma = reciprocal lattice vectors / (2 * pi) * vcell
c     sig0 = initial reciprocal lattice vectors / (2 * pi) * vcell
c     rat(j,na,nt) = atomic positions in lattice coordinates
c     car(j,na,nt) = atomic positions in cartesian coordinates
c
c     output:
c     rat(j,na,nt) = atomic positions in lattice coordinates
c     ratd(j,na,nt) = atomic velocities in lattice coordinates
c     rat2d(j,na,nt) = atomic accelerations in lattice coordinates
c     rat2di(j,na,nt) = intermediary atomic accelerations in lattice coordinates
c     car(j,na,nt) = atomic position after random displacements (cart. coord.)
c     v(i,na,nt) = initial velocity of atom na of type nt (cart. coord.)
c     a(i,na,nt) = acceleration of atom na of type nt (cart. coord.)
c     avec(3,3) = lattice vectors
c     avecd(3,3) = 1st lattice vectors derivatives
c     avec2d(3,3) = 2nd lattice vectors derivatives
c     avec2di(3,3) = intermediary 2nd lattice vectors derivatives
c                    t                      t
c     gd(3,3) = avecd * avec + avecd * avec
c     ut = total potential energy
c     ekint = total kinetic energy
c     etot = total energy
c     the atomic and lattice contributions uta,eka,eta,utl,ekl,etl
c
      implicit real*8 (a-h,o-z)
c
      dimension avec0(3,3),sig0(3,3)
      dimension avec(3,3),avecd(3,3),avec2d(3,3),avec2di(3,3)
      dimension g(3,3),gm1(3,3),gd(3,3),sigma(3,3),pim(3,3)
      dimension gmgd(3,3),sigav(3,3)
      dimension natom(mxdtyp),atmass(mxdtyp)
      dimension ncell(3)
      dimension rat(3,mxdatm,mxdtyp)
      dimension ratd(3,mxdatm,mxdtyp)
      dimension rat2d(3,mxdatm,mxdtyp),rat2di(3,mxdatm,mxdtyp)
      dimension car(3,mxdatm,mxdtyp),v(3,mxdatm,mxdtyp)
      dimension a(3,mxdatm,mxdtyp)
      dimension vint(3,mxdatm,mxdtyp),aint(3,mxdatm,mxdtyp)
      dimension rand(3,mxdatm,mxdtyp)
      dimension rms(mxdtyp),vmean(mxdtyp),ekin(mxdtyp)
      dimension vx2(mxdtyp),vy2(mxdtyp),vz2(mxdtyp)
      dimension carn(3,mxdnat), card(3,mxdnat)
      dimension inda(mxdnat),indt(mxdnat)
      dimension f(3,mxdnat),fint(3,mxdnat),frr(3,3)
      dimension fs(3,mxdnat),fsint(3,mxdnat)
      dimension ran(3,mxdnat),rad(3,mxdnat)
      dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
      dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
      dimension d2(3,3)
c
       character*1 ic
       character*2 calc
c
      data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
c
      ptmass = .5d0 * 1837.36
c
c     information for later plotting
c
      write(20,301) ntype
      write(20,301) (natom(nt),nt=1,ntype)
 301  format(1x,5i5)
c
c     establish maxwellian distribution of velocities if starting calculation
c
      if ((ic .eq. 's') .or. (ic .eq. 'c')) then
c
        call ranv(ntype,natom,atmass,
     1            mxdtyp,mxdatm,
     2            temp,ekint,
     3            v,rand,vmean,rms,vx2,vy2,vz2,ekin)
c
c       initialize avecd and avec2di
c
        ivzero = 1
        if (ivzero .eq. 1) then
          do 50 j=1,3
          do 50 i=1,3
            avecd(i,j) = zero
            avec2d(i,j) = zero
            avec2di(i,j) = zero
            gd(i,j) = zero
 50       continue
        else
c
c     initialize avecd at your own taste if not read from unit 7
c          .
c          .
c          .
        endif
      endif
      do 25 j = 1,3
      do 25 i = 1,3
        gd(i,j) = zero
 25   continue
      do 30 k=1,3
      do 30 j = 1,3
      do 30 i = 1,3
        gd(i,j) = gd(i,j) + avec(k,i) * avecd(k,j) +
     +            avecd(k,i) * avec(k,j)
 30   continue
c
c     update atomic positions using the randomly generated #'s
c
c
c     warning!!!!!!!!! the following initialization has to
c     be carefully planned. the development of the simulation
c     depends strongly on the initialization.
c     the following lines should not be used "as are"!!
c
      iran = 0
      if (iran .eq. 1) then
      do 240 nt = 1,ntype
        do 230 na = 1,natom(nt)
          natot = natot + 1
          rat(1,na,nt) =  rat(1,na,nt) + rand(1,na,nt) / 100.d0
          rat(2,na,nt) =  rat(2,na,nt) + rand(2,na,nt) / 100.d0
          rat(3,na,nt) =  rat(3,na,nt) + rand(3,na,nt) / 100.d0
 230    continue
 240  continue
      endif
c
c     transform position to cartesian coordinates and velocities 
c     to lattice coordinates
c
      natot = 0
      do 700 nt=1,ntype
        natot = natot + natom(nt)
        do 700 na=1,natom(nt)
        do 250 k=1,3
          car(k,na,nt) = zero
          ratd(k,na,nt) = zero
          rat2d(k,na,nt) = zero
 250    continue
        do 500 l=1,3
          do 500 k=1,3
            car(l,na,nt) = car(l,na,nt) + avec(l,k) * rat(k,na,nt)
            ratd(l,na,nt) = v(k,na,nt) * sigma(k,l) / vcell +
     +                      ratd(l,na,nt)
 500    continue
 700  continue
c
c     calculate forces and accelerations in cart. and latt. coord.'s
c
      ilj = 1
      if (ilj .eq. 1) then
      call forclj(avec,sigma,vcell,calc,
     1          ntype,natom,natot,atmass,
     2          mxdtyp,mxdatm,mxdnat,
     3          rcut,ncell,
     4          car,v,a,rat,ratd,rat2d,
     5          card,carn,rad,ran,
     6          indt,inda,f,fs,fint,fsint,frr,ut,w)
      else
      call forc(avec,sigma,vcell,calc,
     1          ntype,natom,natot,atmass,
     2          mxdtyp,mxdatm,mxdnat,mxdnr,
     3          rcut,ncell,
     4          car,v,a,rat,ratd,rat2d,
     5          card,carn,rad,ran,
     6          indt,inda,f,fs,fint,fsint,frr,ut,w,
     7          vpp,fpp,r0,dr,nr)
      endif
c
c      initialize rat2di
c
      do 710 nt = 1,ntype
      do 710 na = 1,natom(nt)
        rat2di(1,na,nt) = zero
        rat2di(2,na,nt) = zero
        rat2di(3,na,nt) = zero
 710  continue
c
c
       if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
c                _1  .
c     calculate g * g    ( = gmgd)
c
      do 760 j = 1,3
      do 760 i = 1,3
        gmgd(i,j) = zero
        pim(i,j) = zero
        do 750 m = 1,3
          gmgd(i,j) = gmgd(i,j) + gm1(i,m) * gd(m,j)
 750    continue
 760  continue
c
c     change forces on particles
c
      do 830 nt = 1,ntype
      do 830 na = 1,natom(nt)
        do 800 k = 1,3
          do 800 m = 1,3
            rat2d(k,na,nt) = rat2d(k,na,nt) -
     -                       gmgd(k,m) * ratd(m,na,nt)
 800      continue
c
c     initialize forces on basis vectors
c
        do 820 j = 1,3
          do 820 i = 1,3
          pim(i,j) = pim(i,j) + atmass(nt) * v(i,na,nt) *
     *               v(j,na,nt)
 820    continue
 830  continue
      do 840 j = 1,3
        do 840 i = 1,3
        pim(i,j) = (pim(i,j) + frr(i,j)) / vcell
c
c       diagonal
c
        if (i .eq. j) then
          pim(i,j) = pim(i,j) - press
        endif
        avec2d(i,j) = zero
 840  continue
c
      do 860 j = 1,3
        do 860 i = 1,3
          do 850 k = 1,3
          avec2d(i,j) = avec2d(i,j) + pim(i,k) *
     *                  sigma(k,j)
 850      continue
        avec2d(i,j) = avec2d(i,j) / cmass
 860  continue
c
c     if new cell dynamics...
c
      if ((calc .eq. 'nd') .or. (calc .eq. 'nm')) then
        call sigp(avec,avecd,avec2d,sigma,vcell)
      endif
c
      if ((calc .eq. 'sd') .or. (calc .eq. 'sm')) then
        call sigs(avec0,avec2d,sig0,v0)
      endif
c
c     strain symmetrization
c
      do 862 i=1,3
      do 862 j=1,3
        d2(i,j) = zero
        do 861 k=1,3
          d2(i,j) = d2(i,j) + avec2d(i,k) * sig0(j,k)
 861    continue
        d2(i,j) = d2(i,j) / v0
 862  continue
c
      do 863 i=1,2
      do 863 j=i+1,3
      d2(i,j) = (d2(i,j) + d2(j,i)) / dois
      d2(j,i) = d2(i,j)
 863  continue
c
c     write(6,*) d2(2,1),d2(3,1),d2(3,2)
c     write(6,*) d2(1,2),d2(1,3),d2(2,3)
c
      do 865 i=1,3
      do 865 j=1,3
        avec2d(i,j) = zero
        do 864 k=1,3
          avec2d(i,j) = avec2d(i,j) + d2(i,k) * avec0(k,j)
 864    continue
 865    continue
c     write(6,*) avec2d(2,1),avec2d(3,1),avec2d(3,2)
c
      endif
c
c     computes initial total energy
c
      uta = ut
      eka = ekint
      eta = uta + eka
      etot = ut + ekint
c
c     lattice contribution
c
      ekl = zero
c
      if ((calc .ne. 'md') .and. (calc .ne. 'mm')) then
        if ((calc .eq. 'nd') .or. (calc .eq. 'nm')) then
        do 890 j = 1,3
          do 890 i = 1,3
            sigav(i,j) = zero
            do 890 l = 1,3
              sigav(i,j) = sigav(i,j) + sigma(l,i) * avecd(l,j)
 890    continue
        do 950 k = 1,3
          tr = zero
          do 900 m = 1,3
            tr = tr + sigav(m,k) * sigav(m,k)
 900      continue
          ekl = ekl + tr
 950    continue
c
      endif
c
        if ((calc .eq. 'sd') .or. (calc .eq. 'sm')) then
        do 960 j = 1,3
          do 960 i = 1,3
            sigav(i,j) = zero
            do 960 l = 1,3
              sigav(i,j) = sigav(i,j) + sig0(l,i) * avecd(l,j)
 960    continue
        do 980 k = 1,3
          tr = zero
          do 970 m = 1,3
            tr = tr + sigav(m,k) * sigav(m,k)
 970      continue
          ekl = ekl + tr
 980    continue
c
        else
c
        do 995 k = 1,3
          tr = zero
          do 990 m = 1,3
            tr = tr + avecd(m,k) * avecd(m,k)
 990      continue
          ekl = ekl + tr
 995    continue
      endif
c
      ekl = ekl * cmass / dois
      utl = press * vcell
      ekint = eka + ekl
      ut = uta + utl
      etot = ekint + ut
      endif
c
      return
      end
c*
c*
      subroutine ranv(ntype,natom,atmass,
     1          mxdtyp,mxdatm,
     2          temp,ekint,
     3          v,rand,vmean,rms,vx2,vy2,vz2,ekin)
c
c     sets up random velocities with maxwellian distribution
c     at temperature t. total linear momentum components are zero
c     rewritten on 1/31/90 by rmw
c
c     input:
c     mxdtyp = array dimension for type of atoms
c     mxdatm = array dimension for atoms of a given type
c     mxdnat = array dimension for total number of atoms
c     ntype = number of types of atoms
c     natom(i) = number of atoms of type i
c     atmass(i) = atomic masses for atoms of type i (in proton masses)
c     temp = temperature in k
c
c     output:
c     v(i,na,nt) = initial velocity of atom na of type nt
c     rand(1,na,nt) = randon # to be used to initialize atomic positions
c     vmean(nt), rms(nt),vx2(nt),vy2(nt),vz2(nt)
c
      implicit real*8 (a-h,o-z)
c
      dimension atmass(mxdtyp),natom(mxdtyp)
      dimension v(3,mxdatm,mxdtyp), rand(3,mxdatm,mxdtyp), p(3)
      dimension vx2(mxdtyp), vy2(mxdtyp), vz2(mxdtyp)
      dimension rms(mxdtyp), vmean(mxdtyp), ekin(mxdtyp)
c
      data b0, b1, c0, c1 / 2.30753d0, 0.27061d0, 0.99229d0, 0.04481d0 /
      data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
c
c     boltz = boltzman constant in ry/kelvin
c     ptmass = proton (neuton) mass in a. u.
c
      factem = 1.57889d5
      boltz = um  / factem
      ptmass = .5d0 * 1837.36
c
c     example run
c
      ntot = 0
      do 50 nt = 1, ntype
        ntot = ntot + natom(nt)
        ekin(nt) = zero
 50   continue
      ekint = zero
c
c
      if (ntot .ne. 1) then
c
c     assign random velocities
c
      t = temp
      if (temp .lt. 1.d-14) t = 1.d-14
       iseed = - 119
       eps = ran3(iseed)
c
c     establish gaussian distribution for each atom kind
c
      do 900 nt = 1,ntype
        vfac = dsqrt(boltz * t / atmass(nt))
        write(6,901)
        write(6,*) 'vfac = ',vfac
        iseed = iseed + 382
        do 200 na = 1,natom(nt)
          do 200 j = 1,3
            eps = ran3(iseed)
            rand(j,na,nt) = eps
            if (eps .lt. 1.d-10) eps = 1.d-10
            if (eps .le. 0.5d0) go to 100
            eps = eps - um
            if (eps .gt. -1.d-10) eps = -1.d-10
100         sig = dsqrt( log(um / (eps * eps)))
            vr = sig - (b0 + b1 * sig) /
     /           (um + c0 * sig + c1 * sig * sig)
            vr = vr * vfac
            if (eps .lt. zero) vr = - vr
            v(j,na,nt) = vr
 200    continue
c
        p(1) = zero
        p(2) = zero
        p(3) = zero
c
c       calculate linear-momentum.
c
        do 400 na=1,natom(nt)
          p(1) = p(1) + v(1,na,nt)
          p(2) = p(2) + v(2,na,nt)
          p(3) = p(3) + v(3,na,nt)
 400    continue
        p(1) = p(1) / dfloat(natom(nt))
        p(2) = p(2) / dfloat(natom(nt))
        p(3) = p(3) / dfloat(natom(nt))
c
c       zero linear momentum for atom type nt
c
        do 600 na = 1,natom(nt)
          v(1,na,nt) = v(1,na,nt) - p(1)
          v(2,na,nt) = v(2,na,nt) - p(2)
          v(3,na,nt) = v(3,na,nt) - p(3)
 600    continue
        do 800 na = 1,natom(nt)
          ekin(nt) = ekin(nt) + (v(1,na,nt) * v(1,na,nt)
     +                        + v(2,na,nt) * v(2,na,nt)
     +                        + v(3,na,nt) * v(3,na,nt)) / dois
 800    continue
        write(6,*) 'ekin(nt)',ekin(nt)
        ekin(nt) = atmass(nt) * ekin(nt)
        ekint = ekint + ekin(nt)
 900  continue
c
c     rescale velocities to give correct temperature
c
      atemp = dois * ekint /
     /            tres / dfloat(ntot - 1) / boltz
      tfac = dsqrt(t / atemp)
      if (temp .lt. 1d-14) tfac = zero
      write(6,*) 'atemp = ',atemp,' k'
      write(6,*) 'tfac = ',tfac
      do 1200 nt = 1,ntype
        vmean(nt) = zero
        rms(nt) = zero
        vx2(nt) = zero
        vy2(nt) = zero
        vz2(nt) = zero
        do 1000 na = 1,natom(nt)
          v(1,na,nt) = v(1,na,nt) * tfac
          v(2,na,nt) = v(2,na,nt) * tfac
          v(3,na,nt) = v(3,na,nt) * tfac
          vmean(nt) = vmean(nt) +
     +      dsqrt(v(1,na,nt) ** 2 + v(2,na,nt) ** 2 + v(3,na,nt) **2)
          vx2(nt) = vx2(nt) + v(1,na,nt) ** 2
          vy2(nt) = vy2(nt) + v(2,na,nt) ** 2
          vz2(nt) = vz2(nt) + v(3,na,nt) ** 2
1000  continue
        vmean(nt) = vmean(nt) / dfloat(natom(nt))
        rms(nt) = dsqrt((vx2(nt) + vy2(nt) + vz2(nt)) /
     /            dfloat(natom(nt)))
        vx2(nt) = dsqrt(vx2(nt) / dfloat(natom(nt)))
        vy2(nt) = dsqrt(vy2(nt) / dfloat(natom(nt)))
        vz2(nt) = dsqrt(vz2(nt) / dfloat(natom(nt)))
1200  continue
      ekint = ekint * tfac * tfac
      else
      ekint = zero
      do 1300 k = 1,3
        rand(k,1,1) = zero
        v(k,1,1) = zero
 1300 continue
      vmean(1) = zero
      rms(1) = zero
      vx2(1) = zero
      vy2(1) = zero
      ekin(1) = zero
      endif
      return
 801    format(1x,5f14.10)
 901    format(/,10x, 'initial conditions',/)
1999    format(1x,//)
      end
c*
c*
      subroutine forclj(avec,sigma,vcell,calc,
     1          ntype,natom,natot,atmass,
     2          mxdtyp,mxdatm,mxdnat,
     3          rcut,ncell,
     4          car,v,a,rat,ratd,rat2d,
     5          card,carn,rad,ran,
     6          indt,inda,f,fs,fint,fsint,frr,ut,w)
c
c     subroutine written by rmw on 2/5/90
c     calculates total energy and force for a number of
c     particles arranged with periodic boundary condition
c     lennard-jones potential (1 atom type) assumed temporarily
c
c
c     input:
c     mxdtyp = array dimension for type of atoms
c     mxdatm = array dimension for atoms of a given type
c     mxdnat = array dimension for total number of atoms
c     ntype = number of types of atoms
c     natom(nt) = number of atoms of type nt
c     atmass(nt) = atomic masses for atoms of type nt (in proton masses)
c     car(j,na,nt) = atomic positions in cartesian coordinates
c     v(j,na,nt) = atomic velocities
c     ncell(k) = # of periodicly reapeating cells along the 3
c                  primitive cell direction
c     avec(3,3) = lattice vectors
c     rcut = cutoff radius for force calculation
c
c     output:
c     a(i,na,nt) = acceleration on atom nat
c     rat2d(i,na,nt) = acceleration on atom nat in latt. coord.
c     ut = total potential energy
c     w = 3 * potential contribution to the virial pressure
c     frr(3,3) = "         "   to the cell forces
c     f(3) = total force in cartesian coordinated
c     fs(3) = "     "     "   lattice  "
c
      implicit real*8 (a-h,o-z)
c
      dimension avec(3,3),sigma(3,3)
      dimension natom(mxdtyp),atmass(mxdtyp)
      dimension car(3,mxdatm,mxdtyp),v(3,mxdatm,mxdtyp)
      dimension a(3,mxdatm,mxdtyp)
      dimension rat(3,mxdatm,mxdtyp),ratd(3,mxdatm,mxdtyp)
      dimension rat2d(3,mxdatm,mxdtyp)
      dimension inda(mxdnat),indt(mxdnat)
      dimension f(3,mxdnat),fint(3,mxdnat),frr(3,3)
      dimension rd(3),fo(3),fr(3,3)
      dimension ncell(3)
      dimension carn(3,mxdnat), card(3,mxdnat)
      dimension sd(3),fs(3,mxdnat),fos(3),fsint(3,mxdnat)
      dimension ran(3,mxdnat), rad(3,mxdnat)
c
       character*2 calc
c
      data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
c
      ptmass = .5d0 * 1837.36
c
      rcutsq = rcut ** 2
c
c     # of cells
c
      n1 = ncell(1)
      n2 = ncell(2)
      n3 = ncell(3)
c
c     construct indices arrays
c
      nato = 0
      do 100 nt = 1,ntype
      do 100 na = 1,natom(nt)
      nato = nato + 1
      do 90 k=1,3
      carn(k,nato) =  car(k,na,nt)
      ran(k,nato) =  rat(k,na,nt)
 90   continue
      inda(nato) = na
      indt(nato) = nt
 100  continue
c
c
       natot = nato
c
c
c     the lenard-jones calculation works for only one atom type
c
      ut = zero
      w = zero
      do 200 k = 1,3
      do 200 l = 1,3
      frr(l,k) = zero
 200  continue
      do 250 nat = 1,natot
      do 250 k = 1,3
      f(k,nat) = zero
      fint(k,nat) = zero
      fs(k,nat) = zero
      fsint(k,nat) = zero
 250  continue
c
c     start loop over particles
c
      do 1000 i1 = -n1,n1
      do 1000 i2 = -n2,n2
      do 1000 i3 = -n3,n3
      if (((i1 .eq. 0) .and. (i2 .eq. 0)).and. (i3 .eq. 0)) then
c
c     particles i and j in the same cell
c
      do 400 i = 1, natot - 1
      do 400 j = i + 1, natot
      do 300 k=1,3
      rd(k) = carn(k,j) - carn(k,i)
      sd(k) = ran(k,j) - ran(k,i)
 300  continue
      rdsq = rd(1) ** 2 +rd(2) ** 2 + rd(3) ** 2
c
      if (rdsq  .le. rcutsq) then
      call elj(rdsq,rd,sd,ulj,fo,fos,fr,wlj,calc)
      w = w + wlj
      ut = ut + ulj
      do 340 jj = 1,3
      do 340 ii = 1,3
      frr(ii,jj) = frr(ii,jj) + fr(ii,jj)
 340  continue
      do 350 k=1,3
      f(k,i) = f(k,i) - fo(k)
      f(k,j) = f(k,j) + fo(k)
      fint(k,i) = fint(k,i) - fo(k)
      fint(k,j) = fint(k,j) + fo(k)
      fs(k,i) = fs(k,i) - fos(k)
      fs(k,j) = fs(k,j) + fos(k)
      fsint(k,i) = fsint(k,i) - fos(k)
      fsint(k,j) = fsint(k,j) + fos(k)
 350  continue
      endif
 400  continue
      else
c
c     particles i and j in different cells
c
      do 700 nat = 1,natot
      do 600 i = 1,3
      rad(i,nat) = zero
      card(i,nat) = carn(i,nat) + dfloat(i1) * avec(i,1) +
     +           dfloat(i2) * avec(i,2) +
     +           dfloat(i3) * avec(i,3)
 600  continue
      do 650 k=1,3
      do 650 l=1,3
      rad(k,nat) = rad(k,nat) + sigma(l,k) / vcell * card(l,nat)
 650  continue
 700  continue
      do 900 i = 1,natot
      do 900 j = 1,natot
      do 750 k=1,3
      rd(k) = card(k,j) - carn(k,i)
      sd(k) = rad(k,j) - ran(k,i)
 750  continue
      rdsq = rd(1) ** 2 + rd(2) ** 2 + rd(3) ** 2
      if (rdsq .le. rcutsq) then
      call elj(rdsq,rd,sd,ulj,fo,fos,fr,wlj,calc)
      w = w + wlj / dois
      ut = ut + ulj / dois
      do 760 jj = 1,3
      do 760 ii = 1,3
      frr(ii,jj) = frr(ii,jj) + fr(ii,jj) / dois
 760  continue
      do 800 k=1,3
      f(k,i) = f(k,i) - fo(k)
      fs(k,i) = fs(k,i) - fos(k)
 800  continue
      endif
 900  continue
      endif
1000  continue
c
c     calculate forces w.r.t center of mass (just to check)
c
      fx = zero
      fy = zero
      fz = zero
      do 1300 nat=1,natot
      fx = fx + f(1,nat)
      fy = fy + f(2,nat)
      fz = fz + f(3,nat)
1300  continue
c
c     calculate acceleration
c
      do 1400 nat = 1,natot
      do 1400 k=1,3
      a(k,inda(nat),indt(nat)) = f(k,nat) / atmass(indt(nat))
      rat2d(k,inda(nat),indt(nat)) = fs(k,nat) / atmass(indt(nat))
1400  continue
      return
      end
c*
c*
      subroutine elj(rdsq,rd,sd,ulj,fo,fos,fr,wlj,calc)
c
c     computes the potential energy and force between 2 particles
c     interacting through a lenard-jones potential separated by
c     distance r
c
c     input:
c     rdsq = square of distance between the particles
c     rd(3) = relative position of particles in cartesian coordinates
c     sig = parameter for lj potential in a.u.
c     epsilon = parameter for lj potential in a.u.
c
c     output:
c     ulj = pair potential energy
c     wlj = pair contribution to the virial pressure
c     fr(3,3) = "     "       to the cell forces
c     fo(3) = force in cartesian coordinated
c     fos(3) = "     "   lattice  "

      implicit real*8 (a-h,o-z)
c                                                                              x
      data zero, um, dois, tres, quatro /0.d0, 1.d0, 2.d0, 3.d0, 4.d0/
c
      dimension rd(3), fo(3),fr(3,3)
      dimension sd(3), fos(3)
c
       character*2 calc
c
      a0 = .529177d0
      ry = 13.6058d0
      factem = 1.57889d5
      boltz = um / factem
c
c     parameters for he
c
c      sig =  2.57d0 / a0
c      epsilon = 10.8d0 / factem
c
c     parameters for ar (ashcroft pg 398)
c
      sig =  3.4 / a0
      epsilon = 1.04d-2 / ry
c
c     parameters for ne (ashcroft pg 398)
c
c      sig =  2.74 / a0
c      epsilon = .31d-2 / ry
c
c     potential energy
c
      sigsq = sig ** 2
      rrsqi = sigsq / rdsq
      ri6 = rrsqi * rrsqi * rrsqi
      ri12 = ri6 * ri6
      ulj = quatro * epsilon * (ri12 - ri6)
c
c     force
c
      factor = 24.d0 * epsilon / sigsq
      ri8 = ri6 * rrsqi
      ri14 = ri12 * rrsqi
      factor = factor * (dois * ri14 - ri8)
      wlj = zero
      do 100 k=1,3
      fo(k) =  factor * rd(k)
      wlj = wlj + fo(k) * rd(k)
      fos(k) = factor * sd(k)
 100  continue
c
c     for cell shape md
c
       if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
      do 300 j = 1,3
      do 200 i = 1,3
      fr(i,j) = fo(i) * rd(j)
 200  continue
 300  continue
      endif
      return
      end
c*
c*
      subroutine forc(avec,sigma,vcell,calc,
     1          ntype,natom,natot,atmass,
     2          mxdtyp,mxdatm,mxdnat,mxdnr,
     3          rcut,ncell,
     4          car,v,a,rat,ratd,rat2d,
     5          card,carn,rad,ran,
     6          indt,inda,f,fs,fint,fsint,frr,ut,w,
     7          vpp,fpp,r0,dr,nr)
c
c     subroutine written by rmw on 5/3/90
c     calculates total energy and force for a number of
c     particles arranged with periodic boundary condition
c
c
c     input:
c     mxdtyp = array dimension for type of atoms
c     mxdatm = array dimension for atoms of a given type
c     mxdnat = array dimension for total number of atoms
c     mxdnr = array dimension for # of points in the potential mesh
c     ntype = number of types of atoms
c     natom(nt) = number of atoms of type nt
c     atmass(nt) = atomic masses for atoms of type nt (in proton masses)
c     car(j,na,nt) = atomic positions in cartesian coordinates
c     v(j,na,nt) = atomic velocities
c     ncell(k) = # of periodicly reapeating cells along the 3
c                  primitive cell direction
c     avec(3,3) = lattice vectors
c     rcut = cutoff radius for force calculation
c     v = interacting potential
c     f = force
c     r0 = 1st point in the potential mesh
c     dr = potential mesh
c
c     output:
c     a(i,na,nt) = acceleration on atom nat
c     rat2d(i,na,nt) = acceleration on atom nat in latt. coord.
c     ut = total potential energy
c     w = 3 * potential contribution to the virial pressure
c     frr(3,3) = "         "   to the cell forces
c     f(3) = total force in cartesian coordinated
c     fs(3) = "     "     "   lattice  "
c
      implicit real*8 (a-h,o-z)
c
      dimension avec(3,3),sigma(3,3)
      dimension natom(mxdtyp),atmass(mxdtyp)
      dimension car(3,mxdatm,mxdtyp),v(3,mxdatm,mxdtyp)
      dimension a(3,mxdatm,mxdtyp)
      dimension rat(3,mxdatm,mxdtyp),ratd(3,mxdatm,mxdtyp)
      dimension rat2d(3,mxdatm,mxdtyp)
      dimension inda(mxdnat),indt(mxdnat)
      dimension f(3,mxdnat),fint(3,mxdnat),frr(3,3)
      dimension rd(3),fo(3),fr(3,3)
      dimension ncell(3)
      dimension carn(3,mxdnat), card(3,mxdnat)
      dimension sd(3),fs(3,mxdnat),fos(3),fsint(3,mxdnat)
      dimension ran(3,mxdnat), rad(3,mxdnat)
      dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
      dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
c
       character*2 calc
c
      data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
c
      ptmass = .5d0 * 1837.36
c
      rcutsq = rcut ** 2
c
c     # of cells
c
      n1 = ncell(1)
      n2 = ncell(2)
      n3 = ncell(3)
c
c     construct indices arrays
c
      nato = 0
      do 100 nt = 1,ntype
      do 100 na = 1,natom(nt)
      nato = nato + 1
      do 90 k=1,3
      carn(k,nato) =  car(k,na,nt)
      ran(k,nato) =  rat(k,na,nt)
 90   continue
      inda(nato) = na
      indt(nato) = nt
 100  continue
c
c
       natot = nato
c
      ut = zero
      w = zero
      do 200 k = 1,3
      do 200 l = 1,3
      frr(l,k) = zero
 200  continue
      do 250 nat = 1,natot
      do 250 k = 1,3
      f(k,nat) = zero
      fint(k,nat) = zero
      fs(k,nat) = zero
      fsint(k,nat) = zero
 250  continue
c
c     start loop over particles
c
      do 1000 i1 = -n1,n1
      do 1000 i2 = -n2,n2
      do 1000 i3 = -n3,n3
      if (((i1 .eq. 0) .and. (i2 .eq. 0)).and. (i3 .eq. 0)) then
c
c     particles i and j in the same cell
c
      do 400 i = 1, natot - 1
      do 400 j = i + 1, natot
      do 300 k=1,3
      rd(k) = carn(k,j) - carn(k,i)
      sd(k) = ran(k,j) - ran(k,i)
 300  continue
      rdsq = rd(1) ** 2 +rd(2) ** 2 + rd(3) ** 2
c
      if (rdsq  .le. rcutsq) then
      nt1 = indt(i)
      nt2 = indt(j)
      call vpair(rdsq,rd,sd,upp,fo,fos,fr,wpp,
     1                 calc,nt1,nt2,mxdtyp,mxdnr,
     2                 vpp,fpp,nr,dr,r0)
      w = w + wpp
      ut = ut + upp
      do 340 jj = 1,3
      do 340 ii = 1,3
      frr(ii,jj) = frr(ii,jj) + fr(ii,jj)
 340  continue
      do 350 k=1,3
      f(k,i) = f(k,i) - fo(k)
      f(k,j) = f(k,j) + fo(k)
      fint(k,i) = fint(k,i) - fo(k)
      fint(k,j) = fint(k,j) + fo(k)
      fs(k,i) = fs(k,i) - fos(k)
      fs(k,j) = fs(k,j) + fos(k)
      fsint(k,i) = fsint(k,i) - fos(k)
      fsint(k,j) = fsint(k,j) + fos(k)
 350  continue
      endif
 400  continue
      else
c
c     particles i and j in different cells
c
      do 700 nat = 1,natot
      do 600 i = 1,3
      rad(i,nat) = zero
      card(i,nat) = carn(i,nat) + dfloat(i1) * avec(i,1) +
     +           dfloat(i2) * avec(i,2) +
     +           dfloat(i3) * avec(i,3)
 600  continue
      do 650 k=1,3
      do 650 l=1,3
      rad(k,nat) = rad(k,nat) + sigma(l,k) / vcell * card(l,nat)
 650  continue
 700  continue
      do 900 i = 1,natot
      do 900 j = 1,natot
      do 750 k=1,3
      rd(k) = card(k,j) - carn(k,i)
      sd(k) = rad(k,j) - ran(k,i)
 750  continue
      rdsq = rd(1) ** 2 + rd(2) ** 2 + rd(3) ** 2
      if (rdsq .le. rcutsq) then
      nt1 = indt(i)
      nt2 = indt(j)
      call vpair(rdsq,rd,sd,upp,fo,fos,fr,wpp,
     1                 calc,nt1,nt2,mxdtyp,mxdnr,
     2                 vpp,fpp,nr,dr,r0)
      w = w + wpp / dois
      ut = ut + upp / dois
      do 760 jj = 1,3
      do 760 ii = 1,3
      frr(ii,jj) = frr(ii,jj) + fr(ii,jj) / dois
 760  continue
      do 800 k=1,3
      f(k,i) = f(k,i) - fo(k)
      fs(k,i) = fs(k,i) - fos(k)
 800  continue
      endif
 900  continue
      endif
1000  continue
c
c     calculate forces w.r.t center of mass (just to check)
c
      fx = zero
      fy = zero
      fz = zero
      do 1300 nat=1,natot
      fx = fx + f(1,nat)
      fy = fy + f(2,nat)
      fz = fz + f(3,nat)
1300  continue
c      write(6,*) 'fx,fy,fz',fx,fy,fz
c
c     calculate acceleration
c
      do 1400 nat = 1,natot
      do 1400 k=1,3
      a(k,inda(nat),indt(nat)) = f(k,nat) / atmass(indt(nat))
      rat2d(k,inda(nat),indt(nat)) = fs(k,nat) / atmass(indt(nat))
1400  continue
      return
      end
c*
c*
      subroutine vpair(rdsq,rd,sd,upp,fo,fos,fr,wpp,
     1                 calc,nt1,nt2,mxdtyp,mxdnr,
     2                 vpp,fpp,nr,dr,r0)
c
c     computes the potential energy and forces between 2 interacting
c     particles separated by distance r. interpolates
c     quadratically 2 functions vpp and dvpp/dr given on a regular spaced grid.
c     written by rmw on 5/3/90
c
c     input:
c     rdsq = square of distance between the particles
c     rdm = distance between the particles
c     rd(3) = relative position of particles in cartesian coordinates
c     sd(3) = relative position of particles in lattice coordinates
c     nt1 and nt2 = type of interacting particles
c     vpp(nt1,nt2,nr)=potential energy between atoms of type nt1 and nt2
c     fpp(nt1,nt2,nr)=force between atoms of type nt1 and nt2
c     nr(nt1,nt2)=# of points in the radial mesh of f and v
c     dr(nt1,nt2)= radial spacing for f and v
c     r0(nt1,nt2)=1st point in the radial mesh for f and v
c
c     output:
c     upp = pair potential energy
c     wpp = pair contribution to the virial pressure
c     fr(3,3) = "     "       to the cell forces
c     fo(3) = force in cartesian coordinated
c     fos(3) = "     "   lattice  "

      implicit real*8 (a-h,o-z)
c
      data zero, um, dois, tres, quatro /0.d0, 1.d0, 2.d0, 3.d0, 4.d0/
c
      dimension rd(3), fo(3),fr(3,3)
      dimension sd(3), fos(3)
      dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
      dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
c
       character*2 calc
c
      a0 = .529177d0
      ry = 13.6058d0
c
c     potential energy
c
c
      rdm = dsqrt(rdsq)
      if (rdm .le. r0(nt1,nt2))
     ~    call fatal(41,rdm,idum)
      if (rdm .le. (r0(nt1,nt2) + dr(nt1,nt2)
     *              * (nr(nt1,nt2) - 1))) then
      rn = (rdm - r0(nt1,nt2))/dr(nt1,nt2)
      n = rn + 0.5d0 + 1
      if (n .eq. 1) n = n + 1
      if (n .eq. nr(nt1,nt2)) n = n-1
      x = rdm - (r0(nt1,nt2) + (n - 1) * dr(nt1,nt2))
c
c     calculate potential energy
c
      c = vpp(nt1,nt2,n)
      b = (vpp(nt1,nt2,n+1) - vpp(nt1,nt2,n-1))/(dois * dr(nt1,nt2))
      a = (vpp(nt1,nt2,n+1) + vpp(nt1,nt2,n-1) - dois * c)/
     /    (dois * dr(nt1,nt2) * dr(nt1,nt2))
      upp = a * x * x + b * x + c
c
c     force
c
      c = fpp(nt1,nt2,n)
      b = (fpp(nt1,nt2,n+1) - fpp(nt1,nt2,n-1))/(dois * dr(nt1,nt2))
      a = (fpp(nt1,nt2,n+1) + fpp(nt1,nt2,n-1) - dois * c)/
     /    (dois * dr(nt1,nt2) * dr(nt1,nt2))
      factor = a * x * x + b * x + c
      wpp = zero
      do 100 k=1,3
      fo(k) = - factor * rd(k) / rdm
      wpp = wpp + fo(k) * rd(k)
      fos(k) = - factor * sd(k) / rdm
 100  continue
      else
      upp = zero
      wpp = zero
      do 110 k = 1,3
      fo(k) = zero
      fos(k) = zero
 110  continue
      endif
c
c     for cell shape md
c
       if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
      do 300 j = 1,3
      do 200 i = 1,3
      fr(i,j) = fo(i) * rd(j)
 200  continue
 300  continue
      endif
      return
      end
c*
c*
      subroutine move(alatt,vcell,
     1          ntype,natom,natot,atmass,cmass,press,
     2          mxdtyp,mxdatm,mxdnat,mxdnr,
     3          rcut,ncell,nst,nstep,calc,dt,
     4          avec,avecd,avec2d,avec2di,avmod,theta,scp,
     5          avec0,sig0,v0,g,gd,gm1,sigma,pim,
     6          rat,ratd,rat2d,rati,rat2di,car,v,a,
     7          card,carn,ran,rad,vpp,fpp,r0,dr,nr,
     8          indt,inda,f,fint,fs,fsint,frr,vint,aint,
     9          uta,eka,eta,utl,ekl,etl,ut,p,ekint,etot)
c
c     subroutine written by rmw on 2/5/90
c     integrates the eq.'s of motion according to the Beeman algorythm
c     assumes the initial quantities q, dq/dt, and d2q/dt2 are given
c     in lattice coordinates
c
c     input:
c     mxdtyp = array dimension for type of atoms
c     mxdatm = array dimension for atoms of a given type
c     mxdnat = array dimension for total number of atoms
c     ntype = number of types of atoms
c     natom(nt) = number of atoms of type nt
c     atmass(nt) = atomic masses for atoms of type nt (in proton masses)
c     rat(j,na,nt) = atomic positions in lattice coordinates
c     ratd(j,na,nt) = atomic velocities     "        "
c     rat2d(i,na,nt) =   "   acceleration    "        "
c     avec(3,3) = lattice vectors
c     avecd(3,3) = 1st lattice vectors derivatives
c     avec2d(3,3) = 2nd lattice vectors derivatives
c     avec0(3,3) = initial lattice vectors
c     sigma(3,3) = reciprocal lattice vectors *  vcell / 2 pi
c     sig0(3,3) = = initial reciprocal lattice vectors *  vcell / 2 pi
c
c     output:
c     rat(j,na,nt) = atomic positions in lattice coordinates
c     ratd(j,na,nt) = atomic velocities     "        "
c     rat2d(i,na,nt) =   "   acceleration    "        "
c     avec(3,3) = lattice vectors
c     avecd(3,3) = 1st lattice vectors derivatives
c     avec2d(3,3) = 2nd lattice vectors derivatives
c     p = internal (virial) pressure
c     ut = new total potential energy
c     ekin = new total kinetic energy
c     etot = total energy
c     we also obtain the same quantities for atomic and lattice components
c     uta,eka,eta,utl,ekl,etl
c     theta(3,3) = angle between lattice vectors
c     avmod(3) = lattice vectors moduli
c     scp(3) = scalar product betwen lattice vectors / respective moduli
c
      implicit real*8 (a-h,o-z)
c
      dimension avec(3,3),avecd(3,3),avec2d(3,3)
      dimension avec0(3,3),sig0(3,3)
      dimension aveci(3,3),avec2di(3,3)
      dimension g(3,3),gm1(3,3),gd(3,3),gmgd(3,3)
      dimension sigma(3,3),pim(3,3),sigav(3,3)
      dimension natom(mxdtyp),atmass(mxdtyp)
      dimension ncell(3)
      dimension car(3,mxdatm,mxdtyp),v(3,mxdatm,mxdtyp)
      dimension a(3,mxdatm,mxdtyp)
      dimension vint(3,mxdatm,mxdtyp),aint(3,mxdatm,mxdtyp)
      dimension rat(3,mxdatm,mxdtyp)
      dimension ratd(3,mxdatm,mxdtyp),rat2d(3,mxdatm,mxdtyp)
      dimension rati(3,mxdatm,mxdtyp),rat2di(3,mxdatm,mxdtyp)
      dimension carn(3,mxdnat), card(3,mxdnat)
      dimension ran(3,mxdnat),rad(3,mxdnat)
      dimension inda(mxdnat),indt(mxdnat)
      dimension f(3,mxdnat),fint(3,mxdnat),frr(3,3)
      dimension fs(3,mxdnat),fsint(3,mxdnat)
      dimension theta(3,3),avmod(3),scp(3)
      dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
      dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
      dimension d2(3,3)
c
       parameter ( pi = 3.141592653589793d0, twopi = 2.0d0*pi)
       parameter ( quatro = 4.0d0, seis = 6.0d0)
c
       character*2 calc
c
      data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
c
c     zero energy components
c
      ut = zero
      ekint = zero
      etot = zero
      uta = zero
      eka = zero
      eta = zero
      utl = zero
      ekl = zero
      etl = zero
      p = zero
c
c     update atomic positions and calculate intermediate velocities
c     and accelerations
c
      do 100 nt = 1,ntype
      do 100 na = 1,natom(nt)
      do 100 k = 1, 3
      rati(k,na,nt) = rat(k,na,nt)
      rat(k,na,nt) = rat(k,na,nt) + dt * ratd(k,na,nt) +
     +               dt * dt * (quatro * rat2d(k,na,nt) - 
     -               rat2di(k,na,nt)) / seis
      rati(k,na,nt) = rat(k,na,nt) - rati(k,na,nt)
      rat2di(k,na,nt) = rat2d(k,na,nt)
      ratd(k,na,nt) = ratd(k,na,nt) + dt * rat2d(k,na,nt)
c
c     move particles back to the box
c
      if (rat(k,na,nt) .gt. um) then
      rat(k,na,nt) = dmod(rat(k,na,nt),um)
      endif
      if (rat(k,na,nt) .lt. zero) then
      x = dabs(rat(k,na,nt)) 
      y = dmod(x,um)
      rat(k,na,nt) = um - y
      endif
 100  continue
c
c     update lattice vectors if cell dynamics
c
       if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
      do 140 j = 1,3
      do 140 i = 1,3
      aveci(i,j) = avec(i,j)
      avec(i,j) = avec(i,j) + dt * avecd(i,j) +
     +            dt*dt * (quatro * avec2d(i,j) -
     -            avec2di(i,j)) / seis
      avec2di(i,j) = avec2d(i,j)
      avecd(i,j) = avecd(i,j) + dt * avec2d(i,j)
 140  continue
c
      itg = 1
      call updg(itg,avec,avecd,
     1          g,gd,gm1,gmgd,sigma,vcell)
       endif
c
c     back to cartesian coordinates
c
      do 250 nt=1,ntype
      do 250 na=1,natom(nt)
      do 200 k=1,3
      car(k,na,nt) = zero
      v(k,na,nt) = zero
      do 200 l=1,3
      car(k,na,nt) = car(k,na,nt) +
     +               avec(k,l) * rat(l,na,nt)
      v(k,na,nt) = v(k,na,nt) +
     +             avec(k,l) * ratd(l,na,nt)
 200  continue
c
c    write atomic positions on tape for later plotting
c
c     ra1 = rat(1,na,nt) - rat(1,1,1)
c     ra2 = rat(2,na,nt) - rat(2,1,1)
c     ra3 = rat(3,na,nt) - rat(3,1,1)
c     write(20,302) ra1,ra2,ra3
c
 250  continue
c
c    calculate new accelerations from forces
c
      ilj = 1
      if (ilj .eq. 1) then
      call forclj(avec,sigma,vcell,calc,
     1          ntype,natom,natot,atmass,
     2          mxdtyp,mxdatm,mxdnat,
     3          rcut,ncell,
     4          car,v,a,rat,ratd,rat2d,
     5          card,carn,rad,ran,
     6          indt,inda,f,fs,fint,fsint,frr,ut,w)
      else
      call forc(avec,sigma,vcell,calc,
     1          ntype,natom,natot,atmass,
     2          mxdtyp,mxdatm,mxdnat,mxdnr,
     3          rcut,ncell,
     4          car,v,a,rat,ratd,rat2d,
     5          card,carn,rad,ran,
     6          indt,inda,f,fs,fint,fsint,frr,ut,w,
     7          vpp,fpp,r0,dr,nr)
      endif
c
c     update cell related quantities
c
      if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
      iloop = 1
c
c     zero pim
c
 700  do 750 j =1,3
      do 750 i =1,3
      pim(i,j) = zero
 750  continue
c
c     obtain new forces
c
      do 830 nt = 1,ntype
      do 830 na = 1,natom(nt)
      do 800 k = 1,3
      do 800 m = 1,3
      rat2d(k,na,nt) = rat2d(k,na,nt) -
     -                 gmgd(k,m) * ratd(m,na,nt)
 800  continue
c
c     initialize forces on basis vectors
c
      do 820 j = 1,3
      do 820 i = 1,3
      pim(i,j) = pim(i,j) + atmass(nt) * v(i,na,nt) *
     *           v(j,na,nt)
 820  continue
 830  continue
      do 840 j = 1,3
      do 840 i = 1,3
      pim(i,j) = (pim(i,j) + frr(i,j)) / vcell
c
c     diagonal
c
      if (i .eq. j) then
      pim(i,j) = pim(i,j) - press
      endif
      avec2d(i,j) = zero
 840  continue
c
      do 860 j = 1,3
      do 860 i = 1,3
      do 850 k = 1,3
      avec2d(i,j) = avec2d(i,j) + pim(i,k) *
     *              sigma(k,j)
 850  continue
      avec2d(i,j) = avec2d(i,j) / cmass
 860  continue
c
c     if new cell dynamics...
c
      if ((calc .eq. 'nd') .or. (calc .eq. 'nm')) then
      call sigp(avec,avecd,avec2d,sigma,vcell)
      endif
c
      if ((calc .eq. 'sd') .or. (calc .eq. 'sm')) then
      call sigs(avec0,avec2d,sig0,v0)
      endif
c
c     strain symmetrization
c
      do 862 i=1,3
      do 862 j=1,3
        d2(i,j) = zero
        do 861 k=1,3
          d2(i,j) = d2(i,j) + avec2d(i,k) * sig0(j,k)
 861    continue
        d2(i,j) = d2(i,j) / v0
 862  continue
c
c
      do 863 i=1,2
      do 863 j=i+1,3
      d2(i,j) = (d2(i,j) + d2(j,i)) / dois
      d2(j,i) = d2(i,j)
 863  continue
c
c     write(6,*) d2(2,1),d2(3,1),d2(3,2)
c     write(6,*) d2(1,2),d2(1,3),d2(2,3)
c
      do 865 i=1,3
      do 865 j=1,3
        avec2d(i,j) = zero
        do 864 k=1,3
          avec2d(i,j) = avec2d(i,j) + d2(i,k) * avec0(k,j)
 864    continue
 865    continue
c     write(6,*) avec2d(2,1),avec2d(3,1),avec2d(3,2)
c
c     correct atomic velocities...
c
      do 940 nt=1,ntype
      do 940 na=1,natom(nt)
      do 900 k=1,3
      ratd(k,na,nt) = rati(k,na,nt) / dt + 
     +                dt * (dois * rat2d(k,na,nt) +
     +                rat2di(k,na,nt)) / seis
 900  continue
      do 930 k = 1,3
      do 920 l = 1,3
      v(k,na,nt) = v(k,na,nt) +
     +             avec(k,l) * ratd(l,na,nt)
 920  continue
 930  continue
 940  continue
c
c     ... and lattice velocities
c
      do 950 j = 1,3
      do 950 i = 1,3
      avecd(i,j) = (avec(i,j) - aveci(i,j)) / dt +
     +            dt* (dois * avec2d(i,j) + avec2di(i,j))/ seis
 950  continue
      iloop = iloop - 1
      if (iloop .gt. 0) then
      itg = 0
      call updg(itg,avec,avecd,
     1          g,gd,gm1,gmgd,sigma,vcell)
      goto 700
      endif
      else
      do 990 nt=1,ntype
      do 990 na=1,natom(nt)
      do 960 k=1,3
      ratd(k,na,nt) = rati(k,na,nt) / dt + 
     +                dt * (dois * rat2d(k,na,nt) +
     +                rat2di(k,na,nt)) / seis
 960  continue
      do 980 k = 1,3
      do 980 l = 1,3
      v(k,na,nt) = v(k,na,nt) +
     +             avec(k,l) * ratd(l,na,nt)
 980  continue
 990  continue
      endif
c
c     calculate basis vectors moduli and angles
c
      if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
      do 1100 k = 1,3
      avmod(k) = zero
      do 1050 l = 1,3
      theta(l,k) = zero
      avmod(k) = avmod(k) + avec(l,k) * avec(l,k)
      do 1000 m = 1,3
      theta(l,k) = theta(l,k) + avec(m,l) * avec(m,k)
 1000  continue
 1050 continue
      avmod(k) = dsqrt(avmod(k))
 1100 continue
      do 1150 k = 1,3
      do 1150 l = 1,3
      x = theta(l,k) / avmod(l) / avmod(k)
      if ((l .eq. 1) .and. (k .eq. 2)) scp(1) = x
      if ((l .eq. 2) .and. (k .eq. 3)) scp(2) = x
      if ((l .eq. 3) .and. (k .eq. 1)) scp(3) = x
      if (x .ge. 0.d0) then
      x = dmin1(1.d0,x)
cray      x = min(1.d0,x)
      else
      x = dmax1(-1.d0,x)
cray      x = max(-1.d0,x)
      endif
      theta(l,k) = dacos(x) * 180.d0 / pi
 1150 continue
      endif
c
c     compute kinetic energy
c
      ekint = zero
      do 1250 nt = 1,ntype
      do 1250 na = 1,natom(nt)
      do 1250 i=1,3
      ek = zero
      do 1200 j=1,3
      ek = ek + ratd(i,na,nt) * g(i,j) * ratd(j,na,nt)
 1200 continue
      ekint = ekint + ek * atmass(nt) / dois
 1250 continue
      etot = ekint + ut
      uta = ut
      eka = ekint
      eta = etot
c
c     lattice contribution
c
      ekl = zero
      if ((calc .ne. 'md') .and. (calc .ne. 'mm'))then
c
      if ((calc .eq. 'nd') .or. (calc .eq. 'nm'))then
      do 1300 j = 1,3
      do 1300 i = 1,3
      sigav(i,j) = zero
      do 1300 l = 1,3
      sigav(i,j) = sigav(i,j) + sigma(l,i) * avecd(l,j)
 1300 continue
      do 1400 k = 1,3
      tr = zero
      do 1350 m = 1,3
      tr = tr + sigav(m,k) * sigav(m,k)
 1350 continue
      ekl = ekl + tr
 1400 continue
      endif
c
      if ((calc .eq. 'sd') .or. (calc .eq. 'sm'))then
      do 1550 j = 1,3
      do 1550 i = 1,3
      sigav(i,j) = zero
      do 1550 l = 1,3
      sigav(i,j) = sigav(i,j) + sig0(l,i) * avecd(l,j)
 1550 continue
      do 1650 k = 1,3
      tr = zero
      do 1600 m = 1,3
      tr = tr + sigav(m,k) * sigav(m,k)
 1600 continue
      ekl = ekl + tr
 1650 continue
      endif
c
      if ((calc .eq. 'cd') .or. (calc .eq. 'cm'))then
      do 1750 k = 1,3
      tr = zero
      do 1700 m = 1,3
      tr = tr + avecd(m,k) * avecd(m,k)
1700  continue
      ekl = ekl + tr
1750  continue
      endif
c
      utl = + press * vcell
      ekl = ekl * cmass / dois
      etl = utl + ekl
c
c     total energy
c
      ekint = eka + ekl
      ut = uta + utl
      etot = ekint + ut
      endif
c
c     calculate "internal (virial) pressure"
c
      p = ( dois * eka + w) / tres /  vcell
c
c
      return
 302  format(1x,3f10.4)
      end
c*
c*
      subroutine updg(itg,avec,avecd,
     1                g,gd,gm1,gmgd,sigma,vcell)
c
c
c     calculates cell related quantities in m.d.
c     calculation involving the cell shape.
c
c     written by r. m. wentzcovitch on 4/11/90
c
c      input:
c      avec(3,3) = lattice vectors
c      avecd(3,3) = derivative of lattice vectors
c      sigma(3,3) = reciprocal lattice vectors / twopi
c      vcell = cell volume
c      itg = tells if sigma and vcell have to be calculated
c                    1 , yes
c             itg =
c                    0 , no
c
c      output:      t
c      g(3,3) = avec * avec
c                    t                      t
c      gd(3,3) = avecd * avec + avecd * avec
c                  _1
c      gm1(3,3) = g
c                   _1
c      gmgd(3,3) = g * gd
c      sigma(3,3) = reciprocal lattice vectors / twopi
c      vcell = cell volume
c
c
       implicit double precision (a-h,o-z)
       parameter ( eps = 1.0d-14)
       parameter (zero = 0.0d0,dois = 2.0d0,um = 1.0d0,tres = 3.0d0)
c
       dimension avec(3,3),avecd(3,3),sigma(3,3)
       dimension g(3,3),gd(3,3),gmgd(3,3),gm1(3,3)
c
c
       if (itg .eq. 1) then
c
c      compute the lattice wave-vectors/twopi and the cell volume
c
       sigma(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
       sigma(2,1) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
       sigma(3,1) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
       sigma(1,2) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
       sigma(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
       sigma(3,2) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
       sigma(1,3) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
       sigma(2,3) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
       sigma(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)
c
c      cell volume
c
       vcell = sigma(1,1)*avec(1,1) + sigma(2,1)*avec(2,1) +
     +         sigma(3,1)*avec(3,1)
       if(vcell .lt. eps) call fatal(50,vcell,idum)
       vcell = dabs(vcell)
       endif
c
c      calculate g, gd, and gm1 matrices
c
       do 20 j = 1,3
       do 20 i = 1,3
       g(i,j) = zero
       gm1(i,j) = zero
       gd(i,j) = zero
 20    continue
       do 26 j = 1,3
       do 26 i = 1,3
       do 25 m = 1,3
       g(i,j) = g(i,j) + avec(m,i) * avec(m,j)
       gm1(i,j) = gm1(i,j) + sigma(m,i) * sigma(m,j)
       gd(i,j) = gd(i,j) + avec(m,i) * avecd(m,j) +
     +           avecd(m,i) * avec(m,j)
 25    continue
       gm1(i,j) = gm1(i,j) / vcell / vcell
 26    continue
c                _1  .
c     calculate g * g    ( = gmgd)
c
      do 40 j = 1,3
      do 40 i = 1,3
      gmgd(i,j) = zero
      do 40 m = 1,3
      gmgd(i,j) = gmgd(i,j) + gm1(i,m) * gd(m,j)
 40   continue
      return
      end
c*
c*
       subroutine rdpp(mxdtyp,mxdnr,ntype,vpp,fpp,r0,dr,nr)
c
c      subroutine reads in ntype * (ntype + 1) pair potentials
c      and the corresponding forces.
c      assumes potentials given in Ry, distances in a^o and
c      forces in Ry/a^o
c
      implicit double precision (a-h,o-z)
      parameter ( dois = 2.0d0)
      parameter ( zero = 0.0d0, um = 1.0d0, tres = 3.0d0)
c
      dimension vpp(mxdtyp,mxdtyp,mxdnr),fpp(mxdtyp,mxdtyp,mxdnr)
      dimension nr(mxdtyp,mxdtyp),dr(mxdtyp,mxdtyp),r0(mxdtyp,mxdtyp)
c
      a0 = .529177d0
      ry = 13.6058d0
c
c     read potentials in the order (1,1), (1,2),...(ntype-1,ntype)
c
      npot = 0
      do 200 nt = 1,ntype
        do 100 ntp = nt,ntype
          npot = npot + 1
          read(30+npot,1001) nr(nt,ntp)
          if (nr(nt,ntp) .gt. mxdnr) call fatal(31,xdum,nr(nt,ntp))
          write(6,2000) nt,ntp
          write(6,2001) nr(nt,ntp)
          nr(ntp,nt) = nr(nt,ntp)
          read(30+npot,1002) r0(nt,ntp),dr(nt,ntp)
c         r0(nt,ntp) = r0(nt,ntp)/ a0 
c         dr(nt,ntp) = dr(nt,ntp)/ a0 
          write(6,2002) r0(nt,ntp),dr(nt,ntp)
          r0(ntp,nt) = r0(nt,ntp)
          dr(ntp,nt) = dr(nt,ntp)
          do 50 mr = 1,nr(nt,ntp)
            read(30+npot,1003) vpp(nt,ntp,mr),fpp(nt,ntp,mr),x
c           vpp(nt,ntp,mr) = vpp(nt,ntp,mr) / ry
c           fpp(nt,ntp,mr) = fpp(nt,ntp,mr) / ry * a0
            vpp(ntp,nt,mr) = vpp(nt,ntp,mr)
            fpp(ntp,nt,mr) = fpp(nt,ntp,mr)
  50      continue
 100    continue
 200  continue
c
c
      return
1001  format(i5)
1002  format(2d15.5)
1003  format(4d20.8)
2000  format(/,'potential between atoms type',i3,' and',i3)
2001  format(/,1x,i3,' mesh points')
2002  format(/,1x,'r0 = ',d13.5,' dr = ',d13.5)
2003  format(1x,3d20.8)
      end
c*
c*
       subroutine sigp(avec,avecd,avec2d,sigma,vcell)
c
c      calculates sigmap matrices and avec2d for
c      new dynamics(rmw 5/30/90)
c
c      input:
c            avec = lattice vectors
c            avecd = time derivative of lattice vectors
c            avec2d = 2nd time derivative of lattice vectors
c            sigma = volume * rec. latt. vectors / 2 pi
c            vcell = cell volume
c
c      output:
c            avec2d = new 2nd time derivative of lattice vectors
c
       implicit double precision (a-h,o-z)
c
       dimension avec(3,3),avecd(3,3),avec2d(3,3)
       dimension sigmap(3,3,3,3),sigmad(3,3),sigma(3,3)
       dimension e(3,3),fp(3,3,3,3),fd(3,3),fm1(3,3)
       dimension fm(3,3),sm(3,3),avint(3,3)
c
      data zero , dois / 0.d0,  2.d0/
c
      sigmap(1,1,1,1) = zero
      sigmap(1,1,1,2) = zero
      sigmap(1,1,1,3) = zero
      sigmap(1,1,2,1) = zero
      sigmap(1,1,2,2) = avec(3,3)
      sigmap(1,1,2,3) = - avec(3,2)
      sigmap(1,1,3,1) = zero
      sigmap(1,1,3,2) = - avec(2,3)
      sigmap(1,1,3,3) = avec(2,2)
c
      sigmap(1,2,1,1) = zero
      sigmap(1,2,1,2) = zero
      sigmap(1,2,1,3) = zero
      sigmap(1,2,2,1) = - avec(3,3)
      sigmap(1,2,2,2) = zero
      sigmap(1,2,2,3) = avec(3,1)
      sigmap(1,2,3,1) = avec(2,3)
      sigmap(1,2,3,2) = zero
      sigmap(1,2,3,3) = - avec(2,1)
c
      sigmap(1,3,1,1) = zero
      sigmap(1,3,1,2) = zero
      sigmap(1,3,1,3) = zero
      sigmap(1,3,2,1) = avec(3,2)
      sigmap(1,3,2,2) = - avec(3,1)
      sigmap(1,3,2,3) = zero
      sigmap(1,3,3,1) = - avec(2,2)
      sigmap(1,3,3,2) = avec(2,1)
      sigmap(1,3,3,3) = zero
c
      sigmap(2,1,1,1) = zero
      sigmap(2,1,1,2) = - avec(3,3)
      sigmap(2,1,1,3) = avec(3,2)
      sigmap(2,1,2,1) = zero
      sigmap(2,1,2,2) = zero
      sigmap(2,1,2,3) = zero
      sigmap(2,1,3,1) = zero
      sigmap(2,1,3,2) = avec(1,3)
      sigmap(2,1,3,3) = - avec(1,2)
c
      sigmap(2,2,1,1) = avec(3,3)
      sigmap(2,2,1,2) = zero
      sigmap(2,2,1,3) = - avec(3,1)
      sigmap(2,2,2,1) = zero
      sigmap(2,2,2,2) = zero
      sigmap(2,2,2,3) = zero
      sigmap(2,2,3,1) = - avec(1,3)
      sigmap(2,2,3,2) = zero
      sigmap(2,2,3,3) = avec(1,1)
c
      sigmap(2,3,1,1) = - avec(3,2)
      sigmap(2,3,1,2) = avec(3,1)
      sigmap(2,3,1,3) = zero
      sigmap(2,3,2,1) = zero
      sigmap(2,3,2,2) = zero
      sigmap(2,3,2,3) = zero
      sigmap(2,3,3,1) = avec(1,2)
      sigmap(2,3,3,2) = - avec(1,1)
      sigmap(2,3,3,3) = zero
c
      sigmap(3,1,1,1) = zero
      sigmap(3,1,1,2) = avec(2,3)
      sigmap(3,1,1,3) = - avec(2,2)
      sigmap(3,1,2,1) = zero
      sigmap(3,1,2,2) = - avec(1,3)
      sigmap(3,1,2,3) = avec(1,2)
      sigmap(3,1,3,1) = zero
      sigmap(3,1,3,2) = zero
      sigmap(3,1,3,3) = zero
c
      sigmap(3,2,1,1) = - avec(2,3)
      sigmap(3,2,1,2) = zero
      sigmap(3,2,1,3) = avec(2,1)
      sigmap(3,2,2,1) = avec(1,3)
      sigmap(3,2,2,2) = zero
      sigmap(3,2,2,3) = - avec(1,1)
      sigmap(3,2,3,1) = zero
      sigmap(3,2,3,2) = zero
      sigmap(3,2,3,3) = zero
c
      sigmap(3,3,1,1) = avec(2,2)
      sigmap(3,3,1,2) = - avec(2,1)
      sigmap(3,3,1,3) = zero
      sigmap(3,3,2,1) = - avec(1,2)
      sigmap(3,3,2,2) = avec(1,1)
      sigmap(3,3,2,3) = zero
      sigmap(3,3,3,1) = zero
      sigmap(3,3,3,2) = zero
      sigmap(3,3,3,3) = zero
c                _1  t           2
c     calculate f = h * h / vcell
c
      do 50 j = 1,3
      do 50 i = 1,3
      fm1(i,j) = zero
      do 45 l = 1,3
      fm1(i,j) = fm1(i,j) + avec(l,i) * avec(l,j)
 45   continue
      fm1(i,j) = fm1(i,j) / vcell / vcell
 50   continue
c                   .t  .
c     calculate e = h * h
c
      do 100 j = 1,3
      do 100 i = 1,3
      e(i,j) = zero
      do 100 m =1,3
      e(i,j) = e(i,j) + avecd(m,i) * avecd(m,j)
 100  continue
c                          ij t             t       ij
c     calculate f' = sigma'  * sigma + sigma * sigma'
c
      do 200 n = 1,3
      do 200 m = 1,3
      do 200 j = 1,3
      do 200 i = 1,3
      fp(i,j,m,n) = zero
      do 200 l = 1,3
      fp(i,j,m,n) = fp(i,j,m,n) + sigmap(i,j,l,m) *
     *              sigma(l,n) + sigma(l,m) * sigmap(i,j,l,n)
 200  continue
c
c     calculate sigmad
c
      do 300 n = 1,3
      do 300 m = 1,3
      sigmad(m,n) = zero
      do 300 j = 1,3
      do 300 i = 1,3
      sigmad(m,n) = sigmad(m,n) + sigmap(i,j,m,n) *
     *              avecd(i,j)
 300  continue
c               .
c     calculate f
c
      do 400 j = 1,3
      do 400 i = 1,3
      fd(i,j) = zero
      do 400 l = 1,3
      fd(i,j) = fd(i,j) + sigmad(l,i) * sigma(l,j) +
     +          sigma(l,i) * sigmad(l,j)
 400  continue
c
c     calculate fm
c
      do 500 j = 1,3
      do 500 i = 1,3
      fm(i,j) = zero
      do 450 l = 1,3
      do 450 k = 1,3
      fm(i,j) = fm(i,j) + e(l,k) * fp(i,j,k,l)
 450  continue
      fm(i,j) = fm(i,j) / dois
 500  continue
c
c     calculate sm
c
      do 600 j = 1,3
      do 600 i = 1,3
      sm(i,j) = zero
      do 600 l = 1,3
      sm(i,j) = sm(i,j) +  avecd(i,l) * fd(l,j)
 600  continue
c
c     calculate new avec2d
c
      do 650 j = 1,3
      do 650 i = 1,3
      avint(i,j) = avec2d(i,j) + fm(i,j) - sm(i,j)
 650  continue
c
c
      do 800 j = 1,3
      do 800 i = 1,3
      avec2d(i,j) = zero
      do 800 m = 1,3
      avec2d(i,j) = avec2d(i,j) + avint(i,m) * fm1(m,j)
 800  continue
c
      return
      end
c*
c*
c*
       subroutine sigs(avec0,avec2d,sig0,v0)
c
c      calculates avec2d for strain dynamics (rmw 9/26/90)
c
c      input:
c            avec0 = initial lattice vectors
c            avec2d = 2nd time derivative of lattice vectors
c            sig0 = volume * rec. latt. vectors / 2 pi
c            v0 = cell volume
c
c      output:
c            avec2d = new 2nd time derivative of lattice vectors
c
       implicit double precision (a-h,o-z)
c
       dimension avec0(3,3),avec2d(3,3),sig0(3,3)
       dimension fm1(3,3),avint(3,3)
c
      data zero , dois / 0.d0,  2.d0/
c
c                _1   t          2
c     calculate f = ho * ho / vo
c
      do 50 j = 1,3
      do 50 i = 1,3
      fm1(i,j) = zero
      do 45 l = 1,3
      fm1(i,j) = fm1(i,j) + avec0(l,i) * avec0(l,j)
 45   continue
      fm1(i,j) = fm1(i,j) / v0 / v0
 50   continue
c
c     calculate new avec2d
c
      do 650 j = 1,3
      do 650 i = 1,3
      avint(i,j) = avec2d(i,j)
 650  continue
c
c
      do 800 j = 1,3
      do 800 i = 1,3
      avec2d(i,j) = zero
      do 800 m = 1,3
      avec2d(i,j) = avec2d(i,j) + avint(i,m) * fm1(m,j)
 800  continue
c
      return
      end
c*
c*
       subroutine fatal(i,xarg,iarg)
c
c      handles the fatal errors
c      and stops the execution of the program
c
       implicit double precision (a-h,o-z)
c
       write(6,1000)
       i10 = i/10
       ir = i - i10*10
       if(i10 .eq. 3) then
c      rdpp
       if (ir .eq. 1) then
          write(6,1031) iarg
       endif
       endif
c
       if(i10 .eq. 4) then
c      vpair
       if (ir .eq. 1) then
          write(6,1041) xarg
       endif
       endif
       if(i10 .eq. 5) then
c        crstl or updg
         if(ir . eq. 0) then
           write(6,1050) xarg
         else if(ir .eq. 1) then
           write(6,1051) iarg
         else if(ir .eq. 2) then
           write(6,1052) iarg
         else if(ir .eq. 3) then
           write(6,1053) iarg
         endif
       endif
       call abort()
 1000  format('  ***fatal error***')
 1031  format('  number of mesh points too large = ',i5)
 1041  format('  interatomic distance too small = ',d12.5)
 1050  format('  cell volume = ',d12.5)
 1051  format('  types of atoms = ',i5)
 1052  format('  number of atoms = ',i5)
 1053  format('  number of steps = ',i5)
       end
c*
c*
      double precision function ran3(idum)
      implicit double precision (a-h,o-z)
      save
c         implicit real*4(m)
c         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d-9)
      dimension ma(55)
c     common /ranz/ ma,inext,inextp
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end
c
c
      double precision function ran3b(idum)
      implicit double precision (a-h,o-z)
c         implicit real*4(m)
c         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d-9)
      dimension ma(55)
      common /ranz/ ma,inext,inextp
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3b=mj*fac
      return
      end
c*
c*
