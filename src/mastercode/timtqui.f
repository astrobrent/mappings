cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******FINDS THE EQUILIBRIUM TEMPERATURE AND THE
c	CORRESPONDING IONISING STATE OF THE GAS AFTER TIME TSTEP
c	TEI : INITIAL GUESS FOR TEMPERATURE
c	TEF,EDENS: FINAL TEMPERATURE AND ELECTRONIC DENSITY
c	CALL SUBR. COOL,TIMION
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine timtqui(tei, tef, edens, hdens, tstep)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision a,b,c1,c2,cc,cof,csub,delmin
      double precision dimax,dinc,dl1,emul,ex
      double precision t2,tl1,tl2,tm1,tmin
      double precision tsub,tti,u,xhyf
      double precision popp0(mxion, mxelem), tei, tef, edens, hdens
     &                 , tstep
c
      integer*4 iel,ion,km
      integer*4 luop,n,nf
c
      character ntest*4
c
      double precision arctanh
c
      arctanh(u) = dlog((1.d0+u)/(1.d0-u))/2.d0
      ntest = 'NO'
      if (tef.lt.0.0d0) ntest = 'Y'
      luop = 23
      delmin = 1.d-3
      nf = 12
      dimax = 1.4d0
      ex = 1.0d0
      tmin = 500.d0
      emul = 4.0d0
      cof = 0.8d0
c
c    ***USING 2 TEMP. (TM1,T2)  AND 2 FRACT. RESID. (DLOS1,DLOS2)
c	IT FITS THE FUNCTION : ARCTANH(DLOS)=B+A*LN(T)
c	NEXT VALUE OF T=DEXP(-B/A)
c
      km = 0
      te0 = tei
      if (te0.le.0.d0) te0 = 2.d4
      if (te0.lt.(3.d0*tmin)) te0 = 3.d0*tmin
      tm1 = te0
      tl1 = dlog(tm1)
      do 78 iel = 1, atypes
      do 77 ion = 1, maxion(iel)
      popp0(ion,iel) = pop(ion,iel)
   77 continue
   78 continue
c
      call timion(tm1, edens, hdens, xhyf, tstep)
c
      call cool(tm1, edens, hdens)
c     if (dabs(dlos).lt.(delmin/4.d0)) goto 200
      if (dabs(dabs(dlos)-1.d0).lt.1.d-6)
     & dlos = (1.d0-1.d-6)*dsign(1.0d0,dlos)
      c1 = arctanh(dlos)
      dinc = dimax*(dabs(dlos) ** ex)
      tl2 = tl1-(dinc*dsign(1.d0,dlos))
      t2 = dexp(tl2)
      if (t2.lt.tmin) t2 = tmin
      tl2 = dlog(t2)
      tti = t2
c
c
      dl1 = dlos
      do 100 n = 1, nf
      do 23 iel = 1, atypes
      do 22 ion = 1, maxion(iel)
      pop(ion,iel) = popp0(ion,iel)
   22 continue
   23 continue
c
      call timion(t2, edens, hdens, xhyf, tstep)
c
      call cool(t2, edens, hdens)
      if (dabs(dabs(dlos)-1.d0).lt.1.d-6)
     & dlos = (1.0d0-1.d-6)*dsign(1.d0,dlos)
      c2 = arctanh(dlos)
      csub = c2-c1
      tsub = tl2-tl1
      if (dabs(tsub).lt.1.d-36) tsub = 1.d-36*dsign(1.0d0,tsub)
      if (tsub.eq.0.d0) tsub = 1.d-36*dsign(1.0d0,csub)
      a = csub/tsub
      b = c1-(a*tl1)
      if ((80.d0*dabs(a)).lt.dabs(b)) a =(dabs(b)/80.d0)*dsign(1.0d0,a)
      b = c1-(a*tl1)
      cc = c1
      c1 = c2
      tl1 = tl2
      tm1 = t2
      tl2 =-(b/a)
      t2 = dexp(tl2)
      if ((t2.le.(tm1*emul)).and.(t2.ge.(tm1/emul))) goto 83
      if (t2.gt.(tm1*emul)) t2 = tm1*emul
      if (t2.lt.(tm1/emul)) t2 = tm1/emul
      km = km+1
      tl2 = dlog(t2)
      emul = 1.d0+(cof*(emul-1.d0))
   83 continue
c
c      *PRINT OUT WHEN TESTED BY SUBR. TESTI
c
      if (t2.lt.tmin) goto 200
      if (ntest .ne. 'Y') goto 50
      if (n .ne. 1) goto 30
      write(luop, 43) 
      write(*, 43) 
   43 format(1h ,t4,2hT1,t14,4hDLOS,t27,2hT2,t39,1hA,t49,1hB,t58,2hC1
     &,t68,2hC2)
      write(*, 33) te0, dl1, tti, dinc
      write(luop, 33) te0, dl1, tti, dinc
   33 format(1h ,1pg12.5,g10.3,g12.5,5x,6h(DINC:,g10.3,1h))
   30 write(*, 37) tm1, dlos, t2, a, b, cc, c2, km
      write(luop, 37) tm1, dlos, t2, a, b, cc, c2, km, emul
c
   37 format(1h ,1pg12.5,g10.3,g12.5,4g10.3,i4,1pg10.3)
   50 if ((dabs(dlos).lt.(delmin*5.d0)).or.(km.ge.3)) goto 200
  100 continue
c
      write(*, 120) dlos, t2
  120 format(44h CONVERGENCE FOR EQUIL. TEMP. IS TOO SLOW   ,3hDL:
     &,1pg9.2,5h  TE:,1pg9.3)
  200 continue
      if (t2.lt.tmin) t2 = tmin
      tef = t2
      do 25 iel = 1, atypes
      do 24 ion = 1, maxion(iel)
      pop(ion,iel) = popp0(ion,iel)
   24 continue
   25 continue
c
      call timion(tef, edens, hdens, xhyf, tstep)
c
      call cool(tef, edens, hdens)
c
      return 
      end
