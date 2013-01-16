cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******SIMPLE PHOTOIONISATION MODEL
c	SPACE STEPS ARE ALL EQUAL , ON THE SPOT APPROXIMATION
c	VARIABLES ENDING BY "LO" ARE DEXPRESSED IN NEPERIAN LOG
c	VARIABLES ENDING BY "LOG" ARE DEXPRESSED IN BASE 10 LOG
c	VARIABLES CONTAINING "UNI" ARE DEXPRESSED IN RUNITS
c	CALL SUBR.  TEEQUI,TOTPHOT,PHOTSOU,SUMDATA,AVRDATA
c	            WRSPPOP,ZETAEFF,TAUDIST
c
c
c
      subroutine photo1()
c
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision de,dh,dr,dtau,dv,dvoluni,epotmi,fi
      double precision qhlo,rav
      double precision rdist,reclo,rfinal,rin,rout,rstep,rstsun,runit
      double precision t,ti,tstep,vfinalo,vflog,vstrlo
c
      integer*4 i,j,luop,m,nmax,npr,ifil
c
      logical iexi
c
      character imod*4, lmod*4, mmod*4
      character nmod*4,where*16
      character filna*20
c
c           Functions
c
      double precision fdilu
c
c
      luop = 21
      limph = 6
      jcon = 'NO'
      jspot = 'YES'
      where = 'Photo 1'
c
      xhmin = 0.01d0
      write(*, 300) 
c
  300 format(///' PHOTOIONISATION MODEL SELECTED : P1'/
     &' CRUDE MODEL : EQUAL SPACE STEPS'/
     &' ELEMENTS : C,Mg,Si,S,Cl  ARE SINGLY IONISED')
      epotmi = 2.1775d-11
      do 5 i = 1, atypes
      do 6 j = 3, maxion(i)
    6 pop(j,i) = 0.0d0
      pop(1,i) = 1.0d0
      if (epot(1,i).lt.epotmi) pop(1,i) = 0.0d0
      pop(2,i) = 1.0d0-pop(1,i)
    5 continue
      do 1118 i = 3, atypes
      arad(2,i) = dabs(arad(2,i))
      if ((pop(2,i).eq.1.0d0).and.(epot(1,i).lt.epotmi)) arad(2,i)
     & =-arad(2,i)
 1118 continue
c
      call photsou(where)
  810 write(*, 820) 
  820 format(//' GIVE THE PHOTOIONISING SOURCE RADIUS IN ',
     &'SOLAR UNITS : ')
      read(*, *, err=810) rstsun
      if (rstsun.le.0.0) goto 810
c
      rstar = 6.9599d10*rstsun
  830 write(*, 840) 
  840 format(/' GIVE FILLING FACTOR AND PEAK HYDROGEN ','DENSITY : ')
c
      read(*, *, err=830) fi, dh
  843 write(*, 847) 
  847 format(/' APPROXIMATE NUMBER OF SPACE STEPS : ')
      read(*, *, err=843) npr
      if (npr.lt.1) goto 843
c
      nmax = 5*npr
      reclo = dlog(2.6d-13)
      qhlo = dlog(qht)
      qtlog = dlog10(qht)+1.09921d0+2.0d0*dlog10(rstar)
      vstrlo = 2.531d0+qhlo+2.0d0*dlog(rstar)
     &-(2.0*dlog(dh))-dlog(fi)-reclo
      rmax = dexp((vstrlo-1.43241d0)/3.0d0)
      rstep = rmax/npr
      remp = 0.01d0*rmax
      if (remp.lt.(1000.0d0*rstar)) remp=1000.0d0*rstar
      t = teff/2.0d0
      runit = rstep
      vunilog = 3.0d0*dlog10(runit)
      rin = remp
      lmod = 'SO'
      imod = 'ALL'
      mmod = 'TAU'
      nmod = 'EQUI'
      dr = 0.0d0
      dv = 0.0d0
      tstep = 0.0d0
c
c    ***WRITE INITIAL PARAMETERS ON OUPUT FILE
c
c       choose new file name photn****.ph1
c
c
      ifil = 0
 1700 ifil = ifil+1
      filna(1:17) = ' '
      write(filna,1710) ifil
 1710 format(i4.4)
      filna = 'photn'//filna
      filna(10:17)='.ph1'
      inquire(file=filna, exist=iexi)
      if (iexi) goto 1700
      open(luop, file=filna, status='NEW')
c
      write(luop, 2000) 
c
 2000 format(1h ,t40,28hSIMPLE PHOTOIONISATION MODEL,t80,
     &41h(EQUAL SPACE STEPS , ON THE SPOT APPROX.)/t40,
     &28h----------------------------///)
      write(luop, 2800) 
 2800 format(/4h MOD,t7,5hTEMP.,t14,5hALPHA,t23,7hCUT-OFF,t34,5hZSTAR
     &,t45,3hQHI,t54,4hQHEI,t64,5hQHEII)
      write(luop, 2850) iso, teff, alnth, cut, zstar, qhi, qhei
     &, qheii
c
 2850 format(1h ,a2,2(1pg9.2),5(1pg10.3))
      write(*, 2620) 
      write(luop, 2620) 
 2620 format(//7h  RSTAR,t15,4hREMP,t27,8hRSTROMG.,t39,5hRSTEP,t51,
     &10hFILLING F.,t63,10hPEAK DENS.)
      write(*, 2622) rstar, remp, rmax, rstep, fi, dh
      write(luop, 2622) rstar, remp, rmax, rstep, fi, dh
 2622 format(1h ,6(1pg10.3,2x))
c
      write(*, 2610) 
 2610 format(///3h  *,t5,11hTEMPERATURE,t17,5hRDIST,t28,3hFHI,t39,4hFHII
     &,t50,8hEL.DENS.,t61,4hDLOS,t72,5hZETAE)
      write(luop, 2627) 
c
c    ***MODEL COMPUTATION STARTS HERE
c
 2627 format(///3h  *,t5,11hTEMPERATURE,t17,5hRDIST,t28,3hFHI,t39,4hFHII
     &,t50,8hEL.DENS.,t61,4hDLOS,t72,5hZETAE,t83,4hDTAU)
      do 400 m = 1, nmax
      rout = rin+rstep
      rdist = ((m-0.5d0)*rstep)+remp
      rav = (rin+rout)/2.0d0
      wdil = fdilu(rstar,rav)
      dvoluni = 4.18879d0*((rout/runit)**3.d0-(rin/runit)**3.d0)
      ti = t
c
      call totphot(t, dh, fi, rav, dr, dv, wdil, lmod)
      call teequi(ti, t, de, dh, tstep, nmod)
      call sumdata(t, de, dh, fi, dvoluni, rstep, rav, imod)
      call taudist(dh, fi, dtau, rstep, rin, pop, mmod)
      call zetaeff(fi, dh)
c
      write(*, 107) m, t, rdist, pop(1,1), pop(2,1), de, 
     &dlos, zetae,dtau      
c      write(*, 108) (pop(1,jj), pop(2,jj),pop(3,jj)
c     &,pop(4,jj),pop(5,jj),pop(6,jj),jj=1,11)
      write(luop, 107) m, t, rdist, pop(1,1), pop(2,1), de, 
     &dlos, zetae, dtau
c      write(luop, 108) (pop(1,jj), pop(2,jj),pop(3,jj)
c     &,pop(4,jj),pop(5,jj),pop(6,jj),jj=1,11)
c
  107 format(1x,i3,1pf9.0,8(1pg11.3))
c  108 format(6(1pg11.3),/6(1pg11.3),/6(1pg11.3),/6(1pg11.3),/
c     &6(1pg11.3),/6(1pg11.3),/6(1pg11.3),/6(1pg11.3),/
c     &6(1pg11.3),/6(1pg11.3),/6(1pg11.3),/6(1pg11.3))
      rin = rout
      if (pop(2,1).lt.xhmin) goto 444
c
c    ***MODEL ENDED, TRANSFER TO DATA OUTPUT PHASE
c
  400 continue
  444 continue
      rfinal = rdist
      rmax = rfinal
      vfinalo = 1.43241d0+(3.0d0*dlog(rfinal))
c
      vflog = vfinalo/dlog(1.0d1)
  788 continue
      write(luop, 2817) 
 2817 format(///41h MODEL ENDED, PRESENTATION OF THE RESULTS/
     &41h ----------------------------------------////)
      write(luop, 2802) 
 2802 format(4h MOD,t7,5hTEMP.,t14,5hALPHA,t23,7hCUT-OFF,t34,5hZSTAR
     &,t45,3hQHI,t54,4hQHEI,t64,5hQHEII)
c
      write(luop, 2850) iso, teff, alnth, cut, zstar, qhi, qhei
     &, qheii
      write(luop, 2620) 
c
      write(luop, 2622) rstar, remp, rmax, rstep, fi, dh
      write(luop, 2110) 
 2110 format(/8h  RFINAL,t15,8hVFINDLOG,t27,6hFINFHI,t39,7hQTOTLOG)
      write(luop, 2113) rfinal, vflog, xhmin, qtlog
 2113 format(1h ,1pg10.3,2x,0pf10.4,2x,1pg10.3,2x,0pf10.4///)
c
c
c    ***PERFORM AVERAGES OF RELEVANT QUANTITIES
c
      call avrdata
c
c
c    ***USES SUBROUTINE WRSPPOP  TO OUTPUT SPECTRUM,TEMP. AND
c	IONIC POPULATIONS
c
c
      call wrsppop(luop)
c
      close(luop) 
      write(*, 2005) 
 2005 format(//38h0OUTPUT CREATED &&&&&&&&&&&&&&&&&&&&&&)
      return 
      end
