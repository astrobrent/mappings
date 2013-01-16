cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******SIMPLE PHOTOIONISATION MODEL
c	SPACE STEPS DERIVED FROM DTAU
c	VARIABLES ENDING BY "LO" ARE DEXPRESSED IN NEPERIAN LOG
c	VARIABLES ENDING BY "LOG" ARE DEXPRESSED IN BASE 10 LOG
c	VARIABLES CONTAINING "UNI" ARE DEXPRESSED IN RUNITS
c	CALL SUBR.  TEEQUI,TOTPHOT,PHOTSOU,SUMDATA,AVRDATA
c	            WRSPPOP,ZETAEFF,TAUDIST
c
c
c
      subroutine photo2()
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision chma,chma0,de,dh,dr,dtau,dtc,dv
      double precision dvoluni,epotmi,fi
      double precision qhlo,rav,rdist,reclo,rfinal,rin,rnext
      double precision rout,rpre,rstep,rstsun,runit,t,ti0,tstep
      double precision vfinalo,vflog,vstrlo
c
      integer*4 i,j,luop,m,nmax,npr
      integer*4 ifil
c
      character imod*4, lmod*4, mmod*4
      character nmod*4, filna*20, where*16
c
      logical iexi
c
c           Functions
c
      double precision fcrit,fdilu
c
c
      luop = 21
c
      limph = 6
      jcon = 'NO'
      jspot = 'YES'
      where = 'Photo 2'
      xhmin = 0.01
      chma0 = 3.0
c
c
      npr = 20
c
      write(*, 300) 
c
c
  300 format(///36h PHOTOIONISATION MODEL SELECTED : P2/
     &34h CRUDE MODEL : ON THE SPOT APPROX./
     &44h ELEMENTS : C,MG,SI,S,CL  ARE SINGLY IONISED)
c
      epotmi = 2.1775d-11
      do 5 i = 1, atypes
      do 6 j = 3, maxion(i)
    6 pop(j,i) = 0.0
      pop(1,i) = 1.0
      if (epot(1,i).lt.epotmi) pop(1,i) = 0.0
      pop(2,i) = 1.0-pop(1,i)
    5 continue
      do 1118 i = 3, atypes
      arad(2,i) = dabs(arad(2,i))
      if ((pop(2,i) .eq. 1.0).and.(epot(1,i).lt.epotmi)) arad(2,i)
     & =-arad(2,i)
c
c    ***INPUT OF PHYSICAL CONDITIONS
c
c
 1118 continue
c
c
      call photsou(where)
c
  810 write(*, 820) 
  820 format(//38h GIVE THE PHOTOIONISING SOURCE RADIUS ,
     &14h(RSUN=6.96E10)/34h$(IN SOLAR UNITS (<1.E8) OR IN CM ,
     &10h(>1.E8) : )
      read(*, *, err=810) rstsun
      if (rstsun.le.0.0) goto 810
      if (rstsun.lt.1.d8) then
      rstar = 6.96d10*rstsun
      else
      rstar = rstsun
      rstsun = rstar/6.96d10
c
c
      end if
c
  830 write(*, 840) 
  840 format(/39h$GIVE FILLING FACTOR AND PEAK HYDROGEN ,10hDENSITY : )
c
c
      read(*, *, err=830) fi, dh
c
  843 write(*, 847) 
  847 format(/44h$STEP VALUE OF THE OPTICAL DEPTH (<=0.05) : )
      read(*, *, err=843) dtau
      if ((dtau.le.0.0).or.(dtau.gt.0.05)) goto 843
      dtau0 = dtau
c
c
      nmax = 20*npr
c
      reclo = dlog(2.6d-13)
      qhlo = dlog(qht)
      qtlog = (dlog10(qht)+1.09921)+(2.0*dlog10(rstar))
      vstrlo = ((((2.531+qhlo)+(2.0*dlog(rstar)))-(2.0*dlog(dh
     &)))-dlog(fi))-reclo
c
c
      rmax = dexp((vstrlo-1.43241)/3.0)
c
  310 write(*, 320) rmax
  320 format(/36h HYDROGEN STROMGREN SPHERE RADIUS : ,1pg10.3/
     &27h GIVE RADIUS OF EMPTY ZONE /
     &54h$(IN CM (>1.D6) OR IN FRACTION OF STROMGREN RADIUS) : )
      read(*, *, err=310) remp
      if (remp.lt.0.0) goto 310
      if (remp.le.1.d6) remp = remp*rmax
      if (remp.lt.rstar) remp = rstar
c
c
      rmax = dexp(((vstrlo+dlog(1.d0+((remp/rmax) ** 3)))-
     &1.43241)/3.0)
      rstep = dmax1(rmax-remp,1.d-6*rmax)/npr
      rnext = rstep
      runit = rstep
      vunilog = 3.0*dlog10(runit)
      rin = remp
      lmod = 'SO'
      imod = 'ALL'
      mmod = 'DIS'
      nmod = 'EQUI'
      dr = 0.0d0
      dv = 0.0d0
      tstep = 0.0d0
      ti0 = 0.0d0
c
c    ***WRITE INITIAL PARAMETERS ON OUPUT FILE
c

c
c       choose new file name photn****.ph2
c
c
      ifil = 0
 1700 ifil = ifil+1
      filna(1:17) = ' '
      write(filna,1710) ifil
 1710 format(i4.4)
      filna = 'photn'//filna
      filna(10:17)='.ph2'
      inquire(file=filna, exist=iexi)
      if (iexi) goto 1700
      open(luop, file=filna, status='NEW')
c
      write(luop, 2000) 
c
c
 2000 format(1h ,t40,28hSIMPLE PHOTOIONISATION MODEL,t80,
     &21h(ON THE SPOT APPROX.)/t40,28h----------------------------///)
c
      write(luop, 2800) 
 2800 format(/4h MOD,t7,5hTEMP.,t16,5hALPHA,t22,7hTURN-ON,t30,7hCUT-OFF
     &,t38,5hZSTAR,t48,3hQHI,t57,4hQHEI,t67,5hQHEII)
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii
c
c
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
c
      write(*, 2620) 
      write(luop, 2620) 
 2620 format(//7h  RSTAR,t15,4hREMP,t27,8hRSTROMG.,t39,4hDTAU,t51,
     &10hFILLING F.,t63,10hPEAK DENS.)
      write(*, 2622) rstar, remp, rmax, dtau, fi, dh
      write(luop, 2622) rstar, remp, rmax, dtau, fi, dh
c
c
 2622 format(1h ,6(1pg10.3,2x))
c
      write(*, 2610) 
 2610 format(///3h  *,t5,11hTEMPERATURE,t17,5hRDIST,t28,3hFHI,t39,4hFHII
     &,t50,8hEL.DENS.,t61,4hDLOS,t72,5hRSTEP)
      write(luop, 2627) 
c
c    ***MODEL COMPUTATION STARTS HERE
c
c
 2627 format(///3h  *,t5,11hTEMPERATURE,t17,5hRDIST,t28,3hFHI,t39,4hFHII
     &,t50,8hEL.DENS.,t61,4hDLOS,t72,5hRSTEP,t83,5hZETAE,t94,4hQHDH)
c
      do 400 m = 1, nmax
      rstep = dmin1(rnext,chma0*rstep)
  450 continue
      rout = rin+rstep
      rav = (rin+rout)/2.0
      rdist = rav
      wdil = fdilu(rstar,rav)
c
c
c
      dvoluni = 4.18879*(((rout/runit) ** 3)-(dble(rin/runit
     &) ** 3))
c
      call totphot(ti0, dh, fi, rav, dr, dv, wdil, lmod)
      call teequi(ti0, t, de, dh, tstep, nmod)
c
c
      call taudist(dh, fi, dtau, rnext, rin, pop, mmod)
c
      if (m .eq. 1) then
      if (rstep.gt.(chma0*rnext)) then
      rstep = chma0*rnext
      goto 450
      else if (rstep.lt.(rnext/5.0)) then
      rstep = rnext/5.0
      goto 450
      end if
      else
      dtc = dabs(t-ti0)/dmax1(t,ti0)
      chma = chma0*(1.d0-(0.6d0*(1.d0-fcrit(dtc,1.1d2))))
      if (rstep.gt.(chma*rpre)) then
      rstep = chma*rpre
      goto 450
      end if
      end if
      rpre = rstep
c
c
      ti0 = t
c
      call sumdata(t, de, dh, fi, dvoluni, rstep, rav, imod)
c
c
c
      call zetaeff(fi, dh)
c
      write(*, 107) m, t, rdist, pop(1,1), pop(2,1), de, 
     &dlos, rstep
c
      write(luop, 107) m, t, rdist, pop(1,1), pop(2,1), de, 
     &dlos, rstep, zetae, qhdh
c
c
  107 format(1x,i3,0pf9.0,9(1pg11.3))
c
      rin = rout
      if ((pop(2,1).lt.xhmin).and.(t.lt.200.0)) goto 444
c
c    ***MODEL ENDED, TRANSFER TO DATA OUTPUT PHASE
c
c
  400 continue
c
  444 continue
      rfinal = rdist
      rmax = rfinal
      vfinalo = 1.43241+(3.0*dlog(rfinal))
c
c
      vflog = vfinalo/dlog(1.0d1)
c
  788 continue
      write(luop, 2817) 
 2817 format(///41h MODEL ENDED, PRESENTATION OF THE RESULTS/
     &41h ----------------------------------------////)
      write(luop, 2800) 
c
c
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii
c
      write(luop, 2620) 
c
c
      write(luop, 2622) rstar, remp, rmax, dtau, fi, dh
c
      write(luop, 2110) 
 2110 format(/8h  RFINAL,t15,8hVFINDLOG,t27,6hFINFHI,t39,7hQTOTLOG)
c
      write(luop, 2113) rfinal, vflog, xhmin, qtlog
c
c
c    ***PERFORM AVERAGES OF RELEVANT QUANTITIES
c
c
 2113 format(1h ,1pg10.3,2x,0pf10.4,2x,1pg10.3,2x,0pf10.4///)
c
c
c    ***USES SUBROUTINE WRSPPOP  TO OUTPUT SPECTRUM,TEMP. AND
c	IONIC POPULATIONS
c
c
      call avrdata
c
c
c
      call wrsppop(luop)
c
      close(luop) 
      write(*, 2005) 
 2005 format(//38h0OUTPUT CREATED &&&&&&&&&&&&&&&&&&&&&&)
      return 
      end
