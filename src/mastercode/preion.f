cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO DERIVE PREIONISATION CONDITIONS IN SHOCKS
c	OUTPUT IONISATION STATES IN POP(6,11) /ELABU/
c	FI : FILLING FACTOR , LTERM : OUTPUT DEVICE NUMBER
c	DH : PEAK H DENSITY , SHOCK VELOCITY : VS
c	TSMAX : MAXIMUM TIME STEP ALLOWED FOR IONISATION
c	USES PHOTON FIELD IN ARRAY : TPHOT
c	OUPUT TEMPERATURE : TEF ,  EL.DENS : DEF ,  AND
c	NUMBER OF PHOTOIONISING PHOTONS : QDIFT
c	CALL SUBR. INTVEC,TEEQUI,TAUDIST
c
c	NB.   SUBR. TOTPHOT NEEDS TO BE PREVIOUSLY CALLED
c
c
      subroutine preion(lterm, luop, tsmax, vs, fi, dh, def, tef, qtot
     &                 , drta)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision ct,dch,def,dh,dq,dq0,dq1,dqmi
      double precision dqmm,dqp,dqpp,dqva,drta,dt,fhi,fi
      double precision fma,frmm,frpp,frta,frti,frtma,frtp
      double precision q1,q2,q3,qavail,qtot,qused,rad,t,tef,tmi
      double precision ts0,ts1,tsm,tsmax,tstep,u,vs,wei
      double precision xhyf,zrr,ztr
c
      integer*4 i,j,l,lterm,luop,m,nf
c
      character nmod*4, mmod*4
      character jsp*4, lpm*4
c
c           Functions
c
      double precision ff,fr
c
      ff(u,ct) = dmin1(0.99999d0,(ct*u)+1.d-5)
      fr(u,ct) = (1.d0-((1.d0-ff(u,ct)) ** ct))/ff(u,ct)
c
      rad = 0.0d0
      ct = 0.6321d0
      frtma = 1.0d1/ct
      jsp = jspot
      jspot = 'YES'
      tem = -1.0d0
      nf = 16
      dqmi = 1.0d-2
      tmi = 0.075d0
      tsm = tsmax
      if (tsm.le.0.d0) tsm = 1.d36
      if (tef.lt.100.d0) then
      tef = 1.d4
      lpm = 'FIXT'
      else
      lpm = 'EQ'
      end if
c
c
c    ***SET POPI(6,11) TO PROTO-IONIC POULATIONS
c	AND FINDS NUMBER OF PHOTONS AVAILABLE : QTOT
c
      call copinto(propop, popi)
c
      call intvec(tphot, q1, q2, q3, qtot)
      frti = 0.07d0
      mmod = 'DIS'
      nmod = 'TIM'
      l = 0
      if (lterm.gt.0) write(lterm, 510) dh, vs, qtot, fi
      if (luop.gt.0) write(luop, 510) dh, vs, qtot, fi
  510 format(/33h HDENS.,VSHOC, QTOT AND FIL.F. : ,4(1pg11.3))
      if (lterm.gt.0) write(lterm, 520) 
      if (luop.gt.0) write(luop, 520) 
  520 format(/9h #   <TE>,t15,6hQAVAIL,t25,5hQUSED,t34
     &,3hFHI,t43,5hTSTEP,t52,4hDRTA,t62,3hTAU,t72,2hDQ/)
c
c    ***ITERATES ON OPTICAL DEPTH : FRTA   UNTIL THE NUMBER OF
c	ABSORBED PHOTONS : QUSED  IS EQUAL TO THE NUMBER OF
c	PHOTONS AVAILABLE : QAVAIL
c
      ztr = 0.5d0
      dqpp = 2.d0
      dqmm = -2.d0
      do 300 m = 1, nf
      if (m .eq. 1) then
      frta = frti
      frtp = frta
      frmm = 1.d38
      frpp = 0.d0
      else
      zrr = dmin1(ztr,dmax1(0.03d0,ztr*(dt/tmi)))
      if ((dq1.gt.zrr).and.(dq1.lt.dqpp)) then
      dqpp = dq1
      frpp = frta
      end if
      if ((dq1.lt.(- zrr)).and.(dq1.gt.dqmm)) then
      dqmm = dq1
      frmm = frta
      end if
c
      fma = 5.d0
      if (((((dabs(dq1)-dabs(dqp))/((dabs(dq1)+dabs(dqp))+1.d-36))
     &.lt.0.05d0).and.(dabs(dq1).gt.0.5d0)).and.(m.gt.2)) fma = 
     &12.0d0
      if (((dqp/(dq1+epsilon)).lt.0.d0).and.(m.gt.2)) then
      wei = 0.5d0*(((dabs(dqp)/((dabs(dqp)+dabs(dq1))+1.d-36))/0.5d0)
     & ** 0.75d0)
      dch = dexp((wei*dlog(frta))+((1.d0-wei)*dlog(frtp)))
      frtp = frta
      frta = dch
      else
      dch = (qavail/qused)*frta
      if (dch.gt.frmm) then
      wei = 0.5d0*(((dabs(dqmm)/((dabs(dqmm)+dabs(dq1))+1.d-36))/
     &0.5d0) ** 0.75d0)
      dch = dexp((wei*dlog(frta))+((1.0d0-wei)*dlog(frmm)))
      else if (dch.lt.frpp) then
      wei = 0.5d0*(((dabs(dqpp)/((dabs(dqpp)+dabs(dq1))+1.d-36))/
     &0.5d0) ** 0.75d0)
      dch = dexp((wei*dlog(frta))+((1.0d0-wei)*dlog(frpp)))
      end if
      frtp = frta
      frta = dmin1(fma*frta,dmax1(frta/fma,dch))
      end if
      end if
      dqp = dq1
c
      call copinto(popi, popf)
c
      do 305 l = 1, 10
      u = frta
      wei = fr(u,ct)
      call averinto(wei, popf, popi, pop)
      call taudist(dh, fi, frta, drta, rad, pop, mmod)
      call copinto(popi, pop)
      tstep = dmin1(drta/vs,tsm)
      if (l.gt.1) tstep = dmax1(ts1/2.d0,dmin1(2.d0*ts1,tstep,tsm))
      ts0 = ts1
      ts1 = tstep
      dt = dabs(ts1-ts0)/(ts1+ts0)
c
      if (lpm .eq. 'EQ') then
      call teequi(tef, t, def, dh, tstep, nmod)
      else
      t = tef
      call timion(t, def, dh, xhyf, tstep)
      end if
c
      call copinto(pop, popf)
c
      qused = 0.d0
      do 411 i = 1, atypes
      do 410 j = 1, maxion(i)
      qused = qused+((zion(i)*0.5d0)*dabs(popf(j,i)-popi(j,i)))
  410 continue
  411 continue
c
      qused = (((dh*fi)*drta)*qused)*fr(u,ct)
      qavail = qtot*tstep
      dq = (qavail-qused)/((qavail+qused)+1.d-36)
      dq0 = dq1
      dq1 = dq
      dqva = dabs(dq0-dq1)
      fhi = popf(1,1)
      tef = t
c
      if (lterm.gt.0) write(lterm, 230) m, tef, qavail, qused
     &, fhi, tstep, drta, frta, dq
      if (luop.gt.0) write(luop, 230) m, tef, qavail, qused, 
     &fhi, tstep, drta, frta, dq
c
  230 format(1h ,i2,f9.0,2(1pg10.3),3(1pg9.2),1pg10.3,0pf9.4)
      if ((((dabs(dq).lt.dqmi).or.(frta.gt.frtma)).or.((tstep
     &.ge.tsm).and.(dq.ge.dqmi))).and.(l.gt.1)) goto 308
      if (((l.gt.1).and.(dt.lt.tmi)).and.(dabs(dq).gt.(2.0*
     &tmi))) goto 300
      if ((l.gt.1).and.(dqva.lt.(dqmi/2.0))) goto 300
      if ((l.gt.3).and.(dt.lt.(2.0*tmi))) goto 300
  305 continue
  300 continue
c
      if (dabs(dq).gt.(4.0*dqmi)) then
      if (luop.gt.0) write(luop, 210) 
      if (lterm.gt.0) write(lterm, 210) 
  210 format(/45h$$$$$$$$$$$$$$$$ CONVERGENCE WAS NOT REACHED ,
     &37hDURING DETERMINATION OF PREIONISATION/)
      end if
c
  308 if (tstep.lt.tsm) goto 903
      if (luop.gt.0) write(luop, 417) tsm
      if (lterm.gt.0) write(lterm, 417) tsm
c
  417 format(/46h [[[[[[[[[[[[[[[ TIME STEP LIMITED BY MAXIMUM ,
     &7hVALUE :,1pg10.3)
  903 continue
      jspot = jsp
c
      tem = -1.d0
c
      return 
      end

