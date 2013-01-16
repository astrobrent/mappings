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
c*******FINDS FINAL IONIC ABUNDANCES AFTER TIME :TSTEP
c	FOR ALL ELEMENTS ; OUTPUT ELECTRONIC DENSITY : DE
c	AND THE FRACTIONAL ABUNDANCE OF THE DIFFERENT SPECIES
c	OF EACH ELEMENT  :  POP(6,11)  IN COMMON BLOCK /ELABU/
c	FHIIF : FINAL FRACTIONAL IONISATION OF HYDROGEN
c	CALL SUBROUTINES IOHYD,IOBAL,COPINTO,DIFPOP
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine timion(t, de, dh, fhiif, tstep)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision t, de, dh, fhiif, tstep, dt
      double precision dcopl,dif,difh,dift,dq,fhii,fhii0,fhiiav
      double precision trea,xhy
      double precision popp0(mxion, mxelem), popp2(mxion, mxelem)
c
      integer*4 i,num,numax
c
      character mod*4, nel*4
c
c           Functions
c
      double precision feldens
c
      trea = 1.d-7
      mod = 'TIM'
      nel = 'ALL'
      xhy = -1.d0
      numax = 8
      fhii0 = pop(2,1)
      de0 = feldens(dh,pop)
      dt = tstep
c
c    ***EVOLVES IONIC POPULATIONS IN ONE STEP AND CHECKS
c	THE AMOUNT OF CHANGE 
c
      call copinto(pop, popp0)
      call iohyd(dh, xhy, t, dt, de, fhii,mod)
      fhiiav = (0.3d0*fhii0)+(0.7d0*fhii)
      call iobal(mod, nel, de, dh, fhiiav, t, dt)
      call copinto(pop, popp2)
      pop(1,1) = popp0(1,1)
      pop(2,1) = popp0(2,1)
c
c
      call iohyd(dh, xhy, t, dt, de, fhii,mod)
      call difpop(pop, popp0, trea, atypes, dift)
      call difpop(pop, popp0, trea, 1, difh)
      call difpop(popp2, pop, trea, 1, dcopl)
      dq = dabs(dlog(de+epsilon)-dlog(de0+epsilon))
      dif = dmax1(dq,100.d0*dcopl,0.6d0*dsqrt(difh),
     &0.4d0*dsqrt(dift))
c
c	WRITE(6,10) DQ,DCOPL,DIFH,DIFT,DIF,NUM
c	WRITE(25,10) DQ,DCOPL,DIFH,DIFT,DIF,NUM
c10	FORMAT (/' ',5(0PF8.4),0PI4)
c	WRITE (6,*) NUM,DT,POP(1,1),POP(2,1)
c
c      *DIVIDES TIMESTEP INTO NUM MINI TIMESTEPS (IF NUM > 1)
c      *DERIVES ABUNDANCES AFTER EACH DT , NUM TIMES
c
      num = min0(1+idint((4.0d0*dif)+0.4d0),numax)
      if (num.le.1) goto 200
c
      call copinto(popp0, pop)
c
      dt = tstep/(dble(num))
      do 100 i = 1, num
      fhii0 = pop(2,1)
      call iohyd(dh, xhy, t, dt, de, fhii,mod)
      fhiiav = (0.3d0*fhii0)+(0.7d0*fhii)
c	WRITE (6,*) DT,POP(1,1),POP(2,1),POPP0(1,1),POPP0(2,1)
      call iobal(mod, nel, de, dh, fhiiav, t, dt)
  100 continue
  200 continue
c
      fhiif = fhii
      return 
      end
