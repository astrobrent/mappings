cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES TOTAL COOLING RATE OF PLASMA , TEMP : T
c	ELECTRON DENSITY : DE , HYDROGEN DENSITY : DH
c	RETURNS TLOSS : TOTAL COOLING RATE (ERG.CM-3.S-1)
c	CALL SUBR.  HYDRO,RESON,INTER,FORBID,FREFRE,COLOSS,
c	            PHEAT,NETGAIN,ALLRATES
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine cool(t, de, dh)
c
      include 'cblocks.inc'
c
c
c
      double precision tll, ell, egg, ett
      character jjmod*4
      double precision t, de, dh
c       
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      jjmod = 'ALL'
      call allrates(t, jjmod)
c
c    ***New hydrogen + helium II line cooling
c    ***Also calls 2 photon and He I calcs
c
      call hydro(t, de, dh)
c
      tll = hloss
c
c
c    ***RESONANCE,INTERCOMBINATION AND FORBIDDEN LINES
c
      call inter(t,de,dh)
      call inter2(t,de,dh)
      call fivelevel(t,de,dh)
      call sixlevel(t,de,dh)
      call ninelevel(t,de,dh)
      call fine3(t,de,dh)
      call ironii(t,de,dh)
      call reson(t,de,dh)
      call reson2(t,de,dh)
      call helif(t,de,dh)
c
      tll = tll+rloss+fslos+floss+f6loss+f9loss+xrloss
c
      tll = tll+xiloss+xheloss+feloss+f3loss
c
c
c    ***FREE-FREE COOLING
c
      call frefre(t, de, dh)
      tll = tll+fflos
c
      cmplos = 0.d0
c
c
c     *** non-relativistic compton cooling/heating ( +ve cooling)
c
c     needs tphot to be up to date
c     Only include compton if completely sure
c
c      call compton(t,de,dh)
c      tll = tll+cmpcool-cmpheat
c
c    ***COLLISIONAL IONISATION LOSSES
c
      call coloss(t, de, dh)
      tll = tll+colos
c
c    ***RECOMB. COOLING (PLUS ESTIMATE OF ON THE SPOT HEATING IF APPL)
c
      call netgain(t, de, dh)
      tll = tll-rngain
c
c    ***PHOTOIONISATION HEATING
c
      call pheat(de, dh)
      tll = tll-pgain
c
c    ***PHOTOELECTRIC GRAIN HEATING AND COLLISIONAL GRAIN LOSSES
c
      call hgrains(t, de, dh)
c
      call hpahs(t,de,dh)
c
      tll = tll-gheat+gcool-paheat
c
c
c    ***CHARGE EXCHANGE HEATING
c
      call cheat(t, de, dh)
      tll = tll-chgain
c
c    ***COSMIC RAY HEATING
c
      call cosmic(t, de, dh)
      tll = tll-cosgain
c
c
c    ***(SEPARATELY) EFFECTIVE LOSS AND GAIN  :  ELOSS,EGAIN
c    Note only include compton heating and cooling if complete sure 
c    Do so at own risk
c
      ell = hloss+rloss+fslos+floss+xrloss+xiloss+xheloss+fflos+colos
     &     +feloss+f3loss+gcool+f6loss+f9loss
c	 +cmpcool
      egg = pgain+cosgain+paheat+gheat
c	  +cmpheat
c
      if (rngain.lt.0.0d0) then
         ell = ell-rngain
      else
         egg = egg+rngain
      end if
c
      if (chgain.lt.0.0d0) then
         ell = ell-chgain
      else
         egg = egg+chgain
      end if
c
      ett = ell+egg
      if (ett.gt.0.d0) then
         dlos = (ell-egg)/ett
      else
         dlos = 1.0d0
      end if
c     
      tloss = tll
      egain = egg
      eloss = ell
c
      if (alphacoolmode.eq.1) then
             tloss = de*de*alphaC0*((1.0d-6*t)**alphaClaw)
      	     egain = 0.d0
             eloss = tloss
      endif
c     
      if (expertmode.gt.0) then
         write(*,*) 'Cool:',tloss,egain,eloss
      endif
c
      return 
      end





