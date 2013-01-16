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
c*******TO FIND IONISATION EQUILIBRIUM AT A GIVEN TEMP.& DENS.
c	FOR ALL ELEMENTS  ;  OUTPUT ELECTRONIC DENSITY  :  DE
c	AND THE FRACTIONAL ABUNDANCE OF THE DIFFERENT SPECIES
c	OF EACH ELEMENT  :  POP(6,11)  IN COMMON BLOCK /ELABU/
c	CALL SUBROUTINES IOHYD,IOBAL
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine equion(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision popzero(mxion, mxelem)
      double precision t, de, dh, tstep, xhy, difma
      double precision difm1, difm2
      double precision treh, trea, zm, xhyf, dif1
      double precision dia, dift, dfh, dft
      integer*4 nf, n, m, iel,inttemp
c
      character mod*4, nff*4, nel*4
c
      tstep = 1.d37
      mod = 'EQUI'
      xhy = -1.0d0
      nf = 6
      difma = 1.d-3
      difm1 = 1.d-4
      difm2 = 0.08d0
      treh = 1.d-8
c
c
      trea = 1.d-5
      nff = 'ALL'
      zm = 0.0d0
      do 30 iel = 3, atypes
 30      zm = zm+zion(iel)
c
c
c    ***ITERATES TO FIND IONIC POPULATIONS AT EQUILIBRIUM
c
      do 222 m = 1, 4
         call copinto(pop, popzero)
         do 10 n = 1, nf
c
            call iohyd(dh, xhy, t, tstep, de, xhyf,mod)
c
            call difpop(pop, popzero, treh, 1, dif1)
c
            dia = dmax1(0.0d0,dlog(5.0d0/(zm+1.d-5))/20.0d0)*pop(2,1)
            if (dif1.ge.dia) then
               nel = 'ALL'
            else if ((m .eq. 1).and.(n.lt.nf)) then
               nel = 'HE'
            else if (m.gt.2) then
               nel = 'ALL'
            end if
c
            if (de.lt.1.d-6) de = 1.d-6
            call iobal(mod, nel, de, dh, xhyf, t, tstep)
            call difpop(pop, popzero, trea, atypes, dift)
c
            call copinto(pop, popzero)
c
            if ((dift.lt.difm1).and.(nel .eq. 'HE')) goto 80
            if ((dift.lt.difm2).and.(nel .eq. 'ALL')) goto 80
c
 10      continue
c
 80      if (((nel .eq. 'ALL').and.(dif1.lt.difma)).and.(dift.lt.
     &        difm2)) goto 300
         call iohyd(dh, xhy, t, tstep, de, xhyf,mod)
         call iobal(mod, nff, de, dh, xhyf, t, tstep)
         inttemp = 1
         call difpop(pop, popzero, treh, inttemp, dif1)
c
         call difpop(pop, popzero, trea, atypes, dift)
         if ((dif1.le.difma).and.(dift.le.difm2)) goto 300
c
 222  continue
c
      if ((dif1.le.difma).and.(dift.le.difm2)) goto 300
      dfh = dif1/difma
      dft = dift/difm2
      write(*, 56) dfh, dft
 56   format('Convergence for Equil. ionisation too slow :','DFH:'
     &     ,1pg9.2,'   DFT:',1pg9.2)
c
 300  continue
      return 
      end


