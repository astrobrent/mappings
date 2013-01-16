cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	COMPUTES COLLISIONAL IONISATION RATES FOR
c	ALL IONS ; RATES RETURNED AS VECTOR COL(6,11)
c	IN COMMON BLOCK /COLION/
c
c	This is a new routine to calculate collisional ionisation
c	rates using an algorithm based on Arnaud and Rothenflug 1985
c	and Younger 1981,1983.  It uses a five parameter fit to the 
c	collision crossection and derives an expresion for the integral
c	over a maxwellian velocity distrbution in electron energies.
c	
c	Numerical procedures sugested by A&R have been revised,
c	and extended to more rigorous precision.
c       The exponential integrals
c	are evaluated in the function farint.  The basis for the second
c	integral comes from Hummer 1983 and earlier references thererin.
c	The Hummer expression is replaced when x is small with a 
c       more accurate convergent series.
c
c	RSS 7/90
c
c	added selector for different calculation methods
c	so the user can choose by setting the mode in ATDAT
c
c	RSS 8/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine collion(telec)
c
      include 'cblocks.inc'
c
c
      double precision telec, f,f0, f1, f2, a, b, c, d
c      double precision ec, bz, r, g, en, e21, t, fo1
      double precision ktev, x, ac, tc
      double precision cdi,cai,fautoi
      integer*4 ion ,atom, shell,is
c
      double precision t,sig
c
c     External Functions
c
      double precision farint
c
      if (mapmode.eq.1) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	original mappings code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c         ec = 1.602193d-12
c         bz = 1.38054d-16
c         r = 1.04d-7
c         t = telec+1.d-1
c         f = dsqrt(t)
c         do 10 atom = 1, atypes
c            do 20 ion = 1, maxion(atom)-1
c               a = epot(ion,atom)/(bz*t)
c               b = beta(ion,atom)
c               c = a*(b-1.0d0)
c               fo1 = fue1(c,0.0d0)
c               g = (1.0d0+a*(b-2.0d0))-((a*(((c*(b-2.0d0))+(2.0d0*b))
c     &             -3.0d0))*fo1)
c               g = (g*c)/b
c               en = effn(ion,atom)
c               e21 = (ec/epot(ion,atom)) ** 2
c               cdi= ((((r*en)*f)*e21)*g)*dexp(- a)
c               if (cdi.lt.0.0d0) cdi= 0.0d0
c               col(ion,atom) = cdi
c 20         continue
c            col(maxion(atom),atom) = 0.0d0
c 10      continue
c
      else
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New mappings code, A&R or L&M colliosonal ionisation
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
         if (collmode.eq.1) then
c     
c	Landini & Monsignori Fossi - Shull rates
c	faster but probably less accurate.....
c       will use A&R for H sequence.
c
            t = telec+1.d-1
            f = dsqrt(t)
            ktev = t*rkb/ev
            do 30 atom = 1, atypes
               do 40 ion = 1, maxion(atom)-1
                  
                  is = mapz(atom)-ion+1
                  if (is.gt.1) then 
                     ac = acol(ion,atom)
                     tc = tcol(ion,atom)
                     cdi = ac*f/(1.d0+0.1d0*(t/tc))*dexp(-tc/t)
                  else
c
c     Use A&R  for Hydrogenic sequence
c
                     a = aar(1,ion,atom)
                     b = bar(1,ion,atom)
                     c = car(1,ion,atom)
                     d = dar(1,ion,atom)
                     x = ipot(1,ion,atom)/ktev
c
c     to prevent numerical errors if x is too large...
c
                     if (x.ge.200.d0) x = 200.d0
c
                     f1 = farint(1,x)
                     f2 = farint(2,x)
                     f0 = a*(1.d0-x*f1)+b*(1.d0+x-x*(2.d0+x)*f1)+c*f1
     &                    +d*x*f2
                     sig = (dexp(-x)/x)*f0
                     cdi = ((6.69d-7)/(ktev**1.5d0))*sig
c
                  endif
                  col(ion,atom)= cdi
 40            continue
               col(maxion(atom),atom) = 0.0d0
 30         continue
         end if
c     
         if (collmode.eq.0) then
c     
c 	Arnaud and Rothenflug rates (default)
c
            ktev = telec*rkb/ev
c     
c     for each type of atom
c
            do 110 atom = 1, atypes
               do 100 ion = 1, maxion(atom)-1
                  sig = 0.d0
                  do 200 shell = 1 , nshells(ion,atom)
                     a = aar(shell,ion,atom)
                     b = bar(shell,ion,atom)
                     c = car(shell,ion,atom)
                     d = dar(shell,ion,atom)
                     x = ipot(shell,ion,atom)/ktev
c
c     to prevent numerical errors if x is too large...
c
                     if (x.ge.200.d0) x = 200.d0
c     
                     f1 = farint(1,x)
                     f2 = farint(2,x)
                     f0 = a*(1.d0-x*f1)+b*(1.d0+x-x*(2.d0+x)*f1)+c*f1
     &                   +d*x*f2
                     sig = sig+(dexp(-x)/x)*f0
c
 200              continue
c
c     direct ionisation
c
                  cdi = ((6.69d-7)/(ktev**1.5d0))*sig
c
c     auto ionisation
c
                  cai = fautoi(telec,atom,ion)
c                  cai = 0.d0
c
c
c                  if (atom.eq.1) then
c                  write(*,*) 'col:',telec,' ',elem(atom),rom(ion),
c     &                 cdi,cai
c                  endif
c
                  col(ion,atom) = cdi+cai
c
                  if (col(ion,atom).lt.0.d0) then
                  col(ion,atom) = 0.d0
                  endif
c
c
 100           continue
               col(maxion(atom),atom) = 0.0d0
 110        continue
c
         endif
c
c     end new mappings code
c
      endif
c
      return
c 
      end
