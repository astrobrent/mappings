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
c     Subroutine to calculate the transitions
c     from metatsable levels of highly ionised ions.
c     long lambda and low ionisations are handled in inter.
c
c	Uses the collision strength/statwei provided
c
c
c	refs: Landini & Monsignori Fosse 1990 A.A.Suppl.Ser. 82:229
c		  Mewe 1985 A.A.Suppl.Ser. 45:11
c
c	RSS 9/90
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine inter2(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision t, de, dh
      double precision egj,ejk,rr
      double precision cgjb,omegab
      double precision abup,y
c
      integer*4 i,atom,ion
c
      xiloss = 0.0d0
      t = dmax1(0.1d0,t)
c
c
      do i = 1,xilines
c     
c	get data from arrays for line # i
c     
c     Note: only line for possible species are read in
c     so no need to check maxion etc...
c     
         atom = xiat(i)
         ion = xiion(i)
c     
c     only calculate for abundant species
c     
         xrbri(i) = 0.d0
c     
         abup = zion(atom)*pop(ion,atom)
c
         if (abup.ge.pzlimit) then
c     
c     note that Egj does equal Ejk
c     
            egj = xiejk(i)
            ejk = xiejk(i)
            omegab = xomeg(i)
c     
c     get scaled energy gap to level j from ground (not k)
c     
            y = egj*ev/(rkb*t)
c     
c     transition power rate rr
c     
            rr = 0.d0
c     
            cgjb = 0.d0
            if (y.lt.loghuge) cgjb = rka*(dexp(-y)/dsqrt(t))*omegab
c     
            rr = ejk*ev*cgjb
c     
c     number to transition abup
c 
c     
c     total power in line
c     
            rr = de*dh*rr*abup
c     
            xiloss = xiloss+rr
c
            xibri(i) = rr/(4.0d0*pi)
c     
c     end population limited loop
c
         endif
c     
c     
      enddo
c     
c     
      return 
c
      end
