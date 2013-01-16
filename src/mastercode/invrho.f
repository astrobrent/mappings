cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******FUNCTION TO COMPUTE  NUMBER DENSITY (de,dh) USING
c       density in g/cc and population in popul
c	USES FUNCTION DENSTOT (amus)
c
c
      subroutine invrho(rho, de, dh, popul)
c
      include 'const.inc'
c
      double precision rho, de, dh
      double precision dt,f
      double precision popul(mxion,mxelem)
      double precision denstot,feldens
c
      dt = denstot(1.0d0)
      f = feldens(1.0d0,popul)
      dh = rho/(amu*dt+me*f)
      de = f*dh
c
      return 
      end
