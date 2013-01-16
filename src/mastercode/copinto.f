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
c*******TO COPY MATRIX OF IONIC POPULATIONS : POPUL(mxion, mxelem)
c	INTO MATRIX POPOUT(mxion, mxelem)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copinto(popul, popout)
c
      include 'cblocks.inc'
c
      double precision popul(mxion, mxelem), popout(mxion, mxelem)
      integer*4 iel, ion
c
      do 100 iel = 1, atypes
         do 100 ion = 1, maxion(iel)
 100        popout(ion,iel) = popul(ion,iel)
c
      return 
      end
