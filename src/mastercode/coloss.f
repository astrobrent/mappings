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
c*******TO FIND COLLISIONAL IONISATION COOLING RATE
c	COLOS (ERG.CM-3.S-1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine coloss(t, de, dh)
c
      include 'cblocks.inc'
c
c
c
      integer*4 jel, jjj
      double precision t, de, dh
c
      colos = 0.0d0
c
      do 10 jel = 1, atypes
         do 20 jjj = 1, maxion(jel)-1
            colos = colos+col(jjj,jel)*epot(jjj,jel)*zion(jel)
     &           *pop(jjj,jel)*de*dh
 20      continue
 10   continue
c
      return 
      end
