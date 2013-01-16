cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******CALCULATES POLYNOME VALUE FOR A GIVEN TEFF
c	COEFFICIENTS IN /STAR/
c
c
c
      subroutine fpol(kstr, tm1, i)
c
      include 'cblocks.inc'
c
c
c
      double precision tm1,te4,xt
      character kstr*4
      integer*4 i
c
c
      te4 = teff/1.d4
c
      if (kstr .ne. 'TNO4') goto 60
c
      xt = (tno4(i,10)/(te4*te4))+tno4(i,11)
c
      tm1 = tno4(i,3)+xt*(tno4(i,4)+xt*(tno4(i,5)+xt*(tno4(i,6)+
     &      xt*(tno4(i,7)+xt*(tno4(i,8)+xt*tno4(i,9))))))
c
      goto 200

   60 if (kstr .ne. 'TNO5') goto 100
c
      xt = (tno5(i,10)/(te4*te4))+tno5(i,11)
c
      tm1 = tno5(i,3)+xt*(tno5(i,4)+xt*(tno5(i,5)+xt*(tno5(i,6)+
     &      xt*(tno5(i,7)+xt*(tno5(i,8)+xt*tno5(i,9))))))
c
      goto 200
c
  100 if (kstr .ne. 'TMET') goto 200
c
      xt = (tmet(i,10)/(te4*te4))+tmet(i,11)
c
      tm1 = tmet(i,3)+xt*(tmet(i,4)+xt*(tmet(i,5)+xt*(tmet(i,6)+
     &      xt*(tmet(i,7)+xt*(tmet(i,8)+xt*tmet(i,9))))))
c
 200  return
c
      end
