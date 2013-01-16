cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO OUTPUT AVERAGE OF POPW & POPCO IN POPOUT
c	THE WEIGHT OFF POPW IS : WEI (0<=WEI<=1)
c
c
      subroutine averinto(wei, popw, popco, popout)
c     
      include 'cblocks.inc'
c     
      double precision popw(mxion, mxelem), popco(mxion, mxelem) 
     &                ,popout(mxion, mxelem)
      double precision wei, weiw, weico
c     
      integer*4 i, j
c     
      weiw = dmax1(0.0d0,dmin1(1.0d0,wei))
      weico = dmax1(0.d0,1.d0-weiw)
c     
      do 100 i = 1, atypes
         do 100 j = 1, maxion(i)
            popout(j,i) = (weiw*popw(j,i))+(weico*popco(j,i))
 100  continue
c
      return 
      end
