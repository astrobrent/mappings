cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Simple 4th order RK integration step for 3 variables
c	  in spherical flow models.  Uses shpderivs
c      x = ???  y1 = r y2 = v y3 = P
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rk4 (x, y1, y2, y3, dy1dx, dy2dx, dy3dx, h,
     &     yout1,yout2,yout3)
c
      include 'const.inc'
c     
c     parameters
c     
      real*8  x, y1, y2, y3, dy1dx, dy2dx, dy3dx
      real*8  xh,h,yout1,yout2,yout3
c     
c     local vars
c     
      real*8 h6, hh
      real*8 dym1, dyt1, yt1
      real*8 dym2, dyt2, yt2
      real*8 dym3, dyt3, yt3
c     
      hh = h * 0.5
      h6 = h / 6.0

      xh = x + hh

      yt1 = y1 + hh * dy1dx
      yt2 = y2 + hh * dy2dx
      yt3 = y3 + hh * dy3dx

      call sphderivs(xh, yt1, yt2, yt3, dyt1, dyt2, dyt3)

      yt1 = y1 + hh * dyt1
      yt2 = y2 + hh * dyt2
      yt3 = y3 + hh * dyt3

      call sphderivs(xh, yt1, yt2, yt3, dym1, dym2, dym3)

      yt1 = y1 + h * dym1
      yt2 = y2 + h * dym2
      yt3 = y3 + h * dym3

      dym1 = dyt1 + dym1
      dym2 = dyt2 + dym2
      dym3 = dyt3 + dym3

      call sphderivs(x + h, yt1, yt2, yt3, dyt1, dyt2, dyt3)
c
c     put answers back into input vars
c
      yout1 = y1 + h6 * (dy1dx + dyt1 + 2.0 * dym1)
      yout2 = y2 + h6 * (dy2dx + dyt2 + 2.0 * dym2)
      yout3 = y3 + h6 * (dy3dx + dyt3 + 2.0 * dym3)

      return
      
      end
