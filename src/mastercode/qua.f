cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c***************************************************
c     A quadratic regression routine.  At present it 
c     accecpts 3 points and so does a perfect fit
c     but it can be extended to fit more points by
c     increasing n and dimensioning larger arrays
c     RSS 1990.
c***************************************************
c     WARNING: The following must be true
c
c     x(1)<>x(2)<>x(3)
c
c     However ordering is not important
c
c***************************************************
      subroutine qua(x,y,a,b,c)
c
      include 'const.inc'
c
      double precision x(3),y(3),a,b,c
      double precision sxx,sxy,sxx2,sx2y,sx2x2
      double precision sigx,sigx2,sigx3,sigx4
      double precision sigy,sigxy,sigx2y,det
      integer*4 n,i
c
c     initialize, don't assume zeros in variables
c
      n = 3
      sigx = 0.0d0
      sigx2 = 0.0d0
      sigx3 = 0.0d0
      sigx4 = 0.0d0
      sigy = 0.0d0
      sigxy = 0.0d0
      sigx2y = 0.0d0
c
      do 100 , i = 1,n
         sigx = sigx+x(i)
         sigx2 = sigx2+(x(i)*x(i))
         sigx3 = sigx3+(x(i)*x(i)*x(i))
         sigx4 = sigx4+(x(i)*x(i)*x(i)*x(i))
         sigy = sigy+y(i)
         sigxy = sigxy+(x(i)*y(i))
         sigx2y = sigx2y+(x(i)*x(i)*y(i)) 
 100  continue
c
      sxx = sigx2 - (sigx*sigx)/n
      sxy = sigxy - (sigx*sigy)/n
      sxx2 = sigx3 - (sigx*sigx2)/n
      sx2y = sigx2y - (sigx2*sigy)/n
      sx2x2 = sigx4 - (sigx2*sigx2)/n
c
      det = (sxx*sx2x2)-(sxx2*sxx2)
c
      a = ((sx2y*sxx)-(sxy*sxx2))/det
      b = ((sxy*sx2x2)-(sx2y*sxx2))/det
      c = (sigy/n)-(b*(sigx/n))-(a*(sigx2/n))
c
      return
c
      end
