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
c     Special Purpose Newtons Method for solving a
c
c     well behaved quartic.
c
c     modified RSS 10.2009
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine quart(a,x1,x2,root,rmod)
c
      include 'const.inc'
c
      double precision a(5),x1,x2,root
      double precision d(5),x,f,df,dx,eps
      integer*4 i,itmax
      character rmod*4
c
c     find -ve curve near root wanted
c
      x  = 0.0
      if (rmod.eq.'FAR') then 
      	dx = (x2)/20.d0
      	do i = 1,40
          x = x+dx
          f = a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
          if (f.gt.0) goto 50
      	enddo
      else
      	dx = (x1)/20.d0
      	do i = 1,40
          x = x+dx
          f = a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
          if (f.lt.0) goto 50
      	enddo
      endif
c     
 50   x = x - (dx*0.5)
c     
c     get derivative
c     
      d(5) = 0.d0
      do i = 1,4
         d(i) = a(i+1)*dble(i)
      enddo
c     
      eps = x*1.d-8
      itmax = 20
c     
      do i = 1,itmax
         f = a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
         df = d(1)+x*(d(2)+x*(d(3)+x*d(4)))
         dx = f/df
         x = x-dx
         if (dabs(dx).lt.(eps)) goto 100
      enddo
c     
 100  root = x
c
      return
c
      end
