cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES RESONANCE LINE COOLING
c       NB.  SUBR. HYDRO SHOULD BE CALLED PREVIOUSLY
c
c	RETURNS RLOSS = RESONANCE COOLING RATE(ERG CM-3 S-1)
c	RBRI(I,J) = BRIGTHNESS OF EACH LINE (ERG CM-3 S-1 STER-1)
c	AVAILABLE IN COMMOM BLOCK /RLINE/
c	CALL SUBR. RATEC, USES FUNCTION FGAUNT
c
      subroutine reson(t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision t, de, dh
      double precision aa,abup,e,f
      double precision po,r12,rr,u,zi
c
      integer*4 jel,jjj,m
c
c
c           Functions
c
      double precision fgaunt
c
c	internal functions
c
      double precision a , x1, fr, lg,x
c
      a(x) = (((5.232d0*x)-6.652d0)*x)+2.878d0
      x1(x) = 0.7179d0*(x**(-0.0584d0))
      fr(x) = 0.588d0*(x**(-0.234d0))
      lg(x) = 0.07d0*(x**(-0.085d0))
c
      rloss = 0.0d0
c
      f = dsqrt(1.0d0/(t+1.d-2))
      u = dmin1(10.0d0,dmax1(0.1d0,t/1.d4))
c
c    ***GET R12 , HENCE COOLING RATE
c
      do m = 1, nlines
c      
         rr = 0.0d0
         jel = ielr(m)
         jjj = ionr(m)
         rbri(m) = 0.d0
         zi = zion(jel)
         po = pop(jjj,jel)
         if (zi*po.ge.pzlimit) then
         e = e12r(m)
         abup = de*dh*zi*po
c
         aa = e/(rkb*t)
         r12 = ((rka*f)*omres(m))*dexp(- aa)
         rr = ((r12*abup)*e)*fgaunt(jjj,0,aa)
         rr = rr/(0.413497d0*a(1.0d0/dble(jjj)))
c
         rloss = rloss+rr
c     
         rbri(m) = rr/(4.d0*pi)
c
         endif
c
      enddo
c     
      return 
c
c
      end
      
