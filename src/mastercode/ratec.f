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
c*******TO COMPUTE COLLISIONAL EXCITATION RATES
c	AT TEMP.: TE  FOR THE N-1 TRANSITION OF
c	HYDROGEN OR HELIUM IONS
c
c	REF.: JOHNSON (1972) APJ 174, 227.
c	      MEWE (1972) A&A 20, 215.
c
c	IONIC CHARGE : NZ     LEVEL : N
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ratec(nz, n, te, rate, ey)
c
      include 'const.inc'
c
      double precision te,rate,ey,ra
      double precision a,an,b,bn,bra,c,coeff,const
      double precision d,ey22,ez21,fij,fn,fra,gag,pbra
      double precision t,u,y,z
      integer*4 nz,n
c
      double precision fue1
c
      ra = 5.2917d-9
      t = te+1.d-5
      u = 1.0d0-(1.0d0/(n**2))
      ey = (ryde*u)*(nz**2)
      y = ey/(rkb*t)
c
c   ****HELIUM ION CASE
      z = (u*0.45d0)+y
      if (nz .eq. 2) then
c
c     Helium Case......
c
         if (n .eq. 2) then
            a = 0.08d0
            b = -0.16d0
            c = 0.11d0
            d = 0.0d0
            a = a+0.04d0
            b = b+0.04
            c = c+0.02d0
            d = d+0.28d0
            fij = 0.8324d0
         else if (n .eq. 3) then
            a = 0.2d0
            b = 0.06d0
            c = 0.0d0
            d = 0.28d0
            fij = 0.1582d0
         else if (n .eq. 4) then
            a = 0.25d0
            b = 0.04d0
            c = 0.0d0
            d = 0.28d0
            fij = 0.0580d0
         else if (n .eq. 5) then
            a = 0.27d0
            b = 0.03d0
            c = 0.0d0
            d = 0.28d0
            fij = 0.0279d0
         else if (n .eq. 6) then
            a = 0.28d0
            b = 0.02d0
            c = 0.0d0
            d = 0.28d0
            fij = 0.0156d0
         else
            a = 0.29d0
            b = 0.02d0
            c = 0.0d0
            d = 0.28d0
            fij = 3.2d0/(n ** 3)
         end if
c     
         gag = (a+((((b*y)-((c*y)*y))+d)*fue1(y,0.d0)))+(c*y)
         rate = (((2.724d-15*dexp(-y))*gag)*fij)/(dsqrt(t)*ey)
c
c
c
      else
c
c   ****HYDROGEN CASE
c
         bn = 4.0d0*((1.0d0+(4.d0/(3.d0*u)))-(0.603d0/(u*u)))
         bn = bn/((n**3)*(u**2))
         fn = (1.133d0-(0.4059d0/u))+(0.07014d0/(u*u))
         fn = (1.96028d0*fn)/((n*u)**3)
         an = (2.0d0*fn)/u
         const = (4.0d0*(ra**2))*dsqrt((((2.0d0*rkb)*t)*pi)/me)
         coeff = (const*(y**2))/u
         pbra = (((1.0d0/z)+0.5d0)*fue1(z,z))*an
         c = an*((1.0d0/y)+0.5d0)
         ez21 = dexp(- z)-(z*fue1(z,z))
         bra = (c*fue1(y,y))-pbra
         ey22 = dexp(- y)-(y*fue1(y,y))
         bra = bra+((bn-(an*dlog(2.0d0/u)))*((ey22/y)-(ez21/z)))
         rate = coeff*bra
c
c     ***COMPENSATE FOR RESONANCE AT LEVEL*3 OF HYDROGEN
c
         if (n .eq. 3) then
            fra = 1.8d0*rate+5.47d-11/dsqrt(t)*
     &           (0.34d0*(t+15.8d4)*dexp(-15.8d4/t)
     &           -0.40d0*(t+25.d4)*dexp(-25.d4/t)
     &           +0.06d0*(t+84.d4)*dexp(-84.d4/t))
            if ((t.le.2.d4).and.(t.ge.5.d3)) fra = (3.5d0*rate+fra)/2.d0
            rate = fra/3.5d0
         end if
c
c
c
      end if
c
c     
      return
c 
      end
