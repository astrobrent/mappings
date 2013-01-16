cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c     
c     copyright 1994 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c
c       Version 1.0.0r
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure to return the derivatives of r,v,P
c	  at any point in a spherical flow 
c     wrt dummy variable u
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sphderivs(u, r, v, P, drdu, dvdu, dPdu)
c
c     
c     Global consts and vars
c     
      include 'cblocks.inc'
c
c     parameters
c     
      real*8 u, r, v, P, drdu, dvdu, dPdu
c     
c     local variables
c     
      real*8 a2, rho, term1, v2, Phi2,h,c
c
c	simple spherical potential
c
      Phi2 = GM / (r * r)
c     
c     read in from file
c
      rho = 1.6e-24 
c     
      a2 = gamma * P / rho
      v2 = v * v
c     
c     
      C = 0.d0
      H = 0.d0
c
c      write(*,*) ' derivs: ',a2,phi2,r
c     
      term1 = ((a2 + a2) / r) - phi2 + ((gamma - 1.0) * (H-C)) / v
c     
      drdu = (v2 - a2)
      dvdu = (v * term1)
      dPdu = -rho * (Phi2 * drdu + (v * dvdu))
c
c      write(*,*)'Derivs: drdu dvdu dPdu',drdu,dvdu,dPdu
c
      return
c
      end

