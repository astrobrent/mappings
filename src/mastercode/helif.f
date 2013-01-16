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
c
c     subroutine to compute collisionally coupled
c     intercombination-forbiden transitions in He like
c     ions (6<z<26).
c
c     RSS 11/90
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine helif(t,de,dh)
c
c
c
      include 'cblocks.inc'
c
c
c
c
      double precision t, de, dh,abup
      double precision betc,rrf,rri,rf,ri
      double precision gbarf,gbari,yi,yf
      double precision ff,fi,bra,z,theta
      double precision brf,cf,cgjbf,cgjbi,ci,egf,egi,omegf
      double precision omegi,pirt38,rat
c
      integer*4 i,j,code,atom,ion
c
c           Functions
c
      double precision fresga
c
      xheloss = 0.d0
c
      pirt38 = pi*8.d0/dsqrt(3.d0)
      t = dmax1(0.1d0,t)
c
      do i = 1, xhelines/2
c
c       get data from arrays for lines # 2i-1 and 2i
c
c     Note: only line for possible species are read in
c     so no need to check maxion etc...
c
         j = 2*i -1
c
         atom = xheat(j)
         ion = xheion(j)
c
c
         xhebri(j) = 0.d0
         xhebri(j+1) = 0.d0
c
         abup = zion(atom)*pop(ion,atom)
c
c     only calculate abundant species
c
         if (abup.ge.pzlimit) then
c
c
c
            egi = xhejk(j)
            egf = xhejk(j+1)
            fi = xhef(j)
            ff = xhef(j+1)
c
c     get mean gaunt factor
c
            gbari = fresga(atom,2,5.d0,egi,t,code)
            gbarf = fresga(atom,2,6.d0,egf,t,code)
c
c     get scaled energy gap to level j from ground (not k)
c
            yi = egi*ev/(rkb*t)
            yf = egf*ev/(rkb*t)
c
c     transition rate rr
c
            rri = 0.d0
            rrf = 0.d0
c
            if ((yi.lt.huge).and.(yf.lt.huge)) then
c
c     casacade sums
c
            ci = 1.065d0
            rat = 7.d0*(yf+0.6d0)/(yf+1.15)
            cf = 1.0d0+(0.4d0*rat*dexp(-0.21d0*yf))
c
c     assuming br=1 and a = 1
c
            omegi = pirt38*fi*ci*gbari*ryd/egi
            omegf = pirt38*ff*cf*gbarf*ryd/egf
c
            cgjbi = rka*(dexp(-yi)/dsqrt(t))
            cgjbf = rka*(dexp(-yf)/dsqrt(t))
c
            rri = cgjbi*omegi
            rrf = cgjbf*omegf
c
c     now couple rates using inter branching ratio betc
c
            bra = xhebra(atom)
            z = mapz(atom)
            theta = t / (z*z*1.d5)
            betc = 330.d0*(z**(-14.57d0))*(theta)**(-0.73d0*
     &           (z**(-0.415)))
            brf = 1/(1+(betc*de*bra))
c
            rf = (rrf+(1-bra)*rri)*brf
            ri = (rri*bra)+(rrf+(1-bra)*rri)*(1-brf)
c
         endif
c
c     number to transition is in abup
c
c
c     photons cm^-3 s^-1
c
            rri = de*dh*abup*ri
            rrf = de*dh*abup*rf
c
c     ergs...
c
            rri = ev*egi*rri
            rrf = ev*egf*rrf
c
            xheloss = xheloss+rri
            xheloss = xheloss+rrf
c
            xhebri(j) = rri/(4.0d0*pi)
            xhebri(j+1) = rrf/(4.0d0*pi)
c
c     end pop limited loop
c
        endif
c
      enddo
c
      return
c
      end
