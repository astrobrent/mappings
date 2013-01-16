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
c     Cosmic Ray Heating, based on simple event rate
c     and thermalisation of mean ejected electrons.
c     
c     Uses Shull 1979 H & He results. Assuming Eo = 35 eV.
c
c     heating = rate * dh(neutral) * Eh(Eo)
c     implicitly includes He contribution
c
c     #secondaries = phi(Eo)
c     adds to primary rate for each of HI, HeI and HeII
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cosmic(t,de,dh)
c
c
      include 'cblocks.inc'
c
c
      double precision t,de,dh,eh
      double precision xe,lxe,a(0:3),ph(0:3)
c
      data a/3.4992857d+1,-4.8392857d+0,-7.4285714e+0,-1.1250000d+0/
      data ph/1.5714286d-3,1.0297619d-1,2.0928571d-1,3.4166667d-2/
c
      cosgain = 0.d0
c
      if (crate.gt.0.d0) then
c
      xe = pop(2,1)/(pop(1,1)+epsilon)
      lxe = dlog10(xe+epsilon)
c
c     cubic fit for eh to log (xe) at Eo = 35 ev
c
      eh = a(0)+lxe*(a(1)+lxe*(a(2)+lxe*a(3)))
c
c     sanity
c
      if (eh.gt.35.d0) eh = 35.d0
      if (lxe.lt.-4.d0) eh = 7.d0
      if (eh.lt.7.d0) eh = 7.d0
c
      cosgain = dh*pop(1,1)*eh*crate*ev
c
c     cubic fit for phi to log (xe) at Eo = 35 ev
c
      cosphi = ph(0)+lxe*(ph(1)+lxe*(ph(2)+lxe*ph(3)))
c
c     sanity
c
      if (cosphi.lt.0.d0) cosphi = 0.d0
      if (lxe.lt.-4.d0) cosphi = 0.75
c
      endif
c
      return
c
      end

