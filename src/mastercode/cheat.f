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
c     Charge exchange heating.  Uses the energy defect from
c     charge exchange reactions.
c
c     chgain is a heating term so -ve is cooling
c
c     v2.5
c
c     RSS 4/91
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine cheat(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision t, de, dh, delta, ab
      integer*4 i,ion,at,rat
c
      chgain = 0.d0
c
      if (chargemode.eq.0) then
c
c     new reactions being used
c
         do i = 1,nchxi
c
            if (chxi(i).gt.0.d0) then
               at = chxiat(i)
               rat = chxix(i)
               ion = chxiio(i)
               delta = chxicos(4,i)
               ab = dh*zion(at)*pop(ion,at)*dh*pop(2,rat)
               if (ab.ge.pzlimit) then
                  chgain = chgain+chxi(i)*ab*ev*delta
               endif
            endif
         enddo
      endif
c
      return
c
      end
