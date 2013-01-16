cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    This will become the standarised method for writing out
c    element abundances and population data.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wionabal(luop,po)
c
c
      include 'cblocks.inc'
c
      integer*4 luop,ionmax,i,j
      double precision po(mxion, mxelem)
      double precision chcksum
c
c    if luop > 0 then writing to a file
c
c      if (luop.gt.0) then
      write(luop, 2050) (elem(i),i = 1, atypes)
 2050 format(1h ,t8,16(4x,a2,4x)/)
      if (.NOT.grainmode) then
      write(luop, 2060) (zion(i),i = 1, atypes)
 2060 format(' Abund. :',16(1pg10.3)//)
      endif
      if (grainmode) then
      write(luop, 2065) (zion0(i),i = 1, atypes)
      write(luop, 2066) (dion(i),i = 1, atypes)
      write(luop, 2067) (zion(i),i = 1, atypes)
 2065 format(' Tot Ab.:',16(1pg10.3)//)
 2066 format(' Depl.  :',16(1pg10.3)//)
 2067 format(' Gas Ab.:',16(1pg10.3)//)
      endif
      write(luop,*) ' '
c
c    find highest ionisation
c
      ionmax = 0
      do i = 1,atypes
          if (maxion(i).ge.ionmax) ionmax = maxion(i)
      enddo
      
      if (ionmax.gt.29) ionmax = 29
c
c
c
 2070 format(a6,2x,16(1pg10.3))
      do j = 1,ionmax
         chcksum = 0.0d0
         do i = 1,atypes
            chcksum = chcksum + po(j,i)
         enddo
         if (chcksum.ne.0.0d0) then
         write(luop, 2070) rom(j), (po(j,i),i = 1, atypes)
         endif
      enddo
c
c
c      else
c
c    writing to standard output, only write abundances
c
c      write(*, 300)
c      write(*, 310) (elem(j), zion(j),j = 1, atypes)
c  300 format(//'Elemental abundances are:'/)
c  310 format(4(4x,a2,1pe10.3))
c
c
c      endif
c
c
c
      return
c
      end
