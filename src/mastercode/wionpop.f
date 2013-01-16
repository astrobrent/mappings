cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    This will become the standarised method for writing out
c    just population data.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wionpop(luop,po)
c
c
      include 'cblocks.inc'
c
      integer*4 luop,ionmax,i,j
      double precision po(mxion, mxelem),chcksum
c
c    if luop > 0 then writing to a file
c
      write(luop, 2050) (elem(i),i = 1, atypes)
 2050 format(1h ,t8,16(4x,a2,4x)/)
c
c    find highest ionisation
c
      ionmax = 0
c
      do i = 1,atypes
          if (maxion(i).ge.ionmax) ionmax = maxion(i)
      enddo
c
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
c
      return
c
      end
