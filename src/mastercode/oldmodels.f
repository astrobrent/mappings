cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******Unsupported models.
c
c
      subroutine oldmodels()
c
c
      include 'cblocks.inc'
c
c
c
c
      character ilgg*4
c
   50 write(*, 100) 
  100 format(///' Superceded models, (not supported)'/
     &          ':::::::::::::::::::::::::::::::::::::::::')
c

      write(*, 340) 
 340  format('    S1  :  Shock, fixed grid'/
     &       '    P1  :  Photoionisation, fixed'/
     &       '    P2  :  Photoionisation, diffuse field'/
     &       '        :  '/
     &       '    E,X :  Exit'/
     &       ' :: ',$)
      read(*, 350 ) ilgg
 350  format(a)
c
c     fiddle lower case entries..
c
      if (ilgg(1:1) .eq.'s') ilgg(1:1) ='S'
      if (ilgg(1:1) .eq.'p') ilgg(1:1) ='P'
      if (ilgg(1:1) .eq.'e') ilgg(1:1) ='E'
c
c     X for exit as well
c
      if (ilgg(1:1) .eq.'X') ilgg(1:1) ='E'
      if (ilgg(1:1) .eq.'x') ilgg(1:1) ='E'
c
c    ***ZERO BUFFER ARRAYS AND RESET COUNTERS
c
      call zer
c
c
      runname = 'Old Run'
c
      if (ilgg(1:2) .eq. 'P1') call photo1
      if (ilgg(1:2) .eq. 'P2') call photo2
      if (ilgg(1:2) .eq. 'S1') call shock1
      if (ilgg .eq. 'X') return
c
      return
c
      end
