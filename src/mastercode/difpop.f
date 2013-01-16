cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******TO DERIVE RELATIVE CHANGE IN IONIC POPULATIONS : DIF
c	LIM : NUMBER OF LAST ELEMENTS TO BE INCLUDED
c	TRE : TRESHOLD FOR MEASURING CHANGE
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine difpop(popin, popfi, tre, lim, dif)
c
      include 'cblocks.inc'
c
c
      double precision popin(mxion, mxelem), popfi(mxion, mxelem)
      double precision tre, dif, thresh, maxdif, tmdif
      integer*4 num, iel, ion, lim, matom, mion
c
      dif = 0.0d0
      num = 0
      thresh = tre
c
      matom = 0
      mion = 0
      maxdif = 0.d0
c
      if (thresh.le.epsilon) thresh = epsilon
      do iel = 1, lim
         do ion = 1, maxion(iel)
            if ((popin(ion,iel).ge.thresh).and.
     &           (popfi(ion,iel).ge.thresh)) then
               tmdif =  dabs(dlog10(popin(ion,iel))
     &              -dlog10(popfi(ion,iel)))
               if (tmdif.ge.maxdif) then
                  maxdif = tmdif
                  matom = iel
                  mion = ion
               endif
               dif = dif + tmdif
               num = num+1
            endif
         enddo
      enddo
c
      dif = dif/(num+epsilon)
c
      if (expertmode.gt.0) then
      write(*,*) 'Max difference: ',elem(matom),rom(mion),'=', maxdif
      endif
c
      return 
c
      end
