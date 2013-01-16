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
c*******TO OBTAIN THE AVERAGE TEMP. , IONIC POP.,DENSITIES AND
c       LINE INTENSITIES USING SUM DATA IN  : /SLINE/,/TONLIN/
c       AND /DELIN/
c       CALL SUBR. FINDTDE
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine avrdata()
c
      include 'cblocks.inc'
c
c
c
      double precision rno
      integer*4 i, j, line,series
c
c    ***AVERAGES IONIC POP,TEMP.,DISTANCES AND DENSITIES
c
      rno = 1.d17
      do 2041 j = 1, atypes
      do 2040 i = 1, maxion(j)
         teav(i,j) = teav(i,j)/(epsilon+pam(i,j))
         deam(i,j) = deam(i,j)/(epsilon+pam(i,j))
         rdisa(i,j) = (rdisa(i,j)/(epsilon+pam(i,j)))*rno
         pam(i,j) = pam(i,j)/(epsilon+deav)
 2040 continue
 2041 continue
c
      deav = dhav/deav
      qhdha = qhdha/(epsilon+dhavv)
c
      zetaeav = zetaeav/(epsilon+dhavv)
      roiii = roiii/(epsilon+weoiii)
      deoiii = deoiii/(epsilon+weoiii)
      rnii = rnii/(epsilon+wenii)
      denii = denii/(epsilon+wenii)
      roii = roii/(epsilon+weoii)
      toii = toii/(epsilon+weoii)
      rsii = rsii/(epsilon+wesii)
      tsii = tsii/(epsilon+wesii)
c
c    ***FIND TNII,TOIII,DESII,DEOII
c
      call findtde
c
c    ***FORM RATIO OF EMISSION LINES RELATIVE TO  H-BETA
c
         do  i = 1, nfions
            do  j = 1, nftrans
              fluxf(j,i) = fluxf(j,i)/(epsilon+fhbeta)
            enddo
         enddo
c
         do  j = 1, nfetrans
            fluxfe(j) = fluxfe(j)/(epsilon+fhbeta)
         enddo
c
         do i = 1,nf3ions
            do  j = 1, nf3trans
               fluxf3(j,i) = fluxf3(j,i)/(epsilon+fhbeta)
            enddo
         enddo
c
         do i = 1,n6ions
            do  j = 1, n6trans
               fluxf6(j,i) = fluxf6(j,i)/(epsilon+fhbeta)
            enddo
         enddo
c
         do i = 1,n9ions
            do  j = 1, n9trans
               fluxf9(j,i) = fluxf9(j,i)/(epsilon+fhbeta)
            enddo
         enddo
c     
      do j = 1, 10
         fluxh(j) = fluxh(j)/(epsilon+fhbeta)
      enddo
      do j = 1, nheilines
         fluxhei(j) = fluxhei(j)/(epsilon+fhbeta)
      enddo
c
        heiioiiibfsum = heiioiiibfsum/(epsilon+fhbeta)
c
      do series = 1, 6
          do line = 1,10 
            hydroflux(line,series) = 
     &          hydroflux(line,series)/(epsilon+fhbeta)
            heliflux(line,series) = 
     &          heliflux(line,series)/(epsilon+fhbeta)
          enddo
      enddo
c
      do  i = 1, atypes
      do series = 1, 2
      do line = 1,10 
            xhydroflux(line,series,i) = 
     &          xhydroflux(line,series,i)/(epsilon+fhbeta)
          enddo
      enddo
      enddo
c
      do j = 1, nlines
         fluxr(j) = fluxr(j)/(epsilon+fhbeta)
      enddo
c
      do 2396 j = 1, mlines
 2396    fluxi(j) = fluxi(j)/(epsilon+fhbeta)
c
      do 2398 j = 1, xlines
 2398    fluxx(j) = fluxx(j)/(epsilon+fhbeta)
c
      do 2399 j = 1, xilines
 2399    fluxxi(j) = fluxxi(j)/(epsilon+fhbeta)
c
      do 2400 j = 1, xhelines
 2400    fheif(j) = fheif(j)/(epsilon+fhbeta)
c
      fheiilr = fheiilr/(epsilon+fhbeta)
      h2qav   = h2qav/(epsilon+fhbeta)
      hei2qa  = hei2qa/(epsilon+fhbeta)
      heii2qa = heii2qa/(epsilon+fhbeta)
c
      do i= 1, atypes
         h2q(i) = h2q(i)/(epsilon+fhbeta)
      enddo
c
      return 
      end
