cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c*******TO ZERO BUFFER VECTORS ONLY
c
c
c
      subroutine zerbuf()
c
      include 'cblocks.inc'
c
c
c
      integer*4 ion, j, i, series, line
c      integer*4 atom
c
c    ***ZERO SUMMED LINE DATA
c
      ion = atypes
      do 10 i = 1, nfions
         do 10 j = 1, nftrans
 10         fluxf(j,i) = 0.0d0
 
      do i = 1, 10
         fluxh(i) = 0.0d0
      enddo
      
      do i = 1, nheilines
         fluxhei(i) = 0.0d0
      enddo
      
	  heiioiiibfsum = 0.d0
	  heiioiiibf    = 0.d0
      do line = 1,10 
      	do series = 1, 6
         hydroflux(line,series) = 0.0d0
         heliflux(line,series) = 0.0d0
         
         hydlin(1,line,series) = 0.0d0
         hydlin(2,line,series) = 1.0d0
         hyduplin(1,line,series) = 0.0d0
         hyduplin(2,line,series) = 1.0d0
         hyddwlin(1,line,series) = 0.0d0
         hyddwlin(2,line,series) = 1.0d0

         hellin(1,line,series) = 0.0d0
         hellin(2,line,series) = 1.0d0
         heluplin(1,line,series) = 0.0d0
         heluplin(2,line,series) = 1.0d0
         heldwlin(1,line,series) = 0.0d0
         heldwlin(2,line,series) = 1.0d0

        enddo
      enddo
c
      do  i = 1, atypes
      do line = 1,10 
          do series = 1, 2
          
            xhydroflux(line,series,i) = 0.d0
            xhydrobri(line,series,i)  = 0.d0

            xhydlin(1,line,series,i)  = 0.d0
            xhydlin(2,line,series,i)  = 1.d0
            xhyduplin(1,line,series,i) = 0.d0
            xhyduplin(2,line,series,i) = 1.d0
            xhyddwlin(1,line,series,i) = 0.d0
            xhyddwlin(2,line,series,i) = 1.d0
            
          enddo
      enddo
      enddo
c     
      do i = 1, mlines
         fluxi(i) = 0.0d0
      enddo
c
      do i = 1, nlines
         fluxr(i) = 0.0d0
      enddo
c
      do i = 1, xlines
          fluxx(i) = 0.0d0
      enddo
c
      do i = 1, xilines
         fluxxi(i) = 0.0d0
      enddo
c
      do i = 1, xhelines
         fheif(i) = 0.0d0
      enddo
c      
      do i = 1, nfetrans
         fluxfe(i) = 0.0d0
      enddo
c
      do i = 1,nf3ions
         do j = 1, nf3trans
            fluxf3(j,i) = 0.d0
         enddo
      enddo
c
      do i = 1,n6ions
         do j = 1, n6trans
            fluxf6(j,i) = 0.d0
         enddo
      enddo
c
      do i = 1,n9ions
         do j = 1, n9trans
            fluxf9(j,i) = 0.d0
         enddo
      enddo
c
      h2qav = 0.0d0
      hei2qa = 0.0d0
      heii2qa = 0.0d0
c
      do i = 1, atypes
         h2q(i) = 0.0d0
      enddo
c
 50   fhbeta = 0.0d0
      deav = 0.0d0
      dhav = 0.0d0
      dhavv = 0.0d0
      qtosoh = 0.0d0
      qhdha = 0.0d0
      zetaeav = 0.0d0
      tlosac = 0.0d0
      tforbi = 0.0d0
      fheiilr = 0.0d0
      vas111 = 0.0d0
      vas222 = 0.0d0
      vas333 = 0.0d0
      vas444 = 0.0d0
      vas555 = 0.0d0
c
c    ***ZERO IONIC WEIGHTED QUANTITIES
c
      do i = 1, atypes
         do j = 1, maxion(i)
            teav(j,i) = 0.0d0
            popint(j,i) = 0.0d0
            deam(j,i) = 0.0d0
            rdisa(j,i) = 0.0d0
            pam(j,i) = 0.0d0
         enddo
c         popint(maxion(i)+1,i) = 0.0d0
      enddo
c
      weoiii = 0.0d0
      wenii = 0.0d0
      wesii = 0.0d0
      weoii = 0.0d0
      roiii = 0.0d0
      rnii = 0.0d0
      rsii = 0.0d0
      roii = 0.0d0
      deoiii = 0.0d0
      denii = 0.0d0
      desii = 0.0d0
      deoii = 0.0d0
      toiii = 0.0d0
      tnii = 0.0d0
      tsii = 0.0d0
      toii = 0.0d0
c
c
c     get ion entry numbers for forbidden OI,OIII,NII,SII,OII
c
         ox3 = 0
         ni2 = 0
         su2 = 0
         ox2 = 0
         ox1 = 0
c
c         
         do i = 1, n6ions
            if ((f6ion(i).eq.3).and.(f6atom(i).eq.zmap(8))) ox3 = i
            if ((f6ion(i).eq.2).and.(f6atom(i).eq.zmap(7))) ni2 = i
         enddo
c
         do i = 1,nfions
            if ((fion(i).eq.1).and.(fatom(i).eq.zmap(8))) ox1 = i
            if ((fion(i).eq.2).and.(fatom(i).eq.zmap(16))) su2 = i
            if ((fion(i).eq.2).and.(fatom(i).eq.zmap(8))) ox2 = i
         enddo
c
c    ***ZERO PHOTON FIELD ARRAYS : UPDIF,DWDIF,UPLIN,DWLIN....
c	BUT NOT THE SOURCE VECTOR : SOUPHO
c     
      do 200 i = 1, infph
         dwdif(i) = 0.0d0
         updif(i) = 0.0d0
         emidif(i) = 0.0d0
         tphot(i) = 0.0d0
c         cspec(i) = 0.0d0
 200  continue
c     
      do 207 i = 1, nlines
         emilin(1,i) = 0.0d0
         emilin(2,i) = 1.0d0
         uplin(1,i) = 0.0d0
         uplin(2,i) = 1.0d0
         dwlin(1,i) = 0.0d0
         dwlin(2,i) = 1.0d0
 207  continue
c     
      return 
c     
      end
      
