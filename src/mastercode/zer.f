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
c*******TO ZERO BUFFER ARRAYS,RESET COUNTERS AND
c	SET DEFAULT VALUES
c
c
c
c
      subroutine zer()
c
      include 'cblocks.inc'
c
      integer*4 ion, j, i,k, series, line,atom
c
c    ***ZERO SUMMED LINE DATA
c
      ion = atypes
      do i = 1, nfions
         do  j = 1, nftrans
            fluxf(j,i) = 0.0d0
         enddo
      enddo
c
      do i = 1, 10
         fluxh(i) = 0.0d0
      enddo
      do i = 1, nheilines
         fluxhei(i) = 0.0d0
      enddo
c      

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
c
      h2qav = 0.0d0
      hei2qa = 0.0d0
      heii2qa = 0.0d0
      do i = 1,atypes
         h2q(i) = 0.d0
      enddo
c
c    ***MAKE RECOMB. COEFFICIENTS EFFECTIVE FOR ALL ELEMENTS
c
   50 fhbeta = 0.0d0
      do  i = 1, atypes
      do  j = 1, maxion(i)
         arad(j,i) = dabs(arad(j,i))
      enddo
      enddo
c
c    ***SET DEFAULT POP. AND ZERO THE WEIGHTED IONIC QUANTITIES
c
      tem = -1
      do  i = 1, atypes
      do  j = 1, maxion(i)
      pop(j,i) = 0.0d0
      teav(j,i) = 0.0d0
      popint(j,i) = 0.0d0
      deam(j,I) = 0.0d0
      rdisa(j,i) = 0.0d0
      pam(j,i) = 0.0d0
      enddo
      pop(1,i) = 1.0d0
      enddo

      deav = 0.0d0
      dhav = 0.0d0
      dhavv = 0.0d0
      vunilog = 0.0d0
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
c    ***ZERO PHOTON FIELD VECTORS
c
      do 200 i = 1, infph
         skipbin(i) = .true.
         soupho(i) = 0.0d0
         dwdif(i) = 0.0d0
         updif(i) = 0.0d0
         emidif(i) = 0.0d0
         tphot(i) = 0.0d0
         IRphot(i) = 0.d0
c         cspec(i) = 0.0d0
 200  continue
c
c     zero grain heating
c
      gheat = 0.d0
      gcool = 0.d0
      paheat = 0.d0
c
c     Dust Column
c
      if(grainmode) then
        dustint = 0.d0
        if(pahfrac.gt.0.d0) then
          do i=1,pahi
            pahint(i) =0.d0
          enddo
   	  pahmode=.FALSE.
          if(ClinPAH) then
	        zion(zmap(6)) = zion0(zmap(6))          
          else
  	        zion(zmap(6))=zion0(zmap(6))* 
     &	      		  (dion(zmap(6))+(1.d0-dion(zmap(6)))*pahCfrac)
          
          endif
        endif    
      endif
c
      do 207 i = 1, xlines
         emilin(1,i) = 0.0d0
         emilin(2,i) = 1.0d0
         uplin(1,i) = 0.0d0
         uplin(2,i) = 1.0d0
         dwlin(1,i) = 0.0d0
         dwlin(2,i) = 1.0d0
  207 continue
c
      qhi = 0.0d0
      qhei = 0.0d0
      qheii = 0.0d0
      zstar = 0.0d0
      alnth = 0.0d0
      teff = 0.0d0
      crate = 0.d0
      cosphi = 0.d0
      do k=1,numtypes
       avgpot(k) = 0.d0
      enddo
c
c     flow globals
c
      vel0 = 0.d0
      vel1 = 0.d0
      rho0 = 0.d0
      rho1 = 0.d0
      pr0 = 0.d0
      pr1 = 0.d0
      dh0 = 0.d0
      dh1 = 0.d0
      de0 = 0.d0
      de1 = 0.d0
      te0 = 0.d0
      te1 = 0.d0
c
c
c
      cut = 0.0d0
      iphom = -1
      ipho = 0
c
c    ***SET DEFAULT VALUES FOR "ON THE SPOT APPROX" : JSPOT
c    ***FOR CASE A,B (0=CASE A  1=CASE B)  :  CASEAB
c    ***NUMBER OF ABSORBING ATOMIC ELEMENTS : LIMPH
c
      jspot = 'NO'
      jcon = 'YES'
      jspec = 'NO'
      limph = atypes
      do atom = 1,atypes
         caseab(atom) = 0.0d0
      enddo
c
      return 
      end

