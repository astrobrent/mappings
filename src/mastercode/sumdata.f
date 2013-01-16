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
c
c
c*******TO INTEGRATE THE LINE AND TEMPERATURE DATA AND THE
c       COLUMN DENSITY
c       USES TEMPERATURE (T), EL. DENS. (DE), HYDR. DENS (DH),
c       FILLING FACT. (FI), THE VOLUME OF THE SLAB CONSIDERED (DVOLUNI)
c       (IN ARBITRARY UNITS (VUNIT)) AND FINALLY THE RADIUS STEP
c       SIZE (DR) (IN CM)
c       THE EQUIVALENT IN CC OF THE VOLUME UNIT : VUNIT  ARE IN
c       VUNILOG (LOG10)
c
c       IMOD = 'ALL'  : ADS ALL QUANTITIES IN BUFFER ARRAYS
c       IMOD = 'COLD' : DERIVES ONLY COLUMN DENSITY AND IMPACT PAR.
c       IMOD = 'REST' : ADS ALL THE OTHER QUANTITIES
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sumdata(t, de, dh, fi, dvoluni, dr, rdis, imod)
c
      
      include 'cblocks.inc'
c
c           Variables
c
      double precision rno,wei,wei2
      double precision t, de, dh, fi, dvoluni, dr, rdis
c      double precision fl,ft
      integer*4 series, line,i,j,atom
      integer*4 luty
c
      character imod*4
c
      rno = 1.d17
      luty = 6
c     
      if ((imod.ne.'ALL').and.(imod.ne.'COLD')
     &     .and.(imod.ne.'REST')) then
         write(luty, 212) imod
 212     format(/' MODE IMPROPERLY SET FOR SUBR. SUMDATA  : ',a4)
         stop 
      end if
c     
c     ***INTEGRATION FOR THE AVERAGE IONISATION : PAM(6,11) ,
c     WEIGHT IN : DEAV  ;  VOLUME STEP : DVOLUNI
c     AVERAGE ELECTRONIC DENSITIES IN DEAM(6,11)
c     AVERAGE DISTANCES IN RDISA(6,11)
c     
      if (imod .ne. 'COLD') then
         wei = ((de*dh)*fi)*dvoluni
         wei2 = wei*de
         dhav = dhav+wei2
         deav = deav+wei
         do 421 j = 1, atypes
            do 420 i = 1, maxion(j)
               deam(i,j) = deam(i,j)+(wei2*pop(i,j))
               rdisa(i,j) = rdisa(i,j)+((wei*pop(i,j))*(rdis/rno))
               pam(i,j) = pam(i,j)+(wei*pop(i,j))
 420        continue
 421     continue
c     
c     ***INTEGRATION FOR THE AVERAGE TEMPERATURES : TEAV(6,11) ,
c     WEIGHT IN : PAM(6,11)  ;  VOLUME STEP : DVOLUNI
c     
         do 426 j = 1, atypes
            do 425 i = 1, maxion(j)
               teav(i,j) = teav(i,j)+((t*wei)*pop(i,j))
 425        continue
 426     continue
c     
c     ***INTEGRATION FOR THE AVERAGE TEMPERATURES GIVEN BY LINE
c     RATIOS : TOIII,TNII  (INTEG. EL. DENS. IN : DEOIII,DENII)
c     WEIGHT IN : WEOIII,WENII  ;  VOLUME STEP : DVOLUNI
c
      if (ox3.ne.0) then
      wei = (f6bri(7,ox3)+f6bri(10,ox3))
      weoiii = weoiii+((wei*dvoluni)*fi)
      roiii = roiii+(((wei*dvoluni)*fi)*(wei/(epsilon+f6bri(13,ox3))))
      deoiii = deoiii+(((de*wei)*dvoluni)*fi)
      endif
      if (ni2.ne.0) then
      wei = (f6bri(7,ni2)+f6bri(10,ni2))
      wenii = wenii+((wei*dvoluni)*fi)
      rnii = rnii+(((wei*dvoluni)*fi)*(wei/(epsilon+f6bri(13,ni2))))
      denii = denii+(((de*wei)*dvoluni)*fi)
      endif
c     
c     ***INTEGRATION FOR THE AVERAGE DENSITIES GIVEN BY LINE
c     RATIOS : TOII,TSII  (INTEG. TEMP. IN : TOII,TSII)
c     WEIGHT IN : WEOII,WESII  ;  VOLUME STEP : DVOLUNI
c     
      if (ox2.ne.0) then 
      wei = fbri(1,ox2)
      weoii = weoii+((wei*dvoluni)*fi)
      roii = roii+(((wei*dvoluni)*fi)*(wei/(epsilon+fbri(2,ox2))))
      toii = toii+(((t*wei)*dvoluni)*fi)
      endif
      if (su2.ne.0) then 
      wei = fbri(2,su2)
      wesii = wesii+((wei*dvoluni)*fi)
      rsii = rsii+(((wei*dvoluni)*fi)*(wei/(epsilon+fbri(1,su2))))
      tsii = tsii+(((t*wei)*dvoluni)*fi)
      endif
c     
c     ***INTEGRATION OF LINE INTENSITIES IN /SLINE/
c     WEIGHT IN : FHBETA  ;  VOLUME STEP : DVOLUNI
c     
      wei = dvoluni*fi
c
         do  j = 1, xlines
            fluxx(j) = fluxx(j)+(xrbri(j)*wei)
         enddo
c
         do  j = 1, xilines
            fluxxi(j) = fluxxi(j)+(xibri(j)*wei)
         enddo
c
         do  j = 1, xhelines
            fheif(j) = fheif(j)+(xhebri(j)*wei)
         enddo
c
         do 200 j = 1, nlines
            fluxr(j) = fluxr(j)+(rbri(j)*wei)
 200     continue
c
         do 210 j = 1, mlines
            fluxi(j) = fluxi(j)+(fsbri(j)*wei)
 210     continue
c
         do j = 1, nfions
            do  i = 1, nftrans
               fluxf(i,j) = fluxf(i,j)+(fbri(i,j)*wei)
            enddo
         enddo
c
        do j = 1, nfetrans
           fluxfe(j) = fluxfe(j)+(febri(j)*wei)
        enddo
c
        do j = 1,nf3ions
           do i = 1, nf3trans
              fluxf3(i,j) = fluxf3(i,j)+(f3bri(i,j)*wei)
           enddo
        enddo
c     
        do j = 1,n6ions
            do i = 1, n6trans
               fluxf6(i,j) = fluxf6(i,j)+(f6bri(i,j)*wei)
            enddo
         enddo
c     
        do j = 1,n9ions
            do i = 1, n9trans
               fluxf9(i,j) = fluxf9(i,j)+(f9bri(i,j)*wei)
            enddo
         enddo
c     
         do  j = 1, 10
            fluxh(j) = fluxh(j)+(hbri(j)*wei)
         enddo
c
         do  j = 1, nheilines
            fluxhei(j) = fluxhei(j)+(heibri(j)*wei)
         enddo
c      
         heiioiiibfsum = heiioiiibfsum+(heiioiiibf*helibri(1,1)*wei)
c
         do line = 1,10 
            do series = 1, 6
              hydroflux(line,series) = hydroflux(line,series)
     &             +(hydrobri(line,series)*wei)
              heliflux(line,series) = heliflux(line,series)
     &             +(helibri(line,series)*wei)
             enddo
         enddo
c
        do  i = 1, atypes
        do series = 1, 2
        do line = 1,10 
        xhydroflux(line,series,i) = xhydroflux(line,series,i)
     &             +(xhydrobri(line,series,i)*wei)
        enddo
        enddo
        enddo
c     
c
         fheiilr = fheiilr+(heiilr*wei)
         h2qav = h2qav+(h2ql*wei)
         hei2qa = hei2qa+(hein(2)*wei)
         heii2qa = heii2qa+(heii2ql*wei)
         do atom = 1,atypes
            h2q(atom) = h2q(atom)+(h2qbri(atom)*wei)
         enddo
         fhbeta = fhbeta+(hbeta*wei)
         tlosac = tlosac+(eloss*wei)
         tforbi = tforbi+(floss*wei)
c
      end if
c     
c     ***INTEGRATION OF THE COLUMN DENSITY IN POPINT(6,11)
c     INTEGRATION OF PHOTON IMPACT PARAMETERS
c     INTEGRATION OF Hydrogen COLUMN (for dust)
c     INTEGRATION OF PAH COLUMN
c     RADIUS STEP SIZE : DR
c     
      if (imod .ne. 'REST') then
         wei2 = ((de*dh)*fi)*dvoluni
         dhavv = dhavv+wei2
         qhdha = qhdha+(wei2*qhdh)
         zetaeav = zetaeav+(wei2*zetae)
         wei = (dh*fi)*dr
         do j = 1, atypes
            do  i = 1, maxion(j)
               popint(i,j) = popint(i,j)+((wei*pop(i,j))*zion(j))
            enddo
c            popint(maxion(j)+1,j) = popint(maxion(j)+1,j)+wei*Zion(j)
         enddo
         if(grainmode) then       
           dustint = dustint+wei
           if(pahmode) then
            do i=1,pahi
             pahint(i) = pahint(i)+wei*pahZ(i)
            enddo
           endif
         endif
      end if
c     
      return 
      end


