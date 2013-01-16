cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******OUTPUT EMISSION SPECTRUM ON DEVICE*LUOP
c	USES SUMMED LINE DATA IN /SLINE/
c	ASSUMED INTENSITY IS RELATIVE TO H-BETA 
c	WHEN MODE='REL'
c
c
      subroutine spectrum(luop, mode)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision chcksum,chklim
c
      integer*4 i,j,jl,k,max,atom
      integer*4 luop
c
      character mode*4
c
c     old
c      double precision eloslog,fhbtlog,forblog
c
c	Formats:
c
c
 1    format(/,t6,'Lambda (A)',t18,5(1pg12.6))
 2    format(/,t6,'Lambda (A)',t18,5(1pg12.5))
 3    format(/,t6,'Lambda (m)',t18,5(1pg12.5))
 4    format( t5,'Int./H-beta',t18,5(1pg12.4))
 5    format(' ',1pg12.6,1x,1pg12.6,1x,1pg12.4,1x,a3,a6)
c
c	5, 6, 9 level atom model formats and fe II (11 & 12)
c
 10   format(/,t2,'Ion: ',a3,a6,$)
 11   format(/,t6,  'Transition',t18,5a12,$)
 12   format(/,t4,'Ionic specie',t18,5(2x,a4,a6),$)
c
c     chklim is a lower limit for line fluxes to print out.
c     1e-4 saves a lot of paper!
c
c      chklim = 0.d0
c
      chklim = epsilon
      if (mode .eq. 'REL') then
         chklim = 1.d-6
      endif
c     
c      if (mode .eq. 'REL') then
c
c    ***CORRECTS TOTAL H-BETA VALUE FOR THE CHOICE OF VUNIT
c	 EMISSIVITY CHANGED INTO LUMINOSITY : JNU*4PI*VOL
c
c      fhbtlog = dlog10((fhbeta*(4.d0*pi))+epsilon)+vunilog
c      forblog = dlog10(tforbi+epsilon)+vunilog
c      eloslog = dlog10(tlosac+epsilon)+vunilog
cc
c 2400 format(/' Flux Scale (H Beta = 1.000) :'/
c     &        ' -----------------------------'/
c     &' H-beta :',f8.4,
c     &' (log(ergs/cm^2/s))')
c
c 2401 format(/' Flux Scale (H Beta = 1.000) :'/
c     &        ' -----------------------------'/
c     &' H-beta  :',f8.4,
c     &' (log(ergs/s))')
c
c      if (jgeo.eq.'S') then
c      write(luop, 2401) fhbtlog
c      else
c      write(luop, 2400) fhbtlog
c      endif
c
c      else
c
c    'ABS'
c    ** NO CORRECTION FOR FACTOR 4PI WITH THIS MODE
c
c      fhbtlog = dlog10(fhbeta+epsilon)+vunilog
c      write(luop, 2425) fhbtlog
c 2425 format(/' Flux Scale (Absolute units) :'/
c     &        ' -----------------------------'/
c     &        ' Units (Log) :',f8.4)
c      end if
c
c     Top Twenty Five Lines - ie rough spectrum, include HBeta printout now too.
c
      call spec2(luop, 'TTWN', mode)
c
c     The full spectrum....
c
      write(luop,2402)
 2402 format(//' Detailed Emission Spectrum :'/
     &         ' ----------------------------')
c
      write(luop, 100) h2qav,h2q(1)
100   format(//,' Hydrogen Spectrum (recombination + collisions)'/
     &          ' ----------------------------------------------'/
     &          ' Old  HI Two-Photon emission :',1pg10.3/
     &          ' New  HI Two-Photon emission :',1pg10.3)
c
101   format(/,' Lyman Series'/
     &         ' ------------')
      write(luop, 101) 
c
      do i = 1,10
      	if (hydroflux(i,1).lt.chklim) hydroflux(i,1) = 0.d0
      enddo
c      
      write(luop,10) elem(1), rom(1)
      write(luop, 1) (hlambda(j,1),j=1,5)
      write(luop, 4) (hydroflux(j,1),j = 1, 5)
c
      write(luop, 1) (hlambda(j,1),j=6,10)
      write(luop, 4) (hydroflux(j,1),j = 6, 10)
c
      
102   format(/,' Balmer Series'/
     &         ' -------------')
      write(luop, 102) 
c
      write(luop,10) elem(1), rom(1)
      write(luop, 1) (hlambda(j,2),j=1,5)
      write(luop, 4) (hydroflux(j,2),j = 1, 5)
c
      write(luop, 1) (hlambda(j,2),j=6,10)
      write(luop, 4) (hydroflux(j,2),j = 6, 10)
c
      
103   format(/,' Paschen Series'/
     &         ' --------------')
      write(luop, 103) 
c
      write(luop,10) elem(1), rom(1)
      write(luop, 1) (hlambda(j,3),j=1,5)
      write(luop, 4) (hydroflux(j,3),j = 1, 5)
c
      write(luop, 1) (hlambda(j,3),j=6,10)
      write(luop, 4) (hydroflux(j,3),j = 6, 10)
c
      
104   format(/,' Brackett Series'/
     &         ' ---------------')
      write(luop, 104) 
c
      write(luop,10) elem(1), rom(1)
      write(luop, 1) (hlambda(j,4),j=1,5)
      write(luop, 4) (hydroflux(j,4),j = 1, 5)
c
      write(luop, 1) (hlambda(j,4),j=6,10)
      write(luop, 4) (hydroflux(j,4),j = 6, 10)
c
      
105   format(/,' Pfund Series'/
     &         ' ------------')
      write(luop, 105) 
c
      write(luop,10) elem(1), rom(1)
      write(luop, 1) (hlambda(j,5),j=1,5)
      write(luop, 4) (hydroflux(j,5),j = 1, 5)
c
      write(luop, 1) (hlambda(j,5),j=6,10)
      write(luop, 4) (hydroflux(j,5),j = 6, 10)
c
      
106   format(/,' Humphreys Series'/
     &         ' ----------------')
      write(luop, 106) 
c
      write(luop,10) elem(1), rom(1)
      write(luop, 1) (hlambda(j,6),j=1,5)
      write(luop, 4) (hydroflux(j,6),j = 1, 5)
c
      write(luop, 1) (hlambda(j,6),j=6,10)
      write(luop, 4) (hydroflux(j,6),j = 6, 10)
c
      if (zmap(2).ne.0) then
c
      write(luop, 300) heii2qa,h2q(2)
300   format(//,' Helium II Spectrum (recombination + collisions)'/
     &          ' -----------------------------------------------'/
     &          ' Old  HeII Two-Photon emission :',1pg10.3/
     &          ' New  HeII Two-Photon emission :',1pg10.3)
c
301   format(/,' n = 1 Series'/
     &         ' -------------'/
     &         ' Bowen Total :', 1pg10.3)
      write(luop, 301) heiioiiibfsum
c
      do i = 1,10
      	if (heliflux(i,1).lt.chklim) heliflux(i,1) = 0.d0
      enddo
c      
      write(luop,10) elem(2), rom(2)
      write(luop, 1) (helambda(j,1),j=1,5)
      write(luop, 4)  (heliflux(j,1),j = 1, 5)
c
      write(luop, 1) (helambda(j,1),j=6,10)
      write(luop, 4)  (heliflux(j,1),j = 6, 10)
c
      
302   format(/,' n = 2 Series'/
     &         ' ------------')
      write(luop, 302) 
c
      write(luop,10) elem(2), rom(2)
      write(luop, 1) (helambda(j,2),j=1,5)
      write(luop, 4)  (heliflux(j,2),j = 1, 5)
c
      write(luop, 1) (helambda(j,2),j=6,10)
      write(luop, 4)  (heliflux(j,2),j = 6, 10)
c
      
303   format(/,' n = 3 Series'/
     &         ' ------------')
      write(luop, 303) 
c
      write(luop,10) elem(2), rom(2)
      write(luop, 1) (helambda(j,3),j=1,5)
      write(luop, 4)  (heliflux(j,3),j = 1, 5)
c
      write(luop, 1) (helambda(j,3),j=6,10)
      write(luop, 4)  (heliflux(j,3),j = 6, 10)
c
      
304   format(/,' Pickering Series'/
     &         ' ----------------')
      write(luop, 304) 
c
      write(luop,10) elem(2), rom(2)
      write(luop, 1) (helambda(j,4),j=1,5)
      write(luop, 4)  (heliflux(j,4),j = 1, 5)
c
      write(luop, 1) (helambda(j,4),j=6,10)
      write(luop, 4)  (heliflux(j,4),j = 6, 10)
c
      
305   format(/,' n = 5 Series'/
     &         ' ------------')
      write(luop, 305) 
c
      write(luop,10) elem(2), rom(2)
      write(luop, 1) (helambda(j,5),j=1,5)
      write(luop, 4)  (heliflux(j,5),j = 1, 5)
c
      write(luop, 1) (helambda(j,5),j=6,10)
      write(luop, 4)  (heliflux(j,5),j = 6, 10)
c
      
306   format(/,' n = 6 Series'/
     &         ' ------------')
      write(luop, 306) 
c
      write(luop,10) elem(2), rom(2)
      write(luop, 1) (helambda(j,6),j=1,5)
      write(luop, 4)  (heliflux(j,6),j = 1, 5)
c
      write(luop, 1) (helambda(j,6),j=6,10)
      write(luop, 4)  (heliflux(j,6),j = 6, 10)     
c
      write(luop, 2495) hei2qa
2495  format(//,' Helium I Spectrum'/
     &          ' ------------------------------------'//
     &          ' HeI Two-Photon emission :',1pg10.3)
c
      write(luop, 10) elem(2), rom(1)
      write(luop,2) ((heilam(j)*1.d8),j = 1 , 5)
      write(luop,4) (fluxhei(j),j = 1, 5)
      write(luop,2) ((heilam(j)*1.d8),j = 6 , 10)
      write(luop,4) (fluxhei(j),j = 6, 10)
      write(luop,2) ((heilam(j)*1.d8),j = 11 , 15)
      write(luop,4) (fluxhei(j),j = 11, 15)
c      
      endif
c
      if (atypes.gt.2) then
         write(luop, 400) 
 400     format(//,' Heavy Atom Two Photon Continuum'/
     &             ' -------------------------------')
 401     format(//,' ',a2,' Two-Photon emission :',1pg10.3)
         do atom = 3,atypes
            if (h2q(atom).gt.chklim) then
               write(luop, 401) elem(atom),h2q(atom)
            endif
         enddo
c
c
      write(luop, 500)
500   format(//,' Heavy Atom Hydrogenic Spectrum ',
     &          '(recombination + collisions)'/
     &          ' -------------------------------',
     &          '----------------------------')
c
501   format(/,a2,' n = 1 Series'/
     &         ' -------------')
502   format(/,a2,' n = 2 Series'/
     &         ' ------------')
c
      do atom = 3 , atypes
c
      if (xhydroflux(1,1,atom).gt.chklim) then
c
      write(luop, 501) elem(atom)
c
      do i = 1,10
      	if (xhydroflux(i,1,atom).lt.chklim) xhydroflux(i,1,atom) = 0.d0
      enddo
c      
      write(luop,10) elem(atom), rom(mapz(atom))
      write(luop, 1) (xhlambda(j,1,atom),j=1,5)
      write(luop, 4) (xhydroflux(j,1,atom),j = 1, 5)
c
      write(luop, 1) (xhlambda(j,1,atom),j=6,10)
      write(luop, 4) (xhydroflux(j,1,atom),j = 6, 10)
c
      
      write(luop, 502) elem(atom)
c
      write(luop,10) elem(atom), rom(mapz(atom))
      write(luop, 1) (xhlambda(j,2,atom),j=1,5)
      write(luop, 4) (xhydroflux(j,2,atom),j = 1, 5)
c
      write(luop, 1) (xhlambda(j,2,atom),j=6,10)
      write(luop, 4) (xhydroflux(j,2,atom),j = 6, 10)
c
      endif
c
c     end atom loop
c
      enddo
c
c     end atypes > 2
c
      endif
c     
      if (n6ions.ne.0) then
c
c
      write(luop, 3427) 
 3427 format(//,' Six Level Atom Models'/
     &          ' ---------------------')
      chcksum = 0.0d0
      do i = 1, n6ions
         do j = 1, n6trans
            chcksum = chcksum + fluxf6(j,i)
         enddo
      enddo
 3428 format(/,' Total Six Level Line Int/H-Beta:',1pg14.6)
      write(luop,3428) chcksum
      chcksum = 0.0d0
      do i = 1, n6ions
         do j = 1, n6trans
            if (f6lam(j,i).lt.(rydlambda*1.d-8)) then
               chcksum = chcksum + fluxf6(j,i)
            endif
         enddo
      enddo
 3429 format(/,' Six Level Lines < 911A Int/H-Beta:',1pg14.6)
      write(luop,3429) chcksum
      do i = 1, n6ions
         chcksum = 0.0d0
         do j = 1, n6trans
            chcksum = chcksum + fluxf6(j,i)
         enddo
         if (chcksum.ge.chklim) then 

            write(luop,10) elem(f6atom(i)), rom(f6ion(i))
            write(luop, 2) ((f6lam(j,i)*1.d8),j=1,5)
            write(luop, 4)  (fluxf6(j,i),j = 1, 5)

            write(luop, 2) ((f6lam(j,i)*1.d8),j=6,10)
            write(luop, 4)  (fluxf6(j,i),j = 6,10)
            
            write(luop, 2) ((f6lam(j,i)*1.d8),j=11,15)
            write(luop, 4)  (fluxf6(j,i),j = 11,15)

        endif
      enddo
c
c
      endif
c
      if (nfions.ne.0) then
c
      write(luop, 2427) 
 2427 format(//,' Five Level Atom Models'/
     &          ' ----------------------')
      chcksum = 0.0d0
      do i = 1, nfions
         do j = 1, nftrans
            chcksum = chcksum + fluxf(j,i)
         enddo
      enddo
 2428 format(/,' Total Five Level Line Int/H-Beta:',1pg14.6)
      write(luop,2428) chcksum
      chcksum = 0.0d0
      do i = 1, nfions
         do j = 1, nftrans
            if (flam(j,i).lt.(rydlambda*1.d-8)) then
               chcksum = chcksum + fluxf(j,i)
            endif
         enddo
      enddo
 2429 format(/,' Five Level Lines < 911A Int :',1pg14.6)
      write(luop,2429) chcksum
      do i = 1, nfions
         chcksum = 0.0d0
         do j = 1, nftrans
            chcksum = chcksum + fluxf(j,i)
         enddo
         if (chcksum.ge.chklim) then 

            write(luop,10) elem(fatom(i)), rom(fion(i))
            write(luop, 2) ((flam(j,i)*1.d8),j=1,5)
            write(luop, 4)  (fluxf(j,i),j = 1, 5)

            write(luop, 2) ((flam(j,i)*1.d8),j=6,10)
            write(luop, 4)  (fluxf(j,i),j = 6,10)

         endif
      enddo
c
      endif
c
      if (nf3ions.ne.0) then
c
      write(luop, 2327) 
 2327 format(//,' Three Level Atom Models'/
     &          ' -----------------------'//)
c
      chcksum = 0.0d0
      do i = 1, nf3ions
         do j = 1, nf3trans
            chcksum = chcksum + fluxf3(j,i)
         enddo
      enddo
 2328 format(/,' Total 3LA Lines Int/H-Beta:',1pg14.6)
      write(luop,2328) chcksum
c
      do i = 1, nf3ions
         chcksum = 0.0d0
         do j = 1, nf3trans
            chcksum = chcksum + fluxf3(j,i)
         enddo
         if (chcksum.ge.chklim) then 
            write(luop, 10) elem(f3atom(i)), rom(f3ion(i))
            write(luop, 3) (f3lam(j,i),j=1,nf3trans)
            write(luop, 4) (fluxf3(j,i),j = 1, nf3trans)
         endif
      enddo
c
      endif
c
      if ((mlines+xilines).ne.0) then
c
c
      max = (mlines+4)/5
      write(luop, 2503) 
 2503 format(//,' Two Level FS and Semi-Forbidden models'/
     &          ' --------------------------------------')
       chcksum = 0.0d0
      do i = 1, mlines
         chcksum = chcksum + fluxi(i)
      enddo
      do i = 1, xilines
         chcksum = chcksum + fluxxi(i)
      enddo
 2504 format(/,' Total Semi-Forbidden Line Int/H-Beta:',1pg14.6)
      write(luop,2504) chcksum
       chcksum = 0.0d0
      do i = 1, mlines
         if (fslam(i).lt.(rydlambda*1.d-8)) then 
            chcksum= chcksum + fluxi(i)
         endif
      enddo
      do i = 1, xilines
         if (xilam(i).lt.rydlambda) then
            chcksum = chcksum + fluxxi(i)
         endif
      enddo
 2505 format(/,' Semi-Forbidden Lines < 911A Int/H-Beta:',1pg14.6)
      write(luop,2505) chcksum
      do i = 1, max
         jl = 5
         k = (i-1)*5
         if ((i*5).gt.mlines) jl = mlines-k

      chcksum = 0.0d0
      do j = 1, jl
         chcksum = chcksum + fluxi(j+k)
      enddo

      if (chcksum.ge.chklim) then 

         write(luop,12) (elem(ielfs(k+j)),rom(ionfs(k+j)),j=1,jl)
         write(luop,2) ((fslam(k+j)*1.d8),j = 1,jl)
         write(luop,4) (fluxi(j+k),j = 1, jl)

      endif
      enddo
c
c
c
      max = (xilines+4)/5
      write(luop, 2553) 
 2553 format(//)
      do 2555 i = 1, max
      jl = 6
      k = (i-1)*5
      if ((i*5).gt.xilines) jl = xilines-k
      chcksum = 0.0d0
      do j = 1, jl
         chcksum = chcksum + fluxxi(j+k)
      enddo
      if (chcksum.ge.chklim) then 
         write(luop,12)(elem(xiat(k+j)),rom(xiion(k+j)),j = 1,jl)
         write(luop,2) ((xilam(k+j)),j = 1,jl)
         write(luop,4) (fluxxi(j+k),j = 1, jl)
      endif
 2555 continue
c
      endif
c
c
c
      if ((nlines+xlines).ne.0) then
c
      max = (nlines+4)/5
      write(luop, 2603) 
 2603 format(//' Two Level Resonance Lines'/
     &         ' -------------------------')
      chcksum = 0.0d0
      do i = 1, nlines
         chcksum = chcksum + fluxr(i)
      enddo
      do i = 1, xlines
         chcksum = chcksum + fluxx(i)
      enddo
 2604 format(/,' Total Resonance Line Int/H-Beta:',1pg14.6)
      write(luop,2604) chcksum
      chcksum = 0.0d0
      do i = 1, nlines
         if (rlam(i).lt.(rydlambda*1.d-8)) then
            chcksum = chcksum + fluxr(i)
         endif
      enddo
      do i = 1, xlines
         if (xrlam(i).lt.rydlambda) then
            chcksum = chcksum + fluxx(i)
         endif
      enddo
 2605 format(/,' Resonance Lines < 911A Int/H-Beta:',1pg14.6)
      write(luop,2605) chcksum

      do i = 1, max
         jl = 5
         k = (i-1)*5
         if ((i*5).gt.nlines) jl = nlines-k
         chcksum = 0.0d0
         do j = 1, jl
            chcksum = chcksum + fluxr(j+k)
         enddo

         if (chcksum.ge.chklim) then 

            write(luop,12)(elem(ielr(k+j)),rom(ionr(k+j)),j = 1,jl)
            write(luop,2) ((rlam(k+j)*1.d8),j = 1,jl)
            write(luop,4) (fluxr(j+k),j = 1, jl)

         endif
      enddo
c     
c
c
      max = (xlines+4)/5
      write(luop, 2703) 
 2703 format(//)
      do 2705 i = 1, max
         jl = 5
         k = (i-1)*5
         if ((i*5).gt.xlines) jl = xlines-k
         chcksum = 0.0d0
         do j = 1, jl
            chcksum = chcksum + fluxx(j+k)
         enddo
         if (chcksum.ge.chklim) then 
            write(luop,12) (elem(xrat(k+j)),rom(xion(k+j)),j = 1,jl)
            write(luop,2) ((xrlam(k+j)),j = 1,jl)
            write(luop,4) (fluxx(j+k),j = 1, jl)
         endif
 2705 continue
c
      endif
c
c
      if (xhelines.ne.0) then
c
      max = (xhelines+4)/5
      write(luop, 2803) 
 2803 format(//' He like ion Intersystem and Forbidden Lines'/
     &         ' -------------------------------------------')
      do i = 1, xhelines
         chcksum = chcksum + fheif(i)
      enddo
 2804 format(/,' Total He X-UV Line Int/H-Beta:',1pg14.6)
      write(luop,2804) chcksum
      do i = 1, max
         jl = 4
         k = (i-1)*4
         if ((i*4).gt.xhelines) jl = xhelines-k
         chcksum = 0.0d0
         do j = 1, jl
            chcksum = chcksum + fheif(j+k)
         enddo
         if (chcksum.ge.chklim) then 
            write(luop,12) (elem(xheat(k+j)),rom(xheion(k+j)),j = 1,jl)
            write(luop,2) ((xhelam(k+j)),j = 1,jl)
            write(luop,4) (fheif(j+k),j = 1, jl)
         endif
      enddo
c
      endif
c
      if (n9ions.ne.0) then
c
      write(luop, 3527) 
 3527 format(//' Fe VII Multiplet Lines'/
     &         ' ----------------------')
      chcksum = 0.0d0
      do i = 1, n9ions
         do j = 1, n9trans
            chcksum = chcksum + fluxf9(j,i)
         enddo
      enddo
 3528 format(/,' Total Fe VII Line Int/H-Beta:',1pg14.6)
      write(luop,3528) chcksum

      do i = 1, n9ions
      
         chcksum = 0.0d0
         do j = 1, n9trans
            chcksum = chcksum + fluxf9(j,i)
         enddo
         if (chcksum.ge.chklim) then 

            write(luop, 10) elem(f9atom(i)), rom(f9ion(i))
            write(luop, 2) ((f9lam(j,i)*1.d8),j=1,5)
            write(luop, 4) (fluxf9(j,i),j = 1, 5)
            
            write(luop, 2) ((f9lam(j,i)*1.d8),j=6,10)
            write(luop, 4) (fluxf9(j,i),j=6,10)

            write(luop, 2) ((f9lam(j,i)*1.d8),j=11,15)
            write(luop, 4) (fluxf9(j,i),j=11,15)

            write(luop, 2) ((f9lam(j,i)*1.d8),j=16,20)
            write(luop, 4) (fluxf9(j,i),j=16,20)

            write(luop, 2) ((f9lam(j,i)*1.d8),j=21,25)
            write(luop, 4) (fluxf9(j,i),j=21,25)

            write(luop, 2) ((f9lam(j,i)*1.d8),j=26,30)
            write(luop, 4) (fluxf9(j,i),j=26,30)
            
            write(luop, 2) ((f9lam(j,i)*1.d8),j=31,35)
            write(luop, 4) (fluxf9(j,i),j=31,35)

            write(luop, 2) (f9lam(36,i)*1.d8)
            write(luop, 4) fluxf9(36,i)
            
        endif
      enddo
c
      endif
c
      if (zmap(26).ne.0) then
c
      max = (nfetrans+4)/5
      write(luop, 2903) 
 2903 format(//' Fe II Multiplet Lines'/
     &         ' ---------------------')
      chcksum = 0.0d0
      do i = 1, nfetrans
         chcksum = chcksum + fluxfe(i)
      enddo
 2904 format(/,' Total FeII Int/H-Beta:',1pg14.6)
      write(luop,2904) chcksum
      chcksum = 0.0d0
      do i = 1, nfetrans
         if (felam(i).lt.0.3) then
            chcksum = chcksum + fluxfe(i)
         endif
      enddo
 2905 format(/,' UV FeII lines < 3000A Int/H-Beta:',1pg14.6)
      write(luop,2905) chcksum
      chcksum = 0.0d0
      do i = 1, nfetrans
         if (felam(i).gt.0.9) then
            chcksum = chcksum + fluxfe(i)
         endif
      enddo
 2906 format(/,' IR FeII lines > 9000A Int/H-Beta:',1pg14.6)
      write(luop,2906) chcksum
c
      do  i = 1, max
         jl = 5
         k = (i-1)*5
         if ((i*5).gt.nfetrans) jl = nfetrans-k
         chcksum = 0.0d0
         do j = 1, jl
            chcksum = chcksum + fluxfe(felamap(j+k))
         enddo
         if (chcksum.ge.chklim) then 
            write(luop,11 ) 
     &           (feid(felamap(k+j)),j = 1,jl)
            write(luop,3) 
     &           ((felam(felamap(k+j))),j = 1,jl)
            write(luop,4)
     &           (fluxfe(felamap(j+k)),j = 1, jl)
         endif
      enddo
c
      endif
c
      write(luop, 2475) 
 2475 format(//' Old Version Hydrogen Lines'/
     &         ' --------------------------'/)     
      write(luop, 10) elem(1), rom(1)
      write(luop,2) ((hlam(j)*1.d8),j = 1,5)
      write(luop,4) (fluxh(j),j = 1, 5)
c
c
      return 
      end
