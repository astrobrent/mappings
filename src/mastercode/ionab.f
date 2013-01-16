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
c     Version 2.0, expanded arrays to cope with higher ionisation
c     stages.
c
c     RSS 8/90
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******CALCULATES IONIC ABUNDANCES AFTER TIME : TSTEP
c	DEPOPULATION RATES : RECO(5) ; #OF IONIC SPECIES : NDE
c	REPOPULATION RATES : PION(5) ,PIONAU(4)
c	INITIAL ABUNDANCES IN AB(6)
c	RETURNS ABUNDANCES IN AB(6) AND DERIVATIVES IN ADNDT(6)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ionab(reco,pion,pionau,ab,adndt,nde,tstep)
c
      include 'cblocks.inc'
c     
c     
c     Variables
c
      double precision astem,cyn,dyn
      double precision repoprate
      double precision rm,ryn,tr,tstep
c     
      double precision reco(mxion), pion(mxion), pionau(mxion)
      double precision ab(mxion), b(mxion, mxion+1), c(mxion, mxion)
      double precision s(mxion, mxion), as(mxion, 2)
      double precision cm(mxion, mxion), adndt(mxion)
      double precision ratetotal
      double precision dtotal, da, dma, delt, tim, dun, rn
c     
      integer*4 i, j, m1, nde
      integer*4 nit, n, nn, nitmi, nitma, kk
      integer*4 j1, j2, k, l, m, maxio
c
      character ndyn*4
c
c    ***CALCULATES TOTAL ABUND. : DTOTAL
c
      dun = 1.00000001d0
c
c     find total abundance
c
      dtotal = 0.d0
      ratetotal=0.d0
      do k = 1, nde
         dtotal = ab(k)+dtotal
         ratetotal=reco(k)+pion(k)+pionau(k)
      enddo
c      if (ratetotal.eq.0.d0) then
c     rates are 0 so nothing changes, leave ab the same
c        do k = 1, nde
c          adndt(k) = 0.0d0
c        enddo
c        write(*,1069) nde,dtotal,ratetotal
c        return
c      endif
      maxio = nde
c     
c     ***ZERO MATRICES C,B,S,AS
c     
      do l = 1, maxio
         as(l,1) = 0.d0
         as(l,2) = 0.d0
         do k = 1, maxio
            c(k,l) = 0.d0
            cm(k,l) = 0.0d0
            s(k,l) = 0.d0
            b(k,l) = 0.d0
         enddo
      enddo
c     
c     ***PUT INITIAL ABUNDANCES IN AS(K,1)
c     
      do  k = 1, nde
         as(k,1) = ab(k)
      enddo
c     
c     ***FORM MATRIX OF IONISATION RATES IN C
c     
      ndyn = 'TEST'
 424  do k = 1, nde-1
         c(k,k+1) = reco(k)
         c(k+1,k+1) =-c(k,k+1)
      enddo
      
      c(1,1) = 0.d0
      do  k = 1, nde-1
         c(k+1,k) = pion(k)
         c(k,k) = c(k,k)-c(k+1,k)
      enddo
c     
c     ***FIND MAX. OF ABS( C )
c     
      do  k = 1, nde-2
         c(k+2,k) = pionau(k)
         c(k,k) = c(k,k)-c(k+2,k)
      enddo
      da = 0.d0
      do k = 1, nde
         da = dmax1(da,dabs(c(k,k)))
      enddo
c     
c     ***FIND INTEGRATION TIME : DELT  (PER STEP)
c     THE POWER TO WHICH MATRIX  S  WILL BE RISEN : (NITMI*10)
c     AND THE NUMBER OF ITERATIONS : NIT
      dyn = 1.d13
c     dyn = 1.d14
      tim = tstep
      nitma = 100
      if (da.lt.(1.d-38*nitma)) da = 1.d-38*nitma
      delt = 1.d0/da
      if (tim.le.delt) goto 470
      if (tim.le.(nitma*delt)) goto 480
      nitma = 20
c     
c     **COMPRESS IONISING RATES IF DYNAMIC RANGE EXCEEDED
c     
      if ((((dyn*1.d-38)*delt)*nitma).gt.1.d0) goto 403
      tr = tim/(dyn*delt)
      if (tr.le.1.0d0) goto 403
      tim = dyn*delt
      if (ndyn .eq. 'OK') goto 403
      j = 0
      rm = 1.d38
      do 493 i = 1, nde-1
         tr = dlog(pion(i)+1.d-37)-dlog(reco(i)+1.d-37)
         if (dabs(tr).ge.rm) goto 493
         j = i
         rm = tr
 493  continue
c     
      if (j.lt.3) goto 403
      ryn = (reco(j)*dyn)/100.d0
      repoprate = pion(j)*1.d4
      do 417 i = 1, j-2
         tr = pion(i)
         cyn = reco(i)*1.d7
         if (((tr.le.ryn).or.(tr.lt.repoprate)).or.(tr.lt.cyn)) goto 417
         rm = ryn/pion(i)
         pion(i) = pion(i)*rm
         reco(i) = reco(i)*rm
         if (i.gt.1) pionau(i-1) = pionau(i-1)*rm
 417  continue
      ndyn = 'OK'
      goto 424
c     
 403  continue
c     write (*,*) '403',tim,nitma,delt,dun
      delt = tim/(nitma*dint((tim/(nitma*delt))+dun))
      da = dlog10(tim/(nitma*delt))
      dma = da-dint(da)
c     write (*,*) da,dma
      dma = nitma*(10.0d0 ** dma)
      delt = delt*(dma/dint(dma+dun))
      da = dlog10(tim/(nitma*delt))
c     write(*,*) delt,dma,da
      dma = da-dint(da)
      nitmi = idint(da)
      nit = idnint(nitma*(10.0d0**dma))
c     write(*,*) '403',tim,delt,da,dma,nitmi,nit
      goto 450
 470  nitmi = 0
      delt = tim
      nit = 1
      goto 450
 480  nitmi = 0
      delt = tim/dint((tim/delt)+dun)
      nit = idint(tim/delt)+1
c     write(*,*) '480',delt,da,dma,nitmi,nit
 450  continue
c     
c     ***MULTIPLY MATRIX BY : DELT
c     
      do l = 1, nde
         do k = 1, nde
            cm(k,l) = c(k,l)
            c(k,l) = c(k,l)*delt
         enddo
      enddo
c     
c     ***SET S AND B AS IDENTITY MATRICES
c     
      do k = 1, nde
          b(k,k) = 1.d0
          s(k,k) = 1.d0
       enddo
c     
c     ***EACH TERM OF THE DEXPANSION IS FORMED IN B
c     AND ADDED TO S
c     
      rn = 0.d0
      da = 0.d0
      m1 = maxio+1
 510  do m = 1, 10
         rn = rn+1.d0
         do k = 1, nde
            do l = 1, nde
               b(l,m1) = 0.d0
               do n = 1, nde
                  b(l,m1) = b(l,m1)+(c(l,n)*b(n,k))
               enddo
            enddo
            do l = 1, nde
               b(l,k) = b(l,m1)/rn
               s(l,k) = s(l,k)+b(l,k)
            enddo
         enddo
      enddo
c     
c     
      do k = 1, nde
         do l = 1, nde
            dma = dmax1(dabs(b(l,k)),0.d0)
         enddo
      enddo
      da = dmax1(da,dma)
      if (dma.gt.0.d0) goto 510
c     
c     ***RAISE MATRIX S TO POWER : (NITMI*10)
c     
      if (nitmi.lt.1) goto 680
      do  nn = 1, nitmi
         do  k = 1, nde
            do  l = 1, nde
               c(l,k) = s(l,k)
            enddo
         enddo
         do  m = 2, 10
            do  k = 1, nde
               do  l = 1, nde
                  b(l,m1) = 0.d0
                  do  n = 1, nde
                     b(l,m1) = b(l,m1)+(c(l,n)*s(n,k))
                  enddo
               enddo
               do l = 1, nde
                  s(l,k) = b(l,m1)
               enddo
            enddo
         enddo
      enddo
c     
c     ***MULTIPLY VECTOR  AS  BY MATRIX  S   (NIT  TIMES)
c     
c     write (*,*) nit
 680  l = 0
      do 745 n = 1, nit
         l = 1-l
         j1 = 1+l
         j2 = 2-l
         do 750 k = 1, nde
            astem = 0.0d0
            do 755 kk = 1 ,maxio
               astem = astem+s(k,kk)*as(kk,j2)
 755        continue
            as(k,j1) = astem
 750     continue
 745  continue
c     
c     FIND NORMALISATION FACTOR : DA
c     
 950  da = 0.d0
      do 900 k = 1, nde
         if (as(k,j1).gt.0.0d0) da = da+as(k,j1)
 900  continue
      da = da/dtotal
      dma = 100.0d0*dabs(1.d0-da)
c     
c     ***NORMALISATION AND WARNING MESSAGE IF ABUNDANCES
c     DO NOT CONSERVE
c     
      if (dma.lt.1.d0) goto 930
      write(*, 940) dma
c     
 940  format(//' WARNING (ionab)!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     &     'ABUNDANCES DID NOT CONSERVE BEFORE NORMALISATION :',g9.3,
     &     ' %')
      write (*,*) 'Input fractions:'
      write (*,*) nde,ratetotal
      write (*,*) 'dtotal,ion1,ion2...'
      write (*,941) dtotal,(ab(k), k = 1,nde)
      write (*,*) 'Input recom rates:'
      write (*,941) (reco(k), k = 1,nde)
      write (*,*) 'Input pion rates:'
      write (*,941) (pion(k), k = 1,nde)
      write (*,*) 'calculated fractions:'
      write (*,*) 'da,ion1,ion2...'
      write (*,*) da,(as(k,j1), k = 1,nde)
 941  format (f12.8,1x,29(f12.8,1x))
 942  format (f12.8)
c     
      stop
c     
 930  do k = 1, nde
         ab(k) = as(k,j1)/da
         if (ab(k).lt.1.0d-38) ab(k) = 0.0d0
      enddo
c     
c     ***CALCULATES DN/DT FOR EACH IONIC SPECIES
c     
      do k = 1, nde
         adndt(k) = 0.0d0
         do  n = 1, nde
            adndt(k) = adndt(k)+(ab(n)*cm(k,n))
         enddo
      enddo
c
c      write(*,1069) nde,dtotal,ratetotal
c 1069 format(i3,1x,2(e12.4E2,1x))
      return 
      end

