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
c*******COMPUTES THE AVERAGE DISTANCE : DRTA  TO OBTAIN
c	A GIVEN ABSORBTION : TAUAV (MMOD='DIS','LIN')
c	OR DERIVES EFFECTIVE TAU  :  TAUAV  FOR A GIVEN 
c	DISTANCE : DRTA  (MMOD='TAU')
c
c	MMOD='DIS' IS HIGHLY NON LINEAR AND PRESUPPOSE EQUILIBRIUM
c	IONISATION , TAKES INTO EFFECT GEOMETRICAL DILUTION .
c
c	RAD : CURVATURE RADIUS AT THE POINT CONSIDERED
c	FI : FILLING FACTOR  ;  DH : PEAK H DENSITY
c	USES PHOTON FLUX : TPHOT  ,  IONIC POPULATIONS : POPUL(J,I)
c	AND PHOTOIONISATION CROSS SECTIONS OF ALL ELEMENTS
c	UP TO*LIMPH
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine taudist(dh, fi, tauav, drta, rad, popul, mmod)
c
      
      include 'cblocks.inc'
c
c           Variables
c
      double precision a,abio,adto,ar,at,bet,cebin,cmpos
      double precision cof,convg,crosec,ctpos,divsi,dri,ecr,eph
      double precision expo,fg,gam,gammax,geo
      double precision phma,pos,relma,relp
      double precision rewe,rr,sc,sca,sca2,scaw,scaw2,se
      double precision sig,sigm,sigmgeo,simax,steep,taunew,testab,u
      double precision ulo,wei,weiad,weim2,weima,wev,wlo
      double precision radius, sigwei, sigtot
      double precision sigen(mxinfph), popul(mxion, mxelem)
      double precision  relwei(mxinfph), rela(16)
      double precision dh, fi, tauav, drta, rad, far, acrs
c
      integer*4 i,ie,inl
      integer*4 j,k,lim,n
c
      character mmod*4
c
      far(u,a) = (a*(dexp(dmin1(24.0d0,(2.0d0*u)/a))-1.0d0))/(dexp(
     &dmin1(24.0d0,(2.0d0*u)/a))+1.0d0)
      acrs(eph,at,bet,se) = (at*(bet+((1.0d0-bet)/eph)))*(eph
     &**(-se))
c
      if (((((mmod .eq. 'DIS').or.(mmod .eq. 'LIN')).or.(mmod .eq. 
     &'TAU')).or.(mmod .eq. 'DPRI')).or.(mmod .eq. 'LPRI')) goto 60
      write(*, 67) mmod
   67 format(//'MODE IMPROPERLY DEFINED IN TAUDIST : ',a4)
      stop 
   60 continue
c
      fg = 0.66d0
      convg = 0.005d0
      sc = 1.d-34
      geo = 1.0d0
c
      if ((mmod .ne. 'DIS').and.(mmod .ne. 'DPRI')) then
         steep = 2.0d0
         expo = 1.0d0
         divsi = 1.0d0
         ctpos = 0.0d0
         cmpos = 0.4d0
         radius = 0.d0
      else
         steep = 150.0d0
         expo = 0.7d0
         divsi = 8.0d0
         ctpos = 0.15d0
         cmpos = 0.005d0
         radius = dabs(rad)
      end if
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (mmod .ne. 'TAU') drta = 0.0d0
      phma = 0.0d0
      relma = 0.0d0
      simax = 1.d-37
      do 500 inl = 1, infph-1
         if (skipbin(inl)) goto 500
         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
         sig = 0.0d0
         rewe = 0.0d0
         do 375 i = 1, limph
 375        rela(i) = 0.0d0
         do 290 i = 1, ionum
            ie = atpho(i)
            if (ie.gt.limph) goto 290
            j = ionpho(i)
            abio = zion(ie)*(popul(j,ie)**expo)
            eph = cebin/ipotpho(i)
            if (eph.lt.1.0d0) goto 293
            crosec = acrs(eph,sigpho(i),betpho(i),spho(i))
            sig = sig+(abio*crosec)
            if (arad(j+1,ie).gt.0.0d0) rela(ie) = 1.0d0
 290     continue
 293     continue
c
         sig = (dh*fi)*sig
c
         do 310 i = 1, limph
 310        rewe = rewe+(rela(i)*zion(i))
c
         phma = dmax1(phma,tphot(inl))
         simax = dmax1(simax,sig*(rewe**2.d0))
         relma = dmax1(relma,rewe)
         sigen(inl) = sig
         relwei(inl) = rewe
 500  continue
c
c    ***AVERAGES SIGMA OVER PHOTON SPECTRUM
c
      simax = simax/(relma ** 2.d0)
      sca = 1.0d0/dmax1(1.d-37,dsqrt(phma)*sc)
      sca2 = 1.0d0/dmax1(sc,(dsqrt(phma)*sca)*sc)
      weima = 0.0d0
      weim2 = 0.0d0
      do 780 inl = 1, infph-1
         if (skipbin(inl)) goto 780
         testab = (relwei(inl)/relma) ** 2.d0
         if (testab.gt.0.001d0) then
            wev = (ephot(inl+1)-ephot(inl))/(ephot(inl)+ephot(inl+1))
            ecr = dsqrt(sigen(inl)/simax)
            wei = ((dsqrt(tphot(inl))*sca)*sca2)*dsqrt(wev)
            relp = wei*testab
            relwei(inl) = relp
            weima = dmax1(weima,ecr*relp)
            weim2 = dmax1(weim2,relp)
         else
            relwei(inl) = 0.0d0
         end if
  780 continue
c
      scaw = 1.0d0/dmax1(1.d-35,1.d-35*weima)
      scaw2 = 1.0d0/dmax1(1.d-35,1.d-35*weim2)
      do 800 k = 1, 14
      pos = dmin1(1.d2,cmpos+(ctpos/(tauav+1.d-20)))
      sigtot = 0.d0
      sigwei = 0.d0
      gammax = 0.d0
c
      do 810 inl = 1, infph-1
         if (skipbin(inl)) goto 810
         ecr = dsqrt(sigen(inl)/simax)
         ar = drta*sigen(inl)
         gam = far(pos*ar,steep)
         gammax = dmax1(gam,gammax)     
         wei = relwei(inl)
         if ((wei.gt.0.0d0).and.(ecr.gt.0.0d0)) then
            adto = dexp(((dlog(wei)+dlog(ecr))+dlog(scaw))-gam)
            weiad = dexp((dlog(wei)+dlog(scaw2))-gam)
         else
            adto = 0.0d0
            weiad = 0.0d0
         end if
         sigwei = sigwei+weiad
         sigtot = sigtot+adto
  810 continue
c
c    ***DETERMINE AVERAGE TAUNEW
c
      sigwei = dmax1(1.d-37,sigwei)
      if (sigwei.gt.0.d0) then
         sigm = dmax1(1.d-37,simax*(dexp(((dlog(sigtot+1.d-37)-dlog(
     &   sigwei))-dlog(scaw))+dlog(scaw2)) ** 2))/divsi
      else
         sigm = 1.d-37
      end if
      taunew = (sigm*geo)*drta
c
c    ***IMPROVES VALUE OF DRTA IF MMOD='DIS','LIN'
c    TAKES INTO ACCOUNT THE AMPLIFICATION FACTOR DUE TO GEOM.
c
      if (mmod .eq. 'TAU') goto 1200
      lim = 5
  720 lim = lim+5
      if (radius .eq. 0.d0) lim = 2
      drta = tauav/sigm
      do 727 n = 1, lim
         dri = drta
         geo = 1.0d0
         if (radius .eq. 0.d0) goto 127
         wlo = (2.d0*dlog(radius))-(2.d0*dlog(radius+(fg*drta)))
         ulo = (2.d0*dlog(radius+(fg*drta)))-(2.d0*dlog(radius))
         if (rad.lt.0.0d0) geo = dexp(wlo)
         if (rad.gt.0.0d0) geo = dexp(ulo)
 127     continue
c
         u = dble(n)*(1.0d0/(dble(lim)-1.0d0))
         if (n .eq. lim) u = 1.0d0
         sigmgeo = sigm*(geo ** u)
         drta = tauav/sigmgeo
         cof = dabs(drta-dri)/dmax1(drta,dri,1.d-8)
  727 continue
c
      if ((cof.gt.convg).and.(lim.lt.26)) goto 720
      if ((mmod .eq. 'DPRI').or.(mmod .eq. 'LPRI'))
     & write(*,10) tauav, taunew, drta, geo,
     & gammax, lim, mmod
c
   10 format(1h ,5(1pg12.4),i4,1x,a4)
      rr = taunew/(tauav+1.d-20)
      if ((dabs(1.0d0-rr).lt.convg).and.(k.gt.2)) goto 900
      goto 800
c
c    ***REPLACES VALUE OF TAUAV USING NEW TAU IF MMOD='TAU'
c
 1200 continue
      tauav = taunew
c
c     *NOW EXIT
c
      goto 900
c
c     *LOOP BACK
c
  800 continue
c
  900 continue
c
      return 
      end

