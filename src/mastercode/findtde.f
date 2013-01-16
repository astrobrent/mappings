cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******FIND TEMPERATURE (TOIII,TNII) USING LINE RATIOS (ROIII,RNII)
c	AND EL. DENS.(DEOIII,DENII) IN COMMON BLOCK /TONLIN/
c
c	FIND DENSITIES (DESII,DEOII) USING LINE RATIOS (RSII,ROII)
c	AND TEMPERATURE (TSII,TOII) IN COMMON BLOCK /DENLIN/
c
c
c
      subroutine findtde()
      
c
c
c
      include 'cblocks.inc'
c
c
c
      double precision ao, bo, co, an, bn, cn, ratio
c
c
c
      double precision a, b, c, dens, tg, tgi
      double precision ra
      integer*4 k, kf
      double precision ftr, fde
c
c
      ftr(ratio,a,b) = b/(dlog(ratio)-dlog(a))
      fde(dens,tg,c) = 1.d0+((c*dens)/dsqrt((tg*(0.5d0+dsign(0.5d0
     &,tg)))+1.d-20))
c
      kf = 20
      ao = 7.7321d0
      bo = 3.2966d4
      co = 4.3764d-4
      an = 6.8662d0
      bn = 2.4995d4
      cn = 2.4037d-3
      toiii = 0.0d0
c
c    *** OIII
c
      tnii = 0.0d0
      ratio = roiii
      if (ratio.le.0.0d0) goto 500
      dens = deoiii
c
      tg = ftr(ratio,ao,bo)
c
      do 100 k = 1, kf
      tgi = tg
      ra = fde(dens,tg,co)*ratio
      tg = ftr(ra,ao,bo)
      if ((dabs(tg-tgi)/tg).lt.0.001d0) goto 110
  100 continue
  110 toiii = tg
c
c    *** NII
c
  500 continue
      ratio = rnii
      if (ratio.le.0.0d0) goto 1000
      dens = denii
c
      tg = ftr(ratio,an,bn)
      do 600 k = 1, kf
      tgi = tg
      ra = fde(dens,tg,cn)*ratio
      tg = ftr(ra,an,bn)
      if ((dabs(tg-tgi)/tg).lt.0.001d0) goto 610
  600 continue
  610 tnii = tg
c
 1000 continue
      return 
      end
