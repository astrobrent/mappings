cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE HEATING DUE TO PHOTOIONISATION
c	USING IONISATION FRACTION OF HYDROGEN, IT INTERPOLATES 
c	THE RATES HEAPH(5,5,11) CALCULATED IN SUBR. PHION
c
c
c
      subroutine pheat(de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision aa,bb,cc,del,fhii,fhmi
      double precision qval,temp,val
      double precision x1,xc,xlast,xlo,y1,y2,y3,ztot
c
      double precision de, dh , quad
c
      integer*4 i,i1,i2,i3,ic,iel,ion
c
      quad(del,y1,y2,y3) = y1+((del*((((4.0d0*y2)-(3.0d0*y1))-y3
     &)+(((y3-(2.0d0*y2))+y1)*del)))/2.0d0)
c
      pgain = 0.0d0
      fhmi = 2.d-5
      ztot = 1.0d0
      do 10 i = 2, atypes
 10      ztot = ztot+zion(i)
c
c    ***USES INTERPOLATION , FHII < 1.0
c
      fhii = dmin1(1.0d0,dmax1((de/dh)/ztot,fhmi))
      if (fhii.lt.1.0d0) then
         xlo = dlog10(fhii+epsilon)
c         xlast = -6.0d0
         xlast = -4.0d0
         x1 = dint(xlo)
         xc = x1
         if (xc.lt.xlast) xc = xlast
         if (x1.lt.(xlast+2.0d0)) x1 = xlast+2.0d0
         del = x1-xlo
         i1 = idnint(1.0d0+dabs(x1))
         i2 = i1+1
         i3 = i1+2
         ic = idnint(1.0d0+dabs(xc))
c
      do 101 iel = 1, limph
         do 100 ion = 1, maxion(iel)-1
            aa = heaph(i1,ion,iel)
            bb = heaph(i2,ion,iel)
            cc = heaph(i3,ion,iel)
            if (cc.gt.0.0d0) then
               qval = dexp(quad(del,dlog(aa),dlog(bb),dlog(cc)))
               val = dmax1(0.0d0,dmin1(qval,heaph(ic,ion,iel)))
               temp = ((zion(iel)*pop(ion,iel))*val)
               pgain = pgain+temp
            else if (aa.gt.0.0d0) then
               qval = quad(del,aa,bb,cc)
               val = dmax1(0.0d0,dmin1(qval,heaph(ic,ion,iel)))
               temp = ((zion(iel)*pop(ion,iel))*val)
               pgain = pgain+temp
            end if
 100     continue
 101  continue
      else
c
c    ***NO INTERPOLATION BECAUSE FHII = 1.0
c
         do 301 iel = 1, limph
            do 300 ion = 1, maxion(iel)-1
               pgain = pgain+(zion(iel)*pop(ion,iel)*heaph(1,ion,iel))
 300        continue
 301     continue
      end if
c     
      pgain = pgain*dh
c
      return 
      end

      
