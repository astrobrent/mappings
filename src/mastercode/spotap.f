cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******CALCULATES FRACTION OF PHOTONS PRODUCED BY RECOMB.
c	OF HEII AND HEIII THAT IONISES HI
c	(ON THE SPOT APPROX. FOR THE DIFFUSE FIELD PHOTONS)
c
c     **REFERENCE : OSTERBROCK,D.E. (1974)
c	            ASTROPHYSICS OF GAZEOUS NEBULAE (FREEMAN)
c
c	            HUMMER,D.G. AND SEATON,M.J. (1964)
c	               MNRAS,V.127,P.230
c
c
c
      subroutine spotap(de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision a21,ahe0,alr,ceffa,ceffb,de,de2,dep
      double precision dey,dh,dm,ds0,ds1,ds2,ds3,ds4
      double precision ds5,fa,fhi,fp,fs,henu2,hnu2
      double precision ph,ph2,qp1,qp3,qs1,t,u,w
      double precision x,x1,xx,y1,y2,y3,yh,z
      double precision zz
      integer*4 n2
c
      double precision pol5,qss,qsp,qsp3,rcf
c
      pol5(x) = dmax1(0.0d0,((((-6.5289d0*x+41.54554d0)*x-
     &97.135778d0)*x+97.0517d0)*x-32.02831d0)*x-0.045645d0)
      qss(u) = 1.d-8*pol5(u ** 0.33333d0)
      qsp(u) = 5.73d-9*pol5((0.175d0+(0.2d0*u)) ** 0.33333d0)
      qsp3(u) = 4.21d-7*pol5((0.165d0+(0.2d0*u)) ** 0.33333d0)
      rcf(u,z) = 0.6412d0*((u/(z*z)) ** (-0.1146d0))
c
      yh = 0.0d0
      ph = 0.0d0
      ph2 = 0.0d0
      dey = 0.0d0
      dep = 0.0d0
      de2 = 0.0d0
c
      u = t/1.d4
      if (u.lt.0.1d0) u = 0.1d0
      if (u.gt.10.d0) u = 10.d0
      w = t/4.d4
      if (w.lt.0.1d0) w = 0.1d0
      if (w.gt.10.d0) w = 10.d0
      ds0 = 13.595d0
      dm = epsilon
      if (zion(2).lt.1.d-4) goto 1000
c
c    ***PH : FRACTION OF RECOMB. OF HEII TO EXCITED STATES
c	ABSORBED BY HI
c
      ahe0 = (dh*zion(2))*pop(1,2)
      a21 = 1.27d-4
      qs1 = qss(u)
      qp1 = qsp(u)
      qp3 = qsp3(u)
      n2 = 1.0d0/(a21+(de*(qs1+qp1)))
      fa = n2*a21
      fs = (de*n2)*qs1
      fp = (de*n2)*qp1
c
      ph = ((1.d0/6.d0)+(0.75d0*(fa+fp)))+(0.56d0*((1.d0/12.d0)+(
     &0.75d0*fs)))
      ds1 = 19.72d0-ds0
      ds2 = (16.0d0-0.92d0)-ds0
      ds3 = (21.21d0-1.49d0)-ds0
      ds4 = 21.21d0-ds0
      ds5 = 16.0d0-ds0
      dep = ev*(((0.75d0)*(((fa*ds1)+((0.56d0*fs)*ds2))+(fp
     &*ds3)))+((0.25d0)*(((2.d0/3.d0)*ds4)+((0.56d0/3.d0)*ds5))
     &))
  100 continue
c
c    ***YH : FRACTION OF RECOMB. OF HEII TO GROUND STATE
c	ABSORBED BY HI
c
      hnu2 = 1.234d20
      henu2 = 7.83d20
      y1 = (dh*fhi)*hnu2
      y2 = ahe0*henu2
      y3 = y1+y2
      if (y3.lt.dm) goto 200
      yh = y1/y3
      ds1 = 24.46d0-ds0
      dey = (ev*ds1)*yh
  200 continue
c
c    ***PH2 : FRACTION OF RECOMB. OF HEIII TO EXCITED STATES
c	ABSORBED BY HI
c
c     **EFFICIENCY OF PROCESS A & B
      ceffa = 0.1d0
      ceffb = 0.85d0
c
c      *A : HEII LYMAN PHOTONS
      x1 = 0.7179d0*(u ** (-0.0584d0))
c      *B : HEII TWO-QUANTA EMISSION
      ph2 = x1*ceffa
      xx = (1.0d0-x1)*1.425d0
c      *C : HEII BALMER CONTINUUM
      ph2 = ph2+(xx*ceffb)
      alr = 0.2946d0*(w ** 0.1972d0)
      ph2 = ph2+alr
      zz = 4.0d0
      ds1 = 30.d0-ds0
      ds2 = 15.d0-ds0
      ds3 = rcf(u,zz)
      de2 = ev*((((x1*ceffa)*ds1)+((xx*ceffb)*ds2))+(alr
     &*ds3))
c
c
  300 continue
 1000 continue
c
      return 
      end
