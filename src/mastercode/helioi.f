cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******DERIVES INTENSITY OF IMPORTANT HEI LINES AND THE RESULTANT
c	USES FUNCTION FGAUNT
c	COOLING RATE : HELOS
c
c
      subroutine helioi(t, de, dh, helos)
c
      include 'cblocks.inc'
c
      double precision t, de, dh, helos, x, u
      double precision abup,ara,cop1,cos1,f,frp1,frs1
      double precision omep1,omes1,pi4,ra21,ra3,rap1,rap3,ras1
      double precision reche,tns3,xp1,xs1
      double precision a2S31S1
      double precision a2P31S1
      double precision at1s
      double precision a2S11S1
c      
      double precision ome1s12p3
      double precision ome1s12s3
      double precision ome1s12s1
      integer*4 atom
c
c           Functions
c
      double precision fgaunt,pol5
      double precision qss,qsp,qsp3,frs3,fusp
c
      pol5(x) = dmax1(0.0d0,(((((((((-(6.5289d0*x))+41.54554d0)*x)-
     &97.135778d0)*x)+97.0517d0)*x)-32.02831d0)*x)-0.045645d0)
      qss(u) = 1.d-8*pol5(u ** 0.3333333d0)
      qsp(u) = 5.73d-9*pol5((0.175d0+(0.2d0*u)) ** 0.3333333d0)
      qsp3(u) = 3.70d-7*pol5((0.218d0+(0.1d0*u)) ** 0.3333333d0)
      frs3(u) = 0.775d0*(u ** 0.0213d0)
      fusp(u) = 0.42d0*(u ** 0.23d0)
c
      u = t/1.d4
      if (u.lt.0.1d0) u = 0.1d0
      if (u.gt.10.d0) u = 10.d0
c
      f     = dsqrt(1.0d0/(t+epsilon))
c
      helos = 0.0d0
c      
      atom  = zmap(2)
      pi4   = 4.d0*pi
c      
      a2S31S1   = 1.13d-4
      a2P31S1   = 1.76d+2
      at1s      = 1.76d+2+1.13d-4
      a2S11S1   = 5.13d+1
c      
      ome1s12p3 = 2.27d-2
      ome1s12s3 = 6.87d-2
      ome1s12s1 = 3.61d-2
      omep1 = 2.57d0
      omes1 = 0.55d0*omep1
c
      reche = de*dh*zion(atom)*pop(2,atom)*rec(4,atom)
      frs1  = reche*fusp(u)*(1.0d0-frs3(u))
      frp1  = reche*(1.0d0-fusp(u))*(1.0d0-frs3(u))
c
c    ***COLLISIONAL RATES FROM GROUND STATE 1S
c
      abup = de*dh*zion(atom)*pop(1,atom)
      ara  = abup*rka*f
c      
      xs1  = heien(2)/(rkb*t)
      xp1  = heien(3)/(rkb*t)
c      
      cos1 = (ara*omes1*dexp(-xs1))*fgaunt(1,1,xs1)
      cop1 = (ara*omep1*dexp(-xp1))*fgaunt(1,1,xp1)
c
c    ***RATES VIA TRIPLET STATES
c
      ra3  = reche*frs3(u)*0.33d0
      tns3 = ra3/(a2S31S1+(de*(qss(u)+qsp(u))))
      ra21 = a2S31S1*tns3
c
      ras1 = de*tns3*qss(u)
      rap1 = de*tns3*qsp(u)
      rap3 = de*tns3*qsp3(u)
c
c    ***TOTAL RATES AND APPROPRIATE COOLING
c
c	2S(T) > 1S(S)
c
      hein(1)   = (heien(1)*ra21)/pi4
      heibri(1) = hein(1)
c
c	2S(T) > 2S(S)
c 
      helos     = helos+(ras1*qren(2))+(cos1*heien(2))
      hein(2)   = (heien(2)*((ras1+frs1)+cos1))/pi4
c
c	2S(T) > 2P(S)
c
      helos     = (helos+(rap1*qren(3)))+(cop1*heien(3))
      hein(3)   = (heien(3)*((rap1+frp1)+cop1))/pi4
      heibri(2) = hein(3)
c
c	2S(T) > 2P(T)
c
      helos     = helos+(rap3*qren(4))
      hein(4)   = ((heien(4)*rap3)/pi4)+hel0(6,j108m)
      heibri(3) = hein(4)
c
      return 
      end
