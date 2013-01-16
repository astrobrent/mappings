cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c*******CALCULATES NET GAIN ARISING FROM HEATING DUE TO THE
c	ON THE SPOT APPROX. MINUS THE COOLING DUE TO RECOMB.
c	CALL SUBROUTINE SPOTAP
c
c     **REFERENCE : SEATON,M.J. (1959) MNRAS V.119,P.81
c	            HUMMER,D.G. AND SEATON,M.J. (1964)
c	               MNRAS,V.127,P.230
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine netgain(t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision t, de, dh, rcfa,rcfb,z
      double precision aa,ab,betr,de2,dep,dey,en,fhi
      double precision ph,ph2,rat,rc,rg,spo,u,yh
c
      integer*4 i,j
c
      rcfa(u,z) = 0.7714d0*((u/(z*z))**(-0.0390d0))
      rcfb(u,z) = 0.6412d0*((u/(z*z))**(-0.1146d0))
c      
      u = t/1.0d4
      if (u.lt.0.025d0) u = 0.025d0
c      if (u.gt.10.0d0) u = 10.0d0
      en = rkb*t
      betr = 0.0d0
      rg = 0.0d0
c
c    ***IF NECESSARY,ADS HEATING DUE TO THE ON THE SPOT APPROX.
c
      if (jspot .eq. 'NO') goto 113
c
      fhi = pop(1,1)
      call spotap(de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
      aa = (dh*zion(1))*pop(2,1)
      rc = aa*rec(3,1)
      aa = (dh*zion(2))*pop(2,2)
      rg = rg+((aa*de)*((dep*rec(4,2))+(dey*(rec(2,2)-rec(4,2)))))
      spo = aa*((yh*(rec(2,2)-rec(4,2)))+(ph*rec(4,2)))
      aa = (dh*zion(2))*pop(3,2)
      spo = spo+((aa*ph2)*rec(5,2))
      rat = spo/(rc+1.d-30)
      if (rat.le.1.0d0) goto 117
      spo = spo/rat
      de2 = de2/rat
  117 rg = rg+(((aa*de)*de2)*rec(5,2))
  113 continue
c
c    ***FINDS COOLING DUE TO RECOMBINATION
c
      do 101 i = 1, atypes
      do 100 j = 1, maxion(i)-1
      	z = dble(j)
      	ab = (dh*zion(i))*pop(j+1,i)
      	if ((i.le.2).and.(jspot .eq. 'YES')) then
      	betr = betr+(((de*ab)*rcfb(u,z))*rec(j+maxion(i),i))
      	else
      	betr = betr+(((de*ab)*rcfa(u,z))*rec(j+1,i))
      	end if
  100 continue
  101 continue
c
c
      rngain = (-(en*betr))+rg
c
      return 
      end

