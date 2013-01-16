cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO FIND INTERPOLATION VALUES USING NEIGHBORING POINTS TO
c	POINT*INL IN VECTOR BUFPHO .
c	BOUNDARIES OF BINS IN VECTOR EPHOT (EV) .
c	ASSUMES POWER LAW FOR INTERPOLATION : 
c		DLOG(Y)=CFLO+ALPHA*DLOG(NU)
c
c
c
      subroutine interpol(inl, cflo, alpha, bufpho)
c
      include 'cblocks.inc'
c
c
c
      double precision xx(5), yy(5), al(5), cc(5), bufpho(mxinfph)
      double precision cflo,alpha, he, hemi, yc, yg, yd
      integer*4 inl, n, n1, modd, modg, k, ini, ifi
c
c
c   ***FIND INTERPOL. VALUES IN THE INTERVAL EPHOT(INL),EPHOT(INL+1)
c	FOR THE MEAN INTENSITY IN BUFPHO(INL) ASSUMING A POWER LAW
c	FUNCTION OF ENERGY .
c
      he = bufpho(inl)
      al(5) = 0.0d0
      cc(5) = -1.d20
      if (he.le.0.0d0) goto 200
      k = 5
c
      cc(5) = dlog(he)
      if (inl.gt.2) then
      ini = 1
      modg = 1
      else
      ini = 3
      modg = 3
      end if
      if (inl.lt.(infph-2)) then
      ifi = 5
      modd = 4
      else
      ifi = 3
      modd = 2
      end if
c     
      if ((inl.le.1)) goto 200 
      if (((ifi-ini).ge.4) .and. 
     &((bufpho(inl-1)+bufpho(inl+1)).lt.(he/7.0d0))) goto 200
c
      do 60 n = ini, ifi
      n1 = (inl+n)-3
      hemi = bufpho(n1)
      if ((hemi.gt.(he/200.0d0)).and.(hemi.lt.(he*200.0d0))) then
      yy(n) = dlog(hemi)
      xx(n) = dlog(0.5d0*(ephot(n1)+ephot(n1+1)))
      else
      if (n.lt.3) modg = 3
      if (n.gt.3) modd = 2
      end if
   60 continue
c
      if (modg.gt.modd) goto 200
      do 120 n = modg, modd
      al(n) = (yy(n+1)-yy(n))/(xx(n+1)-xx(n))
      cc(n) = ((yy(n)*xx(n+1))-(yy(n+1)*xx(n)))/(xx(n+1)-xx(n))
  120 continue
c
      if ((modg .eq. 3).or.(modd .eq. 2)) goto 150
      yc = cc(5)
      yg = cc(1)+(al(1)*xx(3))
      yd = cc(4)+(al(4)*xx(3))
      if (dabs(yg-yc).gt.dabs(yd-yc)) then
      modg = 3
      else
      modd = 2
      end if
c
  150 if (modg .eq. 3) goto 180
      k = 2
      goto 200
  180 if (modd .eq. 2) goto 200
      k = 3
  200 continue
c
c	CONVERTS INTO FUNCTION OF FREQUENCY INSTEAD OF 2XENERGY (EV)
c
      alpha = al(k)
c
      cflo = cc(k)-(33.1191d0*alpha)
c
      return 

      end
