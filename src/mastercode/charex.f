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
c*******CALCULATES CHARGE EXCHANGE RATES WITH HYDROGEN
c	RETURNS RATES IN CHARTE(3,X) AND CHARTE(4,X)
c	ELEMENT NUMBER IS IN CHARTE(1,X)
c	AND THE LOWEST IONIC SPECIE FOR THAT ELEMENT (1<=N<=5)
c	IN CHARTE(2,X)   (EX:N <> N+ CORRESPONDS TO SPECIE 1 )
c	NB. IF CHARCO(1,X) IS PRESET NEGATIVE , REACTION IS DISABLED
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine charex(t)
c
      include 'cblocks.inc'
c
      double precision t, tt, tl, tlim, alph, sigma,tl0,tl1,t4
      double precision x,a,ktev,lowT
      integer*4 i, j, jn
c
c  internal functions
c
      double precision fonct1, fonct2
      fonct1(x,a) = sigma*((x/1.0d4)**a)
      fonct2(x,a) = sigma*dexp(-(a/x))
c
      tlim = 6.0d4
      tt = t+1.d0
      t4 = tt*1.d-4
      tl = dmin1(tt,tlim)
c
      if (chargemode.eq.1) then
c
c     old calcs
c
      do 200 j = 1, nchxold
      do 202 jn = 1, 5
  202 charte(jn,j) = 0.0d0
c
      if (charco(1,j).le.0.0d0) goto 200
      charte(5,j) = charco(7,j)
      do 201 i = 1, 2
         charte(i,j)= charco(i,j)
         sigma = charco(2+i,j)
         alph = charco(4+i,j)
         if (alph.lt.20.0d0) charte(2+i,j) = fonct1(tl,alph)
         if (alph.ge.20.0d0) charte(2+i,j) = fonct2(tt,alph)
 201  continue
 200  continue
c
c     end old mappings exchange reactions
c
      endif
c
c
      if (chargemode.eq.0) then
c
c     new calcs, returns rates in chxr and chxi
c
c     first recombination rates
c
         do i = 1,nchxr
c
            tl0 = dlog10(chxrtemp(1,i))-0.5
            tl1 = dlog10(chxrtemp(2,i))+0.5
            tl = dlog10(t)
c
            if ((tl.gt.tl0).and.(tl.lt.tl1)) then
               chxr(i) = (chxrcos(1,i)*(T4)**chxrcos(2,i))*
     &                (1.d0 + chxrcos(3,i)*dexp(t4*chxrcos(4,i)))
c
c Allow for hard fields by entering constant adiabatic limit
c
			else if (tl.lt.tl0) then
			   lowT=10**(tl0-4)
			   chxr(i) = (chxrcos(1,i)*(lowT)**chxrcos(2,i))*
     &                (1.d0 + chxrcos(3,i)*dexp(lowT*chxrcos(4,i)))
c
            else
               chxr(i) = 0.0d0
            endif
c
         enddo
c
c     then ionising rates
c
         do i = 1,nchxi
c
c
            tl0 = dlog10(chxitemp(1,i))-0.5
            tl1 = dlog10(chxitemp(2,i))+0.5
            tl = dlog10(tt)
c
c
            if ((tl.gt.tl0).and.(tl.lt.tl1)) then
c
                ktev = rkb*tt/ev
c
c     special case for OI-HII reaction
c
               if (chxiat(i).eq.zmap(8)) then
                  chxi(i) = chxicos(1,i)*dexp(-1.d0*chxicos(4,i)/ktev)*
     &                      (1.d0-(0.93d0*dexp(-1.d0*chxicos(3,i)*t4)))
c
               else
c
c     the rest of the reactions
c
                  chxi(i) = (chxicos(1,i)*((t4)**chxicos(2,i)))*
     &                      dexp(-1.d0*chxicos(3,i)*t4)*
     &                      dexp(-1.d0*chxicos(4,i)/ktev)
c
               endif

            else
               chxi(i) = 0.0d0
            endif
c
         enddo
c
c
c
      endif
c
      return 
      end


