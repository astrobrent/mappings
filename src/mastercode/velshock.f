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
c
c*******COMPUTES SHOCK FRONT VELOCITY USING PRESHOCK DENSITY:DHPR
c	PRESHOCK AND POSTSHOCK TEMPERATURES : TEPR,TEPO
c	MAGNETIC FIELD : HMAG  , IONIC POPULATIONS : POP
c	AND FRACTIONAL IONISATION OF HYDROGEN : XHPR
c	NB.   ALL THESE QUANTITIES NEED TO BE DEFINED 
c
c	OUTPUT INVARIANTS AND VELOCITIES IN COMMON BLOCK /VELSHO/
c
c
      subroutine velshock(dhpr, xhpr, tepr, tepo, hmag)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision aa(10), a1(10)
      double precision mu, muf, dhpr, xhpr, tepr, tepo, hmag 
      double precision al,alp,b1,b2,b3,be,bet,co1
      double precision co2,co3,eps,ga,gam,gamf,humag,rho
      double precision root,tmin,u1,u2,v1
c
      integer*4 ncoef
      character rmod*4
c
c           Functions
c
      double precision feldens,fmua,frho
c
      if (xhpr.le.0.d0) xhpr = pop(2,1)
      if (pop(2,1) .ne. xhpr) pop(2,1) = xhpr
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
      if ((pop(1,1)+xhpr).ne.1.0d0) pop(1,1) = 1.0d0-xhpr
      depr = feldens(dhpr,pop)
      mu = fmua(depr,dhpr)
      rho = frho(depr,dhpr)
c  magnetic energy dens. u=B^2/8pi
      humag = (hmag*hmag)/25.13274d0 
      hmago = humag/rho
c
      alp = ((5.d0*rgas)*(tepo-tepr))/mu
      bet = (rgas*tepr)/mu
      gam = (rgas*tepo)/mu
c
      co1 = alp-(2.0d0*(gam-bet))
      al = alp/co1
      be = bet/co1
      ga = gam/co1
c
      co2 =-((2.0d0*al*be)+((ga-al) ** 2)-(be ** 2))
      co3 =-((al*be) ** 2)
      co2 = co2*co1
      co3 = co3*co1
      co1 = 1.0d0
      u1 = ((-co2)+dsqrt((co2*co2)
     &     -((4.d0*co1)*co3)))/(2.d0*co1)
      vshoc = dsqrt(u1)
 100  continue
c     
c     ***APPROXIMATE SHOCK VELOCITY COMPUTED : VSHOC
c     NOW TAKE MAGNETIC FIELD INTO ACCOUNT IN ITERATIVE FASHION
c     
      u1 = vshoc ** 2
      b1 = bet/u1
      b2 = hmag/u1
      b3 = 0.0d0
      aa(1) = 1.0d0
      aa(2) =-((5.0d0*((b1+b2)+1.0d0))/4.0d0)
      aa(3) = ((((5.0d0*b1)/4.0d0)+b2)+0.25d0)-(b3/2.0d0)
      aa(4) = b2
      ncoef = 4
      root = 1.0d0
c     
      call roots(aa, ncoef, root, rmod)
c     
      a1(1) = 1.0d0
      a1(2) = (aa(1)*root)+aa(2)
      a1(3) = (((aa(1)*root)*root)+(aa(2)*root))+aa(3)
      ncoef = 3
      root = 0.25d0
c     
      call roots(a1, ncoef, root, rmod)
c     
      vpo = root*vshoc
      u2 = vpo*vpo
      v1 = vshoc
      vshoc = 
     &  0.5d0*dsqrt(5.d0*(gam-bet)+(4.0d0*hmag)*((vshoc/vpo)-1.0d0)+u2)
     &  +(0.5d0*vshoc)
      eps = dabs((vshoc-v1)/vshoc)
      if (eps.gt.0.001) goto 100
c     
c     ***PHYSICAL CONDITIONS , POST-SHOCK
c     
      dhpo = (dhpr*vshoc)/vpo
      depo = feldens(dhpo,pop)
      mu = fmua(depo,dhpo)
      rho = frho(depo,dhpo)
      hmago = (hmago*vshoc)/vpo
      gam = (rgas*tepo)/mu
      ww = ((5.d0*gam)+(vpo*vpo)+(4.d0*hmago))/2.d0
      pres = rho*(gam+(vpo*vpo)+hmago)
      fm = rho*vpo
      qst = (2.d0*ww)-((5.d0*gam)+(4.d0*hmago)+(vpo*vpo))  !=0?!   
c     
c     ***COMPUTES APPROXIMATE FINAL STATE OF GAS
c     
      tmin = 100.0d0
      muf = fmua(0.0d0,1.0d0)
      gamf = (rgas*tmin)/muf
      vfin = ((tmin*mu)*vshoc)/(tepo*muf)
c     
      qfin = (2.0d0*ww)-(((5.0d0*gamf)+(vfin*vfin))+(((4.0d0*
     &     hmago)*vshoc)/vfin))
c     
      return 
      end
      
