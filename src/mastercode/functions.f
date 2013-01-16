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
c*******TO DETERMINE CASE A,B  (CASEAB  0=CASE A , 1=CASE B)
c   SUBROUTINE  RESTRANS  NEEDS TO BE PREVIOUSLY CALLED
c   USES ESCAPE PROBABILITY OF LYMAN-GAMMA PHOTONS
c   for hydrogen or helium
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      double precision function casab(nz,line,series)
c
      include 'cblocks.inc'
c
      double precision trg,  x_sec, abfrac
      integer*4 line,series,nz
c
c     currently only line=3,series=1 considered
c
      trg = 0.4233d0
c
      x_sec = 0.d0
      if (nz.eq.1) then
          x_sec = hydlin(2,line,series)-1.0d0
      endif
      if (nz.eq.2) then
          x_sec = hellin(2,line,series)-1.0d0
      endif
c      
      if (x_sec.lt.0.0d0) x_sec = 0.0d0
      abfrac = trg ** x_sec
c
      casab = 1.0d0-abfrac
c
      return 
      end
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
c*******GIVES TOTAL NUMBER DENSITY INCLUDING ALL ELEMENTS  .
c   USES RELATIVE ABUNDANCES  :  ZION(I)
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function densnum(dh)
c
      include 'cblocks.inc'
c
      integer*4 i
      double precision dh
c
      densnum = 0.0d0
c
      do 100 i = 1, atypes
         densnum = densnum+zion(i)
  100 continue
c
      densnum = dh*densnum
c
      return 
      end
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
c*******GIVES TOTAL DENSITY IN atomic mass units / cc
c   USES ATOMIC MASSES  :  ATWEI   AND RELAT. ABUND.  :  ZION
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function denstot(dh)
c
      include 'cblocks.inc'
c
      integer*4 i
      double precision dh,d
c
      d = 0.0d0

      do 100 i = 1, atypes
 100     d = d+(zion(i)*atwei(i))
c
      denstot = dh*d
c
      return 
      end
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
c*******COMPUTES INTEGRAL FOR THE PHOTOIONISATION RATES
c   ARGUMENTS B1 AND B2 IN DOUBLE PRECIS.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function dfuint(bet, sa, b1, b2, alo)
c
      include 'const.inc'
c
      double precision da, db, dc, dd, dl, dk, b2, b1
      double precision bet, sa, alo, sa1,df
c
      sa1 = sa+1.0d0
      dfuint = 0.d0
      df = 0.0d0
c
      da = 0.d0
      db = 0.d0
      dc = 0.d0
      dd = 0.d0
      dl = 0.d0
      dk = 0.d0
c
      if (sa1.lt.12.0d0) then
         if (sa .eq. 0.0d0) then
            da = bet*(((b2 ** sa1)/sa1)-((b1 ** sa1)/sa1))
            db = (1.0d0-bet)*(dlog(b2)-dlog(b1))
         else if (sa1 .eq. 0.0d0) then
            da = bet*(dlog(b2)-dlog(b1))
            db = (1.0d0-bet)*(((b2 ** sa)/sa)-((b1 ** sa)/sa))
         else
            da = bet*(((b2 ** sa1)/sa1)-((b1 ** sa1)/sa1))
            db = (1.0d0-bet)*(((b2 ** sa)/sa)-((b1 ** sa)/sa))
         end if
c
         dl = da+db
c
         if (dl.le.0.d0) then
            df = 0.0d0
         else
            if (dabs(alo+dlog(dl)).ge.150.d0) then
c               write (*,*) 'dfuint error:alo+dlog(dl):',(alo+dlog(dl))
c               write (*,*) 'dfuint error:alo,dlog(dl):',alo,dlog(dl)
            else
               df = dexp(alo+dlog(dl))
            endif
         end if
      else
         da = (alo+(sa1*dlog(b2)))-dlog(dabs((sa1)))
         db = (alo+(sa1*dlog(b1)))-dlog(dabs((sa1)))
         if ((dabs(da).lt.loghuge).and.(dabs(db).lt.loghuge)) then
            dl = (bet*(dexp(da)-dexp(db)))*dsign(1.0d0,sa1)
         endif
         dc = (alo+(sa*dlog(b2)))-dlog(dabs((sa)))
         dd = (alo+(sa*dlog(b1)))-dlog(dabs((sa)))
         if ((dabs(dc).lt.loghuge).and.(dabs(dd).lt.loghuge)) then
            dk = ((1.0d0-bet)*(dexp(dc)-dexp(dd)))*dsign(1.0d0,sa)
         endif
c
         df = dmax1(0.d0,dl+dk)
c
      end if
c
      dfuint = df
c
      return 
      end
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
c     This subroutine will evaluate the exponential integrals
c     used by Arnaud and Rothenflug 1985 to evaluate collisional
c     ionisation rates.
c     
c     ref: Arnaud, M., Rothenflug, R. 1985 A. A. Suppl. Ser. 62:425
c     particularly pages 435-436.
c
c     NOTE: there are some algebraic errors in the expressions published
c     and the function for f2 is only valid for large X.  I have used
c     my own series expansion for small x < 3 and AR's one for large X.
c
c     RSS 7/90
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function farint(sub,x)
c
      implicit none

      double precision x
      integer*4 sub
c
c     sub is an index to the two different integrals
c
      double precision xp,p,q,xinv
      double precision fint, f
      double precision a(0:5)
      double precision b(1:4)
      double precision c(1:4)
      double precision pj(0:14), qj(0:14)
      integer*4 j
c
c
c     data for coefficients
c
      data a / -0.57721566490d0,0.99999193d0,-0.24991055d0,
     &0.05519968d0,-0.00976004d0,0.00107857d0/
      data b / 8.5733287401d0,18.059016973d0,8.6347608925d0,
     &0.2677737343d0/
      data c / 9.5733223454d0,25.6329561486d0,21.0996530827d0,
     &3.9584969228d0/

      data pj / 1.d0,2.1658d2,2.0336d4,1.0911d6,3.7114d7,8.3963d8,
     &1.2889d10,1.3449d11,9.4002d11,4.2571d12,1.1743d13,1.7549d13,
     &1.0806d13,4.9776d11,0.0d0 /
c
      data qj / 1.0d0,2.1958d2,2.0984d4,1.1517d6,4.0349d7,9.4900d8,
     &1.5345d10,1.7182d11,1.3249d12,6.9071d12,2.3531d13,4.9432d13,
     &5.7760d13,3.0225d13,3.3641d12 /
c
c     begin calculations
c
      fint = 0.d0
      if (sub.eq.1) then
c
c     evaluate f1(x) = e^x integral(1 to inf)(e^-xt)/t dt
c
c     Eqn 5.1.53 pg231
c
         if (x.lt.1.d0) then
            f = a(0)+x*(a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*a(5)))))-dlog(x)
            fint = dexp(x)*f
         else
c
c     Eqn 5.1.56 pg 231
c
            fint = ( ( b(4)+x*(b(3)+x*(b(2)+x*(b(1)+x))))/
     &             (c(4)+x*(c(3)+x*(c(2)+x*(c(1)+x)))))/x
         endif
c
      else
c
c     sub = 2
c
c     evaluate f2(x) = e^x integral(1 to inf) (e^-xt)ln(t)/t dt (Note: *not* E_2 RSS2006)
c
            p = 0.d0
            q = 0.d0
            xp  = 1.d0
            xinv = 1.d0/x
            do j = 0,14
               p = p + xp * pj(j)
               q = q + xp * qj(j)
               xp = xp*xinv
            enddo
            fint = xinv*xinv*p/q
c
      end if
c
      farint = fint
c
      return
      end

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
c     function returns collisional autoionisation 
c     rate for atom,ion. ref Arnaud&Rothenflug 1985 A.A.Supp.Ser.60:425
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fautoi(t,atom,ion)
c
      include 'cblocks.inc'
c
      double precision t,f1,g,y,a,b,zr,cea
      double precision term1,term2,ktev,iea
      integer*4 atom,ion,isos,z
c
c           Functions
c
      double precision farint
c
      fautoi = 0.d0
      cea = 0.d0
      z = mapz(atom)
      zr = dble(z)
      isos = z - ion + 1
      ktev = rkb*t/ev
c
c
      if (isos.eq.3) then
c
c     Lithium iso series
c
         term1 = (zr-0.835d0)*(zr-0.835d0)
         term2 = (zr-1.62d0)*(zr-1.62d0)
         iea = ryd*term1-0.25d0*term2
         b = 1.d0/(1.d0+2.d-4*zr*zr*zr)
c
c     an effective z due to screening?
c
         zr = zr-0.43d0
c
c     yet another gaunt factor
c
         y = iea/ktev
c
         cea = 0.d0
         if (y.le.150.d0) then 
c
         f1 = farint(1,y)
         g = 2.22d0*f1+0.67d0*(1.d0-y*f1)+0.49d0*y*f1+1.2d0*
     &        y*(1.d0-y*f1)
c
c     see ref for explanation of 1.2 factor (Appendix A)
c
         cea = 1.2d0*(1.60d-7*b*dexp(-y)*g/(zr*zr*dsqrt(ktev)))
c
         endif
c
c     correct for some particular species
c
         if (z.eq.6) cea=cea*0.6d0
         if (z.eq.7) cea=cea*0.8d0
         if (z.eq.8) cea=cea*1.25d0
c
      endif
c
c
      if (isos.eq.11) then
c
c     Sodium iso series
c
         if (z.le.16) then
c
            iea = 26.d0*(zr-10.d0)
            y = iea/ktev
            cea = 0.d0
            if (y.le.150.d0) then 
c
            f1 = farint(1,y)
c
            a = 2.8d-17*(zr-10.d0)**(-0.7d0)
c
            cea = 6.69d7*a*iea*dexp(-y)*(1.d0-y*f1)/dsqrt(ktev)
c
            if (cea.lt.0.d0) cea = 0.d0
            endif
c
         else
c     (16 <z <28)
c
            iea = 11.d0*(zr-10.d0)**1.5d0
            y = iea/ktev
            cea = 0.d0
            if (y.le.150.d0) then 
c
            f1 = farint(1,y)
c
            a = 1.3d-14*(zr-10.d0)**(-3.73d0)
c
            cea = 6.69d7*a*iea*dexp(-y)*(1.d0-(y-y*y+y*y*y*f1)/2.d0)
     &           /dsqrt(ktev)
c
            endif
c
         endif
      endif
c
c
c
      if ((z.ge.18).and.(isos.gt.11).and.(isos.lt.17)) then
c
c     magnesium through sulphur series for heavy elements
c
         if (isos.eq.12) iea = 10.3d0*(zr-10.d0)**1.52d0
         if (isos.eq.13) iea = 18.0d0*(zr-11.d0)**1.33d0
         if (isos.eq.14) iea = 18.4d0*(zr-12.d0)**1.36d0
         if (isos.eq.15) iea = 23.7d0*(zr-13.d0)**1.29d0
         if (isos.eq.16) iea = 40.1d0*(zr-14.d0)**1.10d0
c
         y = iea/ktev
         cea = 0.d0
         if (y.le.150.d0) then 
c
         f1 = farint(1,y)
c
         a = (4.d-13/(zr*zr))/iea
c
c
         cea = 6.69d7*a*iea*dexp(-y)*(1.d0-(y-y*y+y*y*y*f1)/2.d0)
     &        /dsqrt(ktev)
c
         endif
c
      endif
c
c
c
      if ((z.eq.20).and.(ion.le.2)) then
c
c     calcium I and II
c
         iea = 25.d0+4.d0*dble(ion-1)
c
         a = 6.0d-17+3.8d-17*dble(ion-1)
         b = 1.12d0
c
         y = iea/ktev
c
         cea = 0.d0
         if (y.le.150.d0) then 
c
         f1 = farint(1,y)
c     
         cea = 6.69d7*a*iea*dexp(-y)*(1.d0+b*f1)/dsqrt(ktev)
c
         endif
c
      endif
c
      if ((z.eq.26).and.((ion.eq.4).or.(ion.eq.5))) then
c
c     Fe IV and V
c
         iea = 60.d0+13.d0*dble(ion-4)
c
         a = 1.8d-17 - 1.3d-17*dble(ion-4)
         b = 1.00d0
c
         y = iea/ktev
c
         cea = 0.d0
         if (y.le.150.d0) then 
c
         f1 = farint(1,y)
c     
         cea = 6.69d7*a*iea*dexp(-y)*(1.d0+b*f1)/dsqrt(ktev)
c
         endif
c
      endif
c
c     pass back result
c
c      if (z.eq.17) then
c         write(*,*) t,' ',elem(atom),rom(ion),cea
c      endif
c
      fautoi = cea
c
      return
c
      end
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
c*******TO DERIVE AVERAGE POSITIVE CHARGE PER NUCLEON
c   LINEAR:NOR=1   ;   RMS:NOR=2
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function favcha(popul, nor)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision fne,wei,weito
      double precision popul(mxion, mxelem)
      integer*4 i,j,nor,norder
c
      norder = max0(1,min0(4,nor))
      fne = 0.0d0
      weito = 0.0d0
c
      do i = 1, atypes
         wei = zion(i)
         weito = weito+wei
         do j = 2, maxion(i)
            fne = fne+((wei*((j-1)**norder))*popul(j,i))
         enddo
      enddo
c
      favcha = (fne/weito)**(1.0d0/norder)
c
      return 
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     function to evaluate the modified bessel function of 
c     the first kind Io(x)
c
c     This is required for the exact free free gaunt
c     factor calculation
c
c     ref Mewe et al  1986 A.A.Suppl.Ser. 65:511
c     Abramowitz & Stegun Handbook of mathematical tables
c      (eds) 
c
c
c     Io(x) = 1/pi * integral(0-inf) cosh(x * cos t) dt
c
c
c     errors:
c       x <=3.75 : error < 2e-7
c       x >3.75 : error < 2e-7
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fbessi(x)
c
      include 'const.inc'
c
      double precision x,i0,x3,x32
      double precision a1,a2,a3,a4,a5,a6
      double precision b1,b2,b3,b4,b5,b6,b7,b8,b9
c
c     assign polynomial constants
c
      a1 = 3.5156229d0
      a2 = 3.0899424d0
      a3 = 1.2067492d0
      a4 = 0.2659732d0
      a5 = 0.0360768d0
      a6 = 0.0045813d0
c
      b1 = 0.39894228d0
      b2 = 0.01328592d0
      b3 = 0.00225319d0
      b4 = -0.00157565d0
      b5 = 0.00916281d0
      b6 = -0.02057706d0
      b7 = 0.02635537d0
      b8 = -0.01647633d0
      b9 = 0.00392377d0
c
      fbessi = 0.0d0 
      i0 = 0.0d0
c
c     this is only valid for  x > -3.75
c
      if (x.ge.-3.75d0) then
         if (x.le.3.75d0) then
            x3 = x/3.75d0
            x32 = x3*x3
            i0 =1.d0+x32*(a1+x32*(a2+x32*(a3+x32*(a4+x32*(a5+x32*a6)))))
         else
            x3 = 3.75d0/x
            i0 = b1+x3*(b2+x3*(b3+x3*(b4+x3*(b5+x3*(b6+x3*(b7+x3*
     &           (b8+x3*b9)))))))
            i0 = i0/(dsqrt(x)*dexp(-x))
         endif
      else
         write (*,*) 'Error in fbessi, x argument too small'
         write (*,*) '!!!!!!!!!!!!!!!!!!!!!!x = ',x
         stop
      endif
c     
      fbessi = i0
c     
      return
c     
      end
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
c     function to evaluate the modified bessel function of 
c     the second kind exp(x)*Ko(x)
c
c     This is required for the exact free free gaunt
c     factor calculation
c
c     ref Mewe et al  1986 A.A.Suppl.Ser. 65:511
c     Abramowitz & Stegun Handbook of mathematical tables
c      (eds) 
c
c
c     Ko(x) = integral(0 to infinity) cos(xt)/sqrt(t^2+1) dt
c
c     {x > 0}
c
c
c     errors:
c       x <2 : error < 1e-8
c       x >2 : error < 2e-7
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fbessk(x)
c
      include 'const.inc'
c
      double precision fbessi
      double precision x,k0,x2,x22
      double precision a1,a2,a3,a4,a5,a6
      double precision b1,b2,b3,b4,b5,b6,b7
c
c     assign polynomial constants
c
      a1 = 0.42278420d0
      a2 = 0.23069756d0
      a3 = 0.03488590d0
      a4 = 0.00262698d0
      a5 = 0.00010750d0
      a6 = 0.00000740d0
c
      b1 = 1.25331414d0
      b2 = -0.07832358d0
      b3 = 0.02189568d0
      b4 = -0.01062446d0
      b5 = 0.00587872d0
      b6 = -0.00251540d0
      b7 = 0.00053208d0
c
      fbessk = 0.0d0 
c
c     this is only valid for positive x
c
      if (x.ge.0.0d0) then
         if (x.le.2.0d0) then
            x2 = x/2.0d0
            x22 = x2*x2
            k0 = -dlog(x2)*fbessi(x) - euler
           k0 = k0+(x22*(a1+x22*(a2+x22*(a3+x22*(a4+x22*(a5+x22*a6))))))
           k0 = k0*dexp(x)
         else
            x2 = 2.0d0/x
            k0 = b1+x2*(b2+x2*(b3+x2*(b4+x2*(b5+x2*(b6+x2*b7)))))
            k0 = k0/dsqrt(x)
         endif
      else
         write (*,*) 'Error in fbessk, negative x argument'
         write (*,*) '!!!!!!!!!!!!!!!!!!!!!!x = ',x
         stop
      endif
c     
      fbessk = k0
c     
      return
c     
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO DERIVE AVERAGE collsionial ionisation equilibrium 
c       TIME SCALE
c   ASSUMING CONSTANT ELECTRONIC DENSITY : DE
c       Weighted by population fractions and abundances
c
c
      double precision function fcietim(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision t, de, dh,rcol,totn,wei
      integer*4 i,j
c
      rcol = 0.d0
      totn = 0.d0
      fcietim = 0.d0
c
      do i = 1 ,  atypes
         do j = 1 , maxion(i)-1
            totn = totn+pop(j,i)*zion(i)
         enddo
      enddo
c
      do 100 i = 1, atypes
c
      do 201 j = 1 , maxion(i)-1
         wei = pop(j,i)*zion(i)
         rcol = rcol+(wei*(col(j,i)*de+epsilon))
 201  continue
c     
 100  continue
c
      fcietim = totn/(rcol+epsilon)
c
      return 
      end
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
c       Calculates Weighted Collisional Ionisation Timescales
c
c       Weighted by population fractions and abundances
c       Assumes that the rates are all up to date. call
c       Allrates otherwise.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fcolltim(de)
c
      include 'cblocks.inc'
c
      double precision de,rcol,totn,wei
      integer*4 i,j
c
      rcol = 0.d0
      totn = 0.d0
      fcolltim = 0.d0
c
      do i = 1 ,  atypes
            totn = totn+zion(i)
      enddo
c
      do i = 1, atypes
c
      do j = 1 , maxion(i)-1
         wei = pop(j,i)*zion(i)
         rcol = rcol+(wei*col(j,i))
      enddo
c     
      enddo
c
      fcolltim = totn/(de*rcol+epsilon)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******DETERMINE THE CONVERGENCE CRITERION USING THE
c   SINC FUNCTION TO THE POWER ZETA
c   DLOS  MUST  SATISFY  :   -1<=DLOS <=1
c   OUTPUT  :  0<=FCRIT<=1.0
c
      double precision function fcrit(dloss, zeta)
c
      include 'const.inc'
c
      double precision dloss, zeta
c
      fcrit = ((dabs(dsin((pi*dloss)+1.d-10))
     &        /dabs((pi*dloss)+1.d-10))+1.d-10) ** zeta
c
      return 
      end

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
c*******COMPUTES GEOMETRICAL DILUTION FACTOR
c   USING PHOTON SOURCE RADIUS : RSOU  AND RADIUS AT WHICH
c   IT IS CALCULATED : RAD
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fdilu(rsou, rad)
c
      include 'cblocks.inc'
c
      double precision rio,gdilf,phi,cphi
      double precision rsou, rad, rsou2, rad2
c
c     default plane || factor (jgeo=P)
c
      fdilu = 0.5d0
      gdilf = 1.0d0
c
      if (jgeo.eq.'S') then
c
c     spherical geometry
c
         if (rad.le.rsou) then
            gdilf = 1.d0
         else
            rio = (rsou)/rad
            if (rio.gt.4.0d-4) then
c
c     double precision trig functions work
c     Still need to check if this is the correct function
c
               rsou2 = rsou*rsou
               rad2 = rad*rad
               phi = dasin(rio)
               cphi = dcos(phi)
               gdilf = (rsou2-rad2*(1.0d0-cphi))/(rad2*cphi)
            else
               gdilf = (0.5d0*rio*rio)
            endif
         endif
      end if
c
      if (jgeo.eq.'F') then
c
c     finite plane parrallel
c
         if (rad.le.0.d0) then
            gdilf = 1.d0
         else
            rio = (rsou)/rad
            if (rio.gt.4.0d-4) then
c
c     double precision trig functions work
c
               phi = datan(rio)
               cphi = dcos(phi)
               gdilf = (1.0d0-cphi)
            else
               gdilf = (0.5d0*rio*rio)
            endif
         endif
      end if
c
      fdilu = 0.5d0*gdilf
c
      return 
c
      end

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
c     This is the new function to replace restrans. It will calculate
c     the radiative line transfer distance multiplier for the new
c     XLINDAT short wavelength resonance lines <2000 Angstoms.
c
c     Old resonance lines with Lambda longer than this will be simplly
c     added to the diffuse field like the inter-combination lines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c*******COMPUTES THE RESONANCE LINE TRANSFER USING ESCAPE
c   PROBABILITY FORMULATION
c   DV IS THE CHANGE IN VELOCITY OF THE GAS FROM R TO R+DR
c   RETURNS EFFECTIVE STEP LENGTH MULTIPLIER FOR THE LINE.
c     **REF. : CAPRIOTTI,E.R.,1965 APJ.142,1101 (EQ* 88A,B)
c          (SIGN ERRORS IN 88A HAVE BEEN CORRECTED)
c
      double precision function fdismul(t, dh, dr, dv, atom, ion, ejk
     &                                 , fab)
c
c           ejk in ergs as usual
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision abu,con,conte,cte,dreff,dveff,rkte,tau
      double precision tauther,vra,vther,x,z
      double precision t, dh, dr, dv, ejk, fab , telc
      double precision e7a,e7b,argu,erf
      double precision a(5),f(3),e
      integer*4 atom,ion
c
c
c
      e7a(x) = (1.0d0+((0.71d0*x)*dlog((2.0d0*x)+1.d-37)))
     &        +(((((((((a(1)*x)+a(2))*x)+a(3))*x)+a(4))*x)+a(5))*x)
      e7b(x) = ((0.14d0+(0.25d0/dsqrt(dlog(2.0d0*x))))
     &         +dsqrt(dlog(2.0d0*x)))/(3.5449077d0*x)
      argu(z) = 1.0d0/(1.0d0+(0.47047d0*z))
      erf(z) = 1.0d0-(((f(1)+((f(2)+f(3)*argu(z))
     &         *argu(z)))*argu(z))*dexp(-(z*z)))
c
c     constants for some functions
c
      data a/0.00048d0,-0.019876d0,0.08333d0,-0.3849d0,-0.83d0/
      data f/0.3480242d0,-0.0958798d0,0.7478556d0/
c
      dveff = dabs(dv)
      dreff = dabs(dr)
      telc = dmax1(t,10.0d0)
      rkte = rkb*telc
      cte = rkte/mp
      con = 3.0d-18
c
c      fab = xfef(line)/xbr(line)
c      ion = xion(line)
c      atom = xrat(line)
c      ejk = xejk(line)*ev
c
      if ((ejk/rkte).lt.loghuge) then
         conte = con*(1.0d0-dexp(-ejk/rkte))
      else
         conte = con
      endif
c     
c     
      vther = dsqrt((2.0d0*cte)/atwei(atom))
      abu = dh*zion(atom)*pop(ion,atom)*conte*fab*dreff/ejk
c
      tauther = abu/vther
      vra = dveff/vther
      if (vra.le.0.01d0) tau = tauther
      if (vra.gt.0.01d0) tau = (tauther*0.886227d0)*(erf(vra)/vra)
      if (tau.lt.1.5d0) e = e7a(tau)
      if (tau.ge.1.5d0) e = e7b(tau)
c     
      fdismul = dmax1(1.0d0,1.0/(e+1.0d-8))
c     
c     
      return
c 
      end
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
c	Approximate Bowen Flourescence losses to HeII303  tau_HeLa>>>Tau_c
c   Kallman & MacCray 1980
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      double precision function fbowen(t, dv)
c
c           ejk in ergs as usual
c
c     Gives the fration of He II 303A that is converted to Bowen lines
c
      include 'cblocks.inc'
c
      double precision t, dv
c
      double precision aba, aboiii, abheii, abhei
      double precision dveff,vra,vther
      double precision telc
      double precision rkte
      double precision cte
      double precision con,fb,fb2,fb3
c      
      double precision taucr, taula
c      
      integer*4 atom
c
c     constants for some functions
c
c
      fb = 0.d0
c
      dveff = dabs(dv)
      telc  = dmax1(t,10.0d0)
      rkte  = rkb*t
      cte   = rkte/mp
      con   = 2.0d-3
c     
      
      atom   = zmap(8)
      if (atom.gt.0) then
      vther  = dsqrt((2.0d0*cte)/atwei(atom))
      aboiii = zion(atom)*pop(3,atom)
      atom   = zmap(2)
      abheii = zion(atom)*pop(2,atom)
c
      taula = 2.0d6*abheii
c
      abhei  = zion(atom)*pop(1,atom)
      atom   = zmap(1)
      aba    = zion(atom)*pop(1,atom)+9.3d0*abhei
c      
      taucr  = 80.d0*abheii/aboiii
c
      fb2 = dmax1(0.0d0,1.0d0/(1.d0+con*(aba/aboiii)))
      fb3 = dmax1(0.0d0,1.0d0/(1.d0+(taucr/taula)))
c
c	degrade efficiency for velocity gradients.
c     RSS2001
c
      vra = dveff/vther
c      fb  = (fb2*fb3)/(1.d0+(vra**2))
      fb  = (fb2*fb3)*exp(-vra)
c     
c      write(*,*) "Bowen Yield", fb, vra
c
      endif
      
      fbowen = fb
c     
      return
c 
      end
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
c*******DERIVES THE ELECTRONIC DENSITY USING IONIC POPULATIONS
c   IN MATRIX POPUL(6,11)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function feldens(dh, popul)
c
      include 'cblocks.inc'
c
      double precision dh, popul(mxion, mxelem), deelec
c
c
      integer*4 i, j
c
      deelec = 0.0d0
      do 101 i = 1, atypes
         do 100 j = 1, maxion(i)-1
            deelec = deelec+((j*zion(i))*popul(j+1,i))
 100     continue
 101  continue
c
      feldens = dh*deelec
c
      return 
      end
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
c*******TO COMPUTE ELECTRON TO NEUTRAL DENSITY RATIO (NEIO=1)
c   OR ELECTRON TO ION DENSITY RATIO
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function felneur(popul, neio)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision ael,ane
      double precision popul(mxion, mxelem)
      integer*4 i,j,mm,neio
c
      mm = min0(2,max0(1,neio))
      ane = 0.0d0
      ael = 0.0d0
      do 101 i = 1, atypes
         do 100 j = mm, maxion(i)
            ane = ane+(zion(i)*popul(j,i))
            ael = ael+((zion(i)*popul(j,i))*(j-1.0d0))
 100     continue
 101  continue
      if (mm .eq. 1) then
         felneur = ael/ane
      else if ((ael+ane).gt.0.0d0) then
         felneur = ael/ane
      else
         felneur = 1.d0
      end if
c
      return 
      end

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
c   interpolates gfb from table using 2D cubics in log gfb space.
c   based on Numerical Recipes II
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fgfblog(tz,en,il,jl)
c
c     Assymetric order interpolation, assumes linear in en, cubic in tz 
c
      include 'cblocks.inc'
c     
      INTEGER*4 m,n,MMAX,NMAX
      REAL*8 en,tz,y,dy
      REAL*8 x1,x2
      PARAMETER (NMAX=20,MMAX=20)
c
      INTEGER*4 i,j,il,jl,ktz,ken,ken2,n2,m2,ngfbt
CU    USES polint
      REAL*8 ymtmp(MMAX),yntmp(MMAX)
      REAL*8 x1tmp(MMAX),x2tmp(MMAX)
c
      fgfblog = 0.d0
c
      m = 2
      m2 = m+2
c
      ktz  = min(max(il-(m2-1)/2,1),ngfbtz+1-m2)
c
      ngfbt = (ngfbe/2)+1
      ken  = min(max(jl-(m-1)/2,1),ngfbt+1-m)
      ken2 = (ken*2)-1
c
      n = m-1
      n2 = m2-1
c      
      x1 =en
      x2 =tz
c
      do i=ktz,ktz+n2
        x1tmp(i-ktz+1) = gfbtz(i)
      enddo
c
      do j=ken,ken+n
        x2tmp(j-ken+1) = gfbe(j)
      enddo
c
c 4 Linear interps for En between edges
c
      do i=ktz,ktz+n2
        do j=ken2,ken2+n
          yntmp(j-ken2+1) = gfb(i,j)
        enddo
        call polint(x2tmp,yntmp,m,x1,ymtmp(i-ktz+1),dy)
      enddo
c
c 1 Cubic interp for tz
c
      call polint(x1tmp,ymtmp,m2,x2,y,dy)
c
      fgfblog = 10.d0**y
c
      return
      END
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
c   interpolates gff from table using 2D cubics in lin gff space.
c   based on Numerical Recipes II
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fgfflin(m,g2,u,il,jl)
c
      include 'cblocks.inc'
c     
      INTEGER*4 m,n,MMAX,NMAX
      REAL*8 g2,u,y,dy
      REAL*8 x1,x2
      PARAMETER (NMAX=20,MMAX=20)
c
      INTEGER*4 i,j,il,jl,kg2,ku
CU    USES polint
      REAL*8 ymtmp(MMAX),yntmp(MMAX)
      REAL*8 x1tmp(MMAX),x2tmp(MMAX)
c
      fgfflin = 1.0
c
      kg2 = min(max(il-(m-1)/2,1),ngffg2+1-m)
      ku = min(max(jl-(m-1)/2,1),ngffu+1-m)
c
      n = m-1
c      
      x1 =g2
      x2 =u
c
      do i=kg2,kg2+n
        x1tmp(i-kg2+1) = gffg2(i)
      enddo
      do j=ku,ku+n
        x2tmp(j-ku+1) = gffu(j)
      enddo
c
      do i=kg2,kg2+n
        do j=ku,ku+n
          yntmp(j-ku+1) = 10.d0**gff(i,j)
        enddo
        call polint(x2tmp,yntmp,m,x2,ymtmp(i-kg2+1),dy)
      enddo
c
      call polint(x1tmp,ymtmp,m,x1,y,dy)
c
      fgfflin = y
c
      return
      END
c
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
c   interpolates gff from table using 2D cubics in log gff space.
c   based on Numerical Recipes II
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fgfflog(m,g2,u,il,jl)
c
      include 'cblocks.inc'
c     
      INTEGER*4 m,n,MMAX,NMAX
      REAL*8 g2,u,y,dy
      REAL*8 x1,x2
      PARAMETER (NMAX=20,MMAX=20)
c
      INTEGER*4 i,j,il,jl,kg2,ku
CU    USES polint
      REAL*8 ymtmp(MMAX),yntmp(MMAX)
      REAL*8 x1tmp(MMAX),x2tmp(MMAX)
c
      fgfflog = 1.0
c
      kg2 = min(max(il-(m-1)/2,1),ngffg2+1-m)
      ku = min(max(jl-(m-1)/2,1),ngffu+1-m)
c
      n = m-1
c      
      x1 =g2
      x2 =u
c
      do i=kg2,kg2+n
        x1tmp(i-kg2+1) = gffg2(i)
      enddo
      do j=ku,ku+n
        x2tmp(j-ku+1) = gffu(j)
      enddo
c
      do i=kg2,kg2+n
        do j=ku,ku+n
          yntmp(j-ku+1) = gff(i,j)
        enddo
        call polint(x2tmp,yntmp,m,x2,ymtmp(i-kg2+1),dy)
      enddo
c
      call polint(x1tmp,ymtmp,m,x1,y,dy)
c
      y = 10.d0**y
      fgfflog = y
c
      return
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  (C) Copr. 1986-92 Numerical Recipes Software $2'%1&&")#+`[.

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then 
           write(*,*) 'failure in polint'
           stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $2'%1&&")#+`[.
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
c*******TO COMPUTE THE GAUNT FACTOR IN COLLISIONAL
c   EXCITATION FOR ARBITRARY XSI (EN.GAP/KT) AND IONIC SPECIE
c   NDEL IS THE  CHANGE IN THE PRINCIPAL QUANTUM NUMBER
c   
c     **REFERENCES : BELY,O. 1966 PROC. PHYS. SOC. 88,587
c            TARTER,C.B. 1969 AP J SUPPL. 18,1
c                VAN REGEMORTER,H 1962 AP J 136,P 906
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fgaunt(ion, ndel, xsi)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision u,y
      double precision xsi
      integer*4 ion, ndel
c
c           Functions
c
      double precision a,aa,b,bb,c,cc,fue1
c
      a(y) = (((5.232d0*y)-6.652d0)*y)+2.878d0
      b(y) = (((2.79d0*y)+1.764d0)*y)+0.898d0
      c(y) = (((0.7917d0*y)+0.2308d0)*y)+0.10741d0
      aa(y) = (((((- (1.9347d0*y))+3.0267d0)*y)-1.5369d0)*y)+0.3122d0
      bb(y) = (((10.574d0*y)+68.378d0)*y)+3.496d0
      cc(y) = (((- (0.2359d0*y))-0.07154d0)*y)+0.2667d0
c
      fgaunt = 0.2d0
c
c    ***IONIC SPECIES
c
      if (ion.le.1) goto 100
      y = 1.0d0/dble(ion)
c
c     **NO CHANGE IN PRINCIPAL QUANTUM NUMBER
c
      if (ndel.gt.0) goto 10
      u = xsi/c(y)
      fgaunt = (0.413497d0*a(y))*(1.d0+((b(y)*fue1(u,0.d0))/2.3026d0))
      goto 200
   10 continue
c
c     **WITH CHANGE IN PRINCIPAL QUANTUM NUMBER
c
      u = xsi/cc(y)
      fgaunt = (0.413497d0*aa(y))*(1.d0+((bb(y)*fue1(u,0.d0))/2.3026d0))
      goto 200
  100 continue
c
c    ***NEUTRAL SPECIES (ION=1)
c
      if (xsi.lt.0.5d0) then
      fgaunt = 0.276d0*fue1(xsi,xsi)
      else
      fgaunt = (0.066d0/dsqrt(xsi))+(0.033d0/xsi)
      end if
c
      if (ndel.le.0) fgaunt = 2.5d0*fgaunt
c
  200 continue
      return 
      end
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
c	Simple heavyside function
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      double precision function fheavyside(x)
c
      include 'cblocks.inc'
c
      double precision x,  fh
c
      fh = 1.d0
      if (x.lt.0.d0) fh = 0.d0
      fheavyside = fh
c
      return 
      end
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
c*******TO COMPUTE AVERAGE CHARGE OF IONS
c   LINEAR:IORD=1   ;   RMS:IORD=2
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fioncha(popul, iord)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision weito
      double precision popul(mxion, mxelem)
      integer*4 i,iek,j,iord
c
      iek = max(1,iord)
      weito = 0.0d0
c
      fioncha = 0.0d0
      do 101 i = 1, atypes
         do 100 j = 2, maxion(i)
            weito = weito+(zion(i)*popul(j,i))
            fioncha = fioncha+((zion(i)*popul(j,i))*dble((j-1) ** iek))
 100     continue
 101  continue
      if (weito.gt.0.0d0) then
         fioncha = (fioncha/weito) ** (1.d0/dble(iek))
      else
         fioncha = 1.0d0
      end if
c
      return 
      end
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
c     This function calculates the radiative recombination rates
c     for hydrogenic ions.  The Seaton fit to the Asymptotic 
c     expansion for the Kramer-Gaunt term is used at low t (<Z^2x10^6K) 
c
c     ref: Arnaud and Rothenflug 1985 A.A.Suppl.Ser. 60:425
c          Seaton 1959 MNRAS 119
c          Burgess 1958 MNRAS 118
c
c     
c
c     RSS 8/90
c
c     The rates for all except the ground term for H and He are left in
c     rec(ion+maxion(atom),atom) for use in the On-The-Spot Approximation.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New calcualtions, use Seaton for log(t)<5, cubic splines
c     fits to explicit sums of hydrogenic wavefunction overlap intergrals
c     convolved with maxwellian electron distribution.
c
c       Total Radiative Hydrogenic Recombination. 10 < Te/Z^2 <10^9 K Z^2
c
c
c       Cubic Spline Fit Coefficients: x,a,b,c,d [y=a+bx+cx^2+dx^3]
c                                      in the interval x - x+1
c
c       x = Log(Te), y = log(rate) cm^3 s^-1
c
c       Total calculated: Log(T)<5   : Seaton 1959. MNRAS 119:81, eq 36
c                                      Accurate < 0.5%
c                         Log(T)>4.8 : Explicit bound-free wavefunction
c                                      overlap integrals, convolved with
c                                      a Maxwellian. For 210 levels up to
c                                      n = 20, l = 19. Accurate < 1%
c
c       n = 1 calculated: All explicit bound-free wavefunction
c                         overlap integrals, convolved with a Maxwellian.
c                         for to n = 1, l = 0. Accurate < alpha^2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      double precision function fkramer(ion,atom,t)
c
      include 'cblocks.inc'
c
      double precision t,tz2,h,z,z2,h3,rate,lt,lr,lx
      double precision h12,h43,h53,d,rate2,x1
c     
      integer*4 atom,ion,i,j
c
      fkramer = 0.d0
      rate = 0.d0
      rate2 = 0.d0
c
c     Check that the ion is really hydrogenic
c
      if ((mapz(atom)-ion+1).ne.0) then
         write (*,*) 'Warning, incorrect ion in fkramer'
         write (*,*) elem(atom),rom(ion)
         stop
      endif
c
c     OK, we are recombining to a H ion 
c     
c
c     get lambda (alias h)
c
      z = dble(mapz(atom))
      z2 = z*z
      tz2 = t/z2
      h = 157890.d0*z2/t
      h12 = dsqrt(h)
      h3 = h**(-1.d0/3.d0)
      h43 = 1.d0/(h3*h3*h3*h3)
      h53 = h43/h3
      d = 5.197d-14

      if (tz2.lt.1.d3) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Low temperature fit if t<z^2 10^4 K
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate the rate from asymptote....
c
      rate = d*z*h12*(0.4288d0+dlog(h)/2.d0+0.469d0*h3)
c      write(*,*) "low",z,t/z2,rate
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     High temperatures > 10^4.8 Z^2 K
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New spline fit to total rate in hspline(1,x,y) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      lt = dlog10(tz2)
      x1 = hspline(1,1,1)
      i = 1
      do j = 1,18
         if (lt.ge.hspline(1,j,1)) then
            x1 = hspline(1,j,1)
            i = j
         endif
      enddo
c
      lx = lt - hspline(1,i,1)
      lr = lx*hspline(1,i,5)
      do j = 4,3,-1
         lr = lx*(hspline(1,i,j) + lr)
      enddo
      lr = lr+hspline(1,i,2)
c
      rate = (z*(10.d0**lr))
c      write(*,*)"high",z,t/z2,rate
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      endif
c
c     Check to see if we have H or He and calculate excited state
c     recombination rates.
c
      if (mapz(atom).le.2) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     New spline fit to N=1 level in hspline(3,x,y)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      lt = dlog10(tz2)
      x1 = hspline(3,1,1)
      i = 1
      do j = 1,18
         if (lt.ge.hspline(3,j,1)) then
            x1 = hspline(3,j,1)
            i = j
         endif
      enddo
c
      lx = lt - hspline(3,i,1)
      lr = lx*hspline(3,i,5)
      do j = 4,3,-1
         lr = lx*(hspline(3,i,j) + lr)
      enddo
      lr = lr+hspline(3,i,2)
c
      rate2 = rate - (z*(10.d0**lr))
      if (rate2.lt.0.d0) rate2 = 0.d0
c
      if (mapz(atom).eq.1) rec(3,1) = rate2
      if (mapz(atom).eq.2) rec(5,2) = rate2
c
      endif
c
c     pass back the result
c
      if (rate.ge.0.d0) then
         fkramer = rate
      else 
         fkramer = 0.d0
      endif
c
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES  MU  INCLUDING ALL ELEMENTS .
c   USES ATOMIC MASSES  :  ATWEI   AND RELAT. ABUND.  :  ZION
c
c
c
      double precision function fmua(de, dh)
c
      include 'cblocks.inc'
c
      double precision de,dh
      double precision denstot,densnum
c
      fmua = (denstot(dh)+(de*me/amu))/(de+densnum(dh))
c
      return 
      end
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
c     function to evaluate the normalised two photon spectral probability 
c     function of Spitzer and Greenstein 1951 ApJ 114, 407-409
c
c     psi(x)/int(0-1) psi(x) dx
c
c     where int(0-1) psi dx = 3.770
c
c     This is required for the two photon gaunt
c     factor calculation
c
c     Note:
c
c     A(x) = 9 alpha^6 c Rhy/2^10 psi(x)   (for hydrogen)
c          = 4.363932296875 psi(x) in modern constants
c
c        x = nu/nu_2p 
c
c        0 <= x <= 1.0
c
c        function is only symmetrical in photons/unit interval.
c
c     I fit a constrained polynomial of 6th order to range 0 - 0.5 and
c     use symmetry for the rest.  This is significantly more accurate 
c     than the Mewe method.
c
c     errors: max error ~ 0.09 percent
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fpsiy(x)
c
      include 'const.inc'
c
      double precision x,p0,x2
      double precision a1,a2,a3,a4,a5,a6
c
c     assign polynomial constants
c
      a1 = 44.021893d0
      a2 = -228.498800d0
      a3 = 842.899670d0
      a4 = -1991.608300d0
      a5 = 2593.050000d0
      a6 = -1401.556300d0
c
      fpsiy = 0.0d0
      p0 = 0.d0
c
c     this is only valid for 0<= x <=1.0
c
      if ((x.ge.0.0).and.(x.le.1.0)) then
        if (x.gt.0.5d0) then
           x2 = 1.0-x
           p0 = (x2*(a1+x2*(a2+x2*(a3+x2*(a4+x2*(a5+x2*a6))))))
        else
           p0 = (x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6))))))
        endif
      else
         write (*,*) 'Error in fpsiy, x argument out of range 0-1'
         write (*,*) '!!!!!!!!!!!!!!!!!!!!!!x = ',x
         stop
      endif
c     
      fpsiy = p0*0.26525199d0
c     
      return
c     
      end
c
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
c       Calculates Weighted Photo-Ionisation Timescales
c
c       Weighted by population fractions and abundances
c       Assumes that the rates are all up to date. call
c       Allrates otherwise.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fphotim()
c
      include 'cblocks.inc'
c
      double precision rpho,totn,wei
      integer*4 i,j
c
      rpho = 0.d0
      totn = 0.d0
      fphotim = 0.d0
c
      do i = 1 ,  atypes
         totn = totn+zion(i)
      enddo
c
      do i = 1, atypes
c
      do j = 1 , maxion(i)-1
         wei = pop(j,i)*zion(i)
         rpho = rpho+wei*rphot(j,i)
      enddo
c     
      enddo
c
      fphotim = totn/(rpho+epsilon)
c
      return 
      end
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
c*******COMPUTES RADIATION PRESSURE Phi Sigma/c
c
c   Assumes tphot is up to date.  Tphot is in 1/4pi sr flux units
c   Also assumes xsec is up to date from totphot, includes dust 
c   radiation pressure.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fradpress(dr,dh)
c
      include 'cblocks.inc'
c
      double precision wid, blum, dr, tau, tlum
      double precision dustf, grainf
      double precision dustradP, ionradP,dh, dgrad
      integer*4 i,k,dtype
c
      blum = 0.d0
      tlum = 0.d0
      do i = IonStartBin,infph-1
         wid = (ephot(i+1)-ephot(i))*evplk
         tau = dabs(dr*xsec(i))
         if (tau.gt.1.d-8) then
         	blum = blum + tphot(i)*wid*(1.d0-dexp(-tau))
         else 
         	blum = blum + tphot(i)*wid*tau
         endif
         tlum = tlum + tphot(i)*wid
      enddo
c
      blum = blum*4.d0*pi
      tlum = tlum*4.d0*pi

      if(grainmode)then
       dustf=0.d0
       do dtype=1,numtypes
        do k=mindust(dtype),maxdust(dtype)
         do i=1,infph-1
           if(skipbin(i)) goto 1010
           wid = (ephot(i+1)-ephot(i))*evplk
           dgrad = gradedge(k+1) - gradedge(k)
           grainf=(absorp(i,k,dtype)+(1-Gcos(i,k,dtype))
     &           *scatter(i,k,dtype))*dustsig(k,dtype)*dgrad 
           dustf=dustf + tphot(i)*wid*(grainf)
 1010      continue
         enddo
        enddo
       enddo 
       dustradP=4.d0*pi*dh*dustf*dr/cls
      else
       dustradP=0.d0
      endif

      ionradP=blum/cls

c      write(*,'(" Dust, Ion  Rad Pressure: ",2(1pg12.5))') 
c     &  dustradP, ionradP
c
      fradpress = ionradP+dustradP
c
      return 
      end

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
c*******COMPUTES PRESSURE USING FUNCTION FELDENS
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function fpressu(t, dh, popul)
c
      include 'cblocks.inc'
c
      double precision t, dh
      double precision popul(mxion, mxelem),de
      double precision en
c
c           Functions
c
      double precision feldens
c
      de = feldens(dh,popul)
c
      en = zen*dh+de
c
      fpressu = en*rkb*t
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO DERIVE AVERAGE RECOMBINATION TIME SCALE
c   ASSUMING CONSTANT ELECTRONIC DENSITY : DE
c
c
c
      double precision function frectim(t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision t, de, dh
      double precision recab(mxion)
      double precision a,ad,adlt,ar,b,bd,c,d
      double precision et,f,rain,ta,tb,u,u4,wei,weito
c
      integer*4 i,ion,j,jj,mif,mij
c
c           Functions
c
      double precision adf,adltf,arf
c
      arf(u,ar,et) = ar*((u/1.0d4)**(- et))
      adf(u,ad,bd,ta,tb) = 
     & (((1.0+(bd*dexp(-(tb/u))))*dexp(-(ta/u)))*ad)*(u**(-1.5d0))
      adltf(u4,a,b,c,d,f) = 
     & ((1.0d-12*((((a/u4)+b)+(c*u4))+((d*u4)*u4)))*(u4**(-1.5)))
     & *dexp(-(f/u4))
c
      u = dmax1(0.01d0,t)
      u4 = u/1.d4
      ion = atypes
      rain = 0.0d0
      weito = 0.0d0
c
      do 100 i = 1, atypes
         mif = 0
         mij = 2
         do 120 j = mij, maxion(i)
            recab(j) = 0.0d0
            ar = arad(j,i)
            if (ar.gt.0.0d0) then
               ad = adi(j,i)
               et = xrad(j,i)
               ta = t0(j,i)
               bd = bdi(j,i)
               tb = t1(j,i)
               a = adilt(j,i)
               b = bdilt(j,i)
               c = cdilt(j,i)
               d = ddilt(j,i)
               f = fdilt(j,i)
               adlt = adltf(u4,a,b,c,d,f)
               if (adlt.lt.0.0) adlt = 0.0d0
               recab(j) = (arf(u,ar,et)+adf(u,ad,bd,ta,tb))+adlt
               mif = j
            else
               mij = j+1
            end if
 120     continue
c     
         if (mif.ge.mij) then
            wei = 1.0d0
            weito = weito+(zion(i)*wei)
            do  j = mij, mif
               do jj = mij, j
                  rain = rain+(wei/(recab(jj)+1.d-36))
               enddo
            enddo
         end if
 100  continue
c
      frectim = rain/(((de*dh)*weito)+(1.d-36*rain))
c
      return 
      end
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
c       Calculates Weighted Recombination timescale
c
c       Weighted by population fractions and abundances
c       Assumes that the rates are all up to date. call
c       Allrates otherwise.
c
c   ASSUMES CONSTANT ELECTRONIC DENSITY : DE
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function frectim2(de)
c
      include 'cblocks.inc'
c
      double precision de,rrec,totn,wei
      integer*4 i,j
c
      rrec = 0.d0
      totn = 0.d0
      frectim2 = 0.d0
c
      do i = 1 ,  atypes
         totn = totn+zion(i)
      enddo
c
      do i = 1, atypes
c
      do j = 2 , maxion(i)
         wei = pop(j,i)*zion(i)
         rrec = rrec+(wei*rec(j,i)) 
      enddo
c     
      enddo
c
      frectim2 = totn/(de*rrec+epsilon)
c
      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c***************************************************
c
c   This routine is to calculate the gaunt factors
c   for the new XLINLIST resonance line data.
c
c   The five parameter fit from Landini&Monsignori Fosse 1990
c   is used with extentions for He & Li like ions from Mewe 1981.
c
c   A new (rational) lookup table code has been implemented
c   to replace the unusual one used by L&M...
c
c     NOTE: Helio will still call the old fgaunt
c     sign errors for C and y^2D corrected
c
c
c
c     RSS 9/1990.
c***************************************************
c
      double precision function fresga(atom,isoe,trans,ejk,telec, code)
c
      include 'cblocks.inc'
c
      double precision f1,f2,f3
      integer*4 isoe,atom,fclookup,code,trid
      double precision ejk,telec,y,e1y
      double precision A,B,C,D,E,trans,z
c
c           Functions
c
      double precision farint
c
c   internal functions (ref Mewe 1981 A.A.Suppl.Ser. 45:11  pg 20)
c   f4 not used since the y>>1 approx not used
c
      f1(z) = 1.d0 - 1.3d0/z
      f2(z) = 1.d0 + 1.5d0/z
      f3(z) = 1.d0 + 7.0d0/z
c
      fresga = 0.2d0
c
c   find y and the first exponential integral (see farint.f)
c
      y = ejk*ev/(rkb*telec)
      e1y = farint(1,y)
c
c
c   get code for transition
c
      code = fclookup(isoe,trans)
c
c
      if (code.gt.-1) then
c
c
c
      if (code.lt.1) code = 1
      if (code.gt.37) code = 1
c
c
c   get parameters from data arrays
c
      A = arg(code)
      B = brg(code)
      C = crg(code)
      D = drg(code)
      E = erg(code)
c
c
c
c      if (mapz(atom).eq.2) then
c         write (*,*) atom,isoe,trans,ejk,telec
c         write (*,*) code,A,B,C,D,E,e1y
c      endif
c
c
      if (isoe.ne.2) then 
          fresga = A+(E+y*(B+y*(y*D-C)))*e1y+y*(C+D)-y*y*D
      else
         trid = idint(trans)
         z = mapz(atom)
          if (trans.le.4) then
              A = A*f1(z)
              D = D*f1(z)
          endif
          if (trans.eq.5) then
              C = C*f2(z)
              D = D*f3(z)
          endif
          if (trans.eq.6) then
              C = C*f2(z)
              D = D*f2(z)
          endif
          fresga =  A+(E+y*(B+y*(y*D-C)))*e1y+y*(C+D)-y*y*D
      endif
c
      else
c
c     transition from meta-stable level fef = omega/statwei already
c
         fresga = 1.d0
c
      endif
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******FUNCTION TO COMPUTE DENSITY IN G/CC USING NUMBER DENSITY
c   USES FUNCTION DENSTOT (amus)
c
c
      double precision function frho(de, dh)
c
      include 'cblocks.inc'
c
      double precision de, dh
      double precision denstot
c
      frho = (amu*denstot(dh))+(me*de)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES INTEGRAL DEXP(-U)/U (FROM ZERO TO INFINITY)
c   TIMES DEXP(U)*EXP(-REXP)
c
c     **REF:  HANDBOOK OF MATH. FUNCTIONS,ABRAMOVITZ AND STEGUN,
c          1970,(AMS 55),P231
c
c
c
      double precision function fue1(u, rexp)
c
      include 'const.inc'
c
      double precision u, rexp

      if (u.gt.1.0d0) goto 100
c
      fue1 = (((((((-57.72157d0+99.99919d0*u)
     &-(24.99106d0*u*u))+(5.51997d0*u*u*u))
     &-(0.976d0*u*u*u*u))
     &+(0.10786d0*u*u*u*u*u))
     &-(100.d0*dlog(u)))/100.d0)
     &*dexp(u)
      goto 500
c
  100 fue1 = (((u*u)+(2.334733d0*u))+0.250621d0)
     &/((((u*u)+(3.330657d0*u))+1.681534d0)*u)
c
  500 fue1 = fue1*dexp(-rexp)
c
      return 
      end
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
c     Function to calculate the coulomb integral at frequency nu
c
c     ref: Krolick and McKee 1978 ApJS 37,459
c         
c   fromt: Johnston and Dawson 1973, Phys Fluids 16 722 (not sighted)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function lambda(t,de,dh,nu)
c
      include 'cblocks.inc'
c
      double precision t,de,dh,nu
      double precision con,x,y
      double precision nuplasma,ear
      double precision zp,nu0,zav
c
c           Functions
c
      double precision favcha
c
      lambda = 0.d0
c
      con = (3.d0/16.d0)
c
c
      x = dsqrt((ryde)/(rkb*t))
      y = 1.d0/x
      ear = 1.d0/(de*bohr**3.d0)
      nuplasma = dsqrt((pi*de*eesu)/me)/pi
      nu0 = nuplasma/nu
      zav = favcha(pop,1)
      zp = rt2*zav*x*x
c
      lambda = con*y*ear*(dmin1(1.d0,nu0)/dmax1(zp,x))
c
      return
c
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE BNU (ERGS.S-1.CM-2.HZ-1.SR-1)
c   FROM PLANK ENERGY DISTRIBUTION
c   TS : TEMP.    RNUH : FREQUENCY IN RYDBERGS (NU/NUHYD)
c
c
c
      double precision function fplank(ts, rnuh)
c
      include 'const.inc'
c
      double precision a,b,c
      double precision ts, rnuh
c
c   1 Ryd = 13.59844eV
c
c     F = pi B(nu) = 2 pi h c^-2 nu^3 /(exp(-h nu/k t)-1)
c
c     2(Rydberg ev/h)^3 h c^-2 = 0.5241732391691488 
c
c     Ryd ev/rkb = 1.578023832938005e+5
c
c     this routine returns B (1/pi Inu units)
c
c
      a = 0.5241732391691488*(rnuh*rnuh*rnuh)
      b = (1.578023832938005d+5*rnuh)/ts
      if (b.gt.1.0d2) goto 50
      c = 1.0d0/(dexp(b)-1.0d0)
      goto 100
   50 c = dexp(- b)
  100 fplank = a*c
c
      return 
      end
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
c     fclookup: converts the L&M transition code + iso series
c     to an ordinal value for fresgau.
c
c     NOTE: transitions such as 14c have been entered as 14.3
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      integer*4 function fclookup(isos,tran)
c
      include 'cblocks.inc'
c
      integer*4 isos, trid,code
      double precision tran,trd
c
c    assign defaults.
c
      fclookup = 1
      code = 0
c
c    Sorry, this is going to be a bit dirty...
c
c    metastable transitions are handled elsewhere
c     (inter2.f)
c
      trid = idint(tran)
      trd = dble(trid)
c
c    find if we have a whole number
c
      if (trd.eq.tran) then
c
c    we have a whole number, use the code array
c
         code = clkup(isos,trid)
      else
c     
c      handle case by case, assuming prior knowlege of allowed (sic)
c     transition numbers in XLINDAT....
c
c        isos  tran  rcode gcode
c          4   14.1    14   100  
c          4   14.2    20   200  
c          4   14.3    20   200  
c         10    1.1     9    53  
c         10    3.1    12    58  
c         10    3.2    13    61  
c         12    6.1    17   103  
c
         if (isos.eq.4) then
            code = 20
            if (tran.eq.14.1d0) code = 14
         endif
         if (isos.eq.10) then
            code = 9
            if (tran.eq.3.1d0) code = 12
            if (tran.eq.3.2d0) code = 13
         endif
         if (isos.eq.12) then
            code = 17
         endif
      endif
c     
c    pass back result
c
      fclookup = code
c
c
c
      return
c
      end
