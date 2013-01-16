cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES BALMER RECOMBINATION LINES FOR HI (CASE B)
c   AND HEI (CASE B), HEII 4686 (CASE B)
c
c     **REFERENCES : PENGELLY,R.M.(1964) MNRAS,V.127,P.145
c                    BROCKLEHURST,M.(1971) MNRAS,V.153,P.471
c
c
      subroutine hydrec(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision t, telc, de, dh, hbet
      double precision a1,a2,a3,a4,a5,ab1,ab2,abde
      double precision amul,argu,b1,b2,bex,cormax,delc,delmu
      double precision delo10,deloef,delomi,dlo,efden,fd,fmden,fnorm
      double precision hbt0,he4471
      double precision palph,raln,rmm,t5
      double precision temul,tinv,tlo,tlolo,tmax,tmin,tspw,x
c
      integer*4 k,l
c
c
c           Functions
c
      double precision fcor,fhb,fhb0,fli,fun
c
      fhb(dlo,t5,tlo) = 1.d-25*(dexp(3.8328d0+0.083585d0*tlo
     &-0.05177d0*tlo*tlo))*(0.995d0*(1.d0+(dexp(-6.67d0
     &+0.25516d0*dlo)*(t5 ** (-1.2d0)))))
c
      fhb0(dlo,t5,tlo,fd) = (0.2392d-25*dexp((3.8328d0+
     &(0.083585d0*tlo))-((0.05177d0*tlo)*tlo)))*(0.995d0*
     &(1.0+(fd*(dexp((-6.67d0)+(0.25516d0*dlo))*(t5**(-1.2d0))))))
c
      fli(x,a1,a2,a3,a4,a5) = (((a5*x+a4)*x+a3)*x+a2)*x+a1
      fcor(x,a1,a2,a3,a4) = (((((a4*x)+a3)*x)+a2)*x)+a1
      fun(x,b1,b2) = b1+(b2*x)
c
      tmin = 500.d0
      tmax = 2.d5
      delomi = -1.d0
c
      cormax = 0.25d0
c
      tlo = dlog(t+1.d-3)
      dlo = dlog(de+1.d-6)
      telc = dmax1(tmin,dmin1(tmax,t))
      t5 = telc/5.d3
      tinv = 1.d4/telc
      tlolo = dlog(dlog(telc))
      ab1 = -19.10033d0
      ab2 = 8.705114d0
      argu = fun(tlolo,ab1,ab2)
      delo10 = dmax1(delomi,dlog10(de+1.d-3))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***FINDS ABSOLUTE FLUX FOR H-BETA
c
      abde = ((de*dh)*zion(1))*pop(2,1)
      hbet = (abde*fhb(dlo,t5,tlo))/(4.d0*pi)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***DERIVES THE OTHER LINES FLUX USING THEIR PREDICTED
c   RATIOS RELATIVE TO H-BETA
c
      do l = 1, 5
          a1 = balcoe(l,1)
          a2 = balcoe(l,2)
          a3 = balcoe(l,3)
          a4 = balcoe(l,4)
          a5 = balcoe(l,5)
          raln = fli(argu,a1,a2,a3,a4,a5)
          a1 = balcoe(l,6)
          a2 = balcoe(l,7)
          a3 = balcoe(l,8)
          a4 = balcoe(l,9)
          b1 = balcoe(l,10)
          b2 = balcoe(l,11)
          deloef = delo10*fun(tinv,b1,b2)
          delc = fcor(deloef,a1,a2,a3,a4)
          if (dabs(delc).gt.cormax) delc = dsign(cormax,delc)
          hbri(l) = hbet*(dexp(raln)*(1.d0+delc))
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***HEII 4686 , COMPUTES FIRST P-ALPHA AT TE/4 (USING H-ALPHA)
c
      tlo = dlog((t/4.d0)+1.d-3)
      telc = dmax1(tmin,dmin1(tmax,t/4.d0))
      t5 = telc/5.d3
      tinv = 1.d4/(telc/2.d0)
      tlolo = dlog(dlog(telc/2.d0))
      abde = ((de*dh)*zion(2))*pop(3,2)
      hbt0 = fhb(dlo,t5,tlo)/(4.d0*pi)
c
      argu = fun(tlolo,ab1,ab2)
      a1 = balcoe(1,1)
      a2 = balcoe(1,2)
      a3 = balcoe(1,3)
      a4 = balcoe(1,4)
      a5 = balcoe(1,5)
      raln = fli(argu,a1,a2,a3,a4,a5)
c
      a1 = balcoe(1,6)
      a2 = balcoe(1,7)
      a3 = balcoe(1,8)
      a4 = balcoe(1,9)
      b1 = balcoe(1,10)
      b2 = balcoe(1,11)
      tspw = 4.d0*dsqrt(tinv/16.d0)
      deloef = delo10*fun(tspw,b1,b2)
      delc = 3.d0*fcor(deloef,a1,a2,a3,a4)
      delc = dmin1(cormax,dmax1(-cormax,delc))
      palph =((0.102d0*(tinv**0.1716d0))*hbt0)*(dexp(raln)*(1.d0+delc))
c
      heiilr = (8.d0*abde)*palph
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***FINDS FLUX OF HEI 4471 USING H-BETA FUNCTION
c
      tlo = dlog(1.d4*((((t/2.d0)+1.d-3)/1.d4)**1.264d0))
      telc = dmax1(tmin,dmin1(tmax,t/2.d0))
      t5 = telc/5.d3
      tinv = 1.d4/telc
      abde = ((de*dh)*zion(2))*pop(2,2)
      fmden = zion(2)
      hbt0 = fhb0(dlo,t5,tlo,fmden)/(4.d0*pi)
      he4471 = abde*hbt0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***DERIVES OTHER HEI LINE FLUXES USING RELATIVE RATIO TO H-ALPHA
c
      a1 = balcoe(1,6)
      a2 = balcoe(1,7)
      a3 = balcoe(1,8)
      a4 = balcoe(1,9)
      b1 = balcoe(1,10)
      b2 = balcoe(1,11)
      deloef = delo10*fun(tinv,b1,b2)
      delc = fcor(deloef,a1,a2,a3,a4)
      delc = dmin1(cormax,dmax1(- cormax,delc))
c
      a1 = balcoe(1,1)
      a2 = balcoe(1,2)
      a3 = balcoe(1,3)
      a4 = balcoe(1,4)
      a5 = balcoe(1,5)
c
      k = 6
      do l = 1, 12
c      if ((hel0(2,l).gt.0.d0).or.(j108m.eq.l)) then
      temul = dabs(hel0(2,l))
      tlolo = dlog(dlog((temul*telc)*2.d0))
      argu = fun(tlolo,ab1,ab2)
      if (temul .eq. 0.5d0) then
      fnorm = 3.04804d0
      else if (temul.eq.0.25d0) then
      fnorm = 3.30346d0
      else if (temul .eq. 0.1d0) then
      fnorm = 3.74000d0
      else
      fnorm = 3.04804d0*((temul/0.5d0)**(-0.12712d0))
      end if
      raln = fli(argu,a1,a2,a3,a4,a5)
c
      rmm = dexp(raln)/fnorm
      amul = hel0(3,l)
      bex = hel0(4,l)
      delmu = hel0(5,l)
      efden = (1.d0+(delmu*delc))*(1.d0+(4.3d-3*delmu))
      hel0(6,l) = (he4471*(amul*(rmm ** bex)))*efden
      heibri(l+3) = hel0(6,l)
      
c      if ((j108m .ne. l).and.(hel0(2,l).gt.0.d0)) then
c      hbri(k) = hel0(6,l)
c      k = min0(10,k+1)
c      end if
      
c      end if
      enddo
c
      return 
      end
