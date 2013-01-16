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
c     v 2.2 new charge exchange reactions method.
c
c
c     RSS 4/91
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******DERIVES IONIC POPULATIONS FOR HYDROGEN AFTER TIME : TSTEP
c	USES INITIAL FRACTIONAL IONISATION  :  XHY (OR POP(2,1))
c	DEFAULT VALUE (USING POP(2,1)) ASSIGNED TO XHY IF SET TO -1
c	FINAL ELECTR. DENSITY : DEF ; FINAL FRACTIONAL IONIS.: XHYF
c	FINAL IONIC POPULATIONS IN POP(6,11) AND DERIVATIVES IN DNDT
c	CALL SUBR. ALLRATES,SDIFEQ,IONSEC
c
c
      subroutine iohyd(dh, xhy, t, tstep, def, xhyf,mod)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision co,rc,phi,rkn,eps,abu
      double precision dh, xhy, t, tstep, def, xhyf
      double precision a1,a2,aa,b1,b2,c1,c2
      double precision da,de,de2,del,dep,dey,dma,fhi
      double precision fhii,ph,ph2,spo,yh,dehyd
c
c     Old charge exchange  
c      double precision chex
c
      integer*4 at,ion,j,iel,ies
c
      character jjmod*4, mod*4
c
c           Functions
c
      double precision feldens
c
      if (tstep.ge.0.0d0) goto 300
      write(*, 200) 
  200 format(' STOP !!!!!!!!!!!!!!!!!! NEGATIVE TIME STEP'//)
      stop 
  300 continue
c
      if ((jspot .eq. 'YES').or.(jspot .eq. 'NO')) goto 60
      write(*, 67) jspot
   67 format(//' FLAG "JSPOT" IMPROPERLY DEFINED : ',a4)
      stop 
   60 continue
c
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      jjmod = 'ALL'
c
c    ***INITIAL CONDITIONS
c
      call allrates(t, jjmod)
      fhii = pop(2,1)
      fhi = pop(1,1)
      if ((xhy.lt.0.0d0).or.(xhy.gt.1.0d0)) goto 327
      fhii = xhy
      fhi = 1.d0-xhy
  327 continue
      dehyd = feldens(dh,pop)
c
c     **COL- , PHOTO- IONISATION AND RECOMBINATION RATES
c
      del = dmax1(0.d0,dehyd-((dh*zion(1))*fhii))/dh
      co = dh*col(1,1)
      rc = dh*rec(2,1)
      phi = rphot(1,1)
c     
c    **IONISATION RATE CORRECTION DUE TO THE ON THE SPOT APPROX.
c
      spo = 0.0d0
c
      if (jspot .ne. 'YES') goto 113
c
      call spotap(de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
      rc = dh*rec(3,1)
      aa = (dh*zion(2))*pop(2,2)
      spo = aa*((yh*(rec(2,2)-rec(4,2)))+(ph*rec(4,2)))
      aa = (dh*zion(2))*pop(3,2)
      spo = spo+((aa*ph2)*rec(5,2))
  113 continue
c
c     **IONISATION RATE CORRECTION DUE TO NON-THERMAL ELECTRONS
c
      call ionsec
      phi = phi+rasec(1,1)
c
c    ***ADD CHARGE EXCHANGE RATES BETWEEN HYDROGEN AND HEAVIER ELEMENTS
c
      rkn = 0.0d0
      eps = 0.0d0
c
      if (chargemode.eq.0) then 
c
c     MAPPINGS 2.2 charge rates
c
c
c     First, recombination reactions
c     with neutral H, ionise H
c
         do j = 1, nchxr
c
            if (chxr(j).gt.0.d0) then
c
c     have a non zero rate
c
               at = chxrat(j)
               ies = chxrx(j)
c
               if (ies.eq.zmap(1)) then
c
c     have H reaction
c
                  ion = chxrio(j)
                  abu = dh*zion(at)
c
c
                  if (abu*pop(ion,at).ge.pzlimit) then
c     
c     significant ions only
c
                     eps = eps+((chxr(j)*abu)*pop(ion,at))

                  endif
               endif
            endif
c
c        end j, ionising neutral H reactions
c
         enddo
c
c     Ionising reactions
c     with ionised H
c
c
         do j = 1, nchxi
c
            if (chxi(j).gt.0.d0) then
c
c     have a non zero rate
c
               at = chxiat(j)
               ies = chxix(j)
c
               if (ies.eq.zmap(1)) then
c
c     have H reaction
c
                  ion = chxiio(j)
                  abu = dh*zion(at)
c
                  if (abu*pop(ion,at).ge.pzlimit) then
c     
c     significant ions only
c
                     rkn = rkn+((chxi(j)*abu)*pop(ion,at))
c
                  endif
               endif
            endif
c
c        end j, recom ionised H reactions
c
         enddo
c
c     end char exchange rates
c
      endif
c
      if (chargemode.eq.1) then
c
c     old mappings charge exchanges
c
      do 40 j = 1, nchxold
         iel = idint(charte(1,j))
         if (iel.lt.2) goto 40
         ion = idint(charte(2,j))
         if (rec(ion+1,iel).le.0.0d0) goto 40
         ies = idint(charte(5,j))
         if (ies .ne. 1) goto 40
         abu = dh*zion(iel)
         if (abu*pop(ion,iel).lt.pzlimit) goto 40
         if (abu*pop(ion+1,iel).lt.pzlimit) goto 40
         rkn = rkn+((charte(3,j)*abu)*pop(ion,iel))
         eps = eps+((charte(4,j)*abu)*pop(ion+1,iel))
   40 continue
c
c     end old charge reactions
c
      endif
c
c
c     for equilibrium we need only the rate ratios
c
c      if (mod.eq.'EQUI') then
c
c         rkn = rkn*fhii
c         eps = eps*fhi
c
c         fhi = (rc+rkn)/(co+spo+phi+eps+rc+rkn)
c         fhii = (co+spo+phi+eps)/(co+spo+phi+eps+rc+rkn)
c
c    
c      else
c     
c     ***SOLVE DIFFERENTIAL EQUATION
c     USING CHARGE EXCHANGE AS ONE COMPONENT
c
c      eps = (eps/(fhii+epsilon))
c      rkn = (rkn/(fhi+epsilon))
c     
c        chex = eps-rkn
c
c     
c     ***SET COEFFICIENTS A,B,C FOR THE DIFFERENTIAL EQUATION
c     
c         a1 = (del*(co+spo))+phi
c         b1 = ((1.d0-del)*co)-phi+spo-(del*rc)+chex
c         c1 =-(rc+co+chex)
c
c         a2 = (1.d0+del)*(rc-spo)
c         b2 = ((((-(co*(del+1.d0)))-phi)+spo)-(rc*(del+2.d0)))-chex
c         c2 =-c1
c
c     ***SOLVE DIFFERENTIAL EQUATION
c     USING CHARGE EXCHANGE AS TWO SEPERATE PARTS
c
c     ***SET COEFFICIENTS A,B,C FOR THE DIFFERENTIAL EQUATION
c     


         a1 = (del*(co+spo))+phi+eps
         b1 = ((1.d0-del)*co)-phi+spo-(del*rc)-(eps+rkn)
         c1 =-(rc+co)

         a2 = (1.d0+del)*(rc-spo)+rkn
         b2 = (((-(co*(del+1.d0)))-phi)+spo)-rc*(del+2.d0)-(eps+rkn)
         c2 =-c1
c     
         call sdifeq(a1, b1, c1, fhii, tstep)
c     
         call sdifeq(a2, b2, c2, fhi, tstep)
c     
c      endif
c     
c
c    ***WARNING MESSAGE AND NORMALISATION IF ABUNDANCES
c	DO NOT CONSERVE
c
      da = fhi+fhii
      dma = 1.d2*dabs(1.d0-da)
      if (dma.ge.1.d0) then
       write(*, 940) dma
  940  format(' WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     & ' ABUND. OF H DID NOT CONSERVE BEFORE NORMALISATION :',g9.3,2h %
     & //)
       write(*,*) rc,co, phi
       write(*,*) rphot(1,1),rasec(1,1) 
       write(*,*) del, dehyd, dh
       write(*,*) pop(2,1),pop(1,1)
       fhi = fhi/da
       fhii = fhii/da
      endif
c      if (fhii.lt.1.d-1) fhi = 1.d0-fhii
c      if (fhi.lt.1.d-1) fhii = 1.d0-fhi
c
c    ***SET OUTPUT VARIABLES AND DETERMINE TIME DERIVATIVE
c
      pop(1,1) = fhi
      pop(2,1) = fhii
      xhyf = fhii
      def = feldens(dh,pop)
      dndt(1,1) = dh*((((c2*fhi)+b2)*fhi)+a2)
c
      dndt(2,1) = dh*((((c1*fhii)+b1)*fhii)+a1)
      return 
      end

