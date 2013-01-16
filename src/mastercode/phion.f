cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO CALCULATE PHOTOIONISATION RATES IN RPHOT(mxion, mxelem) ,
c	PHOTOHEATING RATES IN HEAPH(5,28,16) AND
c	BASIS FOR SECONDARY IONISATION IN ANR(2,28,16),WNR(2,28,16)
c
c	USES LOCAL MEAN INTENSITY OF RADIATION JNU CONTAINED
c	IN VECTOR : TPHOT     ( NUMBER OF CROSS SECTIONS : IONUM )
c	IT MULTIPLIES JNU BY 4PI AND INTEGRATES THE NUMBER OF PHOTONS
c
c	USES DOUBLE PRECIS. FUNCTION : DFUINT  TO INTEGRATE
c	CALL SUBR. INTERPOL
c
c
c
      subroutine phion()
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision alh,alo,alpha,augen,augerg,augt,bet,cflo
      double precision daiau,daip,dait,e1,e2,eau2,eau2ev,ed,etr
      double precision piplklo,rnutr,rr1,rr2,sa,satmp
      double precision Cedge,Oedge
c
      integer*4 i,ie,inl,j,n
c
c           Functions
c
      double precision dfuint
c
c   **  EFF(A,B,ED)=B+1.0/(ED/A+1.0/(1.0-B+1.0D-12))
c   **  COEFFICIENTS A,B FROM SHULL,M.J. APJ.234,P761 (1979)
c
      double precision dain, daih, eff1, eff2, eff3, eff4, eff5
c
c
      eff1(ed) = 1.0d0
      eff2(ed) = 0.708d0+(1.0d0/((ed/31.0d0)+(1.0d0/0.292d0)))
      eff3(ed) = 0.383d0+(1.0d0/((ed/16.2d0)+(1.0d0/0.617d0)))
      eff4(ed) = 0.194d0+(1.0d0/((ed/7.20d0)+(1.0d0/0.806d0)))
      eff5(ed) = 0.113d0+(1.0d0/((ed/3.25d0)+(1.0d0/0.887d0)))
c
c
c     ln(4pi/h)
c
      piplklo = 62.8098090594892475d0
c     
      Cedge=280.d0
      Oedge=533.d0
      eau2 = 30.0d0
      eau2ev = eau2*ev
      qtosoh = 0.0d0
c
      do i = 1, atypes
         do j = 1, 3
            rasec(j,i) = 0.0d0
         enddo
         do j = 1, 4
 5          auphot(j,i) = 0.0d0
         enddo
         do j = 1, maxion(i)-1
            rphot(j,i) = 0.0d0
            do n = 1, 2
               anr(n,j,i) = 0.0d0
               wnr(n,j,i) = 0.0d0
               heaph(n,j,i) = 0.0d0
            enddo
            do n = 3, 5
               heaph(n,j,i) = 0.0d0
            enddo
         enddo
      enddo
c     
c
c    If there is dust include atoms contained in dust above Auger limit
c
c
      if (.NOT.grainmode) then !(grainmode = False)
     
       do 500 inl = IonStartBin, infph-1
c
c note test for zero is a problem in double precis
c tphot bins should be truncated to zero in totphot now
c
c now includes a comparison between photons/cm2/s and ions/cm3
c and if this is less than 1 cm/s then skipbin is true
c this is tested and set in totphot
c
c      if (tphot(inl).le.0.0d0) goto 500
         if (skipbin(inl)) goto 500
         e1 = ephot(inl)
         e2 = ephot(inl+1)
c     
c     ***FIND INTERPOLATION VALUES IN THE INTERVAL E2-E1 (IN EV)
c     FOR THE MEAN INTENSITY JNU ASSUMING A POWER LAW FUNCTION
c     OF FREQUENCY
c     
         call interpol(inl, cflo, alpha, tphot)
c     
c     
c     ***INTEGRATES IN FREQUENCY FOR EACH ENERGY BIN
c     
         alh = alpha-1.0d0
         rr2 = e2/ryd
         if (rr2.lt.1.0001d0) goto 333
         rr1 = e1/ryd
         if (rr1.lt.1.0d0) rr1 = 1.0d0
         bet = 1.0d0
         alo = (piplklo+cflo)+((1.0+alh)*dlog(ryd*evplk))
         qtosoh = qtosoh+dfuint(bet,alh,rr1,rr2,alo)
 333     continue
c     
         alh = alpha-1.0d0
         do 400 i = 1, ionum
            etr = ipotpho(i)
            rr2 = e2/etr
            if (rr2.lt.1.0001d0) goto 500
            rr1 = e1/etr
            if (rr1.lt.1.0d0) rr1 = 1.0d0
            ie = atpho(i)
            j = ionpho(i)
            augen = augpho(i)
            bet = betpho(i)
            sa = alh-spho(i)
            rnutr = etr*evplk
c     
c     ***DERIVES PHOTOIONISATION RATE
c     
            alo = ((piplklo+dlog(sigpho(i)))+cflo)+((1.d0+alh)*dlog
     &           (rnutr))
            dain = dfuint(bet,sa,rr1,rr2,alo)
            if (dain.le.0.d0) goto 400
            if (augen.le.0.0d0) then
               rphot(j,ie) = rphot(j,ie)+dain
            else
               auphot(j,ie) = auphot(j,ie)+dain
            end if
c     
c     ***INTEGRATES HEAPH TO BE USED IN SUBR. : PHEAT
c     AND ANR,WNR  FOR THE IONISATION CORRECTION .
c     
            if (ie.gt.limph) goto 400
            satmp = sa+1.d0
            daih = dfuint(bet,satmp,rr1,rr2,alo)
            daip = ((daih-dain)*etr)*ev
            if (daip.gt.0.0d0) then
               ed = dmax1(1.d-3,(((rr1*etr)+e2)/2.d0)-(etr+10.2d0))
               heaph(1,j,ie) = heaph(1,j,ie)+(eff1(ed)*daip)
               heaph(2,j,ie) = heaph(2,j,ie)+(eff2(ed)*daip)
               heaph(3,j,ie) = heaph(3,j,ie)+(eff3(ed)*daip)
               heaph(4,j,ie) = heaph(4,j,ie)+(eff4(ed)*daip)
               heaph(5,j,ie) = heaph(5,j,ie)+(eff5(ed)*daip)
               dait = daip-(dain*ryde)
               if (dait.ge.0.0d0) then
                  anr(1,j,ie) = anr(1,j,ie)+dait
                  wnr(1,j,ie) = wnr(1,j,ie)+dain
                  dait = daip-(dain*eau2ev)
                  if (dait.ge.0.0d0) then
                     anr(2,j,ie) = anr(2,j,ie)+dait
                     wnr(2,j,ie) = wnr(2,j,ie)+dain
                  end if
               end if
            end if
c     
            if (augen.gt.ryd) then
               augt = augen-10.2d0
               augerg = augen*ev
               daiau = dain*augerg
               heaph(1,j,ie) = heaph(1,j,ie)+(eff1(augt)*daiau)
               heaph(2,j,ie) = heaph(2,j,ie)+(eff2(augt)*daiau)
               heaph(3,j,ie) = heaph(3,j,ie)+(eff3(augt)*daiau)
               heaph(4,j,ie) = heaph(4,j,ie)+(eff4(augt)*daiau)
               heaph(5,j,ie) = heaph(5,j,ie)+(eff5(augt)*daiau)
               dait = dain*(augerg-ryde)
               if (dait.ge.0.0d0) then
                  anr(1,j,ie) = anr(1,j,ie)+dait
                  wnr(1,j,ie) = wnr(1,j,ie)+dain
                  dait = dain*(augerg-eau2ev)
                  if (dait.ge.0.0d0) then
                     anr(2,j,ie) = anr(2,j,ie)+dait
                     wnr(2,j,ie) = wnr(2,j,ie)+dain
                  end if
               end if
            end if
c     
 400     continue
c     
 500   continue
 
      else !(grainmode = True)
      
       do 1000 inl = IonStartBin, infph-1
c
c note test for zero is a problem in double precis
c tphot bins should be truncated to zero in totphot now
c
c now includes a comparison between photons/cm2/s and ions/cm3
c and if this is less than 1 cm/s then skipbin is true
c this is tested and set in totphot
c
c      if (tphot(inl).le.0.0d0) goto 500
         if (skipbin(inl)) goto 1000
         e1 = ephot(inl)
         e2 = ephot(inl+1)
c     
c     ***FIND INTERPOLATION VALUES IN THE INTERVAL E2-E1 (IN EV)
c     FOR THE MEAN INTENSITY JNU ASSUMING A POWER LAW FUNCTION
c     OF FREQUENCY
c     
         call interpol(inl, cflo, alpha, tphot)
c     
c     
c     ***INTEGRATES IN FREQUENCY FOR EACH ENERGY BIN
c     
         alh = alpha-1.0d0
         rr2 = e2/ryd
         if (rr2.lt.1.0001d0) goto 666
         rr1 = e1/ryd
         if (rr1.lt.1.0d0) rr1 = 1.0d0
         bet = 1.0d0
         alo = (piplklo+cflo)+((1.0+alh)*dlog(ryd*evplk))
         qtosoh = qtosoh+dfuint(bet,alh,rr1,rr2,alo)
 666     continue
c
         do 800 i = 1, ionum
            etr = ipotpho(i)
            rr2 = e2/etr
            if (rr2.lt.1.0001d0) goto 1000
            rr1 = e1/etr
            if (rr1.lt.1.0d0) rr1 = 1.0d0
            ie = atpho(i)
            j = ionpho(i)
            augen = augpho(i)
            bet = betpho(i)
            sa = alh-spho(i)
            rnutr = etr*evplk
c     
c     ***DERIVES PHOTOIONISATION RATE
c     
            alo = ((piplklo+dlog(sigpho(i)))+cflo)+((1.d0+alh)*dlog
     &           (rnutr))
            dain = dfuint(bet,sa,rr1,rr2,alo)
            if (e2.gt.Oedge) then 
               dain=dain*invdion(ie)
            elseif ((ie.eq.6).AND.(e2.gt.Cedge)) then
               dain=dain*invdion(ie)
            endif
            if (dain.le.0.d0) goto 800
            if (augen.le.0.0d0) then
               rphot(j,ie) = rphot(j,ie)+dain
            else
               auphot(j,ie) = auphot(j,ie)+dain
            end if
c     
c     ***INTEGRATES HEAPH TO BE USED IN SUBR. : PHEAT
c     AND ANR,WNR  FOR THE IONISATION CORRECTION .
c     
            if (ie.gt.limph) goto 800
            satmp = sa+1.d0
            daih = dfuint(bet,satmp,rr1,rr2,alo)
            if (e2.gt.Oedge) then 
               daih=daih*invdion(ie)
            elseif ((ie.eq.6).AND.(e2.gt.Cedge)) then
               daih=daih*invdion(ie)
            endif
            daip = ((daih-dain)*etr)*ev
            if (daip.gt.0.0d0) then
               ed = dmax1(1.d-3,(((rr1*etr)+e2)/2.d0)-(etr+10.2d0))
               heaph(1,j,ie) = heaph(1,j,ie)+(eff1(ed)*daip)
               heaph(2,j,ie) = heaph(2,j,ie)+(eff2(ed)*daip)
               heaph(3,j,ie) = heaph(3,j,ie)+(eff3(ed)*daip)
               heaph(4,j,ie) = heaph(4,j,ie)+(eff4(ed)*daip)
               heaph(5,j,ie) = heaph(5,j,ie)+(eff5(ed)*daip)
               dait = daip-(dain*ryde)
               if (dait.ge.0.0d0) then
                  anr(1,j,ie) = anr(1,j,ie)+dait
                  wnr(1,j,ie) = wnr(1,j,ie)+dain
                  dait = daip-(dain*eau2ev)
                  if (dait.ge.0.0d0) then
                     anr(2,j,ie) = anr(2,j,ie)+dait
                     wnr(2,j,ie) = wnr(2,j,ie)+dain
                  end if
               end if
            end if
c     
            if (augen.gt.ryd) then
               augt = augen-10.2d0
               augerg = augen*ev
               daiau = dain*augerg
               heaph(1,j,ie) = heaph(1,j,ie)+(eff1(augt)*daiau)
               heaph(2,j,ie) = heaph(2,j,ie)+(eff2(augt)*daiau)
               heaph(3,j,ie) = heaph(3,j,ie)+(eff3(augt)*daiau)
               heaph(4,j,ie) = heaph(4,j,ie)+(eff4(augt)*daiau)
               heaph(5,j,ie) = heaph(5,j,ie)+(eff5(augt)*daiau)
               dait = dain*(augerg-ryde)
               if (dait.ge.0.0d0) then
                  anr(1,j,ie) = anr(1,j,ie)+dait
                  wnr(1,j,ie) = wnr(1,j,ie)+dain
                  dait = dain*(augerg-eau2ev)
                  if (dait.ge.0.0d0) then
                     anr(2,j,ie) = anr(2,j,ie)+dait
                     wnr(2,j,ie) = wnr(2,j,ie)+dain
                  end if
               end if
            end if
c     
 800     continue
c     
 1000   continue

       endif
c     
c      call intvec(tphot, q1, q2, q3, q4)
c
c      write(*,*) 'qtosoh q4 ',qtosoh , q4
c
c     add cosmic ray ionization rate, crate, to neutral hydrogen
c     and Helium
c     crate of order 10^-17, cosphi 0-0.75
c
      rphot(1,1) = rphot(1,1)+crate*(1.d0+cosphi)
      rphot(1,2) = rphot(1,2)+crate*(1.d0+cosphi)
c
c      do 600 i = 1,atypes
c         do 600 j = 1,maxion(i)-1
c               write (*,*) 'rphot:',elem(i),rom(j),rphot(j,i)
c 600     continue
c 
      return 
      end
      
