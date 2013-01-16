cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO EVOLVE IONIC POPULATIONS,TEMPERATURE AND DENSITY WITH
c	TIME (TIME STEP : TIM)
c	THE PHOTON FIELD MUST HAVE BEEN SPECIFIED BEFORE (TOTPHOT)
c	INITIAL AND FINAL IONIC POPUL. IN POP(6,11) /ELABU/
c	JDEN : DENSITY BEHAVIOR ('C' OR 'B')
c	PRESCC : NORMALISED PRESSURE (USED ONLY WHEN JDEN='B')
c		(IT IS USED TO DETERMINE DENSITY BY USING TEMP.)
c
c	NB. EXIT IF RATIO OF EL. DENSITY TO TOTAL DENSITY <=EXL
c	   (SET TO -1.0 WHEN STARTING WITH INITIALLY NEUTRAL GAS)
c	   EXIT ALSO WHEN TEMPERATURE DROPS BELOW  TEX
c	   (TEX SHOULD BE SMALLER THEN TEI)
c	CALL SUBROUTINES COOL,TIMION,COPINTO
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine evoltem(tei, tef, de, dh, prescc, tim, exl, tex, lut)
c
      include 'cblocks.inc'
c
      double precision popj(mxion, mxelem), dleq(4), tmuav(4), cave(4)
      double precision tei, tef, de, dh, prescc, tim, exl, tex,u
      double precision tu, tj, tuf, tuco, dtuj, dtuf, tura, adex
c
      double precision a,adte,ampma
      double precision cavx,cc,ch0,chan,chanp,chxx
      double precision convg,cop0,copeq,cri,cv0
      double precision dchan,dcop,dea,dedhma,dedhmi
      double precision dej,deu,dfco,dfdf,dfm00,dfmr
      double precision dfmu,dfr,dfr00,dfri,dfrpr
      double precision dha,dhj,dhu,div,dlexi
      double precision dljok,dlosj,dlosp,dlosu,dlsk,dlsl,dlto
      double precision dmu,dtim,dtimp,eca,eps
      double precision fnre,froo,ftu,ftu0
      double precision gaj,gamma0,gammap,gammav,gau,gav,gm0,gma,gmi
      double precision prej1,preu1,ra0,radmu
      double precision rain,rdt,rga,rgami,rmuj,rmuu
      double precision rrmu,rxx,ta,tmi,tmin,tmu,tmu0
      double precision tmuma,tmuvar,ttm
      double precision varf,wei,weis,xhy
c
      integer*4 i,ical,ichh,iconsis,idiv,iflip,imu,j,jok
      integer*4 l,lmax,ltt,mma
      integer*4 lut
c
c     External Functions
c
      double precision feldens,fmua,fpressu
c
c
c
c     Internal Functions
c
      double precision far
c
      far(u,a) = (a*(dexp(dmin1(24.0d0,(2.0d0*u)/a))-1.d0))
     &/(dexp(dmin1(24.0d0,(2.0d0*u)/a))+1.0d0)
c
      tmin = 1.0d0
      convg = 0.005d0
      dlexi = 7.d-6
      dfr00 = 0.40d0
      dfm00 = 0.20d0
      ampma = 5.0d0*dfr00
      ltt = 0
      lmax = 12
      mma = 4
      tmuma = 8.0d0
      do 123 j = 1, mma
         cave(j) = 2.0d0
         tmuav(j) = tmuma/3.0d0
 123     dleq(j) = 1.0d0
c
      eps = 0.0d0
      tj = dmax1(dabs(tei),1.1d0*tmin)
      texi = dmax1(tmin,dmin1(tex,0.95d0*tei))
      tef = tj
      if (jden .ne. 'B') then
      dh = dh
      else
      dh = (dh*prescc)/fpressu(tj,dh,pop)
      end if
      de = feldens(dh,pop)
c
      if (tim.le.0.0d0) return 
      dedhmi = 0.0d0
      dedhma = 0.0d0
      do 127 j = 1, atypes
         dedhma = dedhma+zion(j)
         if (arad(2,j).le.0.0d0) dedhmi = dedhmi+zion(j)
 127  continue
c
c
      call copinto(pop, popj)
c
c
      dhj = dh
      dej = de
      prej1 = fpressu(1.0d0,dhj,popj)
      rmuj = fmua(dej,dhj)
c
      call cool(tj, dej, dhj)
c
      dlosj = dlos
      gaj = tloss/prej1
c
      if (jden .eq. 'B') gaj = gaj/tj
      eca = 1.0d0
      dlosp = ((5.0d-3*dlosj)+dsign(1.d-5,dlosj))/2.0d0
      if ((jspot .eq. 'YES').and.(dlosj.gt.0.1)) dlosp = 5.d-3*
     &dsqrt(dlosj)
      wei = 0.975d0
      gammav =-(((2.0d0/3.0d0)*((wei*tloss)
     &+((1.0d0-wei)*eloss)))/(prej1*tj))
c
      if (jden .eq. 'B') gammav = 0.60d0*gammav
      if (lut.gt.0) write(lut, 9) tj,eps,0.0d0, eloss, 
     &floss, dhj, popj(1,1), popj(2,1), rmuj, dlosj
c
      if (lut.gt.6) write(*, 9) tj, eps, 0.0, eloss, floss, 
     &dhj, popj(1,1), popj(2,1), rmuj, dlosj
c
      dfr = dfr00
      dfmr = dfm00/5.0d0
      dfmu = 1.0d0
      tmu = tmuma
      dtimp = 1.d36
      jok = 0
      dljok = 2.d-4
  703 continue
c
c
      ichh = 0
      imu = 0
      idiv = 0
      gamma0 = gammav
      dfmr = dsqrt(dfmr*dfm00)
      dfri = dfr00
      dfr = dmin1(dfri,dsqrt((dfri*dfr)*((dabs(dlosp)+1.d-8) ** eca)
     &))
      dtim = dmin1(dabs(dfr/(gammav+1.d-37)),dmax1(0.0d0,tim-eps))
      dfr = dmin1(dfr,((1.0d-2*tim)/(((eps/50.0)+dtim)+1.d-20)
     &)+(dfr/2.0d0))
  702 dtim = dmin1(dabs((dfr*dfmu)/(gammav+1.d-37)),dmax1(0.0d0,tim
     &-eps),tmu*dtimp)
      dfrpr = (dabs(gammav)*dtim)/dfmu
      dfr = dsqrt(dfr*dmin1(dfr,(4.6d0*tmu)*dfrpr))
  701 continue
      gammap = gammav
      rdt = dtim/dtimp
      tu = tj*dexp((dmin1(gammav*dtim,ampma)))
      varf = (dfr/4.0d0)+(1.1d0*convg)
      ftu0 = 0.0d0
c
      iflip = 0
      do 707 l = 1, lmax
         ta = (tj+tu)/2.d0
         if (jden .ne. 'B') then
            dha = dhj
         else
            dha = (dhj*tj)/ta
         end if
c
         call copinto(popj, pop)
         call timion(ta, dea, dha, xhy, dtim)
c
         if (jden .ne. 'B') then
            dhu = dhj
         else
            dhu = (dh*prescc)/fpressu(tu,dh,pop)
         end if
c
         deu = feldens(dhu,pop)
         rmuu = fmua(deu,dhu)
         preu1 = fpressu(1.0d0,dhu,pop)
         call cool(tu, deu, dhu)
c
         dlosu = dlos
         gau = tloss/preu1
         if (jden .eq. 'B') gau = gau/tu
c
         dmu = ((2.d0/(rmuj+rmuu))*(rmuj-rmuu))/dtim
         if (dabs(gau).ge.dabs(gaj)) then
            gm0 = gaj
            gma = gau
         else
            gm0 = gau
            gma = gaj
         end if
         gmi = (1.d-4*dabs(gma))+epsilon
         rxx = (dej*dhu)/(((deu*dhj)+epsilon)+((1.d-20*dej)*
     &        dhu))
         rga = dabs(gm0)/dmax1(dabs(gma),gmi)
         if (((dlosu.gt.0.9d0).and.(dabs(gau).lt.dabs(gaj))).and.(rxx
     &        .gt.1.0d0)) rga = dmin1(rga*far(rxx ** 0.4d0,10.d0),1.d0)
         copeq = dsign(dmax1(dabs(gma),gmi),gma)
         if ((gm0*dsign(1.d0,gma)).lt.0.d0) then
            cop0 = copeq
            dcop = dabs(cop0-gm0)
         else
            cop0 = copeq-gm0
            dcop = dabs(cop0)
         end if
         weis = 0.50d0
c
         gav = dsign(dexp((weis*dlog(dabs(gma+cop0)))+((1.d0-weis)*
     &        dlog(dabs(gm0+cop0)+gmi)))-dcop,gma)
         tura = tu/tj
         adex = 1.d0
c
         if ((dabs(tura)-1.d0).gt.1.d-3) adex = ((0.5d0*(tura+1.d0)
     &        )/(tura-1.d0))*dlog(tura)
         if (jden .ne. 'B') then
            adte =-((adex*((2.d0/3.d0)*gav))/ta)
            gammav = adte-dmu
         else
            adte =-((2.d0/3.d0)*gav)
            gammav = (3.d0/5.d0)*(adte-dmu)
         end if
c
c
         rrmu = dabs(dmu)/((dabs(adte)+(1.d-10*dabs(dmu)))+epsilon)
         wei = 4.0d0
         radmu = far(rrmu/wei,1.d0) ** wei
         if (radmu.lt.0.05d0) radmu = 0.0d0
         rga = dexp((1.d0-radmu)*dlog(rga+1.d-10))
         dtuj = (tj*dexp((dmin1(gammav*dtim,ampma))))-tj
         tuco = tj+dsign(dmin1(0.8d0*tj,dabs(tu-tj)),dtuj)
         dtuf = (tj*dexp((dmin1(gammav*dtim,ampma))))-tuco
         tuf = tuco+dsign(dmin1(varf*tuco,dabs(dtuf)),dtuf)
         ftu = (tuf-tu)/dmax1(tuf,tu,1.d0)
         cri = dsqrt((dabs(ftu)*dabs(tuf-tu))/(dabs(tuf-tj)+(1.d-5*tj)))
         chan = dabs(dlog(tuf/tj))
         chanp = dabs(dlog(tu/tj))
         iconsis = idint((1.1d0*ftu0)*ftu)/(1.d-20+dabs(ftu0*ftu))
         ftu0 = ftu
         dfdf = dtim*dabs(dmu)
c
c
         if (dfdf.gt.dfmr) then
            imu = imu+1
            if (ltt.gt.6) write(lut, 734) gammav, adte, dmu, 
     &           rga, rdt, dtim, dfr, chan, chanp, cri, dfrpr, 
     &           iconsis, iflip, ichh, imu, idiv
            div = dmax1(1.25d0,dmin1(100.0d0,(dfdf/dfmr) ** 1.5d0))
            dtim = dtim/div
            wei = dmin1(0.5d0,dfr/(chan+1.d-12))
            gammav = (wei*gammav)+((1.0d0-wei)*gamma0)
            goto 701
         end if
c
c
         chxx = chan ** 0.30d0
         rgami = far(chxx-dmin1(0.1d0,0.75d0*chxx,1.0d-3/
     &        (chan+1.d-5)),1.0d0)
         ra0 = rgami/(rga+1.d-6)
         rain = 1.0d0/(ra0+1.d-6)
         if (rga.lt.rgami) then
            idiv = idiv+1
            if (ltt.gt.6) then
               write(lut, 734) gammav, adte, dmu, 
     &              rga, rdt, dtim, dfr, chan, chanp, cri, dfrpr,
     &              iconsis, iflip, ichh, imu, idiv
            end if
            if (idiv .eq. 1) then
               div = dmax1(1.3d0,far(1.1d0*(ra0 ** dmin1(1.d0,0.6d0
     &              +((1.d0/ra0) **0.5d0))),4.d0))
            else if (idiv .eq. 2) then
               div = dmin1(20.d0,dmax1(1.4d0,1.3d0*(ra0**0.65d0)))
            else
               div = dmin1(50.d0,2.d0*(ra0 ** (0.6d0+(rgami/1.5d0))))
            end if
            dtim = dtim/div
            wei = 0.66d0
            gammav = (wei*gammav)+((1.d0-wei)*gammap)
            goto 701
         end if
c
c
         dfco = (dfri*((dfr/dfri) ** 0.8))+3.d-4
         dchan = (chan-dfco)/dmax1(chan,dfco)
         cc = (0.2d0+(0.2d0*((ichh*ichh)/4.d0)))+
     &        (0.2d0*(((l-1)*(l-1))/4.d0))
         cc = far(cc,1.1d0)
         if ((dchan.gt.cc).or.(chan.gt.(1.15*dfri))) then
            ichh = ichh+idint(dsign(1.d0,dchan))
            if ((dabs(dble(ichh)/2.0d0).eq.dint(dabs(dble(ichh)/2.0d0)))
     &           .or.(chan.gt.(1.20d0*dfri))) then
               if (ltt.gt.6) write(lut, 734) gammav, adte,
     &              dmu, rga, rdt, dtim, dfr, chan, chanp, cri, dfrpr,
     &              iconsis, iflip, ichh, imu, idiv
               gammav = (dfmu*gammav)*((chan/dfr) ** 0.20d0)
               dfr = 0.75d0*dfr
               goto 702
            end if     
         end if
c
c
         if (iconsis.lt.0) then
            if (((tu-tj)*(tuf-tj)).lt.0.d0) then
               iflip = iflip+1
               if (iflip.gt.1) then
                  if (ltt.gt.6) write(lut, 734) gammav, adte,
     &                 dmu, rga, rdt, dtim, dfr, chan, chanp, cri, 
     &                 dfrpr, iconsis, iflip, ichh, imu, idiv
                  dfr = dfr/3.d0
                  wei = dmin1(0.5d0,dfr/(chan+1.d-12))
                  gammav = (wei*gammav)+((1.d0-wei)*gamma0)
                  dfmu = dsqrt(dfmu)
                  goto 702
               end if
            end if
            varf = varf/2.0
            tuf = tuco+dsign(dmin1(varf*tuco,dabs(dtuf)),dtuf)
            wei = 0.66
c     
            tuf = (wei*tuf)+((1.0-wei)*tuco)
         else if (iconsis.gt.0) then
            varf = dmin1(dfr/2.0,1.25*varf)
            tuf = tuco+dsign(dmin1(varf*tuco,dabs(dtuf)),dtuf)
            ftu = (tuf-tu)/dmax1(tuf,tu,1.d0)
            cri = dsqrt((dabs(ftu)*dabs(tuf-tu))/(dabs(tuf-tj)+
     &           (1.d-5*tj)))
         end if
c
c
         if (ltt.gt.6) write(lut, 734) gammav, adte, dmu, 
     &        rga, rdt, dtim, dfr, chan, chanp, cri, dfrpr, iconsis,
     &        iflip, ichh, imu, idiv
 734     format(t5,6(1pg9.2),5(0pf7.4),5(i3))
c
c
         ical = l
         chan = dabs(dlog(tuf/tj))
         cv0 = 1.0+(3.0*(dabs(dlosu) ** 0.50))
         ch0 = dmin1(cv0,dmax1(1.2/cv0,dfr/(chan+1.d-36)))
         ftu = (tuf-tu)/dmax1(tuf,tu,1.d0)
         cri = dsqrt((dabs(ftu)*dabs(tuf-tu))/(dabs(tuf-tj)+(1.d-5*tj)))
         tu = tuf
         if ((cri.lt.convg).and.(dchan.lt.cc)) goto 787
 707  continue
c
      l = l-1
      if (lut.gt.0) write(lut, 783) cri
 783  format(' Convergence too sloww in subr. evoltem :',f5.2)
c
 787  continue
      do 370 j = 1, mma-1
         cave(j) = cave(j+1)
 370  continue
      cave(mma)=ical
      wei = 0.0d0
      cavx = 0.0d0
      do 378 j = 1, mma
         wei = wei+j
         cavx = cavx+(j*cave(j))
 378  continue
      cavx = cavx/wei
      dlosp = dmin1(4.d0,dmax1(0.2d0,(5.d-3+dabs(dlosu))/(5.d-3+dabs(
     &     dlosj))))*dlosu
      if ((cavx.lt.1.2d0).and.(rain.gt.10.2d0)) then
         eca = 1.d0-((0.07d0*far(rain/100.d0,1.d0))/(cavx ** 2))
      else if (ical.gt.1) then
         eca = (eca) ** 0.3d0
      end if
c
      eps = eps+dtim
      dhj = dhu
      ta = (ta*tuf)/tj
      tj = tuf
      dej = deu
      rmuj = rmuu
      prej1 = preu1
      dlosj = dlosu
      if ((gaj*dsign(1.0d0,gau)).gt.0.0d0) then
         if (dabs(dlosu).gt.0.02d0) then
            gav = (gau+(gau-((gav*dhu)/dha)))*((dha/dhu) ** 2)
         end if
      end if
      gaj = gau
      dtimp = dtim
      wei = 0.80d0-(0.50d0*(dabs(dlosu) ** 0.50d0))
      dfmu = (wei*ch0)+((1.d0-wei)*dfmu)
      froo =-(dmu/(adte+1.d-36))
      if (jden .ne. 'B') then
         adte =-((adex*((2.d0/3.d0)*gav))/ta)
         gammav = (1.d0+froo)*adte
      else
         adte =-((2.d0/3.d0)*gav)
         gammav = ((3.d0/5.d0)*(1.0+froo))*adte
      end if
c     
      call copinto(pop, popj)
c
      if (lut.gt.0) write(lut, 9) tuf, eps, dtim, eloss
     &     , floss, dhu, pop(1,1), pop(2,1), rmuu, dlosu, dfr
      if (lut.gt.6) write(*, 9) tuf, eps, dtim, eloss, floss
     &     , dhu, pop(1,1), pop(2,1), rmuu, dlosu, dfr
c
      call sumdata(tuf, deu, dhu, 1.0d0, dtim, dtim, dtim, 'REST')
c
      if (jspec.eq.'YES') then
      	  call spec2(lut,'TTWN','ABS')
      endif
c
 9    format(t10,(0pf10.2),7(1pg9.2),3(1pg11.4),0pf4.2)
c
 10   format(' ',0pf10.2,1x,1pg9.2,9(1x,1pg12.5))
c
      do 103 i = 1, mma-1
         tmuav(i) = tmuav(i+1)
         dleq(i) = dleq(i+1)
 103  continue
c     
      tmuav(mma) = dmin1(tmuma,dmax1(0.1d0,rdt))
      dleq(mma) = dlosu*dlosu
      if (dabs(dlosu).gt.6.d-3) then
         dljok = 2.d-4
         jok = 0
      end if
      dlto = 0.0d0
      tmu0 = 0.0d0
      wei = 0.0d0
      do 107 i = 1, mma
         tmu0 = tmu0+(dsqrt(dble(i)) * tmuav(i))
         wei = wei + dsqrt(dble(i))
         dlto = dlto+dleq(i)
 107  continue
      tmu0 = tmu0/wei
      dlto = dsqrt(dlto/mma)
      dlsl = (dleq(mma-2)*dleq(mma-3))/((dleq(mma)*dleq(mma-
     &     1))+1.d-36)
      dlsk = dleq(mma-1)/(dleq(mma)+1.d-36)
      if ((dlto.lt.dljok).and.(dlsl.lt.1.0)) then
         jok = jok+1
         dljok = dmin1(1.4*dljok,dmax1(1.7*dlto,2.d-4))
      end if
c     
      fnre = ((deu/dhu)-dedhmi)/dedhma
      tmuvar = 0.0d0
      wei = 0.0d0
      do 101 i = 1, mma
         wei = wei+dsqrt(dble(i))
         tmuvar = tmuvar+(dsqrt(dble(i))*dabs((tmuav(i)-tmu0)
     &        /tmu0))
 101  continue
      tmuvar = tmuvar/wei
      ttm = 1.0d0
      do 223 i = 2, mma
         ttm = ttm*dmin1(tmuav(i)/dmax1(tmuav(i-1),1.0d0),1.6d0)
 223  continue
      ttm = ttm ** (1.0d0/(mma-1.0d0))
      ttm = dmax1(1.10d0,dmin1(1.65d0,ttm))
      tmu = dmin1(tmuma,tmu0*dmin1(2.0d0,ttm
     &     *((1.0d0+tmuvar)**0.5d0)))
      tmi = 1.4d0
      if ((cavx.lt.1.32d0).and.(rain.gt.5.2d0)) then
         tmi = dmax1(1.4d0,(3.d0*rga)/cavx)
      end if
      tmu = dmax1(tmi,tmu)
      if (((((fnre.le.exl).and.(dlosu.gt.0.d0)).or.((dlto.le.
     &     dlexi).and.(dmin1(dlsl,dlsk).ge.1.0d0))).or.(((jok.gt.3)
     &     .and.(dleq(mma).lt.((2.d0*dlexi)*jok))).and.(rrmu.lt.
     &     0.25d0))).or.(tuf.le.texi)) goto 767
c     
      if (eps.lt.(0.99999d0*tim)) goto 703
 767  continue
      dh = dhu
      de = deu
c     
      tef = dmax1(tuf,tmin)
      if (lut.gt.0) write(lut, 9) tef, tim
      if (lut.gt.6) write(*, 9) tef, tim
c
c
      return
c
      end
