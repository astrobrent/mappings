cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADYFLOWPLANE PARALLEL RADIATING SHOCK
c    TREATES THE ELECTRON GAS AND THE IONS AS A TWO-PHASE MEDIUM
c    INCLUDES DIFFUSE FIELD , USES CASE B
c    CALL SUBROUTINES  ROOTS,COOL,TIMION,INTVEC,TAUDIST,EQUION,
c                      WRSPPOP,SUMDATA,AVRDATA,VELSHOCK,ZETAEFF,
c                      LOCALEM,NEWDIF,TOTPHOT,ZERBUF,PREION
c              AVERINTO,COPINTO
c
      subroutine compsh3(luter, dhpr, xhpr, tepr, tepo, hmag, banfil)
c
      include 'cblocks.inc'
c
      character imod*4, jmod*4, lmod*4, mmod*4
c
      double precision teel(nstmax), teion(nstmax)
      double precision u, usq, qf, qs, dq, upr, gam, qfmo, qsmo
      double precision mu2,eac
c
      double precision adj,ag0,agdiv,agmu,agpro,agq,canor,cavx
      double precision cfm,chasl,chate,cofcrit,continuesa,convg
      double precision d,de,de2,dedhma,dedhmi,def,dexi,dfa0,dfaa
      double precision dfv1,dh,dh2,dhpr,distav,dl0,dlf,dli
      double precision dll0,dlla,dlli,dlosa,dlosav,dqro
      double precision dr,drdw,drta,drup,dstep
      double precision dte,dtee,dtef,dtel0,dtema,dtemi,dteto
      double precision dv,dva,dvdw,dvi,dvup
      double precision egairr,el0,elf,eli,elosrr
      double precision fa,fel0,fel1,felio,felnf,ff,ff0,fgae
      double precision fi,floe,fmi,fmi0,fmui,fnorm,frdv,frdv0
      double precision frdw,frta,frto,ftim,ga0,gaf,gai,hmag,humag
      double precision poppro,prec1,prec2,prec3,pros,qdift,qdl,rad
      double precision rdis,rfvl,rho,rla,rsl,t,taa,tdw,tee
      double precision tefin,tei,teio,tel0,tel1,tel2,tela
      double precision telf,telim,tepo,tepr,tioa,tiof,tl0,tlf,tli
      double precision tlosa0,tlosrr,tmag,tmp107,tsa,tscal,tsfin
      double precision tsmax,tsrsl,tst,tstep,tstep0,tstepi,tsterr
      double precision tta,tup,vcons,velav,vshockm,vsi
      double precision wdt,web,wec,wef,wei,xhav,xhpr,xhyf,zia,zsq
c
      integer*4 i,i107,ic107,icpl,iflag,iii,irsl
      integer*4 j,kk,l,m,mi,mma,ndi
      integer*4 luter,luop,lupf,lusl,lucty,np
      integer*4 nstma2
c
      character filna*20, filnb*20, filnc*20,banfil*20, fnam*20,ispo*4
      character caller*4, wmod*4, pfx*8
c
c     External Functions
c
      double precision fcrit,feldens,felneur,fioncha,fmua,frho
c
c
      nstma2 = nstmax-50
c
      luter = 6
      mma = 300
      jcon = 'YES'
      jspot = 'NO'
      ispo = 'SII'
      limph = 1
c
      do i = limph, atypes
        if (zion(i).gt.1.d-2) limph = i
      enddo
c
      dedhmi = 0.0d0
      dedhma = 0.0d0
      do j = 1, atypes
        dedhma = dedhma+zion(j)
        if (arad(2,j).le.0.d0) dedhmi = dedhmi+zion(j)
      enddo
c
      luop = 21
      lupf = 27
      lusl = 45
      lucty = 5
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
c
c
      call s3fheader(fnam,filna,filnb,filnc,luop,lupf,lusl
     &                 ,dhpr, xhpr, tepr, tepo,tmag, banfil)
c
c    ***SHOCK INITIALISATION FINISHED
c    COMPUTATION STARTS HERE ; ITERATION INDEX : MI
c
      do 888 mi = 1, iter
      rad = 1.d38
      wdil = diluf
      fi = 1.0d0
c
      qdift = 0.0d0
c
      if (mi .eq. 1) goto 921
      dva = 0.0d0
      tei = tepr
      dr = 1.d0
c
c    ***SETTING UP OF PREIONISATION CONDITIONS
c
c
      if (jfin .eq. 'A') tsmax = 1.d2*tsmax
      open(luop, file=filna, status='OLD', access='APPEND') 
      write(luop, 677) jpre
c
  677 format(//' SETTING UP OF THE NEXT PREIONISATION CONDITIONS :' 
     &,a)
      lmod = 'UP'
      call localem(tei, depr, dhpr)
c      call contcont(tei, depr, dhpr)
c
c
      call totphot(tei, dhpr, fi, rad, dr, dva, wdil, lmod)
      if (jpre .ne. 'U') goto 927
      do 937 l = 1, 4
c
      vsi = vshoc
      call preion(lucty, luop, tsmax, vshoc, fi, dhpr, depr, tepr,qdift
     &, d)
      xhpr = pop(2,1)
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
      if ((dabs(vsi-vshoc)/dmax1(vsi,vshoc)).lt.0.05) goto 927
  937 continue
  927 continue
c
      if (jpre .ne. 'F') goto 963
      call copinto(popf, pop)
      xhpr = pop(2,1)
      tepr = tefin
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
c
  963 continue
c
      if (jpre .ne. 'T') goto 903
      call equion(tepr, depr, dhpr)
      xhpr = pop(2,1)
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
c
  903 continue

      if (jpre .ne. 'C') goto 961
c
      call copinto(popini, pop)
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
  961 continue
c
c
      close(luop) 
c
      call zerbuf
  921 continue
c
c
      call copinto(pop, popi)
      call copinto(pop, popf)
      call copinto(pop, popini)
c
      xh(1) = xhpr
      xh(2) = xh(1)
      veloc(1) = vshoc
      veloc(2) = vpo
      vshockm = vshoc/1.d5
      te(1) = tepr
      te(2) = tepo
      teel(1) = te(1)
      teel(2) = teel(1)
      teion(1) = te(1)
      teion(2) = te(2)+(felneur(pop,1)*(te(2)-teel(2)))
      dhy(1) = dhpr
      dhy(2) = dhpo
      dist(1) = 0.0d0
      dist(2) = 0.0d0
      deel(1) = depr
      deel(2) = depo
c
      timlps(1) = 0.0d0
      qdl = 1.5+((0.2*(tepo/1.d5))*(nstma2/(loop)))
      fnorm = qdl/(1.0-dexp(- qdl))
      vcons =-dabs((vpo-vfin)/(loop-4.0))
      upr = vpo
      u = upr
      usq = u*u
      humag = hmago*(vpo/u)
      gam = (((pres*u)/fm)-usq)-humag
      qf = (2.d0*ww)-(((5.d0*gam)+(4.d0*humag))+usq)
      qfmo =-(((5.d0*gam)+(4.d0*humag))+usq)
c
      open(luop, file=filna, status='OLD', access='APPEND') 
      write(luop, 2023) loop, mi, tepr, qdift
 2023 format(//i4,' STEPS',t40,' ITERATION #',i2,10x,'TEPR :' 
     &     ,f8.0,10x,'PRESHOCK QTOT : ',1pg10.3/t40,' -----------'/)
c     write(*, 2027) mi, tepr, qdift, fnam
 2027 format(/' ITER. #',i2,5x,'TEPR :',f8.0,'K',5x,'PRESHOCK QTOT :', 
     &     1pg10.3,5x,a12)
      write(luop, 2047) jpre
 2047 format(1h ,t9,28hINITIAL STATE OF IONISATION ,
     &     26hFOR ALL ELEMENTS  (JPRE : ,a1,5h )  :/)
c
      call wionpop(luop,popi)
c
      write(*, 108) 
 108  format(/' STEP NUMBERS AS THEY ARE COMPLETED :'/' L#',t8,'<Te>'
     &,t18,'<TLOSS>',t27,'<DLFR>',t36,'TSTEP',t45,'<DENS.>',t55
     &,'<EL.DENS>',t69,'DV')
c
      write(luop, 1220) elem(1),rom(1), elem(1),rom(2),
     &(elem(iel1), rom(j),j = 1, 6), (elem(iel2),rom(j),j = 1,6)
 1220 format(//' STEP NUMBERS AS THEY ARE COMPLETED :'//' L#',t7,
     &'<LOSS.R.>',t18,'TSTEP',t26,14(1x,a2,a3,1x)/)
      close(luop) 
c
c
c    DOMWNSTREAM COMPUTATION ; STEP INDEX : M
c
      drta = 0.0d0
      frta = 0.03d0
      eac = 100.d0
      canor = 3.5d0
      frdv = 1.d-2
      ff = 1.0d0
      frto = 0.d0
      fmi = fnorm*dexp(-((qdl*frto)/loop))
      convg = 0.004d0
      do 183 i = 1, 5
  183 cav(i) = canor
      cavx = canor
c
c
c
c
      m = 1
    1 m = m+1
      iflag = 1
      icpl = 0
      eac = dsqrt(50.d0*dmax1(50.d0,eac))
      irsl = 0
c
c    ***IT OVERWRITES STORED QUANTITIES (FROM 100 UP) WHEN NUMBER
c    OF STEPS IS BIGGER THEN 250
c
      if (m .eq. mma) then
      m = 100
      do 23 i = 0, 1
         teel(m-i) = teel(mma-i)
         te(m-i) = te(mma-i)
         teion(m-i) = teion(mma-i)
         dhy(m-i) = dhy(mma-i)
         xh(m-i) = xh(mma-i)
         deel(m-i) = deel(mma-i)
         dist(m-i) = dist(mma-i)
         veloc(m-i) = veloc(mma-i)
         timlps(m-i) = timlps(mma-i)
         hoiii(m-i) = hoiii(mma-i)
         hoii(m-i) = hoii(mma-i)
         hnii(m-i) = hnii(mma-i)
         hsii(m-i) = hsii(mma-i)
   23 continue
      end if
      qs = qf
c
      qsmo = qfmo
      dtema = 0.25d0
      dtemi = 0.033d0
      if ((m.lt.(loop/3)).and.((((deel(m)/dhy(m))-dedhmi)/
     &dedhma).lt.0.2d0)) dtemi = dtemi/dabs(dlog(dmax1(((deel(m)/dhy(
     &m))-dedhmi)/dedhma,1.d-10)))
      fmi0 = fmi
      fmi = fnorm*dexp(- ((qdl*frto)/loop))
      frdv = (frdv*(fmi0/fmi))*(((cavx+0.05)/canor) ** (-0.5))
      if (m.gt.2) then
      if (((te(m).lt.200.0).and.(dexi.lt.7.d-2)).or.(dexi.lt.
     &2.d-2)) dtema = 0.36
      if (wef.gt.0.75) ff = dmin1(1.d0,ff**dmin1(1.d0,0.5d0+((dteto
     &/dtema)*(5.d0**((1.d0/wef)**3.d0)))))
      if (dteto.lt.(1.05*dtemi)) then
      fmui = 0.97*dsqrt(dtemi/(dteto+1.d-9))
      frdv = frdv*dmin1(2.d0,fmui)
      else if ((dte.lt.(2.0*dtemi)).and.(dteto.gt.(2.5*dtemi))
     &) then
      agq = 1.60
      frdv = frdv*dmin1(agq-((agq-1.d0)*((dmin1(dteto,dtema)/
     &dtema)**0.5d0)),dsqrt((2.d0*dtemi)/dte))
      end if
      end if
      frdv = dmin1(frdv,1.d0)
      dvi = fmi*vcons
      if (dabs(dvi/veloc(m)).gt.(0.75d0*dtema))
     &       dvi =-((veloc(m)*0.75d0)*dtema)
  107 continue
c
c    ***DERIVES GUESSED VALUES FOR : DE,DH,TLOSS,DR,TSTEP AND T
c
      ic107 = iflag
c
c      *STEP INDEX M=2
      icpl = icpl+1
      if (m .eq. 2) then
      xh(2) = popi(2,1)
      dh = dhy(2)
      de1 = feldens(dh,popi)
      de = de1
      deel(2) = de
      t = te(2)
      tela = teel(2)
      tstep = 0.0d0
      dr = 1.d5
c
c
c
      lmod = 'DW'
      call localem(tela, de, dh)
c      call contcont(tela, de, dh)
      call totphot(tela, dh, fi, rad, dr, dva, wdil, lmod)
c
c
      call cool(tela, de, dh)
c
c
c
      el0 = eloss
      ga0 = egain
      dl0 = dlos
      tl0 = tloss
c
c      **STEP INDEX M=3
c
      else if (m .eq. 3) then
      dfaa = dabs((fmi*frdv)/dfv1)
      tloss = tl0
      dh = dhy(3)
      de = deel(3)
      t = te(3)+(((0.3*(te(3)-te(2)))*te(3))/te(2))
      tela = teel(3)+(((0.3*(teel(3)-teel(2)))*teel(3))/teel(2))
      tela = dmax1((1.0-(dtema/1.9))*teel(m),dmin1((1.0+((dtema
     &*(1.0+dtema))/1.9))*teel(m),tela))
      dlosa = dl0
c
c      *STEP INDEX M>3
c
      else
      if (dabs(dlf).gt.0.7d0) then
      dfaa = dabs((fmi*frdv)/dfv1)
      adj = dfaa/dfa0
      tsterr = (tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36
     &))))*adj
c
c
      else if ((dabs(dlf).gt.0.3d0).and.(dabs(dl0).gt.0.3d0)) then
      agmu = 5.d0*(1.d0+((3.d0*dabs(dl0)) ** 2.d0))
      dfaa = dabs((fmi*frdv)/dfv1)
      dfaa = dmin1(agmu,dmax1(1.d0/agmu,dfaa))
      adj = dfaa/dfa0
      adj = dmin1(agmu,dmax1(1.d0/agmu,adj))
      tsterr = (tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36
     &))))*adj
      tsterr = dmin1(agmu*tstep0,dmax1(tstep0/agmu,tsterr))
      else
      adj = 1.d0
      dfaa = 1.d0
      agmu = 5.d0
      tsterr = (tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36
     &))))*adj
      tsterr = dmin1(agmu*tstep0,dmax1(tstep0/agmu,tsterr))
      end if
c
c
c
      if (irsl .eq. 0) then
      dlosa = dl0+((0.5d0*(dl0-dli))*dfaa)
      if (dabs(dlosa).gt.1.d0) dlosa = dsign(1.d0,dlosa)
      if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.(ff*tta))) tsterr
     & = ff*tta
      if (iflag.gt.1) tsterr = dmin1(tsterr,1.10*tsa)
      web = fcrit(dlosa,eac)
      web = fcrit(dlosa/(1.d0+(1.d1*web)),eac)
      tstep = dexp(((1.d0-web)*dlog(tsterr+1.d-30))+(web*((0.3d0
     &*dlog(ff*tta))+(0.7d0*dlog(tsterr+1.d-30)))))
      else
      dlosa = dl0+(0.25d0*((dlf-dl0)+((dl0-dli)*dfaa)))
      if (dabs(dlosa).gt.1.d0) dlosa = dsign(1.d0,dlosa)
      if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.(ff*tta))) tsterr
     & = ff*tta
      web = fcrit(dlosa,eac)
      web = fcrit(dlosa/(1.d0+(1.d1*web)),eac)
      tstep = dexp(((1.d0-web)*dlog(dmin1(tsterr,tsa)+1.d-30))+(
     &web*dlog(tsa+1.d-30)))
      end if
c
c
c
      tloss = tl0+((0.5d0*(tl0-tli))*dfaa)
      tlf = tloss+(tloss-tl0)
c
c
c
      if ((tli/tl0).gt.0.d0) then
      dh = dhy(m)+(((0.5d0*(dhy(m)-dhy(m-1)))*(dhy(m)/dhy(m-
     &1)))*dfaa)
      t = te(m)+(((0.5d0*(te(m)-te(m-1)))*(te(m)/te(m-1)))
     &*dfaa)
      tel1 = teel(m)+(((0.5*(teel(m)-teel(m-1)))*(teel(m)/
     &teel(m-1)))*dfaa)
      else
      dh = dhy(m)-((0.5*(dhy(m)-dhy(m-1)))*dfaa)
      t = te(m)-((0.5*(te(m)-te(m-1)))*dfaa)
      tel1 = teel(m)-((0.5*(teel(m)-teel(m-1)))*dfaa)
      end if
c
c
c
      t = dmax1((1.d0-(dtema/2.d0))*te(m),dmin1((1.d0+((dtema*(
     &1.d0+dtema))/2.d0))*te(m),t))
      dlla = dsign(dmin1(dabs(t-teion(m)),dmin1(2.d0,dfaa)*dabs(dll0)),
     &dll0)
      tioa = teion(m)+(0.5d0*dlla)
      tel2 = t+((t-tioa)/(felneur(popi,1)+1.d-10))
c
c
c
      if (((tstep.gt.(12.d0*tscal)).or.((dabs(te(m)-teel(m))/te(
     &m)).lt.1.d-2)).or.((tstep.gt.(2.5d0*tscal)).and.((dabs(t-
     &tel2)/t).lt.0.05d0))) then
      tela = tel1
      else
      wei = dmin1(0.6d0,(0.2d0*tstep)/tscal)
      tela = (wei*tel1)+((1.0-wei)*tel2)
      end if
c
c
c
      tela = dmax1((1.d0-(dtema/1.9d0))*teel(m),dmin1((1.d0+((dtema
     &*(1.d0+dtema))/1.9d0))*teel(m),tela))
      de = (deel(m)*dh)/dhy(m)
      end if
c
c
c    ***DERIVES AVERAGES FOR DH,DE,MU,RHO AND TE,
c    COMPUTES CONDITIONS AT THE END OF TIMESTEP
c
  102 frdv0 = frdv
      ff0 = ff
      dv = dvi*dsign(frdv,tl0)
      dva = dabs(dv)
      u = upr+(dv)
      usq = u*u
      humag = hmago*(vpo/u)
      gam = (((pres*u)/fm)-usq)-humag
      qf = (2.d0*ww)-(((5.d0*gam)+(4.d0*humag))+usq)
      qfmo =-(((5.d0*gam)+(4.d0*humag))+usq)
      dq = (qfmo-qsmo)/2.d0
      veloc(m+1) = u
      velav = veloc(m)+(dv/2.d0)
      dstep = ((veloc(m+1)+veloc(m))*tstep)/2.d0
      dr = 0.5d0*dstep
c
c
c
      if (iflag.gt.ic107) goto 1257
      open(lupf, file=filnb, status='OLD', access='APPEND') 
      write(*, 1230) m, tela, tloss, 
     &dlosa, tstep, dh, de, dv
 1230 format(1h ,i3,f10.0,1pg11.3,0pf7.3,4(1pg11.3))
      write(lupf, 1243) m, iflag, ic107, teel(m), tl0, tlf, 
     &dlosa, dv, tstep, dh, de, frdv, ff
 1243 format(1h ,3(i3),0pf9.0,9x,7(1pg11.3),2(1pg9.2))
      close(lupf) 
 1257 continue
c
c
c
      lmod = 'DW'
      wei = 0.5d0
      if (iflag .eq. ic107) then 
        call localem(tela, de, dh)
c        call contcont(tela,de,dh)
      endif
c
      call totphot(tela, dh, fi, rad, dr, dva, wdil, lmod)
      call copinto(popi, pop)
      call timion(tela, def, dh, xhyf, tstep)
      call copinto(pop, popf)
c
c
c
      call averinto(wei, popi, popf, pop)
c
c
c
      xh(m+1) = xhyf
      xhav = (xh(m)+xh(m+1))/2.d0
      dhy(m+1) = (dhpo*vpo)/veloc(m+1)
      dh1 = dhy(m)
      dh2 = dhy(m+1)
      dh = (dh1+dh2)/2.d0
      de1 = feldens(dh1,popi)
      de2 = feldens(dh2,popf)
      deel(m+1) = de2
      de = (de1+de2)/2.d0
      mu2 = fmua(de2,dh2)
      rho = frho(de,dh)
      te(m+1) = (mu2*gam)/rgas
      te(m+1) = dmax1(0.15d0*te(m),dmin1(te(m)/0.15d0,te(m+1)))
      t = (te(m+1)+te(m))/2.d0
      t = dmax1((1.d0-(dtema/1.8d0))*te(m),dmin1((1.d0+((dtema*(
     &1.d0+dtema))/1.8d0))*te(m),t))
c
c
c
      ndi = 15
      fel0 = felneur(popi,1)
      fel1 = felneur(popf,1)
      felio = felneur(pop,2)
      zsq = fioncha(pop,2) ** 2
      zia = fioncha(pop,1)
      tst = tstep/ndi
      do 12 iii = 1, 10
      telim = teel(m)
      taa = te(m)
      tiof = teion(m)
      do 33 kk = 1, ndi
      rla = ((9.42d0+(1.5d0*dlog(telim)))-(0.5d0*dlog(de+1.d-26)))
     &-dlog(zia+1.d-20)
      tscal = ((3.197d-3*fmua(0.d0,dh))*(((telim*1836.d0)+(tiof/
     &fmua(0.d0,dh)))**1.5d0))/(((((1.d0+(1.d0/felio))*rla)*de)
     &*zsq)+1.d-20)
      dlli = (taa-tiof)-((taa-tiof)*dexp(- (tst/tscal)))
      taa = te(m)+((kk*(te(m+1)-te(m)))/ndi)
      dlli = dsign(dmin1(dabs(taa-tiof),dabs(dlli)),dlli)
      tiof = tiof+dlli
      if (tiof.lt.((1.d0-dtema)*teion(m))) tiof = dmax1(tiof,0.2d0*
     &teion(m))
      felnf = fel0+((kk*(fel1-fel0))/ndi)
      tee = taa+((taa-tiof)/felnf)
      telim = dmax1((1.d0-((4.d0*dtema)/ndi))*telim,dmin1((1.d0+(
     &((4.0*dtema)*(1.d0+dtema))/ndi))*telim,tee))
      telim = dmax1((1.d0-(dtema/0.9d0))*teel(m),dmin1((1.d0+((dtema
     &*(1.d0+dtema))/0.9d0))*teel(m),telim))
   33 continue
c
c
c
      teion(m+1) = tiof
      if (teion(m+1).lt.((1.d0-dtema)*teion(m))) teion(m+1) = 
     &dmax1(teion(m+1),0.2d0*teion(m))
      tel0 = tela
      teel(m+1) = te(m+1)+(((te(m+1))-teion(m+1))/(
     &felneur(popf,1)+1.d-10))
      teel(m+1) = dmax1(0.1d0*teel(m),dmin1(teel(m)/0.1d0,teel(m+1
     &)))
      telf = teel(m+1)
      telf = dmax1((1.d0-(dtema/0.9d0))*teel(m),dmin1((1.d0+((dtema
     &*(1.d0+dtema))/0.9d0))*teel(m),telf))
      tela = (teel(m)+teel(m+1))/2.d0
      tela = dmax1((1.d0-(dtema/1.9d0))*teel(m),dmin1((1.d0+((dtema
     &*(1.0+dtema))/1.9d0))*teel(m),tela))
      if (iii.gt.4) tela = (tela+tel0)/2.d0
      dtel0 = dabs(tela-tel0)/dmax1(tela,tel0)
      if (dtel0.lt.4.d-5) goto 14
   12 continue
   14 continue
c
c
c
      wdt = dmax1(1.d0-(0.4d0*dabs(dlog10(1.d-30+(teel(m)/te(m)))))
     &,dmin1(0.75d0,0.5d0+(iflag-ic107)))
      dte = (dabs(te(m+1)-te(m))/dmax1(te(m),te(m+1)))+1.d-10
      dtef = dabs(teion(m)-teion(m+1))/dmax1(teion(m),teion(m+1))
      dtee = dabs(teel(m+1)-teel(m))/dmax1(teel(m),teel(m+1))
      dteto = dmax1(0.8d0*dte,wdt*dtee,0.25d0*dtef)
c
      dqro = dq*rho
      tlosrr = tloss
      elosrr = eloss
      egairr = egain
      tsterr = tstep
c
c    ***FINDS FINAL AND AVERAGE COOLING RATE
c
      wei = 0.5d0
      lmod = 'DW'
      mmod = 'DIS'
      call localem(tela, de, dh)
c      call contcont(tela,de,dh)
c
      if ((dabs(dlosa).lt.0.99d0).or.(dabs(dl0).lt.0.99d0)) then
      call totphot(tela, dh, fi, rad, dstep, dva, wdil, lmod)
      endif
      call taudist(dh, fi, frta, drta, rad, pop, mmod)
      call copinto(popf, pop)
      call cool(telf, de2, dh2)
      call averinto(wei, popi, popf, pop)
c
c
c
      elf = eloss
      eloss = (el0+elf)/2.d0
      gaf = egain
      egain = (ga0+gaf)/2.d0
      dlf = dlos
      dlosa = (dl0+dlf)/2.d0
      tlf = tloss
      tloss = (tl0+tlf)/2.d0
      rsl = tl0/dsign(dmax1(dabs(tlf),0.01d0*dabs(tl0),1.d-38),tlf)
      wef = fcrit(dlosa,eac)
      wef = fcrit(dlosa/(1.d0+(1.d1*wef)),eac)
c
c
c
      tsfin = dabs((2.d0*dqro)/(tl0+tlf))
      tsrsl = dmin1((tsfin),dabs((2.d0*dqro)/dmax1(dabs(tl0),dabs(
     &tlf))))
      tta = (drta/((veloc(m)+veloc(m+1))/2.d0))+1.d-20
      iflag = iflag+1
      i107 = iflag-ic107
c
c    test convergence
c
      pros = dmax1(1.d0,((i107*fcrit(dlosa,1.d1))+1.d-4)**0.35d0)
      cofcrit = dmin1(1.001d0,convg+(1.5d0*fcrit(dlosa/pros,3.d1)))
      if ((((rsl.lt.0.0d0).and.(dabs(rsl).lt.0.5d0)).and.(irsl.lt.2
     &)).and.(dteto.gt.((1.d0+((0.1d0*(irsl+1))*(irsl+1)))*
     &dtemi))) goto 104
c
c
c
      if (dteto.gt.dtema) goto 104
      floe = dabs((eloss-elosrr)/dmax1(dabs(eloss),dabs(elosrr)))
      fgae = dabs((egain-egairr)/dmax1(egain,egairr,1.d-33))
      ftim = dabs((tsfin-tsterr)/dmax1(tsfin,tsterr))
c
c    ***ITERATES IF NECESSARY ON TIME STEP TO IMPROVE CONVERGENCE
c
      if ((ftim.lt.convg).or.((i107.gt.1).and.(ftim.lt.cofcrit
     &))) goto 110
  104 continue
      tmp107 = dble(i107)
      eac = dmin1(2.d3,eac*dmax1(1.0d0,2.0d0*dsqrt(tmp107/8.d0)))
      chate = (0.7d0*dtemi)/dteto
      if (rsl.gt.0.d0) then
      if (((((rsl.le.3.0d0).or.(irsl.ge.2)).or.(dteto.le.dtemi))
     &.or.(i107.ge.2)).or.(dteto.gt.dtema)) goto 536
      chasl = 0.5d0
      irsl = irsl+1
      frdv = frdv0*dmax1(chate,chasl)
  536 continue
c
c
c
     tsa = dexp(((1.d0-wef)*dlog(tsfin+1.d-30))+((wef*0.5d0)*(
     &dlog(ff*tta)+dlog(dmin1(tsterr,ff*tta)+1.d-30))))
      agdiv = dmax1(10.d0,((tsterr/(tsfin+1.d-30))+0.01d0)**(0.25
     &d0+(0.025d0*dble(i107))))
      tsa = dmax1(tsa,(fcrit(dlosav,1.d0)*tsterr)/agdiv)
      if ((i107.gt.4).and.(dabs(dlosa).gt.0.35d0)) then
      wei = dmin1(0.5d0,1.d0*((i107/1.d1)**2))*(1.d0-fcrit(dlosa,
     &1.d0))
      tsa = ((1.d0-wei)*tsa)+(wei*dmin1(1.45d0*tsa,dmax1(0.7d0*
     &tsa,tsterr)))
      end if
c
c
c
      else
      if (((dabs(rsl).ge.0.5d0).or.(irsl.ge.2)).or.(dteto.le.
     &dtemi)) goto 151
      chasl = dabs(rsl)/5.d0
      irsl = irsl+1
      frdv = frdv0*dmax1(chate,chasl)
  151 continue
c
c
c
      fa = dabs(tsterr-tsfin)/dmax1(tsterr,tsfin)
      tsa = (tsfin+tsterr)/2.d0
      cfm = 0.5d0-(0.35d0*fcrit(dlosa,2.5d1))
      if (fa.gt.cfm) tsa = tsterr*(1.0d0+dsign(cfm,dlog(1.d-38+(
     &tsfin/(tsterr+1.d-30)))))
      tsa = dmin1(ff*tta,dmax1(tsa,tsrsl/2.0))
      end if
c
c
c
      tstep = tsa
      open(lupf, file=filnb, status='OLD', access='APPEND') 
      write(*, 1233) tela, tloss, dlosa, 
     &tsterr, tsfin, tta, rsl
 1233 format(4h    ,f10.0,1pg11.3,0pf7.3,4(1pg11.3))
      write(lupf, 1237) m, iflag, ic107, teel(m), teel(m+1), 
     &tl0, tlf, dlosa, dv, tsterr, tsfin, tta, tstep, rsl
 1237 format(1h ,3(i3),2(0pf9.0),8(1pg11.3),0pf8.3)
      close(lupf) 
c
c
c
      ag0 = 1.d0
      if ((dteto.gt.dtema).and.((i107.gt.1).or.(dabs(dlosa).gt.
     &0.3d0))) then
      ag0 = dmin1(0.6d0*(1.d0-dtema),0.7d0*((dtema/dteto) ** 2.d0))
      eac = dmin1(2.d3,1.5*eac)
      else if ((i107.gt.7).and.((((deel(m)/dhy(m))-dedhmi)/
     &dedhma).lt.0.2d0)) then
      ag0 = 0.125d0
      else if (i107.gt.7) then
      ag0 = 0.3d0
      end if
c
c
c
      if (ag0.lt.1.d0) then
      wec = wef
      frdv = frdv*dexp((1.0-wec)*dlog(ag0))
      ff = ff*dexp(wec*dlog(ag0))
      end if
c
c
c
      if (iflag.gt.30) goto 341
      if ((frdv .ne. frdv0).or.(ff .ne. ff0)) goto 107
      goto 102
c
c    ***IF CONVERGENCE UNACHIEVABLE : EXIT GRACEFULLY
c
  341 write(*, 376) m
  376 format(/' <<<<<<<<<<<< NO CONVERGENCE AT STEP#',i4/)
  110 continue
c
c    ***WHEN CONVERGENCE ACHIEVED , PREPARES NEXT TIME STEP
c
      do 197 j = 1, 4
  197 cav(j) = cav(j+1)
      agpro = i107+dmin1(1.5d0,icpl-1.d0)
      cav(5) = (fcrit(dlosa,8.d0)*canor)+((1.d0-fcrit(dlosa,8.d0))
     &*agpro)
      cavx = 0.d0
      i = 5
      do 191 j = 5, 6-i, -1
  191 cavx = cavx+(j*cav(j))
      cavx = cavx/((dble(i)+1.d0)*(dble(i)/2.d0))
c
c
c
      dist(m+1) = dist(m)+dstep
      timlps(m) = timlps(m-1)+tsterr
      tsmax = 5.d0*timlps(m)
      tefin = te(m+1)
c
c    ***INTEGRATES DIFFUSE FIELD,COMPUTES AVERAGE SPECTRUM AND
c    SUMS UP RELEVANT QUANTITIES
c
      dexi = ((de/dh)-dedhmi)/dedhma
      rdis = (dist(m)+dist(m+1))/2.d0
      tdw = tela
      tup = ((tela*te(2))*teel(2)) ** (1.d0/3.d0)
      drdw = dstep
      dvdw = dva
      drup = dist(m)
      dvup = dsqrt(vpo*u)
      frdw = 0.5d0
      jmod = 'LODW'
      imod = 'ALL'
c
c
c
      call cool(tela, de, dh)
      call zetaeff(fi, dh)
      call newdif(tdw,tup,dh,fi,rad,drdw,dvdw,drup,dvup,frdw,jmod)
c
c
c    ***OUTPUT PART OF THE DATA AND TEST ENDING
c
      call sumdata(tela, de, dh, fi, drdw, drdw, rdis, imod)
c
      if (ox3.ne.0) then 
        hoiii(m) = (fluxf6(7,ox3)+fluxf6(10,ox3))/(fluxh(2)+epsilon)
      endif
      if (ox2.ne.0) then 
        hoii(m) = (fluxf(1,ox2)+fluxf(2,ox2))/(fluxh(2)+epsilon)
      endif
      if (ni2.ne.0) then 
        hnii(m) = (fluxf6(7,ni2)+fluxf6(10,ni2))/(fluxh(2)+epsilon)
      endif
      if (su2.ne.0) then 
        hsii(m) = (fluxf(1,su2)+fluxf(2,su2))/(fluxh(2)+epsilon)
      endif
      if (ox1.ne.0) then 
      if (ispo .eq. ' OI') hsii(m) = fluxf(3,ox1)/(fluxh(1)+epsilon)
      endif
      if (jfin .ne. 'E') goto 739
      if (istf .eq. 1) rfvl = hoiii(m)
      if (istf .eq. 2) rfvl = hoii(m)
      if (istf .eq. 3) rfvl = hnii(m)
      if (istf .eq. 4) rfvl = hsii(m)
  739 continue
c
c
c
      open(luop, file=filna, status='OLD', access='APPEND') 
      open(lupf, file=filnb, status='OLD', access='APPEND') 
      open(lusl, file=filnc, status='OLD', access='APPEND') 
c
c
c
      write(*, 1233) tela, tloss, dlosa, 
     &tstep, tsfin, tta, rsl
      write(luop, 1200) m, tloss, tstep, pop(1,1), pop(2,1), (
     &pop(j,iel1),j = 1, 6), (pop(i,iel2),i = 1, 6)
 1200 format(1h ,i3,1pg10.2,1pg9.2,1x,14(0pf7.4))
      write(lupf, 1237) m, iflag, ic107, teel(m), teel(m+1), 
     &tl0, tlf, dlosa, dv, tsterr, tsfin, tta, tstep, rsl
      write(lupf, 1247) cavx, caseab(1), dexi, fgae, floe, ftim, 
     &cofcrit, t, tela, dteto, wef
 1247 format(7h CAVX :,0pf6.3,3x,10hCASE A,B :,0pf4.2,3x,7hFEXIT :
     &,1pg10.3,3x,4(0pf7.4),2(0pf10.0),2(0pf7.4))
c
      write(lusl,*)'It#   distance         dr            Te    ',
     &'   Loss rate        de            dh   '   
      write(lusl,1210) m,dist(m),dstep,tela,tloss,de,dh
 1210 format(i3.3,4(1pg14.7))
      call wionabal(lusl,pop)
c
      close(luop) 
      close(lupf) 
      close(lusl) 
c
c
c
      if (iflag.gt.30) goto 788
      if ((t.lt.3000.0d0).and.(dexi.lt.xhmin)) goto 788
      if ((jfin .eq. 'B').and.(te(m+1).le.texi)) goto 788
      if (((jfin .eq. 'C').and.(m.gt.5)).and.(fval.lt.timlps(m))
     &) goto 788
      if (((((jfin .eq. 'D').and.(m.gt.5)).and.(pop(istf+1,iel1)
     &.lt.prec2)).and.(pop(istf,iel1).lt.prec1)).and.(pop(istf,
     &iel1).lt.fval)) goto 788
      if ((((jfin .eq. 'E').and.(m.gt.(loop/5.0))).and.(rfvl
     &.lt.fval)).and.(prec3.gt.fval)) goto 788
      prec2 = pop(istf+1,iel1)
      prec1 = pop(istf,iel1)
      prec3 = rfvl
c
c
c
      call copinto(popf, popi)
      call copinto(popf, pop)
c
c
c
      upr = u
      tli = tl0
      tl0 = tlf
      eli = el0
      el0 = elf
      gai = ga0
      ga0 = gaf
      dli = dl0
      dl0 = dlf
      tstepi = tstep0
      tstep0 = tsterr
      if (dlosa.gt.0.5d0) tstep0 = tsfin
      tlosa0 = tloss
      dfv1 = fmi*frdv
      dfa0 = dfaa
      dll0 = dlf
      frto = frto+frdv
      goto 1
c
c    ***OUTPUT QUANTITIES STORED IN ARRAYS
c
  788 continue
c
c
c
      if (mi .eq. iter) then
      open(luop, file=filna, status='OLD', access='APPEND') 
      write(luop, 1250) 
 1250 format(//16h LIST OF STORED ,13hQUANTITIES : //3h L*,t10,4h<TE>
     &,t18,6h<TEEL>,t27,7h<TEION>,t37,7h<DENS.>,t46,9h<EL.DENS>,t56,
     &6h<DIST>,t66,3h<V>,t73,9h<ELPS.TM>,t83,8hOIII/H-B,t92,7hOII/H-B
     &,t101,7hNII/H-A,t110,7hSII/H-A/)
c
c
c
      do 1240 l = 2, m
      distav = (dist(l)+dist(l+1))/2.0
      velav = (veloc(l+1)+veloc(l))/2.0
      t = (te(l)+te(l+1))/2.0
      tela = (teel(l)+teel(l+1))/2.0
      teio = (teion(l)+teion(l+1))/2.0
      dh = (dhy(l)+dhy(l+1))/2.0
      de = (deel(l)+deel(l+1))/2.0
      write(luop, 1280) l, t, tela, teio, dh, de, distav, velav
     &, timlps(l), hoiii(l), hoii(l), hnii(l), hsii(l)
 1280 format(1h ,i3,3(0pf10.0),2(1pg10.3),7(1pg9.2))
 1240 continue
      close(luop) 
c
c
c
      end if
c
c
c
      dva = dsqrt(vpo*u)
      dr = 1.d0
      lmod = 'UP'
c
      call totphot(tepr, dhpr, fi, rad, dr, dva, wdil, lmod)
c
      caller = 'S3'
      pfx = 'psoup'
      np = 5
      wmod = 'REAL'
c
      call wpsou(caller,pfx,np,wmod,
     & tepr,depr,dhpr,
     & veloc(m),vpo,vshoc,
     & dist(m),dstep,
     & timlps(m),tsterr,
     & (2.d0/4.d0),tphot)
c
c***    MODEL TERMINATED
c
  888 continue
c
c
c
      open(luop, file=filna, status='OLD', access='APPEND') 
c
c
c
      if (iter.le.1) goto 531
      write(luop, 533) 
  533 format(///' INITIAL PARAMETERS CHARACTERISING THE ',
     &'LAST ITERATION :'/)
c
c
c
 2020 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' --------------------',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write(luop, 2043) zgas
 2043  format(1h ,t9,'Abundances of the elements relative',
     &' to hydrogen :','  (Zgas=',f7.4,' Zsun)'/
     &' and the proto-ionisation used.')
c
      call wionabal(luop,poppro)
c
 2807  format(/t9,' PHOTOIONISATION SOURCE AT THE SHOCK FRONT :')
      write(luop, 2800)
 2800  format(/' MOD',t7,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t47,'FQHI',t56,'FQHEI',t66,'FQHEII',t76,'DILUF')
c
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii, diluf
c
 2850  format(1h ,a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))
c
      vshockm = vshoc/1.d5
      write(luop, 2100) tepo, vshockm
c
 2100  format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' ---------------------')
c
      write(luop, 2020) dhpr, xhpr, tepr, tmag
      write(luop, 2047) jpre
c
      call wionpop(luop,popini)
c
      write(luop, 2807) 
      write(luop, 2800) 
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii, diluf
      write(luop, 2813) qdift
 2813 format(/1h ,t10,35h PRESHOCK UPSTREAM IONISING FLUX : ,1pg10.3/)
      write(luop, 2100) tepo, vshockm
c
c    ***OUTPUT SPECTRUM AND AVERAGE IONIC POP. AND TEMP.
c
c
c
  531 continue
c
c
c
      call avrdata
c
      call wrsppop(luop)
      close(luop) 
c
c
c
      write(*, 2005) fnam
 2005 format(/' OUTPUT CREATED IN &&&&&& FILE : ',a/)
c
c
c
      return 
      end
