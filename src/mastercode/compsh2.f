cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADY FLOW PLANE PARALLEL RADIATING SHOCK
c	INCLUDES DIFFUSE FIELD , USES CASE B
c	CALL SUBROUTINES  ROOTS,COOL,TIMION,INTVEC,TAUDIST,EQUION,
c	                  WRSPPOP,SUMDATA,AVRDATA,VELSHOCK,ZETAEFF,
c	                  LOCALEM,NEWDIF,TOTPHOT,ZERBUF,PREION
c			  AVERINTO,COPINTO
c
c
c
      subroutine compsh2(luter, dhpr, xhpr, tepr, tepo, hmag, banfil)
c
      include 'cblocks.inc'
c
      double precision mu2,dhpr, xhpr, tepr, tmag,t,de,dh,dr
      double precision cspd,tnloss,press,en,ue,rhotot,wmol,fmui
c
      double precision adj,canor,cavx,cfm,chasl,chate,cofcrit,convg
      double precision d,de2,dedhma,dedhmi,def,dexi,dfa0,dfaa,dfv1,dh2
      double precision distav,dl0,dlf,dli,dlosa,dlosav,dq,dqro
      double precision drdw,drta,drup,dstep,dte,dtema,dtemi
      double precision dv,dva,dvdw,dvi,dvup,eac,egairr
      double precision el0,elf,eli,elosrr,fa,ff,fgae,fi,floe,fmi,fmi0
      double precision fnorm,frdv,frdv0,frdw,frta,frto,ftim
      double precision ga0,gaf,gai,gam,hmag,humag
      double precision prec1,prec2,prec3,pros,qdift,qdl,qf,qs
      double precision rad,rdis,rfvl,rho,rsl
      double precision tdw,tefin,tei,tepo,tl0,tlf,tli
      double precision tlosa0,tlosrr,tsa,tsfin,tsmax,tsrsl
      double precision tstau,tstep,tstep0,tstepi,tsterr
      double precision tup,u,u2,vcons,velav,vshockm,vsi
      double precision wei,xhav,xhy,xhyf
c
      integer*4 i,i107,ic107,iflag,irsl,j,l,m,mi
      integer*4 luop,lupf,lusl,np
      integer*4 luter,luions(4),lsp,nstma2
c
      character imod*4, jmod*4, lmod*4, mmod*4
      character carac*46, pfx*8, tab*4
      character filna*20, filnb*20, filnc*20,filn(4)*20
      character banfil*20, fnam*20,caller*4,wmod*4,ispo*4
c
c     External Functions
c
      double precision fcrit,feldens,fmua,fpressu,frho
c
      jspot = 'NO'
      jcon = 'YES'
      ispo = 'SII'
c
      nstma2 = nstmax-50
c
      if ((iel1 .eq. 5).and.(iel2 .eq. 5)) ispo = ' OI'
c
      do 713 i = limph, atypes
      if (zion(i).gt.0.01) limph = i
  713 continue
c
      dedhmi = 0.0d0
      dedhma = 0.0d0
      do 127 j = 1, atypes
      dedhma = dedhma+zion(j)
      if (arad(2,j).le.0.d00) dedhmi = dedhmi+zion(j)
c
  127 continue
c
      lsp = 26
      luop = 21
      lupf = 27
      lusl = 45
      luter = 6
c
      do i = 1,ieln
         luions(i) = 21+i
      enddo
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
      call s2fheader(fnam,filna,filnb,filnc,filn,
     &luop,lupf,lusl,luions,
     &dhpr, xhpr, tepr, tepo,tmag, banfil)
c
c
c
      do 888 mi = 1, iter
      rad = 1.d38
      wdil = diluf
      fi = 1.0d0
c
      qdift = 0.0d0
      if (mi .eq. 1) goto 921
      dva = 0.0d0
      tei = tepr
      dr = 1.d0
c
c    ***SETTING UP OF PREIONISATION CONDITIONS
c
c
      if (jfin .eq. 'A') tsmax = 1.d2*tsmax
c
      open(luop, file=filna, status='OLD', access='APPEND') 
      write(luop, 677) jpre
c
  677 format(//51h SETTING UP OF THE NEXT PREIONISATION CONDITIONS : 
     &,a1)
c
      lmod = 'UP'
      call localem(tei, depr, dhpr)
c      call contcont(tei, depr, dhpr)
c
c
      call totphot(tei, dhpr, fi, rad, dr, dva, wdil, lmod)
c
      if (jpre .ne. 'U') goto 927
      do 937 l = 1, 4
c
      vsi = vshoc
c
      call preion(luter, luop, tsmax, vshoc, fi, dhpr, depr, tepr,qdift
     &, d)
c
      xhpr = pop(2,1)
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
      if ((dabs(vsi-vshoc)/dmax1(vsi,vshoc)).lt.0.05) goto 927
  937 continue
c
  927 continue
c
      if (jpre .ne. 'F') goto 963
      call copinto(popf, pop)
      xhpr = pop(2,1)
c
      tepr = tefin
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
c
  963 continue
c
      if (jpre .ne. 'T') goto 903
c
      call equion(tepr, depr, dhpr)
      xhpr = pop(2,1)
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
c
  903 continue
c
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
c
c
  921 continue
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
      dhy(1) = dhpr
      dhy(2) = dhpo
      dist(1) = 0.0d0
      dist(2) = 0.0d0
      deel(1) = depr
      deel(2) = depo
      timlps(1) = 0.0d0
c
      qdl = 1.5d0+((0.2d0*(tepo/1.d5))*(nstma2/(loop)))
      fnorm = qdl/(1.0d0-dexp(- qdl))
      vcons =-dabs((vpo-vfin)/(loop-4.0d0))
c
c
c
      qf = qst
      open(luop, file=filna, status='OLD', access='APPEND') 
c
      write(luop, 2023) loop, mi, tepr, qdift
 2023 format(//i4,6h STEPS,t40,12h ITERATION #,i2,10x,7hTEPR : 
     &,f8.0,10x,16hPRESHOCK QTOT : ,1pg10.3/t40,12h -----------/)
c
      write(*, 2027) mi, tepr, qdift, fnam
 2027 format(/8h ITER. #,i2,5x,6hTEPR :,f8.0,1hK,5x,16hPRESHOCK QTOT : 
     &,1pg10.3,5x,a12)
      write(luop, 2047) jpre
 2047 format(1h ,t9,28hINITIAL STATE OF IONISATION ,
     &26hFOR ALL ELEMENTS  (JPRE : ,a1,5h )  :/)
c
      call wionabal(luop,popi)
c
      carac(1:33) = 'I  II IIIIV V  VI '
      if (luter.gt.0) write(luter, 108) 
  108 format(/37h STEP NUMBERS AS THEY ARE COMPLETED :/3h L*,t8,4h<TE>
     &,t18,7h<TLOSS>,t27,6h<DLFR>,t36,5hTSTEP,t45,7h<DENS.>,t55,
     &9h<EL.DENS>,t69,2hDV)
      write(luop, 1220) elem(1), carac(1:3), elem(1), carac(4:7
     &), (elem(iel1), carac(j:j+3),j = 1, 16, 3), (elem(iel2), carac(j
     &:j+3),j = 1, 16, 3)
 1220 format(//37h STEP NUMBERS AS THEY ARE COMPLETED ://3h L*,t7,
     &9h<LOSS.R.>,t18,5hTSTEP,t26,14(1x,a2,a3,1x)/)
c
      close(luop) 
c
c
c	DOMWNSTREAM COMPUTATION ; STEP INDEX : M
c
      drta = 0.0d0
      frta = 0.007d0
      eac = 100.0d0
      canor = 4.0d0
      frdv = 1.0d0
      ff = 1.0d0
      frto = 0.0d0
      fmi = fnorm*dexp(-((qdl*frto)/loop))
      convg = 0.004d0
      do 183 i = 1, 5
  183 cav(i) = canor
      cavx = canor
c
c
      do 1 m = 2, 299
      iflag = 1
      irsl = 0
c
      if (luter.gt.0) write(luter, *) 
c
      qs = qf
      dtema = 0.25
      dtemi = 0.033
      if ((m.lt.(loop/3.0)).and.((((deel(m)/dhy(m))-dedhmi)/
     &dedhma).lt.0.2)) dtemi = dtemi/dabs(dlog(dmax1(((deel(m)/dhy(
     &m))-dedhmi)/dedhma,1.d-10)))
      fmi0 = fmi
      fmi = fnorm*dexp(- ((qdl*frto)/loop))
      dte = dabs(te(m)-te(m-1))/dmax1(te(m-1),te(m))
      frdv = (frdv*(fmi0/fmi))*(((cavx+0.05)/canor) ** (-0.5))
      fmui = 0.97*dsqrt(dtemi/(dte+1.d-9))
      if (dte.lt.(1.05d0*dtemi)) frdv = frdv*dmin1(2.0d0,fmui)
      if (frdv.gt.1.0) frdv = 1.0d0
      dvi = fmi*vcons
c
c    ***DERIVES GUESSED VALUES FOR : DE,DH,TLOSS,DR,TSTEP AND T
c
      if (dabs(dvi/veloc(m)).gt.0.2) dvi =-(veloc(m)*0.2)
  107 continue
c
      ic107 = iflag
c
c      *STEP INDEX M=2
c
      if (m .eq. 2) then
      xh(2) = popi(2,1)
      dh = dhy(2)
      de1 = feldens(dh,popi)
      de = de1
      deel(2) = de
      t = tepo
      tstep = 0.0d0
c
c
      dr = 1.d5
      lmod = 'DW'
      call localem(t, de, dh)
c      call contcont (t,de,dh)
      call totphot(t, dh, fi, rad, dr, dva, wdil, lmod)
      call cool(t, de, dh)
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
c
      dfaa = dabs((fmi*frdv)/dfv1)
      tloss = tl0
      dh = dhy(3)
      de = deel(3)
      t = te(3)
      dlosa = dl0
      te(4) = 0.0
c
      else
c
c      *STEP INDEX M>3
c
      if (((te(m).lt.200.d0).and.(dexi.lt.7.d-2)).or.(dexi.lt.
     &2.d-2)) dtema = 0.4d0
      dfaa = dabs((fmi*frdv)/dfv1)
      adj = dfaa/dfa0
c
      if ((dl0.lt.0.4d0).and.(dli.lt.0.4d0)) adj = 1.0d0
c
      if (irsl .eq. 0) then
      dlosa = dl0+((0.5d0*(dl0-dli))*dfaa)
      if (dabs(dlosa).gt.1.d0) dlosa = dsign(1.d0,dlosa)
      tsterr = (tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36
     &))))*adj
c
      if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.tstau)) tsterr = 
     &tstau
c
      wei = fcrit(dlosa,eac)
      wei = fcrit(dlosa/(1.d0+(10.d0*wei)),eac)
      tstep = dexp(((1.d0-wei)*dlog(tsterr+1.d-30))+(wei*((0.2d0
     &*dlog(tstau))+(0.8d0*dlog(tsterr+1.d-30)))))
      else
      dlosa = dl0+(0.25d0*((dlf-dl0)+((dl0-dli)*dfaa)))
      if (dabs(dlosa).gt.1.d0) dlosa = dsign(1.d0,dlosa)
      tsterr = (tstep0+((tstep0-tstepi)*(tstep0/(tstepi+1.d-36
     &))))*adj
c
      if ((dabs(dlosa).lt.0.7d0).and.(tsterr.gt.tstau)) tsterr = 
     &tstau
c
      wei = fcrit(dlosa,eac)
      wei = fcrit(dlosa/(1.d0+(10.d0*wei)),eac)
      tstep = dexp(((1.d0-wei)*dlog(tsterr+1.d-30))+(wei*dlog(
     &tsa+1.d-30)))
      end if
c
c
c
      tloss = tl0+((0.5d0*(tl0-tli))*dfaa)
      tlf = tloss+(tloss-tl0)
      if ((tli/tl0).gt.0.d0) then
      dh = dhy(m)+(((0.5d0*(dhy(m)-dhy(m-1)))*(dhy(m)/dhy(m-
     &1)))*dfaa)
      t = te(m)+(((0.5d0*(te(m)-te(m-1)))*(te(m)/te(m-1)))
     &*dfaa)
c
      else
c
      dh = dhy(m)-((0.5d0*(dhy(m)-dhy(m-1)))*dfaa)
      t = te(m)-((0.5d0*(te(m)-te(m-1)))*dfaa)
      end if
      t = dmax1((1.d0-(dtema/2.d0))*te(m),dmin1((1.d0+((dtema*(
     &1.d0+dtema))/2.d0))*te(m),t))
      de = (deel(m)*dh)/dhy(m)
      te(m+1) = 0.0d0
c
      end if
c
c    ***DERIVES AVERAGES FOR DH,DE,MU,RHO AND TE,
c	COMPUTES CONDITIONS AT THE END OF TIMESTEP
c
  102 frdv0 = frdv
      dv = dvi*dsign(frdv,tl0)
      dva = dabs(dv)
      veloc(m+1) = veloc(m)+dv
      velav = veloc(m)+(dv/2.d0)
      u = veloc(m+1)
      u2 = u*u
      dstep = ((veloc(m+1)+veloc(m))*tstep)/2.d0
      dr = 0.5d0*dstep
c
c
c
      if (iflag.gt.ic107) goto 1257
      open(lupf, file=filnb, status='OLD', access='APPEND') 
c
      if (luter.gt.0) write(luter, 1230) m, t, tloss, dlosa, 
     &tstep, dh, de, dv
 1230 format(1h ,i3,f10.0,1pg11.3,0pf7.3,4(1pg11.3))
      write(lupf, 1243) m, iflag, ic107, te(m), tl0, tlf, dlosa
     &, dv, tstep, dh, de, frdv, ff
 1243 format(1h ,3(i3),0pf9.0,9x,7(1pg11.3),2(1pg9.2))
      close(lupf) 
c
 1257 continue
c
c
      lmod = 'DW'
      if (iflag .eq. ic107) then 
         call localem(t, de, dh)
c         call contcont(t,de,dh)
      endif
c
      call totphot(t, dh, fi, rad, dr, dva, wdil, lmod)
      call copinto(popi, pop)
      call timion(t, def, dh, xhyf, tstep)
      call copinto(pop, popf)
c
c
c
      xh(m+1) = xhyf
      xhav = (xh(m)+xh(m+1))/2.d0
      humag = (hmago*vpo)/u
      gam = (((pres*u)/fm)-u2)-humag
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
      if (te(m+1).lt.((1.d0-dtema)*te(m))) 
     &   te(m+1)=dmax1(te(m+1),0.15d0*te(m))
      t = (te(m+1)+te(m))/2.d0
      t = dmax1((1.d0-(dtema/1.8d0))*te(m),dmin1((1.d0+((dtema*(
     &1.d0+dtema))/1.8d0))*te(m),t))
      dte = dabs(te(m+1)-te(m))/dmax1(te(m),te(m+1))
      qf = (2.d0*ww)-(((5.d0*gam)+(4.d0*humag))+u2)
      dq = (qf-qs)/2.d0
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
      call averinto(wei, popi, popf, pop)
      call localem(t, de, dh)
c      call contcont(t,de,dh)
c
      if ((dabs(dlosa).lt.0.99d0).or.(dabs(dl0).lt.0.99d0)) then
         call totphot(t, dh, fi, rad, dstep, dva, wdil, lmod)
      endif
c
      call taudist(dh, fi, frta, drta, rad, pop, mmod)
      call copinto(popf, pop)
      call cool(te(m+1), de2, dh2)
      call averinto(wei, popi, popf, pop)
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
      rsl = tl0/tlf
c
      tsfin = dabs((2.d0*dqro)/(tl0+tlf))
      tsrsl = dmin1(tsfin,dabs((2.d0*dqro)/dmax1(dabs(tl0),dabs(tlf))))
c
c      if photon field is very small, then skipbin may be true for all
c     bins, in this case drta will be exactly zero, so tstau will be
c     set equal to the dynamic timescale from the Rankine-Hugenoit
c     eqns.  Note; once the step size increases then the field will 
c     build up and drta should return a non zero result.
c
c      if (drta.gt.0.d0) then 
c         tstau = (drta/((veloc(m)+veloc(m+1))/2.d0))+1.d-20
c      else
         tstau = tsfin
c      endif
c
      iflag = iflag+1
      i107 = iflag-ic107
c
c    ***TEST CONVERGENCE
c
      pros = dmax1(1.d0,((i107*fcrit(dlosa,10.d0))+1.d-4)**0.35d0)
      cofcrit = dmin1(1.001d0,convg+(1.5d0*fcrit(dlosa/pros,3.d1)))
      if ((((rsl.lt.0.d0).and.(dabs(rsl).lt.0.5d0)).and.(dte.gt.
     &dtemi)).and.(irsl.lt.2)) goto 104
      if (dte.gt.dtema) goto 104
      floe = dabs((eloss-elosrr)/dmax1(dabs(eloss),dabs(elosrr)))
      fgae = dabs((egain-egairr)/dmax1(egain,egairr,1.d-33))
      ftim = dabs((tsfin-tsterr)/dmax1(tsfin,tsterr))
      if ((ftim.lt.convg).or.((i107.gt.1).and.(ftim.lt.cofcrit
     &))) goto 110
c
c    ***ITERATES IF NECESSARY ON TIME STEP TO IMPROVE CONVERGENCE
c
  104 continue
      chate = (0.8d0*dtemi)/dte
      if (rsl.gt.0.d0) then
      if (((((rsl.le.3.d0).or.(irsl.ge.2)).or.(dte.le.dtemi))
     &.or.(i107.ge.2)).or.(dte.gt.dtema)) goto 536
      chasl = 0.5d0
      irsl = irsl+1
      frdv = frdv0*dmax1(chate,chasl)
  536 wei = fcrit(dlosa,eac)
      wei = fcrit(dlosa/(1.d0+(1.d1*wei)),eac)
      tsa = dexp(((1.d0-wei)*dlog(tsfin+1.d-30))+((wei*0.5d0)*(
     &dlog(tstau)+dlog(tsterr+1.d-30))))
      tsa = dmax1(tsa,(fcrit(dlosav,1.d0)*tsterr)/1.d1)
      if ((i107.gt.3).and.(dabs(dlosa).gt.0.25d0)) then
      wei = dmin1(0.5d0,1.6d0*((i107/10.d0)**2))*(1.d0-fcrit(dlosa,
     &eac/4.d0))
c
      tsa = ((1.d0-wei)*tsa)+(wei*dmin1(1.45d0*tsa,dmax1(0.7d0*
     &tsa,tsterr)))
c
      end if
      else
      if (((dabs(rsl).ge.0.5d0).or.(irsl.ge.2)).or.(dte.le.dtemi))
     & goto 151
      chasl = dabs(rsl)/5.0d0
      irsl = irsl+1
      frdv = frdv0*dmax1(chate,chasl)
  151 continue
      fa = dabs(tsterr-tsfin)/dmax1(tsterr,tsfin)
      tsa = (tsfin+tsterr)/2.d0
      cfm = 0.5d0-(0.35d0*fcrit(dlosa,2.5d1))
      if (fa.gt.cfm) tsa = tsterr*(1.d0+dsign(cfm,dlog(epsilon
     &+(tsfin/(tsterr+1.d-30)))))
      tsa = dmin1(tstau,dmax1(tsa,tsrsl/2.d0))
c
      end if
c
c
      tstep = tsa
      if (dte.gt.dtema) then
      frdv = frdv0*dmin1(0.6d0*(1.d0-dtema),0.7d0*((dtema/dte)** 
     &1.5d0))
      else if ((i107.gt.7).and.((((deel(m)/dhy(m))-dedhmi)/
     &dedhma).lt.0.2d0)) then
      frdv = frdv0/8.0d0
      else if (i107.gt.7) then
      frdv = frdv0/3.0d0
c
      end if
c
c
c
      open(lupf, file=filnb, status='OLD', access='APPEND') 
c
      if (luter.gt.0) write(luter, 1233) t, tloss, dlosa, 
     &tsterr, tsfin, tstau, rsl
 1233 format(4h    ,f10.0,1pg11.3,0pf7.3,4(1pg11.3))
      write(lupf, 1237) m, iflag, ic107, te(m), te(m+1), tl0
     &, tlf, dlosa, dv, tsterr, tsfin, tstau, tstep, rsl
 1237 format(1h ,3(i3),2(0pf9.0),8(1pg11.3),0pf8.3)
      close(lupf) 
      if (iflag.gt.30) goto 341
      if (frdv .ne. frdv0) goto 107
c
c    ***IF CONVERGENCE UNACHIEVABLE : EXIT GRACEFULLY
c
      goto 102
  341 continue
      write(*, 376) m
  376 format(/' <<<<<<<<<<<< NO CONVERGENCE AT STEP#',i4/)
c
  110 continue
c
c    ***WHEN CONVERGENCE ACHIEVED , PREPARES NEXT TIME STEP
c
      do 197 j = 1, 4
  197 cav(j)=cav(j+1)
      cav(5)=(fcrit(dlosa,8.d0)*canor)+((1.d0-fcrit(dlosa,8.d0))
     &*i107)
      cavx = 0.0d0
      i = 5
      do 191 j = 5, 6-i, -1
  191 cavx = cavx+(j*cav(j))
c
      cavx = cavx/((i+1.d0)*(i/2.d0))
      dist(m+1) = dist(m)+dstep
      timlps(m) = timlps(m-1)+tsterr
      tsmax = 5.d0*timlps(m)
      tefin = te(m+1)
c
c    ***INTEGRATES DIFFUSE FIELD,COMPUTES AVERAGE SPECTRUM AND
c	SUMS UP RELEVANT QUANTITIES
c
      dexi = ((de/dh)-dedhmi)/dedhma
      rdis = (dist(m)+dist(m+1))/2.d0
      tdw = t
      tup = dsqrt(tepo*t)
      drdw = dstep
      dvdw = dva
      drup = dist(m)
      dvup = dsqrt(vpo*u)
      frdw = 0.5d0
c
c
c
      jmod = 'LODW'
      imod = 'ALL'
      call cool(t, de, dh)
      call zetaeff(fi, dh)
      call newdif(tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup, frdw, 
     &jmod)
c
c
c    ***OUTPUT PART OF THE DATA AND TEST ENDING
c
      call sumdata(t, de, dh, fi, drdw, drdw, rdis, imod)
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
        hsii(m) = (fluxf(1,su2)+fluxf(2,su2))/(fluxh(1)+epsilon)
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
      open(luop, file=filna, status='OLD', access='APPEND') 
      open(lupf, file=filnb, status='OLD', access='APPEND') 
      open(lusl, file=filnc, status='OLD', access='APPEND') 
c
      tab = char(9)
      do i = 1 ,ieln 
      open(luions(i),file=filn(i),status='OLD',access='APPEND')
c
 3000 format(I3.3,A1,1pg12.5,A1,29(1pg11.4,A1))
      if (i.eq.1) then
      write(luions(i),3000) m,tab,T,tab,(pop(j,iel1),tab, 
     &j = 1,maxion(iel1))
      endif
      if (i.eq.2) then
      write(luions(i),3000) m,tab,T,tab,(pop(j,iel2),tab,
     &j = 1,maxion(iel2))
      endif
      if (i.eq.3) then
      write(luions(i),3000) m,tab,T,tab,(pop(j,iel3),tab,
     &j = 1,maxion(iel3))
      endif
      if (i.eq.4) then
      write(luions(i),3000) m,tab,T,tab,(pop(j,iel4),tab,
     &j = 1,maxion(iel4))
      endif
c
      close(luions(i))
c
      enddo
c
      if (luter.gt.0) write(luter, 1233) t, tloss, dlosa, 
     &tstep, tsfin, tstau, rsl
      write(luop, 1200) m,tloss,tstep,pop(1,1),pop(2,1),(
     &pop(j,iel1),j = 1, 6), (pop(i,iel2),i = 1, 6)
 1200 format(1h ,i3,1pg10.2,1pg9.2,1x,14(0pf7.4))
      write(lupf, 1237) m, iflag, ic107, te(m), te(m+1),tl0
     &, tlf, dlosa, dv, tsterr, tsfin, tstau, tstep, rsl
      write(lupf, 1247) cavx, caseab(1), dexi, fgae, floe, ftim, 
     &cofcrit, t, dte
 1247 format(7h CAVX :,0pf6.3,3x,10hCASE A,B :,0pf4.2,3x,7hFEXIT :
     &,1pg10.3,3x,4(0pf7.4),1(0pf10.0),1(0pf7.4))
c
c
      press = fpressu(t,dh,pop)
c
      en = zen*dh
c
      ue = 3/2*(en+de)*rkb*t
      tnloss = tloss/(en*de)
c
      rhotot = frho(de,dh)
      cspd = dsqrt(5/3*press/rhotot)
      wmol = rhotot/(en+de)
c
      tab = char(9)
      write(lusl,1210) m,tab,dist(m),tab,dstep,tab,t,tab,de,tab,
     &dh,tab,en,tab,tloss,tab,tnloss,tab,rhotot,tab,press,tab,ue,tab,
     &cspd,tab,wmol,tab
c
 1210 format(i3.3,A1,13(1pg12.5,A1))
c
c
c      wmod = 'FILE'
c      call wmodel(luter,t,de,dh,dstep,wmod)
c
c      call wionabal(luter,pop)
c
      close(luop) 
      close(lupf) 
      close(lusl) 
c
c
c      wmod = 'SCRN'
c      call wmodel(luter,t,de,dh,dstep,wmod)
c
c first upstream integrated field
c
c     diagnostics turned off
c
      if (.false.) then 
      dva = dsqrt(vpo*u)
      lmod = 'UP'
c
      call totphot(t, dh, fi, rad, dstep, dva, wdil, lmod)
c
      caller = 'S2'
      pfx = 'psoup'
      np = 5
      wmod = 'REAL'
c
      call wpsou(caller,pfx,np,wmod,
     & t,de,dh,
     & veloc(m+1),dva,vshoc,
     & dist(m+1),dstep,
     & timlps(m),tsterr,
     & (2.d0/4.d0),tphot)
c
c      pfx = 'psupskp'
c      np = 7
c      sfx = 'sou'
c      ns = 3
c      sf = ''
c      call newfile(pfx,np,sfx,ns,sf)
 100  format(1pf14.7,1x,a)
c      open(lsp, file=sf, status='NEW')
c      do i = 1,infph-1
c         if (skipbin(i)) then
c            write(lsp,100) ephot(i),'1'
c         else
c            write(lsp,100) ephot(i),'0'
c         endif
c      enddo
c      close(lsp)
c
c     then local downstream
c
      dva = 0.0d0
      lmod = 'DW'
c
      call totphot(t, dh, fi, rad, dstep, dva, wdil, lmod)
c
      caller = 'S2'
      pfx = 'psodw'
      np = 5
      wmod = 'REAL'
c
      dv = dvi*dsign(frdv,tl0)
c
      call wpsou(caller,pfx,np,wmod,
     & t,de,dh,
     & veloc(m+1),dv,vshoc,
     & dist(m+1),dstep,
     & timlps(m),tsterr,
     & (2.d0/4.d0),tphot)
c
c      pfx = 'psdwskp'
c      np = 7
c      sfx = 'sou'
c      ns = 3
c      sf = ''
c      call newfile(pfx,np,sfx,ns,sf)
c      open(lsp, file=sf, status='NEW')
c      do i = 1,infph-1
c         if (skipbin(i)) then
c            write(lsp,100) ephot(i),'1'
c         else
c            write(lsp,100) ephot(i),'0'
c         endif
c      enddo
c      close(lsp)

      endif
c
c
      if (iflag.gt.30) goto 788
      if ((t.lt.1500.0).and.(dexi.lt.xhmin)) goto 788
      if ((t.lt.200.0)) goto 788
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
      call copinto(popf, popi)
      call copinto(popf, pop)
c
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
      if (dlosa.gt.0.5) tstep0 = tsfin
      tlosa0 = tloss
      dfv1 = fmi*frdv
      dfa0 = dfaa
      frto = frto+frdv
c
c
c     END MAIN LOOP IN M
c
c
    1 continue
c
      m = m-1
  788 continue
c
c    ***OUTPUT QUANTITIES STORED IN ARRAYS
c
      open(luop, file=filna, status='OLD', access='APPEND') 
      write(luop, 1250) ispo
 1250 format(//16h LIST OF STORED ,13hQUANTITIES : //3h L*,t9,4h<TE>
     &,t16,6h<FHII>,t24,7h<DENS.>,t34,9h<EL.DENS>,t43,6h<DIST>,t53,3h<V>
     &,t63,9h<ELPS.TM>,t74,8hOIII/H-B,t84,7hOII/H-B,t94,7hNII/H-A
     &,t104,a3,4h/H-A/)
      do 1240 l = 2, m
      distav = (dist(l)+dist(l+1))/2.0
      velav = (veloc(l+1)+veloc(l))/2.0
      t = (te(l)+te(l+1))/2.0
      dh = (dhy(l)+dhy(l+1))/2.0
      de = (deel(l)+deel(l+1))/2.0
      xhy = (xh(l)+xh(l+1))/2.0
      write(luop, 1280) l, t, xhy, dh, de, distav, velav, 
     &timlps(l), hoiii(l), hoii(l), hnii(l), hsii(l)
 1280 format(' ',i3,f10.0,f7.4,9(1pg10.3))
 1240 continue
c
c
      close(luop) 
c
c     output upstream field photon source file
c
      dva = dsqrt(vpo*u)
      dr = 1.d0
      lmod = 'UP'
c
      call totphot(tepr, dhpr, fi, rad, dr, dva, wdil, lmod)
c
      caller = 'S2'
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
c     main loop
c
  888 continue
c
      open(luop, file=filna, status='OLD', access='APPEND') 
c
c     formats
c
 2020 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' --------------------',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
c
      if (iter.le.1) goto 531
      write(luop, 533) 
  533 format(///38h INITIAL PARAMETERS CARACTERISING THE ,
     &16hLAST ITERATION :/)

      write(luop, 2020) dhpr, xhpr, tepr, tmag
      write(luop, 2047) jpre
c
      call wionpop(luop,popini)
c
c
 2807 format(/t9,' PHOTOIONISATION SOURCE AT THE SHOCK FRONT :')
 2800 format(/' MOD',t7,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t47,'FQHI',t56,'FQHEI',t66,'FQHEII',t76,'DILUF')
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))
c
      write(luop, 2807) 
      write(luop, 2800) 
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii, diluf
      write(luop, 2813) qdift
 2813 format(/1h ,t10,35h PRESHOCK UPSTREAM IONISING FLUX : ,1pg10.3/)
 2100  format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' ---------------------')
      write(luop, 2100) tepo, vshockm
  531 continue
c
c    ***OUTPUT SPECTRUM AND AVERAGE IONIC POP. AND TEMP.
c
c
      call avrdata
c
c
c
      call wrsppop(luop)
c
      close(luop) 
      write(*, 2005) fnam
c
 2005 format(/30h OUTPUT PRINTED &&&&&& FILE : ,a12/)
      return 
      end
