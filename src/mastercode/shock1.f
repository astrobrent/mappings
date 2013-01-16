cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADY FLOW PLANE PARALLEL RADIATING SHOCK
c    CALL SUBROUTINES  POPCHA,ROOTS,COOL,TIMION,
c                      WRSPPOP,SUMDATA,AVRDATA
c                      VELSHOCK,TOTPHOT
c
      subroutine shock1()
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision mu2, t, de, dh, epotmi
      double precision de2,dh2,dhpr,distav,dq,dr,dstep,dte
      double precision dtemp1,dtemp2,dtemp3,dv,fcool,fi,fnorm,frdv
      double precision frdv0,ftim,gam,hmag,humag
      double precision prec1,prec2,prec3,qdl,qf,qs,rdis,rfvl,rho
      double precision rho2,tepo,tepr,tlosa0,tlosa1,tlosav,tlosrr,tmag
      double precision tstep,tstep0,tstep1,tstepa,tstepf,u,u2,vcons
      double precision velav,vshockm,xhav,xhpr,xhy,xhyf
c
      integer*4 iflag,l
      integer*4 luop,nstma2,ibell
      integer*4 i,j,m,ifil
c
      character imod*4, lmod*4
      character carac*46, filna*20,ispo*4
      character ilgg*4, caract*4,where*16
c
      logical iexi
c
c           Functions
c
      double precision feldens,fmua,frho
c
      nstma2 = nstmax-50
c
      luop = 21
      ispo = 'SII'
c
 1400 write(*, 1130) 
 1130 format(///'SHOCK MODEL SELECTED : S1'/
     &'PHOTOIONISATION NOT TAKEN INTO ACCOUNT')
c
      ibell = 117901063
      jcon = 'NO'
      jspot = 'YES'
      where = 'Shock 1'
      epotmi = ryde
      limph = 2
      do 713 i = limph, atypes
      if (zion(i).gt.0.01d0) limph = i
  713 continue
c
c    ***SET UP INITIAL IONISATION CONDITIONS
c
      where = 'pre-ionisation'
      call popcha(where)
      where = 'Shock 1'
c
      carac(1:33) = '   '
      do 1120 i = 3, atypes
      j = 1+index(carac(1:33),'   ')
      caract = '   '
 1010 format(a2)
 1011 format(a)
      if ((epot(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) then
           write(caract, 1010) elem(i)
      endif
      carac(j:j+2) = caract // '   '
 1120 continue
      i = j+2
 1110 write(*, 1112) carac(1:i)
 1112 format(//' Do you allow the elements : ',a/
     &' to recombine freely to neutral state (Y/N) ?')
      read(*, 1011, err=1110) ilgg
      ilgg = ilgg(1:1)
      ilgg = ilgg(1:1)
      if (ilgg .eq. 'y') ilgg = 'Y'
      if (ilgg .eq. 'n') ilgg = 'N'
      if (ilgg .eq. 'Y') goto 1116
      if (ilgg .ne. 'N') goto 1110
      do 1118 i = 3, atypes
      arad(2,i) = dabs(arad(2,i))
      if ((epot(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) arad(2,i)
     & =-arad(2,i)
 1118 continue
 1116 continue
c
c    ***NOW INPUT OF PHYSICAL CONDITIONS
c
 1018 write(*, 1020) 
 1020 format(//' Select shock conditions :')
 1025 write(*, 1030) 
 1030 format(//' Give initial hydrogen density and ionised fraction:'/)
      read(*, *, err=1025) dhpr, xhpr
      if (((xhpr.lt.0.0d0).or.(dhpr.le.0.0d0)).or.(xhpr.gt.1.d0)) 
     &goto 1025
c
 1042 write(*, 1040) 
 1040 format(//' Give initial and post-shock temperatures : ')
      read(*, *, err=1042) tepr, tepo
      if (tepr.le.10) tepr = 10**tepr
      if (tepo.le.10) tepo = 10**tepo
c
      if ((tepo.le.tepr).or.(tepr.le.100.0d0)) then
      write(*, 29) ibell
   29 format(50h UNSUITABLE ANSWER(S) ****************************,a1)
      goto 1042
      endif
c
 1047 write(*, 1051) 
 1051 format(//' Give magnetic field (micro-gauss) : ')
      read(*, *, err=1047) tmag
      if (tmag.lt.0.0d0) goto 1047
      hmag = tmag*1.0d-6
c
 1300 write(*, 1310) nstma2
 1310 format(//' GIVE NUMBER OF STEPS (MAX:',i4,') '/
     &' AND CHOOSE TWO ELEMENTS TO BE FOLLOWED  :'/
     &' ( #steps  Z1  Z2):' )
      read(*, *, err=1300) loop, iel1, iel2
c
      if (iel1.gt.28) goto 1300
      if (iel2.gt.28) goto 1300
c      
c
      if (iel1.le.0) goto 1300
      if (iel2.le.0) goto 1300
c
      iel1 = zmap(iel1)
      iel2 = zmap(iel2)
c      
      if (iel1.eq.0) goto 1300
      if (iel2.eq.0) goto 1300
c
      if ((loop.gt.nstma2).or.(loop.lt.8)) goto 1300
c
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
c
      vshockm = vshoc/1.d5
      qdl = 1.5d0+((0.2d0*(tepo/1.d5))*(nstma2/(loop)))
      fnorm = qdl/(1.0d0-dexp(-qdl))
c
c    ***SET TYPE OF ENDING FOR THE PROGRAM
c
      xhmin = 0.01d0
      write(*, 1800) vshockm, qdl
 1800 format(//' SHOCK VELOCITY :',f8.1,' KM/SEC',7x,'Q FINE-TUNING :'
     &,f7.2)
 1810 write(*, 1820) xhmin, elem(iel1)
 1820 format(//' CHOOSE TYPE OF ENDING FOR THE PROGRAM :'/t6,
     &' A-  NORMAL ENDING  (FINAL H+/H :',2pf4.1,'% )'/t6,
     &' B-  ACCORDING TO FINAL TEMPERATURE'/t6,
     &' C-  ACCORDING TO ELAPSED TIME'/t6,
     &' D-  ACCORDING TO THE IONISED FRACTION OF : ',a2/t6,
     &' E-  ACCORDING TO LINE INTENSITIES'/t6,' R-  RE-INITIALIZE'/' ' 
     &,t55,':: ',$)
      read(*, 1821, err=1810) jfin
 1821 format(a)
      jfin = jfin(1:1)
c
      if ((jfin.lt.'A').or.((jfin.gt.'E').and.(jfin .ne. 'R'))) 
     &goto 1810
      if (jfin .eq. 'R') return 
      if (jfin .eq. 'A') goto 1895
      if (jfin .ne. 'D') goto 1877
c
 1830 write(*, 1840) 
 1840 format(/' WHICH IONISING STAGE :' )
      read(*, *, err=1810) istf
      if ((istf.lt.2).or.(istf.gt.maxion(iel1))) goto 1810
c
 1877 continue
c
      if (jfin .ne. 'E') goto 1838
 1863 write(*, 1868) 
 1868 format(/' CHOOSE LINE NUMBER : ' ,
     &' [OIII]   2 [OII]   3 [NII]   4 [SII]  ::')
      read(*, *, err=1810) istf
      if ((istf.lt.1).or.(istf.gt.4)) goto 1810
 1838 continue
c
 1823 write(*, 1825) 
 1825 format(/' GIVE FINAL VALUE : ')
      read(*, *, err=1810) fval
      if (fval.lt.0.0d0) goto 1810
c
      if (jfin .eq. 'B') texi = fval
 1895 continue
c
      do 912 i = 1, atypes
      do 911 j = 1, maxion(i)
      popi(j,i) = pop(j,i)
  911 continue
  912 continue
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
      xh(1) = xhpr
      xh(2) = xh(1)
      veloc(1) = vshoc
      veloc(2) = vpo
      te(1) = tepr
      te(2) = tepo
      dhy(1) = dhpr
      dhy(2) = dhpo
      dist(1) = 0.0d0
      dist(2) = 0.0d0
      deel(1) = depr
      deel(2) = depo
      timlps(1) = 0.0d0
      vcons =-dabs((vpo-vfin)/(dble(loop)-4.0d0))
      qf = qst
c
c
c    choose new file name SHCKN****.SH1
c
c
      ifil = 0
 1700 ifil = ifil+1
      filna(1:17) = ' '
      write(filna,1710) ifil
 1710 format(i4.4)
      filna = 'shckn'//filna
      filna(10:17)='.sh1'
      inquire(file=filna, exist=iexi)
      if (iexi) goto 1700
      open(luop, file=filna, status='NEW')
c
c
      write(luop, 2000) 
 2000 format(' ',t40,'PLANE PARALLEL SHOCK MODEL  (NO PHOTOIONISATION)'
     &/' ',t40,'--------------------------'//)
      write(luop, 2020) dhpr, xhpr, tepr, tmag
 2020 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' --------------------',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write(luop, 2043) zgas
 2043 format(' ',t9,'ABUNDANCES OF THE ELEMENTS RELATIVE',
     &' TO HYDROGEN AND THEIR INITIAL STATE OF IONISATION : (ZGAS=',
     &f7.4,' ZSUN )')
c
      call wionabal(luop,popi)
c
      write(luop, 2100) tepo, vshockm
 2100 format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' ---------------------')
      write(*, 108) 
  108 format(//' STEP NUMBERS AS THEY ARE COMPLETED :'/' L#',t8,'<TE>'
     &,t20,'<TLOSS>',t32,'TSTEP',t44,'<DENS.>',t56,'<EL.DENS>',t68,'<V>'
     &/)
      write(luop, 1220) elem(1), rom(1), elem(1), rom(2)
     &, (elem(iel1),rom(j),j=1,6),(elem(iel2),rom(j),j=1,6)
 1220 format(//' STEP NUMBERS AS THEY ARE COMPLETED :'//' # ',t7,
     &'<LOSS.R.>',t18,'TSTEP',t26,14(1x,a2,a3,1x)/)
c
c    ***SHOCK INITIALISATION FINISHED
c    DOWNSTREAM COMPUTATION STARTS HERE , STEP INDEX : M
c
c
      frdv = 1.0d0
      do 1 m = 2, 249
c
c    ***GUESS ON THE AVERAGE COOLING RATE
c
      iflag = 0
      if (m.gt.2) goto 77
      dh = dhy(m)
      de1 = feldens(dh,popi)
      de = de1
      t = te(m)
      tstep = 0.0d0
      call cool(t, de, dh)
   77 continue
      qs = qf
      if (m.le.3) goto 102
      tloss = tlosav+((tlosav-tlosa0)*(tlosav/(tlosa0+epsilon)))
      if ((tloss*dsign(1.0d0,tlosav)).lt.(0.3d0*dabs(tlosav))) tloss = 
     &0.3d0*tlosav
      tstep = tstepf+((tstepf-tstep0)*(tstepf/(tstep0+epsilon)))
      dh = dhy(m)+((0.5d0*(dhy(m)-dhy(m-1)))*(dhy(m)/dhy(m-1)))
      t = te(m)+((0.5d0*(te(m)-te(m-1)))*(te(m)/te(m-1)))
  102 continue
c
c    ***DERIVES AVERAGES FOR DH,DE,MU,RHO AND TE,
c    COMPUTES CONDITIONS AT THE END OF TIMESTEP
c
      call timion(t, de, dh, xhyf, tstep)
      frdv0 = frdv
      dv = vcons*dsign(frdv,tloss)
      dv = (dv*fnorm)*dexp(-((qdl*(dble(m)-2.0d0))/(dble(loop)-2.0d0)))
      if (dabs(dv/veloc(m)).gt.0.3d0) dv =-(veloc(m)*0.3d0)
      veloc(m+1) = veloc(m)+dv
      velav = veloc(m)+(dv/2.0d0)
      u = veloc(m+1)
      u2 = u*u
      write(*, 1230) m, t, tloss, tstep, dh, de, velav
 1230 format(1h ,i3,2x,f11.1,1x,5(1pg10.3,2x))
c
      xh(m+1) = xhyf
      xhav = (xh(m)+xh(m+1))/2.0d0
      do 51 j = 1, atypes
      do 50 i = 1, maxion(j)
      popf(i,j) = pop(i,j)
      pop(i,j) = (popi(i,j)+popf(i,j))/2.0d0
   50 continue
   51 continue
c
      humag = (hmago*veloc(2))/u
      gam = (((pres*u)/fm)-u2)-humag
      dhy(m+1) = (dhy(2)*veloc(2))/veloc(m+1)
      dh1 = dhy(m)
      dh2 = dhy(m+1)
      dh = (dh1+dh2)/2.0d0
      de1 = feldens(dh1,popi)
      de2 = feldens(dh2,popf)
      deel(m+1) = de2
      de = (de1+de2)/2.0d0
      mu2 = fmua(de2,dh2)
      rho1 = frho(de1,dh1)
      rho2 = frho(de2,dh2)
      rho = (rho2+rho1)/2.0d0
      te(m+1) = (mu2*gam)/rgas
      if (te(m+1).lt.1.0d0) te(m+1) = 1.0d0
      t = (te(m+1)+te(m))/2.0d0
      dte = dabs(te(m+1)-te(m))/dmax1(te(m),te(m+1))
      if (dte.gt.0.92d0) frdv = frdv/2.0d0
      qf = (2.0d0*ww)-(((5.0d0*gam)+(4.0d0*humag))+u2)
      dq = (qf-qs)/2.0d0
      iflag = iflag+1
      tlosrr = tloss
      tstepa = tstep
      dr = ((veloc(m+1)+veloc(m))*tstep)/4.0d0
c
c   ***DETERMINE CASE A,B AND TEST CONVERGENCE OF COOLING RATE
c
      lmod = 'CAB'
      dtemp1 = 1.0d0
      dtemp2 = 0.0d0
      dtemp3 = 0.5d0
      call totphot(t, dh, dtemp1, dtemp2, dr, dv, dtemp3, lmod)
      call cool(t, de, dh)
c
      if (dlos.lt.0.6d0) tloss = 0.3d0*eloss
      tlosav = tloss
      tstep = (dq*rho)/tloss
      tstepf = tstep
      fcool = dabs((tloss-tlosrr)/(tloss+1.d-38))
      ftim = dabs((tstep-tstepa)/(tstep+1.d-30))
c
      if (((fcool.lt.0.005).and.(ftim.lt.0.005)).and.(frdv .eq. 
     &frdv0)) goto 110
      do 61 j = 1, atypes
      do 60 i = 1, maxion(j)
      pop(i,j) = popi(i,j)
   60 continue
   61 continue
      goto 102
  110 continue
c
c
c    ***COMPUTES DISTANCES AND FLUXES
c    PREPARES NEXT TIMESTEP
c
      dstep = ((veloc(m+1)+veloc(m))*tstep)/2.0d0
      dist(m+1) = dist(m)+dstep
      timlps(m) = timlps(m-1)+tstep
      tlosa0 = tlosa1
      tlosa1 = tlosav
      tstep0 = tstep1
      tstep1 = tstepf
c
      rdis = (dist(m)+dist(m+1))/2.0d0
      fi = 1.0d0
      imod = 'ALL'
c
      call sumdata(t, de, dh, fi, dstep, dstep, rdis, imod)
c
c
      if (ox3.ne.0) then 
        hoiii(m) = (fluxf6(7,ox3)+fluxf6(7,ox3))/(fluxh(2)+epsilon)
      endif
      if (ox2.ne.0) then 
        hoii(m) = (fluxf(1,ox2)+fluxf(2,ox2))/(fluxh(2)+epsilon)
      endif
      if (ni2.ne.0) then 
        hnii(m) = (fluxf6(7,ni2)+fluxf6(10,ni2))/(fluxh(2)+epsilon)
      endif
      if (su2.ne.0) then 
        hsii(m) =(fluxf(1,su2)+fluxf(2,su2))/(fluxh(2)+epsilon)
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
      write(luop, 1200) m, tlosav, tstep, pop(1,1), pop(2,1), (
     &pop(j,iel1),j = 1, 6), (pop(i,iel2),i = 1, 6)
 1200 format(1h ,i3,1pg10.2,1pg9.2,1x,14(0pf7.4))
      write(*, 1235) t, tloss, tstep, dh, de, caseab(1)
c
 1235 format(1h ,5x,f11.1,1x,4(1pg10.3,2x),5hCAB :,0pf5.2/)
      if ((real(m).gt.(real(loop)/2.0))
     &.and.(xhav.lt.(xhmin/1.d3))) goto 788
      if ((t.lt.1000.d0).and.((de/dh).lt.
     &(xhmin*dsqrt(1000.d0/(t+1.d0))))) goto 788
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
      do 66 j = 1, atypes
      do 65 i = 1, maxion(j)
      popi(i,j) = popf(i,j)
      pop(i,j) = popf(i,j)
   65 continue
   66 continue
    1 continue
      m = m-1
  788 continue
c
c    ***MODEL FINISHED, TRANSFER TO DATA OUTPUT PHASE
c    PERFORM AVERAGES OF RELEVANT QUANTITIES
c
c
c    ***AVERAGES AND OUTPUT STORED QUANTITIES
c
      call avrdata
c
      write(luop, 1250) 
 1250 format(//' PROGRAM ENDED , LIST OF STORED QUANTITIES : '//
     &' L#',t9,'<TE>',t16,'<FHII>',t24,'<DENS.>',t34,'<EL.DENS>',t43,
     &'<DIST>',t53,'<V>',t63,'<ELPS.TM>',t74,'OIII/H-B',t84,'OII/H-B'
     &,t94,'NII/H-A',t104,'SII/H-A'/)
      do 1240 l = 2, m
      distav = (dist(l)+dist(l+1))/2.0d0
      velav = (veloc(l+1)+veloc(l))/2.0d0
      t = (te(l)+te(l+1))/2.0d0
      dh = (dhy(l)+dhy(l+1))/2.0d0
      de = (deel(l)+deel(l+1))/2.0d0
      xhy = (xh(l)+xh(l+1))/2.0d0
      write(luop, 1280) l, t, xhy, dh, de, distav, velav, 
     &timlps(l), hoiii(l), hoii(l), hnii(l), hsii(l)
 1280 format(1h ,i3,f10.0,f7.4,9(1pg10.3))
 1240 continue
c
c
c    ***OUTPUT SPECTRUM AND AVERAGE IONIC POP. AND TEMP.
c
c
      call wrsppop(luop)
c
      close(luop) 
      write(*, 2005) 
c
 2005 format(//'OUTPUT CREATED &&&&&&&&&&&&&&&&&&&&&&'/)
      return 
      end
