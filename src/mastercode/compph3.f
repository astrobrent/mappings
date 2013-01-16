cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       PHOTOIONISATION MODEL
c       DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c       SPACE STEPS DERIVED FROM CHANGE IN OPTICAL DEPTH  :  DTAU
c
c       CALL SUBR.: TEEQUI,COPINTO,TOTPHOT,NEWDIF,LOCALEM,EVOLTEM,
c                   SUMDATA,AVRDATA,WRSPPOP,ZETAEFF,TAUDIST
c
c
c
      subroutine compph3(dhn, fin, banfil)
c
      include 'cblocks.inc'
c
c
      double precision aadn,agmax,agmu,aper,aperi,apero
      double precision a,b,c,x,fx,fp,n0,r0
      double precision de,dedhma,dedhmi,deeq,dh,dhf,dhl,dpw,dxp
      double precision dhlma,dhn,dhp,difg,difma,difmi,difp0
      double precision dift,disin,dison,dlos0,dlos1,dlosav
      double precision dr,drmu,drp,drxp,drzf,drzi,dtau,dtauon
      double precision dtaux,dtco,dtco0,dtco1,dtl,dtlma,dtoff
      double precision dton,dton0,dton1,durmax,dv,dvoluni
      double precision eqscal,exl,exla
      double precision fhi,fhii,fi,fi0,fi1,fidhcc,fin,frdw
      double precision frx,fthick
      double precision poppre(mxion, mxelem)
      double precision popend(mxion, mxelem), ppre(mxion, mxelem)
      double precision prescc,prf,protim
      double precision rad,radi,recscal,rmm,rww,runit
      double precision t,taux,te1a,tef,telast
      double precision temp1,temp2,tempre,tep,tex,texm
      double precision tii0,tprop,trea,treach,treach0,treach1,tspr
      double precision wd,wdil0,wdil1,wei,weiph
c
      integer*4 m,maxio,mma,n,nidh,niter
      integer*4 i,j,luop,lupf,lusl,lut0
      integer*4 ia_sjm(mxelem),ia(6)
c
      character imod*4, jmod*4, lmod*4, mmod*4, nmod*4,jjd*4
      character linemod*4, spmod*4
      character newfil*20, filna*20, banfil*20, fnam*20, filnb*20
c
c     External Functions
c
      double precision fcrit,fdilu, feldens
      double precision fpressu,frectim
c
c
c     internal functions
c
      double precision dif, frad, fsight, fring
c
      dif(a,b) = dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
      frad(x,n0,fx,fp,a,b,c,r0)  = n0*(fx*dexp((x-a)/r0)
     &        + fp*((x/a)**b))+c
      fsight(aper,radi) = 1.d0-(dcos(dasin(dmin1(1.d0,aper/(radi+(
     & 1.d-38*aper)))))*(1.d0-(dmin1(1.d0,aper/(radi+(1.d-38*aper)))
     & **2)))
      fring(aperi,apero,radi) = dmax1(0.d0,fsight(apero,radi)-fsight(
     & aperi,radi))
c
      jcon = 'NO'
      jspot = 'NO'
      limph = atypes

c
      mma = 7/dtau0
      exla = 5.d-5
      trea = 1.d-5
      drmu = 3.0d0
      difma = 0.05d0
      dtlma = 0.03d0
      dhlma = 0.010d0
      dton0 = 0.0d0
      dton1 = 0.0d0
      dtoff = 0.0d0
      scalen = 1.0d16
c
      luop = 21
      lupf = 27
      lusl = 45
      lut0 = 0
c
      ia(1) = 1
      ia(2) = 2
c
      ia_sjm(1) = 1
      ia_sjm(2) = 2
c
c***    CHOICE OF ATOMIC ELEMENTS TO OUTPUT IN APN FILE
c
      ia(3) = 3
      ia(4) = 4
      ia(5) = 5
      ia(6) = 9
c
      ia_sjm(3) = 3
      ia_sjm(4) = 4
      ia_sjm(5) = 5
      ia_sjm(6) = 6
      ia_sjm(7) = 7
      ia_sjm(8) = 8
      ia_sjm(9) = 9
      ia_sjm(10) = 10
      ia_sjm(11) = 11
c
c    ***DERIVE FILENAME : PHN**.PH3
c
      call p3fheader(newfil,filna,banfil,fnam,filnb,luop,lupf,lusl,
     &               ia,dhn,fin)
c
c
c*********************************************************************
c
c    ***COMPUTATION STARTS HERE
c
c*********************************************************************
c
      fi = fin
      fi0 = fin
      fi1 = fin
      fidhcc = fin*dhn
      dv = 0.0d0
      runit = (rmax-remp)/300.d0
      dtau = dtau0/2.0d0
      if (jgeo .eq. 'S') then
      vunilog = 3.0d0*dlog10(runit)
      else
      vunilog = dlog10(runit)
      end if
      if ((jeq .eq. 'E').or.(jeq .eq. 'P')) then
      nmod = 'EQUI'
      else
      nmod = 'TIM'
      end if
      if (jden .ne. 'B') then
      agmax = 0.05d0
      aadn = 2.5d0
      else
      agmax = 0.80d0
      aadn = 5.0d0
      end if
      agmu = agmax
      m = 0
c
c***********************************
c
c     ***ITERATE, M = STEP NUMBER
c
c***********************************
c
 1    m = m+1
c
      jspot = 'NO'
      niter = 0
      nidh = 0
      difmi = 1.d30
      wei = fcrit(dtau/dtau0,1.0d0)
      dtau = (dtau0+((3.d0*wei)*dtau))/(1.0d0+(3.d0*wei))
      wei = 3.d0
      agmu = (agmax+(wei*agmu))/(1.d0+wei)
c
c****************************************************************
c
c    ***IF M=1 , CHECK CONVERGENCE FOR THE DENSITY AND IONIC POP.
c       AT THE INNER BOUNDARY
c
c****************************************************************
      if (m .eq. 1) then
         te0 = 1.d4
         dr = 0.d0
         drp = dr
         if (jden .ne. 'F') then
            dh0 = dhn
         else
            dh0 = frad(remp,dhn,xfac,pfac,afac,bfac,cfac,scalen)
         end if
         de0 = feldens(dh0,pop0)
         dedhmi = 0.d0
         dedhma = 0.d0
         do 127 j = 1, 11
            dedhma = dedhma+zion(j)
            if (arad(2,j).le.0.d0) dedhmi = dedhmi+zion(j)
 127     continue
         call copinto(pop0, pop)
 12      continue
c
c
         if (jgeo .eq. 'S') then
            rad0 = remp
            dis0 = rad0
            wdil0 = fdilu(rstar,rad0)
         else
            rad0 = 0.d0
            dis0 = remp
            wdil0 = wdpl
         end if
c
c
         if ((fin.ne.1.d0).and.(jden.eq.'B')) fi0=dmin1(1.d0,fidhcc/dh0)
         lmod = 'DW'
         call localem(te0, de0, dh0)
c
         call totphot(te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
c
         if ((jeq .eq. 'E').or.(jeq .eq. 'P')) then
            dton0 = 1.d33
            call teequi(te0, te0, de0, dh0, dton0, nmod)
         else
            distim(1,0) = dis0
            distim(2,0) = dis0/cls
            treach0 = distim(2,0)
            dton0 = dmax1(0.d0,dmin1(tlife,telap-((2.d0*dis0)/cls)))
            jjd = jden
            jden = 'C'
            exl = -1.d0
            tex = 0.95d0*tm00
            prescc = fpressu(1.d4,dhn,pop)
c
            call copinto(pop0, pop)
            call evoltem(tm00,te0,de0,dh0,prescc,dton0,exl,tex,lut0)
c
            jden = jjd
c
            if (dton0.le.0.d0) then
              write(*, 305) dton0
 305          format(/,' AGE WAS WRONGLY SET SINCE IT IS SHORTER THAN',
     &         /,' NEEDED AT THE FIRST SPACE STEP , DTON0 :',1pg10.3,/)
              return 
            end if
            dtco0 = (1.5d0*fpressu(te0,dh0,pop))/(eloss+1.d-36)
         end if
c
         dhp = dh0
         tep = te0
         dlos0 = dlos
c
         if (jden .eq. 'C') then
            dh0 = dhn
         else if (jden .eq. 'F') then
            dh0 = frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen)
         else
            prescc = fpressu(1.d4,dhn,pop)
            dh0 = (dhn*prescc)/fpressu(te0,dhn,pop)
         end if
         dhl = dif(dhp,dh0)
         write(*, 472) te0, dhl
         if (dhl.ge.dhlma) goto 12
         tspr = dton0
c
         tempre = te0
         call copinto(pop, poppre)
         call copinto(pop, ppre)
c
c***********************************************
c     END OF INNER BOUNDRY FOR M = 1
C
      end if
c
c***********************************************
c
c
c***********************************************
c
c    ***ESTIMATE OF VALUES AT THE OUTER BOUNDARY
c
c***********************************************
c
      call copinto(poppre, pop)
      call copinto(poppre, popend)
c
      te1 = te0+(((te0-tep)*te0)/tep)
      te1 = dmax1(te0*agmu,dmin1(te0/agmu,te1))
      te1a = te1
      t = (te1+te0)/2.d0
      if (jden .eq. 'C') then
         dh = dhn
      else if (jden .eq. 'F') then
         dh1 = frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen)
      else
         dh = (dhn*prescc)/fpressu(t,dhn,pop)
      end if
c
c    ***DETERMINE SPACE STEP : DR
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
      mmod = 'DIS'
      lmod = 'SO'
c
      call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
      call taudist(dh, fi, dtau, dr, rad0, pop, mmod)
c
   30 continue
c
      if ((dr.gt.(drmu*drp)).and.(m.gt.1)) dr = drmu*drp
      if (jgeo .eq. 'S') then
         if (jden .eq. 'F') then
            if ((frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen))
     &         .gt.ofac) then
	            dpw = dr
	            if (pfac.gt.0.0d0) then
		            if (bfac.gt.0.0d0) then
		               dpw = dis0*((1.1d0**(1/bfac))-1.d0)
		            else if (bfac.lt.0.0d0) then
		               dpw = dis0*((0.9d0**(1/bfac))-1.d0)
		            endif
	            endif
	            dxp = dr
	            if (xfac.gt.0.d0) then
	               dxp = dabs(0.05d0*scalen)
	            endif
	            dr = 1.d0/((1.d0/dpw)+(1.d0/dxp)+(1.d0/dr))
            endif
         endif
         rad1 = rad0+dr
         dis1 = rad1
         wdil1 = fdilu(rstar,rad1)
         if (abc9.lt.huge) then
            dvoluni = ftpi*((fring(abc0,abc9,rad1)*((rad1/runit) ** 3))
     &      -(fring(abc0,abc9,rad0)*((rad0/runit) ** 3)))
         else
            dvoluni = ftpi*(((rad1/runit) ** 3)-((rad0/runit) ** 3))
         end if
      else
         dis1 = dis0+dr
         rad1 = 0.d0
         wdil1 = wdpl
         dvoluni = dr/runit
      end if
c
c     *** dr done
c
c************************************************************
c
c    ***CALCULATION OF IONIC POPULATION AT THE OUTER BOUNDARY
c
c************************************************************
      if (jden .eq. 'C') then
         dh1 = dhn
      else if (jden .eq. 'F') then
         dh1 = frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen)
      else
         dh1 = (dhn*prescc)/fpressu(te1,dhn,popend)
      end if
      de1 = feldens(dh1,popend)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi1=dmin1(1.d0,fidhcc/dh1)
      rww = 1.d0
      t = (te0+(rww*te1))/(1.d0+rww)
      dh = (dh0+(rww*dh1))/(1.d0+rww)
      de = (de0+(rww*de1))/(1.d0+rww)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
      lmod = 'DW'
      wei = 1.d0/(1.d0+rww)
c
      call averinto(wei, poppre, popend, pop)
      call localem(t, de, dh)
      call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)
c
c     PRINT OUT TPHOT AS ERGS AND PHOTONS
c
c      write (*,*) '**************************************************'
c      write (*,*) m,'th tphot array by bin number and eV:'
c      write (*,*) 'As Photons/cm2/s & skipbin'
c      rhoion = densnum(dh)
c      write (*,*) 'ion density = ',rhoion
c      write (*,*) '**************************************************'
c      do 3306 itph = 1, 230
c         teng = 0.5d0*(ephot(itph)+ephot(itph+1))
c         tphots = ((tphot(itph)/ev)/teng)*cls
c         write (*,*) itph,teng,tphots,skipbin(itph)
c 3306 continue
c      write (*,*) '**************************************************'
c      write (*,*) ' '
c
      call zetaeff(fi, dh)
c
c***********************************************
c    *EQUILIBRIUM IONISATION AND TEMPERATURE .
c***********************************************
c
      if ((jeq .eq. 'E').or.(jeq .eq. 'P')) then
         dton1 = 1.d33
         call copinto(popend, pop)
         call teequi(te1, tef, de1, dh1, dton1, nmod)
      else
c
c      *FINITE AGE ,NON EQUILIBRIUM IONISATION , DETERMINE TIME STEP
c      TAKES INTO ACCOUNT IONISING FRONT VELOCITY .
c
         distim(1,m) = dis1
         distim(2,m) = distim(2,m-1)+(dr/(viofr+(1.d-35*dr)))
         tprop = distim(2,m)
         dtauon = 1.0d0
         mmod = 'DIS'
         call taudist(dh, fi, dtauon, disin, 0.d0, pop0, mmod)
         fthick = dmin1(disin,dis1-distim(1,0))
         dison = dis1-fthick
         protim = distim(2,0)
         do 330 j = m-1, 0, -1
            if (distim(1,j).gt.dison) goto 330
            prf = (dison-distim(1,j))/((distim(1,j+1)-distim(1,j))+
     &      1.d-38)
            protim = distim(2,j)+(prf*(distim(2,j+1)-distim(2,j)))
            goto 340
 330     continue
 340     continue
c
c
         treach1 = dmax1(dis1/cls,protim+(fthick/cls))
         durmax = dmax1(0.d0,tlife+((dis1/cls)-treach1))
         dton1 = dmax1(0.d0,dmin1(durmax,(telap-(dis1/cls))-treach1))
         deeq = dmax1(de1/3.d0,exla*dh1)
         eqscal = 5.d0*frectim(2.d0*te1,deeq,dh1)
         if (dton1.lt.eqscal) then
c
c      *TIME STEP TOO SHORT , EQUILIBRIUM TEMPERATURE NOT REACHED
c
            jjd = jden
            jden = 'C'
            exl = -1.d0
            tex = 0.95d0*tm00
            call copinto(pop0, pop)
            call evoltem(tm00,tef,de1,dh1,prescc,dton1,exl,tex,lut0)
            jden = jjd
         else
c
c      *TIME STEP LONG ENOUGH TO ASSUME EQUILIBRIUM TEMPERATURE
c
            call copinto(pop0, pop)
            call teequi(te1, tef, de1, dh1, dton1, nmod)
         end if
         dtco1 = (1.5d0*fpressu(tef,dh1,pop))/(eloss+1.d-36)
      end if
c
c***********************************************
c  END EQUILIBRIUM IONISATION AND TEMPERATURE .
c***********************************************
c
      if (jden .eq. 'C') then
         dhf = dhn
      else if (jden .eq. 'F') then
         dhf = frad(dis1,dhn, xfac,pfac,afac,bfac,cfac,scalen)
      else
         dhf = (dhn*prescc)/fpressu(tef,dhn,pop)
      end if
c
c************************************************
c    ***COMPARES END VALUES WITH PREVIOUS STEP
c************************************************
c
      call difpop(pop, popend, trea, 11, dift)
      call difpop(pop, poppre, trea, 11, difp0)
      call copinto(pop, popend)
c
      drxp = dr
      dtl = dif(tef,te1)
      dhl = dif(dh1,dhf)
      difg = dmax1(dift/difma,dtl/dtlma,dhl/dhlma)
      niter = niter+1
      if (difg.le.difmi) then
         texm = te1
         dtaux = dtau
         difmi = difg
      end if
c
      write(*, 472) tef, dhl, dtl, dift, difg, difp0, dtau, 
     &dr
  472 format(t7,0pf8.0,5(0pf8.3),3(1pg9.2))
c
c***************************************************************
c    ***DERIVE AVERAGE QUANTITIES FOR THE SPACE STEP CONSIDERED
c***************************************************************
c
      rww = 1.d0
      t = (te0+(rww*tef))/(1.d0+rww)
      dh = (dh0+(rww*dhf))/(1.d0+rww)
      de = (de0+(rww*de1))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
      disav = (dis0+(rww*dis1))/(1.d0+rww)
      wei = 1.d0/(1.d0+rww)
      call averinto(wei, poppre, popend, pop)
      if (jgeo .eq. 'S') then
         rad = disav
         wdil = fdilu(rstar,rad)
      else
         wdil = wdpl
         rad = 0.d0
      end if
c
c**********************************************************
c    ***CHECK CONVERGENCE OF END VALUES WITH PREVIOUS STEP
c**********************************************************
c
      temp1 = ((de0/dh0)-dedhmi)/dedhma
      temp2 = (de1/(de0+1.d-7))*((dh0/dhf) ** 2)
      if (difg.ge.1.d0) then
         weiph = 0.5d0
         rmm = 1.d5
         if ((niter.gt.11).and.(te1 .eq. texm)) then
            write(*, 172) dhl, dtl, dift, difg
 172        format(' NO SPATIAL CONVERGENCE :',4(0pf8.3))
            goto 31
         else if (niter .eq. 12) then
            te1 = texm
            dtau = dtaux
            difmi = 0.d0
         else if (((temp1.gt.0.01d0).and.(temp2.lt.0.2d0)).or.(dhl
     &      .ge.0.07d0)) then
            agmu = dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-((1.d0-agmu)*
     &      0.6d0))
            nidh = nidh+1
            if (nidh.lt.3) niter = niter-1
            te1 = (te0+dmax1(agmu*te0,dmin1(te0/agmu,te1)))/2.d0
            weiph = 0.85d0
            rmm = 0.25d0
            dtau = (0.8d0*rmm)*dtau
         else if (((niter.eq.3).or.(niter.eq.6)).and.((jden.ne. 
     &      'B').or.(difg .eq. (dift/difma)))) then
            agmu = dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-
     &      ((1.d0-agmu)*0.65d0))
            nidh = max0(0,min0(nidh-2,idnint(1.0d0+dint(niter/3.0d0))))
            te1 = (((te1+te0)+te1a)+te0)/4.0
            te1 = dmax1(te0*agmu,dmin1(te0/agmu,te1))
            weiph = 0.75d0
            rmm = 0.4d0
            dtau = (0.8d0*rmm)*dtau
         endif
            te1 = tef
            goto 120
 34         dtau = dtau/2.d0
 120        continue
c
            mmod = 'DIS'
            lmod = 'SO'
            call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
            call taudist(dh, fi, dtau, drzi, rad0, poppre, mmod)
            call taudist(dh, fi, dtau, drzf, rad0, popend, mmod)
c
            dr = dexp((weiph*dlog(drzi+.1d0))+((1.d0-weiph)*dlog(drzf
     &      +.1d0)))
            if ((dr.gt.(rmm*drxp)).and.(dtau.gt.1.d-16)) goto 34
            wei = dmin1(1.d0,(1.2d0*dr)/drxp)
            call averinto(wei, popend, poppre, popend)
            goto 30
 31     continue
      end if
c
      te1 = tef
      dlos1 = dlos
      telast = t
      frx = 1.d0-pop(jpoen,ielen)
c
c*****************************************************************
c    ***END  CHECK CONVERGENCE OF END VALUES WITH PREVIOUS STEP
c*****************************************************************
c
c
c*******************************************************
c    ***INTEGRATES COLUMN DENSITY AND END DIFFUSE FIELD
c*******************************************************
      if (jgeo .eq. 'S') then
         jmod = 'OUTW'
         frdw = 1.d0
      else
         jmod = 'LODW'
         frdw = 0.5d0
      end if
      imod = 'COLD'
      jspot = 'NO'
c
c
c
      call zetaeff(fi, dh)
      call localem(t, de, dh)
      call newdif(t, t, dh, fi, rad0, dr, dv, 0.d0, dv, frdw, jmod)
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
c
c*****************************************************************
c    ***CALCULATION OF AVERAGE QUANTITIES GIVING EQUAL WEIGHTS TO
c       EQUAL VOLUMES
c*****************************************************************
c
      if (jgeo .eq. 'S') then
         rww = (1.d0+((0.75d0*dr)/(rad0+dr))) ** 2
      else
         rww = 1.d0
      end if
c
      t = (te0+(rww*te1))/(1.d0+rww)
      dh = (dh0+(rww*dh1))/(1.d0+rww)
      de = (de0+(rww*de1))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav = (dis0+(rww*dis1))/(1.d0+rww)
      dlosav = (dlos0+(rww*dlos1))/(1.d0+rww)
      dton = (dton0+(rww*dton1))/(1.d0+rww)
      treach = (treach0+(rww*treach1))/(1.d0+rww)
      dtco = (dtco0+(rww*dtco1))/(1.d0+rww)
      wei = 1.d0/(1.d0+rww)
c
      call averinto(wei, poppre, popend, pop)
c
      if (jgeo .eq. 'S') then
      rad = disav
      wdil = fdilu(rstar,rad)
      else
      wdil = wdpl
      rad = 0.0d0
      end if
c
c*************************************************************
c    ***CALCULATION OF IONIC POPULATION AND TEMPERATURE WHEN
c       OR IF SOURCE SWITCHED OFF
c*************************************************************
c
      dtoff = dmax1(0.d0,(telap-((2.d0*disav)/cls))-tlife)
      if ((jeq .ne. 'E').and.(dtoff.gt.0.0)) then
         if (jeq .eq. 'P') recscal = frectim(t,0.8d0*de,dh)
         wd = abc3*wdil
         lmod = 'SO'
         jspot = 'YES'
         tii0 = t
c
         call totphot(t, dh, fi, rad, dr, dv, wd, lmod)
c
         call evoltem(tii0, t, de, dh, prescc, dtoff, exla, tm00, lut0)
c
         if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
      end if
c
c*********************************************************
c
c    ***COMPUTES SPECTRUM OF THE REGION CONSIDERED
c
c*********************************************************
c
      imod = 'REST'
      call cool(t, de, dh)
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
c
c    ***OUTPUT QUANTITIES RELATED TO SPACE STEP
c
      fhi = pop(1,1)
      fhii = pop(2,1)
      open(unit=luop, file=filna, status='OLD', access='APPEND')
      open(unit=lupf, file=filnb, status='OLD', access='APPEND')
      open(unit=lusl, file=newfil, status='OLD', access='APPEND')
      write(luop, 461) m, t, dlosav, dis1, disav, dr, fi, dh, 
     &de, fhi, fhii, zetae, qhdh
      if (jeq .eq. 'F') write(luop, 462) fthick, treach, dtco, 
     &dton, dtoff
  462 format(t71,5(1pg9.2))
      if (jeq .eq. 'P') write(luop, 473) recscal, dtoff
  473 format(t98,5(1pg9.2))
  460 format(i4,0pf8.0,1pg9.1,12(1pg9.2))
  461 format(i4,0pf8.0,1pg9.1,11(1pg9.2))
c
      write(*, 460) m,t,dlosav,dton,dtoff,disav,dr,dh,fhi
      write(lupf, 478) m, ((dmax1(0.d0,dmin1(0.9999d0,1.d-4*
     &idnint(1.d4*pop(j,ia(i))))),j = 1, maxion(ia(i))-1),i = 1, 6)
  478 format(1h ,i3,26(0pf5.4))
      write(lusl, 678) m, ((pop(j,ia_sjm(i)),j = 1, maxion(
     &ia_sjm(i))),i = 1, 11), t, de, dh, fi, dvoluni, dr, disav
  678 format(i3,66(e16.5),/)
      close(luop) 
      close(lupf) 
      close(lusl) 
c
c     ***END OUTPUT
c
c
c**********************************************************************
c
c       CHECKPOINT BATCHMAP
c      call cpt
c
c**********************************************************************
c
c
c     ***TEST ENDING CONDITIONS
c
      do 513 n = 1, ionum
      if ((atpho(n).eq.ielen).and.
     &(ionpho(n).eq.jpoen)) goto 517
  513 continue
  517 taux = sigpho(n)*popint(jpoen,ielen)
      if (((jend .eq. 'A').or.(jend .eq. 'B')).and.(frx.le.fren)) 
     &goto 505
      if ((jend .eq. 'C').and.(telast.le.tend)) goto 505
      if ((jend .eq. 'D').and.(taux.ge.tauen)) goto 505
      if ((jend .eq. 'E').and.(dis1.ge.diend)) goto 505
c
c     removed until it`s purpose is clear..... 8/90
c
c      if ((((de1/dh1)-dedhmi)/dedhma).le.exla) goto 505
c
c
c    ***RESET INNER BOUNDARY QUANTITIES FOR NEXT SPACE STEP
c
      call copinto(popend, poppre)
      tep = te0
      te0 = te1
      de0 = de1
      dh0 = dh1
      dton0 = dton1
      treach0 = treach1
      dtco0 = dtco1
      fi0 = fi1
      rad0 = rad1
      dis0 = dis1
      wdil0 = wdil1
      drp = dr
      dlos0 = dlos1
c
c**********************************************************************
c
c     LOOP BACK AND INCREMENT M, DO NEXT STEP.....
c
c**********************************************************************
c
      goto  1
c
c
c
  505 continue
c***********************************************************
c
c    ***MODEL ENDED ; OUTPUT RESULTS ***********************
c
c***********************************************************
c
      open(luop, file=filna, status='OLD', access='APPEND')
      open(lupf, file=filnb, status='OLD', access='APPEND')
      write(luop, 530)
 530  format(//' Model ended',t20,'Tfinal',t28,'DISfin',t38,'FHXF',t48,
     &'TAUXF',t58,'Ending')
      write(luop, 510) t, dis1, frx, taux, jend
 510  format(12h -----------,t18,0pf8.0,3(1pg10.3),4x,a4)
      if ((jeq .ne. 'E').and.(abc3.gt.0.0)) write(luop, 519
     &) abc3
 519   format(/' Final fractional luminosity of source after ',
     &'turn off :',1pg10.3)
      if ((abc9.lt.huge).and.(jgeo .eq. 'S')) write(luop,
     &      967) abc0, abc9
 967  format(/,' NB. :::::::::::::::Integeration through the line',
     & 'of sight for a ring aperture of radii :',2(1pg10.3))
      write(lupf, 2047) tempre, tspr
 2047 format(//t4,'Preionisation conditions for step#0 ',
     &'for all elements   (TEpr :',0pf8.0,' Time :',1pg9.2,' )  :'/)
      write(lupf, 2050) (elem(i),i = 1, atypes)
 2050 format(1h ,t8,16(4x,a2,4x)/)
      maxio = 0
      do 2200 j = 1,atypes
         if (maxio.le.maxion(j)) maxio = maxion(j)
 2200 continue
      do 2100 j = 1,maxio
      write(lupf, 2070) rom(j), (ppre(j,i),i = 1
     &, atypes)
 2070 format(a7,t8,16(1pg10.3))
 2100 continue
c
c
c    ***COMPUTE AVERAGE SPECTRUM AND WRITES IT IN FILE PHN
c
      call avrdata
      call wrsppop(luop)
      linemod = 'LAMB'
      spmod   = 'REL'
      call spec2(luop,linemod,spmod)
      close(lupf)
      close(luop)
      write(*, 2005) fnam
 2005  format(//' Output created &&&&&& File : ',a/)
c
      return 
c
c
      end
