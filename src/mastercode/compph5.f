cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       PHOTOIONISATION MODEL
c       DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c
c       RSS 1999
c
c       PHOTO5:  As for P4 but adds integrated radiation pressure
c               (altered for dust pressure)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compph5(dhn, fin, banfil,difma,dtlma,dhlma)
c
      include 'cblocks.inc'
c
      double precision poppre(mxion, mxelem), popend(mxion, mxelem)
      double precision ppre(mxion, mxelem),difma,dtlma,dhlma
      double precision thist(3),dpw,dxp
      double precision popz(mxion, mxelem),difp,treap
      double precision ctime,rtime,ptime,eqtime,stime
c
      double precision aadn,agmax,agmu,aper,aperi,apero
      double precision a,b,c,x,fx,fp,n0,r0,fv,e,f
      double precision de,dedhma,dedhmi,deeq
      double precision dh,dhf,dhl,dhn,dhp,frp
      double precision difg,difmi,difp0,dift
      double precision disin,dison,dlos0,dlos1,dlosav
      double precision dr,drp,drxp,drzf,drzi
      double precision dtau,dtauon,dtaux,dtco,dtco0,dtco1
      double precision dtl,dtoff,dton,dton0,dton1,durmax
      double precision dv,dva,dvoluni,eqscal,exl,exla
      double precision fhi,fhii,fi,fi0,fi1,fidhcc,fin,frdw,frx,fthick
      double precision prescc,prf,protim,rad,radi,recscal,rmm,rww
      double precision t,taux,te1a,tef,telast,temp1,temp2,tempre
      double precision tep,tex,texm,tii0,tprop,trea,runit 
      double precision treach,treach0,treach1,tspr
      double precision wd,wdil0,wdil1,wei,weiph
      double precision ponk, radpint, pres0, pres1
      double precision popinttot
      double precision A_v,Ffuv,dnu,eng,Hab,pahtest(mxinfph)
      double precision pahsum, pahQ0, pahQion, pahQfac
      double precision q1,q2,q3,q4
      double precision def
c
      integer*4 np,luop,lups,lusp,lusl,lupf,IR1
      integer*4 lut0,m,maxio,mma,n,nidh,niter
      integer*4 luions(4),i,j,k
      integer*4 fuvmin,fuvmax,pahabsmax
c
c     Old computer time calculations
c      real tarray(2)
c      real dtarr,dtime
c      integer*4 nt0,nt1, time
c
      character jjd*4, pollfile*12,tab*4
      character newfil*20, filna*20, banfil*20, fnam*20, filnb*20
      character imod*4, jmod*4, lmod*4, mmod*4, nmod*4,wmod*4,ispo*4
      character linemod*4, spmod*4
      character filnc*20, filnd*20,filn(4)*20, fspl*20
      character pfx*8,caller*4,sfx*4
      character IRfile*20,chargefile*20
      logical iexi
c
c     External Functions
c
      double precision fcolltim,fcrit,fdilu,feldens,fradpress
      double precision fphotim,fpressu,frectim,frectim2
c
c     internal functions
c
      double precision dif,frad,fsight,fring
c
      dif(a,b) = dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
      frad(x,n0,fx,fp,a,b,c,r0,fv,e,f)  = 
     &    n0*(fx*dexp((x-a)/r0)+fp*((x/a)**b)+fv*(((r0-x)/e)**(-f)))+c
      fsight(aper,radi) = 1.d0-(dcos(dasin(dmin1(1.d0,aper/(radi+(
     &    1.d-38*aper)))))*(1.d0-(dmin1(1.d0,aper/(radi+(1.d-38*aper)))
     &    **2)))
      fring(aperi,apero,radi) = dmax1(0.d0,fsight(apero,radi)-fsight(
     &    aperi,radi))
c
      jspot = 'NO'
      jcon = 'YES'              ! Enable full continuum calculation
      ispo = 'SII'
      IRfile=' '
      tab = char(9)
c
      limph = atypes
c
c      nt0 = time()
c      nt1 = 0

      mma = nstmax-1
c      scalen = 1.0d16
      exla = 5.d-5
      trea = 1.d-5
      dton0 = 0.0d0
      dton1 = 0.0d0
      dtoff = 0.0d0
c
c     Define logical unit numbers
      luop = 21
      lups = 22
      lusp = 23
      lupf = 27
      lusl = 45
      lut0 = 0
      IR1 = 91
c
      ieln = 4
      do i = 1,ieln
         luions(i) = 22+i
      enddo
c
c
c    ***DERIVE FILENAME : PHN**.PH5
c
      call p5fheader(newfil,banfil,fnam,filna,filnb,filnc,filnd,filn,
     &     luop,lupf,lusl,lusp,lups,luions,dhn,fin)

      chargefile = ' '
c
c  Set up for PAH existence
c
      if (grainmode.AND.(pahCfrac.ne.0.d0)) then
      	k=1
	do while (ephot(k).le.6.0)
	  k=k+1
	enddo
c     FUVMIN = photon bin corresponding to 6 eV energy
	fuvmin=k-1
	do while (ephot(k).le.1.359907296E+1)
	  k=k+1
	enddo
c     FUVMAX = photon bin corresponding to 13.6 eV energy
	fuvmax=k-1
	do while (ephot(k).le.262.012)
	  k=k+1
	enddo
c     PAHABSMAX = photon bin corresponding to 262.0 eV energy
c       maximal PAH absorption bin
	pahabsmax=k-1
c
c Calculate photon flux in FUV range using Weingartner & 
c   Draine (2001) ISRF 
c 
	Ffuv=0.d0
	pahsum=0.d0
	do k=fuvmin,fuvmax
	  eng=(ephot(k+1)+ephot(k))/2.d0
	  dnu=(ephot(k+1)-ephot(k))*evplk
	  pahsum=pahsum+pahnabs(k)*dnu
	  if (eng.lt.9.26) then
	    Ffuv=Ffuv+2.548E-18*eng**(-1.3322)/eV*pahnabs(k)*dnu
	  elseif (eng.lt.11.2) then
	    Ffuv=Ffuv+1.049E-16*eng**(-3.d0)/eV*pahnabs(k)*dnu	    
	  else !11.2<eng<13.6
	    Ffuv=Ffuv+4.126E-13*eng**(-6.4172)/eV*pahnabs(k)*dnu
	  endif
	enddo
	pahQ0=Ffuv/pahsum
	write(*,520) pahQ0
 520	format(/'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/
     &        'ISRF photon flux :',1pg12.5,/
     &	       '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      endif
c
c*********************************************************************
c
c   The routine does some initial setup, guesses at the parameters at 
c   the inner boundary of the first spatial step, then for each spatial
c   step it guesses at the parameters at the outer boundary of the step,
c   iterates till these parameters converge, and then iterates over 
c   spatial steps till the completion criteria are achieved.
c
c      write(*,*)'***COMPUTATION STARTS HERE'
c
c*********************************************************************

      fi = fin                  ! Average filling factor in space step
      fi0 = fin                 ! Filling factor at inner step boundary
      fi1 = fin                 ! Filling factor at outer step boundary
      fidhcc = fin*dhn          !
      dv = 0.0d0
      runit = (rmax-remp)/300.d0
      dtau = dtau0/2.d0
      if (jgeo .eq. 'S') then
        vunilog = 3.0d0*dlog10(runit)
      else
        vunilog = dlog10(runit)
      end if
      if ((jeq .eq. 'E').or.(jeq .eq. 'P').or.(jeq.eq.'C')) then
c     Equilibrium structure (E=equilibrium, P=post-equilibrium decay, 
c     C=collisional ionization equilibrium)
         nmod  = 'EQUI'
      else
c     Equilibrium structure (F=finite source)
         nmod  = 'TIM'
      end if
c
c  set convergence criterion parameters
c  agmax = ? ; aadn = ?  
c
      if (jden .ne. 'B') then
c        Density structure (B=isobaric)
         agmax = 0.05d0
         aadn  = 2.5d0
      else
c     Density structure (C=isochroic, F=function, I=input file)
         agmax = 0.80d0
         aadn  = 5.0d0
      end if
      agmu  = agmax
      
      ponk  = 0.d0
      if (jden.eq.'B') ponk  = dhn*1e4
      radPint = 0.d0
c
c***********************************************
c
c      write(*,*)'***ITERATE, M = STEP NUMBER'
c
c***********************************************
c
      m = 1        ! Spatial step counter
c
 1    continue     ! Loop over steps through Nebula
c
      jspot = 'NO'
      niter = 0                 ! Convergence loop counter
      nidh  = 0
      difmi = 1.d30
      wei   = fcrit(dtau/dtau0,1.0d0) !Sinc function --WHY?
      dtau  = (dtau0+((19.d0*wei)*dtau))/(1.0d0+(19.d0*wei))
      wei   = 3.d0
      agmu  = (agmax+(wei*agmu))/(1.d0+wei)
c
c****************************************************************
c
c      write(*,*)'IF M=1 , CHECK CONVERGENCE FOR THE DENSITY'
c      write(*,*)'AND IONIC POP. AT THE INNER BOUNDARY'
c
c****************************************************************
      if (m .eq. 1) then
c
         if (jeq.eq.'C') then
c     Equilibrium structure (C=collisional ionization equilibrium)
c     TE0 = Electron temperature at inner spatial step boundary
            te0 = frad(remp,tofac,txfac, tpfac,tafac,tbfac,tcfac,tscalen
     &		,0.d0,1.d0,0.d0)
         else
            te0 = 1.d4
         endif
c
         if (jthm.eq.'T') te0 = Fixtemp
c
         dr = 0.d0
         drp = dr
         if (jden .ne. 'F') then
c     Density structure (C=isochroic, B=isobaric, I=input file)
c     DH0 = Hydrogen density at inner spatial step boundary
            dh0 = dhn
         else
c     Density structure (F=function)
            dh0 = frad(remp, dhn, xfac, pfac, afac, bfac, cfac,
     &           scalen, vfac, efac, ffac)
         end if

c     Derive electron density using ionic populations in POP0.
c     DE0 = Electron density at the inner spatial step boundary
c     DEDHMA = total metallicity of the gas (relative to hydrogen)
c     DEDHMI = total metallicity of gas for all species with 
c              -ve or 0 recombination rates from 2nd ion state??? 
         de0 = feldens(dh0, pop0)
         dedhmi = 0.d0
         dedhma = 0.d0
         do j = 1, atypes
            dedhma = dedhma+zion(j)
            if (arad(2,j).le.0.d0) dedhmi = dedhmi+zion(j)
         enddo
         call copinto(pop0, pop)
c
c     Derive geometrical dilution factor using photon source radius RSTAR
c     and radial distance from source of the inner step boundary RAD0. 
         if (jgeo .eq. 'S') then
            rad0 = remp
            dis0 = rad0
            wdil0 = fdilu(rstar, rad0)
         endif
         if (jgeo .eq. 'F') then
            rad0 = 0.d0
            dis0 = remp
            wdil0 = fdilu(rstar,dis0)
         end if
         if (jgeo .eq. 'P') then
            rad0 = 0.d0
            dis0 = remp
            wdil0 = wdpl
         end if
c
c  Loop back here if hydrogen density change is too great
 12      continue   
c
         if ((fin.ne.1.d0).and.(jden.eq.'B')) fi0=dmin1(1.d0,fidhcc/dh0)
         lmod = 'DW'
c     Compute local emissivity in vectors EMIDIF (diffuse field) and
c     EMILIN (line field) used by routines TOTPHOT and NEWDIF
         call localem(te0, de0, dh0)
c     Compute mean intensity of radiation in vector TPHOT using
c     source vector SOUPHO
         call totphot(te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
c
         if ((jeq .eq. 'E').or.(jeq .eq. 'P').and.(jthm.eq.'S')) then
            dton0 = 1.d33
            call teequi(te0, te0, de0, dh0, dton0, nmod)
         else if ((jeq .eq. 'C') .or. (jthm .eq. 'T')) then
c     Equilibrium structure (C=collisional ionization equilibrium)
c     Thermal structure (T=isothermal)

c     Compute ionization equilibrium at initial temperature and 
c     density for all elements. Outputs are electron density DE 
c     and the fractional abundances of the different species of 
c     each element POP(6,11)
            call equion(te0, de0, dh0)
c
            call localem(te0,de0,dh0)
            call totphot(te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
            call equion(te0, de0, dh0)
c
            call localem(te0,de0,dh0)
            call totphot(te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
            call equion(te0, de0, dh0)
c
         else
            distim(1,0) = dis0
            distim(2,0) = dis0/cls
            treach0 = distim(2,0)
            dton0 = dmax1(0.d0,dmin1(tlife,telap-((2.d0*dis0)/cls)))
            jjd = jden
            jden = 'C'
            exl = -1.d0
            tex = 0.95d0*tm00
            prescc = dhn*rkb*1.d4+fradpress(0.0d0,dhn)
c
            call copinto(pop0, pop)
            call evoltem(tm00,te0,de0,dh0,prescc,dton0,exl,tex,lut0)
c
            jden = jjd
c
            if (dton0.le.0.d0) then
              write(*, 305) dton0
 305          format(/,' Age was incorrectly set, it is shorter than',/,
     &            ' needed at the first space step , DTon0 :',1pg10.3,/)
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
            dh0 = frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac
     &                ,efac,ffac)
         else
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c

            frp    = fradpress(dr,dhn)
            prescc = dhn*rkb*1.d4 + radPint + frp
            ponk   = (prescc/rkb)
            pres0  = fpressu(te0,dhn,pop)
            dh0    = dhn *prescc/(pres0)
         end if
c
         dhl = dif(dhp,dh0)
         write(*, 472) te0, dhl
         if (dhl.ge.dhlma) goto 12
         tspr = dton0
c
         tempre = te0
         call copinto(pop, poppre)
         call copinto(pop, ppre)
c
c*****************************************************
c      write(*,*)'END OF INNER BOUNDRY FOR M = 1'
c*****************************************************
c
      else !m ne 1
c
c     calculate qhdn into zone
c
        q1 = 0.0d0
        q2 = 0.0d0
        q3 = 0.0d0
        q4 = 0.0d0
c
c     Integrate local mean intensities to derive number of H/He 
c     ionizing photons. Returns q1=QAH1, q2=QAHEI, q3=QAHEII, 
c     q4=QATOT.
        call intvec(tphot, q1, q2, q3, q4)
c
        qhdin=q4/(dh1*zen)
c
      end if
c
c*****************************************************************
c
c      write(*,*)'***ESTIMATE OF VALUES AT THE OUTER BOUNDARY'
c
c*****************************************************************
c
      call copinto(poppre, pop)
      call copinto(poppre, popend)
c
      if (m.lt.4) then
c
c     linear extrap from previous 2 steps
c
         te1 = te0+(((te0-tep)*te0)/tep)
         te1 = dmax1(te0*agmu,dmin1(te0/agmu,te1))
         te1a = te1
      else
         te1a = te1
         te1 = thist(3)
      end if
c
      if (jthm.eq.'T') te1 = fixtemp
c
      if (te1.lt.0.d0) te1 = 290.d0
      if (jden .eq. 'C') then
         dh1 = dhn
      else if (jden .eq. 'F') then
         dh1 = frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac
     &             ,efac,ffac)
      else
c     error if pop is very different from equilibrium at te1
c         frp    = fradpress
c         dh1 = (dhn*prescc)/(fpressu(te1,dhn,pop)+frp)
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
            frp    = fradpress(dr,dh)
            prescc = dhn*rkb*1.d4 + radPint + frp
            ponk   = (prescc/rkb)
            pres1  = fpressu(te1,dhn,pop)
            dh1    = dhn *prescc/(pres1)
      end if
      t = (te1+te0)/2.d0
      dh = (dh0+dh1)/2.d0
c
c***********************************************
c      write(*,*)'***DETERMINE SPACE STEP : DR'
c***********************************************
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
      mmod = 'DIS'
      lmod = 'SO'
c
c     get new ionisation at t predicted
c     equion calculates a new de
c     
      call equion(t,de,dh)
c
      call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
c     Compute distance DR to obtain a given total photon absorption 
c     fraction DTAU
      call absdis(t, dh, fi, dtau, dr, rad0, pop)

 30   continue                  ! Loop to achieve temperature convergence

      if (jgeo .eq. 'S') then
         if (jden .eq. 'F') then
            if ((frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac
     &           ,efac,ffac)).gt.ofac) then
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
            dvoluni = ftpi*((fring(abc0,abc9,rad1)*((rad1/runit)**3))
     &      -(fring(abc0,abc9,rad0)*((rad0/runit) ** 3)))
         else
            dvoluni = ftpi*(((rad1/runit)**3)-((rad0/runit)**3))
         end if
      endif
      if (jgeo.eq. 'P') then
         dis1 = dis0+dr
         rad1 = 0.d0
         wdil1 = wdpl
         dvoluni = dr/runit
      end if
      if (jgeo.eq. 'F') then
         dr = 1.d15
c        dtau = dtau0
         call absdis(t, dh, fi, dtau, dr, rad0, pop)
         dis1 = dis0+dr
         rad1 = 0.d0
         wdil1 = fdilu(rstar,dis1)
         dvoluni = dr/runit
      end if
c
c**********************************************
c       write(*,*)'*** dr done',dr
c**********************************************
c



c******************************************************************
c
c      write(*,*)'CALCULATION OF IONIC POPULATION AT THE OUTER'
c     of current step
c
c******************************************************************
      if (jden .eq. 'C') then
         dh1 = dhn
      else if (jden .eq. 'F') then
         dh1 = frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac
     &             ,efac,ffac)
      else
c         frp    = fradpress
c         dh1 = dhn*prescc/(fpressu(te1,dhn,popend)+frp)
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
c        write(*,'(" P Pr",2(1pg12.5))') fpressu(te1,dhn,popend),frp

            frp    = fradpress(dr,dh)
            prescc = dhn*rkb*1.d4 + radPint + frp
            ponk   = (prescc/rkb)
            pres1  = fpressu(te1,dhn,popend)
            dh1    = dhn *prescc/(pres1)
      end if
      de1 = feldens(dh1,popend)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi1=dmin1(1.d0,fidhcc/dh1)
      rww = 1.d0
      if (jeq.eq.'C') te1 = frad(dis1,tofac,txfac,tpfac,tafac,
     & 				tbfac,tcfac,tscalen,0.d0,1.d0,0.d0)
      t = (te0+(rww*te1))/(1.d0+rww)
      dh = (dh0+(rww*dh1))/(1.d0+rww)
      de = (de0+(rww*de1))/(1.d0+rww)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
c
c  Determine dust continuum emission for this region
c
c      if ((grainmode).AND.(IRmode.ne.0)) then
c      lmod = 'ALL'
c       call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)
c       call dusttemp(t,dh,de, fi, dr,m,phot0)
c      endif
c
c
      lmod = 'DW'
      wei = 1.d0/(1.d0+rww)
c     Form weighted average of POPPRE and POPEND in POP
      call averinto(wei, poppre, popend, pop)
      call localem(t, de, dh)
      call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)

c     Calculate effective density of photons per particle
      call zetaeff(fi, dh)
c
c*****************************************************************
c      write(*,*)'*EQUILIBRIUM IONISATION AND TEMPERATURE .'
c*****************************************************************
c
      if ((jeq .eq. 'E').or.((jeq .eq. 'P').and.(jthm.eq.'S'))) then
         dton1 = 1.d33
         call copinto(popend, pop)
         call teequi(te1, tef, de1, dh1, dton1, nmod)
      else if ((jeq.eq.'C').or.(jthm.eq.'T')) then
c
         te1 = frad(dis1,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen,
     &    	    0.d0,1.d0,0.d0)
         call copinto(popend, pop)
c
         call equion(te1, de1, dh1)
c
         i = 0
c
c     PJM Aug 2008: Changed IF/GOTO loop to DO WHILE. Loop was testing 
c     value of DIFT that was not set by DIFPOP. It now triggers on the 
c     value of DIFP that is set by DIFPOP. It is also unnecessary to call
c     COPINTO(POP,POPZ) at the beginning and the end of the loop.
         i = 0
         treap = 1.d-12
         difp = 1.d0
         do while ((difp .ge. 0.001) .and. (i .le. 4))
 
c     Store current ionic populations POP in POPZ.
            call copinto(pop,popz)
c     Compute local emissivity in vectors EMIDIF and EMILIN.
            call localem(te1,de1,dh1)
c     Compute mean intensity of radiation in vector TPHOT.
            call totphot(te1,dh1,fi,rad,dr,dv,wdil,lmod)
c     Calculate effective density of photons per particle.
            call zetaeff(fi,dh1)
            call equion(te1,de1,dh1)
c     Determine average logarithmic change in ionic populations.
            call difpop(pop,popz,treap,atypes,difp)

            i = i+1

         end do
c
         tef = te1
c         
      else
c
c     Finite age, non-equilibrium ionisation. Determine time step.
c     Takes into account ionization front velocity.
c
         distim(1,m) = dis1
         distim(2,m) = distim(2,m-1)+(dr/(viofr+(1.d-35*dr)))
         tprop = distim(2,m)
         dtauon = 1.0d0
         mmod = 'DIS'
c
         call absdis(t, dh, fi, dtauon, disin, 0.0d0,pop0)
         fthick = dmin1(disin,dis1-distim(1,0))
         dison = dis1-fthick
         protim = distim(2,0)
         do 330 j = m-1, 0, -1
            if (distim(1,j).gt.dison) goto 330
            prf = (dison-distim(1,j))/((distim(1,j+1)-distim(1,j))
     &            +epsilon)
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
c     Compute equilibrium temperature and ionization state of the gas
            call teequi(te1, tef, de1, dh1, dton1, nmod)
         end if
         dtco1 = (1.5d0*fpressu(tef,dh1,pop))/(eloss+epsilon)
      end if
c
c*****************************************************************
c       write(*,*)'END EQUILIBRIUM IONISATION AND TEMPERATURE .'
c       write(*,*) t,te1,tef
c*****************************************************************
c
      if (jthm.eq.'T') tef = fixtemp

c
      if (jden .eq. 'C') then
         dhf = dhn
      else if (jden .eq. 'F') then
         dhf = frad(dis1,dhn,xfac,pfac,afac, bfac,cfac,scalen,vfac
     &             ,efac,ffac)
      else
c         frp    = fradpress
c         dhf = dhn*prescc/(fpressu(tef,dhn,pop)+frp)
c         write(*,'(" P Pr",2(1pg12.5))') fpressu(tef,dhn,pop),frp
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
            frp    = fradpress(dr,dh)
            prescc = dhn*rkb*1.d4 + radPint + frp
            ponk   = (prescc/rkb)
            pres1  = fpressu(tef,dhn,pop)
            dhf    = dhn *prescc/(pres1)
c            write(*,'(" P/k Pr/k dPr/k",3(1pg12.5))') 
c     &               ponk, radPint/rkb, frp/rkb
      end if
      def = feldens(dhf,pop)

c
c************************************************************
c   write(*,*)'***COMPARES END VALUES WITH PREVIOUS STEP'
c************************************************************
c
      call difpop(pop, popend, trea, atypes, dift)
      call difpop(pop, poppre, trea, atypes, difp0)
      call copinto(pop, popend)
c
      drxp = dr
      dtl = dif(tef,te1)
      dhl = dif(dh1,dhf)
c
c     strict criteria
c
c     difg = dmax1(dift/difma,dtl/dtlma,dhl/dhlma)
c
c     relaxed cirteria
c     
c     difg = dhl/dhlma
c     
c     lenient
c     
      difg = dmax1(dtl/dtlma,dhl/dhlma)
c     
c     no criteria
c
c      difg = 0.d0
c     
c
      niter = niter+1
      if (difg.le.difmi) then
         texm = te1
         dtaux = dtau
         difmi = difg
      end if
c
      write(*, 472) tef, dhl, dtl, dift, difg, difp0, dtau,dr
  472 format(t7,0pf8.0,5(0pf8.3),3(1pg9.2))
c
c***************************************************************
c  write(*,*)'AVERAGE QUANTITIES FOR THE SPACE STEP CONSIDERED'
c***************************************************************
c
      rww = 1.d0
      t = (te0+(rww*tef))/(1.d0+rww)
      dh = (dh0+(rww*dhf))/(1.d0+rww)
      de = (de0+(rww*def))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
      disav = (dis0+(rww*dis1))/(1.d0+rww)
      wei = 1.d0/(1.d0+rww)
      call averinto(wei, poppre, popend, pop)
c
      if (jgeo .eq. 'S') then
         rad = disav
         wdil = fdilu(rstar,rad)
      endif
      if (jgeo .eq. 'F') then
         rad = 0.d0
         wdil = fdilu(rstar,disav)
      endif
      if (jgeo .eq. 'P') then
         wdil = wdpl
         rad = 0.d0
      end if
c
c****************************************************************
c  write(*,*)'CHECK CONVERGENCE OF END VALUES WITH PREVIOUS STEP'
c****************************************************************
c
      temp1 = ((de0/dh0)-dedhmi)/dedhma  !ionization vs 1st level ionization
      temp2 = (de1/(de0+1.d-7))*((dh0/dhf) ** 2) !change in ionization state
      if ((difg.ge.1.d0).and.(qhdh.gt.1d3)) then
         weiph = 0.5d0
         rmm = 1.d5
         if ((niter.gt.11).and.(te1 .eq. texm)) then
c     No convergence of T or  H dens after 12 iterations
            write(*, 172) dhl, dtl, dift, difg
 172        format(' NO SPATIAL CONVERGENCE :',4(0pf8.3))
            goto 31
         else if (niter .eq. 12) then
c     Convergence of T and H dens after 12 iterations
            te1 = texm
            dtau = dtaux
            difmi = 0.d0
         else if (((temp1.gt.0.01d0).and.(temp2.lt.0.2d0)).or.(dhl
     &      .ge.0.07d0)) then
c     changes in ionization state or Hdens are big
            agmu = dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-((1.d0-agmu)*
     &      0.6d0))
            nidh = nidh+1
            if (nidh.lt.3) niter = niter-1
            te1 = (te0+dmax1(agmu*te0,dmin1(te0/agmu,te1)))/2.d0
            weiph = 0.85d0
            rmm = 0.25d0
            dtau = (0.2d0*rmm)*dtau
         else if (((niter.eq.3).or.(niter.eq.6)).and.((jden.ne. 
     &      'B').or.(difg .eq. (dift/difma)))) then
c at 3 or 6 iterations and not isobaric or convergence in T slower than H dens 
            agmu = dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-
     &      ((1.d0-agmu)*0.65d0))
            nidh = max0(0,min0(nidh-2,idnint(1.0d0+dint(niter/3.0d0))))
            te1 = (((te1+te0)+te1a)+te0)/4.0
            te1 = dmax1(te0*agmu,dmin1(te0/agmu,te1))
            weiph = 0.75d0
            rmm = 0.4d0
            dtau = (0.2d0*rmm)*dtau
         else
c  None of the above criteria hold
            te1 = tef
            weiph = 0.85d0
            rmm = 0.25d0
            dtau = (0.2d0*rmm)*dtau
         endif
         
 34      dtau = dtau/2.d0       ! Loop back to here if dr > rmm*drxp

         mmod = 'DIS'
         lmod = 'SO'
         drzi = dr
         drzf = dr
         call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
         call absdis(t, dh, fi, dtau, drzi, rad0, poppre)
         call absdis(t, dh, fi, dtau, drzf, rad0, popend)
c
c
         dr = dexp((weiph*dlog(drzi+.1d0))+((1.d0-weiph)
     &                *dlog(drzf+.1d0)))
         if ((dr.gt.(rmm*drxp)).and.(dtau.gt.1.d-16)) goto 34
         if (dr.lt.drp*1.d-2) then
            dr = drp*1.d-2
            niter = niter+11
            te1 = texm
         endif
         if (dif(dr,drxp).le.0.01d0) then
            niter = niter+11
            te1 = texm
         endif
         wei = dmin1(1.d0,(1.2d0*dr)/drxp)
         call averinto(wei, popend, poppre, popend)
c         
c     Loop back until space step (???) converges.        
         goto 30
c
c     Jump here if no spatial convergence.        
 31      continue

      else if ((difg.ge.1.d0).and.(qhdh.gt.1d4)) then
c     THIS IS THE SAME TEST AS FOR THE MAIN TEST ABOVE!!!????
         dtau=dtau/4.d0
      end if
c
c
      te1 = tef
      dlos1 = dlos
      telast = t
      frx = 1.d0-pop(jpoen,ielen)
c
c*****************************************************************
c      write(*,*)'END CHECK CONVERGENCE STEP'
c*****************************************************************
c
c****************************************************************
c      write(*,*)'INTEGRATES COLUMN DENSITY AND END DIFFUSE FIELD'
c****************************************************************
      if (jgeo .eq. 'S') then
         jmod = 'OUTW'
         frdw = 1.d0
      endif
      if (jgeo .eq. 'P') then
c         jmod = 'LOUP'
c         frdw = 0.0d0
          jmod = 'OUTW'
          frdw = 1.0d0
      end if
      imod = 'COLD'
c
c
c
c
c  Determine dust continuum emission for this region
c
      lmod = 'DW'

      if ((grainmode).AND.(IRmode.ne.0)) then
       call localem(t, de, dh)
       call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)
       call dusttemp(t,dh,de, fi, dr,m)
      endif
c
c
      call localem(t, de, dh)
c      for diagnostic purposes
c      call totphot(t, dh, fi, rad0, dr, dv, wdil1, lmod)
c 
      call zetaeff(fi, dh)
      call newdif(t, t, dh, fi, rad0, dr, dv, 0.d0, dv, frdw, jmod)
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
c
c*****************************************************************
c      write(*,*)'CALCULATION OF AVERAGE QUANTITIES '
c       EQUAL VOLUMES
c*****************************************************************
c
      if (jgeo .eq. 'S') then
         rww = (1.d0+((0.75d0*dr)/(rad0+dr))) ** 2
      else
         rww = 1.d0
      end if
c
      t      = (te0+(rww*tef))/(1.d0+rww)
      dh     = (dh0+(rww*dhf))/(1.d0+rww)
      de     = (de0+(rww*def))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav  = (dis0+(rww*dis1))/(1.d0+rww)
      dlosav = (dlos0+(rww*dlos1))/(1.d0+rww)
      dton   = (dton0+(rww*dton1))/(1.d0+rww)
      treach = (treach0+(rww*treach1))/(1.d0+rww)
      dtco   = (dtco0+(rww*dtco1))/(1.d0+rww)
      wei    = 1.d0/(1.d0+rww)
c
      call averinto(wei, poppre, popend, pop)
c
	  frp     = fradpress(dr,dh)
	  radPint = radPint + frp
c
      if (jgeo .eq. 'S') then
         rad = disav
         wdil = fdilu(rstar,rad)
      endif
      if (jgeo .eq. 'F') then
         rad = 0.d0
         wdil = fdilu(rstar,disav)
      endif
      if (jgeo .eq. 'P') then
         wdil = wdpl
         rad = 0.d0
      end if
c
c*************************************************************
c    ***CALCULATION OF IONIC POPULATION AND TEMPERATURE WHEN
c       OR IF SOURCE SWITCHED OFF
c*************************************************************
c
      dtoff = dmax1(0.d0,(telap-((2.d0*disav)/cls))-tlife)
      if ((jeq .ne. 'E').and.(jeq.ne.'C').and.(dtoff.gt.0.0)) then
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
c      write(*,*)'COMPUTES SPECTRUM OF THE REGION CONSIDERED'
c
c*********************************************************
c
      imod = 'REST'
c     Calculate total cooling rate of the plasma
      call cool(t, de, dh)
      
c     Integrate the line and temperature data and column density
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
c
c     record line ratios
c     
c         hoii(m) = (fbri(1,ox2)+fbri(2,ox1))/(hbri(2)+epsilon)
c
      if (ox3.ne.0) then
         hoiii(m) = (f6bri(7,ox3)+f6bri(10,ox3))/(hbri(2)+epsilon)
      endif
      if (ox2.ne.0) then
         hoii(m) = (fbri(3,ox1))/(hbri(2)+epsilon)
      endif
      if (ni2.ne.0) then
         hnii(m) = (f6bri(7,ni2)+f6bri(10,ni2))/(hbri(2)+epsilon)
      endif
      if (su2.ne.0) then
         hsii(m) = fbri(1,su2)+fbri(2,su2)/(hbri(2)+epsilon)
      endif
      if (ox1.ne.0) then
        if (ispo .eq. ' OI') hsii(m) = fbri(3,ox1)/(hbri(2)+epsilon)
      endif
c
c    Test for PAH existence (if grainmode true)
c
      if (grainmode.AND.(pahfrac.gt.0.d0)) then
c
c  Calculate habing photodissociation parameter
c      
        Ffuv=0.d0
        do i=fuvmin,fuvmax-1
           dnu=(ephot(i+1)-ephot(i))*evplk
           eng=0.5*(ephot(i+1)+ephot(i))*eV
c    photons s-1 cm-2
           Ffuv=Ffuv+4*pi*tphot(i)*dnu/eng
        enddo
        Hab=Ffuv/(cls*dh)
c
c Calcualte PAH ionization/Local ISRF ratio (Q factor)
c
	pahsum=0.d0
	Ffuv=0.d0
	do i=fuvmin,pahabsmax-1
          dnu=(ephot(i+1)-ephot(i))*evplk
          eng=0.5*(ephot(i+1)+ephot(i))*eV
	  pahsum=pahsum+pahnabs(i)*dnu	  
	  Ffuv=Ffuv+4*pi*tphot(i)/eng*pahnabs(i)*dnu
	enddo
	pahQion=Ffuv/pahsum	
	pahQfac=pahQion/pahQ0
c
        write(*,210) Hab,pahQfac
 210	format('Habing parameter =',1pg12.5,', PAH ion/ISRF =',1pg12.5)
c
c Determine if PAHs exist
c If so, deplete gas of Carbon in PAHs
c
	if (.NOT.pahmode) then
	  if (pahend.eq.'H') then !Habing
	    if (hab.lt.pahlimit) then
	      pahmode=.TRUE.
	    endif
	  elseif (pahend.eq.'Q') then !QHDH
	    if (qhdh.lt.pahlimit) then
	      pahmode=.TRUE.
	    endif
	  elseif (pahend.eq.'I') then !Q factor
	    if (pahQfac.lt.pahlimit) then
	      pahmode=.TRUE.
	    endif
	  endif
	  if(pahmode) then
	    if (ClinPAH) then !C linked to PAHs
	       zion(zmap(6)) = zion0(zmap(6))*dion(zmap(6))
	    else
	        zion(zmap(6))=zion(zmap(6)) -
     &	      		   zion0(zmap(6))*(1.d0-dion(zmap(6)))*pahCfrac
            endif
	  endif
	endif
      endif
c
c
c    ***OUTPUT QUANTITIES RELATED TO SPACE STEP
c
c
c     write balance for first step or if poll file exists
c
c     
      fhi  = pop(1,1)
      fhii = pop(2,1)
c
      pollfile = 'balance'
      inquire(file=pollfile, exist=iexi)
      if (((jbal.eq.'YES').and.(m.eq.1)).or.(iexi)) then
         caller = 'P5'
         pfx = 'p4bal'
         pfx = jbfx
         np = 5
         wmod = 'NONE'
         call wbal(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0)
      endif
c
      pollfile = 'photons'
      inquire(file=pollfile, exist=iexi)
      if ((iexi)) then
         caller = 'P5'
         pfx = 'psou'
         np = 4
         wmod = 'REAL'
         dva = 0
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
      endif
c
      pollfile = 'nebphot'
      inquire(file=pollfile, exist=iexi)
      if ((iexi)) then
         lmod='NEBL'
         call totphot(t, dh, fi, rad, dr, dv, wdil1, lmod)
         caller = 'P5'
         pfx = 'nsou'
         np = 4
         wmod = 'REAL'
         dva = 0
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
         lmod='ALL'
         call totphot(t, dh, fi, rad, dr, dv, wdil1,lmod)
      endif
c
      pollfile = 'sophot'
      inquire(file=pollfile, exist=iexi)
      if ((iexi)) then
         lmod='SO'
         call totphot(t, dh, fi, rad, dr, dv, wdil1,lmod)
         caller = 'P5'
         pfx = 'ssou'
         np = 4
         wmod = 'REAL'
         dva = 0
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
         lmod='ALL'
         call totphot(t, dh, fi, rad, dr, dv, wdil1,lmod)
      endif
c
      pollfile = 'loclphot'
      inquire(file=pollfile, exist=iexi)
      if ((iexi)) then
         lmod='LOCL'
         call totphot(t, dh, fi, rad, dr, dv, wdil1,lmod)
         caller = 'P5'
         pfx = 'lsou'
         np = 4
         wmod = 'REAL'
         dva = 0
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
         lmod='ALL'
         call totphot(t, dh, fi, rad, dr, dv, wdil1,lmod)
      endif
c
      pollfile = 'IRphot'
      inquire(file=pollfile, exist=iexi)
      if ((iexi).AND.(IRmode.gt.0)) then
         caller = 'P5'
         pfx = 'irsou'
         np = 5
         wmod = 'REAL'
         dva = 0
         do i=1,infph-1
            pahtest(i)=0.5d0*(ephot(i+1)+ephot(i))*eV*IRphot(i)*dr
         enddo
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,pahtest)
      endif
c
      pollfile = 'pahphot'
      inquire(file=pollfile, exist=iexi)
      if ((iexi).AND.(pahmode)) then
         caller = 'P5'
         pfx = 'pahs'
         np = 4
         wmod = 'REAL'
         dva = 0
         do i=1,infph
            pahtest(i)=paheng*pahflux(i)
         enddo
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,pahtest)
      endif
      
      pollfile = 'speclocal'
      inquire(file=pollfile, exist=iexi)
      if (iexi) then
         caller = 'P5'
         pfx = 'local'
         np = 5
      	 sfx = 'ph5'
         call newfile(pfx,np,sfx,3,fspl)
      
         open(lusp, file=fspl, status='NEW')

         linemod = 'LAMB'
         spmod   = 'REL'
         call speclocal(lusp
     &   ,tloss,eloss,egain,dlos
     &   ,t,dh,de, pop(1,1)
     &   ,disav,dr
     &   ,linemod,spmod)
         
         close(lusp) 
         
      endif
c
c   Infrared output (every 1 step)
c
      pollfile='IRlocal'
      inquire(file=pollfile, exist=iexi)
       if ((IRmode.ne.0).AND.iexi) then
        if (IRfile.EQ.' ') then
           caller = 'P5'
           np = 6
           pfx = 'IRflux'
           sfx = 'sou'
           call newfile(pfx,np,sfx,3,IRfile)
           open(IR1, file=IRfile, status='NEW')
           write(IR1,2500)
 2500      format('%  Infrared flux per region ',/,
     &            '%  MAPPINGS III v1.2.0md',
     &            '%  given as:',/,
     &            '%    energy edge (eV) ,,continuumflux, ', 
     &            'IRflux Fnu(erg s-1 cm-2 Hz-1Sr-1)',/,
     &            '%')
           close(IR1)
         endif
         do i=1,infph-1
            cnphot(i)=0.d0
         enddo
         call freefree(t,de,dh)
         call freebound(t,de,dh)
         call twophoton(t,de,dh)
         open(IR1, file=IRfile, status='OLD',access='APPEND')
         write(IR1,2501) m
         write(IR1,*) infph
 2501    format(/,' Region ',i4.4)
         do i=1,infph
          write(IR1,2502) ephot(i),cnphot(i)*dr,IRphot(i)*dr
         enddo
 2502    format(1pg14.7,' ',1pg14.7,' ',1pg14.7)
         close(IR1)
c         do i=1,mxinfph
c            storephot(i)=tphot(i)
c         enddo
c         call totphot(t, dh, fi, rad, dr, dv, wdil, 'LOCL')
c         call wpsoufile(caller,IRneb,
c     &        t,de,dh,
c     &        0.d0,0.d0,qhdh,
c     &        disav,dr,
c     &        0.d0,0.d0,
c     &        1.d0,tphot)
c         do i=1,mxinfph
c            tphot(i)=storephot(i)
c         enddo        
      endif
c
c     graincharge output
c
      if (grainmode) then
        pollfile = 'graincharge'
        inquire(file=pollfile, exist=iexi)
        if (iexi) then
          if (chargefile.eq.' ') then
           caller = 'P5'
           np = 5
           pfx = 'grpot'
           sfx = 'ph5'
           call newfile(pfx,np,sfx,3,chargefile)
           open(IR1, file=chargefile, status='NEW')
           write(IR1,2510)
 2510      format('%  grain charge in each region ',/,
     &            '%  MAPPINGS III v1.2.0md',
     &            '%  given as:',/,
     &            '%  m,Av. distance,Temp,dh,de,',/,
     &            '%  grain charge(allsizes) (graphite),',/,
     &            '%  grain charge(allsizes) silicate',/,
     &            '%')
           write(IR1,2511) (grainrad(i),i=mindust(1),maxdust(1))
 2511      format('gra radii',10(1pg11.4))
           write(IR1,*)
           write(IR1,2512) (grainrad(i),i=mindust(2),maxdust(2))
 2512      format('sil radii',10(1pg11.4))
           close(IR1)
          endif
          open(IR1, file=chargefile,status='OLD',access='APPEND')
           write(IR1,2515) m,disav,t,dh,de
          write(IR1,2516)(grainpot(i,1),i=mindust(1),maxdust(1))
          write(IR1,2517)(grainpot(i,2),i=mindust(2),maxdust(2))
          close(IR1)
 2515      format(/,i3,4(1pg14.5))
 2516      format('gra ',30(1pg11.4))
 2517      format('sil ',30(1pg11.4))
        endif 
      endif

c
      open(luop, file=filna, status='OLD', access='APPEND')
      open(lupf, file=filnb, status='OLD', access='APPEND')
      open(lups, file=filnc, status='OLD', access='APPEND')
      if (jall.eq.'YES') then
        open(lusl, file=newfil, status='OLD', access='APPEND')
      endif
c
  461 format(i4,1pg11.4,1x,1pg10.3,5(1pg11.4),2(1pg10.3),2(1pg11.4)
     &     ,2(1pg10.3))
 1461 format(t4,10(1pg10.3))
      write(luop, 461) m, t, dlosav, disav, dis1, dr, dh, de,fhi, fhii,
     &       qhdin, qhdh,caseab(1), caseab(2)
c
c      write(luop, 461) m, t, dlosav, dis1, disav, dr, fi, dh, de
c      write(luop, 1461) fhi, fhii, zetae, qhdh,hoiii(m),hoii(m),
c     &                  hnii(m),hsii(m),hbeta,caseab(1)
c
      if (jeq .eq. 'F') write(luop, 462) fthick, treach, dtco, 
     &  dton, dtoff
  462 format(t71,5(1pg9.2))
c
  473 format(t98,5(1pg9.2))
      if (jeq .eq. 'P') write(luop, 473) recscal, dtoff
c
 465  format('    Case A-B (HI, HeII): ',2(1pg11.3))
      write(*,465) caseab(1), caseab(2)

 460  format(' #',t7,'Teav',t15,'DLOSav',t24,'Rfin',t33,'dr',t42,
     &'nHav',t51,'FHI',t60,'QHDNf',t69,'[OIII]',/
     & i4,0pf8.0,1pg9.1,6(1pg9.2))
c
      write(*, 460) m,t,dlosav,dis1,dr,dh,fhi,qhdh,hoiii(m)
c
c
c
      if (jiel.eq.'YES') then
       do i = 1 ,ieln 
        open(luions(i),file=filn(i),status='OLD',access='APPEND')
c
 3000   format(29(1pg14.7,1x))
        if (i.eq.1) then
          write(luions(i),3000) disav,dr,t,
     &        (pop(j,iel1), j = 1,maxion(iel1))
        endif
        if (i.eq.2) then
          write(luions(i),3000) disav,dr,t,
     &        (pop(j,iel2), j = 1,maxion(iel2))
        endif
        if (i.eq.3) then
          write(luions(i),3000) disav,dr,t,
     &        (pop(j,iel3), j = 1,maxion(iel3))
        endif
        if (i.eq.4) then
          write(luions(i),3000) disav,dr,t,
     &        (pop(j,iel4), j = 1,maxion(iel4))
        endif
c
        close(luions(i))
       enddo
c
c     end .tsr files
      endif
c
c
c      nt1 = time()
c
c      dtarr = dtime(tarray)
c
c 100  format('Step time (map3,sys): ',2(f8.1,1x),/) 
c      write(*,100) tarray(1),tarray(2)
c
c      nt0 = nt1
c
      if (jall.eq.'YES') then
         write(lusl, 1250)
         write(lusl, 678) m, t, de, dh, fi, 
     &              dvoluni, dr, disav, (disav-remp)
 1250    format(//' Step   Te Ave.(K)   ',
     &                '  ne(cm^-3)   ',
     &                '  nH(cm^-3)   ',
     &                '  fill. Fact. ',
     &                '  dVol Unit.  ', 
     &                '    dR (cm)   ',
     &                ' Dist.Ave.(cm)',
     &                ' D - Remp.(cm)')
 678     format(i3,8(1pg14.7)/)
      endif
c
      wmod = 'PROP'
      call wmodel(lupf,t,de,dh,disav,wmod)
      wmod = 'LOSS'
      call wmodel(lups,t,de,dh,disav,wmod)
c
c     experiment with timescales
c
      ctime = fcolltim(de)
      rtime = frectim2(de)
      ptime = fphotim()
c
      eqtime = (1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
      eqtime = dabs(1.d0/eqtime)
c
c     experiment with timescales
c
      stime = (1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      stime = dabs(1.d0/stime)
c
 700  format(' Timescales: ',5(1pg11.4,1x))
      write(*,700) ctime,rtime,ptime,eqtime,stime
c
      if (jall.eq.'YES') then
        call wionabal(lusl,pop)
        call wionabal(lusl,popint)
      endif
c
      close(luop) 
      close(lupf) 
      if (jall.eq.'YES') close(lusl) 
      close(lups) 
c
c     ***END OUTPUT
c
c
c**********************************************************************
c
c       CHECKPOINT BATCHMAP
c
      pollfile = 'save'
      inquire(file=pollfile, exist=iexi)
      if (iexi) then
         call dumpp4
      endif
c
c
c
c**********************************************************************
c
c
c     ***TEST ENDING CONDITIONS
c
      do 513 n = 1, ionum
      if ((atpho(n).eq.ielen).and.
     &(ionpho(n).eq.jpoen)) goto 517
c
c  Note: loop out of bound for Column measurements of HII
c
  513 continue
  517 taux = sigpho(n)*popint(jpoen,ielen)
      if (((jend.eq.'A').or.(jend.eq.'B')).and.(frx.le.fren)) 
     &goto 505
      if ((jend.eq.'C').and.(telast.le.tend)) goto 505
      if ((jend.eq.'D').and.(taux.ge.tauen)) goto 505
      if ((jend.eq.'E').and.(dis1.ge.diend)) goto 505
      if (jend.eq.'F') then
       if (jpoen.gt.maxion(ielen)) then
        popinttot=0
        do i=1,maxion(ielen)
          popinttot=popinttot+popint(i,ielen)
        enddo
        if (popinttot.ge.tauen) goto 505
       else
         if (popint(jpoen,ielen).ge.tauen) goto 505
       endif
      endif
c
c Check Visual extinction condition or end when HII<1%
c

      if (jend.eq.'H') then
         A_v = dustAv*dustint
         if ((A_v.gt.A_vend).OR.(de/dh.lt.1e-5)) goto 505
      endif
c
      if (((((de1/dh1)-dedhmi)/dedhma).le.exla).and.(.NOT.grainmode))
     &   goto 505
      pollfile = 'terminate'
      inquire(file=pollfile, exist=iexi)
      if (iexi) goto 505
c
c
c    ***RESET INNER BOUNDARY QUANTITIES FOR NEXT SPACE STEP
c
      call copinto(popend, poppre)
c
      tep = te0
      te0 = tef

c     Record electron temperatures in THIST for quadratic extrapolation
c     when m > 3.

      thist(1) = thist(2)
      thist(2) = tep
      thist(3) = tef
c
      de0 = def
      dh0 = dhf
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
c      write(*,*)'LOOP BACK AND INCREMENT M, DO NEXT STEP.....'
c
c**********************************************************************
c
      m = m+1
      write(*,'(/)')

c
c     stop things getting silly....
c
      if (m.gt.mma) goto 505
c
c     loop back
c     
      goto 1
c
c
c
  505 continue
c
c********************************************************************
c
c     write(*,*)'MODEL ENDED ; OUTPUT RESULTS **************'
c
c********************************************************************
c
      open(luop, file=filna, status='OLD', access='APPEND')
      open(lupf, file=filnb, status='OLD', access='APPEND')
      write(luop, 530) 
  530 format(//' Model ended',t20,'Tfinal',t28,'DISfin',t38,'FHXF',t48,
     &'TAUXF',t58,'Ending')
      write(luop, 510) t, dis1, frx, taux, jend
  510 format(12h -----------,t18,0pf8.0,3(1pg10.3),4x,a4)
      if ((jeq .ne. 'E').and.(abc3.gt.0.0)) write(luop, 519
     &) abc3
  519 format(/' Final fractional luminosity of source after ',
     &'turn off :',1pg10.3)
      if ((abc9.lt.huge).and.(jgeo .eq. 'S')) write(luop, 
     &967) abc0, abc9
  967 format(/' NB. :::::::::::::::Integration through the line', 
     &'of sight for a ring aperture of radii :',2(1pg10.3))
      write(lupf, 2047) tempre, tspr
 2047 format(//t4,'Preionisation conditions for step#0 ',
     &'for all elements   (TEpr :',0pf8.0,' Time :',1pg9.2,' )  :'/)
      write(lupf, 2050) (elem(i),i = 1, atypes)
 2050 format(1h ,t8,16(4x,a2,4x)/)
      maxio = 0
      do j = 1,atypes
         if (maxio.le.maxion(j)) maxio = maxion(j)
      enddo
      do j = 1,maxio
         write(lupf, 2070) rom(j), (ppre(j,i),i = 1
     &        , atypes)
 2070    format(a7,t8,16(1pg10.3))
      enddo
c
      if (jsou.eq.'YES') then
c
	      caller = 'P4'
	      pfx = jpfx
	      np = 5
	      wmod = 'NFNU'
	      dva = 0
	      call wpsou(caller,pfx,np,wmod,
     &     t,de,dh,
     &     0.d0,0.d0,qhdh,
     &     disav,dr,
     &     0.d0,0.d0,
     &     1.d0,tphot)

c
c     Down stream photon field
c
	      caller = 'P4'
	      pfx = jpfx
	      np = 5
	      wmod = 'REAL'
	      dva = 0
	      call wpsou(caller,pfx,np,wmod,
     &     t,de,dh,
     &     0.d0,0.d0,qhdh,
     &     disav,dr,
     &     0.d0,0.d0,
     &     1.d0,tphot)
c     
      endif
c
c
      if (jspec.eq.'YES') then
c
c     write nebula spectrum in nuFnu (ergs/s/cm2/sr)
c     
c
          lmod = 'NEBL'
          call totphot(te1, dh1, fi, dis0, dr, 0.0d0, wdil0, lmod)
	      caller = 'P4'
	      pfx = jpfx
	      np = 5
	      wmod = 'NFNU'
	      dva = 0
	      call wpsou(caller,pfx,np,wmod,
     &     t,de,dh,
     &     0.d0,0.d0,qhdh,
     &     disav,dr,
     &     0.d0,0.d0,
     &     1.d0,tphot)
c
c     Down stream nebula only photon field
c
	      caller = 'P4'
	      pfx = jpfx
	      np = 5
	      wmod = 'REAL'
	      dva = 0
	      call wpsou(caller,pfx,np,wmod,
     &     t,de,dh,
     &     0.d0,0.d0,qhdh,
     &     disav,dr,
     &     0.d0,0.d0,
     &     1.d0,tphot)
c     
c
      endif
c
c    ***COMPUTE AVERAGE SPECTRUM AND WRITES IT IN FILE PHN
c
      close(lupf) 
c
      call avrdata
      call wrsppop(luop)
      close(luop) 
c  
      open(lusp, file=filnd, status='OLD', access='APPEND')
      linemod = 'LAMB'
      spmod   = 'REL'
      call spec2(lusp,linemod,spmod)
      close(lusp) 
c
      write(*, 2005) fnam
 2005 format(//' Output created &&&&&& File : ',a/)
c
      return 
c
      end
      
