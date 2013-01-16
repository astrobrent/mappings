cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       PHOTOIONISATION MODEL
c       DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c
c       RSS 1991/97
c
c       PHOTO 4: step length based on crossection weighted photon
c       absorption fraction, in absdis.  This removes the need for
c       most of the convergence tests and associated problems
c       in earlier models.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compph4(dhn, fin, banfil,difma,dtlma,dhlma)
c
      include 'cblocks.inc'
c
c
c
      double precision poppre(mxion, mxelem), popend(mxion, mxelem)
      double precision ppre(mxion, mxelem),difma,dtlma,dhlma
      double precision thist(3),steps(3),dpw,dxp
      double precision popz(mxion, mxelem),difp,treap
      double precision ctime,rtime,ptime,eqtime,stime
c
      double precision aadn,agmax,agmu,aper,aperi,apero
      double precision a,b,c,x,fx,fp,n0,r0
      double precision de,dedhma,dedhmi,deeq
      double precision dh,dhf,dhl,dhn,dhp
      double precision difg,difmi,difp0,dift
      double precision disin,dison,dlos0,dlos1,dlosav
      double precision dr,drmd,drmu,drp,drxp,drzf,drzi
      double precision dtau,dtauon,dtaux,dtco,dtco0,dtco1
      double precision dtl,dtoff,dton,dton0,dton1,durmax
      double precision dv,dva,dvoluni,eqscal,exl,exla
      double precision fhi,fhii,fi,fi0,fi1,fidhcc,fin,frdw,frx,fthick
      double precision prescc,prf,protim,rad,radi,recscal,rmm,rww
      double precision t,taux,te1a,tef,telast,temp1,temp2,tempre
      double precision tep,tex,texm,tii0,tprop,trea,runit 
      double precision treach,treach0,treach1,tspr
      double precision wd,wdil0,wdil1,wei,weiph
      double precision def
c
      integer*4 np,luop,lups,lusp,lusl,lupf
      integer*4 j,lut0,m,maxio,mma,n,nidh,niter
      integer*4 luions(4),i
c      integer*4 nt0,nt1,time
c      real tarray(2)
c      real dtarr,dtime
c      
      character jjd*4, pollfile*12,tab*4
      character newfil*20, filna*20, banfil*20, fnam*20, filnb*20
      character imod*4, jmod*4, lmod*4, mmod*4, nmod*4,wmod*4,ispo*4
      character linemod*4, spmod*4
      character filnc*20, filnd*20,filn(4)*20, fspl*20
      character pfx*8,caller*4,sfx*4
      
      logical iexi
c
c     External Functions
c
      double precision fcolltim,fcrit,fdilu,feldens
      double precision fphotim,fpressu,frectim,frectim2
c
c     internal functions
c
      double precision dif,frad,fsight,fring
c
      dif(a,b) = dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
      frad(x,n0,fx,fp,a,b,c,r0) = n0*(fx*dexp((x-a)/r0)
     &                            +fp*((x/a)**b))+c
      fsight(aper,radi) = 1.d0-(dcos(dasin(dmin1(1.d0,aper/(radi+
     & (1.d-38*aper)))))*(1.d0-(dmin1(1.d0,aper/(radi+(1.d-38*aper)))
     & **2)))
      fring(aperi,apero,radi) = dmax1(0.d0,fsight(apero,radi)-fsight(
     & aperi,radi))
c
      jspot = 'NO'
      jcon = 'YES'
      ispo = 'SII'
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
      drmu = 3.0d0
      drmd = 3.0d-2
      dton0 = 0.0d0
      dton1 = 0.0d0
      dtoff = 0.0d0
      steps(1) = -3
      steps(2) = -2
      steps(3) = -1
c
      luop = 21
      lups = 22
      lusp = 23
      lupf = 27
      lusl = 45
      lut0 = 0
c
      ieln = 4
      do i = 1,ieln
         luions(i) = 22+i
      enddo
c
c
c    ***DERIVE FILENAME : PHN**.PH4
c
      call p4fheader(newfil,banfil,fnam,
     &filna,filnb,filnc,filnd,filn,
     &luop,lupf,lusl,lusp,lups,luions,dhn,fin)
c
c
c*********************************************************************
c
c      write(*,*)'***COMPUTATION STARTS HERE'
c
c*********************************************************************
c
      fi = fin
      fi0 = fin
      fi1 = fin
      fidhcc = fin*dhn
      dv = 0.0d0
      runit = (rmax-remp)/300.d0
      dtau = dtau0/2.d0
      if (jgeo .eq. 'S') then
      vunilog = 3.0d0*dlog10(runit)
      else
      vunilog = dlog10(runit)
      end if
      if ((jeq .eq. 'E').or.(jeq .eq. 'P').or.(jeq.eq.'C')) then
      nmod  = 'EQUI'
      else
      nmod  = 'TIM'
      end if
      if (jden .ne. 'B') then
      agmax = 0.05d0
      aadn  = 2.5d0
      else
      agmax = 0.80d0
      aadn  = 5.0d0
      end if
      agmu  = agmax
c
c***********************************************
c
c      write(*,*)'***ITERATE, M = STEP NUMBER'
c
c***********************************************
c
      m = 1
c
 1    continue
c
      jspot = 'NO'
      niter = 0
      nidh  = 0
      difmi = 1.d30
      wei   = fcrit(dtau/dtau0,1.0d0)
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
            te0 = tofac*frad(remp,tofac,txfac,tpfac,tafac,tbfac,tcfac
     &            ,tscalen)
         else
            te0 = 1.d4
         endif
c
         if (jthm.eq.'T') te0 = Fixtemp
c
         dr = 0.d0
         drp = dr
         if (jden .ne. 'F') then
            dh0 = dhn
         else
            dh0 = frad(remp,dhn,xfac,pfac,afac,bfac,cfac,scalen)
         endif
c
         de0 = feldens(dh0,pop0)
         dedhmi = 0.d0
         dedhma = 0.d0
         do 127 j = 1, atypes
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
c
         if ((fin.ne.1.d0).and.(jden.eq.'B')) fi0=dmin1(1.d0,fidhcc/dh0)
         lmod = 'DW'
         call localem(te0,de0,dh0)
         call totphot(te0, dh0, fi0, rad0, dr, dv, wdil0, lmod)
c
         if ((jeq .eq. 'E').or.(jeq .eq. 'P').and.(jthm.eq.'S')) then
            dton0 = 1.d33
            call teequi(te0, te0, de0, dh0, dton0, nmod)
         else if((jeq.eq.'C').or.(jthm.eq.'T')) then
c
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
            prescc = dhn*rkb*1.d4
c
            call copinto(pop0, pop)
            call evoltem(tm00,te0,de0,dh0,prescc,dton0,exl,tex,lut0)
c
            jden = jjd
c
            if (dton0.le.0.d0) then
               write(*, 305) dton0
 305           format(/,' Age was incorrectly set, it is shorter than',/
     &           ,' needed at the first space step , DTon0 :',1pg10.3,/)
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
            prescc = dhn*rkb*1.d4
            dh0 = dhn*prescc/(fpressu(te0,dhn,pop))
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
      end if
c
c
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
         dh1 = frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen)
      else
c     error if pop is very different from equilibrium at te1
         dh1 = (dhn*prescc)/(fpressu(te1,dhn,pop))
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
c      call equion(t,de,dh)
c
      call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
      call absdis(t, dh, fi, dtau, dr, rad0, pop)
c
   30 continue
c
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
         dh1 = frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen)
      else
         dh1 = dhn*prescc/(fpressu(te1,dhn,popend))
      end if
      de1 = feldens(dh1,popend)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi1=dmin1(1.d0,fidhcc/dh1)
      rww = 1.d0
      if (jeq.eq.'C') 
     &     te1 = frad(dis1,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen)
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
c
      call zetaeff(fi, dh)
c
c*****************************************************************
c      write(*,*)'*EQUILIBRIUM IONISATION AND TEMPERATURE .'
c*****************************************************************
c
      if ((jeq .eq. 'E').or.(jeq .eq. 'P').and.(jthm.eq.'S')) then
         dton1 = 1.d33
         call copinto(popend, pop)
         call teequi(te1, tef, de1, dh1, dton1, nmod)
      else if ((jeq.eq.'C').or.(jthm.eq.'T')) then
c
         te1 = frad(dis1,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen)
         call copinto(popend, pop)
c
         call equion(te1, de1, dh1)
c
         i = 0
c
         treap = 1.d-12
         difp = 0.d0
c
 80      continue
c
         call copinto(pop, popz)
c
         call localem(te1,de1,dh1)
         call totphot(te1, dh1, fi, rad, dr, dv, wdil, lmod)
         call zetaeff(fi, dh1)
         call equion(te1, de1, dh1)
c
         call difpop(pop, popz, treap, atypes, difp)
c
         call copinto(pop, popz)
c
         i = i+1
c
         if ((dift.ge.0.001).and.(i.le.4)) goto 80
         tef = te1
c         
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
c
         call absdis(t, dh, fi, dtauon, disin, 0.0d0,pop0)
         fthick = dmin1(disin,dis1-distim(1,0))
         dison = dis1-fthick
         protim = distim(2,0)
         do 330 j = m-1, 0, -1
            if (distim(1,j).gt.dison) goto 330
            prf = (dison-distim(1,j))/((distim(1,j+1)-distim(1,j))
     &            + epsilon)
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
         dhf = frad(dis1,dhn,xfac,pfac,afac, bfac,cfac,scalen)
      else
         dhf = dhn*prescc/(fpressu(tef,dhn,pop))
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
c      difg = dmax1(dift/difma,dtl/dtlma,dhl/dhlma)
c     
c     lenient
c     
      difg = dmax1(dtl/dtlma,dhl/dhlma)
c     
c     relaxed cirteria
c     
c     difg = dhl/dhlma
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
            dtau = (0.2d0*rmm)*dtau
         else if (((niter.eq.3).or.(niter.eq.6)).and.((jden.ne. 
     &      'B').or.(difg .eq. (dift/difma)))) then
            agmu = dmin1(1.d0-((1.d0-agmax)/aadn),1.d0-
     &      ((1.d0-agmu)*0.65d0))
            nidh = max0(0,min0(nidh-2,idnint(1.0d0+dint(niter/3.0d0))))
            te1 = (((te1+te0)+te1a)+te0)/4.0
            te1 = dmax1(te0*agmu,dmin1(te0/agmu,te1))
            weiph = 0.75d0
            rmm = 0.4d0
            dtau = (0.2d0*rmm)*dtau
         else
            te1 = tef
            weiph = 0.85d0
            rmm = 0.25d0
            dtau = rmm*dtau
         endif
 34         dtau = dtau*0.5
c 120        continue !NO idea what this is from
c
            mmod = 'DIS'
            lmod = 'SO'
            call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
            call absdis(t, dh, fi, dtau, drzi, rad0, poppre)
            call absdis(t, dh, fi, dtau, drzf, rad0, popend)
c
c
            dr = dexp((weiph*dlog(drzi+.1d0))+((1.d0-weiph)
     &           *dlog(drzf+.1d0)))
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
            goto 30
 31     continue
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
c  Determine dust continuum emission for this region (save for p5)
c
c      if ((grainmode).AND.(IRmode.ne.0)) then
c        call dusttemp(t,dh,de,dr,m)
c      endif

      call zetaeff(fi, dh)
      call localem(t, de, dh)
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
      t      = (te0+(rww*te1))/(1.d0+rww)
      dh     = (dh0+(rww*dh1))/(1.d0+rww)
      de     = (de0+(rww*de1))/(1.d0+rww)
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
      call cool(t, de, dh)
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
c
c    ***OUTPUT QUANTITIES RELATED TO SPACE STEP
c
c
c     write balance for first step or if poll file exists
c
c     
      fhi = pop(1,1)
      fhii = pop(2,1)
c
      pollfile = 'balance'
      inquire(file=pollfile, exist=iexi)
      if (((jbal.eq.'YES').and.(m.eq.1)).or.(iexi)) then
         caller = 'P4'
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
      if (iexi) then
         caller = 'P4'
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
      
      pollfile = 'speclocal'
      inquire(file=pollfile, exist=iexi)
      if (iexi) then
         caller = 'P5'
         pfx = 'local'
         np = 5
      	 sfx = 'ph4'
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
      open(luop, file=filna, status='OLD', access='APPEND')
      open(lupf, file=filnb, status='OLD', access='APPEND')
      open(lups, file=filnc, status='OLD', access='APPEND')
      if (jall.eq.'YES') then
      open(lusl, file=newfil, status='OLD', access='APPEND')
      endif
c
  461 format(i4,0pf9.2,1x,1pg10.3,6(1pg10.3))
 1461 format(t4,10(1pg10.3))
      write(luop, 461) m, t, dlosav, dis1, disav, dr, fi, dh, de
      write(luop, 1461) fhi, fhii, zetae, qhdh,hoiii(m),hoii(m),
     &                  hnii(m),hsii(m),hbeta,caseab(1)
c
      if (jeq .eq. 'F') write(luop, 462) fthick, treach, dtco, 
     &dton, dtoff
  462 format(t71,5(1pg9.2))
c
  473 format(t98,5(1pg9.2))
      if (jeq .eq. 'P') write(luop, 473) recscal, dtoff
c
 465  format('    Case A-B (HI, HeII): ',2(1pg11.3))
      write(*,465) caseab(1), caseab(2)

 460  format(' #',t7,'Teav',t15,'DLOSav',t24,'Rfin',t33,'dr',t42,
     & 'nHav',t51,'FHI',t60,'QHDNf',t69,'[OIII]',/,
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
 3000 format(29(1pg14.7,1x))
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
c
      enddo
c
c     end .tsr files
c
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
 1250 format(//' Step   Te Ave.(K)   ',
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
  513 continue
  517 taux = sigpho(n)*popint(jpoen,ielen)
      if (((jend.eq.'A').or.(jend.eq.'B')).and.(frx.le.fren)) 
     &goto 505
      if ((jend.eq.'C').and.(telast.le.tend)) goto 505
      if ((jend.eq. 'D').and.(taux.ge.tauen)) goto 505
      if ((jend.eq.'E').and.(dis1.ge.diend)) goto 505
      if ((jend.eq.'F').and.(popint(jpoen,ielen).ge.tauen))
     &goto 505
c
      if ((((de1/dh1)-dedhmi)/dedhma).le.exla) goto 505
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
      te0 = te1
c
c     record more ts in thist for quad etxrapolation
c     when m >3
c
      thist(1) = thist(2)
      thist(2) = tep
      thist(3) = te1
c
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
c      write(*,*)'LOOP BACK AND INCREMENT M, DO NEXT STEP.....'
c
c**********************************************************************
c
      m = m+1
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
  967 format(/' NB. :::::::::::::::Integeration through the line', 
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
     &        , atypes)
 2070    format(a7,t8,16(1pg10.3))
 2100 continue
c
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
c
c
      end
      
