cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c     
c     copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c     Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     PHOTOIONISATION MODEL
c     DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c     
c     RSS 1999
c     
c     PHOTO5:  As for P4 but adds integrated radiation pressure
c     (altered for dust pressure)
c     
c     PHOTO5 Lite - many options removed for clarity and testing
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine compph6(dhn, fin, banfil,difma,dtlma,dhlma)
c     
      include 'cblocks.inc'
c     
c     
c     
      double precision poppre(mxion, mxelem), popend(mxion, mxelem)
      double precision ppre(mxion, mxelem),difma,dtlma,dhlma
      double precision thist(3),steps(3)
c     double precision ctime,rtime,ptime,eqtime,stime
c     
      double precision aadn,agmax,agmu
      double precision a,b,c,x,fx,fp,n0,r0,fv,e,f
      double precision de,dedhma,dedhmi
      double precision dh,dhf,dhl,dhn,dhp,frp
      double precision difg,difmi,difp0,dift
      double precision dlos0,dlos1,dlosav
      double precision dr,drmd,drmu,drp,drxp
      double precision dtau,dtaux,dtco,dtco0,dtco1
      double precision dtl,dtoff,dton,dton0,dton1
      double precision dv,dva,dvoluni,exla
      double precision fhi,fhii,fi,fi0,fi1,fidhcc,fin,frdw,frx,fthick
      double precision prescc,rad,recscal,rww
      double precision t,taux,te1a,tef,telast,tempre
      double precision tep,texm,trea,runit 
      double precision treach,treach0,treach1,tspr
      double precision wdil0,wdil1,wei
      double precision ponk, radpint, pres0, pres1

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
      real tarray(2)
      real dtarr,dtime
c     
      character pollfile*12,tab*4
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
      double precision fcrit,fdilu,feldens,fradpress, fpressu
c     
c     internal functions
c     
      double precision dif,frad
c     
      dif(a,b) = dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
      
      frad(x,n0,fx,fp,a,b,c,r0,fv,e,f)  = 
     &     n0*(fx*dexp((x-a)/r0)+fp*((x/a)**b)+fv*(((r0-x)/e)**(-f)))+c
c     
      jspot = 'NO'
      jcon = 'YES'
      ispo = 'SII'
      IRfile=''
      tab = char(9)
c     
      limph = atypes
c     
      mma = nstmax-1
c     scalen = 1.0d16
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
      IR1 = 91

c     
      ieln = 4
      do i = 1,ieln
         luions(i) = 22+i
      enddo
c     
c     
c     ***DERIVE FILENAME : PHN**.PH5
c     
      call p6fheader(newfil,banfil,fnam,
     &     filna,filnb,filnc,filnd,filn,
     &     luop,lupf,lusl,lusp,lups,luions,dhn,fin)

      chargefile = ' '
c     
c     Set up for PAH existence
c     
      if (grainmode.AND.(pahCfrac.ne.0.d0)) then
         k=1
         do while (ephot(k).le.6)
            k=k+1
         enddo
         fuvmin=k-1
         do while (ephot(k).le.1.359907296E+1)
            k=k+1
         enddo
         fuvmax=k-1
         do while (ephot(k).le.262.012)
            k=k+1
         enddo
         pahabsmax=k-1
c     
c     Calculate photon flux in FUV range using Weingartner & 
c     Draine (2001) ISRF 
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
            else                !11.2<eng<13.6
               Ffuv=Ffuv+4.126E-13*eng**(-6.4172)/eV*pahnabs(k)*dnu
            endif
         enddo
         pahQ0=Ffuv/pahsum
         write(*,520) pahQ0
 520     format(/'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/
     &        'ISRF photon flux :',1pg12.5,/
     &        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      endif
c     
c*********************************************************************
c     
c     write(*,*)'***COMPUTATION STARTS HERE'
c     
c*********************************************************************
c     
      fi     = fin
      fi0    = fin
      fi1    = fin
      fidhcc = fin*dhn
      
      dv = 0.0d0
      
      runit = (rmax-remp)/300.d0
      dtau = dtau0*0.5d0
      
      if (jgeo .eq. 'S') then
         vunilog = 3.0d0*dlog10(runit)
      else
         vunilog = dlog10(runit)
      end if

      nmod  = 'EQUI'

      if (jden .ne. 'B') then
         agmax = 0.05d0
         aadn  = 2.5d0
      else
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
c     write(*,*)'***ITERATE, M = STEP NUMBER'
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
c     write(*,*)'IF M=1 , CHECK CONVERGENCE FOR THE DENSITY'
c     write(*,*)'AND IONIC POP. AT THE INNER BOUNDARY'
c     
c****************************************************************
      if (m .eq. 1) then
c     
         te0 = 1.d4
c     
         dr = 0.d0
         drp = dr
         if (jden .ne. 'F') then
            dh0 = dhn
         else
            dh0 = frad(remp,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,
     &           efac,ffac)
         endif
c     
         de0 = feldens(dh0,pop0)
         dedhmi = 0.d0
         dedhma = 0.d0
         do  j = 1, atypes
            dedhma = dedhma+zion(j)
            if (arad(2,j).le.0.d0) dedhmi = dedhmi+zion(j)
         enddo
         
         call copinto(pop0, pop)
         
 12      continue
c     
c     
         if (jgeo .eq. 'S') then
            rad0 = remp
            dis0 = rad0
            wdil0 = fdilu(rstar,rad0)
         endif
         if (jgeo .eq. 'P') then
            rad0 = 0.d0
            dis0 = remp
            wdil0 = wdpl
         end if

c     
         if ((fin.ne.1.d0).and.(jden.eq.'B')) fi0=dmin1(1.d0,fidhcc/dh0)
c     
c     First field is the inner edge, source term is only one non-zero
c     and dr = 0.d0 to just get the inner zone boundary field.
c     mode could also be 'SO' or source only, but with dr = 0.d0 
c     DW == SO
c     
         lmod = 'DW'
         call localem(te0,de0,dh0)
         call totphot(te0, dh0, fi0, rad0, 0.d0, dv, wdil0, lmod)

         dton0 = 1.d33
         call teequi(te0, te0, de0, dh0, dton0, nmod)

         dhp   = dh0
         tep   = te0
         dlos0 = dlos
c     
         if (jden .eq. 'C') then
            dh0 = dhn
         else
            frp    = fradpress(dr,dhn)
            prescc = dhn*rkb*1.d4 + radPint + frp
            ponk   = (prescc/rkb)
            pres0  = fpressu(te0,dhn,pop)
            dh0    = dhn *prescc/(pres0)
c     write(*,'(" P/k Pr/k dPr/k",3(1pg12.5))') 
c     &               ponk, radPint/rkb, frp/rkb
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
c     write(*,*)'END OF INNER BOUNDRY FOR M = 1'
c*****************************************************
c     
      else                      !m ne 1
c     
c     calculate qhdn into zone
c     
         q1 = 0.0d0
         q2 = 0.0d0
         q3 = 0.0d0
         q4 = 0.0d0
c     
         call intvec(tphot, q1, q2, q3, q4)
c     
         qhdin=q4/(dh1*zen)
c     
      end if
c     
c     
c     
c*****************************************************************
c     
c     write(*,*)'***ESTIMATE OF VALUES AT THE OUTER BOUNDARY'
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
      if (te1.lt.0.d0) te1 = 290.d0
      
      if (jden .eq. 'C') then
         dh1 = dhn
      else if (jden .eq. 'F') then
         dh1 = frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac
     &        ,efac,ffac)
      else
         frp    = fradpress(dr,dh)
         prescc = dhn*rkb*1.d4 + radPint + frp
         ponk   = (prescc/rkb)
         pres1  = fpressu(te1,dhn,pop)
         dh1    = dhn *prescc/(pres1)
c     write(*,'(" P/k Pr/k dPr/k",3(1pg12.5))') 
c     &               ponk, radPint/rkb, frp/rkb
         
      end if

      rww = 1.d0
      
      t  = (te0+(rww*te1))/(1.d0+rww)
      dh = (dh0+(rww*dh1))/(1.d0+rww)
      wei = 1.d0/(1.d0+rww)
      
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)

c     
c***********************************************
c     write(*,*)'***DETERMINE SPACE STEP : DR'
c***********************************************
c     
c     Get field at start of zone, dr = 0.d0, and then guess a step
c     size with the current opacities in pop
c     
c     
      lmod = 'DW'
      call totphot(t, dh, fi, rad0, 0.d0, dv, wdil0, lmod)
      mmod = 'DIS'
      call absdis(t, dh, fi, dtau, dr, rad0, pop)
c     
 30   continue
c     
      if (jgeo .eq. 'S') then
         rad1 = rad0+dr
         dis1 = rad1
         wdil1 = fdilu(rstar,rad1)
         dvoluni = ftpi*(((rad1/runit) ** 3)-((rad0/runit) ** 3))
      endif
      if (jgeo.eq. 'P') then
         dis1 = dis0+dr
         rad1 = 0.d0
         wdil1 = wdpl
         dvoluni = dr/runit
      endif
c     
c**********************************************
c     write(*,*)'*** dr done',dr
c**********************************************
c     
c     
c******************************************************************
c     
c     write(*,*)'CALCULATION OF IONIC POPULATION AT THE OUTER'
c     of current step
c     
c******************************************************************
      if (jden .eq. 'C') then
         dh1 = dhn
      else if (jden .eq. 'F') then
         dh1 = frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac
     &        ,efac,ffac)
      else
         frp    = fradpress(dr,dh)
         prescc = dhn*rkb*1.d4 + radPint + frp
         ponk   = (prescc/rkb)
         pres1  = fpressu(te1,dhn,popend)
         dh1    = dhn *prescc/(pres1)
      end if
      
      de1 = feldens(dh1,popend)
c     
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi1=dmin1(1.d0,fidhcc/dh1)

      if (jgeo .eq. 'S') then
         rww = (1.d0+((0.75d0*dr)/(rad0+dr))) ** 2
      else
         rww = 1.d0
      end if
      
      t  = (te0+(rww*te1))/(1.d0+rww)
      dh = (dh0+(rww*dh1))/(1.d0+rww)
      de = (de0+(rww*de1))/(1.d0+rww)
c     
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
c     
      wei = 1.d0/(1.d0+rww)
      call averinto(wei, poppre, popend, pop)
      

      lmod = 'DW'
      
      call localem(t, de, dh)
      call totphot(t, dh, fi, rad0, dr, dv, wdil0, lmod)
      call zetaeff(fi, dh)
c     
c*****************************************************************
c     write(*,*)'*EQUILIBRIUM IONISATION AND TEMPERATURE .'
c*****************************************************************
c     

      dton1 = 1.d33
      call copinto(popend, pop)
      call teequi(te1, tef, de1, dh1, dton1, nmod)
c     
c*****************************************************************
c     write(*,*)'END EQUILIBRIUM IONISATION AND TEMPERATURE .'
c*****************************************************************

c     
      if (jden .eq. 'C') then
         dhf = dhn
      else if (jden .eq. 'F') then
         dhf = frad(dis1,dhn,xfac,pfac,afac, bfac,cfac,scalen,vfac
     &        ,efac,ffac)
      else
         frp    = fradpress(dr,dh)
         prescc = dhn*rkb*1.d4 + radPint + frp
         ponk   = (prescc/rkb)
         pres1  = fpressu(tef,dhn,pop)
         dhf    = dhn *prescc/(pres1)
         write(*,'(" P/k Pr/k dPr/k",3(1pg12.5))') 
     &        ponk, radPint/rkb, frp/rkb

      end if
      
      def = feldens(dhf,pop)

c     
c     get a final field consistent with equilibrium
c     
      if ((grainmode).AND.(IRmode.ne.0)) then
         call dusttemp(tef,dhf,def, fi, dr, m)
      endif
c     
c************************************************************
c     write(*,*)'***COMPARES END VALUES WITH PREVIOUS STEP'
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
c     difg = 0.d0
c     
c     
      niter = niter+1
      if (difg.le.difmi) then
         texm = te1
         dtaux = dtau
         difmi = difg
      end if
c     
      if ( niter .eq. 1) then
         write(*, 471) 
 471     format(t7,'    Tend     dhl     dtl    dift',
     &        '    difg   difp0     dtau       dr')
      endif
      write(*, 472) niter, tef, dhl, dtl, dift, difg, difp0, dtau,dr
 472  format(t3,i4,0pf8.0,5(0pf8.3),3(1pg9.2))
      
c     c
c     c****************************************************************
c     c      write(*,*)'CHECK CONVERGENCE OF END VALUES WITH PREVIOUS STEP'
c     c****************************************************************
c     c
c     
c     simple, iterate 1 times usually
c     
      if ((niter.lt.2).or.(difg.gt.1.d0)) goto 30
c     
      te1 = tef 
      dlos1 = dlos
      telast = t
      frx = 1.d0-pop(jpoen,ielen)
c     
c*****************************************************************
c     write(*,*)'END CHECK CONVERGENCE STEP'
c*****************************************************************
c     
c     
c***************************************************************
c     write(*,*)'AVERAGE QUANTITIES FOR THE SPACE STEP CONSIDERED'
c***************************************************************
c     

      if (jgeo .eq. 'S') then
         rww = (1.d0+((0.75d0*dr)/(rad0+dr))) ** 2
      else
         rww = 1.d0
      end if

      t = (te0+(rww*tef))/(1.d0+rww)
      dh = (dh0+(rww*dhf))/(1.d0+rww)
      de = (de0+(rww*def))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi = dmin1(1.d0,fidhcc/dh)
      disav = (dis0+(rww*dis1))/(1.d0+rww)
      dlosav = (dlos0+(rww*dlos1))/(1.d0+rww)
      dton   = (dton0+(rww*dton1))/(1.d0+rww)
      treach = (treach0+(rww*treach1))/(1.d0+rww)
      dtco   = (dtco0+(rww*dtco1))/(1.d0+rww)
      
      wei = 1.d0/(1.d0+rww)
      call averinto(wei, poppre, popend, pop)
c     
      if (jgeo .eq. 'S') then
         rad = disav
         wdil = fdilu(rstar,rad)
      endif
      if (jgeo .eq. 'P') then
         wdil = wdpl
         rad = 0.d0
      end if


c****************************************************************
c     write(*,*)'INTEGRATES COLUMN DENSITY AND END DIFFUSE FIELD'
c****************************************************************
      
      if (jgeo .eq. 'S') then
         jmod = 'OUTW'
         frdw = 1.d0
      endif
      if (jgeo .eq. 'P') then
c     jmod = 'LOUP'
c     frdw = 0.0d0
         jmod = 'OUTW'
         frdw = 1.0d0
      end if

c     
c     compute diffuse field for zone implied by average state
c     
      call localem(t, de, dh)
      call zetaeff(fi, dh)
c     
c     sum diffuse field over dr, accumulating fraction that reaches
c     rad0+dr.
c     
      
      jmod = 'OUTW'
      call newdif(t, t, dh, fi, rad0, dr, dv, 0.d0, dv, frdw, jmod)

c     
c     sum ion columns
c     
      imod = 'COLD'
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
      
c     
c*****************************************************************
      lmod = 'DW'
      call totphot(t, dh, fi, rad0+dr, 0.d0, dv, wdil1, lmod)
c
      frp     = fradpress(dr,dh)
      radPint = radPint + frp
c*********************************************************
c     
c     write(*,*)'COMPUTES SPECTRUM OF THE REGION CONSIDERED'
c     
c*********************************************************
c     
      imod = 'REST'
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
c     
c     record line ratios
c     
c     hoii(m) = (fbri(1,ox2)+fbri(2,ox1))/(hbri(2)+epsilon)
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
c     Test for PAH existence (if grainmode true)
c     
      if (grainmode.AND.(pahfrac.gt.0.d0)) then
c     
c     Calculate habing photodissociation parameter
c     
         Ffuv=0.d0
         do i=fuvmin,fuvmax-1
            dnu=(ephot(i+1)-ephot(i))*evplk
            eng=0.5*(ephot(i+1)+ephot(i))*eV
c     photons s-1 cm-2
            Ffuv=Ffuv+4*pi*tphot(i)*dnu/eng
         enddo
         Hab=Ffuv/(cls*dh)
c     
c     Calcualte PAH ionization/Local ISRF ratio (Q factor)
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
 210     format('Habing parameter =',1pg12.5,', PAH ion/ISRF =',1pg12.5)
c     
c     Determine if PAHs exist
c     If so, deplete gas of Carbon in PAHs
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
     &                 zion0(zmap(6))*(1.d0-dion(zmap(6)))*pahCfrac
               endif
            endif
         endif
      endif
c     
c     
c     ***OUTPUT QUANTITIES RELATED TO SPACE STEP
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
         caller = 'P6'
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
         caller = 'P6'
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
         caller = 'P6'
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
         caller = 'P6'
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
         caller = 'P6'
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
         caller = 'P6'
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
         caller = 'P6'
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
         caller = 'P6'
         pfx = 'local'
         np = 5
      	 sfx = 'ph5'
         call newfile(pfx,np,sfx,3,fspl)
         
         open(lusp, file=fspl, status='NEW')

         linemod = 'LAMB'
         spmod   = 'REL'
         call speclocal(lusp
     &        ,tloss,eloss,egain,dlos
     &        ,t,dh,de, pop(1,1)
     &        ,disav,dr
     &        ,linemod,spmod)
         
         close(lusp) 
         
      endif
c     
c     Infrared output (every 1 step)
c     
      pollfile='IRlocal'
      inquire(file=pollfile, exist=iexi)
      if ((IRmode.ne.0).AND.iexi) then
         if (IRfile.EQ.'') then
            caller = 'P6'
            np = 6
            pfx = 'IRflux'
            sfx = 'sou'
            call newfile(pfx,np,sfx,3,IRfile)
            open(IR1, file=IRfile, status='NEW')
            write(IR1,2500)
 2500       format('%  Infrared flux per region ',/,
     &           '%  MAPPINGS III v1.2.0md',
     &           '%  given as:',/,
     &           '%    energy edge (eV) ,,continuumflux, ', 
     &           'IRflux Fnu(erg s-1 cm-2 Hz-1Sr-1)',/,
     &           '%')
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
c     do i=1,mxinfph
c     storephot(i)=tphot(i)
c     enddo
c     call totphot(t, dh, fi, rad, dr, dv, wdil, 'LOCL')
c     call wpsoufile(caller,IRneb,
c     &        t,de,dh,
c     &        0.d0,0.d0,qhdh,
c     &        disav,dr,
c     &        0.d0,0.d0,
c     &        1.d0,tphot)
c     do i=1,mxinfph
c     tphot(i)=storephot(i)
c     enddo        
      endif
c     
c     graincharge output
c     
      if (grainmode) then
         pollfile = 'graincharge'
         inquire(file=pollfile, exist=iexi)
         if (iexi) then
            if (chargefile.eq.' ') then
               caller = 'P6'
               np = 5
               pfx = 'grpot'
               sfx = 'ph5'
               call newfile(pfx,np,sfx,3,chargefile)
               open(IR1, file=chargefile, status='NEW')
               write(IR1,2510)
 2510          format('%  grain charge in each region ',/,
     &              '%  MAPPINGS III v1.2.0md',
     &              '%  given as:',/,
     &              '%  m,Av. distance,Temp,dh,de,',/,
     &              '%  grain charge(allsizes) (graphite),',/,
     &              '%  grain charge(allsizes) silicate',/,
     &              '%')
               write(IR1,2511) (grainrad(i),i=mindust(1),maxdust(1))
 2511          format('gra radii',10(1pg11.4))
               write(IR1,*)
               write(IR1,2512) (grainrad(i),i=mindust(2),maxdust(2))
 2512          format('sil radii',10(1pg11.4))
               close(IR1)
            endif
            open(IR1, file=chargefile,status='OLD',access='APPEND')
            write(IR1,2515) m,disav,t,dh,de
            write(IR1,2516)(grainpot(i,1),i=mindust(1),maxdust(1))
            write(IR1,2517)(grainpot(i,2),i=mindust(2),maxdust(2))
            close(IR1)
 2515       format(/,i3,4(1pg14.5))
 2516       format('gra ',30(1pg11.4))
 2517       format('sil ',30(1pg11.4))
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
 461  format(i4,1pg11.4,x,1pg10.3,5(1pg11.4),2(1pg10.3),2(1pg11.4)
     &     ,2(1pg10.3))
 1461 format(t4,10(1pg10.3))
      write(luop, 461) m, t, dlosav, disav, dis1, dr, dh, de,fhi, fhii,
     &     qhdin, qhdh,caseab(1), caseab(2)
c     
c     write(luop, 461) m, t, dlosav, dis1, disav, dr, fi, dh, de
c     write(luop, 1461) fhi, fhii, zetae, qhdh,hoiii(m),hoii(m),
c     &                  hnii(m),hsii(m),hbeta,caseab(1)
c     
      if (jeq .eq. 'F') write(luop, 462) fthick, treach, dtco, 
     &     dton, dtoff
 462  format(t71,5(1pg9.2))
c     
 473  format(t98,5(1pg9.2))
      if (jeq .eq. 'P') write(luop, 473) recscal, dtoff
c     
c     465  format('    Case A-B (HI, HeII): ',2(1pg11.3))
c     write(*,465) caseab(1), caseab(2)

 460  format(/,' #',t7,'Teav',t15,'DLOSav',t24,'Rfin',t33,'dr',t42,
     &     'nHav',t51,'FHI',t60,'QHDNf',t69,'[OIII]',/
     &     i4,0pf8.0,1pg9.1,6(1pg9.2))
c     
      write(*, 460) m,t,dlosav,dis1,dr,dh,fhi,qhdh,hoiii(m)
c     
c     
c     
      if (jiel.eq.'YES') then
         do i = 1 ,ieln 
            open(luions(i),file=filn(i),status='OLD',access='APPEND')
c     
 3000       format(29(1pg14.7,x))
            if (i.eq.1) then
               write(luions(i),3000) disav,dr,t,
     &              (pop(j,iel1), j = 1,maxion(iel1))
            endif
            if (i.eq.2) then
               write(luions(i),3000) disav,dr,t,
     &              (pop(j,iel2), j = 1,maxion(iel2))
            endif
            if (i.eq.3) then
               write(luions(i),3000) disav,dr,t,
     &              (pop(j,iel3), j = 1,maxion(iel3))
            endif
            if (i.eq.4) then
               write(luions(i),3000) disav,dr,t,
     &              (pop(j,iel4), j = 1,maxion(iel4))
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
c     
      dtarr = dtime(tarray)
c     
 100  format('      Step time (map3,sys): ',2(f8.1,x),/) 
      write(*,100) tarray(1),tarray(2)
c     
c     
      if (jall.eq.'YES') then
         write(lusl, 1250)
         write(lusl, 678) m, t, de, dh, fi, 
     &        dvoluni, dr, disav, (disav-remp)
 1250    format(//' Step   Te Ave.(K)   ',
     &        '  ne(cm^-3)   ',
     &        '  nH(cm^-3)   ',
     &        '  fill. Fact. ',
     &        '  dVol Unit.  ', 
     &        '    dR (cm)   ',
     &        ' Dist.Ave.(cm)',
     &        ' D - Remp.(cm)')
 678     format(i3,8(1pg14.7)/)
      endif
c     
      wmod = 'PROP'
      call wmodel(lupf,t,de,dh,disav,wmod)
      wmod = 'LOSS'
      call wmodel(lups,t,de,dh,disav,wmod)
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
c     CHECKPOINT BATCHMAP
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
     &        (ionpho(n).eq.jpoen)) goto 517
c     
c     Note: loop out of bound for Column measurements of HII
c     
 513  continue
 517  taux = sigpho(n)*popint(jpoen,ielen)
      if (((jend.eq.'A').or.(jend.eq.'B')).and.(frx.le.fren)) 
     &     goto 505
      if ((jend.eq.'C').and.(telast.le.tend)) goto 505
      if ((jend.eq.'D').and.(taux.ge.tauen)) goto 505
      if ((jend.eq.'E').and.(dis1.ge.diend)) goto 505
      if ((jend.eq.'F').and.(popint(jpoen,ielen).ge.tauen))
     &     goto 505
c     
c     Check Visual extinction condition or end when HII<1%
c     

      if (jend.eq.'H') then
         A_v = dustAv*dustint
         if ((A_v.gt.A_vend).OR.(de/dh.lt.1e-5)) goto 505
      endif
c     
      if (((((de1/dh1)-dedhmi)/dedhma).le.exla).and.(.NOT.grainmode))
     &     goto 505
      pollfile = 'terminate'
      inquire(file=pollfile, exist=iexi)
      if (iexi) goto 505
c     
c     
c     ***RESET INNER BOUNDARY QUANTITIES FOR NEXT SPACE STEP
c     
      call copinto(popend, poppre)
c     
      tep = te0
      te0 = tef
c     
c     record more ts in thist for quad etxrapolation
c     when m >3
c     
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
c     write(*,*)'LOOP BACK AND INCREMENT M, DO NEXT STEP.....'
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
 505  continue
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
 530  format(//' Model ended',t20,'Tfinal',t28,'DISfin',t38,'FHXF',t48,
     &     'TAUXF',t58,'Ending')
      write(luop, 510) t, dis1, frx, taux, jend
 510  format(12h -----------,t18,0pf8.0,3(1pg10.3),4x,a4)

      write(lupf, 2047) tempre, tspr
 2047 format(//t4,'Preionisation conditions for step#0 ',
     &     'for all elements   (TEpr :',0pf8.0,' Time :',1pg9.2,' )  :'/)
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
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)

c     
c     Down stream photon field
c     
         caller = 'P4'
         pfx = jpfx
         np = 5
         wmod = 'REAL'
         dva = 0
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
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
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
c     
c     Down stream nebula only photon field
c     
         caller = 'P4'
         pfx = jpfx
         np = 5
         wmod = 'REAL'
         dva = 0
         call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,qhdh,
     &        disav,dr,
     &        0.d0,0.d0,
     &        1.d0,tphot)
c     
c     
      endif
c     
c     ***COMPUTE AVERAGE SPECTRUM AND WRITES IT IN FILE PHN
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
      
