cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c     
c     copyright 1994 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c
c       Version 1.0.0r
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c*******PHOTOIONISATIONMODEL
c     DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c     
c     space steps derived from optical depth of predicted temperature
c     and ionistation, extrapolated from previous steps
c     
c     PHOTO5 LITE - EQ only
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine photo6()
c     
      include 'cblocks.inc'
c     
c           Variables
c
      double precision blum,ilum,dhlma,dhn,dht,difma,dtlma
      double precision epotmi,fin
      double precision qhlo,rechy,reclo
      double precision rstsun,vstrlo,wid, starlum,astar
      double precision presreg,tinner,dh10, scale,qhtav
c
      integer*4 j,luin,i
c
      character where*16
      character ilgg*4,inputlum*4
      character banfil*12
c     
c           Functions
c
      double precision densnum,fdilu,fgetdhfromradius
c     
 1010 format(a)
      luin = 10
      limph = atypes
      jcon =  'YES'
      jspot = 'NO'
      dtau0 = 0.025d0
c     
      write(*, 300) 
c     
c     ***INITIAL IONISATION CONDITIONS
c     
 300  format(///
     &     ' Photoionisation model P6 selected:'/
     &     ' Diffuse Field : Outward integration: Radiation Pressure'/
     &     ' *******************************************************',/)
c     
      epotmi = ryde
c     
c     set up ionisation state
c     
      where = 'protoionisation'
      call popcha(where)
      where = 'Photo 6'
c     
c     
c     ***CHOOSING PHOTON SOURCE AND GEOMETRICAL PARAMETERS
c     
      where = 'Photo 6'
c     
 11   call photsou(where)
c
c     Possible to have no photons     
c     
      qhlo = dlog(qht+epsilon)
      rechy = 2.6d-13
      reclo = dlog(rechy)
c     
 333  continue
c     
 330  write(*, 335) 
 350  format(a)
 335  format(/' Spherical or plane II  geometry (s/p) : ',$)
      read(*, 350) jgeo
      jgeo = jgeo(1:1)
c     
      if ((jgeo.eq.'S').or.(jgeo.eq.'s')) jgeo = 'S'
      if ((jgeo.eq.'P').or.(jgeo.eq.'p')) jgeo = 'P'
c     
      if       ((jgeo .ne. 'S')
     &     .and.(jgeo .ne. 'P')) goto 330
c     
c     ***IF SYMMETRY IS SPHERICAL :
c     
      if (jgeo .eq. 'S') then
c     
c     
c     1st total energy in Inu
c     
c     
         blum = 0.d0
         ilum = 0.d0
c     
         do i = 1,infph-1
            wid = ephot(i+1)-ephot(i)
            blum = blum + soupho(i)*wid*evplk
            if (ephot(i).ge.Ryd) ilum = ilum + soupho(i)*wid*evplk
         enddo
c
c     Plane Parallel x pi
c
         blum = pi*blum
         ilum = pi*ilum
c
         rstar = 1.d0
c
         if (blum.gt.0.d0) then
c     
 810     write(*, 820) 
 820     format(/' Define the source size or luminosity '/
     &        ' ::::::::::::::::::::::::::::::::::::::'/
     &        '    R   :  By Radius'/
     &        '    L   :  By Luminosity'/
     &        '    P   :  By Ionizing Photons'//
     &        ' :: ',$)
 825     read(*,350) inputlum
         inputlum = inputlum(1:1)
c     
         if (inputlum.eq.'r') inputlum = 'R'
         if (inputlum.eq.'l') inputlum = 'L'
         if (inputlum.eq.'p') inputlum = 'P'
c     
         if ((inputlum.ne.'R').and.(inputlum.ne.'L').and.
     &       (inputlum.ne.'P')) goto 810
c     
         If (inputlum.eq.'R') then
c     
 826        format(/,/,' Give photoionisation source radius',
     &           ' (Rsun=6.96e10)',/
     &           ' (in solar units (<1.e8) or in cm (>1.e8) : ',$)
            write(*,826)
            read(*, *) rstsun
            if (rstsun.le.0.0d0) goto 810
            if (rstsun.lt.1.d8) then
               rstar = 6.96d10*rstsun
            else
               rstar = rstsun
               rstsun = rstar/6.96d10
            end if
c
            astar = 4.d0*pi*rstar*rstar
            blum = 4.d0*pi*rstar*rstar*blum
            ilum = 4.d0*pi*rstar*rstar*ilum
c     
 829        format(//' *******************************************'/
     &           '  Source Total Luminosity : ',1pg12.5/
     &           '  Source Ryd+  Luminosity : ',1pg12.5/
     &           '  Source Ryd+  Photons    : ',1pg12.5/
     &           ' *******************************************')
            write(*,829) blum,ilum, qht*astar
c
            mdot = 0.d0
c            
c     end def by radius
c     
         endif
c         
         if (inputlum.eq.'L') then
c     
 1827       format(//' Total or Ionising Luminosity (T/I)',$)
            write(*,1827)
 1828       read(*,350) ilgg
            ilgg = ilgg(1:1)
c     
            if (ilgg.eq.'t') ilgg = 'T'
            if (ilgg.eq.'i') ilgg = 'I'
c
            if ((ilgg.ne.'I').and.(ilgg.ne.'T')) ilgg = 'I'
c            
 827        format(//' Give ionising source luminosity '/
     &               ' (log (<100) or ergs/s (>100)) : ',$)
2827        format(//' Give bolometric source luminosity '/
     &               ' (log (<100) or ergs/s (>100)) : ',$)
            if (ilgg.eq.'T') then
            	write(*,2827)
            else
            	write(*,827)
            endif
            
            read(*, *) starlum
            if (starlum.lt.1.d2) starlum = 10.d0**starlum
c     
c     
            if (ilgg.eq.'I')  then
                 rstsun = dsqrt(starlum/(ilum*4.d0*pi))
            else
                 rstsun = dsqrt(starlum/(blum*4.d0*pi))
            endif
c     
c     
 828        format(//' ******************************************'/
     &           '  Source radius : ',1pg12.5/
     &           ' ******************************************')
            write(*,828) rstsun
c     
            astar = 4.d0*pi*rstsun*rstsun
            blum = astar*blum
            ilum = astar*ilum
            write(*,829) blum,ilum,qht*astar
c     
            rstar = rstsun
            rstsun = rstar/6.96d10
c
c     end def by luminosity
c     
         endif
c
         if (inputlum.eq.'P') then
c            
3827        format(//' Give source ionising photon rate '/
     &               ' (log (<100) or photons/s (>100)) : ',$)
            write(*,3827)
            read(*, *) starlum
            if (starlum.lt.1.d2) starlum = 10.d0**starlum
c
c     P/s = qht*4.d0*pi*rstsun*rstsun
c     
             rstsun = dsqrt(starlum/(qht*4.d0*pi))
c     
c     
            write(*,828) rstsun

            astar = 4.d0*pi*rstsun*rstsun
            blum = astar*blum
            ilum = astar*ilum
            write(*,829) blum,ilum,qht*astar
c     
c     
            rstar = rstsun
            rstsun = rstar/6.96d10
c
c    end source radius with photons
c
         endif
c     
      endif
c     
      endif
c     
c     
 420  format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' Setting the physical structure  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
      write(*,420)
c
c     ***THERMAL STRUCTURE
c
      jthm = 'S'     ! self-consistent Temp
      jden = 'C'     ! Isochoric (const density) default
      jeq  = 'E'     ! Equilibrium calculations
c     
c     
c     ***DENSITY BEHAVIOR
c     
 480  write(*, 485) 
 485  format(/' Choose a density structure '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    C  : isoChoric,  (const volume)'/
     &     '    B  : isoBaric, (const pressure)'/
     &     ' :: ',$)
      read(*, 350 ) jden
c     
      if ((jden(1:1).eq.'C').or.(jden(1:1).eq.'c')) jden = 'C'
      if ((jden(1:1).eq.'B').or.(jden(1:1).eq.'b')) jden = 'B'
c     
      if ((jden.ne.'C').and.(jden.ne.'B')) goto 480
c     
c     
c     
c     Estimate strom radius for nominal density
c     if we have any photons that is!
c     
c     if isobaric, get pressure regime
c     
      if (jden.eq.'B') then
 497     format(/' Give pressure regime (p/k, <10 as log) : ',$)
         write(*, 497) 
         read(*,*) presreg
         if (presreg.le.10.d0) presreg = 10.d0**presreg
c     
 499     format(//' Give estimate of mean temperature ',
     &        ' [T ~1e4, <10 as log]'/
     &        ' (for initial size estimates only, actual temperature'/
     &        ' and structure calculated later in detail, '/
     &        ' if in doubt use 1e4 (or 4))'/
     &        ' : ',$)
         write(*, 499) 
         read(*,*) tinner
         if (tinner.le.10.d0) tinner = 10.d0**tinner
c     
         dhn = presreg/tinner
         dh10 = presreg/1e4
c     
 498     format(/' *******************************************'/
     &           '  H density at 1e4 K : ',1pg12.5/
     &           ' *******************************************')
         write(*, 498) dh10
c
c     jden = B isobaric
c     
      endif
c     
 490  continue
      if (jden.eq.'C') then
 491     format(/' Give nominal hydrogen density  : ',$)
         write(*, 491) 
         read(*, *) dhn
      endif
c
      dht = densnum(dhn)
c      
 492  format(/' ********************************************'/
     &        '   Mean Ion Density (N):',1pg12.5/
     &        '   Mean H  Density (dh):',1pg12.5/
     &        ' ********************************************')
c
      write(*, 492) dht,dhn
c     
 493  format(/' Give a filling factor (0<f<=1) : ',$)
      Write(*, 493) 
      read(*, *) fin
      if (((fin.le.0.0d0).or.(fin.gt.1.0d0)).or.(dhn.le.0.0d0)) 
     &     goto 490
c
      if ((jgeo .eq. 'S')) then
c
c       Spherical Geometry
c
      if (qht.gt.0.d0) then
c      
         qtlog = (dlog10(qht)+1.09921d0+2.0d0*dlog10(rstar))
         vstrlo = ((((2.531d0+qhlo)+(2.d0*dlog(rstar)))-(2.d0*dlog(
     &        dht)))-dlog(fin))-reclo
         rmax = dexp((vstrlo-1.43241d0)/3.d0)
c
c     write (*,*) qht,dht
c
         qhdin = qht/dht
         qhdav = (4.d0*qht/dht)*fdilu(rstar,(rmax+rstar)*0.5d0)
         qhtav = (4.d0*qht/(cls*dhn))*fdilu(rstar,(rmax+rstar)*0.5d0)
c     
 312     write(*, 320) rmax, qhdin, qhdav
 320     format(/' ********************************************'/
     &        '  Estimated hydrogen Stromgren radius of the'/
     &        '  order of ................:',1pg10.3,' cm.'/
     &        '  Impact parameters:QHDN in:',1pg10.3,/
     &        '                    QHDN av:',1pg10.3,/
     &        ' ********************************************'//
     &        ' Give initial radius in terms of distance',
     &        ' or Q(N) or U(H) (d/q/u):',$)
         read(*, 1010) ilgg
         ilgg = ilgg(1:1)
c     
         if (ilgg.eq.'D') ilgg = 'd'
         if (ilgg.eq.'Q') ilgg = 'q'
         if (ilgg.eq.'U') ilgg = 'u'         
         if ((ilgg.ne.'q').and.(ilgg.ne.'d').and.
     &       (ilgg.ne.'u')) goto 312
c     
         if (ilgg.eq.'u') then  
 307        write(*,308)
 308        format(/' Give U at inner radius (<=0 as log) : ',$)
            read(*, *) remp
c     
            if (remp.le.0) remp = 10.d0**remp
c     
c     Assume large distance from source and scale by (rmax+rstar)/2
c     
            remp = ((rmax+rstar)/2.d0)*dsqrt(qhtav/remp)
c     
         endif       
         if (ilgg.eq.'q') then  
 314        write(*,321)
 321        format(/' Give Q at inner radius (<= 100 as log) : ',$)
            read(*, *) remp
c     
            if (remp.le.100) remp = 10.d0**remp
c     
c     Assume large distance from source and scale by (rmax+rstar)/2
c     
            remp = ((rmax+rstar)/2.d0)*dsqrt(qhdav/remp)
c     
         endif
c     
         if (ilgg.eq.'d') then      
c     
c     give R_0 directly
c     
 313        write(*,322)
 322        format(/,' Give initial radius ',
     &           ' (in cm (>1) or in fraction of Stromgren radius) : '
     &           ,$)
            read(*, *) remp
c     
            if (remp.lt.0.0d0) goto 312
            if (remp.le.1) remp = remp*rmax
            if (remp.lt.rstar) remp = rstar
c
            if (jden.eq.'I') then
               dhn = fgetdhfromradius(remp)
            endif
c     
         endif
c     
         rmax = dexp(((vstrlo+dlog(1.d0+((remp/rmax)**3)))-
     &        1.43241d0)/3.0d0)
c     
c
c     end with photons (qht>0)
c
         else
c
c     no photons
c
c     
c     give R_0 directly
c     
 1220        write(*,1230)
 1230        format(/' Give initial radius:',$)
             read(*, *) remp
c     
             if (remp.lt.0.0d0) goto 1220
             if (remp.lt.rstar) remp = rstar            
             rmax = remp*100.0d0
c
         endif
c
         wdpl = fdilu(rstar,remp)
         qhdin = 4.d0*qht/dht*wdpl
         wdpl = fdilu(rstar,(rmax+remp)*0.5d0)
         qhdav = 4.d0*qht/dht*wdpl
 323     format(/' ******************************************'/
     &           '  Empty radius : ',1pg12.5/
     &           '  QHDN in      : ',1pg12.5/
     &           '  QHDN <r>     : ',1pg12.5/
     &           ' ******************************************'//)
         write(*,323) remp,qhdin,qhdav
c
c   end Spherical
c
      endif
c     
c     ***IF GEOMETRY PLANE-PARALLEL :
c     
      if (jgeo .eq. 'P') then
         rstar = 0.0d0
         remp = 0.0d0
c
      if (qht.gt.0.d0) then
c      
c     Determine number of ionising photons
c
c     
c     1st total energy in Inu
c     
c     
         blum = 0.d0
         ilum = 0.d0
c     
         do i = 1,infph-1
            wid = ephot(i+1)-ephot(i)
            blum = blum + soupho(i)*wid*evplk
            if (ephot(i).ge.Ryd) ilum = ilum + soupho(i)*wid*evplk
         enddo
c
c     Plane Parallel x pi
c
         blum = pi*blum
         ilum = pi*ilum
c
c
c
 930     write(*,935)
 935     format(/' Give Ionizing Flux at inner edge by: '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    B  : Bolometric Flux (erg/s/cm^2)'/
     &     '    I  : Ionising Flux (erg/s/cm^2)'/
     &     '    U  : Ionisation parameter (Using dh)'/
     &     '    N  : No change (Use current flux)'/
     &     ' :: ',$)
         read(*, 350 ) ilgg
c     
         if ((ilgg(1:1).eq.'B').or.(ilgg(1:1).eq.'b')) ilgg = 'B'
         if ((ilgg(1:1).eq.'I').or.(ilgg(1:1).eq.'i')) ilgg = 'I'
         if ((ilgg(1:1).eq.'U').or.(ilgg(1:1).eq.'u')) ilgg = 'U'
         if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) ilgg = 'N'
         
         if (ilgg.eq.'B') then
            write(*,940) 
 940        format(/,'Give Bolometric flux at inner edge of cloud:')
            read(*,*) scale
            scale = scale/blum
         else if (ilgg.eq.'I') then 
            write(*,945) 
 945        format(/,'Give Ionising flux at inner edge of cloud:')
            read(*,*) scale
            scale = scale/ilum
         else if (ilgg.eq.'U') then 
            write(*,950) 
 950      format(/,'Give Ionisation parameter at inner edge of cloud:')
            read(*,*) scale
            if (scale.le.0) scale = 10**scale
            scale = scale*dhn*cls/qht
         else if (ilgg.eq.'N') then 
            scale=1.d0
         else
            goto 930
         endif
c
         do j = 1,infph
            soupho(j) = soupho(j)*scale
         enddo
c
c     
         blum=blum*scale
         ilum=ilum*scale
         qht=qht*scale
 955     format(//' *********************************************'/
     &     '  Total Flux at inner edge : ',1pg12.5,' (erg/s/cm^2)',/,
     &     '  Ionising Flux at inner edge: ',1pg12.5,' (erg/s/cm^2)',/,
     &     '  Ionisation parameter at inner edge : ',1pg12.5/
     &     ' *********************************************',/)
         write(*,955) blum,ilum,qht/(dhn*cls)
c
c

         qtlog = dlog10(qht)
         vstrlo = ((qhlo-(2.0d0*dlog(dht)))-dlog(fin))-reclo
         rmax = dexp(vstrlo)
c
c        Dilution = 0.5, qht includes geo dil of 0.5
c
         qhdav = ((2.0d0*qht)/dht)*0.5
c     
 400     write(*, 410) qhdav
 410     format(/' ************************************************'/
     &           '  Impact parameter (g=0.5)  QHDN av:',1pg12.5,/
     &           ' ************************************************'//
     &           ' Give the geometrical dilution factor (<=0.5) : ',$)
         read(*, *) wdpl
c     
         if ((wdpl.le.0.0d0).or.(wdpl.gt.0.5d0)) goto 400
c     
         rmax = rmax*wdpl
         qhdin = ((2.0d0*qht)*wdpl)/dht
         qhdav = qhdin
c
      endif
c     
      endif
c     
c     
c     ***EQUILIBRIUM OR FINITE AGE ASSUMPTION
c
c keep input for script compatibility, but allow only eq models
c
c     
 459  format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Ionisation balance calculations'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
      write(*,459)
 460  write(*, 465) 
 465  format(//' Choose type'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '     E  : Equilibrium balance.'/
c     &     '     F  : Finite source lifetime.'/
c     &     '     P  : Post-Equilibrium decay.'/
c     &     '        :'/
c     &     '     C  : CIE for given Te profile.'/
     &     '  :: ',$)
      read(*, 350) jeq
c
      jeq = 'E'
c     
c      if ((jeq(1:1).eq.'E').or.(jeq(1:1).eq.'e')) jeq = 'E'
c      if ((jeq(1:1).eq.'F').or.(jeq(1:1).eq.'f')) jeq = 'F'
c      if ((jeq(1:1).eq.'P').or.(jeq(1:1).eq.'p')) jeq = 'P'
c      if ((jeq(1:1).eq.'c').or.(jeq(1:1).eq.'c')) jeq = 'C'
c     
c      if ((jeq .ne. 'E') ) goto 460
c     &     .and.(jeq .ne. 'F')
c     &     .and.(jeq .ne. 'P'))
c     &     goto 460
c     
 843  write(*, 847) 
 847  format(/' Step value of the photon absorbsion fraction: ',$)
      read(*, *) dtau0
      if ((dtau0.le.0.0d0)) goto 843
c     
c     standard convergence limits
c     
      difma = 0.05d0
      dtlma = 0.1d0
      dhlma = 0.10d0
c     
      if (expertmode.gt.0) then
 2000    format(/' Do you wish to alter the convergence criteria ?',
     &        ' (not recommended) : ',$)
         write(*,2000)
         read(*, 350) ilgg
c     
         if ((ilgg(1:1).eq.'Y').or.(ilgg(1:1).eq.'y')) then 
 849        write(*,850)
 850        format(/' Set the convergence criteria, the maximum',
     &       ' changes allowed'/
     &       ' for three parameters: ion population, te and hydrogen',
     &       'density.'/ 
     &       ' Convergence criteria: difma (0.05),dtlma (0.03)',
     &       ',dhlma(0.015) :',$)
c     
            read(*,*) difma,dtlma,dhlma
         endif
c     
c     end expertmode
c     
      endif
c     
c     ***PRINT SETTING UP
c     
      
         write(*, 623) dtau0         
 623     format(/' Summary : Thermal and Ionic Equilibrium :'/
     &           ' -----------------------------------------'/
     &           '  dTau        :',0pf7.4)
     
c
 620  format(/' Radiation Field :'/
     &   ' -----------------'/
     &   ' Rsou.',t15,' Remp.',t25,' Rmax',t35,' DILUav'/
     &   ,1pg10.3,t15,1pg10.3,t25,1pg10.3,t35,1pg10.3//
     &   ' Hdens',t15,' Ndens',t25,' F.F.',t35,' QHDNin',t45,' QHDNav'/
     &   ,1pg10.3,t15,1pg10.3,t25,0pf9.6,t35,1pg10.3,t45,1pg10.3)
c
      write(*, 620) rstar, remp, rmax, wdpl,
     &              dhn, dht, fin, qhdin , qhdav
c
 503  continue
c     
c     ***TYPE OF EXIT FROM THE PROGRAM
c     
      ielen = 1
      jpoen = 1
      tend  = 10.0d0
      fren  = 0.01d0
      diend = 0.0d0
      tauen = 0.0d0
c     
c     
 507  format(/' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  Boundry conditions  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write(*,507)
c     
 510  format(/' Choose a model ending : '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Radiation bounded, HII <',2pf5.2,'%'/
     &     '    B  :   Ionisation bounded, XII < Y'/
     &     '    C  :   Temperature bounded, Tmin.'/
     &     '    D  :   Optical depth limited, Tau'/
     &     '    E  :   Density bounded, Distance'/
     &     '    F  :   Column Density limited, atom, ion'/
     &     '       :   '/
     &     '    R  :   (Reinitialise)'/
     &     '    G  :   (Reset geometry only)'/
     &     ' :: ',$)
 505  write(*, 510) fren
      read(*, 350) jend
c     
      if ((jend(1:1).eq.'A').or.(jend(1:1).eq.'a')) jend = 'A'
      if ((jend(1:1).eq.'B').or.(jend(1:1).eq.'b')) jend = 'B'
      if ((jend(1:1).eq.'C').or.(jend(1:1).eq.'c')) jend = 'C'
      if ((jend(1:1).eq.'D').or.(jend(1:1).eq.'d')) jend = 'D'
      if ((jend(1:1).eq.'E').or.(jend(1:1).eq.'e')) jend = 'E'
      if ((jend(1:1).eq.'F').or.(jend(1:1).eq.'f')) jend = 'F'
      if ((jend(1:1).eq.'R').or.(jend(1:1).eq.'r')) jend = 'R'
      if ((jend(1:1).eq.'G').or.(jend(1:1).eq.'g')) jend = 'G'
c     
      if ((jend .ne. 'A').and.(jend .ne. 'B').and.(jend .ne. 
     &    'C').and.(jend .ne. 'D').and.(jend .ne. 'E').and.
     &    (jend .ne. 'F').and.(jend.ne. 'R').and.(jend .ne. 'G'))
     &     goto 503
      if (jend .eq. 'R') return 
      if (jend .eq. 'G') goto 333
      if ((jend .eq. 'B').or.(jend .eq. 'D')) then
 525     write(*, 530) 
 530     format(/' Applies to element (Atomic number) : ',$)
         read(*, *) ielen
         if (zmap(ielen).eq.0) goto 525
         ielen = zmap(ielen)
         jpoen = 1
         if (arad(2,ielen).le.0.0d0) jpoen = 2
      end if
c     
      if (jend .eq. 'B') then
 535     write(*, 540) elem(int(ielen))
 540     format(/' Give final ionisation fraction of ',a2,' : ',$)
         read(*, *) fren
         if ((fren.lt.0.0d0).or.(fren.gt.1.0d0)) goto 535
      end if
c     
      if (jend .eq. 'D') then
 545     write(*, 550) elem(int(ielen))
 550     format(/' Give final optical depth at threshold of ',a2,':',$)
         read(*, *) tauen
         if (tauen.le.0.0d0) goto 545
      end if
c     
      if (jend .eq. 'C') then
 555     write(*, 560) 
 560     format(/' Give final temperature (<10 as log): ',$)
         read(*,*) tend
         if (tend.lt.1.0d0) goto 555
      end if
c     
      if (jend .eq. 'E') then
 565     write(*, 570) 
 570     format(/' Give distance or radius at which the density',
     &        ' drops: '/
     &        ' (in cm (>1E6) or as a fraction of the Stromgren'/
     &        ' radius (<1E6)) : ',$)
         read(*,*) diend
         if (diend.lt.1.d6) diend = diend*rmax
         if (diend.le.remp) goto 565
      end if
c     
      if (jend .eq. 'F') then
 585     write(*, 580) 
 580     format(/' Give final column density (<100 as log): ',$)
         read(*,*) tauen
         if (tend.lt.1.0d2) tauen = 10.d0**tauen
 590     write(*, 595) 
 595     format(/' Applies to element (Atomic number) : ',$)
         read(*, *) ielen
         if (zmap(ielen).eq.0) goto 590
         ielen = zmap(ielen)
 596     write(*, 597) 
 597     format(/' Applies to ion stage : ',$)
         read(*, *) jpoen
         if (jpoen.le.0) goto 596
      end if
c     
 2107 format(//
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' Output Requirements  '/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
      write(*,2107)
c     
c     
c     
 2100 format(//' Choose output settings : '/
     &     ' :::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Standard output (photnxxxx,phapnxxxx).'/
     &     '    B  :   Standard + monitor 4 specific ions.'/
     &     '    C  :   Standard + all ions file.'/
     &     '    D  :   Standard + final source + nu-Fnu spectrum.'/
     &     '    F  :   Standard + first balance.'/
     &     '    G  :   Everything'//
     &     ' :: ',$)
 2101 write(*, 2100)
      read(*, 350) ilgg
c     
      if ((ilgg(1:1).eq.'A').or.(ilgg(1:1).eq.'a')) ilgg = 'A'
      if ((ilgg(1:1).eq.'B').or.(ilgg(1:1).eq.'b')) ilgg = 'B'
      if ((ilgg(1:1).eq.'C').or.(ilgg(1:1).eq.'c')) ilgg = 'C'
      if ((ilgg(1:1).eq.'D').or.(ilgg(1:1).eq.'d')) ilgg = 'D'
      if ((ilgg(1:1).eq.'F').or.(ilgg(1:1).eq.'f')) ilgg = 'F'
      if ((ilgg(1:1).eq.'G').or.(ilgg(1:1).eq.'g')) ilgg = 'G'
c     
      jiel = 'NO'
      if ((ilgg.eq.'B').or.(ilgg.eq.'G')) jiel = 'YES'
c     
c     
c     default monitor elements
c     
      iel1 = zmap(1)
      iel2 = zmap(2)
      iel3 = zmap(6)
      iel4 = zmap(8)
c     
      if (jiel.eq.'YES') then
 146     format(//' Give 4 elements to monitor (atomic numbers): ',$)
         if (runmode.ne.'batchrun') then
            write(*,146)
         endif
c     
         read(*,*) iel1,iel2,iel3,iel4
         iel1 = zmap(iel1)
         iel2 = zmap(iel2)
         iel3 = zmap(iel3)
         iel4 = zmap(iel4)
c     
      endif
c     
      jall = 'NO'
      if ((ilgg.eq.'C').or.(ilgg.eq.'G')) jall = 'YES'
c     
      jspec = 'NO'
      if ((ilgg.eq.'D').or.(ilgg.eq.'G')) jspec = 'YES'
c     
      jsou = 'NO'
      jpfx = 'psend'
c     
      if ((ilgg.eq.'D').or.(ilgg.eq.'G')) then
c     
         jsou = 'YES'
c     
c     get final field file prefix
c     
 120     format (a8)
 110     format(//' Give a prefix for final source file : ',$)
         write(*,110)
         read (*,120) jpfx
      endif
c     
      jbal = 'NO'
      jbfx = 'p6bal'
c     
      if ((ilgg.eq.'F').or.(ilgg.eq.'G')) then
c     
         jbal = 'YES'
c     
c     get final field file prefix
c     
 130     format(//' Give a prefix for first ion balance file : ',$)
         write(*,130)
         read (*,120) jbfx
      endif
c     
c     get runname
c     
 100  format (a64)
      
 2120 format(//' Give a name/code for this run: ',$)
      write(*,2120)
      read (*,100) runname
c     
c     current UNIX batch system...
c     
      if (runmode.ne.'batchstart') then
c     
c     
c     remains of old vax batch system (not used)
c     
         banfil = 'INTERACTIVE'
c     
         if (jden.eq.'B') dhn = dh10
c     
         call compph6(dhn, fin, banfil,difma,dtlma,dhlma)
c     
      endif
c     
      return 
c     
      end

