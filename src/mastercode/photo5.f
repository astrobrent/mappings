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
c     DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION, RADIATION 
c     PRESSURE
c     
c     space steps derived from optical depth of predicted temperature
c     and ionistation, extrapolated from previous steps
c     
c     
c     CHOICE OF EQUILIBRIUM OR FINITE AGE CONDITIONS
c     
c     CALL SUBR.: PHOTSOU,COMPPH5
c     
c     NB. COMPUTATIONS PERFORMED IN SUBROUTINE COMPPH5
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine photo5()
c     
      include 'cblocks.inc'
c     
c           Variables
c
      double precision blum,ilum,dhlma,dhn,dht,difma,dthre,dti,dtlma
      double precision epotmi,fin
      double precision qhlo,rcin,rechy,reclo
      double precision rstsun,vstrlo,wid, starlum,astar
      double precision presreg,tinner,dh10, scale,qhtav
c
      integer*4 j,luin,i
c
      character carac*36,where*16
      character caract*4,ilgg*4,fnam*80,inputlum*4
      character banfil*12, str*80
c
      logical iexi
c     
c           Functions
c
      double precision densnum,fdilu,fgetdhfromradius
c     
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
 300  format(///,
     &     ' Photoionisation model P5 selected:',/,
     &     ' Diffuse Field : Outward integration: Radiation Pressure',/,
     &     ' *******************************************************')
c     
      epotmi = ryde
c     
c     set up ionisation state
c     
      where = 'protoionisation'
      call popcha(where)
      where = 'Photo 4'
c     
c     artificial support for low ionisation species
c     
      if (expertmode.gt.0) then
c     
 1010    format(a)
         carac(1:33) = '   '
         do 1120 i = 3, atypes
            j = 1+index(carac(1:33),'   ')
            caract = '   '
            if (epot(1,i).lt.epotmi) write(caract, 1010) elem(i)
            carac(j:j+2) = caract // '   '
 1120    continue
         i = j+2
 1110    write(*, 1112) carac(1:i)
         ilgg = '    '
c     
 1112    format(//,' Allow the elements : ',a,/,
     &        ' to recombine freely? (y/n) : ',$)
c     
         read(*, 1010) ilgg
         ilgg = ilgg(1:1)
c     
         if (ilgg .eq. 'Y') goto 1116
         if (ilgg .eq. 'y') goto 1116
         if ((ilgg.ne.'N').and.(ilgg.ne.'n')) goto 1110
c     
         do 1118 i = 3, atypes
            arad(2,i) = dabs(arad(2,i))
            if (epot(1,i).lt.epotmi) arad(2,i) =-arad(2,i)
 1118    continue
 1116    continue
c     
      endif
c     
c     
c     ***CHOOSING PHOTON SOURCE AND GEOMETRICAL PARAMETERS
c     
      where = 'Photo 4'
c     
 11   call photsou(where)
c
c     Possible to have no photons     
c     
      qhlo = dlog(qht+epsilon)
      rechy = 2.6d-13
      reclo = dlog(rechy)
c     
c     
c     
 330  write(*, 335) 
 350  format(a)
 335  format(/,' Spherical or plane II  geometry (s/p/f) : ',$)
      read(*, 350) jgeo
      jgeo = jgeo(1:1)
c     
      if ((jgeo.eq.'S').or.(jgeo.eq.'s')) jgeo = 'S'
      if ((jgeo.eq.'P').or.(jgeo.eq.'p')) jgeo = 'P'
      if ((jgeo.eq.'F').or.(jgeo.eq.'f')) jgeo = 'F'
c     
      if       ((jgeo .ne. 'S')
     &     .and.(jgeo .ne. 'P')
     &     .and.(jgeo .ne. 'F')) goto 330
c     
c     ***IF SYMMETRY IS SPHERICAL :
c     
      if ((jgeo .eq. 'S').or.(jgeo.eq.'F')) then
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
 820     format(/,' Define the source size or luminosity ',/,
     &        ' ::::::::::::::::::::::::::::::::::::::',/,
     &        '    R   :  By Radius',/,
     &        '    L   :  By Luminosity',/,
     &        '    P   :  By Ionizing Photons',//,
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
     &           ' (Rsun=6.96e10)',/,
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
 829        format(//' *******************************************',/,
     &           '  Source Total Luminosity : ',1pg12.5,/,
     &           '  Source Ryd+  Luminosity : ',1pg12.5,/,
     &           '  Source Ryd+  Photons    : ',1pg12.5,/,
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
 1827       format(//,' Total or Ionising Luminosity (T/I)',$)
            write(*,1827)
 1828       read(*,350) ilgg
            ilgg = ilgg(1:1)
c     
            if (ilgg.eq.'t') ilgg = 'T'
            if (ilgg.eq.'i') ilgg = 'I'
c
            if ((ilgg.ne.'I').and.(ilgg.ne.'T')) ilgg = 'I'
c            
 827        format(//,' Give ionising source luminosity ',/,
     &               ' (log (<100) or ergs/s (>100)) : ',$)
2827        format(//,' Give bolometric source luminosity ',/,
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
 828        format(//,' ******************************************',/,
     &           '  Source radius : ',1pg12.5,/,
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
3827        format(//,' Give source ionising photon rate ',/,
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
 420  format(//,' :::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     ' Setting the physical structure  '/,
     &     ' :::::::::::::::::::::::::::::::::::::::::::::::::'/,/)
      write(*,420)
c
c     ***THERMAL STRUCTURE
c
      jthm = 'S'
      if (expertmode.gt.0) then
c
 1480 write(*, 1485) 
 1485 format(/,' Choose a thermal structure '/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    S  : Self-consistent,  (recommended)',/,
     &     '    T  : isoThermal, (non-physical, not recommended)',/,
     &     ' :: ',$)
      read(*, 350 ) jthm
c     
      if ((jthm(1:1).eq.'S').or.(jthm(1:1).eq.'s')) jthm = 'S'
      if ((jthm(1:1).eq.'T').or.(jthm(1:1).eq.'t')) jthm = 'T'
      if ((jthm(1:1).eq.'F').or.(jthm(1:1).eq.'f')) jthm = 'F'
      if ((jthm(1:1).eq.'I').or.(jthm(1:1).eq.'i')) jthm = 'I'
c     
      if ((jthm.ne.'S').and.(jthm.ne.'T').and.(jthm.ne.'F').and.
     &    (jthm.ne.'I')) goto 1480
c
      if (jthm.eq.'T') then
c
 1486    format(/,' Give the fixed electron temperature (<=10 as log):'
     &         ,$)
         write(*,1486)
         read(*,*) fixtemp
c
         if (fixtemp.le.10.d0) fixtemp = 10.d0**fixtemp
c
      endif
      endif
c     
c     
c     ***DENSITY BEHAVIOR
c     
 480  write(*, 485) 
 485  format(/,' Choose a density structure '/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    C  : isoChoric,  (const volume)',/,
     &     '    B  : isoBaric, (const pressure)',/,
     &     '    F  : Functional form, F(r)',/,
     &     '    I  : Input file: r dh(r)',/,
     &     ' :: ',$)
      read(*, 350 ) jden
c     
      if ((jden(1:1).eq.'C').or.(jden(1:1).eq.'c')) jden = 'C'
      if ((jden(1:1).eq.'B').or.(jden(1:1).eq.'b')) jden = 'B'
      if ((jden(1:1).eq.'F').or.(jden(1:1).eq.'f')) jden = 'F'
      if ((jden(1:1).eq.'I').or.(jden(1:1).eq.'i')) jden = 'I'
c     
      if ((jden.ne.'C').and.(jden.ne.'B').and.(jden.ne.'F').and.
     &    (jden.ne.'I')) goto 480
c     
      if (jden.eq.'I') then
 425     fnam = 'default'
         if (runmode.ne.'batchrun') then
            write (*,430)
 430        format(/,' Enter file name : ',$)
         endif
         read(*,350) structfile
         inquire(file=structfile, exist=iexi)
c     
         if (iexi) then
c     
c     found the file...
c
            open(luin, file=structfile, status='OLD') 
 435        read(luin,350) str
            if (str(1:1).eq.'%') goto 435
            write(*,350) str
            read(luin,*) nsteps
            write(*,*) nsteps, ' Data points to read.'
            do i = 1,nsteps
               read(luin,*) dist(i),dhy(i)
c               write(*,*) dist(i),dhy(i)
            enddo
            close(luin)
c
c     dummy value for dhn = will set it to dh(r) for rinner later
c
            dhn = dhy(1)
c
            write(*,*) ' Done.'
c
         else
c     
            write(*,*)' File Not Found :',structfile
            goto 480
c     
         endif
c
c end  input file structure
c
      endif
c     
      if (jden .eq. 'F') then
 732     write(*, 735) 
 735     format(/,' Density Function ',/,
     &        ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &        ' n(x) = n0*[fx*exp((x-a)/r0) + fp*((x/a)**b) +', 
     &               ' fv*(((r0-x)/e)**-f)]+c  ',/,
     &        ' Give parameters n0 fx fp fv a b c e f & r0 (cgs): ',$)
         read(*, *) ofac, xfac, pfac, vfac, afac, bfac, cfac, 
     &              efac, ffac, scalen
         dhn = ofac
      end if
c     
c     Estimate strom radius for nominal density
c     if we have any photons that is!
c     
c     if isobaric, get pressure regime
c     
      if (jden.eq.'B') then
 497     format(/,' Give pressure regime (p/k, <10 as log) : ',$)
         write(*, 497) 
         read(*,*) presreg
         if (presreg.le.10.d0) presreg = 10.d0**presreg
c     
 499     format(//,' Give estimate of mean temperature ',
     &        ' [T ~1e4, <10 as log]',/,
     &        ' (for initial size estimates only, actual temperature',/,
     &        ' and structure calculated later in detail, ',/,
     &        ' if in doubt use 1e4 (or 4))',/,
     &        ' : ',$)
         write(*, 499) 
         read(*,*) tinner
         if (tinner.le.10.d0) tinner = 10.d0**tinner
c     
         dhn = presreg/tinner
         dh10 = presreg/1e4
c     
 498     format(/,' *******************************************',/,
     &           '  H density at 1e4 K : ',1pg12.5,/,
     &           ' *******************************************')
         write(*, 498) dh10
c
c     jden = B isobaric
c     
      endif
c     
 490  continue
      if (jden.eq.'C') then
 491     format(/,' Give nominal hydrogen density  : ',$)
         write(*, 491) 
         read(*, *) dhn
      endif
c
      dht = densnum(dhn)
c      
 492  format(/,' ********************************************',/,
     &        '   Mean Ion Density (N):',1pg12.5,/,
     &        '   Mean H  Density (dh):',1pg12.5,/,
     &        ' ********************************************')
c
      write(*, 492) dht,dhn
c     
 493  format(/,' Give a filling factor (0<f<=1) : ',$)
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
 320     format(/,' ********************************************',/,
     &        '  Estimated hydrogen Stromgren radius of the',/,
     &        '  order of ................:',1pg10.3,' cm.',/,
     &        '  Impact parameters:QHDN in:',1pg10.3,/,
     &        '                    QHDN av:',1pg10.3,/,
     &        ' ********************************************',//,
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
 308        format(/,' Give U at inner radius (<=0 as log) : ',$)
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
 321        format(/,' Give Q at inner radius (<= 100 as log) : ',$)
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
        blum = 0.d0
        ilum = 0.d0
c     
         do i = 1,infph-1
           wid = ephot(i+1)-ephot(i)
           blum = blum + soupho(i)*wid*evplk
           if (ephot(i).ge.Ryd) ilum = ilum + soupho(i)*wid*evplk
         enddo
c
c     Plane Parallel (Sr=1/pi) x pi
c
         blum = pi*blum
         ilum = pi*ilum

         wdpl = fdilu(rstar,remp)
         qhdin = 4.d0*qht/dht*wdpl
         blum = blum*wdpl
         ilum = ilum*wdpl
         wdpl = fdilu(rstar,(rmax+remp)*0.5d0)
         qhdav = 4.d0*qht/dht*wdpl
 323     format(/,' ******************************************',/,
     &           '  Empty radius : ',1pg12.5,/,
     &           '  QHDN in      : ',1pg12.5,/,
     &           '  QHDN <r>     : ',1pg12.5,/,
     &           '  Total intensity (erg/s/cm2)   : ',1pg12.5,/,
     &           '  Ionizing intensity (erg/s/cm2): ',1pg12.5,/,
     &           ' ******************************************',//)
         write(*,323) remp,qhdin,qhdav,blum,ilum
c
         abc0 = 0.d0
         abc9 = 1.d38
c     
 60      write(*, 62) 
 62      format(/,' Volume integration over the whole sphere?',
     &        ' (y/n) : ',$)
         read(*, 350) ilgg
c     
         if ((ilgg(1:1).eq.'Y').or.(ilgg(1:1).eq.'y')) ilgg = 'Y'
         if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) ilgg = 'N'
c     
         if ((ilgg .ne. 'Y').and.(ilgg .ne. 'N')) goto 60
c     
         if (ilgg .eq. 'N') then
 921        write(*, 923) 
 923        format(/,' Integration through the line of sight ',
     &           'for a centered ring aperture.',/,
     &           'Give inner and outer radii (cm) :')
            read(*, *) abc0, abc9
            if ((abc0.lt.0.0d0).or.(abc0.ge.(0.999d0*abc9))) goto 921
         endif
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
 935     format(/,' Give Ionizing Flux at inner edge by: ',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    B  : Bolometric Flux (erg/s/cm^2)',/,
     &     '    I  : Ionising Flux (erg/s/cm^2)',/,
     &     '    U  : Ionisation parameter (Using dh)',/,
     &     '    N  : No change (Use current flux)',/,
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
 955     format(//,' *********************************************',/,
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
 410     format(/' ************************************************',/,
     &           '  Impact parameter (g=0.5)  QHDN av:',1pg12.5,/,
     &           ' ************************************************',//,
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
      if (jgeo.eq.'F') then
c
        if (qht.gt.0.d0) then
c      
         qtlog = (dlog10(qht)+1.09921d0)+(2.0d0*dlog10(rstar))
         vstrlo = ((((2.531d0+qhlo)+(2.d0*dlog(rstar)))-(2.d0*dlog(
     &        dht)))-dlog(fin))-reclo
         rmax = dexp((vstrlo-1.43241d0)/3.d0)
         qhdin = qht/dht
         qhdav = ((2.d0*qht)/dht)*fdilu(rstar,(rmax+rstar)*0.5d0)
c     
 512     write(*, 520) rmax, qhdin, qhdav
 520     format(/,' ********************************************',/,
     &        '  Estimated hydrogen Stromgren length of the',/,
     &        '  order of ................:',1pg10.3,' cm.',/,
     &        '  Impact parameters:QHDN in:',1pg10.3,/
     &        '                    QHDN av:',1pg10.3,/
     &        ' ********************************************',//,
     &        ' Give radius of empty zone in front of source ',/,
     &        ' (in cm (0,>100)) or in log(cm) (>0<=100)) : ',$)
         read(*, *) remp
c     
         if (remp.lt.0.0d0) goto 512
         if ((remp.le.1.d2).and.(remp.gt.0.d0)) remp = 10.d0**remp
c     
         rmax = dexp(((vstrlo+dlog(1.d0+((remp/rmax) ** 3)))-
     &        1.43241d0)/3.0d0)
c
        else
c     
c     give R_0 directly
c     
 1240        write(*,1250)
 1250        format(/' Give initial radius:',$)
             read(*, *) remp
c     
             if (remp.lt.0.0d0) goto 1240
             if (remp.lt.rstar) remp = rstar            
             rmax = remp*100.0d0
c
         endif
c     
         write(*,*) rstar,remp,rmax
         wdpl = fdilu(rstar,remp)
         qhdin = ((4.0d0*qht)*wdpl)/dht
         wdpl = fdilu(rstar,(rmax+remp)/2.0d0)
         qhdav = ((4.0d0*qht)*wdpl)/dht
c     
      endif
c     
 333  continue
c     
      abc3 = 0.0d0
      tm00 = 0.0d0
c     
c     ***EQUILIBRIUM OR FINITE AGE ASSUMPTION
c     
 459  format(//,' ::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '  Ionisation balance calculations',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write(*,459)
 460  write(*, 465) 
 465  format(//' Choose type',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '     E  : Equilibrium balance.',/,
     &     '     F  : Finite source lifetime.',/,
     &     '     P  : Post-Equilibrium decay.',/,
c     &     '        :'/
c     &     '     C  : CIE for given Te profile.'/
     &     '  :: ',$)
      read(*, 350) jeq
c     
      if ((jeq(1:1).eq.'E').or.(jeq(1:1).eq.'e')) jeq = 'E'
      if ((jeq(1:1).eq.'F').or.(jeq(1:1).eq.'f')) jeq = 'F'
      if ((jeq(1:1).eq.'P').or.(jeq(1:1).eq.'p')) jeq = 'P'
      if ((jeq(1:1).eq.'c').or.(jeq(1:1).eq.'c')) jeq = 'C'
c     
      if (      (jeq .ne. 'E')
     &     .and.(jeq .ne. 'F')
     &     .and.(jeq .ne. 'P'))
     &     goto 460
c     
c     forced CIE for given tprofile.
c     useful for hot x-ray bubbles.
c     
      if(jeq .eq. 'C') then 
 1200    write(*, 1210) 
 1210    format(/,' Fixed Thermal Profile:',/,
     &       ' :::::::::::::::::::::::::::::::::::::::::::::::::::',//,
     &       ' f(x) = t0*[fx*exp((x-ta)/r0) + fp*((x/ta)**tb)]+tc  ',/,
     &       ' (with r in r0 units , t0<10K as log)',/,
     &       ' Give parameters t0 ta tb tc r0 : ',$)
         read(*, *) tofac, txfac, tpfac, tafac, tbfac, tcfac, tscalen
c
      endif
c     
c     non equilibrium
c     
      if ((jeq .ne. 'E').and.(jeq .ne. 'C')) then
         rcin = 0.9
         if (jgeo .eq. 'S') then
            dti =-(dlog(1.0d0-(rcin ** 3))/(rechy*dht))
         else
            dti =-(dlog(1.0d0-rcin)/(rechy*dht))
         end if
         dti = (rmax/3.d10)+dmax1(dti,rmax/3.d10)
         dthre = (1.2/dht)/3d-13
c     
c     finite lifetime
c     
         if (jeq .eq. 'F') then
 467       write(*, 470) rmax, dti
 470       format(/' Source turn on and finite lifetime:',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     ' (Hydrogen Stromgrem radius of the order of',1pg10.3,
     &     ' cm',/,' Characteristic time scale for ionisation: ',1pg7.1,
     &    ' sec)',/,' Give elapsed time since source turned on (sec)',/,
     &     ' and give life-time of ionising source (sec) : ',$)
           read(*, *) telap, tlife
           if ((telap.le.0.0d0).or.(tlife.le.0.0d0)) goto 467
           tm00 = 100.0d0
 777       write(*, 772) 
 772       format(//' Initial temperature of neutral gas : ',$)
           read(*, *, err=777) tm00
           if (tm00.lt.1.0d0) goto 777
         else
           tlife = 0.0d0
 431       write(*, 433) rmax, dthre
 433       format(/' Source turn-off and final luminosity',/,
     &      ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &      ' (Hydrogen Stromgrem radius of the order of',1pg10.3,' cm'
     &      /,'Characteristic time scale for recombination: ',1pg7.1,
     &   ' sec)',/,' Give elapsed time since source turned off (sec)',/,
     &      ' and fractional final lumnosity of the source (>=0): ',$)
           read(*, *) telap, abc3
           if ((telap.le.0.0d0).or.(abc3.lt.0.0d0)) goto 431
         end if
      end if
c     
 843  write(*, 847) 
 847  format(/,' Step value of the photon absorption fraction: ',$)
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
     &       ' changes allowed',/,
     &       ' for three parameters: ion population, te and hydrogen',
     &       'density.',/ ,
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
      if (jeq .eq. 'E') then
      
         write(*, 623) dtau0         
 623     format(/' Summary : Thermal and Ionic Equilibrium :',/,
     &           ' -----------------------------------------',/,
     &           '  dTau        :',0pf7.4)
     
      else if (jeq .eq. 'F') then
      
         write(*, 624) telap, tlife, dtau0
 624     format(/' Summary : Non-Equilibrium + Finite Source Life :',/,
     &           ' ------------------------------------------------',/,
     &           '  Age         :',1pg10.3,' sec',/,
     &           '  Source Life :',1pg10.3,' sec',/,
     &           '  dTau        :',0pf7.4)

      else
      
         write(*, 625) telap, dtau0
 625     format(/' Summary : Equilibrium + Source Switch off :',/,
     &           ' -------------------------------------------',/,
     &           '  Switch off  :',1pg10.3,' sec',/,
     &           '  dTau        :',0pf7.4)

      end if
c
 620  format(/' Radiation Field :',/,
     &   ' -----------------',/,
     &   ' Rsou.',t15,' Remp.',t25,' Rmax',t35,' DILUav',/
     &   ,1pg10.3,t15,1pg10.3,t25,1pg10.3,t35,1pg10.3,//,
     &   ' Hdens',t15,' Ndens',t25,' F.F.',t35,' QHDNin',t45,' QHDNav',/
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
 507  format(/' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '  Boundry conditions  ',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      write(*,507)
c     
 510  format(/' Choose a model ending : ',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    A  :   Radiation bounded, HII <',2pf5.2,'%',/,
     &     '    B  :   Ionisation bounded, XII < Y',/,
     &     '    C  :   Temperature bounded, Tmin.',/,
     &     '    D  :   Optical depth limited, Tau',/,
     &     '    E  :   Density bounded, Distance',/,
     &     '    F  :   Column Density limited, atom, ion',/,
     &     '       :   ',/,
     &     '    R  :   (Reinitialise)',/,
     &     '    G  :   (Reset geometry only)',/,
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
     &        ' drops: ',/,
     &        ' (in cm (>1E6) or as a fraction of the Stromgren',/,
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
 597     format(/' Applies to ion stage',
     &         ' (>maxion for total abundance column): ',$)
         read(*, *) jpoen
         if (jpoen.le.0) goto 596
         if (jpoen.gt.maxion(ielen)) jpoen=maxion(ielen)+1
      end if
c     
 2107 format(//
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     ' Output Requirements  ',/,
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',//)
      write(*,2107)
c     
c     
c     
 2100 format(//' Choose output settings : ',/,
     &     ' :::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &     '    A  :   Standard output (photnxxxx,phapnxxxx).',/,
     &     '    B  :   Standard + monitor 4 specific ions.',/,
     &     '    C  :   Standard + all ions file.',/,
     &     '    D  :   Standard + final source + nu-Fnu spectrum.',/,
     &     '    F  :   Standard + first balance.',/,
     &     '    G  :   Everything',//,
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
      jbfx = 'p4bal'
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
         call compph5(dhn, fin, banfil,difma,dtlma,dhlma)
c     
      endif
c     
      return 
c     
      end

