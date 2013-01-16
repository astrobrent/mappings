cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORM VECTOR : TPHOT REPRESENTING THE MEAN INTENSITY
c     OF RADIATION : JNU (AT THE POINT CONSIDERED) USING
c     SOURCE VECTOR : SOUPHO(=INU AT RDIST=RSTAR).
c
c     ADDS DIFFUSE FIELD CONTAINED IN VECTORS :
c          UPDIF,DWDIF,UPLIN AND DWLIN
c     ALSO ADDS LOCAL EMISSIVITY VECTORS : EMIDIF,EMILIN
c     
c     NOTE : UNITS FOR THE VECTORS EMIDIF & EMILIN ARE 
c              IN NUMBERS OF PHOTONS
c              (INSTEAD OF ERGS LIKE ALL OTHER VECTORS)
c
c   GENERALLY SUBR. LOCALEM  NEEDS TO BE PREVIOUSLY CALLED
c   GENERALLY SUBR. NEWDIF  NEEDS TO BE PREVIOUSLY CALLED
c
c     ABSORPTION : TAU  IS DETERMINED USING THE INTEGRATED
c     COLUMN DENSITY OF MATTER : POPINT  FOR THE SOURCE VECTOR
c     AND THE COLUMN DENS. WITHIN THE SPACE STEP ITSELF FOR ALL
c     VECTORS (FROM R UP TO R+DR  USING POPUL. IN : POP(6,11))
c
c     THE GEOMETRICAL DILUTION FACTOR : WDIL  NEEDS TO BE GIVEN.
c     GENERALLY : WDIL=(RDIST**2)/(4*RSTAR**2)  (SEE FUNCTION FDILU).
c     IN THE PLANE-PARALLEL CASE, WDIL HAS TO BE EQUAL TO 0.5.
c
c     THE DILUTION FACTOR WITHIN THE SPACE STEP ITSELF IS
c     DETERMINED FROM THE RADIUS OF CURVATURE : RAD  DEFINED
c     AT THE INNER LIMIT OF THE SPACE STEP : DR
c     FOR THE UPSTREAM VECTOR, IT CORRESPONDS (BY CONVENTION FOR
c     AN INWARD FLUX) TO A CONCENTRATION FACTOR
c     (GIVE A VALUE FOR RAD OF ZERO FOR PLANE PARALLEL GEOMETRY)
c     CALL SUBR. RESTRANS,CASAB
c
c
c     LMOD='DW'   : ADDS SOURCE VECTOR AND DOWNSTREAM VECTOR.
c     LMOD='UP'   : ADDS SOURCE VECTOR AND UPSTREAM VECTOR.
c     LMOD='ALL'  : ADDS ALL COMPONENTS.
c     LMOD='SO'   : ADDS ONLY SOURCE VECTOR (ON THE SPOT APPROX.)
c                   (IN THIS CASE SUBR. LOCALEM IS UNNECESSARY).
c     LMOD='CAB'  : THE SUBROUTINE ONLY DETERMINE CASE A,B
c                   (NO PHOTON FIELD AT ALL).
c     LMOD='LOCL' : ONLY BUILD LOCAL DIFFUSE FIELD.
c     LMOD='NEBL' : ONLY INTEGRATED NEBULA EMISSION (IE NO SOURCE)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine totphot(t, dh, fi, rad, dr, dv, wd, lmod)

c     T    = Electron tempeature in K.
c     DH   = Total Hydrogen density
c     FI   = Filling factor
c     RAD  = Distance from photoionization source in cm
c     DR   = Width of space step in cm
c     DV   = Width of velocity step in cm s^-1
c     LMOD = Code for radiation components to include.

      include 'cblocks.inc'
c
c           Variables
c
      double precision t, dh, fi, rad, dr, dv, wd
      character lmod*4
      
      double precision dadw,daup,den,dismul
      double precision drh,dvh,dwf,emf,energ
      double precision swdil
      double precision totem,upf,wid
      double precision tauso, tau, sigmt, curad, wlo, ulo
      double precision tau0, dusttau0, srcf0
      double precision dwdil,rhoion,phots,z,p,f,es,casab,fb,dfb
      
      double precision pathlength
      double precision localflux
      double precision escape
      
      double precision rf0,rf1
      double precision tlin0,tlin1,thlin0,thlin1
      double precision txhlin0,txhlin1,thelin0,thelin1
c
      double precision ti, dustloss
c
      double precision dusttau, dustsigmat 
      double precision dusttrans, dusttotem 
      double precision downflux, upflux, retain
c
      double precision sum1,sum2,sum3,tlocalem
      double precision dsum1,dsum2,dsum3
c
      integer*4 inl,lin,m,n,atom,ion
      integer*4 bincount
      integer*4 line,series
c
c           Functions
c
      double precision densnum, localout, transferout
      double precision fdilu,fdismul,fbowen
c
c     constants
c
      rhoion = densnum(dh)
c
c       write(*,*) 'Totphot:',t, dh, fi, rad, dr, dv, wd, lmod
c
      if ((lmod(1:2) .ne. 'DW').and. 
     &    (lmod(1:2) .ne. 'UP').and.
     &    (lmod(1:3) .ne. 'ALL').and.
     &    (lmod(1:2) .ne. 'SO').and. 
     &    (lmod(1:4) .ne. 'LOCL').and.
     &    (lmod(1:4) .ne. 'NEBL').and.
     &    (lmod(1:3) .ne. 'CAB')) then
          write(*, 67) lmod(1:4)
 67       format(//,"MODE IMPROPERLY DEFINED IN TOTPHOT :",a4)
          stop 
      endif
      
      if ((fi.gt.1.0d0).or.(fi.le.0.0d0)) then
        write(*, 811) fi
 811    format(//,"INCONSISTENT FILLING FACTOR :",1pg10.3)
        stop 
      endif
c     
      if ((wd.gt.0.5d0).or.(wd.lt.0.0d0)) then
        write(*, 837) wd
 837    format(/,"GEOMETRICAL DILUTION FACTOR IS IMPROPER IN",
     &     " SUBR. TOTPHOT  :",1pg10.3)
        stop 
      endif
c
c     Compute distance multipliers for resonance lines and
c     determine case A,B.
c
      drh = dr*0.5d0
      dvh = dv*0.5d0
c
      do line = 1,xlines
        z = zion(xrat(line))
        p = pop(xion(line),xrat(line))
        emilin(2,line) = 1.d0
        if (xbr(line).gt.0.d0) then
         f = xfef(line)/xbr(line)
         es = xejk(line)*ev
         if ((z*p.ge.pzlimit))then
            emilin(2,line) =
     &        fdismul(t, dh, drh, dvh , xrat(line), xion(line), es,f)
         endif
        endif
      enddo
c
c
      fb = 0.d0
      if (zmap(2).gt.0)then 
      	fb = fbowen(t,dvh)
      endif
c
      do atom = 1, atypes
       do line = 1,10
        do series = 1,6
         ion = mapz(atom)
         z = zion(atom)
         p = pop(ion,atom)
         f = hydrogf(line,series)
c
c     Hydrogen.
         if (atom .eq. 1) then 
          hydlin(2,line,series) = 1.d0
          es = (LmeV/hlambda(line,series))*ev
c
          if ((z*p.ge.pzlimit).and.(hbin(line,series).ne.0).and.
     &       (series.eq.1))then
            hydlin(2,line,series) = 
     &       fdismul(t, dh, drh, dvh , atom, ion, es,f)
          endif
         endif
c         
c     Helium.
         if (atom .eq. 2) then 
          hellin(2,line,series) = 1.d0
          es = (LmeV/helambda(line,series))*ev
c
          if ((z*p.ge.pzlimit).and.(hebin(line,series).ne.0).and.
     &       (series.eq.1))then
            hellin(2,line,series) = 
     &       fdismul(t, dh, drh, dvh , atom, ion, es,f)
          endif
         endif
c         
c     Heavier elements.
         if ((atom .gt. 2).and.(series.lt.3)) then 
          xhydlin(2,line,series,atom) = 1.d0
          es = (LmeV/xhlambda(line,series,atom))*ev
c
c         if ((z*p.ge.pzlimit).and.(xhbin(line,series,atom).ne.0))then
            xhydlin(2,line,series,atom) = 
     &        fdismul(t, dh, drh, dvh , atom, ion, es,f)
c         endif
         endif
c         
        enddo
       enddo
      enddo
c
c     call casab with hydrogen Lyman gamma
c
      if (zmap(1) .ne. 0) then
           line = 3
           series = 1
           caseab(1) = casab(zmap(1),line,series)
      endif
c
c     call casab with Helium Lyman gamma
c
      if (zmap(2) .ne. 0) then
           line = 3
           series = 1
           caseab(2) = casab(zmap(2),line,series)
      endif
c
c     Assume heavy series are case A
c
      do atom = 3,atypes
         caseab(atom) = 0.0d0
      enddo
c      
      if (lmod .eq. 'CAB') return 
c
c     Set mode and internal dilution factors.
c
      upf = 1.0d0
      dwf = 1.0d0
      emf = 1.0d0
      if ((lmod .ne. 'UP').and.(lmod .ne. 'ALL')) upf = 0.0d0
      if ((lmod .ne. 'DW').and.(lmod .ne. 'ALL')) dwf = 0.0d0
c
      if((lmod.eq.'LOCL')) dwf = 1.0d0
      if((lmod.eq.'NEBL')) dwf = 1.0d0
c
      if (lmod .eq. 'SO') emf = 0.0d0
      if (lmod .eq. 'NEBL') emf = 0.0d0
c
      curad = rad
      dadw = 1.0d0
      daup = 1.0d0

      if ((curad.gt.0.d0).and.(jgeo.eq.'S')) then 
        wlo = (2.0d0*dlog(curad))-(2.0d0*dlog(curad+dr))
        ulo = (2.0d0*dlog(curad+dr))-(2.0d0*dlog(curad))
        dadw = dexp(wlo)
        daup = dexp(ulo)
      endif
c
      dwdil = 1.d0
      if (jgeo.eq.'F') then
        dwdil = fdilu(rshock,rad-dilradpho(inl))
      endif
c
c   Additional source dilution because sources are Inu (ie 1/pi)
c   and wdil = 0.5 in plane parallel (1/2pi) and phion integrates
c   4 pi, so we need an additional 0.5 for sources on one side
c   due to units being used (so phion effectively integrates Inu*0.5 
c   sources in pp over only 2pi radians not 4). Local 1/4pi fields
c   are integrated normally and lines have their distance mulipliers
c   and are not strictly one sided.  Handled implicitly in spherical models
c   by fdilu.
c
      swdil = 1.d0
      if (jgeo.eq.'P') then
        swdil = 0.5d0
      endif
c         
c         write(*,*) 'dils',wd,dwdil,swdil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Form vector TPHOT bin by bin.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      bincount = 0
	dustloss = 0.0d0
	ti = 0.0d0
        sum1 = 0.d0
        dsum1= 0.d0
        sum2 = 0.d0
        dsum2= 0.d0
        sum3 = 0.d0
        dsum3= 0.d0
        tlocalem=0.d0
c
      do 500 inl = 1, infph-1   ! Main loop over photon energies.

         den = 0.5d0*(ephot(inl)+ephot(inl+1))
         wid = (ephot(inl+1)-ephot(inl))
         
         energ      = den*ev
         tauso      = 0.d0
         sigmt      = 0.d0
         dustsigmat = 0.d0
         
         call crosssections(inl, tauso, sigmt, dustsigmat)
         
         sigmt      = dh*sigmt  ! total crossection, including dust
         dustsigmat = dh*dustsigmat  ! just the dust component
            
         pathlength = dr*fi
         tau0 = pathlength*sigmt
         dusttau0 = pathlength*dustsigmat  

c     XSEC is use to determine "skip bin" for photoionization
c     calculations only, don't include dust in these.

         xsec(inl) = (sigmt-dustsigmat)*fi

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Add diluted source field and integrated diffuse (source) field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Downstream, source attenuated by local opacity.
c
         srcf0 = 0.d0
         rf0   = 0.d0
         rf1   = 0.d0
         
         tau     = tau0       ! tau = 0 if dr = 0.d0
         dusttau = dusttau0
         
         dusttrans = 0.d0
c         
         if ((lmod .ne. 'LOCL')) then 
c
c      Source field if not local and not nebula only
c
             if (lmod .ne. 'NEBL') then
                srcf0 = wd*swdil*soupho(inl)*dexp(-tauso)
             endif
c
c     Integrated diffuse field and source
c
             downflux  = dwf*dadw*dwdil*dwdif(inl)+srcf0
             upflux    = upf*daup*updif(inl)
c             
             escape    = transferout(1.d0, tau)
c             
             rf0       = downflux*escape
             rf1       = upflux
c
             retain    = 1.d0-transferout(1.d0, dusttau)
             dusttrans = dusttrans + downflux*retain
c
c             if (inl.gt.(infph*0.75d0)) then
c              sum3=sum3+(downflux)*(1.d0-escape)*wid*evplk
c              dsum3=dsum3+(downflux)*retain*wid*evplk
c             else if (inl.gt.(infph*0.625d0)) then
c              sum2=sum2+(downflux)*(1.d0-escape)*wid*evplk
c              dsum2=dsum2+(downflux)*retain*wid*evplk
c             else if (inl.gt.(infph*0.5d0)) then
              sum1=sum1+(downflux)*(1.d0-escape)*wid*evplk
              dsum1=dsum1+(downflux)*retain*wid*evplk
c             endif
c               
         endif
c         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Now Add Integrated Resonance Lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Metal Resonance Lines
c     
         tlin0 = 0.0d0
         tlin1 = 0.0d0
         
         if (lmod .ne. 'LOCL') then 
           do n = 1,xlines
            if ((xbin(n).eq.inl)) then
            
               lin    = n

               downflux  = dwf*dadw*dwdil*dwlin(1,lin)
               upflux    = upf*daup*uplin(1,lin)
c
               dismul    = emilin(2,lin)
               tau       = tau0
               escape    = transferout(dismul, tau)
c
               tlin0     = tlin0     + downflux*escape
               tlin1     = tlin1     + upflux
c
               dusttau   =  dusttau0
               retain    = 1.d0-transferout(dismul, dusttau)
               dusttrans = dusttrans + downflux*retain
c               
               sum1 = sum1 +downflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +downflux*retain*wid*evplk
c
            endif
           enddo
         endif
c
c   Hydrogen
c
         thlin0 = 0.0d0
         thlin1 = 0.0d0
c
         if (lmod .ne. 'LOCL') then 
c
          do line = 1,10
           do series = 1,6
            if (hbin(line,series).eq.inl) then
c            
               downflux = dwf*dadw*dwdil*hyddwlin(1,line,series)
               upflux   = upf*daup*hyduplin(1,line,series)
          
               dismul   = hydlin(2,line,series)
               tau      = tau0
               escape   = transferout(dismul, tau)

               thlin0   = thlin0 + downflux*escape
               thlin1   = thlin1 + upflux
c
               dusttau   =  dusttau0
               retain    = 1.d0-transferout(dismul, dusttau)
               dusttrans = dusttrans + downflux*retain
c               
               sum1 = sum1 + downflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +downflux*retain*wid*evplk
c               sum1 = sum1 + downflux*retain*wid*evplk
c               
            endif
           enddo
          enddo
         endif
c
c
c   Helium
c
         thelin0 = 0.0d0
         thelin1 = 0.0d0
c
         if (lmod .ne. 'LOCL') then 
c
          do line = 1,10
           do series = 1,6
            if (hebin(line,series).eq.inl) then
            

               downflux = dwf*dadw*dwdil*heldwlin(1,line,series)
               upflux   = upf*daup*heluplin(1,line,series)
               
               dismul   = hellin(2,line,series)
               tau      = tau0
               escape   = transferout(dismul, tau)

               thelin0  = thelin0 + downflux*escape
               thelin1  = thelin1 + upflux

               dusttau   =  dusttau0
               retain    = 1.d0-transferout(dismul, dusttau)
               dusttrans = dusttrans + downflux*retain
c               
               sum1 = sum1 +downflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +downflux*retain*wid*evplk
c
c
            endif
           enddo
          enddo
         endif
c
c   Hydrogenic heavies
c
         txhlin0     = 0.0d0
         txhlin1     = 0.0d0
c
         if (lmod .ne. 'LOCL') then 
c
          do atom = 3,atypes
           if ((xhyddwlin(1,1,1,atom) .gt. epsilon)
     &      .or.(xhyduplin(1,1,1,atom) .gt. epsilon)) then
            do series = 1,2
             do line = 1,10
              if (xhbin(line,series,atom).eq.inl) then
            
               downflux = dwf*dadw*dwdil*xhyddwlin(1,line,series,atom)
               upflux   = upf*daup*xhyduplin(1,line,series,atom)
               
               dismul   = xhydlin(2,line,series,atom)
               tau      = tau0
               escape   = transferout(dismul, tau)

               txhlin0  = txhlin0 + downflux*escape
               txhlin1  = txhlin1 + upflux
               
               dusttau   = dusttau0
               retain    = 1.d0-transferout(dismul, dusttau)
               dusttrans = dusttrans + downflux*retain
c               
               sum1 = sum1 +downflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +downflux*retain*wid*evplk
c
              
              endif
             enddo
            enddo
           endif
          enddo
         endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Now add local diffuse field if on-the-spot is not used
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***ADS LOCAL EMISSIVITY CONTRIBUTION TO THE NEW DIFFUSE FIELD
c     *NB. ONLY WHEN THE ON THE SPOT APPROX. IS NOT USED : JSPOT
c     
c     Local field intensity is the fraction of the local field that
c     is  passed onto the next zone , ie (1-exp(-tau))/tau
c
c
         totem     = 0.0d0
         dusttotem = 0.0d0
c
         if ((emf .eq. 0.0d0).or.(jspot .eq. 'YES')) goto 490
c
         pathlength = (dr*fi)
         energ      = den*ev
         localflux  = energ*emidif(inl)*pathlength
c   
         tau     = tau0
         escape  = localout(1.d0, tau)
         totem   = totem     + localflux*escape

         dusttau    = dusttau0
         retain     = 1.d0-localout(1.d0, dusttau)
         dusttotem  = dusttotem + localflux*retain
c               
         sum1 = sum1 +localflux*(1.d0-escape)*wid*evplk
         dsum1 = dsum1 +localflux*retain*wid*evplk
c         sum2 = sum2 + localflux*escape*wid*evplk
         tlocalem=tlocalem+localflux*escape*wid*evplk
c
c   Resonance Lines
c
         do m = 1,xlines
            if ((xbin(m).eq.inl).and.(emilin(1,m).gt.epsilon)) then
               lin     = m
               
               energ     = ev*xejk(lin)
               localflux = energ*emilin(1,lin)*pathlength

               tau     = tau0
               dismul  = emilin(2,lin)
               escape  = localout(dismul, tau)
               totem   = totem     + localflux*escape
               
               dusttau    = dusttau0
               retain     = 1.d0-localout(dismul,dusttau)
               dusttotem  = dusttotem + localflux*retain
c               
               sum1 = sum1 +localflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +localflux*retain*wid*evplk
c               sum3 = sum3 + localflux*escape*wid*evplk
               tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
         enddo
c
c   Hydrogen
c
         do series = 1,6
           do line = 1,10
            if ((hbin(line,series).eq.inl).and.
     &          (hydlin(1,line,series).gt.epsilon))then
               energ     = ev*LmeV/hlambda(line,series)
               localflux = energ*hydlin(1,line,series)*pathlength
               
               tau        = tau0
               dismul     = hydlin(2,line,series)
               escape     = localout(dismul, tau)
               totem      = totem     + localflux*escape

               dusttau    =  dusttau0
               retain     = 1.d0-localout(dismul,dusttau)
               dusttotem  = dusttotem + localflux*retain
c               
               sum1 = sum1 +localflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +localflux*retain*wid*evplk
c               dsum2 = dsum2 + localflux*escape*wid*evplk
               tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
           enddo
         enddo
c
c
c   Helium
c
         do series = 1,6
           do line = 1,10
            if ((hebin(line,series).eq.inl).and.
     &          (hellin(1,line,series).gt.epsilon))then
c                    
               dfb = 1.d0
c               if ((line.eq.1).and.(series.eq.1)) then
c                 dfb = 1.d0-fb
c                endif

               energ = ev*LmeV/helambda(line,series)
               localflux = dfb*energ*hellin(1,line,series)*pathlength
              
               tau    = tau0
               dismul = hellin(2,line,series)
               escape = localout(dismul, tau)
               totem  = totem     + localflux*escape

               dusttau   = dusttau0
               retain    = 1.d0-localout(dismul, dusttau)
               dusttotem = dusttotem + localflux*retain
c               
               sum1 = sum1 +localflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +localflux*retain*wid*evplk
c               dsum3 = dsum3 + localflux*escape*wid*evplk
               tlocalem=tlocalem+localflux*escape*wid*evplk
c
            endif
           enddo
         enddo
c
c   Heavy Hydrogen
c
         do atom = 3,atypes
          if ((xhydlin(1,1,1,atom).gt.epsilon)) then
           do series = 1,2
            do line = 1,10
             if (xhbin(line,series,atom).eq.inl) then
               energ     = ev*LmeV/xhlambda(line,series,atom)
               localflux = energ*xhydlin(1,line,series,atom)*pathlength
               
               tau       = tau0
               dismul    = xhydlin(2,line,series,atom)
               escape    = localout(dismul, tau)
               totem     = totem     + localflux*escape

               dusttau   = dusttau0
               retain    = 1.d0-localout(dismul,dusttau)
               dusttotem = dusttotem + localflux*retain
c               
               sum1 = sum1 +localflux*(1.d0-escape)*wid*evplk
               dsum1 = dsum1 +localflux*retain*wid*evplk
c               dsum3 = dsum3 + localflux*escape*wid*evplk
               tlocalem=tlocalem+localflux*escape*wid*evplk
c
             endif
            enddo
           enddo
          endif
         enddo
c
 490     continue  ! Skip to here if EMF = 0 or JSPOT = yes.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
         totem         = totem*dwf
         dusttotem     = dusttotem*dwf
c
c******************************************************************

c     Further attenuate the downstream field by local x-sect.
c     Upstream field is already attenuated by popint in newdif.

c     RF0 is the downward flux(???).
c     RF1 is the upward flux(???).

         tphot(inl) = rf0+rf1
         
c     Add metal resonance lines escaped downward flux (TLIN0) and 
c     upward flux (TLIN1).

         tphot(inl) = tphot(inl)+tlin0+tlin1

c     Add hydrogen lines escaped downward flux (THLIN0) and 
c     upward flux (THLIN1).

         tphot(inl) = tphot(inl)+thlin0+thlin1

c     Add helium lines escaped downward flux (THELIN0) and 
c     upward flux (THELIN1).

         tphot(inl) = tphot(inl)+thelin0+thelin1

c     Add hydrogenic heavies lines escaped downward flux (TXHLIN0)
c     and upward flux (THLIN1).

         tphot(inl) = tphot(inl)+txhlin0+txhlin1

c     Add local emissivity contribution to the diffuse field (TOTEM).

         tphot(inl) = tphot(inl)+totem

c     Add component of the downward flux absorbed by dust(???) and ???.

         dustphot(inl) = dusttrans+dusttotem

c     Put local fraction of crossections into dustphot.

         dustloss = dustloss+dustphot(inl)*wid*evplk

         ti = ti+tphot(inl)*wid*evplk

c******************************************************************
c
         if (lmod .ne. 'LOCL') then 
c
c     ephot is in eV
c     rhoion = densnum(dh)
c     energ = mean energy of bin in ergs
c     
           energ = den*ev
           phots = (tphot(inl)/energ)*evplk*wid
c
           skipbin(inl) = .FALSE.
           bincount = bincount +1
           if (((xsec(inl)/dh)*phots).le.epsilon) then
            skipbin(inl) = .TRUE.
            bincount = bincount-1
           endif
c
           if (photonmode.eq.0) skipbin(inl) = .TRUE.
c         
         endif
c     
 500  enddo  ! End of main loop over photon energies
c     
c     ***FORCE RECALCULATION OF PHOTOIONISING RATES
c     
c 
c      if (grainmode) then 
c          write(*,1000) 4.d0*pi*dustloss, 4.d0*pi*ti, 4.d0*pi*(ti+dustloss)
c 1000     format('Tphot Intensity ',3(1pg14.7))		 
c      endif
c        write(*,1000) fpi*sum1,fpi*dsum1,fpi*tlocalem,fpi*ti
 1000   format('Tphot Intensity ',4(1pg14.7))		 
c        write(*,1010) fpi*sum2,fpi*sum3,fpi*(dsum2),fpi*dsum3
c 1010   format('local Intensity ',4(1pg14.7))		 

      if (lmod .ne. 'LOCL') then 
        ipho = ipho+1
      endif
c
c     Jspec disabled
c
      return 
      end

