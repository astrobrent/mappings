cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO COMPUTE NEW DIFFUSE FIELD VECTORS AT THE END OF SPACE
c     STEP,USING LOCAL EMISSIVITY PRODUCED BY SUBR. LOCALEM)
c     
c     DOWNSTREAM AND UPSTREAM TEMPERATURES : TDW,TUP
c     DOWNSTREAM AND UPSTREAM VARIATION OF DISTANCE :DRDW,DRUP
c     DOWNSTREAM AND UPSTREAM VARIATION OF VELOCITY : DVDW,DVUP
c     REM: FOR THE NON LOCAL COMPONENT,THESE QUANTITIES REFER
c     TO THE TOTAL VARIATION OF T,DR,DV THAT TAKES
c     PLACE BETWEEN THE ORIGIN AND THE POINT CONSIDERED
c     FRDW : FRACTION OF LOCAL EMISSIVITY THAT IS ADDED TO THE
c     DOWNSTREAM COMPONENT;THE COMPLEMENT IS ADDED TO
c     THE UPSTREAM VECTORS
c     
c     THE DILUTION FACTOR IN CASE OF SPHERICAL SYMMETRY IS
c     DERIVED FROM THE CURVATURE RADIUS AT THE INNER LIMIT
c     OF THE SPACESTEP(=LOCAL DR)  FOR PLANE PARALLEL SYMETRY
C     THE RAD (DISTANCE) IS IGNORED, FOR CYLINDERS IT IS USED
C     TO FORM WEIGHTED SUMS FOR DIFFUSE FIELD DILUTIONS
c     BY CONVENTION,THE DOWNSTREAM VECTORS ARE RESERVED FOR
c     THE OUTWARD (DIRECTION OF INCREASING R) FLUX AND THE
c     THE UPSTREAM ONES FOR THE INWARD FLUX
c     
c     NOTE : UNITS FOR THE VECTORS EMIDIF,EMILIN ARE IN NUMBERS
c     OF PHOTONS (INSTEAD OF ERGS LIKE ALL OTHER VECTORS)
c     
c     
c     CALL SUBR.RESTRANS
c     
c     JMOD='LODW' MEANS THAT ONLY THE DOWNSTREAM DIFFUSE FIELD
c     IS LOCAL.THE FRACTION ADDED TO THE UPSTREAM
c     COMPONENT IS FURTHER ATTENUATED USING
c     THE INTEGRATED COLUMN DENSITIES : POPINT
c     JMOD='LOUP' REVERSED SITUATION TO THE PREVIOUS CASE
c     JMOD='OUTW' OUTWARD ONLY APPROXIMATION WHERE ALL THE
c     LOCAL EMISSIVITY IS ADDED TO THE DOWNSTREAM
c     COMPONENT;NO UPSTREAM DIFFUSE FIELD
c     JMOD='DWUP' BOTH STREAMS ARE LOCAL
c     
c     
c     
      subroutine newdif (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup
     &                   ,frdw, jmod)
c     
      include 'cblocks.inc'
c
c           Variables
c
      double precision aph,at,bet,crosec,den,dildw
      double precision dilup,dismul,dr,drem,dremh,dvem,dvemh,dwem
      double precision dwex, dwf, dwliem, energ, eph
      double precision se, sigmul,telc
      double precision upem, upex, upf, upliem, wadw, waup, wedw, weup
      double precision tauso, tau, sigmt, curad, rm, wlo, ulo
      double precision tdw, tup, dh, fi, rad, drdw, dvdw, drup
      double precision dvup,frdw,emis,z,p,f,es,dwsc,upsc,pahsum
      double precision dfbu,dfbd,dfb
      double precision fbu,fbd,fb
      double precision drdwh,dvdwh,drduph,dvduph
c
      integer*4 line,series,i,j, atom,ion
      integer*4 ie,inl,lin,m,n,dtype
c
      character jmod*4
c
c     old
c     double precision a,aa,dtype
c     integer*4 k,nz
c           Functions
c
      double precision acrs,fdismul,fbowen
c     
      acrs(eph,at,bet,se) = (at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c     
c      
c     Check parameters and replace with sane values if necessary   
c      
c
c filling factor 
c
      if ((fi.gt.1.0d0).and.(fi.le.0.0d0)) then
           write(*, 811) fi
 811      format(//' INCONSISTENT FILLING FACTOR IN NEWDIF: ',1pg10.3)
          fi = 1.0
      endif
c 
      if ((tdw.lt.0.d0)
     &     .or.(tup.lt.0.d0)
     &     .or.(rad.lt.0.d0)
     &     .or.(drdw.lt.0.d0)
     &     .or.(drup.lt.0.d0)
     &     .or.(dvdw.lt.0.d0)
     &     .or.(dvup.lt.0.d0)) then
     
        write(*, 873) tdw,tup,rad,drdw,drup,dvdw,dvup
        
 873    format(/,/,'ONE OF THE ARGUMENTS HAS A NEGATIVE VALUE IN NEWDIF'
     &         ,'****',/,7(1pg10.3))
 
 
        tdw = dabs(tdw)
        tup = dabs(tup)
        rad = dabs(rad)
        drdw = dabs(drdw)
        drup = dabs(drup)
        dvdw = dabs(dvdw)
        dvup = dabs(dvup)
        
      endif
c     
      if ((jmod .eq. 'LODW')
     &     .or.(jmod .eq. 'LOUP')
     &     .or.(jmod .eq. 'OUTW')
     &     .or.(jmod .eq. 'DWUP')) goto 60
     
      write(*, 67) jmod
 67   format(//' MODE IMPROPERLY DEFINED IN NEWDIF : ',a4)
      stop 
      
 60   continue
c     
      if ((frdw.gt.1.d0).or.(frdw.lt.0.d0)) then
             write(*, 77) frdw
 77         format(/,/,' INCONSISTENT VALUE FOR THE ARGUMENT FRDW :'
     &             ,1pg10.3)
            frdw = 0.5d0
      endif
c     
c     ***SET MODE AND INTERNAL DILUTION FACTORS
c     
      dwf = frdw
      if (jmod .eq. 'OUTW') dwf = 1.0d0
      upf = 1.d0-dwf
c
      dwex = 1.d0
      if (jmod .eq. 'LOUP') dwex = 0.d0
      upex = 1.d0-dwex
      if (jmod .eq. 'DWUP') upex = 1.d0
c
      if (jmod .eq. 'LOUP') then
        dr = drup
        drem = drup
        dvem = dvup
        telc = tup
      else
        dr = drdw
        drem = drdw
        dvem = dvdw
        telc = tdw
      endif
c
c     Plane Parallel / Finite Cylinder Default
c
      curad = rad
      dildw = 1.d0
      dilup = 1.d0
      
      if (jgeo.eq.'S') then
c
c     Spherical
c      
          rm = curad+drdw
          wlo = (2.d0*dlog(curad))-(2.d0*dlog(rm))
          rm = (curad+dr)-drup
      
          if (rm.le.0.d0) then
              write(*, 123) drup, curad
 123        format(/,/,' DRUP IS LARGER THAN THE RADIUS OF CURVATURE :'
     &              ,2(1pg10.3))
              stop
          endif
           
          ulo = (2.d0*dlog(curad+dr))-(2.d0*dlog(rm))
          dildw = dexp(wlo)
          dilup = dexp(ulo)
c
c     Spherical
c      
      endif
c     
      wadw = dildw
      wedw = 1.d0
      waup = 1.d0
      weup = dilup
      if (jmod .ne. 'LOUP') goto 78
      wadw = 1.d0
      wedw = dildw
      waup = dilup
      weup = 1.d0
 78   if (jmod .ne. 'DWUP') goto 71
      wadw = dildw
      wedw = 1.d0
      waup = dilup
      weup = 1.d0
 71   continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         Begin Routine Proper
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ***COMPUTES RESONANCE-LINES DISTANCE MULTIPLIERS
c     
      dremh = drem/2.d0
      dvemh = dvem/2.d0
c     
      do line = 1,xlines
      
         z = zion(xrat(line))
         p = pop(xion(line),xrat(line))
c         
         uplin (2,line) = 1.d0
         dwlin (2,line) = 1.d0
         emilin(2,line) = 1.d0
c
         if (xbr(line).gt.epsilon) then
         f = xfef(line)/xbr(line)
         es = xejk(line)*ev
c     
         if ((z*p.ge.pzlimit))then
            emilin(2,line) = fdismul(telc, dh, dremh, dvemh, 
     &                       xrat(line), xion(line), es,f)
         endif
         endif
c
       enddo
c     
c   Hydrogenic lines
c
          fb  = 0.0d0
          fbu = 0.0d0
          fbd = 0.0d0
          if (zmap(2).gt.0)then 
            fb  = fbowen(telc,dvemh)
            fbu = fbowen(tup,dvup)
            fbd = fbowen(tdw,dvdw)
          endif
          
          heiioiiibf = fb
c
c	loop in atomic number
c
         do atom = 1, atypes
         do series = 1,6
         do line = 1,10
         
c           
            ion  = mapz(atom)
            z = zion(atom)
            p = pop(ion,atom)
            f = hydrogf(line,series)
c            
            if (mapz(atom) .eq. 1) then
              hyduplin(2,line,series) = 1.d0
              hyddwlin(2,line,series) = 1.d0
              hydlin(2,line,series) = 1.d0
              j = hbin(line,series)
              if (j.ne.0) then
              es = LmeV/hlambda(line,series)*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 hydlin(2,line,series) = fdismul(telc, dh, dremh, 
     &                                   dvemh, atom, ion, es,f)
              endif
              endif
            endif
c            
            if (mapz(atom) .eq. 2) then
              heluplin(2,line,series) = 1.d0
              heldwlin(2,line,series) = 1.d0
              hellin(2,line,series) = 1.d0
              j = hebin(line,series)
              if (j.ne.0) then
              es = LmeV/helambda(line,series)*ev
c
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 hellin(2,line,series) = fdismul(telc, dh, dremh
     &                                   ,dvemh, atom, ion, es,f)
              endif
              endif
            endif
c
            if ((mapz(atom).gt.2).AND.(series.lt.3)) then
              xhyduplin(2,line,series,atom)= 1.d0
              xhyddwlin(2,line,series,atom)= 1.d0
              xhydlin(2,line,series,atom)= 1.d0
              j = xhbin(line,series,atom)
              if (j.ne.0) then
              es = LmeV/xhlambda(line,series,atom)*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 xhydlin(2,line,series,atom) = fdismul(telc, dh
     &                               , dremh, dvemh, atom, ion, es,f)
              endif
              endif
            endif
c
         enddo
         enddo
         enddo
c     
c     do up and down if needed
c     
      if (dwf.gt.0.d0) then
         drdwh=0.5d0*drdw
         dvdwh=0.5d0*dvdw
c
c first the general resonance lines     
c
         do line = 1,xlines
            z = zion(xrat(line))
            p = pop(xion(line),xrat(line))
            dwlin(2,line) = 1.d0

            if (xbr(line).gt.0.0) then
            f = xfef(line)/xbr(line)
            es = xejk(line)*ev
            if ((z*p.ge.pzlimit))then
               dwlin(2,line) = fdismul(tdw, dh, drdwh, dvdwh
     &                         , xrat(line), xion(line), es,f)
            endif
            endif
         enddo
c     
c
c   Hydrogenic lines
c
         do atom = 1,atypes
         do series = 1,6
         do line = 1,10
c         
            ion  = mapz(atom)
            z = zion(atom)
            p = pop(ion,atom)
            f = hydrogf(line,series)
c            
            if (mapz(atom) .eq. 1) then
              hyddwlin(2,line,series) = 1.d0
              j = hbin(line,series)
              if (j.ne.0) then
              es = (LmeV/hlambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 hyddwlin(2,line,series) = fdismul(tdw, dh, drdwh
     &                                 , dvdwh, atom, ion, es,f)
              endif
              endif
            endif
c            
            if (mapz(atom) .eq. 2) then
              heldwlin(2,line,series) = 1.d0
              j = hebin(line,series)
              if (j.ne.0) then
              es = (LmeV/helambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 heldwlin(2,line,series) = fdismul(tdw, dh, drdwh
     &                                     , dvdwh, atom, ion, es,f)
              endif
              endif
            endif
c
            if ((mapz(atom).gt.2).and.(series.lt.3)) then
              xhyddwlin(2,line,series,atom) = 1.d0
              j = xhbin(line,series,atom)
              if (j.ne.0) then
              es = LmeV/xhlambda(line,series,atom)*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 xhyddwlin(2,line,series,atom) = fdismul(tdw, dh
     &                                , drdwh, dvdwh, atom, ion, es,f)
              endif
              endif
            endif
c
         enddo
         enddo
         enddo
c
      endif
c     
      if (upf.gt.0.d0) then
         druph=0.5d0*drup
         dvuph=0.5d0*dvup
c
c first the general resonance lines     
c
         do line = 1,xlines
c         
            z = zion(xrat(line))
            p = pop(xion(line),xrat(line))
            uplin(2,line) = 1.d0
            if (xbr(line).gt.0.0) then
            f = xfef(line)/xbr(line)
            es = xejk(line)*ev
            if ((z*p.ge.pzlimit))then
               uplin(2,line) = fdismul(tup, dh, druph, dvuph, xrat(line)
     &                         , xion(line), es,f)
            endif
            endif
c
         enddo
c
c   Hydrogenic lines
c
         do atom = 1,atypes
         do series = 1,6
         do line = 1,10
c         
            ion  = zmap(atom)
            z = zion(atom)
            p = pop(ion,atom)
            f = hydrogf(line,series)
c            
            if (mapz(atom) .eq. 1) then
              hyduplin(2,line,series) = 1.d0
              j = hbin(line,series)
              if (j.ne.0) then
              es = (LmeV/hlambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 hyduplin(2,line,series) = fdismul(tup, dh, druph, dvuph
     &                      , atom, ion, es,f)
              endif
              endif
            endif
c            
            if (mapz(atom) .eq. 2) then
              heluplin(2,line,series) = 1.d0
              j = hebin(line,series)
              if (j.ne.0) then
              es = (LmeV/helambda(line,series))*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 heluplin(2,line,series) = fdismul(tup, dh, druph, dvuph
     &                      , atom, ion, es,f)
              endif
              endif
            endif
c
            if ((mapz(atom).gt.2).and.(series.lt.3)) then
              xhyduplin(2,line,series,atom) = 1.d0
              j = xhbin(line,series,atom)
              if (j.ne.0) then
              es = LmeV/xhlambda(line,series,atom)*ev
              if ((z*p.ge.pzlimit).and.(series.eq.1))then
                 xhyduplin(2,line,series,atom) = fdismul(tup, dh, druph
     &                                         , dvuph, atom, ion, es,f)
              endif
              endif
            endif
c
         enddo
         enddo
         enddo

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***UPDATES DIFFUSE FIELD VECTORS UPDIF, DWDIF etc BIN BY BIN
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 500 inl = 1, infph-1
c     
         den = 0.5d0*(ephot(inl)+ephot(inl+1))
         energ = den*ev
         tauso = 0.d0
         sigmt = 0.d0
c     
c     ***DERIVES ABSORPTION CROSS-SECTION FOR BIN*INL
c     
         do i = 1, ionum
            ie = atpho(i)
            if (ie.le.limph) then
               j = ionpho(i)
               eph = den/ipotpho(i)
               if (eph.ge.1.d0) then
                  aph = zion(ie)*pop(j,ie)
                  if (aph.gt.pzlimit) then
                     crosec = acrs(eph,sigpho(i),betpho(i),spho(i))
                     tauso = tauso+(popint(j,ie)*crosec)
                     sigmt = sigmt+(aph*crosec)
                  endif
               endif
            endif
         enddo
c     
         if (grainmode) then 
            if (pahmode) then
c     
c     pahs
c     
               do m=1,pahi
                 if (pahion(m).gt.0) then
                   crosec = pahiext(inl)*pahfrac
                 else
                   crosec = pahnext(inl)*pahfrac
                 endif
                 tauso = tauso+pahint(m)*crosec
                 sigmt = sigmt + pahZ(m)*crosec
               enddo
c     
c     end pahs
c     
            endif 
c     
c     add dust next
c
	 if (.NOT.ClinPAH) then !C not linked to PAH or no PAHs
 	   do dtype=1,numtypes
            tauso = tauso+dustint*dcrosec(inl,dtype)
            sigmt = sigmt+dcrosec(inl,dtype)
 	   enddo
	 else
           if (pahmode) then
             do m=1,pahi
               tauso= tauso+pahint(m)*dcrosec(inl,1) !pahint=dustint for C
             enddo
             sigmt = sigmt+dcrosec(inl,1)
           endif
 	   do dtype=2,numtypes
            tauso = tauso+dustint*dcrosec(inl,dtype)
            sigmt = sigmt+dcrosec(inl,dtype)
 	   enddo
	 endif	   
c     
c     end dust
c     
         endif
c     
c     ***ATTENUATION OF THE ORIGINAL INTENSITY OF THE DIFFUSE FIELD
c     
         sigmt = dh*sigmt
         tau   = dr*fi*sigmt
c     
         dwdif(inl) = dwdif(inl)*(wadw*dexp(-(dwex*tau)))
         updif(inl) = updif(inl)*(waup*dexp(-(upex*tau)))
c     
c     distance multipliers attenuate lines for emilin etc
c
         do n = 1,xlines
            if (xbin(n) .eq. inl) then
               lin    = n
               dismul = (dwex*dwlin(2,lin))+(upex*uplin(2,lin))
               sigmul = sigmt*dismul
               tau    = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                              tau = (dr*fi)*sigmul
               dwlin(1,lin) = dwlin(1,lin)*(wadw*dexp(-(dwex*tau)))
               uplin(1,lin) = uplin(1,lin)*(waup*dexp(-(upex*tau)))
            endif
         enddo
c
c   Hydrogenic lines
c
         do series = 1,6
         do line = 1,10
            if (hbin(line,series) .eq. inl) then
               dismul = (dwex*hyddwlin(2,line,series))
     &                  +(upex*hyduplin(2,line,series))
               sigmul = sigmt*dismul
               tau = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                  tau = (dr*fi)*sigmul
               hyddwlin(1,line,series) = hyddwlin(1,line,series)
     &                  *(wadw*dexp(-(dwex*tau)))
               hyduplin(1,line,series) = hyduplin(1,line,series)
     &                  *(waup*dexp(-(upex*tau)))
            endif
         enddo
         enddo
c
c   Helium lines
c
         do series = 1,6
         do line = 1,10
            if (hebin(line,series) .eq. inl) then
            dismul = (dwex*heldwlin(2,line,series))
     &               +(upex*heluplin(2,line,series))
            sigmul = sigmt*dismul
            tau = 1.d10
            if ((tau/dismul).gt.((dr*fi)*sigmt)) tau = (dr*fi)*sigmul
c               
            dfbu = 1.d0
            dfbd = 1.d0
            if ((line.eq.1).and.(series.eq.1)) then
                 dfbu = 1.d0-fbu
                 dfbd = 1.d0-fbd
            endif
c
            heldwlin(1,line,series) = dfbd*heldwlin(1,line,series)
     &                                *(wadw*dexp(-(dwex*tau)))
            heluplin(1,line,series) = dfbu*heluplin(1,line,series)
     &                                *(waup*dexp(-(upex*tau)))
            endif
         enddo
         enddo
c
         do atom = 3,atypes
         do series = 1,2
         do line = 1,10
            if (xhbin(line,series,atom) .eq. inl) then
               dismul = (dwex*xhyddwlin(2,line,series,atom))
     &                  +(upex*xhyduplin(2,line,series,atom))
               sigmul = sigmt*dismul
               tau = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                  tau = (dr*fi)*sigmul
               xhyddwlin(1,line,series,atom) = 
     &           xhyddwlin(1,line,series,atom)*(wadw*dexp(-(dwex*tau)))
               xhyduplin(1,line,series,atom) = 
     &           xhyduplin(1,line,series,atom)*(waup*dexp(-(upex*tau)))
            endif
         enddo
         enddo
         enddo
c     
c     ***ADS LOCAL EMISSIVITY CONTRIBUTION TO THE NEW DIFFUSE FIELD
c     
         tau = (dr*fi)*sigmt
c
         emis =  emidif(inl)*dr*fi*energ
c         
         dwem = (dwf*emis)
         upem = (upf*emis)
c     
c     
c     Add local diffuse field through this zone
c     
         dwdif(inl) = dwdif(inl)+(dwem*(wedw*dexp(-(upex*tau))))
         updif(inl) = updif(inl)+(upem*(weup*dexp(-(dwex*tau))))
c     
c     
c     need to add dust scattering contribution, 
c     
c     
         if (grainmode) then
c     
c     not so good, need to assume that old tphot is up to date
c     
	  dwsc=0.d0
	  upsc=0.d0
	 if (ClinPAH.AND.(.NOT.pahmode)) then	!Cgrains not around  
 	   do dtype=2,numtypes
            dwsc = dwsc + tphot(inl)*dr*fi*dh*dfscacros(inl,dtype)
            upsc = upsc + tphot(inl)*dr*fi*dh*dbscacros(inl,dtype)
           enddo
         else
 	   do dtype=1,numtypes
            dwsc = dwsc + tphot(inl)*dr*fi*dh*dfscacros(inl,dtype)
            upsc = upsc + tphot(inl)*dr*fi*dh*dbscacros(inl,dtype)
           enddo         
         endif
c
c     Add PAH scattering
c           
          if (pahmode) then
           pahsum=((pahZ(4)+pahZ(5))*pahisca(inl)
     &           +(1.d0-(pahZ(4)+pahZ(5)))*pahnsca(inl))
     &           *0.5d0*(1+pahncos(inl))
           dwsc = dwsc + tphot(inl)*dr*fi*dh*pahsum*pahfrac
           pahsum=((pahZ(4)+pahZ(5))*pahisca(inl)
     &              +(1.d0-(pahZ(4)+pahZ(5)))*pahnsca(inl))
     &              *0.5d0*(1-pahncos(inl))
           upsc = upsc + tphot(inl)*dr*fi*dh*pahsum*pahfrac
          endif

c     
c     Add dust scattered source light to diffuse fields
c     
            dwdif(inl) = dwdif(inl)+(dwsc*(wedw*(dexp(-(upex*tau)))))
            updif(inl) = updif(inl)+(upsc*(weup*(dexp(-(dwex*tau))))) 
c     
         endif
c     
c     do emilin resonance lines
c     
         do m = 1,xlines
c
            if (xbin(m).eq.inl) then
               lin = m
               dismul = (dwex*dwlin(2,lin))+(upex*uplin(2,lin))
               sigmul = sigmt*dismul
               tau = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                    tau = (dr*fi)*sigmul
               if (tau.gt.1.d-4) then 
                  emis = (emilin(1,lin)/sigmul)*(1.d0-dexp(-tau))
                  dwliem = (dwf*emis)*energ
                  upliem = (upf*emis)*energ
               else
                  dwliem = ((dwf*emilin(1,lin))*(fi*dr))*energ
                  upliem = ((upf*emilin(1,lin))*(fi*dr))*energ
               endif
               dismul = (upex*dwlin(2,lin))+(dwex*uplin(2,lin))
               tau = 1.d10
               if ((tau/dismul).gt.tauso) tau = tauso*dismul
               dwlin(1,lin) = dwlin(1,lin)
     &                        +(dwliem*wedw*dexp(-(upex*tau)))
               uplin(1,lin) = uplin(1,lin)
     &                        +(upliem*weup*dexp(-(dwex*tau)))
c     
c     end xbin = inl
c     
            endif
c     
c     end xlines do loop
c     
         enddo
c
c Hydrogen
c
         do series = 1,6
         do line = 1,10
c
            if (hbin(line,series) .eq. inl) then
               dismul = (dwex*hyddwlin(2,line,series))
     &                  +(upex*hyduplin(2,line,series))
               sigmul = sigmt*dismul
               tau = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                  tau = (dr*fi)*sigmul
               if (tau.gt.1.d-4) then 
                  emis = (hydlin(1,line,series)/sigmul)
     &                  *(1.d0-dexp(- tau))
                  dwliem = (dwf*emis)*energ
                  upliem = (upf*emis)*energ
               else
                  dwliem = ((dwf*hydlin(1,line,series))*(fi*dr))*energ
                  upliem = ((upf*hydlin(1,line,series))*(fi*dr))*energ
               endif
               dismul = (upex*hyddwlin(2,line,series))
     &                  +(dwex*hyduplin(2,line,series))
               tau = 1.d10
               if ((tau/dismul).gt.tauso) tau = tauso*dismul
c               
               hyddwlin(1,line,series) = hyddwlin(1,line,series)
     &                  +(dwliem*wedw*dexp(-(upex*tau)))
               hyduplin(1,line,series) = hyduplin(1,line,series)
     &                  +(upliem*weup*dexp(-(dwex*tau)))
c
            endif
c
         enddo
         enddo
c
c Helium
c
         do series = 1,6
         do line = 1,10
c
            if (hebin(line,series) .eq. inl) then
               dismul = (dwex*heldwlin(2,line,series))
     &                  +(upex*heluplin(2,line,series))
               sigmul = sigmt*dismul
               tau = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                  tau = (dr*fi)*sigmul
               if (tau.gt.1.d-4) then 
                  emis = (hellin(1,line,series)/sigmul)
     &                  *(1.d0-dexp(- tau))
                  dwliem = (dwf*emis)*energ
                  upliem = (upf*emis)*energ
               else
                  dwliem = ((dwf*hellin(1,line,series))*(fi*dr))*energ
                  upliem = ((upf*hellin(1,line,series))*(fi*dr))*energ
               endif
               dismul = (upex*heldwlin(2,line,series))
     &                  +(dwex*heluplin(2,line,series))
               tau = 1.d10
               if ((tau/dismul).gt.tauso) tau = tauso*dismul
c
               dfb  = 1.d0
               if ((line.eq.1).and.(series.eq.1)) then
                 dfb  = 1.d0-fb
               endif
c               
             heldwlin(1,line,series) = heldwlin(1,line,series)
     &                          +(dfb*dwliem*wedw*dexp(-(upex*tau)))
c
             heluplin(1,line,series) = heluplin(1,line,series)
     &                          +(dfb*upliem*weup*dexp(-(dwex*tau)))
c
            endif
c
         enddo
         enddo
c
c
c Heavy Hydrogenic
c
         do atom = 3,atypes
         if (xhydlin(1,1,1,atom).gt.epsilon) then
         do series = 1,2
         do line = 1,10
c
            if (xhbin(line,series,atom) .eq. inl) then
               dismul = (dwex*xhyddwlin(2,line,series,atom))
     &                  +(upex*xhyduplin(2,line,series,atom))
               sigmul = sigmt*dismul
               tau = 1.d10
               if ((tau/dismul).gt.((dr*fi)*sigmt)) 
     &                  tau = (dr*fi)*sigmul
               if (tau.gt.1.d-4) then 
                  emis = (xhydlin(1,line,series,atom)/sigmul)
     &                  *(1.d0-dexp(- tau))
                  dwliem = (dwf*emis)*energ
                  upliem = (upf*emis)*energ
               else
                  dwliem = ((dwf*xhydlin(1,line,series,atom))
     &                  *(fi*dr))*energ
                  upliem = ((upf*xhydlin(1,line,series,atom))
     &                  *(fi*dr))*energ
               endif
               dismul = (upex*xhyddwlin(2,line,series,atom))
     &                  +(dwex*xhyduplin(2,line,series,atom))
               tau = 1.d10
               if ((tau/dismul).gt.tauso) tau = tauso*dismul
c               
               xhyddwlin(1,line,series,atom) = 
     &                  xhyddwlin(1,line,series,atom)
     &                  +(dwliem*wedw*dexp(-(upex*tau)))
               xhyduplin(1,line,series,atom) = 
     &                  xhyduplin(1,line,series,atom)
     &                  +(upliem*weup*dexp(-(dwex*tau)))
c
            endif
c
         enddo
         enddo
         endif
         enddo
c
c end vector loop inl
c     
 500  continue
c     
      return 
      end









