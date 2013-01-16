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
c     BAG09 : JMOD='DWUP' or 'LOUP" are never used?!?
c     
      subroutine newdif (tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup
     &                   ,frdw, jmod)
c     
      include 'cblocks.inc'
c
c           Input Variables
c
      double precision tdw, tup, dh, fi, rad
      double precision drdw, dvdw, drup, dvup, frdw
      character jmod*4
c
c           Local Variables 
c
c      double precision at,bet,eph,se
      double precision den,dildw
      double precision dilup,dismul,dr,dvem,dvemh,dwem
      double precision dwex, dwf, dwliem, energ,telc
      double precision upem, upex, upf, upliem, wadw, waup, wedw, weup
      double precision tauso, tau, tau0, sigmt
      double precision curad, rm, wlo, ulo
      double precision emis,z,p,f,es
      double precision dfbu,dfbd,dfb
      double precision fbu,fbd,fb
      double precision dustsigmat
      double precision pathlength
      double precision escape
c
c      double precision drdwh,dvdwh,druph,dvuph  !upwind never used?
c
c          External Functions
c
      double precision fdismul,localout, tauline, fbowen
c
      integer*4 line, series, j, atom, ion
      integer*4 inl, lin, m, n
c
c      
c     Check parameters and replace with sane values if necessary      
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
 
 
        tdw  = dabs(tdw)
        tup  = dabs(tup)
        rad  = dabs(rad)
        drdw = dabs(drdw)
        drup = dabs(drup)
        dvdw = dabs(dvdw)
        dvup = dabs(dvup)
        
      endif
c     
      if ((jmod .ne. 'LODW')
     &     .and.(jmod .ne. 'LOUP')
     &     .and.(jmod .ne. 'OUTW')
     &     .and.(jmod .ne. 'DWUP')) then
     
        write(*, 67) jmod
 67     format(//' MODE IMPROPERLY DEFINED IN NEWDIF : ',a4)
        stop 
      
      endif
c     
      if ((frdw.gt.1.d0).or.(frdw.lt.0.d0)) then
        write(*, 77) frdw
 77     format(/,/,' INCONSISTENT VALUE FOR THE ARGUMENT FRDW :'
     &         ,1pg10.3)
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
        dvem = dvup
        telc = tup
      else
        dr = drdw
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
 123      format(/,/,' DRUP IS LARGER THAN THE RADIUS OF CURVATURE :'
     &            ,2(1pg10.3))
          stop
        endif  
        ulo = (2.d0*dlog(curad+dr))-(2.d0*dlog(rm))
        dildw = dexp(wlo)
        dilup = dexp(ulo)
      endif
c     
      if (jmod.eq.'LOUP') then !upwind local
        wadw = 1.d0
        wedw = dildw
        waup = dilup
        weup = 1.d0
      else if (jmod.eq.'DWUP') then !downwind & upwind local
        wadw = dildw
        wedw = 1.d0
        waup = dilup
        weup = 1.d0
      else !downwind local (STANDARD)
        wadw = dildw
        wedw = 1.d0
        waup = 1.d0
        weup = dilup
      endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         Begin Routine Proper
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ***COMPUTES RESONANCE-LINES DISTANCE MULTIPLIERS
c
c     NOTE: single resonance line multipliers never used 
c           all split up into dw and up     
c
      dvemh = dvem*0.5d0
c
c     Downstream
c
      if (dwf.gt.0.d0) then
c
c first the general resonance lines     
c
c         drdwh=drdw*0.5d0  ! old diffuse crosses full zone ie dr
c         dvdwh=dvdw*0.5d0  ! so replaced d?dwh with d?dw BAG09
         do line = 1,xlines
           z = zion(xrat(line))
           p = pop(xion(line),xrat(line))
           dwlin(2,line) = 1.d0

           if (xbr(line).gt.0.0) then
             f = xfef(line)/xbr(line)
             es = xejk(line)*ev
             if ((z*p.ge.pzlimit))then
               dwlin(2,line) = fdismul(tdw, dh, drdw, dvdw
     &                         , xrat(line), xion(line), es,f)
             endif
           endif
         enddo
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
                  hyddwlin(2,line,series) = fdismul(tdw, dh, drdw
     &                                 , dvdw, atom, ion, es,f)
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
                  heldwlin(2,line,series) = fdismul(tdw, dh, drdw
     &                                     , dvdw, atom, ion, es,f)
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
     &                                , drdw, dvdw, atom, ion, es,f)
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
c     Upstream
c  
      if (upf.gt.0.d0) then
c
c first the general resonance lines     
c
         do line = 1,xlines   
           z = zion(xrat(line))
           p = pop(xion(line),xrat(line))
           uplin(2,line) = 1.d0
           if (xbr(line).gt.0.0) then
             f = xfef(line)/xbr(line)
             es = xejk(line)*ev
             if ((z*p.ge.pzlimit))then
               uplin(2,line) = fdismul(tup, dh, drup, dvup, xrat(line)
     &                         , xion(line), es,f)
             endif
           endif
         enddo
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
              hyduplin(2,line,series) = 1.d0
              j = hbin(line,series)
              if (j.ne.0) then
                es = (LmeV/hlambda(line,series))*ev
                if ((z*p.ge.pzlimit).and.(series.eq.1))then
                  hyduplin(2,line,series) = fdismul(tup, dh, drup, dvup
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
                 heluplin(2,line,series) = fdismul(tup, dh, drup, dvup
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
                 xhyduplin(2,line,series,atom) = fdismul(tup, dh, drup
     &                                         , dvup, atom, ion, es,f)
                endif
              endif
            endif
c
           enddo
          enddo
         enddo
c
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***UPDATES DIFFUSE FIELD VECTORS UPDIF, DWDIF etc BIN BY BIN
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Bowen fluorescence fractions for He II 303A
c
      fb  = 0.0d0
      fbu = 0.0d0
      fbd = 0.0d0
      if (zmap(2).gt.0)then 
        fb  = fbowen(telc,dvemh)
        fbu = fbowen(tup,dvup)
        fbd = fbowen(tdw,dvdw)
      endif
c
      pathlength = dr*fi        
      do 500 inl = 1, infph-1
c     
         den        = 0.5d0*(ephot(inl)+ephot(inl+1))
         energ      = den*ev
         tauso      = 0.d0
         sigmt      = 0.d0
         dustsigmat = 0.d0
            
         call crosssections(inl, tauso, sigmt, dustsigmat)
            
         sigmt      = dh*sigmt
         dustsigmat = dh*dustsigmat

c     
c     ***ATTENUATION OF THE ORIGINAL INTENSITY OF THE DIFFUSE FIELD
c     
         tau0      = pathlength*sigmt
         tau       = tau0
c        
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
               tau        =   tauline(dismul,tau0)
c
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
               tau        =   tauline(dismul,tau0)
c
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
            tau        =   tauline(dismul,tau0)
c            
            dfbu = 1.d0
            dfbd = 1.d0
            if ((line.eq.1).and.(series.eq.1)) then
              dfbu = 1.d0-fbu
              dfbd = 1.d0-fbd
            endif
c  
            heldwlin(1,line,series) = heldwlin(1,line,series)
     &                                *(wadw*dexp(-(dwex*tau)))
            heluplin(1,line,series) = heluplin(1,line,series)
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
               tau        =   tauline(dismul,tau0)
c  
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
         tau    = tau0
         escape = localout( 1.d0, tau)         
c
c Factor of x2 for PP geometry fixed RSS 2006
c
         emis   =  energ*emidif(inl)*pathlength
c
c old newdif - RSS 06
c
c         dwem = (dwf*emis)
c         upem = (upf*emis)
c     Ne newdif allows for escape-- consistent with totphot BAG 06
c       
         dwem = (dwf*emis)*escape
         upem = (upf*emis)*escape
c


c     
c     Add local diffuse field through this zone, assuming 
c
c     original had additional extinction... not sure why
c      as already calculated escape fraction.
c         dwdif(inl) = dwdif(inl)+(dwem*(wedw*dexp(-(upex*tau))))
c         updif(inl) = updif(inl)+(upem*(weup*dexp(-(dwex*tau))))
c
c     
c
         dwdif(inl) = dwdif(inl)+(dwem*wedw)   
         updif(inl) = updif(inl)+(upem*weup)   

c     
c     Scattering presumed to be 100% forward
c     set to Zero till I fix it BAG 05/2006
c         if (grainmode) then
c     
c     not so good, need to assume that old tphot is up to date
c     
c	  dwsc=0.d0
c	  upsc=0.d0
c	 if (ClinPAH.AND.(.NOT.pahmode)) then	!Cgrains not around  
c 	   do dtype=2,numtypes
c            dwsc = dwsc + tphot(inl)*dr*fi*dh*dfscacros(inl,dtype)
c            upsc = upsc + tphot(inl)*dr*fi*dh*dbscacros(inl,dtype)
c           enddo
c         else
c 	   do dtype=1,numtypes
c            dwsc = dwsc + tphot(inl)*dr*fi*dh*dfscacros(inl,dtype)
c            upsc = upsc + tphot(inl)*dr*fi*dh*dbscacros(inl,dtype)
c           enddo         
c         endif
c
c     Add PAH scattering
c           
c          if (pahmode) then
c           pahsum=((pahZ(4)+pahZ(5))*pahisca(inl)
c     &           +(1.d0-(pahZ(4)+pahZ(5)))*pahnsca(inl))
c     &           *0.5d0*(1+pahncos(inl))
c           dwsc = dwsc + tphot(inl)*dr*fi*dh*pahsum*pahfrac
c           pahsum=((pahZ(4)+pahZ(5))*pahisca(inl)
c     &              +(1.d0-(pahZ(4)+pahZ(5)))*pahnsca(inl))
c     &              *0.5d0*(1-pahncos(inl))
c           upsc = upsc + tphot(inl)*dr*fi*dh*pahsum*pahfrac
c          endif
c
c     
c     Add dust scattered source light to diffuse fields
c     
c            dwdif(inl) = dwdif(inl)+(dwsc*(wedw*(1.d0)))
c            updif(inl) = updif(inl)+(upsc*(weup*(0.d0))) 
c     
c         endif
c
c     
c     do emilin resonance lines
c     
         do m = 1,xlines
            if (xbin(m).eq.inl) then
            
               lin = m

               dismul = (dwex*dwlin(2,lin))+(upex*uplin(2,lin))
               
               tau        =              tau0
               escape     = localout( dismul, tau)

               energ      =      ev*xejk(lin)               
               emis       = energ*emilin(1,lin)*pathlength
               
               dwliem = (dwf*emis)*escape
               upliem = (upf*emis)*escape
c
c
               dwlin(1,lin) = dwlin(1,lin) +(dwliem*wedw)
               uplin(1,lin) = uplin(1,lin) +(upliem*weup)
c
c I don't understand this bit: How can the local distance multiplier 
c be used for the integrated opacity back to the source?
c
c And what is upex*dwlin(2,lin) mean? (or dwex*uplin(2,lin))
c
c Also, here we have the down stream local escape, being split *again*
c by upex
c                  dwliem*wedw*dexp(-(upex*tau))
c 
c  After dwliem is already split downstream:
c
c                  dwliem = (dwf*emis)*energ
c
c
c  for photo outward models:
c
c                 dwf = 1.0  upf = 0.0
c                 wedw = 1.0 weup = 1.0
c                 dwex = 1.0 upex = 0.0
c
c  so dexp(-(upex*tau)) = 1.0
c  and dexp(-(dwex*tau)) = some small number, but upliem = 0.0
c
c   (upex*dwlin(2,lin)) = 0.0
c   (dwex*uplin(2,lin)) = upstream distance multiplier
c 
ccccccc
c
c               dismul = (upex*dwlin(2,lin))+(dwex*uplin(2,lin))
c               tau = 1.d10
c               if ((tau/dismul).gt.tauso) tau = tauso*dismul
c               dwlin(1,lin) = dwlin(1,lin)
c     &                        +(dwliem*wedw*dexp(-(upex*tau)))
c               uplin(1,lin) = uplin(1,lin)
c     &                        +(upliem*weup*dexp(-(dwex*tau)))
c
ccccccc
c
c
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
     
               tau        =              tau0
               escape   = localout( dismul, tau)

               energ      = ev*LmeV/hlambda(line,series)
               emis      = energ*hydlin(1,line,series)*pathlength
               
               dwliem = (dwf*emis)*escape
               upliem = (upf*emis)*escape
c
c               
               hyddwlin(1,line,series) = hyddwlin(1,line,series)
     &                  +(dwliem*wedw)
               hyduplin(1,line,series) = hyduplin(1,line,series)
     &                  +(upliem*weup)
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
     
               tau        =              tau0
               escape   = localout( dismul, tau)

               energ      = ev*LmeV/helambda(line,series)
               emis       = energ*hellin(1,line,series)*pathlength
               
               dwliem = (dwf*emis)*escape
               upliem = (upf*emis)*escape
c
               dfb  = 1.d0
               if ((line.eq.1).and.(series.eq.1)) then
                 dfb  = 1.d0-fb
               endif
c             
             heldwlin(1,line,series) = heldwlin(1,line,series)
     &                          +(dwliem*wedw)
             heluplin(1,line,series) = heluplin(1,line,series)
     &                          +(upliem*weup)
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
       
               tau        =              tau0
               escape   = localout( dismul, tau)

               energ      = ev*LmeV/xhlambda(line,series,atom)
               emis       = energ*xhydlin(1,line,series,atom)*pathlength
               
               dwliem = (dwf*emis)*escape
               upliem = (upf*emis)*escape
c
               xhyddwlin(1,line,series,atom) = 
     &                  xhyddwlin(1,line,series,atom)
     &                  +(dwliem*wedw)
               xhyduplin(1,line,series,atom) = 
     &                  xhyduplin(1,line,series,atom)
     &                  +(upliem*weup)
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









