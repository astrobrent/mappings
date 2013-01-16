cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    
c     Set grain parameters
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine grainpar()
c     
c     
      include 'cblocks.inc'
c     
c     
      double precision volume,pow,pow1,const,fw
      double precision dgr,ghr,siggrain(mxdtype),dgrad,dustmax
      double precision Nc
      integer*4 i,j,k,inl,dtype
      character ilgg*4
c     
      if (runmode.ne.'batchrun') then
 100     format(/' Include dust? (y/n) :',$)
         write(*,100)
      endif
      read(*,'(a)') ilgg
      ilgg = ilgg(1:1)
      if (ilgg.eq.'y') ilgg = 'Y'
c     
c     defaults
c     
c     Standard Grain parameters
c     
c     
c     Bgrain: threshold parameter
c     Yinf: yield parameter
c     alh: accomodation parameter
c     segrain,spgrain: sticking for electrons and protons
c     siggrain: projected area of grains per H atom
c     graindens : grain  density in g/cm^-3
c     
c     
      siggrain(1) = 1.5e-21
      Bgrain = 8.0d0
      segrain = 0.5
      spgrain = 1.d0
      gheat = 0.d0
      gcool = 0.d0
      Yinf = 0.05d0
      haccom = 0.14d0
      graindens(1) = 2.0d0
      amin(1) = 0.005
      amax(1) = 0.25
      galpha = -3.5d0
      Nc = 468
      mindust(1)=1
      maxdust(1)=dustbinmax
      mindust(2)=1
      maxdust(2)=dustbinmax
c     
      grainmode = .FALSE.
      pahmode = .FALSE.
      pahCfrac=0.d0
      pahfrac=0.d0
      IRmode = 0
      IRtemp = .False.
      do dtype=1,numtypes
      do inl=1,infph-1
        dcrosec(inl,dtype) = 0.d0
        dfscacros(inl,dtype) = 0.d0
        dbscacros(inl,dtype) = 0.d0
      enddo
      enddo

      if (ilgg.eq.'Y') then
         grainmode = .TRUE.
         call depcha
c     
c     
c     Initialize anyway for safety
         numtypes=2
c
c     Detemine grain distribution
c     
 110     write(*,111) 
 111     format(/' Enter the grain distribution model '/
     &        ' ::::::::::::::::::::::::::::::::::::::'/
     &        '    P   :  Powerlaw     N(a) = k a^alpha'/
     &        '    M   :  MRN distribution'/
     &	      '    S   :  Grain Shattering Profile'/
     &        '   (Note: need full DUSTDATA for P or S)'/
c     &        '    E   :  Exponential  N(a) = k exp(a/alpha)'/
     &        /' :: ',$)
         read(*,'(a)') gdist
         gdist = gdist(1:1)
         if (gdist.eq.'p') gdist = 'P'
         if (gdist.eq.'m') gdist = 'M'
         if (gdist.eq.'s') gdist = 'S'
c
         if (gdist.eq.'M') then
c
c     MRN dist N(a) = k a^-3.5 C amin=50A Sil amin=100A,amax=2500 A
c
 115      galpha = -3.5
          amin(1) = 0.0050
          amin(2) = 0.0100
          amax(1) = 0.2500
          amax(2) = 0.2500
c
c  convert to cm
c
           amin(1) = amin(1)*1d-4
           amax(1) = amax(1)*1d-4
           amin(2) = amin(2)*1d-4
           amax(2) = amax(2)*1d-4

         else if (gdist.eq.'P') then
c
c	Power - Law N(a) = A^alpha
c                 
 120      write(*,121)
 121      format(/,'Enter alpha  (MRN=-3.5)::',$)
          read(*,*) galpha

c
c     Enter min and max grain radii (within 1e-3, 10 mu)
c     then find corresponding grain radius bins
c     
c     Graphite first (type=1)
c
 130      if (runmode.ne.'batchrun') then
 131         format(/' Give graphite grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
             write(*,131)
          endif
          read(*,*) amin(1)
c     
          if (amin(1).le.0.d0) goto 130
          if (amin(1).ge.10.d0) amin(1) = amin(1)/1d4
          if ((amin(1).lt.1.0d-3).OR.(amin(1).gt.10.d0)) goto 130
c     
 140      if (runmode.ne.'batchrun') then
 141         format(/' Give graphite grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
             write(*,141)
          endif
          read(*,*) amax(1)
c     
          if (amax(1).le.0.d0) goto 140
          if (amax(1).ge.10.d0) amax(1) = amax(1)/1d4
          if ((amax(1).lt.1.0d-3).OR.(amax(1).gt.10.d0)) goto 140
          if (amax(1).lt.amin(1)) goto 130

c     
c     Then Silicate (type=2)
c
 150      if (runmode.ne.'batchrun') then
 151         format(/' Give silicate grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
             write(*,151)
          endif
          read(*,*) amin(2)
c     
          if (amin(2).le.0.d0) goto 150
          if (amin(2).ge.10.d0) amin(2) = amin(2)/1d4
          if ((amin(2).lt.1.0d-3).OR.(amin(2).gt.10.d0)) goto 150
c     
 160      if (runmode.ne.'batchrun') then
 161         format(/' Give silicate grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
            write(*,161)
          endif
          read(*,*) amax(2)
c     
          if (amax(2).le.0.d0) goto 160
          if (amax(2).ge.10.d0) amax(2) = amax(2)/1d4
          if ((amax(2).lt.1.0d-3).OR.(amax(2).gt.10.d0)) goto 160
          if (amax(2).lt.amin(2)) goto 150
c
c  convert to cm
c
          amin(1) = amin(1)*1d-4
          amax(1) = amax(1)*1d-4
          amin(2) = amin(2)*1d-4
          amax(2) = amax(2)*1d-4

         else if (gdist.eq.'S') then
c
c	Grain shattering distribution based on Jones et al 1996
c	 with additional exponential cutoffs
c	N(a) = k * a^-3.3 * (exp[-(a/amin)^-3])/(exp[(a/amax)^3])
c
310	format(/,' Shattered grain distribution:',/
     &   ,'N(a) = k * a^alpha * (exp[-(a/amin)^-3])/(exp[(a/amax)])',/)
	write(*,310) 
c     Enter min and max grain radii (within 1e-3, 10 mu)
c     then find corresponding grain radius bins
c     
 320      write(*,321)
 321      format(/,'Enter alpha  (Jones=-3.3, MRN=-3.5)::',$)
          read(*,*) galpha

c     Graphite first (type=1)
c
 330      if (runmode.ne.'batchrun') then
 331         format(/' Give graphite grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
             write(*,331)
          endif
          read(*,*) amin(1)
c     
          if (amin(1).le.0.d0) goto 330
          if (amin(1).ge.10.d0) amin(1) = amin(1)/1d4
          if ((amin(1).lt.1.0d-3).OR.(amin(1).gt.10.d0)) goto 330
c     
 340      if (runmode.ne.'batchrun') then
 341         format(/' Give graphite grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
             write(*,341)
          endif
          read(*,*) amax(1)
c     
          if (amax(1).le.0.d0) goto 340
          if (amax(1).ge.10.d0) amax(1) = amax(1)/1d4
          if ((amax(1).lt.1.0d-3).OR.(amax(1).gt.10.d0)) goto 340
          if (amax(1).lt.amin(1)) goto 330

c     
c     Then Silicate (type=2)
c
 350      if (runmode.ne.'batchrun') then
 351         format(/' Give silicate grain min. radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
             write(*,351)
          endif
          read(*,*) amin(2)
c     
          if (amin(2).le.0.d0) goto 350
          if (amin(2).ge.10.d0) amin(2) = amin(2)/1d4
          if ((amin(2).lt.1.0d-3).OR.(amin(2).gt.10.d0)) goto 350
c     
 360      if (runmode.ne.'batchrun') then
 361         format(/' Give silicate grain max. (> grain min) radius ',
     &           '(<10 in mu, >=10 in Angs.) :',$)
            write(*,361)
          endif
          read(*,*) amax(2)
c     
          if (amax(2).le.0.d0) goto 360
          if (amax(2).ge.10.d0) amax(2) = amax(2)/1d4
          if ((amax(2).lt.1.0d-3).OR.(amax(2).gt.10.d0)) goto 360
          if (amax(2).lt.amin(2)) goto 350
c
c  convert to cm
c
          amin(1) = amin(1)*1d-4
          amax(1) = amax(1)*1d-4
          amin(2) = amin(2)*1d-4
          amax(2) = amax(2)*1d-4
         else 
c
c  Wrong Grain size distribution
c         
         
            goto 110
         endif

c
c  get edges
c
	if ((gdist.eq.'M').OR.(gdist.eq.'P')) then
         mindust(1)=1
         maxdust(1)=dustbinmax
         do i=1,dustbinmax
           if((amin(1).ge.gradedge(i)).AND.(amin(1).lt.gradedge(i+1))) 
     &      then
             mindust(1)=i
             amin(1) = grainrad(i)
           endif
           if((amax(1).gt.gradedge(i)).AND.(amax(1).le.gradedge(i+1)))
     &      then
             maxdust(1)= i
             amax(1) = grainrad(i)
           endif
         enddo
           mindust(2)=1
           maxdust(2)=dustbinmax
         do i=1,dustbinmax
           if((amin(2).ge.gradedge(i)).AND.(amin(2).lt.gradedge(i+1)))
     &      then
             mindust(2)=i
             amin(2) = grainrad(i)
           endif
           if((amax(2).gt.gradedge(i)).AND.(amax(2).le.gradedge(i+1)))
     &      then
             maxdust(2) = i
             amax(2) = grainrad(i)
           endif
         enddo
	endif
         write(*,135) amin(1)*1d4
         write(*,145) amax(1)*1d4
 135     format('***************************',/,
     &      'Graphite grain min radius (mu): ',e11.4)
 145     format('Graphite grain max radius (mu): ',e11.4,/,
     &      '***************************')
         write(*,155) amin(2)*1d4
         write(*,165) amax(2)*1d4
 155     format('***************************',/,
     &      'Silicate grain min radius (mu): ',e11.4)
 165     format('Silicate grain max radius mu): ',e11.4,/,
     &      '***************************')

c
c	Photoelectric parameters
c
            Bgrain       = 8.0
            Yinf         = 0.05d0
c            Yinf         = 0.1


         if ((gdist.eq.'M')) then
            graindens(1) = 1.85
            graindens(2) = 2.5
            goto 219
          endif

 170     write(*,171)
 171     format(/' Give graphite grain density (~1.85 g/cm^3) :',$)

         read(*,*) graindens(1)     
         if (graindens(1).le.0.d0) goto 170

 180     write(*,181)
 181     format(/' Give silicate grain density (~2.5 g/cm^3) :',$)
         read(*,*) graindens(2)     
         if (graindens(2).le.0.d0) goto 180
c
c     Brent: Still need to be changed to distribution
c     
c 190     if (runmode.ne.'batchrun') then
c 200        format(/' Give grain photoelectric ',
c     &           'threshold(~8 eV [5-13]) :',$)
c            write(*,200)
c         endif
c         read(*,*) Bgrain
c     
c         if (Bgrain.lt.5.d0) goto 190
c         if (Bgrain.gt.13.d0) goto 190
c     
c 210     if (runmode.ne.'batchrun') then
c 211        format(/' Give grain photoelectric ',
c     &           'yield (~0.1 at inf.) :',$)
c            write(*,211)
c         endif
c         read(*,*) Yinf
c     
c         if (Yinf.lt.0.d0) goto 210
c         if (Yinf.gt.1.d0) goto 210
c     
c
c     
c     get PAHs at this point so carbon depletion is consistent
c     
         ClinPAH=.False.
c
 219     if (runmode.ne.'batchrun') then
 220        format(/' Include PAH molecules? (y/n) :',$)
            write(*,220)
         endif
         read(*,'(a)') ilgg
         ilgg = ilgg(1:1)
         if (ilgg.eq.'y') ilgg = 'Y'
C          ilgg='N'
c     
        if (ilgg.eq.'Y') then
         Nc=468
c            
 230      if (runmode.ne.'batchrun') then
 231       format(/' Give fraction of Carbon Dust',
     &            ' Depletion in PAHs :',$)
           write(*,231)
          endif
          read(*,*) pahCfrac
c     
          if (pahCfrac.lt.0.d0) goto 230
          if (pahCfrac.gt.1.d0) goto 230
c
 234      write(*,235)
 235      format(/,' Choose PAH switch on;',/
     &     ,' :::::::::::::::::::::::::::::::::::::::::::::::::::',/
     &     ,'   H :  Habing photodissociation parameter < Value',/
     &     ,'   I :  PAH ionizing radiation/ local ISRF ratio < Value',/
     &     ,'   Q :  QHDH < Value',/
     &     ,' :: ',$)
	  read(*,*) ilgg
	  pahend=ilgg(1:1)
	  write(*,*) ilgg, pahend
	  if (pahend.eq.'h') pahend='H'
	  if (pahend.eq.'i') pahend='I'
	  if (pahend.eq.'q') pahend='Q'
	  if (.NOT.((pahend.eq.'H')
     &        .OR.(pahend.eq.'I')
     &        .OR.(pahend.eq.'Q'))) goto 234
 440	  write(*,441)
 441	  format(/'Give PAH switch on Value :',$)
 	  read(*,*) PAHlimit
 	  if (PAHlimit.lt.0.d0) goto 440
	  	  
c    
c   Convert to #PAH grains per H atom

          pahfrac=zion(zmap(6))*(1.d0-dion(zmap(6)))*pahCfrac/Nc
          
          write(*,'("PAH per H (ppm): ",1pg11.4)') pahfrac*1.d6

          write(*,442) 
 442      format('Do you wish graphite grains to be cospatial'
     &           ,' with PAHs?',/
     &           ,'(ie PAH destroyed = C grain destroyed):',$)
          read(*,*) ilgg
          ilgg = ilgg(1:1)
          if (ilgg.eq.'y') ilgg = 'Y'
	  if (ilgg.eq.'Y') then
	    ClinPAH =.TRUE.
	  endif         
        endif
 260  write(*,261)
 261  format(/' Evaluate dust temperatures and IR flux? (P5 only)',/,
     &  '(WARNING: will slow down computing time) (y/N):',$)
      read(*,'(a)') ilgg
      ilgg = ilgg(1:1)
      if ((ilgg.eq.'y').OR.(ilgg.eq.'Y')) then
 265     write(*,266)
 266     format(/,' Which IR model? ',/
     &      ,' ::::::::::::::::::::::::::::::::::::::',/
     &      ,'    Q :  Quick (50 bins, ~10% accurate IR)',/
     &      ,'    I :  Intermediate (100 bins, ~2% accurate',/
     &      ,'    S :  Slow (400 bins, accurate Temp. distributions)',/
     &      ,/,' :: ',$)
         read(*,'(a)') ilgg
         ilgg = ilgg(1:1)
         if ((ilgg.eq.'q').OR.(ilgg.eq.'Q')) then
           IRmode=1
         else if ((ilgg.eq.'i').OR.(ilgg.eq.'I')) then
           IRmode=2
         else if ((ilgg.eq.'s').OR.(ilgg.eq.'S')) then
           IRmode=3
         else
           goto 265
         endif
         
c       if (expertmode.gt.0) then
         write(*,270)
  270    format(/,'Output dusttemperature data (large file) (y/n):',$) 
         read(*,'(a)') ilgg
         ilgg = ilgg(1:1)
         if ((ilgg.eq.'y').OR.(ilgg.eq.'Y')) then
           IRtemp=.True.
         else
           IRtemp=.False.
         endif
c       endif
      endif
c
c     Solve for factor k in distribution (num dens wrt H)
c     and determine distribution
	if ((gdist.eq.'M').OR.(gdist.eq.'P')) then
c	Power law distribution (incl MRN)
c
c
c     Graphite grains
c
c     solve for weight of dust
c
         fw = atwei(zmap(6))*zion(zmap(6))*(1.d0-dion(zmap(6)))
     &        *(1.d0-pahCfrac)
c        fw = atwei(zmap(6))*zion(zmap(6))*(1.d0-dion(zmap(6)))
         fw = fw*amu
         dustmass(1)=fw
c
c     Then determine 'volume and solve for k
c
         if (mindust(1).eq.maxdust(1)) then
           volume=ftpi*amax(1)**3.d0
         else
           if(galpha.ne.3.d0) then
             pow = galpha+3.d0
             pow1 = pow+1.d0
             volume = ftpi*(((amax(1)**pow1)/pow1)-
     &                ((amin(1)**pow1)/pow1))
           else
             volume = ftpi*(log(amax(1)) - log(amin(1)))
           endif
         endif
         const = fw/(volume*graindens(1))
c
c     Evaluate graphite grain distribution per H atom
c     (1e4 is diff factor between cubic and square microns in cm)
c
         if (mindust(1).eq.maxdust(1)) then
           dgrad=gradedge(mindust(1)+1)-gradedge(mindust(1))
           dustnum(maxdust(1),1) = const/dgrad
           dustsig(maxdust(1),1) = pi*const*grainrad(maxdust(1))
     &              *grainrad(maxdust(1))/dgrad
         else
          do i=mindust(1),maxdust(1)
            dustnum(i,1) = const * grainrad(i)**galpha
            dustsig(i,1) = pi*const*grainrad(i)**(galpha+2.d0)
          enddo
         endif
c     
c     Silicate grains
c
c     solve for weight of dust
c
         fw = 0.d0
         do i = 1,atypes
            if (mapz(i).ne.6) then
               fw = fw + atwei(i)*zion(i)*(1.d0-dion(i))
            endif
         enddo
         fw = fw*amu
         dustmass(2)=fw
c
c     Then determine volume and solve for k
c
         if (mindust(2).eq.maxdust(2)) then
           volume=ftpi*amax(2)**3.d0
         else
           if(galpha.ne.3.d0) then
             pow = galpha+3.d0
             pow1 = pow+1.d0
             volume = ftpi*(((amax(2)**pow1)/pow1)-
     &                 ((amin(2)**pow1)/pow1))
           else
             volume = ftpi*(log(amax(2)) - log(amin(2)))
           endif
         endif
         const = fw/(volume*graindens(2))
c
c     Evaluate silicate grain distribution per H atom
c
         if (mindust(2).eq.maxdust(2)) then
           dgrad=gradedge(mindust(2)+1)-gradedge(mindust(2))
           dustnum(maxdust(2),2) = const/dgrad
           dustsig(maxdust(2),2) = pi*const*grainrad(maxdust(2))
     &             *grainrad(maxdust(2))/dgrad
         else
           do i=mindust(2),maxdust(2)
             dustnum(i,2) = const * grainrad(i)**galpha 
             dustsig(i,2) = pi*const*grainrad(i)**(galpha+2.d0)
           enddo
         endif
        else if (gdist.eq.'S') then
c
c    Grain shattering distribution
c        
c     Graphite grains
c
c     solve for weight of dust
c
         fw = atwei(zmap(6))*zion(zmap(6))*(1.d0-dion(zmap(6)))
     &        *(1.d0-pahCfrac)
c        fw = atwei(zmap(6))*zion(zmap(6))*(1.d0-dion(zmap(6)))
         fw = fw*amu
c
c     Then determine 'volume and solve for k
c
	 volume = 0.d0
         do i=1,dustbinmax
           dgrad = gradedge(i+1)-gradedge(i)
	   dustnum(i,1) = grainrad(i)**(galpha)*
     &         exp(-(grainrad(i)/amin(1))**(-3) -(grainrad(i)/amax(1)))
     	   volume = volume + dustnum(i,1)*grainrad(i)**3*dgrad
     	 enddo
         const = fw/(ftpi*volume*graindens(1))
c
c     Evaluate graphite grain distribution per H atom
c     (1e4 is diff factor between cubic and square microns in cm)
c
	 dustmax=0.d0
         do i=1,dustbinmax
           dustnum(i,1) = const*dustnum(i,1)
           dustsig(i,1) = pi*grainrad(i)**2*dustnum(i,1)
           if (dustnum(i,1).gt.dustmax) dustmax = dustnum(i,1)
         enddo
         i=1
	 do while(dustnum(i,1).lt.1.0d-10*dustmax)
	  i=i+1
	 enddo
	 mindust(1)=i
	 i=i+1
	 do while((dustnum(i,1).gt.1.0d-10*dustmax).AND.(i.lt.dustbinmax))
	  i=i+1
	 enddo	
	  maxdust(1) = i
c     
c     Silicate grains
c
c     solve for weight of dust
c
         fw = 0.d0
         do i = 1,atypes
            if (mapz(i).ne.6) then
               fw = fw + atwei(i)*zion(i)*(1.d0-dion(i))
            endif
         enddo
         fw = fw*amu
c
c     Then determine volume and solve for k
c
	 volume = 0.d0
         do i=1,dustbinmax
           dgrad = gradedge(i+1)-gradedge(i)
	   dustnum(i,2) = grainrad(i)**(galpha)*
     &        exp(-(grainrad(i)/amin(2))**(-3) 
     &            -(grainrad(i)/amax(2))**(3))
     	   volume = volume + dustnum(i,2)*grainrad(i)**3*dgrad
     	 enddo
         const = fw/(ftpi*volume*graindens(2))
c
c     Evaluate silicate grain distribution per H atom
c
	 dustmax=0.d0
         do i=1,dustbinmax
           dustnum(i,2) = const*dustnum(i,2)
           dustsig(i,2) = pi*grainrad(i)**2*dustnum(i,2)
           if (dustnum(i,2).gt.dustmax) dustmax = dustnum(i,2)
         enddo
         i=1
	 do while(dustnum(i,2).lt.1.0d-10*dustmax)
	  i=i+1
	 enddo
	 mindust(2)=i
	 i=i+1
	 do while((dustnum(i,2).gt.1.0d-10*dustmax).AND.(i.lt.dustbinmax))
	  i=i+1
	 enddo	
	 maxdust(2) = i
   
        endif
     
c     
c     calc ratios
c
         dgr = 0.d0
c     
         do i = 1,atypes
            if (mapz(i).ne.6) then
            dgr = dgr + atwei(i)*zion(i)*(1.d0-dion(i))
            endif
            if (mapz(i).eq.6) then
            dgr = dgr+atwei(i)*zion(i)*(1.d0-dion(i))*(1.d0-pahCfrac)
            endif
         enddo
c     
         ghr = 0.d0
c     
         do i = 1,atypes
            ghr = ghr + atwei(i)*zion(i)*dion(i)
         enddo
c     
         dgr = dgr/ghr
c     
         fw = (atwei(zmap(6))*pahfrac*Nc)/ghr
c
         do i = 1,numtypes
          siggrain(i) = 0.d0
          do j=mindust(i),maxdust(i)
            dgrad=gradedge(j+1)-gradedge(j)
            siggrain(i)=siggrain(i)+dustsig(j,i)*dgrad
C            write(*,*) dustnum(j,i), dustsig(j,i)
          enddo
         enddo
         if (runmode.ne.'batchrun') then
c
 240     format(//'  Projected dust area/H atom (graphite,silicate):'/
     &            ' ****************************'/
     &            '  ',1pg12.4,1x,' (cm^2)',1x,1pg12.4,1x,' (cm^2)' )
         write(*,240) siggrain(1),siggrain(2)
c
 250     format( /'  Composition Ratios:  '/
     &            ' ****************************'/
     &            '   Gas/H     :',1pg12.4/
     &            '   Dust/Gas  :',1pg12.4/
     &            '   PAH/Gas   :',1pg12.4/
     &            ' ****************************'/
     &           ' (ratios by mass) ')
c     
         write(*,250) ghr,dgr,fw
c     
         endif
c     
c     deplete...
c     
         do i = 1,atypes
            if (i.eq.zmap(6)) then !seperate C because of PAHs
	      if (.NOT.ClinPAH) zion(i) = 
     &                    zion0(i)*(dion(i)+(1.d0-dion(i))*pahCfrac)
            else
              zion(i) = zion0(i)*dion(i)
            endif
         enddo
c
c  Calculate extinction & scattering cross-sections for use in totphot 
c     & newdif
c     
         do dtype=1,numtypes
          do k=mindust(dtype),maxdust(dtype)
            do inl=1,infph-1  
              dgrad=gradedge(k+1)-gradedge(k)
              dcrosec(inl,dtype) = dcrosec(inl,dtype) + 
c     &        (absorp(inl,k,dtype)+scatter(inl,k,dtype))
     &        absorp(inl,k,dtype)
     &             *dustsig(k,dtype)*dgrad
c              dfscacros(inl,dtype) = dfscacros(inl,dtype)
c     &             + 0.5d0*(1.d0+gcos(inl,k,dtype))*
c     &             scatter(inl,k,dtype)*dustsig(k,dtype)*dgrad
c              dbscacros(inl,dtype) = dbscacros(inl,dtype)
c     &             + 0.5d0*(1.d0-gcos(inl,k,dtype))*
c     &             scatter(inl,k,dtype)*dustsig(k,dtype)*dgrad
            enddo
          enddo
         enddo
c
c     Calculate visual extinction (per H column) (V=5470A=2.27eV)
c

        do inl=1,infph-1
          if ((ephot(inl).lt.2.2666435).AND.
     &        (ephot(inl+1).gt.2.2666435)) then
            dustAv = 2.5*dlog(dexp(1.d0))
     &               *(dcrosec(inl,1)+dcrosec(inl,2))
            goto 1010
          endif
        enddo
 1010   continue

c
c   Init fixed parameters for new IR temp calcs
c   in dustemp.f
c
      call dustinit()
c
      endif
c     
c     
      return
c     
      end
