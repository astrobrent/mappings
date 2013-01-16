cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	New resonance line data in XLINDAT, requires new 
c	code so that the effective osc strength f' can be used 
c	effectively.
c
c
c	calls the new function fresga to get gaunt factors needed.
c
c
c	refs: Landini & Monsignori Fosse 1990 A.A.Suppl.Ser. 82:229
c		  Mewe 1985 A.A.Suppl.Ser. 45:11
c
c	RSS 8/90
c
c	NOTE: transition energy does not give Egj, it gives Ejk
c	and in some cases this may not be a good approximation.....
c     ie transition H6
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES RESONANCE LINE COOLING
c       NB.  SUBR. HYDRO SHOULD BE CALLED PREVIOUSLY
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine reson2(t, de, dh)
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision t, de, dh
      double precision abup,pirt38,rr,y
      double precision egj,ejk,tr,fef,gbar,fresga
      double precision cgjb,omegab,e1
c
      integer*4 i
      integer*4 atom,ion,isos,code
c
c           Functions
c
      double precision farint
c
      xrloss = 0.0d0
      pirt38 = pi*8.d0/dsqrt(3.d0)
      t = dmax1(0.1d0,t)
c
c
      do i = 1,xlines
c     
c	get data from arrays for line # i
c     
c     Note: only line for possible species are read in
c     so no need to check maxion etc...
c     
         atom = xrat(i)
         ion = xion(i)
         isos = xiso(i)
c     
         xrbri(i) = 0.d0
c
c
c
         abup = zion(atom)*pop(ion,atom)
c
c     only calculate abundant species
c     

         if (abup.ge.pzlimit) then
c     
c     note that Egj does not necesarily equal Ejk
c     
            egj = xegj(i)
            ejk = xejk(i)
            tr = xtrans(i)
            fef = xfef(i)
c
c
c     
c     get scaled energy gap to level j from ground (not k)
c     
            y = egj*ev/(rkb*t)
c
c     
c     
c     get mean gaunt factor
c     
            gbar = fresga(atom,isos,tr,egj,t,code)
c
c
c     
c     transition power rate rr
c     
            rr = 0.d0
c     
            cgjb = 0.d0
c
            if (code.ge.0) then
               omegab = pirt38*fef*gbar*ryd/egj
            else 
               omegab = fef
            endif
c
c     Special cases (sorry)
c
c
c     CIV 1549/51
c
            if ((mapz(atom).eq.6).and.(ion.eq.4)
     &           .and.(egj.lt.10.d0)) then
            e1 = farint(1,y)
            omegab = 0.6790+0.2364*y
            omegab = omegab+((0.4202-0.04195*y-0.2364*y*y)*e1)
            omegab = 9.3685*omegab
            if(xrlam(i).gt.1.550E+03) then
               omegab = (omegab/3.0)
            endif
            if(xrlam(i).lt.1.550E+03) then
               omegab = omegab*(2.0/3.0)
            endif
            endif            
c
c     OVI 1031.9/1037.6
c
            if ((mapz(atom).eq.8).and.(ion.eq.6)
     &           .and.(egj.lt.13.d0)) then
            e1 = farint(1,y)
            omegab = 0.6548+0.1358*y
            omegab = omegab+((0.3272+0.1598*y-0.1358*y*y)*e1)
            omegab = 5.483*omegab
            if(xrlam(i).gt.1.035E+03) then
               omegab = (omegab/3.0)
            endif
            if(xrlam(i).lt.1.035E+03) then
               omegab = omegab*(2.0/3.0)
            endif
            endif            
c
            if (y.lt.loghuge) cgjb = rka*(dexp(-y)/dsqrt(t))
c
            rr = cgjb*omegab
            
c     number to transition is in abup
c     
c     
c     photons cm^-3 s^-1
c     
            rr = de*dh*abup*rr
c
c     ergs...
c
            rr = ev*ejk*rr
c
c
c            if (hellcoolmode.eq.1) then
c
c     cooling for He in opt thin case
c
c               if ((mapz(atom).ne.1).and.(mapz(atom).ne.2)) then
c     
c     cooling contribution from H&He  lines are done in hydro
c
c
c                  xrloss = xrloss+rr
c
c               endif
c            endif
c     
c            if (hellcoolmode.eq.0) then
c
c     cooling for He in opt thick case
c
c               if ((mapz(atom).ne.1)) then
c     
c     cooling contribution from H  lines are done in hydro
c
c
c                  xrloss = xrloss+rr
c
c               endif
c            endif
c
c
            xrloss = xrloss+rr
            xrbri(i) = rr/(4.0d0*pi)
c     
c     
c     end population limited loop
c
         endif
c     
c     
      enddo
c     
c     
      return 
c
      end
