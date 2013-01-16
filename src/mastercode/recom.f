cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c     
c     copyright 1994 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c*******TOCOMPUTE RADIATIVE AND DIELECTRONIC
c     RECOMBINATION RATES ; UNITS : CM**3 SEC**-1
c     RETURNS VECTOR REC(6,11) IN COMMON BLOCK /RECOMB/
c     
c     *FOR H AND HE , RATES IN REC(3,1),REC(4,2) AND REC(5,2)
c     FOR THE ON THE SPOT APPROXIMATION (WHEN USED)
c     
c     **REFERENCE : ALDROVANDI,S.M.V., AND PEQUIGNOT,D. (1973)
c     ASTRON. ASTROPHYS. V.25, 137
c     ASTRON. ASTROPHYS. V.47, 321
c     **REFERENCE : SEATON,M.J. (1959) MNRAS V.119, 81
c     (1962) MNRAS V.125, 437
c     **REFERENCE : NUSSBAUMER,H., AND STOREY,P.J. (1983)
c     SUBMITTED TO ASTRON. ASTROPHYS.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     allow old and new calculations of dielectronic rates
c     dielcmode = 0 means new calc based on
c     Landini & Monsignori Fosse, Arnausd and Rothenflug
c     & Shull
c     
c     mapmode = 1 means use old mappings methods
c     this is limited to the old Z range of atoms
c     
c     RSS 8/90
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     New total rec rates based on OP data for S and Ar. Use Spline
c     fits to rec calculations.  Still require High T dielectronic
c     recombination terms.
c     
c     RSS94
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      subroutine recom(t)
c     
      include 'cblocks.inc'
c     
      double precision t, u, u4, ar, xr, ad, bd, ta, tb
      double precision lx,lt,lr,x1,tp
      double precision a, b, c, d, e, f, r, tp4
      double precision ch,ct,et
      double precision fcor, adlt
c     
      integer*4 i,j,jj
      integer*4 ion,atom, maxio, chh,isos
c     
c
c           Functions
c
      double precision arf, adf, adltf, hifa, flowt, fkramer
c
      arf(u,r,e) = r*((u/1.0d4)**(-e))
c     
      adf(u,ad,bd,ta,tb) =
     &    (((1.0d0+(bd*dexp(-(tb/u))))*dexp(-(ta/u)))*ad)*(u**(-1.5d0))
c     
      adltf(u4,a,b,c,d,f) = ((1.0d-12*((((a/u4)+b)+(c*u4))+(
     &     (d*u4)*u4)))*(u4**(-1.5d0)))*dexp(-(f/u4))
c     
      hifa(u) = ((0.4288d0+(0.5d0*dlog(1.5789d5/u)))+(0.469d0*((
     &     1.5789d5/u)**(-(1.0d0/3.0d0)))))/1.9954d0
c     
      flowt(u,ct,et) = dmax1(1.d0,(((u/1.d4)/dsqrt(ct))
     &    **(-dmax1(0.0d0,et-0.5d0)))/dmax1(1.d0,hifa(u/dsqrt(ct))))
c     
c     keep t positive
c     
      tp = dmax1(0.01d0,t)
c     
      tp4 = tp/1.d4
c     
      do atom = 1, atypes
c   
        rec(1,atom) = 0.0d0
        maxio = maxion(atom)    
        if (mapz(atom).eq.2) maxio = 4
c    
        do ion = 2, maxio     
          isos = mapz(atom)-ion+1     
          rec(ion,atom) = 0.0d0
c   
c   Handle Hydrogenic ions separately
c    skip to end as there are no dielectronic rates 
c    for these ions     
c   
          if (isos.eq.0) then
             rec(ion,atom) = fkramer(ion,atom,t)
             goto 20
          endif
c   
          ar = arad(ion,atom)
          xr = xrad(ion,atom)    
          if ((mapz(atom).eq.2).and.(ion.eq.4)) then
             ar = 2.7d-13
             xr = 0.787d0
          endif
c
          if (ar.lt.0.0d0) goto 20
c   
          ad = adi(ion,atom)
          bd = bdi(ion,atom)
          ta = t0(ion,atom)
          tb = t1(ion,atom)

          if (oprec) then     
c     
c     Handle SII-VI and ArII-VI ions separately
c     
            if (((mapz(atom).eq.16).and.(ion.lt.7)).or.
     &           ((mapz(atom).eq.18).and.(ion.lt.7)))then
c     
c     Get right spline in i - one must exist 
              i = 1
              do j = 1,nrcspl
                if ((atom.eq.recsplat(j)).and.
     &              (ion.eq.recsplio(j))) then
                     i = j
                endif
              enddo
c     
c     Evaluate Spline at log(tp) and convert back to normal units
c     
c     Get nodes
              lt = dlog10(tp)
              x1 = recspl(i,1,1)
              jj = 1
              do j = 1,recspln(i)
                if (lt.ge.recspl(i,j,1)) then
                  x1 = recspl(i,j,1)
                  jj = j
                endif
              enddo
c     
c     interp and evaluate
              lx = lt - recspl(i,jj,1)
              lr = lx*recspl(i,jj,5)
              do j = 4,3,-1
                 lr = lx*(recspl(i,jj,j) + lr)
              enddo
              lr = lr+recspl(i,jj,2)     
              lr = 10.d0**lr
c     
c     Add High T Dielec
              rec(ion,atom) = lr+adf(tp,ad,bd,ta,tb)
c     
c     Sanity Check
              rec(ion,atom) = dmax1(0.d0,rec(ion,atom))
              goto 20
            endif
c
c   End seperate Sulfur and Argon treatment	  
	  endif
c
c     Do remaining Ions
c          
c     1st check for low t corr     
          if (t.lt.1.d5) then
c     
c     use low t corrections...
            a = adilt(ion,atom)
            b = bdilt(ion,atom)
            c = cdilt(ion,atom)
            d = ddilt(ion,atom)
            f = fdilt(ion,atom)
c     
            chh = ion-1
            if (chh.gt.maxion(atom)-1) chh=chh-(maxion(atom)-1)
            ch = dble(chh)
c        
            fcor = flowt(tp,ch,xr)
c     
            adlt = adltf(tp4,a,b,c,d,f)
            if (adlt.lt.0.0d0) adlt = 0.0d0
c        
            rec(ion,atom) = (arf(tp,ar,xr)*fcor)+adf(tp,ad,bd,ta,tb)
     &                      +adlt
c
c
         else
c     
c     no low t correction needed...
c     
            rec(ion,atom) = arf(tp,ar,xr)+adf(tp,ad,bd,ta,tb)
c     
         endif
c     
c     
c     End Recombination Rate Calcs
c     
 20      continue
c
c     Sanity Check
c     
         rec(ion,atom) = dmax1(0.d0,rec(ion,atom))
c     
c     Debug and test
c
c      if (atom.eq.1) then
c     
c         write(*,*) t,' ',rec(ion,atom),
c     &       ' ',elem(atom),rom(ion)
c
c      endif
c     
c     write(*,*) ar,xr
c     write(*,*) a,b,c,d,f
c     write(*,*) ad,bd,ta,tb
c     arfc = arf(tp,ar,xr)
c     adfc = adf(tp,ad,bd,ta,tb)
c     write(*,*) t,' ',elem(atom),rom(1), rec(1,atom)
c     write(*,*) t,' ',elem(atom),rom(2), rec(2,atom)
c     write(*,*) t,' ',elem(atom),rom(3), rec(3,atom)
c     write(*,*) t,' ',elem(atom),rom(ion),
c     &                 rec(ion,atom),arfc,adfc,fcor,adlt
c     endif
c      endif
c
c     end ion
c
        enddo
c     
c	   end atom
c
      enddo
c     
c     
      return
c     
      end
