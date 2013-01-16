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
c
c
c     Subroutine freebound calculates the free-bound continuum contribution
c     to the local diffuse field..  H like atoms are treated with lookuptables
c     with muliple levels.  Other atoms are derived from photo- crossections
c     with an additions excited level contribution summed into the first
c     edge above ground.
c
c     adds the emission to cnphot (in photons/cm3/s/Hz/sr)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine freebound(t,de,dh)
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision t,de,dh,bev,en2,abde
      double precision gfbh,gfbi,enh,zpop
      double precision rydu,rz3,zn2,t12
      double precision z2,z4,fbg
      double precision f1,df,ri,tz,fbconst
      integer*4 i,it,atom,ion,n0,ngfbt,isos
c
      double precision abio,at,bet,cebin,crosec,enu
      double precision eph,excerg,et,phots
      double precision phoi,rkt,se,statwei
      integer*4 io,ie,j,inl
c
c     External Functions
c
      integer*4 GetIndexUpwards
      double precision fgfblog
c
c     Internal Functions
c
      double precision acrs,hir
c
c
c  -138.446097d0 = ln((h^3)/(me^3 c^2)  (me/2)^(1.5)  4/sqrt(¹) /(4¹))
c
c	Osterbrok Appendix 1
c
c
      hir(rkt,enu,excerg) = dexp((-138.446097d0+(2.0d0*dlog(enu)))-
     &                          (1.5d0*dlog(rkt))-(excerg/rkt))
      acrs(eph,at,bet,se) = (at*(bet+((1.0d0-bet)/eph)))/(eph**se)
c
      rkt    = rkb*t
      rydu   = (ryde)/(rkt)
      t12    = 1.d0/dsqrt(t)
c
      do inl = 1,infph
         fbph(inl) = 0.d0
      enddo
c
c     The Reimann zeta function for argument 3
c
      rz3 = 1.202056903159594d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ground level fb continuum from photo X-sections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ground level free-bound continuum. Transplant from localem.
c     This is better because it is based on the actual photo
c     cross-sections in use, rather than an H approximation.
c
         do inl = 1,infph-1
c
c     check energy compared to kt
c
            cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
            enu   = cebin*evplk
            bev   = cebin*ev
            et    = (bev)/(rkt)
c
            if (dabs(et).lt.loghuge) then
c
c     **ADS RECAPTURE TO GROUND STATE FOR ELEMENTS UP TO*LIMPH
c     STATWEI = RATIO OF STATISTICAL WEIGHT OF RECOMB. TO ION
c     EX. HYDR.: STATWEI 2/1=2
c     
              fbg = 0.0d0
c
              do io = 1, ionum
c
               ie      = atpho(io)
               j       = ionpho(io)
               isos    = (mapz(ie) - j + 1)
               abio    = zion(ie)*pop(j+1,ie)
	           statwei = stwtpho(io)
c	           
               if ((abio.gt.pzlimit).and.(statwei.gt.0.0)) then
c               
	          phoi     = ipotpho(io)
	          if (cebin.ge.phoi) then
	             excerg = ev*(cebin-phoi)
	             eph    = cebin/phoi
	             crosec = acrs(eph,sigpho(io),betpho(io),spho(io))
	             fbg    = fbg+abio*crosec*hir(rkt,enu,excerg)
     &                        *statwei
	          endif
c
               endif
c
              enddo
c            
c               
c    not collecting Mewe style gaunt factors normally - place after 
c     phots line when using
c
c        fbgau(inl)  = fbgau(inl) + bev*phots/(ffk*de*de*t12*dexp(-et))
c              
            phots       = de*dh*fbg
            fbph(inl)   = fbph(inl)   + phots
            cnphot(inl) = cnphot(inl) + phots
        endif
c            
        enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Excited level fb continuum
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     loop through all ions present
c
      do atom = 1 , atypes
         do ion = 2 , maxion(atom)
c
            zpop = zion(atom)*pop(ion,atom)
c
            if (zpop.ge.pzlimit) then
c
               abde    = zpop*de*dh
               fbconst = ffk*abde*t12
c 
               isos = (mapz(atom) - ion + 1)
c
c     get secondary ionisation edge
c
               en2  = epot2(ion-1,atom)*ev
c
c     calculate or get effective charge (z) ,n0 and gamma (alias y)
c
               zn2 = zn0(ion-1,atom)
               zn2 = zn2*zn2
c
               z2 = dble(ion)-1.0d0
               z2 = z2*z2
               z4 = z2*z2
c
c     calculate f1 factor (a sum over n)
c
c     term f1 (reimann zeta(3) - sum n^-3 up to no)
c
               f1   = 1.d0
               n0   = fn0(ion-1,atom)
               do i = 2, n0
                  ri = dble(i)
                  f1 = f1+1.0d0/(ri*ri*ri)
               enddo
c
               f1 = rz3 - f1
c
               do inl = 1,infph-1
c
c     check energy compared to kt
c
                  cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
                  bev   = cebin*ev
                  et    = (bev)/(rkt)
c
                  if (dabs(et).lt.loghuge) then
c
                  gfbh = 0.d0
                  if (isos.eq.0)then
c
c     Add higher level h recombination from lookup tables
c              
                     enh  = dlog10(min((bev/epot(ion-1,atom)),1.d0))
c
                     if (enh.gt.gfbe(1))then
                       tz  = dlog10(t/zn2)
                       if (tz.lt.5.69897d0) then
c
c     less than 500,000K/Z^2
c
                         ngfbt = (ngfbe/2)+1
                         ie    = GetIndexUpwards(enh, ngfbt  , gfbe)
                         it    = GetIndexUpwards(tz , ngfbtz , gfbtz)
                         gfbh  = zn2*fgfblog(tz,enh,it,ie)
c
                       endif
                     endif
                  endif
c                 
c     only calculate above 2nd ionisation edge for other elements
c
                  df  = 0.d0
                  if (bev.ge.en2) then
                     df = 2.d0*z4*f1*dexp(en2/(rkb*t))
                  endif
c
c     Correct gfbh by 2  Will fix table soon.
c
                  gfbi = 0.d0
c
                  if (isos.eq.0) then
c
c     Hydrogenic Series with Lookup Table
c
                     gfbi = 2.0*gfbh
c
                  else 
                     gfbi = rydu*df
                  endif
c            
c               
c    not collecting Mewe style gaunt factors normally - place after 
c     phots line when using
c
c       fbgau(inl)  = fbgau(inl) + bev*phots/(ffk*de*de*t12*dexp(-et))
c              
                  phots       = fbconst*gfbi*dexp(-et)/bev
                  fbph(inl)   = fbph(inl)   + phots
                  cnphot(inl) = cnphot(inl) + phots
c
c              end et limit
c
                  endif
c
c              end radiation loop
c
               enddo
c
c     end if loop for zpop > pzlimit
c
         endif
c
c     end ion and atom loops
c
         enddo
c
      enddo
c
      return
c
      end
      
