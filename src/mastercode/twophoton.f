cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine TWOPHOTON calculates the two-photon continuum 
c     contribution to the local diffuse field for H and He-like ions.
c     It adds the emission to CNPHOT (in photons/cm3/s/Hz/sr).
c
c     RSS 2000: Modified to only do He-like 2p emission from helium until 
c     heavies are sorted out.
c
c     Must be called soon after HYDRO so that COLLRATE and RECRATE are 
c     up to date, i.e., from CONTCONT from LOCALEM.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine twophoton(t,de,dh)
c
      include 'cblocks.inc'
c
      double precision t,de,dh,bev,cebin
      double precision fpsiy,invpi4,rkt,t12
      double precision y,nz2,z6,pcoll
      double precision r2q,ab,ab0,ab1,abp,tz
      double precision u,phots,omes1,x
c
      integer*4 atom,nz,inl
c
c     old
c      double precision et,br0,ez,f
c      double precision a21,ara,cos1,frs1
c      double precision omep1,ra3,ras1
c      double precision reche,tns3,es1,fgaunt
c
c     Internal functions for hydrogen.
c
      double precision qpr,qel
c
c     Functions for Helium I
c
      double precision pol5
      double precision qss,qsp,frs3,fusp
c
      qpr(u) = 4.74d-4*(u ** (-0.151d0))
      qel(u) = 0.57d-4*(u ** (-0.373d0))
c     
      pol5(x) = dmax1(0.0d0,(((((((((-(6.5289d0*x))+41.54554d0)*x)-
     &  97.135778d0)*x)+97.0517d0)*x)-32.02831d0)*x)-0.045645d0)
c
      qss(u) = 1.d-8*pol5(u ** 0.3333333d0)
      qsp(u) = 5.73d-9*pol5((0.175d0+(0.2d0*u)) ** 0.3333333d0)
c
      frs3(u) = 0.775d0*(u ** 0.0213d0)
      fusp(u) = 0.42d0*(u ** 0.23d0)
c
      invpi4 = 1/(4.d0*pi)
c
c     get scaled energies wrt Te
c
      rkt = rkb*t
      t12    = 1.d0/dsqrt(t)
c
      u   = bev/(rkt)
c
      do inl = 1,infph
         p2ph(inl) = 0.d0
      enddo
c
      omes1 = 0.332d0
c
      do atom = 1,atypes 
c
c   H like atoms
c     
         nz    = mapz(atom)
         nz2   = 1.d0/dble(nz*nz)
         ab1   = zion(atom)*pop(nz+1,atom)
         ab0   = zion(atom)*pop(nz,atom)
         ab    = 0.5d0*(ab0+ab1)
c
         r2q = 0.d0
c     
         if (ab.gt.pzlimit) then
c     
c     use precomputed rates from hydro.f
c
            r2q   =  ab0*collrate2p(atom) + ab1*recrate2p(atom)
c
c     scale temperature 
c     
            tz   = 1.d-4*t*nz2
c     
c     proton collision
c     
            abp   = dh*zion(1)*pop(2,1)
            pcoll = (abp*qpr(tz)+de*qel(tz))*nz
c     
c     A scales as nz^6
c     
            nz   = nz*nz
            z6   = nz*nz*nz
c
            r2q  = r2q/(1.0d0+(pcoll/(8.227*z6)))

            r2q  = de*dh*r2q*invpi4
c
            do inl = 1,infph-1
               cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
               bev   = cebin*ev
               y     = bev/(0.75d0*epot(mapz(atom),atom))               
               if (y.lt.1.d0) then
                   phots       = 2.d0*plk*y*fpsiy(y)*r2q/bev
                   p2ph(inl)   = p2ph(inl) + phots
                   cnphot(inl) = cnphot(inl) + phots
               endif
            enddo
         endif
      enddo
      
c
c     The following assumes that HELIOI has been called.
c     This includes the de*dh/4pi etc. already - see helioi line 79.
c      
       r2q   = hein(2)/heien(2)

       do inl = 1,infph-1
c
         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
         bev   = cebin*ev
         y     = bev/heien(2)
c         
         if (y.lt.1.d0) then
             phots       = 2.d0*plk*y*fpsiy(y)*r2q/bev
             p2ph(inl)   = p2ph(inl)   + phots
             cnphot(inl) = cnphot(inl) + phots
         endif
       enddo
c
c     He I collision strengths.
c
c      omep1 = 2.57d0
c      omes1 = 0.55d0*omep1
c
c
c      do atom = 2, atypes
c      atom = 2
c
c   He like atoms
c     
c         nz    = mapz(atom)-1
c         nz2   = dble((nz+1)*(nz+1))
c         ab1   = zion(atom)*pop(nz+1,atom)
c         ab0   = zion(atom)*pop(nz,atom)
c         ab    = 0.5d0*(ab0+ab1)
c         ez    = ezz(nz,atom)
c         es1   = ez/rkt
c
c         r2q = 0.d0
c     
c         if (ab.gt.pzlimit) then
c
c            u     = dabs(t+epsilon)/nz2
c            f     = dsqrt(1.0d0/u)
c            
c            if (es1.lt.loghuge) then
c     
c            br0   = (1.d0-1.34d0/mapz(atom))
c            z6    = nz2*nz2*nz2
c
c    ***COLLISIONAL RATES FROM GROUND STATE 1S
c
c            ara   = ab0*rka*f
c         
c            cos1  = ara*br0*(omes1/nz2)*dexp(-es1)*fgaunt(1,1,es1)
c
c   Recombination into triplet and singlet systems
c
c            tz    = 1.d-4*u
c            if (tz.lt.0.1) tz = 0.1d0
c            if (tz.gt.10.0) tz = 10.d0
c            
c            reche = ab1*rec(nz+1,atom)
c            frs1  = (reche*fusp(tz))*(1.0d0-frs3(tz))
c            ra3   = reche*frs3(tz)
c            a21   = 1.27d-4 * z6
c            tns3  = ra3/(1.d0+(de*(qss(tz)+qsp(tz)))/a21)
c            ras1  = tns3*qss(tz)
c 
c            r2q   = de*dh*(ras1+frs1+cos1)*invpi4
c            write(*,*) "HeI atom ras1  frs1  cos1",atom, ras1, frs1, cos1
c            write(*,*) "HeI atom r2q",atom, r2q
c            
c            do inl = 1,infph-1
c
c               cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
c               bev   = cebin*ev
c               y     = bev/ez
c               
c              if (y.lt.1.d0) then
c
c                   phots       = 2.d0*plk*y*fpsiy(y)*r2q
c                   p2ph(inl)   = p2ph(inl) + phots
c                   cnphot(inl) = cnphot(inl) + phots
c
c               endif
c            
c            enddo
c
c           endif
c
c         endif
c
c   He like atoms
c     
c      enddo
c
      return
c
      end
