cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.1.0md
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Program to determine dust temperature distribution using the 
c     algorithm and equations of Draine & Li 2001, except use H for 
c     enthalpy not U (which is grain potential here)
c
c     Output is IRphot in phot s-1 cm-3 Hz-1 Sr-1 (4pi)
c    
c     Brent Groves 5/11/01 
c
      subroutine dustinit()
c
c Call once dust setup is known to save repeat intialisation
c
      include 'cblocks.inc'
c
      double precision massatom(2)
      integer*4 dtype, k,l,inl

c
      massatom(1)=12.d0*amu            ! mass of atom (assume pure C ) (g)
      massatom(2)=60.0855d0*amu        ! mass of atom (assume SiO2) (g) 
      do l=1,dustbinmax
       grainvol(l)=ftpi*(grainrad(l))**3 !grain grainvolume in cm^3
       d_rad(l)=gradedge(l+1)-gradedge(l)
       do dtype=1,numtypes
c    number of atoms in the grain
        atom_no(dtype,l)=Grainvol(l) * graindens(dtype)/massatom(dtype) 
       enddo
      enddo
      
c
c     setup energy bns for quick calculation
c
      do inl=1,infph-1
       dustsigt(inl)=0.d0
       do dtype=1,numtypes
         dustsigt(inl)=dustsigt(inl)+dcrosec(inl,dtype)
       enddo
      enddo
c
      delE = ephot(infph)/ephot(1)
      delE = dlog10(delE)/dble(dinfph-1)
      edphot(1) = ephot(1)
      do k=1,dinfph-1
        edphot(k+1) = ephot(1)*10.d0**(k*delE)
        dengy(k)=0.5d0*(edphot(k+1)+edphot(k))*eV
        ddE(k)=(edphot(k+1)-edphot(k))*eV
        v_e(k)=dsqrt(2.d0*dengy(k)/me)  !electron velocity (cm s-1)
      enddo

c
      return
     
      end

      
      subroutine dusttemp(T_e,Hdens,n_e, fi, dr,IRcount)

      include 'cblocks.inc'
c

      double precision e_kt,wid,den
      double precision delta_T,dT,invkT
      double precision Tmin,Tmax, dTmin, dTmax
      double precision dTmin0, dTmax0
c

      double precision Hdens,grarea,fi,dr
      double precision n_e,T_e, s_f, phi
      double precision newflx, oldflx, flxratio
      double precision flxtest1, flxtest2, flxabs
      double precision Teq, Tlimit
      double precision totaldust,flxabstot
c
      integer*4 dtype,IRcount,Nmax,Nmax0,Nlimit
      integer*4 absmax,dabsmax
      integer*4 i,j,k,l,inl
c
      character tempfile*20
      character pfx*8,sfx*4
c
      logical Tfine, QuickIR
c
c     old
c       double precision dcool(dinfph, MxTempbin)
c
       double precision MBdist
c
c      functions
c
       double precision planck, a, b
c
       planck(a,b) = (kbb*b*b*b*2)/(dexp(b*a)-1.d0)
c
c  Dust parameters 
c
c
c      massatom(1)=12.d0*amu            ! mass of atom (assume pure C ) (g)
c      massatom(2)=60.0855d0*amu        ! mass of atom (assume SiO2) (g) 
c      do l=1,dustbinmax
c       Grainvol(l)=ftpi*(grainrad(l))**3 !grain volume in cm^3
c       d_rad(l)=gradedge(l+1)-gradedge(l)
c       do dtype=1,numtypes
c    number of atoms in the grain
c        atom_no(dtype,l)=Grainvol(l) * graindens(dtype)/massatom(dtype) 
c       enddo
c      enddo
c
c  Determine IR mode
c
      QuickIR=.True.
c
      if (IRmode.eq.1) then
        QuickIR=.True.
        Nlimit=50
      else if (IRmode.eq.2) then
        QuickIR=.True.
        Nlimit=100
      else if (IRmode.eq.3) then
        QuickIR=.False.
      else
        write(*,*) 'IR program incorrect'
        return
      endif
c
      if (IRtemp) then
         pfx='temp'
         sfx='dat'
         call newfile(pfx,4,sfx,3,tempfile)
         open(15,file=tempfile,status='unknown',access='APPEND')
      endif
c
c  Initialise Temperature bins
c			
      Tmin   =    1.d0
      Tmax   = 1201.d0
c
c   force initial temp/grid calcs.
c
        dTmin0 = 0.0d0
        dTmax0 = 0.0d0
c
c     setup energy bns for quick calculation
c
      flxtest1 = 0.d0
      flxabstot   = 0.d0
      do inl=1,infph-1
       wid=(ephot(inl+1)-ephot(inl))*evplk
       phot(inl)= 4.d0*pi*dustphot(inl)*wid     
c  convert from Jnu (1/4pi) to Fnu.dnu (erg s-1 cm-2)
       engy(inl)=0.5d0*(ephot(inl+1)+ephot(inl))
       if (engy(inl).gt.Bgrain) then
         phot(inl)=phot(inl)*(1.d0-Yinf*(1.d0-Bgrain/engy(inl))) 
c  losses due to PE effect
       endif
       photabs(inl) = phot(inl)
       flxtest1=flxtest1 +phot(inl)
       engy(inl)=engy(inl)*eV
      enddo
c      write(*,*)  ' dusttemp flux:',flxtest1
 
c
c      delE = ephot(infph)/ephot(1)
c      delE = dlog10(delE)/dble(dinfph-1)
c      edphot(1) = ephot(1)
c      do k=1,dinfph-1
c        edphot(k+1) = ephot(1)*10**(k*delE)
c        dengy(k)=0.5d0*(edphot(k+1)+edphot(k))*eV
c        ddE(k)=(edphot(k+1)-edphot(k))*eV
c      enddo
c
c
c  set_up collisional heating bins
c
      do k=1,dinfph-1
c       v_e=dsqrt(2.d0*dengy(k)/me)  !electron velocity (cm s-1)
       collheat(k) = MBdist(T_e,dengy(k))*v_e(k)*dengy(k)*ddE(k)
      enddo

c
      if(T_e.ne.0.0)then
        e_kt = eV/(rkb*T_e)
      else
        e_kt = 0.d0
      endif
c
c  Initialise Region Flux output
c



      do k=1,infph-1
        IRFlux(k)=0.d0
      enddo
c
c  Do for both Silicates and Graphite grains
c   Type=1=Silicates, Type=2=Graphite
c
      totaldust=0.d0
      do 3001 dtype=1,numtypes
c
c  Set-up Temp range
       dTmax=Tmax

c
c  determine temperature for grain sizes considered 
c      
       do l=mindust(dtype),maxdust(dtype)
        totaldust=totaldust+dustsig(l,dtype)*d_rad(l)
c
c  Calculate sticking factor
c
        phi = grainpot(l,dtype)*e_kt
        if (phi.ge.0.d0) then
           s_f = 0.5d0*(1.d0+phi)
        else
           s_f = 0.5d0*dexp(phi)
        endif
c
c Areas in sq. microns within Transmatrix
c 
        grarea=dustsig(l,dtype)*d_rad(l)
        do inl=1,infph-1
         Cabs(inl)=grarea*absorp(inl,l,dtype)
         if (Cabs(inl).ne.0.d0) absmax=inl
        enddo
c 
c     Set up absphot bins
c
        flxabs=0.d0
        j=1
        do k=1,dinfph-1
          dabsphot(k)=0.d0
 90       den=0.5d0*(ephot(j+1)+ephot(j))
           if (dustsigt(j).eq.0.d0) goto 95
          if (den.lt.edphot(k+1)) then
c           write(*,*) Cabs(j),dustsigt(j),photabs(j)
c
c  determine energy in scaled dust bin
c  using fraction of absorption due to grain l (ie sigma_l/sigma_tot)
c
            dabsphot(k)=dabsphot(k)
     &       + (Cabs(j)/dustsigt(j))
     &       * (photabs(j)/(dr*fi*Hdens)) !4pi in phot
            flxabs=flxabs
     &       + (Cabs(j)/dustsigt(j))
     &       * (photabs(j)/(dr*fi*Hdens)) !4pi in phot
           if (dabsphot(k).ne.0.d0) dabsmax=k
           j=j+1
           if (j.gt.(infph-1)) goto 95
           goto 90
          endif

        enddo
 95     continue
        flxabstot=flxabstot+flxabs
        Teq = (flxabs/(grarea*stefan))**(0.25)
          if (IRtemp) then
 105              format(i4,2(1pg15.7,x))
           write(15,105) dtype,grainrad(l),Teq
          endif
c         write(*,2000) dtype,l,flxabs,Teq
c 2000   format(2(i4,x),'flxabs= ',1pg12.4,'Teq= ',1pg12.4)
c
c  setup initial Temperature & flxmeasures
c     (maxmimum dTmax in smallest (first) grain)
c 

        Nmax=50
        Nmax0=Nmax
        dTmin=Tmin
        oldflx=1.d2
        newflx=1.d0
        flxratio=1.d2
        Tfine = .True.
        
        
c
c  set P(T) Tlimit
c
C        if (grainrad(l).lt.1d-5) then
          Tlimit=1.0d-15
C        else
C          Tlimit=1.0d-10
C        endif

        do 173 while (Tfine)
c
c     Clear Flux distribution
c
         do k=1,infph-1
           flxdist(k,l)=0.d0
         enddo
c
c recompute Temp grids if either dTmax or dTmin change
c

          if (( dTmin .ne. dTmin0).or.(dTmax .ne. dTmax0)
     &        .or.(Nmax.ne.Nmax0)) then
c
c  setup Temp grid
c 
          delta_T= (dTmax-dTmin)/Nmax ! temp grid bin separation
          dT    = delta_T*0.5d0 ! 0.5 temp grid bin separation

c          write(*,*) 'tmax tmin delt dt', dTmax, dTmin, delta_T, dT

          T_edge(1)=dTmin
          do i=1,Nmax
            T_edge(i+1)=dTmin+i*delta_T
            T_grid(i)=T_edge(i+1)-dT ! make temp grid
          enddo
c
c  Setup grid for Plank spectrum at all dust temp. and photon energies
c
c           open(69,file='dcool.dat')
          do i=1,Nmax 
            invkT = 1.d0/(rkb*T_grid(i))
            do k=1,infph-1
              BBemiss(k,i)=planck(invkT,engy(k))
            enddo
          enddo
c
c record Ts for most recent BB and grid calcs
c

          dTmin0 = dTmin
          dTmax0 = dTmax
c
          endif
c
          call initgrids(dtype,Grainvol(l),atom_no(dtype,l),T_grid,T_edge,
     &       H_grid,Hmin,Hmax,deltaH,Nmax,mxTempbin)
      
          do j=1,Nmax
             do i=1,Nmax
               invH(i,j) = 0.d0
              if ( i .ne. j)invH(i,j) = 1.d0/(H_grid(i) - H_grid(j))
             enddo
          enddo
              

c
c     solve for transition matrix
c
          if (dabsmax.ne.1) then
           call transmatrix(grarea, n_e ,s_f, Nmax, absmax,dabsmax,t_e)

c  With transition matrix obtained now solve for P vector
             call probsolve(Tr_matrix,T_prob,Nmax, mxTempbin)
          endif
c
c     Calculate Flux distribution and total flux
c
          newflx=0.d0
          do j=1,Nmax
           if (T_prob(j).ne.0.d0) then
            do k=1,infph-1
             if(BBemiss(k,j).lt.1.d-50) goto 520    !past exp dropoff
              flxdist(k,l)=flxdist(k,l)+T_prob(j)*BBemiss(k,j)
     &           *absorp(k,l,dtype)*grarea
            enddo
 520        continue
           endif
          enddo 
          do k=1,infph-1
            wid=(ephot(k+1)-ephot(k))*ev
            newflx=newflx+flxdist(k,l)*4*pi*wid
          enddo
c
c
c
c          if (IRtemp) then
c           pfx='temp'
c           sfx='dat'
c           call newfile(pfx,4,sfx,3,tempfile)
c           open(15,file=tempfile,status='unknown',access='APPEND')
c           write(15,*) Nmax
c 105       format(2(1pg15.7,x))
c           do j=1,Nmax
c            write(15,105) T_grid(j),T_prob(j)
c           enddo
c           close(15)
c          endif

c
c     QuickIR using only 50 bins for energy conservation
c
          if (QuickIR) then
           k=Nmax
           do while (T_prob(k).lt.Tlimit)
             dTmax=T_edge(k)
             k=k-1
           enddo
           j=1
           do while (T_prob(j).lt.Tlimit)
            dTmin=T_edge(j+1)
            j=j+1
           enddo
c  In case hotter than expected
           if ((k.eq.Nmax).AND.(dTmax.eq.Tmax)) dTmax=Tmax+500  
           if ((j.lt.4).AND.((Nmax-k).lt.4)) then
              Nmax=Nmax*2.d0
           endif
           if (Nmax.gt.Nlimit) Tfine = .False.
c
c  Not quick do full IR for temp distributions
c
          else  
           if (Nmax.lt.200) then
            k=Nmax
            do while (T_prob(k).lt.Tlimit)
             dTmax=T_edge(k)
             k=k-1
            enddo
            j=1
            do while (T_prob(j).lt.Tlimit)
             dTmin=T_edge(j+1)
             j=j+1
            enddo
c  In case hotter than expected
            if ((k.eq.Nmax).AND.(dTmax.eq.Tmax)) dTmax=Tmax+500  
            if ((j.lt.4).AND.((Nmax-k).lt.4)) Nmax=Nmax*2.0
           else
            k=Nmax
            do while (T_prob(k).lt.1.d-15)
              dTmax=T_edge(k)
              k=k-1
            enddo
            j=1
            do while (T_prob(j).lt.1.d-15)
             dTmin=T_edge(j+1)
             j=j+1
            enddo
            if ((j.lt.(Nmax/20)).AND.
     &          ((Nmax-k).lt.(Nmax/20))) Nmax=Nmax*2.d0
           endif
           if (Nmax.gt.400) Tfine = .False.
          endif
  173    enddo
c
c       flxtest1=flxtest1+flxabs
c
c  Absorbed & Emitted per Grain
c
c       flxtest2=0.d0
c       do k=1,infph-1
c          wid=(ephot(k+1)-ephot(k))*ev
c          flxtest2=flxtest2+flxdist(k,l)*wid
c       enddo
c       flxtest2=flxtest2*4*pi*d_rad(l)
c       write(*,*) grainrad(l),'E IRFlux Diff =', 
c     &    (1.d0 - flxtest2/flxabs)
c
c
      enddo   !dust size loop

c
c  Finally determine radiation field
c   (use trapezoidal integration over dust sizes)
c
       if(mindust(dtype).ne.maxdust(dtype))then
        do l=mindust(dtype),maxdust(dtype)
         do k=1,infph-1
           if (flxdist(k,l).lt.1.0d-50) goto 540    !past exp dropoff
c           IRFlux(k)=IRFlux(k)+0.5d0*d_rad(l)*(flxdist(k,l-1)
c                  +flxdist(k,l))
            IRFlux(k)=IRFlux(k)+flxdist(k,l)
         enddo
 540     continue
        enddo
       else
        l=mindust(dtype)
        do k=1,infph-1
          if(flxdist(k,l).lt.1.d-50) goto 620
            IRFlux(k)=IRFlux(k)+flxdist(k,l)
        enddo
 620    continue
       endif
 3001 continue
      
c
c  Detemine IRphot (photons/s/cm^3/Hz)
c     flux in ergs s-1 Hatom-1 Hz-1 Sr-1 (due to dustsig)
c
c
      flxtest2=0.d0
      do k=1,infph-1
        wid=(ephot(k+1)-ephot(k))*ev
        flxtest2  = flxtest2+IRFlux(k)*4*pi*wid
        IRphot(k) = IRFlux(k)*Hdens*plk/engy(k)
      enddo
      
      write(*,1001) dr*fi*Hdens*flxabstot,flxtest1,dr*fi*Hdens*flxtest2
 1001 format('IRin1=',1pg12.5,' IRin2=',1pg12.5,' IRout=',1pg12.5)
        Teq = (flxtest2/(totaldust*stefan))**(0.25)
      write(*,1002) Teq
 1002 format('Average Dust Temp: ',1pg14.7)

          if (IRtemp) then
           close(15)
          endif
      return
      end

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine intialises Grids and grid values
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initgrids(dtype,V,atom_no,T_grid,T_edge,
     & H_grid,Hmin,Hmax,deltaH,Nmax,MxBin)
      implicit none
      integer*4 dtype,Nmax, i, MxBin

      double precision V,atom_no,T_grid(MxBin),T_edge(MxBin),
     & H_grid(MxBin),Hmin(MxBin),Hmax(MxBin),deltaH(MxBin)
      
      
c     external functions      
      double precision Sil_enth,Gra_enth
      

c  calculate enthalpy for each temperature 
      if(dtype.eq.1)then
       do  i=1,Nmax
         H_grid(i)=Sil_enth(atom_no,v,T_grid(i))         
         Hmin(i)=Sil_enth(atom_no,v,T_edge(i))
         Hmax(i)=Sil_enth(atom_no,v,T_edge(i+1))
         deltaH(i) = Hmax(i) - Hmin(i)
       enddo
      else
       do  i=1,Nmax
         H_grid(i)=Gra_enth(atom_no,v,T_grid(i))         
         Hmin(i)=Gra_enth(atom_no,v,T_edge(i))
         Hmax(i)=Gra_enth(atom_no,v,T_edge(i+1))
         deltaH(i) = Hmax(i) - Hmin(i)
       enddo
      endif

      return
     
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Function integrates Silicate Heat Capacity to get Enthalpy (H-n)
c Heat capacity from experimental values (see G. Mark Voit 1991)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function Sil_enth(n,v,t)
      
      implicit none
      double precision n,v,t
       
      double precision f, fint1, fint2, fint3, fint4
      double precision c1, c2, c3, t2
       
      f = (1.d0 - (2.d0/n))*v
       
      c1 = 5.83333333333333333d+7
      c2 = 8.76008564039770429d+8
      c3 = 8.54862846766113158d+9
       
      fint1 = 0.d0
      fint2 = 0.d0
      fint3 = 0.d0
      fint4 = 0.d0
       
      if (t.le.50.d0) then
           fint1 = 1.40d3*(t*t*t)/3.d0
      else if (t.le.150.d0) then
           fint1 = c1
      	   t2    = 50.d0
       	   fint2 = 2.1647d4*((t**2.3d0) - (t2**2.3d0))/2.3d0
      else if (t.le.500.d0) then
       	   fint1 = c1
       	   fint2 = c2
       	   t2    = 150.d0
      	   fint3 = 4.8369d5*((t**1.68d0) - (t2**1.68d0))/1.68d0
      else
       	   fint1 = c1
    	   fint2 = c2
       	   fint3 = c3
       	   t2    = 500.d0
       	   fint4 = 3.3103d7*(t - t2)
      endif
       
      Sil_enth =f*( fint1 + fint2 + fint3 + fint4)
       
      return
       
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c Function integrates Graphite Heat Capacity to get Enthalpy (H-n)
c Heat capacity fitted (to within 3%) to experimental values (se G&D)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function Gra_enth(n,v,t)
       
      implicit none
      double precision n,v,t
      double precision Hatom

       
      Hatom=(4.15d-22*T**3.3)/
     &  (1 + 6.51d-3*T + 1.5d-6*T**2 + 8.3d-7*T**2.3)
       
      Gra_enth = (1.d0-2.d0/n)* n * Hatom

      return
       
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Function evaluates Planck spectrum given Temp and photon energy
c  ergs s-1 cm-2 Hz-1 Sr-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      double precision function Planck(invkT,E)
c invkT = 1/(kT)
c      
c      include 'const.inc'
c      double precision invkT, E, a, b
c
cc      Planck=(2.d0*E**3.d0/((cls**2)*(plk**3)))/(dexp(E/(rkb*T))-1.d0)
c      a = dexp(E*invkT)-1.d0
c      b = kbb*E*E*E*2
c      Planck=b/a
c      return
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Function evaluates Maxwell-Boltzmann distribution given 
c  electron Temp and energy
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function MBdist(T,E)
      
      include 'const.inc'
      double precision T,E, a, b, c, invpi
      
      parameter (invpi = 0.318309886183791d0)

c      MBdist=2.d0*dsqrt(E/(pi*(rkb*T)**3.d0))*dexp(-E/(rkb*T))
      a = 1.d0/(rkb*T)
      b = dsqrt(E*(a*a*a)*invpi)
      c = dexp(-E*a)
      MBdist=2*b*c

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine Makes Transisition matrix
c  |1 2 0 0 ... 0|      
c  |3 1 2 0 ... 0|
c  |.   .       :|
c  |:     .     :|
c  |: ....3 1 2 0|
c  |3 ......3 1 2|
c  |3 ......3 3 1|
c  2 is Thermal continuous cooling approx. only from level above
c  3 is heating as in Draine & Li
c  1 is the sum of the column representing losses from that level
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine transmatrix(grarea, n_e, s_f, Nmax, absmax, dabsmax, 
     &     t_e)
      
      include 'cblocks.inc'

      double precision grarea
      double precision n_e,s_f,stick,t_e

      integer*4 Nmax
      integer*4 absmax,dabsmax
      
      double precision den,wid
      double precision cool,heat
      double precision W1,W2,W3,W4,W5
      double precision Gfac,invfac
      double precision dht1, invdht1
      double precision Hmn, Hmx

      integer*4 i,j,k

c
c  Initialise
c
      do j=1,Nmax
       do i=1,Nmax
        Tr_matrix(i,j) = 0.d0
       enddo
      enddo
c
c  Cooling (#2)
c  Uses midpoint integration with Qabs energy points
c  Will improve later
c
      do i=2,Nmax   
       cool =0.d0
       do k=1,absmax
        den = 0.5d0*(ephot(k+1)+ephot(k))*eV
        wid = (ephot(k+1)-ephot(k))*eV
        if(den.gt.H_grid(i)) goto 501
        cool = cool+BBemiss(k,i)*Cabs(k)*wid
       enddo
 501   Tr_matrix(i-1,i)=(4.d0*pi)*invH(i,i-1)*cool
      enddo
      
      stick = n_e*grarea*s_f*(segrain*dexp(-t_e/2.0d5))

c      stick=0.d0
      
      do i=2,Nmax-1
c
c  Heating (#3) from l->u (j->i)
c  Allows for finite width of H bins
c  Uses midpoint integration with Qabs energy points
c

      Hmn = Hmin(i)
      Hmx = Hmax(i)

       do j=1,i-1
        W1=Hmn-Hmax(j)
        W2=MIN((Hmn-Hmin(j)),(Hmx-Hmax(j)))
        W3=MAX((Hmn-Hmin(j)),(Hmx-Hmax(j)))
        W4=Hmx-Hmin(j)
        invfac = 1.d0/(deltaH(i)*deltaH(j))
C
c  Photon & Collisional (currently turned off) heating
c
        heat=0.d0
        do k=1,dabsmax
         if(dengy(k).ge.W1) then
           if(dengy(k).lt.W2) then
             Gfac=(dengy(k)-W1)*invfac
           elseif(dengy(k).lt.W3) then
             Gfac=MIN(deltaH(i),deltaH(j))*invfac
           elseif(dengy(k).le.W4) then
             Gfac=(W4-dengy(k))*invfac
           else
             goto 512
           endif
           heat=heat+Gfac*(dabsphot(k)+stick*collheat(k)) 
         endif
        enddo  
 512    continue
        Tr_matrix(i,j)=deltaH(i)*invH(i,j)*heat   
       enddo
c  Allow for intrabin absorbtions
       heat=0.d0
       dht1 = deltaH(i-1)
      invdht1 = 1.d0/dht1
       do k=1,dabsmax
        if(dengy(k).gt.dht1) goto 521
         heat=heat+(1.d0-(dengy(k)*invdht1))*dabsphot(k)
       enddo
 521   continue
       Tr_matrix(i,i-1)=Tr_matrix(i,i-1)+invH(i,i-1)*heat
      enddo
      
      Hmn = Hmin(Nmax)
      do j=1,Nmax-1
       W1=Hmn-Hmax(j)
       W4=Hmn-Hmin(j)
       W5=1.d0/(W4-W1)
       heat=0.d0
c  Photon & collisional heating
       do k=1,dabsmax
        if((dengy(k).gt.W1).AND.(dengy(k).lt.W4)) then
          heat=heat+(dengy(k)-W1)*W5*(dabsphot(k)
     &             +stick*collheat(k))
        elseif(dengy(k).ge.W4) then
          heat=heat+dabsphot(k)+stick*collheat(k)
        endif
       enddo
       Tr_matrix(Nmax,j)=1.d0*invH(Nmax,j)*heat
      enddo
      heat=0.d0
      dht1 = deltaH(Nmax-1)
      invdht1 = 1.d0/dht1
      do k=1,dabsmax
       if(dengy(k).le.dht1) then
        heat=heat+(1.d0-dengy(k)*invdht1)*dabsphot(k)
       else
        goto 531
       endif
      enddo
 531  continue
      Tr_matrix(Nmax,Nmax-1)=Tr_matrix(Nmax,Nmax-1)+invH(Nmax,Nmax-1)
     &                * heat
c
c  Losses from i (Heating and Cooling) (#1)
c
c      open(69,file='Trmatrix.dat')
c      do i=1,Nmax
c       sum=0.d0
c       do j=1,Nmax
c         sum=sum+Tr_matrix(j,i)  !Note: Tr_matrix(i,i)=0 so j!=i is included
c       enddo
C       Tr_matrix(i,i)=-sum
c       write(69,169) (Tr_matrix(j,i),j=1,Nmax)
c      enddo
c 169  format(400(1pg11.4,x))
c      write(69,*)
c      close(69)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This subroutine solves for probability distribution given
c    Transition matrix Tr_matrix(i,j) (see G & D)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine probsolve(Tr_matrix,T_prob,Nmax,MxBin)

      implicit none      
      integer*4 Nmax,submax,MxBin
      parameter (submax=401)      !No. of Temp. bins in subroutine
      
      double precision Tr_matrix(MxBin,MxBin),T_prob(MxBin)
      double precision Bm(submax,submax), Xo(submax)
      double precision sum
      integer*4 i,j
       
      do i=1,Nmax
       Bm(Nmax,i)=Tr_matrix(Nmax,i)
       do j=Nmax-1,i,-1
         Bm(j,i)=Tr_matrix(j,i)+Bm(j+1,i)
       enddo
      enddo

      Xo(1)=1.d0

      do j=2,Nmax
       sum=0.d0
       do i=1,j-1
        sum=sum+Bm(j,i)*Xo(i)
       enddo
       if(Tr_matrix(j-1,j).ne.0.d0) then
        Xo(j)=sum/Tr_matrix(j-1,j)
       else
        print*, 'T(j-1,j)=0. IN trouble! j=',j
        stop
       endif
      enddo
      
      sum=0.d0
      do j=1,Nmax
       sum=sum+Xo(j)
      enddo
      if(sum.eq.0.d0) then
        print*, 'sum of X=0. IN trouble!'
        stop
      endif
c solve for T_prob
      do j=1,Nmax
       T_prob(j)=Xo(j)/sum
       if(T_prob(j).lt.1.0d-20) T_prob(j)=0.d0
      enddo

      return
      end
