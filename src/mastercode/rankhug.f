cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      rankhug solves the Rankine Hugoniot conditions between
c      any two points with full specification of v0, and 
c      loss rate, returns v1.
c
c      uses te0 and te1 :pre and post temps
c      magnetic field Hmag (Gauss), populations in common pop.
c      also need dh0, H num density
c
c     RSS 06/91
c     RSS 11/92
c     RSS 7/94  Fixed imaginary velcities
c     modified RSS 10.2009
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine rankhug(tpr,delpr,dhpr,vpr,hmag,rmod,tl,tstep)
c
      include 'cblocks.inc'
c
c
c
c
      double precision a(5), root
      double precision tpr,delpr,dhpr,vpr,hmag,humag
      double precision tl, tstep,c0,c1,c2
      double precision G,lambda, en0,Delta
      double precision frho,rv0,beta
      double precision cmpf,x1,x2,r1,r2
c      double precision x,y
c      integer*4 i
      character rmod*4
c
c     set globals, may be needed for iterations
c
      vel0 = vpr
      te0 = tpr
      dh0 = dhpr
      de0 = delpr
c
c
c
      if (vel0.lt.0.d0) then
        if (runmode.ne.'batchrun') then
         write(*,*) 'ERROR in RankHug: vel0 <= 0'
         write(*,*) vel0
        endif
        stop
      endif
c
c     adiabatic index gamma = 5/3 for monatomic gas
c
      gamma = 5.d0/3.d0
c
      G = (2.d0*gamma)/(gamma-1.d0)
c
c     ideally:    G = 5
c
c     get pressure and density terms
c     
      en0 = zen*dh0+de0
c
      pr0 = en0*rkb*te0
c
c     Hmag not used yet, but the mag field is calced anyway
c
      rho0 = frho(de0,dh0)
      humag = (hmag*hmag)/(8.d0*pi)
      hmago = humag/rho0
c
c     energy loss over step (assumed constant tl)
c
      Lambda = 2.d0*tl*tstep/rho0
c
c     compile quadratic terms
c
      c2 = G-1.d0
      c1 = -1.d0*G*(pr0/(rho0*vel0)+vel0)
      c0 = (vel0*vel0)+(G*(pr0/rho0))-Lambda
c
      Delta = (c1*c1)-(4.d0*(c2*c0))
c
      if (Delta .lt.0.d0) then 
         if (hmag.eq.0.d0) then
            write(*,*) 'r0,p0,t0,v0 :',rho0,pr0,te0,vel0
            write(*,*) 'c0,c1,c2,Delta :',c0,c1,c2,Delta
            write(*,*) 'Loss rate, timestep, lam:',tl,tstep,lambda
            stop
         else
c            write(*,*) 'r0,p0,t0,v0 :',rho0,pr0,te0,vel0
c            write(*,*) 'c0,c1,c2,Delta :',c0,c1,c2,Delta
c            write(*,*) 'Loss rate, timestep, lam:',tl,tstep,lambda
c
c     Use previous vel if delta  lt 0
c
            vel1 = vel0
            r1 = vel0
            r2 = vel0
c
         endif
      else

         r1 = ((-1.d0*c1)-dsqrt(Delta))/(2.d0*c2)
         r2 = ((-1.d0*c1)+dsqrt(Delta))/(2.d0*c2) 
c
c         write(*,*) 'Roots ',r1,r2
c
         if (rmod.eq.'FAR') then 
            vel1 = r2
         else
            vel1 = r1
         endif

      endif
c
      if (hmag.gt.0.d0) then
c
c	use non-magnetic velocity as initial guess for
c	quartic mag field solution
c
           beta = humag*4.d0
           rv0 = rho0*vel0
           a(5) = (G-1.d0)
           a(4) = -1.d0*G*(pr0/rv0+vel0+(beta/(4.d0*rv0)))
           a(3) = ((vel0*vel0)+(G*(pr0/rho0))-Lambda+beta)
           a(2) = G*(beta*vel0/(4.d0*rho0))
           a(1) = -1.d0*beta*vel0*vel0
c
c           do i = 1,20
c
c              x = r2*0.05*i
c              y = a(1) +x*(a(2)+x*(a(3)+x*(a(4)+x*a(5))))
c
c              write(*,'("x,y ",2(1pg12.5,1x))') x,y
c
c           enddo
c
           x1 = r1
           x2 = r2
           
           call quart(a,x1,x2,root,rmod)
c
c         write(*,*) 'MRoots ',r1,r2,root
c
c
c          if ((root/r1).lt.0.01) then
c             root = vel0
c          endif
c
          vel1 = root
c
      endif	 
c
      cmpf = vel0/vel1
      rho1 = rho0*cmpf
c
      call invrho(rho1,de1,dh1,pop)
c
      pr1 = pr0+(rho0*vel0*vel0*(1.d0-(vel1/vel0)))
     &      +humag*(1.d0-(vel0*vel0)/(vel1*vel1))
c
      te1 = pr1/((zen*dh1+de1)*rkb)
      bm0 = hmag
      hmag = hmag*cmpf
      bm1 = hmag
c
      return 
      end

