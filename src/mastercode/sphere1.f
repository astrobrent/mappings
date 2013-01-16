c     
      subroutine sphere1()
c     
c     
c     { }
c     { subroutine to comput 1D}
c     { isochoric spherical flow models}
c     { }
c     { with gravity , and atomic heating / cooling}
c     { }
c
c     Using a pre-determined density profile an pdV work
c     and given Mdot the velocity profile.
c     use dhy(500) for gram density here.
c     Version 1.0.0r
c     
      include 'cblocks.inc'
c     
c           Variables
c
      double precision bmag,dhn,dtime,gmd
      double precision rhozero
      double precision rstsun,tim0,tim1,ue,xhii
      double precision t, de, dh, G,  rad,  vtemp
      double precision lum,wid,blum,area
      double precision R0,R1,a,a2,dummyteff, gamma_dyn
      double precision dr,dv,dt,rs,as,q,q0,q1
      double precision en, pop1(mxion, mxelem)
      double precision fi,v02,v12,drfraction
      double precision ctime, rtime, ptime, ftime, stime, eqtime
      double precision sonic1, sonic2, sonic_ratio, rat1, rat2, hmc
     &                 , hmcdiff
c
      integer*4 luin, luout
      integer*4 count, step, i
c
      character where*16, title*80 ,nmod*4, wmod*4
      character pfx*8, sfx*4,fn*20, ispo*4, lmod*4, specmode*4
c     
      logical iexi
c
c           Functions
c
      double precision fdilu,fpressu,frho
      double precision fGetdhFromRadius, fGetVelFromRadius, feldens
      double precision fGetxhFromRadius, frectim2, fcolltim, fphotim
c     
c     
c     Formats
c     
c     
 1    format(a)
 95   format(12(1pg10.3,1x))
c 95   format(7(1pg10.3,1x))
c 95   format(7(1pg14.7,1x))
c     
c     init bookkeeping
c     
c     zero arrays and define modes: all ions, no on-thespot,
c     full continuum
c     calculations..
c     
      call zer
      call zerbuf
c
      wgain = 0.0d0
c
      do i = 1,500
         dist(i) = 0.0
         veloc(i) = 0.0
c rho in dhy cgs
         dhy(i) = 0.0
c pressure in xh cgs
         xh(i) = 0.0  
      enddo
c     
      limph = atypes
      jspot = 'NO'
      jcon = 'YES'
      jspec = 'NO'
      ispo = 'SII'
      specmode = 'DW'
      jgeo = 'S'
      jden = 'B'
      where = 'Sphere 1'
      luout = 11
c     
c     
c     Set up for MAPPINGS style output
c     
c     INTRO
c     
c 50   format(
c     &     ' SPHERE 1 model selected:'/
c     &     ' Isochoric Spherical flow in gravitational potential.'/
c     &     ' ****************************************************')
 50   format(
     &     ' ****************************************************'/
     &     ' Model xxx: Iteration xxx;  SPHERE1 output'/
     &     ' ****************************************************')
      if (runmode.ne.'batchrun') then
         write(*,50)
      endif
c
      pfx = 'sph'
      sfx = 'sp1'
      call newfile(pfx,3,sfx,3,fn)
c
      open(unit=luout,file=fn,status='NEW')
      write(luout,50)
      close(unit=luout)
c
c
      gamma = 5.0 / 3.0

c        gamma = 4.0 / 3.0

      G = (gamma/(gamma-1))
c     
c     get ionisation balance to work with...
c     
      where = 'pre-ionisation'
      call popcha(where)
      where = 'Sphere 1'
c     
c     set up current field conditions
c     
      call photsou(where)
c
c
      blum = 0.d0
c     
      do i = 1,infph-1
         wid = ephot(i+1)-ephot(i)
         blum = blum + soupho(i)*wid*ev/plk
      enddo
c     
      luin = 10
 420  title = ' '
      if (runmode.ne.'batchrun') then
         write (*,fmt=430)
 430     format(/' Enter profile file name : ',$)
      endif
      read(*,440) title
 440  format(a)
      inquire(file=title, exist=iexi)

      if (iexi) then 
         open(unit=luin, file=title, status='OLD')
         read(luin,'(a80)') title
         write(*,*) title
         read(luin,*) lum,mdot,dummyteff
         read(luin,*) rs,as
         write(*,'(2(1pg14.7,1x))')lum, mdot,rs,as,dummyteff
         read(luin,*) nsteps
c     
         if (nsteps.lt.2) then
            write(*,*) 'Error:Less than two profiles points...'
            stop
         endif
c     
c     Read dhy as gm/cm^3 here, not h denisty
c
         do i = 1,nsteps
            read(luin,*) dist(i),dhy(i),veloc(i),xh(i)
         enddo
c     
         close(luin)
      else
         write(*,*) 'File: not found',title
         goto 420
      endif
      
c     
c
      if (lum.eq.0.d0) then
         rstar = 0.d0
      else
         rstar = 1/(2*pi)*dsqrt(lum/blum)
c     
      write(*,*) ' ****************************************'
      write(*,*) ' Eff. Star Radius : ',rstar
      write(*,*) ' ****************************************'
c
      endif
c     
      rstsun = rstar/6.96d10
c     
c     Neutron star 1.4 M_0
c     
c      GM = 6.67259e-8 * 1.4 * 1.991e33      
c
 55   format(/'Give the potential mass (M_sun)',$)
      write(*,55)
      read(*,*) GMD
      GM = 6.67259e-8 * GMD * 1.991e33


 57   format(/'Give the gas index used by dynamical flow program',$)
      write(*,57)
      read(*,*) gamma_dyn

c
c     Sphere flow model preferences
c     
 60   format(///
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' Setting the flow conditions'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      if (runmode.ne.'batchrun') then
         write(*,60)
      endif
c     
 70   format(//' Give Initial Radius and step size:',/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '  R_0(cm) (Radius <100 as a log)'/
     &     '  delta R (<0.1)'/
     &     ' :: ',$)
c     
 200  if (runmode.ne.'batchrun') then
         write(*,70)
      endif
      read(*, *) r0,drfraction
c     
      if (r0.le.100.d0) r0 = 10.d0**r0
      if (r0.lt.1.d6) goto 200
      if (drfraction.gt.0.1) drfraction = 0.1
      if (drfraction.le.0.0) drfraction = 0.01
c     
      vel0 = fGetVelFromRadius(r0)
      rho0 = fGetdhFromRadius(r0)
      pr0 = fGetxhFromRadius(r0)
c
       write(*,*) rho0,vel0,pr0
c
      remp = r0
      rad0 = remp
      dis0 = rad0
c     
c     
c     ***DENSITY BEHAVIOR
c     
      jden = 'C'
c     
      call invrho(rho0,de0,dhn,pop)
      write(*,*) de0,dhn
      dh0 = dhn      
c     
      de0 = feldens(dh0,pop)
      rho0 = frho(de0,dh0)
c
c     Get H=C equilibrium temp for initial point
c
c     Ist Radiation:
c
      wdil = fdilu(rstar,r0)
c
c      write(*,*) ' Dilution Factor: ',wdil
c
      dr = drfraction*r0
      t = Pr0/((de0 + dh0*zen)*rkb)
      write(*,*) t,pr0,de0,dh0,zen,rkb
c
c     Hack...
c
      if (lum.eq.0.d0) then
c
 75      format(/'Give the intial gas temperature',$)
         write(*,75)
         read(*,*) t
c     
         if (t.lt.10.0) t = 10.d0**t
c     
      endif
c     
      te0 = t

c     dv = 0.0d0


      vtemp = fGetVelFromRadius(r0-dr)
      dv = vtemp - vel0


      fi = 1.d0
c
      write(*,*) 'Te',t
      write(*,*) 'vel0, vtemp, dv',vel0, vtemp, dv
c
c     Iterate on intial temp
c
 10   call localem(t,de0,dh0)
      lmod = 'ALL'
      call totphot(t, dh0, fi, rad0, dr, dv, wdil, lmod)
      call zetaeff(fi, dh0)
c     
      if (lum.eq.0.d0) then
         call equion(te0, de0, dh0)
      else
         nmod = 'EQUI'
         dt = 1.d33
         call teequi(t, te0, de0, dh0, dt, nmod)
      endif
c
      pr0 = fpressu(te0,dh0,pop)
      pr1 = pr0
c
      write(*,*) te0,dlos
c
      wmod = 'SCRN'
      call wmodel(0,te0,de0,dh0,dr,wmod)
c
      if (abs( (pr1-pr0)/pr1 ).gt.1.d-2) then
         pr1 = pr0
         t = te0
         goto 10
      endif
c
      area = 4*pi*r0*r0
c
      dtime = (dr/vel0)
c
      write(*,*) 'M/Msun =', GMD
      write(*,*) 'L (erg/s) =', lum
      write(*,*) 'Mdot (g/s) =', mdot
      write(*,*) 'Teff (K) =', dummyteff
      write(*,*) 'gamma_dyn =', gamma_dyn
      write(*,*) 'ndata = xxxxx'
      
      write(*,*) ' ****************************************'
c      Write(*,*) ' Initial Sonic Radius: ', rs
      Write(*,*) ' Initial "Sonic" Radius: ', rs
      Write(*,*) ' Initial H density: ', dh0
      Write(*,*) ' Initial electron density: ', de0
      Write(*,*) ' Initial Temperature : ',te0
      Write(*,*) ' Initial Pressure: ',pr0
      write(*,*) ' Initial Velocity: ',vel0
      write(*,*) ' Initial Step Size: ',dr
      write(*,*) ' Initial Dynamical Timescale: ',dtime
      write(*,*) ' ****************************************'
c
c
      open(unit=luout,file=fn,status='OLD',access='APPEND')
c
      write(luout,3456) 'M/Msun =', GMD
3456   format(a9, 1pe10.3)
      write(luout,3457) 'L (erg/s) =', lum
3457   format(a12, 1pe10.3)
      write(luout,3458) 'Mdot (g/s) =', mdot
3458   format(a13, 1pe10.3)
      write(luout,3459) 'Teff (K) =', dummyteff
3459   format(a11, 1pe10.3)
      write(luout,3460) 'gamma_dyn =', gamma_dyn
3460   format(a12, 1pe20.13)
      write(luout,3461) 'ndata = xxxxx'
3461   format(a14)
      write(luout,*) '*************************************************'
c      Write(luout,*) ' Initial Sonic Radius: ', rs
      Write(luout,*) ' Initial "Sonic" Radius: ', rs
      Write(luout,*) ' Initial H density: ', dh0
      Write(luout,*) ' Initial electron density: ', de0
      Write(luout,*) ' Initial Temperature : ',te0
      Write(luout,*) ' Initial Pressure: ',pr0
      write(luout,*) ' Initial Velocity: ',vel0
      write(luout,*) ' Initial Step Size: ',dr
      write(luout,*) ' Initial Dynamical Timescale: ',dtime
      write(luout,*) '*************************************************'
c
      close(luout)
c
c      rs = 1.d12
c
      if (vel0 .eq. 0.0)  vel0 = epsilon
c     
c     fill in other default parameters
c     
      Bmag = 0.d0
      vshoc = vel0
      rhozero = rho0
c     
c     total ions number
c     
      zen = 0.d0
      do i = 1,atypes
         zen = zen+zion(i)
      enddo
c     
      en = dh*zen
c     
c     prepare for time step in ions
c     
c     record initial conditions, t0 and l0 and pop
c     
      call copinto(pop,pop0)
c     
      count = 0
      step = 0

      tim0 = 0.d0

c     
c     Integration loop reentry point
c     
      write(*,*) ' R(cm)   Te(K)   Rho(g/cc)  v/a  (H-C)(erg/cc/s)', 
     &           '(H-C)/(H+C)  t_eq(s)  t_s (s)'

      write(*,*) ' R(cm)   dr(cm)  time(s)    dt(s)     de(/cc)',
     &           '   dh(/cc)  GM/(2 ra^2)  sonic_rat'
      write(*,*) '--------------------------------------------------',
     &  '----------' 
c
      open(unit=luout,file=fn,status='OLD',access='APPEND')
      write(luout,*) '  R(cm)      Te(K)    Rho(g/cc)     v/a',
     &  '   GM/(2ra^2)   son_rat  (H-C)(erg/cc/s) (H-C)/(H+C) t_eq(s)',
     &  '    ratGM    rata2   P/k'
      write(luout,*) '----------------------------------------------',
     &  '-----------------------------'
      close(luout)
c     
 150  step = step +1
c
c      if (r0.lt. 0.1 * rs) goto 1000
c      if (r0.lt. 0.001 * rs) goto 1000
      if (r0.lt. 0.01 * rs) goto 1000
c     
      call copinto(pop0,pop)
c     
      dr = drfraction*r0
c
c      tim0 = 0.d0
c
c      dv = 0.0d0

      vtemp = fGetVelFromRadius(r0-dr)
      dv = vtemp - vel0

c      write(*,*) 'vel0, vtemp, dv',vel0, vtemp, dv

      wdil = fdilu(rstar,r0)
      call localem(te0,de0,dh0)
      lmod = 'ALL'
      call totphot(te0, dh0, fi, r0, dr, dv, wdil, lmod)
      call zetaeff(fi, dh0)
c
      call cool(te0,de0,dh0)
c
      ue = 1.5*(de0 + dh0*zen)*te0*rkb
      dt = dabs(10*ue/tloss)
      dtime = (dr/vel0)
c
      ctime = fcolltim(de0)
      rtime = frectim2(de0)
      ptime = fphotim()
      
c      write(*,*) 'Ctime, Rtime, Ptime'
c      write(*,'(3(1pg12.5,1x))') Ctime, Rtime, Ptime
c
      eqtime = (1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
      eqtime = dabs(1.d0/eqtime)
c
c     experiment with timescales
c
      stime = (1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      stime = dabs(1.d0/stime)
c
      ftime = 1.d0/((1.d0/(stime*0.1))
     &             +(1.d0/(eqtime*0.1)))
c
c      write(*,*) dtime, dt, ftime
c
      dt = abs(1.0/((1.0/dt) + (1.0/dtime)+ (1.0/ftime)))
c
      dr = vel0*dt
c
c      write(*,*) 'dr  dt :',dr, dt
c
      a2 = gamma * pr0 / rho0
      a = sqrt(a2)

c       sonic_ratio = GM/(2.*a2*r0)

      sonic1 = GM/(2.*a2*r0)

      sonic2 = (GM/6.)*((egain-eloss)/Rho0)

      sonic2 = sonic2 / a**5

       sonic_ratio = sonic1*(1.-sonic2)

        hmc = (egain-eloss)/Rho0 

        rat1 = 0.67*(hmc*r0*r0)/(GM*vel0)

        rat2 = 0.33*(hmc*r0)/(a2*vel0)


        hmcdiff = ( (egain-eloss)/(egain+eloss) )

c
      write(*,95) r0, Te0, Rho0, vel0/a, (egain-eloss), hmcdiff,
     &            eqtime, stime
      write(*,95) r0, dr, tim0, dt, de0, dh0, sonic1, sonic_ratio
      open(unit=luout,file=fn,status='OLD',access='APPEND')
c      write(luout,95) r0, dr, tim0, dt, Te0, de0, dh0
      write(luout,95) r0, Te0, Rho0, vel0/a, sonic1, sonic_ratio, 
     &                (egain-eloss), hmcdiff, eqtime, rat1, rat2,
     &                   pr0/1.38e-16
      close(luout)
c
      r1 = r0 - dr
      tim1 = tim0+dt
c
c     conditions at end of step
c
      vel1 = fGetVelFromRadius(r1)
      rho1 = fGetdhFromRadius(r1)

      dv = vel1 - vel0

c
c      write(*,*) 'r1      vel1        rho1     dv'
c      write(*,'(4(1pg14.7,1x))') r1,vel1,rho1, dv
c      write(*,*) 'dh0         de0        pop(6,5)'
c      write(*,'(3(1pg14.7,1x))') dh0,de0,pop(6,5)
c
      call invrho(rho1,de1,dh1,pop)
c
c     inner neq convergence loop
c
      count = 0
      q1 = 0.d0
      q0 = 0.d0
c
 151  continue
c
      count = count+1
c
      call copinto(pop0,pop)
c
      t = (te0+te1)/2.d0
      dh = (dh0+dh1)/2.d0
      de = feldens(dh,pop)
      xhii = pop(2,1)


          rad = 0.5 * ( r0 + r1 )

c
      call timion(t,de,dh,xhii,dt)
c
      call copinto(pop,pop1)
c
c     after timion on the average t,de,dh, now
c     calculate new lossrate at te1,de1,dh1
c
c     calculate the radiation field and atomic rates
c
      call localem(t,de,dh)

c
c               write(*,*) 'before'
c               write(*,*) 'tphot'
c               write(*,*) tphot
c               write(*,*) 't,dh,fi,rad,dr,dv,wdil,specmode'
c               write(*,*) t,dh,fi,rad,dr,dv,wdil,specmode
c


      call totphot(t, dh, fi, rad, dr, dv, wdil, specmode)

c
c               write(*,*) 'after'
c               write(*,*) 'tphot'
c               write(*,*) tphot
c               write(*,*) 't,dh,fi,rad,dr,dv,wdil,specmode'
c               write(*,*) t,dh,fi,rad,dr,dv,wdil,specmode
c

      call zetaeff(fi, dh)
c
      dh = dh1
      de1 = feldens(dh,pop)
      de = de1
c
      call cool(t,de,dh)


c          write(*,*) 't, de, dhi, tloss, egain, eloss'
c          write(*,*) t, de, dhi, tloss, egain, eloss

c          write(*,*) 'vel0, vel1, rho0, rho1'
c          write(*,*) vel0, vel1, rho0, rho1

c
      v02 = vel0*vel0
      v12 = vel1*vel1
c
      Q = 2.d0*(egain-eloss)*dt/(rho0+rho1)
c
c      Q = 0.d0
c
      Pr1 = (rho1/G)*(
     &     (v02-v12)*0.5 
     &     + GM*(r0-r1)/(r0*r1) 
     &     + G*Pr0/rho0
     &     + Q)
c
c      write(*,*) (v02-v12)*0.5, GM*(r0-r1)/(r0*r1), G*Pr0/rho0, Q
c      
      te1 = Pr1 / ((de1 + dh1*zen) * rkb)
c
      write(*,'(5(1pg14.7,1x))') te1,pr1,q,egain,eloss,fflos,cmplos
c
      q0 = q1
      q1 = q

c
c     if we are flip-flopping then
c     if q1 is for lower temperature then force convergence to q1
c     otherwise convergence will fail and a new iteration will
c     occur until q1 <= q0.
c
c
      if (count.gt.10) then
         if (q1.le.q0) q0=q1
c
c   Taking higher temperature solution --- just a test.
c
c       if (q1.ge.q0) q0=q1
      endif
c

c      q = dabs(1.d0-(q1/q0))
      q = dabs( (q0-q1) / q0 )
c
      if ((count.lt.3).or.(q.gt.1.d-3)) goto 151
c      if ((count.lt.3).or.(q.gt.2.d-3)) goto 151
c
      pr0 = pr1
      r0 = r1
      tim0 = tim1
      vel0 = vel1
      rho0 = rho1
      de0 = de1
      dh0 = dh1
      te0 = te1
c     
      call copinto(pop1,pop0)
c

c    Monitor ion pop and losses/gains at each step.

c       if(r0.lt.1e13)then
c      if(r0.ge.1e8)then
c      open(unit=luout,file=fn,status='OLD',access='APPEND')
c      call wionpop(luout,pop0)
c       wmod = 'FILE'
c      call wmodel(luout,te0,de0,dh0,dr,wmod)
c      close(luout)
c      endif


      goto 150
c
 1000 continue
c
      return
c     
c     end sphere
c     
      end
