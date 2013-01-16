cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.2.0md
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Heating due to grains: grain photo-electron ejection.
c                          : +PAH contribution
c                          : +collisional grain cooling
c     
c     Assumes tphot (J_nu) is up to date.
c     
c     Output in gheat, gcool, paheat
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine hgrains(t,de,dh)
c     
      include 'cblocks.inc'
c     
      double precision t,de,dh
      double precision den,wid,phots(mxinfph),B,je,acc_H,yi
      double precision eu,jec,jpr,elmax,prmax,phi
      double precision eng,engmax,englow
      double precision je1,jec1,jpr1,rydmax,Qabs
      double precision pmax
      double precision gam,del,delmax
      double precision sigma,mred,dgrad
      double precision se,sp,spf
c      double precision total
      double precision U, U1, U2, Umid, y1, y2, ymid
c
      integer minbin,maxbin
      integer*4 i,j,nbins,nr,check,Usolv,k,dtype
c     
 100  format(5(1pg11.4,1x))
c     
c     Grain initialization
c
      do dtype=1,numtypes
         avgpot(dtype) = 0.d0
      enddo
      gheat = 0.d0
      gcool = 0.d0
c
c     Init neutral grains.
c
      if (grainmode .and. (Yinf.ne.0.d0)) then     

c     Some global variables are copied to local variables for readibilty.
c     Standard parameters: B (cutoff) and Yinf (yield factor),
c     (B in eV), Qabs as fraction of geometric size,
c     acc_H : accomodation factor.
        acc_H = haccom
        B = Bgrain
        yi = Yinf
c     Sticking factors for electrons (se) and protons (sp)
        se = segrain*dexp(-t/2.0d5)
        sp = spgrain    
c     Electron and proton currents, less phi factor.
c     These need to be multiplied by phi factor for different 
c     potentials.
        elmax = -de*se*dsqrt((8*rkb*t)/(pi*me))
        prmax = pop(2,1)*dh*sp*dsqrt((8*rkb*t)/(pi*mp))
c
c     Cooling collision rates. Separate from grains.
c
        delmax = de*dh*((2.d0*rkb*t)**(1.5d0))/dsqrt(pi*me)
        mred = dsqrt(me/mh)
c
c     Do for all types and grainsize distribution.
c
        do 1000 dtype = 1,numtypes
         do 2000 k = mindust(dtype),maxdust(dtype)
c    
          dgrad = gradedge(k+1)-gradedge(k)
          sigma = dustsig(k,dtype)*dgrad !sigma = grain area per h atom
          rydmax = Ryd
          nbins = 0
          pmax = 0.d0
c      total = 0.d0
c
c     Calculate phots(inl) as same thruout (only diff bins accessed)
c     and get totalphot and maxphot
c
       minbin=0
       maxbin=0
       do i = 1,infph-1

       if ((ephot(i).gt.B)) then
c
         Qabs = absorp(i,k,dtype)
c            
         if (Qabs.gt.0.d0) then
           if (minbin.eq.0) minbin=i
           den = 0.5d0*(ephot(i)+ephot(i+1))
           wid = (ephot(i+1)-ephot(i))*evplk
c
c     only look where absorption will occur
c                
           phots(i) = Qabs*tphot(i)*wid/(den*ev)
           phots(i) = Yi*(1.d0 - B/den)*phots(i)
c
           if (phots(i).gt.pmax) pmax = phots(i)
           maxbin=i  !sets max bin to last bin where Qabs > 0
c           total = total+ phots(i)
c
c     end Qabs>0
c
            endif
c     
         endif
      enddo
c
c     Now solve for charge balance to get grain potential
c     photo-electron rate:
c      
c     Use binary search method to get grain charge such that 
c     f(U)=e-sticking + p-sticking + photoelectric=0
c
c     Changed to allow for all photon energies up to 290 eV.

               check = 0
               nr = rydmax
               nr = 10*((nr/10)+1)
               U1 = -1.d0*dble(nr)
               U2 = 1.d0*dble(nr)

c     Get f(U1), check that its +ve, 0, or -ve.
c     If +ve continue with binary search. If -ve, double U1 and try again.
c     If 0 (ie U=U1), exit search.
c     (Note: f(U) only equals zero at one point)

 200           continue

               U = U1
               eu = ev*U1
               phi = eu/(rkb*t)
c     U is -ve 
               jec1 = elmax*dexp(phi)+epsilon
               jpr1 = prmax*(1.d0-phi)

c     Get photo currents.
c     Loop through whole vector in case vector changes one day.
c     Only want bins below 8 Ryd and above cutoff.

               je1 = 0.d0
               do j = minbin,maxbin-1
c     As U -ve, all photoelectrons escape (spf=1).
                 je1 = je1+phots(j)
               enddo


               y1 = je1+jpr1+jec1
c               write (69,*) 'U1'
c               write (69,100) U, je1, jpr1, jec1, y1 
               if ((y1 .lt. 0) .and. (check .lt. 5)) then
                  U1 = 2.d0*U1
                  check = check+1
                  go to 200
               else if (y1 .eq. 0) then
                  U2 = U1
                  goto 300
               else if (check .ge. 5) then
                  write (*,*) 'Grain Charge error'
                  write (*,100) U2,U1,je1,jpr1,jec1
                  write (*,100) elmax,prmax,phi
                  write (*,100) de,t,se
                  stop
               endif

c
c    Get f(U2), check that its +ve,0, or -ve.
c    If -ve, continue with binary search. If +ve, double U2 and try again.
c    If 0 (ie U=U2), exit search

 250           check = 0
               U = U2
               eu = ev*U2
               phi = eu/(rkb*t)
c  U +ve.
               jec1 = elmax*(1.d0+phi)
               jpr1 = prmax*dexp(-phi)+epsilon

               je1 = 0.d0

                do j = minbin,maxbin-1
                 if (ephot(j) .gt. (B+U)) then
                  eng = 0.5d0*(ephot(j)+ephot(j+1))
                  engmax = (eng-B)
                  englow = U
c   U +ve, so some fraction of photoelectrons can't escape potential.
                  if (englow .gt. (0.5*engmax)) then
                     spf = 2.d0*((engmax-englow)/engmax)**2
                  else       ! 0 < englow < engmax/2 
                     spf = 1.d0-2.d0*(englow/engmax)**2
                  endif
                  je1 = je1+phots(j)*spf
                 endif
                enddo
    
               y2 = je1+jpr1+jec1
c               write (69,*) 'U2'
c               write (69,100) U, je1, jpr1, jec1, y2 
               if ((y2 .gt. 0) .and. (check .lt. 5)) then
                  U2 = 2.d0*U2
                  check = check+1
                  goto 250
               else if (y2 .eq. 0) then
                  U1 = U2
                  goto 300
               else if (check .ge. 5) then
                  print*,'graincharge error2'
                  stop
               endif

c     Now that F(U1) > 0 > F(U2) (ie U1 < U < U2) use binary search to find U
c     to several decimal places (up to 1e-6*original Umid).
c

               do Usolv = 1,20
                  Umid = (U1+U2)/2.d0 
                  U = Umid
                  eu = ev*U
                  phi = eu/(rkb*t)
                  if (phi .ge. 0.d0) then
                     jec1 = elmax*(1.d0+phi)
                     jpr1 = prmax*dexp(-phi)+epsilon
                  else
                     jec1 = elmax*dexp(phi)+epsilon
                     jpr1 = prmax*(1.d0-phi)
                  endif
                  je1 = 0.d0
                   do j = minbin,maxbin-1
c     ephot > B due to Minbin.
                     if ((ephot(j) .gt. (B+U))) then
                      eng = 0.5d0*(ephot(j)+ephot(j+1))
                      engmax = (eng-B)
                      englow = U
                      if (englow .le. 0) then
                         spf = 1.d0
                      else if (englow .lt. (0.5*engmax)) then
                         spf = 1.d0-2.d0*(englow/engmax)**2
                      else
                         spf = 2.d0*((engmax-englow)/engmax)**2
                      endif
                      je1 = je1+phots(j)*spf
                     endif
                   enddo

                  ymid = je1+jpr1+jec1
c                  write (69,100) U, je1, jpr1, jec1, ymid 
                  if (ymid .gt. 0.d0) then      ! Too much PE heating
                     U1 = Umid
                  else if (ymid .lt. 0.d0) then ! Electron current too large
                     U2 = Umid
                  else
                     goto 300
                  endif
               enddo

 300           continue

c     Then grainpot(k,dtype) = U.

               U = (U1+U2)/2.d0

c     Have balance potential in U, get heating.

               eu = ev*U
               phi = eu/(rkb*t)

c               write (*,*) 'grains:',U,eu,phi

               if (phi .ge. 0.d0) then
                  jec = elmax*(1.d0+phi)
                  jpr = prmax*dexp(-phi)+epsilon
               else
                  jec = elmax*dexp(phi)+epsilon
                  jpr = prmax*(1.d0-phi)
               endif
c
c     Get heating term, (ergs not photons).
c
               je = 0.d0
               gam = 0.d0
               do j = minbin,maxbin-1
                if ((ephot(j) .gt. (B+U))) then
                   eng = 0.5d0*(ephot(j)+ephot(j+1))
                   engmax = (eng-(B+U))
                   englow = -U
c
c     Calculate average e- escape energy.
c     (Note: eng=eng/spf, phots=phots*spf but cancelled out here)
                   if (englow .ge. 0) then               ! U -ve
                      spf = 1.d0
                      eng = 0.5d0*(engmax + englow)
                   else if ((engmax+englow) .gt. 0) then ! engmax > abs(englow)
                      spf = 1.d0-2.d0*(englow/(engmax-englow))**2
                      eng = 0.5d0*(-4.d0*englow**3/
     &                     (3.d0*(engmax-englow)**2)+ engmax+englow)
                   else                                  ! engmax < abs(englow)
                      spf = 2.d0*(engmax/(engmax-englow))**2
                      eng = 2.d0*engmax**3/(3.d0*(engmax-englow)**2)
                   endif
                   gam = gam+eng*ev*phots(j)
                   je = je+phots(j)*spf
                endif
               enddo

               gam = gam*dh*sigma*4.d0

c     Collisional cooling.

c     Have phi so...

               del = delmax*sigma
               if (phi .ge. 0.d0) then
                  del = del*(se*(2.d0+phi)+mred*dexp(-phi)*
     &                 (acc_H+acc_H+phi))
               endif
               if (phi .lt. 0.d0) then
                  del = del*((se*(2.d0-phi)*dexp(phi))+
     &                 (mred*acc_H*(2.d0-phi)))
               endif

c     Grain size/type.

               grainpot(k,dtype) = U
               avgpot(dtype) = avgpot(dtype) + U
               gheat = gheat + gam
               gcool = gcool + del

 2000       continue            ! End loop over grain sizes

            avgpot(dtype) = avgpot(dtype)/float(maxdust(dtype)
     &           -mindust(dtype)+1)

c            write (*,100) U

 1000    continue               ! End loop over grain types
         
      endif
      
c      write (*,1069) gcool, gheat
c 1069 format ('graincool: ',1pg14.7,' grainheat: ',1pg14.7)
      
      return
c     
      end

