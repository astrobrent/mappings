cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                       Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     heating due to PAHS  : photo-electron ejection.
c     
c     assumes tphot (J_nu) is up to date.
c     
c     output in paheat (erg s^-1)
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine hpahs(t,de,dh)
c     
      include 'cblocks.inc'
c
      double precision t,de,dh
      double precision den,wid,pahrad,sigma
      double precision phi,ryd4,B
      double precision phots,photsn(mxinfph),photsi(mxinfph)
      double precision Jpe(pahi),Jec(pahi)
      double precision elmax,se,eC,qtmfac
      double precision englow,engmax,engesc
      double precision pahstate(pahi+2,pahi+1),pahsolve(pahi+1)
      double precision pahheat,pahcool
c      double precision sum,beta
      logical pahbin(mxinfph)
c
      integer*4 i,j,Nc
c     

c     
      paheat = 0.d0
      paheng = 0.d0
c      sum=0.d0
c     
      if(pahmode) then
c
c
c     assuming coronene (Nc=24)
c
            Nc = 468
c
c  Use Coronene Nc=24, radius ~ 5.85 A, S.A = (2/pi)* pi*a**2 (2/pi) 
c   due to disc orientation sigma= 6.85E-15
c    from http://ois.nist.gov/pah/sp922_Detail.cfm?ID=362
c     new uses 10A Draine & Li Grain
c     have dropped 4pi from sigma as it seems to fit.
c
            sigma = (10.0)**2*1.d-16 !*4.d0*pi 
c
            pahrad=10.0D-8

c     
            ryd4 = 4.d0*Ryd
c     
            B=pahIP(1)
c
            se=0.5d0
c           
c
c  Maximum electon sticking rate
c
            elmax = de*se*dsqrt((8*rkb*t)/(pi*me))*sigma
c     
c     Note: e^2/C=3.42 eV where C=2a/pi = capacitance of disk
c     
            eC=(eesu**2*pi)/(2*pahrad)
            qtmfac=1.d0/(1.d0+(27d-08/pahrad)**0.75)
            phi=eC/(rkb*t)
            eC=eC/eV
c
            do j = 1,infph-1
c     
               pahbin(j)=.FALSE.
               if ((ephot(j).gt.B).and.(ephot(j).lt.ryd4)) then     
c                  
                 den = 0.5d0*(ephot(j)+ephot(j+1))*ev
                 wid = (ephot(j+1)-ephot(j))*evplk
    
                 phots=tphot(j)*wid*4.d0*pi/den
                 photsn(j)=pahnabs(j)*phots
                 photsi(j)=pahiabs(j)*phots
                 if ((photsn(j).gt.0).OR.(photsi(j).gt.0)) then
                  pahbin(j)=.TRUE.
                endif
               endif
             enddo
c
c  Calculate electon collision rate and photoionisation rates
c
c     Negative PAHs
c
             do i=1,2
               Jpe(i)=0.d0
c
c  Electron capture to that level 
               Jec(i+1)=elmax*dexp(pahion(i+1)*phi)
c
c  Photoionisation from that level
               do j=1,infph-1
                if(pahbin(j)) then
                if(ephot(j).gt.pahIP(i)) then
                  Jpe(i)=Jpe(i)+photsn(j)*pahyield(j,i)
                endif
                endif
               enddo
             enddo
c
c     Neutral PAHs
c
             i=3
               Jpe(i)=0.d0
c
c  Electron capture to that level 
               Jec(i+1)=elmax*(1.d0+(pahion(i+1)*phi))
c
c  Photoionisation from that level
               do j=1,infph-1
                if(pahbin(j)) then
                if(ephot(j).gt.pahIP(i)) then
                  Jpe(i)=Jpe(i)+photsn(j)*pahyield(j,i)
                endif
                endif
               enddo
            
c
c     Positive PAHs
c
             do i=4,4
               Jpe(i)=0.d0
c
c  Electron capture to that level 
               Jec(i+1)=elmax*(1.d0+(pahion(i+1)*phi))
c
c  Photoionisation from that level
               do j=1,infph-1
                if(pahbin(j)) then
                if(ephot(j).gt.pahIP(i)) then
                  Jpe(i)=Jpe(i)+photsi(j)*pahyield(j,i)
                endif
                endif
               enddo
             enddo
c
c     Solve ionisation state of PAH
c               
             do i=2,6
               do j=1,7
                  pahstate(j,i)=0.d0
               enddo
             enddo

c     sum f(Z)=1.d0
             do i=1,6
               pahstate(i,1)=1.d0
             enddo
c     Jpe-Jec=0
             do i=2,5
               pahstate(i-1,i)=Jpe(i-1)
               pahstate(i,i)=-Jec(i)
             enddo
c solve
             call mdiag(5,pahstate,pahsolve)
c Calculate (+ve) ionisation fraction
             do i=1,pahi
                pahZ(i)=pahsolve(i)
             enddo

c     
c Calculate heating & cooling & energy absorbed
c
c Uses both Dopita Sutherland (2000) ApJ 539 742
c and Weingartner & Draine(2001) ApJS 134, 263
c engesc =average escape energy
c        = ef(E)/f(E) from W&D
c

             pahheat=0.d0
             pahcool=0.d0
C Negative PAHs
             do i=1,2
               englow=-(pahion(i)+1)*eC*qtmfac
               do j=1,infph-1
                den = 0.5d0*(ephot(j)+ephot(j+1))
                wid = (ephot(j+1)-ephot(j))*evplk
                paheng=paheng+pahZ(i)*pahnabs(j)*tphot(j)
     &                *(1.d0-pahyield(j,i))*wid*4.d0*pi
c                 sum=sum+pahZ(i)*pahnabs(j)*tphot(j)*wid*4.d0*pi
                if ((pahbin(j)).AND.(ephot(j).gt.pahIP(i))) then
                  engmax=den-pahIP(i)+englow
                  engesc=0.5d0*(engmax+englow)
                  pahheat=pahheat + pahZ(i)*photsn(j)*pahyield(j,i)
     &                              *engesc*eV
                endif
               enddo
               pahcool=pahcool + pahZ(i+1)*elmax*rkb*t
     &                 *(2-pahion(i+1)*phi)*dexp(pahion(i+1)*phi)
              enddo
C Neutral PAHs
             i=3
             englow=-(pahion(i)+1)*eC
             do j=1,infph-1
               den = 0.5d0*(ephot(j)+ephot(j+1))
               wid = (ephot(j+1)-ephot(j))*evplk
               paheng=paheng+pahZ(i)*pahnabs(j)*tphot(j)
     &               *(1.d0-pahyield(j,i))*wid*4.d0*pi
c                 sum=sum+pahZ(i)*pahnabs(j)*tphot(j)*wid*4.d0*pi
              if ((pahbin(j)).AND.(ephot(j).gt.pahIP(i))) then
                 engmax=den-pahIP(i)
c                beta=1.d0/(engmax-englow)**2                 
c                engesc=(engmax+englow-(4.d0*beta*englow**3)/3.d0)
c     &                  /(2.d0*(1.d0-2.d0*beta*englow**2))
                 engesc=0.5d0*engmax*(engmax-2*englow)/(engmax-3*englow)

                 pahheat=pahheat + pahZ(i)*photsn(j)*pahyield(j,i)
     &                              *engesc*eV
               endif
             enddo
             pahcool=pahcool + pahZ(i+1)*elmax*rkb*t
     &                 *(2+pahion(i+1)*phi)
              
C Positive PAHs
             i=4
             englow=-(pahion(i)+1)*eC
             do j=1,infph-1
              den = 0.5d0*(ephot(j)+ephot(j+1))
              wid = (ephot(j+1)-ephot(j))*evplk
              paheng=paheng+pahZ(i)*pahiabs(j)*tphot(j)
     &                *(1.d0-pahyield(j,i))*wid*4.d0*pi
c                sum=sum+pahZ(i)*pahiabs(j)*tphot(j)*wid*4.d0*pi
              if((pahbin(j)).AND.(ephot(j).gt.pahIP(i)))then
                engmax=den-pahIP(i)
c                beta=1.d0/(engmax-englow)**2                 
c                engesc=(engmax+englow-(4.d0*beta*englow**3)/3.d0)
c     &                /(2.d0*(1.d0-2.d0*beta*englow**2))
                engesc=0.5d0*engmax*(engmax-2*englow)/(engmax-3*englow)
 
               pahheat=pahheat + pahZ(i)*photsi(j)*pahyield(j,i)
     &                *engesc*eV
                
              endif
             enddo
             pahcool=pahcool + pahZ(i+1)*elmax*rkb*t
     &                 *(2+pahion(i+1)*phi)                
                   
             i=5
              do j=1,infph-1
                wid = (ephot(j+1)-ephot(j))*evplk
                paheng=paheng+pahZ(i)*pahiabs(j)*tphot(j)
     &                *(1.d0-pahyield(j,i))*wid*4.d0*pi
c                sum=sum+pahZ(i)*pahiabs(j)*tphot(j)*wid*4.d0*pi
              enddo
c
            paheat = (pahheat-pahcool)*pahfrac*dh
c
c  Testing grain charge - print output
c  
c
c 100        format(5(1pg11.4,1x))
c            write(*,*) 'PAHheat'
c            write(*,100) (pahZ(i),i=1,5)
c            write(*,100) (Jpe(i),i=1,5)
c            write(*,100) (Jec(i),i=1,5)            
c            write(*,100) paheng*dh*pahfrac,elmax,pahheat,pahcool,paheat
c           write(*,100) sum, paheng+pahheat, sum-(paheng+pahheat)

c     
      endif
c     
      return
c     
      end


