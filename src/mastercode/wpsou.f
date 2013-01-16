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
c     standard routine to write out a source vector, given the 
c     field in tphot.
c
c     Vector Units: Fnu  = ergs/s/cm2/Hz/sr
c
c     External files are all 3D flux units 1/(4pi) because they are
c     generally written from diffuse field vectors.  2D (1/pi) units
c     are only used internally for source vectors and are converted
c     when files are read in in photsou.f
c
c     wmod = 'REAL' use total emission from slab, don't use dr
c     wmod = 'NORM' then normalise to 1cm slab, use dr
c     wmod = 'NFNU' then write nu v nuFnu (ergs/s/cm2/sr)
c
c     Now puts out two columns in the files: energy in eV and flux.
c
c     RSS 12/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wpsou(caller,pfx,np,wmod,
     &                 t,de,dh,
     &                 vl,dv,vl0,
     &                 di,dr,
     &                 eltim,tst,
     &                 scale,tp)
c
      include 'cblocks.inc'
c
      double precision tp(mxinfph),tl(mxinfph),t,de,dh,dr
      double precision bv,di,eltim,tst,vl,dv,vl0,scale
      double precision dhpr
      double precision blum, ilum, wid
      double precision q0, q1, q2, q3
c
      integer*4 lups,i,j,np
c
      character caller*4,wmod*4,tab*4
      character pfx*8,sfx*4,fps*20
c
c
      tab = char(9)
      lups = 44
c
c     
c     total energy in Inu
c     
c     
      blum = 0.d0
      ilum = 0.d0
c     
c     cvt to 1/4pi units in tl and sum
c
      do i = 1,infph-1
          tl(i) = 0.d0
          if (tp(i).ge.epsilon) then
              tl(i) = scale*tp(i)
           endif
c
         wid = ephot(i+1)-ephot(i)
         blum = blum + tl(i)*wid*evplk

         if (ephot(i).ge.Ryd) ilum = ilum + tl(i)*wid*evplk

      enddo
c
      blum = 4.d0*pi*blum
      ilum = 4.d0*pi*ilum
c
c
c    ***INTEGRATE NUMBER OF EMEGENT PHOTONS TO IONISE H,HE (=PI*JNU)
c    Assuming plane parallel geometry and dilution of 0.5 and intvec
c    assumes 1/4pi units in vector.
c
      q1   = 0.d0
      q2   = 0.d0
      q3   = 0.d0
      q0   = 0.d0
c
      call intvec(tl, q1, q2, q3,q0)
c
      if (wmod.eq.'NFNU') then
c
      fps = ' '
      sfx = 'nfn'
      call newfile(pfx,np,sfx,3,fps)
c
      else
c
c     output upstream field photon source file
c
      fps = ' '
      sfx = 'sou'
      call newfile(pfx,np,sfx,3,fps)
c
      endif
c
      open(lups, file=fps, status='NEW')
      i = infph-1
c
      if (wmod.eq.'NFNU') then
         write(lups, 862) dhpr,runname,wmod
 862     format('%NUFNU SPECTRUM '/
     &        '% DIFFUSE FIELD PLUS SOURCE, DHPR',1pg9.2/
     &        '% RUN:',a64/
     &        '% MODE:',a4/
     &        '% UNITS: nu v nuFnu (Hz v ergs/s/cm2/sr)'/
     &        '% Fnu in 3D units (1/(4pi)) sr not '/
     &        '% 2D source (1/pi) units. '/
     &        '% bin energies are lower edge of bins.'/
     &        '% fluxes are average over bin.')
      else
         write(lups, 860) dhpr,runname,wmod
 860     format('% PHOTON SOURCE FILE '/
     &        '% DIFFUSE FIELD PLUS SOURCE, DHPR',1pg9.2/
     &        '% RUN: ',a64/
     &        '% MODE:',a4/
     &        '% UNITS:eV vs Fnu (ergs/s/cm2/Hz/sr)'/
     &        '% Fnu in 3D units (1/(4pi)) sr not '/
     &        '% 2D source (1/pi) units. '/
     &        '% bin energies are lower edge of bins.'/
     &        '% fluxes are average over bin.')
      endif
c
 302  format('%:::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '% TOTAL LOSS:',1pg12.5,' % EFF. LOSS  :',1pg12.5/
     &     '% EFF.  GAIN:',1pg12.5,' % FRAC. RESID:',1pg12.5/
     &     '% TEMP.     :',1pg12.5,' % H DENSITY  :',1pg12.5/
     &     '% EL. DENS  :',1pg12.5,' % FR. NEUT. H:',1pg12.5/
     &     '% SHOCK VEL.:',1pg12.5,' % FLOW VELOC.:',1pg12.5/
     &     '% DELTA VEL :',1pg12.5,' % DISTANCE.  :',1pg12.5/
     &     '% SLAB DEPTH:',1pg12.5,' % ELAPSE TIME:',1pg12.5/
     &     '% TIME STEP :',1pg12.5,' % SOURCE MOD.:',A12    /
     &     '% DILUT. F. :',1pg12.5,' % SRC. TEMP. :',1pg12.5/
     &     '% ALPHA     :',1pg12.5,' % CUT-OFF    :',1pg12.5)
c     
      write(lups,302) 
     &     tloss,eloss,egain,dlos,t,dh,de,pop(1,1),
     &     vl0,vl,dv,di,dr,eltim,tst,
     &     iso,wdil,teff,alnth,cut
c     
 303  format('%:::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '%: Total intensity   : ',1pg12.5,' (ergs/s/cm^2)      :'/
     &     '%: Ion.  intensity   : ',1pg12.5,' (ergs/s/cm^2)      :'/
     &     '%: FQHI  (1-1.8Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     '%: FQHeI (1.8-4Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     '%: FQHeII   (>4Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     '%: FQ tot   (>1Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     '%:::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c     
      write(lups, 303) blum,ilum,q1,q2,q3,q0
c     
      write (lups, *) 'Produced by ',caller,' :MAPPINGS III v'
     &                ,theVersion
      write (lups, *)  FieldVersion
      write (lups, *)  infph-1
c     
c     **MULTIPLY BY TWO TO CHANGE DIFF. FIELD IN INTENSITY
c     
c     IMPORTANT UPDATE this is only true for split up/down diffuse fields.
c     
c     Assuming diffuse field is split into two equal upstream and 
c     downstream components - ie in shocks.  All routines will now send
c     scaling factor in scale parameter so that outward only and single slab
c     routines don't double the flux by mistake.  Units will always be in /sr
c     (1/4pi) units.
c     
c     
      do j = 1, infph-1
c     
         bv = 0.5d0*(ephot(j)+ephot(j+1))
c     
         if ((tp(j).ge.epsilon).and.(dr.gt.0.d0)) then
            if (wmod.eq.'NORM') tl(j) = (scale*tp(j)/dr)
            if (wmod.eq.'REAL') tl(j) = (scale*tp(j))
            if (wmod.eq.'NFNU') tl(j) = (scale*bv*evplk*tp(j))
         else
            tl(j) = 0.d0
         endif
c     
         if (wmod.eq.'NFNU') then
 865        format(1pg14.7,' ',1pg14.7)
            if (tl(j).gt.(2.d0*epsilon)) then
               write(lups,865) ephot(j)*evplk,tl(j)
            else
               write(lups,865) ephot(j)*evplk,epsilon
            endif
         else
 864        format(1pg14.7,' ',1pg14.7)
            if (tl(j).gt.(2.d0*epsilon)) then
               write(lups,864) ephot(j),tl(j)
            else
               write(lups,864) ephot(j),epsilon
            endif
         endif
c     
      enddo
c     
      close(lups)
c     
      return
c
      end
