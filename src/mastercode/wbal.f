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
c     writes out a machine readable ionisation balance.
c
c     wmod ignored at this stage
c
c     rss 1992
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine wbal(caller,pfx,np,wmod,
     &t,de,dh,
     &vl,dv,vl0,
     &di,dr,
     &eltim,tst)     
c
      include 'cblocks.inc'
c
      double precision t,de,dh,dr
      double precision di,eltim,tst,vl,dv,vl0
c
      character caller*4,wmod*4,tab*4
      character pfx*8,sfx*4,fn*20
      integer*4 lups,i,np,nentries,j
c
      tab = char(9)
      lups = 44
c
c     write out balance file
c
      fn = ' '
      sfx = 'bln'
      call newfile(pfx,np,sfx,3,fn)
c
      open(lups, file=fn, status='NEW')
c
c
c
      write(lups,'("%")')
      write(lups,'("% Ionisation Balance File")')
      write(lups,'("%")')
      write(lups,'("%")')
c
      write(lups, 300) runname,wmod
 300  format('% RUN:',a64,'MODE:',a4) 
      write(lups,'("%")')
c
 302   format('%::::::::::::::::::::::::::'/
     &15h% TOTAL LOSS  :,1pg11.3/15h% EFF. LOSS   :,1pg11.3/
     &15h% EFF.  GAIN  :,1pg11.3/15h% FRAC. RESID.:,1pg11.3/
     &15h% TEMP.       :,1pg11.4/15h% H DENSITY   :,1pg11.3/
     &15h% ELECTR. DENS:,1pg11.3/15h% FR. NEUT. H :,1pg11.3/
     &15h% SHOCK VELOC.:,1pg11.3/15h% FLOW VELOC. :,1pg11.3/
     &15h% DELTA VELOC :,1pg11.3/15h% DISTANCE.   :,1pg11.3/
     &15h% SLAB DEPTH  :,1pg11.3/15h% ELAPSED TIME:,1pg11.3/
     &15h% TIME STEP   :,1pg11.3/15h% SOURCE MOD. :,A11    /
     &15h% DILUT. F.   :,1pg11.3/15h% SRC. TEMP.  :,1pg11.3/
     &15h% ALPHA       :,1pg11.3/15h% CUT-OFF     :,1pg11.3/
     &'%::::::::::::::::::::::::::')
c
       write(lups,302) 
     & tloss,eloss,egain,dlos,t,dh,de,pop(1,1),
     & vl0,vl,dv,di,dr,eltim,tst,
     & iso,wdil,teff,alnth,cut
c
       write(lups,*) 'Ionisation balance by ',caller,
     & ' :MAPPINGS III v',theVersion
c
       nentries = 0
c
       do i = 1 , atypes
          do j = 1,maxion(i)
             nentries = nentries+1
          enddo
       enddo
c
       write(lups,*) nentries
c
      do i = 1 , atypes
         do j = 1,maxion(i)
 445        format(i2,1x,i2,1x,1pg14.6)
            write(lups,445) mapz(i),j,pop(j,i)
         enddo
      enddo
      close(lups)
c
      return
c
      end
