cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Subroutine to get the file header buisness out where
c     it can be worked on, and/or modified for other headers.
c     
c     first used for photo3 and photo4
c     
c     RSS 8/90
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine p4fheader(newfil,banfil,fnam,
     &filna, filnb, filnc, filnd, filn,
     &luop, lupf, lusl, lups, lusp, luions,dhn,fin)
c     
      include 'cblocks.inc'
c     
c           Variables
c
      double precision dhn,fin,fl
      double precision dht
      integer*4 i,ie,j
      integer*4 luop,lups,lupf,lusl,lusp,luions(4)
      character fn*20,pfx*8,sfx*4
      character newfil*20, filna*20, banfil*20, fnam*20, filnb*20
      character filnc*20, filnd*20, filn(4)*20
c      
      double precision densnum
c
      dht = densnum(dhn)
c
      fn = ' '
      pfx = 'photn'
      sfx = 'ph4'
      call newfile(pfx,5,sfx,3,fn)
      filna = fn
      fnam = filna
c     
c     
      fn = ' '
      pfx = 'phapn'
      sfx = 'ph4'
      call newfile(pfx,5,sfx,3,fn)
      filnb = fn
c
      fn = ' '
      pfx = 'phlss'
      sfx = 'ph4'
      call newfile(pfx,5,sfx,3,fn)
      filnc = fn
c     
      fn = ' '
      pfx = 'spec'
      sfx = 'ph4'
      call newfile(pfx,4,sfx,3,fn)
      filnd = fn
c     
      if (jall.eq.'YES') then
         fn = ' '
         pfx = 'allion'
         sfx = 'ph4'
         call newfile(pfx,6,sfx,3,fn)
         newfil = fn
      endif
c     
      open(luop, file=filna, status='NEW') 
      open(lupf, file=filnb, status='NEW') 
      open(lups, file=filnc, status='NEW') 
      open(lusp, file=filnd, status='NEW') 
      if (jall.eq.'YES') open(lusl, file=newfil, status='NEW') 
c     
c     ***WRITES ON FILE INITIAL PARAMETERS
c     
      write(luop, 410) runname,banfil,fnam
      write(lupf, 410) runname,banfil,fnam
      write(lups, 410) runname,banfil,fnam
      write(lusp, 410) runname,banfil,fnam
      
      if (jall.eq.'YES') then
         write(lusl, 410) runname,banfil,fnam
      endif
      
 410  format(' Photoionisation model P4',
     &' (diffuse field calculated ,outward integration)'/ 
     &' ---------------------',
     &' (on-the-spot approximation, when source turned off)'//
     &' Run   : ',a64/
     &' Input : ',a12/
     &' Output: ',a13/)
      
      close(lusp) 

      write(lupf, 412) 
      
 412  format(' This file contains the plasma properties for ',
     &'each model step,'/
     &' and the pre-ionisation conditions used'/ )
 413  format(' This file contains the ionisation fraction for ',
     &'all atomic elements as a function of distance,'/
     &' and the abundances used'// )
      if (jall.eq.'YES') then
         write(lusl, 413) 
         write(lusl, 2222) (elem(i),zion(i),i = 1, atypes)
 2222    format('ABUND:',16(a2,':',1pg10.3,1x)//)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     photn file header
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jeq .eq. 'E') then
      
         write(luop, 623) dtau0         
 623     format(/' Summary : Thermal and Ionic Equilibrium :'/
     &           ' -----------------------------------------'/
     &           '  dTau        :',0pf7.4)
     
      else if (jeq .eq. 'F') then
      
         write(luop, 624) telap, tlife, dtau0
 624     format(/' Summary : Non-Equilibrium + Finite Source Life :'/
     &           ' ------------------------------------------------'/
     &           '  Age         :',1pg10.3,' sec'/
     &           '  Source Life :',1pg10.3,' sec'/
     &           '  dTau        :',0pf7.4)

      else
      
         write(luop, 625) telap, dtau0
 625     format(/' Summary : Equilibrium + Source Switch off :'/
     &           ' -------------------------------------------'/
     &           '  Switch off  :',1pg10.3,' sec'/
     &           '  dTau        :',0pf7.4)

      end if
c
 620  format(/' Radiation Field :'/
     &        ' -----------------'/
     &        ' Rsou.',t15,' Remp.',t25,' Rmax',t35,' DILUav'/
     &         1pg10.3,t15,1pg10.3,t25,1pg10.3,t35,1pg10.3//
     &        ' Hdens',t15,' Ndens',t25,' F.F.',t35,' QHDNin',t45,
     &        ' QHDNav'/
     &         1pg10.3,t15,1pg10.3,t25,0pf9.6,t35,1pg10.3,t45,1pg10.3)
c
      write(luop, 620) rstar, remp, rmax, wdpl,
     &                 dhn, dht, fin, qhdin , qhdav
c
c
c 
      write(luop, 422) 
 422  format(//' Model Parameters:'/
     &         ' -----------------'/
     &         ' ',t2,'Jden',t8,'Jgeo',t14,'Jend',t20,'Ielen',t26,
     &         'Jpoen',t35,'Fren',t41,'Tend',t47,'DIend',t57,
     &         'TAUen',t67,'Jeq',t74,'Teini')
c     
      write(luop, 420) jden, jgeo, jend, ielen,jpoen, fren,
     &                 tend, diend, tauen, jeq, tm00
 420  format('  ',3(a4,2x),2(i2,4x),0pf6.4,
     &            0pf6.0,2(1pg10.3),3x,a4,0pf7.1)
c
c
c
      write(luop, 430) 
 430  format(//' Description of Photon Source :'/
     &         ' ------------------------------')
c
      write(luop, 2800) 
 2800 format(/' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
c
c
      write(luop, 435) 
 435  format(//' Description of each space step:'/
     &         ' -------------------------------')
      write(luop, 440) 
      write(luop, 1440) 
 440  format(/' #',t7,'Teav',t15,'DLOSav',t24,'DISout',t33,'DISav',t42,
     &     'dr',t51,'Fill.F.',t60,'H.Den.',t69,'El.Den.')
 1440 format(t7,'FHI',t15,'FHII',t24,'ZETAef',t33,'QHDN',t42,'[OIII]'
     &     ,t51,'[OII]',t60,'[NII]',t69,'[SII]',t78,'Caseab H',
     &     ' Caseab He')
      if (jeq .eq. 'F') write(luop, 467) 
 467  format(t73,'Fthick',t82,'Treach',t91,'DTCO',t100,'Dton',t109,
     &'DToff')
      if (jeq .eq. 'P') write(luop, 421) 
 421  format(t100,'RECscal',t109,'DToff')
      write(*, 447) fnam
 447  format(//' Output in file : ',a//
     &' #',t7,'Teav',t15,'DLOSav',t24,'DISav',t33,'dr',t42,
     &'Hden',t51,'FHI',t60,'QHDN',t69,'[OIII]')
      close(luop) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phapn file header
c
c     new format: contains plasma properties in line, for 
c     easy plotting later.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 100  format(//'Dist.',t14,'Te',t26,'de',t38,'dh',t50,
     &'en',t62,'FHI',t74,'Density',t86,'Pressure',t98,
     &'Flow',t110,'Ram Press.',t122,'Sound Spd',t134,
     &'Slab depth',t146,'Alfen Spd.',t158,
     &'Mag. Field',t170,'Grain Pot.'//)
      write(lupf,100)
c
      close(lupf) 
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phlss file header
c
c     contains plasma cooling/heating in line, for 
c     easy plotting later.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 200  format(//'Dist.',t14,'Te',t26,'de',t38,'dh',t50,'en',t62,
     &     'Dloss',t74,'Eloss',t86,'Egain',t98,'Coll. H',t110,
     &     'Chrg. Ex.',t122,'Reson',t134,
     &     'Xreson',t146,'Inter/fine',t158,'He int-forb',t170,
     &     'Forbid',t182,'Fe II',t194,'2Photon',t206,
     &     'Compton',t218,'Free-Free',t230,'Coll. Ion',t242,
     &     'Photo.',t254,'Recomb.',t266,'Cosmic',t278,'Grains'//)
      write(lups,200)
c
      close(lups) 
c
      if (jall.eq.'YES') close(lusl) 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     monitor ions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c
      if (jiel.eq.'YES') then
c
c     write ion balance files if requested
c
         do i= 1, ieln
c     
            if (i.eq.1) ie = iel1
            if (i.eq.2) ie = iel2
            if (i.eq.3) ie = iel3
            if (i.eq.4) ie = iel4
            fn = ' '
            pfx = elem(ie)
            sfx = 'ph4'
            if ((mapz(ie).eq.1).or.(mapz(ie).eq.6).or.      
     &           (mapz(ie).eq.7).or.(mapz(ie).eq.8).or.
     &           (mapz(ie).eq.16)) then
               call newfile(pfx,1,sfx,3,fn)
            else
               call newfile(pfx,2,sfx,3,fn)
            endif
            filn(i) = fn
         enddo
c     
         do i = 1,ieln
            open(luions(i),file=filn(i),status='NEW')
         enddo
c     
c     
         do i = 1,ieln 
c     
            write(luions(i), *) fl,runname
c     
 3010       format(A4,29A6)
            
            if (i.eq.1) then
               write(luions(i),3010) 'Te ',(rom(j),j=1,maxion(iel1))
            endif
            if (i.eq.2) then
               write(luions(i),3010) 'Te ',(rom(j),j=1,maxion(iel2))
            endif
            if (i.eq.3) then
               write(luions(i),3010) 'Te ',(rom(j),j=1,maxion(iel3))
            endif
            if (i.eq.4) then
               write(luions(i),3010) 'Te ',(rom(j),j=1,maxion(iel4))
            endif
c     
         enddo
c     
         do i=1,ieln
            close(luions(i))
         enddo
c     
      endif

c
      return
c     
      end
