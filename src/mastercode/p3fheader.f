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
c    Subroutine to get the file header buisness out where
c    it can be worked on, and/or modified for other headers.
c
c    first used for photo3 and photo4
c
c    RSS 8/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine p3fheader(newfil,filna,banfil,fnam,filnb,luop,
     &                     lupf,lusl,ia,dhn,fin)

c
      include 'cblocks.inc'
c
c           Variables
c
      double precision dhn,fin
c
      integer*4 i,j,luop,lupf,lups,lusl
      integer*4 ifil,ia(6)
c
      logical iexi
c
      character carac*20
      character newfil*20, filna*20, banfil*20, fnam*20, filnb*20
c
      ifil = 0
 1737 ifil = ifil+1
      filna(1:17) = ' '
      newfil(1:17) = ' '
      write(filna, 1827) ifil
      write(newfil, 1827) ifil
 1827 format(i4.4)
      filna = 'photn' // filna
      newfil = 'allion' // newfil
      filna(10:17) = '.ph3'
      inquire(file=filna, exist=iexi) 
      if (iexi) goto 1737
      inquire(file=newfil, exist=iexi)
      if (iexi) goto 1737
      fnam = filna(1:16)
      filnb = filna
      filnb(3:5) = 'apn'
      open(luop, file=filna, status='NEW', err=1737) 
      open(lupf, file=filnb, status='NEW') 
      open(lusl, file=newfil, status='NEW') 
c
c    ***WRITES ON FILE INITIAL PARAMETERS
c
c     
c     ***WRITES ON FILE INITIAL PARAMETERS
c     
      write(luop, 410) runname,banfil,fnam
      write(lupf, 410) runname,banfil,fnam
      write(lups, 410) runname,banfil,fnam
      write(lusl, 410) runname,banfil,fnam
      
 410  format(' Photoionisation model P3',
     &' (diffuse field calculated ,outward integration)'/ 
     &' ---------------------',
     &' (on-the-spot approximation, when source turned off)'//
     &' Run   : ',a64/
     &' Input : ',a12/
     &' Output: ',a13/)
      
      write(lusl, 2222) (zion(i),i = 1, atypes)
 2222 format(' ',11(1pg10.3))
c
      write(lupf, 412) 
  412 format(' This file contains the ionisation fraction for ',
     &'various atomic elements as a function of step number,'/
     &' and the pre-ionisation conditions used'/ )
      carac(1:20) = 'I  II IIIIV V  VI '
      write(lupf, 1220) 
     &((elem(ia(i)), rom(j),j = 1, maxion(ia(i))),i = 1, 6)
c
 1220 format(/'  #',t5,26(a2,a3))
      if (jeq .eq. 'E') then
      write(luop, 625) dtau0
  625 format(/' ',t5,'Summary ; Thermal and ionic equilibrium',5x,
     &'dtau :',0pf7.4)
      else if (jeq .eq. 'F') then
      write(luop, 627) telap, tlife, dtau0
  627 format(/' ',t5,'Summary  ;  Age :',1pg10.3,4x,'Source Life :'
     &,1pg10.3,' sec',4x,'dtau :',0pf7.4)
      else
      write(luop, 623) telap, dtau0
  623 format(/' ',t5,'Summary  ;  Equil. plus Switching off for :'
     &,1pg10.3,' sec',4x,'dtau :',0pf7.4)
      end if
c
      write(luop, 620) 
  620 format(' ',t5,'Rsou.',t15,'Remp.',t25,'Rmax',t35,'Diluav',t45,
     &'Fill.F.',t54,'Hdens',t64,'QHDHin',t74,'QHDHav')
      write(luop, 615) rstar, remp, rmax, wdpl, fin, dhn, qhdin
     &, qhdav
  615 format(1x,4(1pg10.3),0pf9.6,3(1pg10.3))
      write(luop, 422) 
  422 format(//' ',t2,'Jden',t8,'Jgeo',t14,'Jend',t20,'Ielen',t26,
     &'Jpoen',t35,'Fren',t41,'Tend',t47,'DIend',t57,'TAUen',t67,'Jeq'
     &,t74,'Teini')
      write(luop, 420) jden, jgeo, jend, ielen, jpoen, fren, 
     &tend, diend, tauen, jeq, tm00
  420 format('  '  
     &,3(a4,2x),2(i2,4x),0pf6.4,0pf6.0,2(1pg10.3),3x,a4,0pf7.1)
      write(luop, 430) 
  430 format(//' Description of Photon Source :')
      write(luop, 2800) 
 2800 format(/' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii
c
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
      write(luop, 435) 
  435 format(//' Description of each space step as computation',
     &' proceeds :')
      write(luop, 440) 
  440 format(/' #',t7,'Teav',t15,'DLOSav',t24,'DISout',t33,'DISav',t42,
     &'dr',t51,'Fill.F.',t60,'H.Den.',t69,'El.Den.',t78,'FHI',t87,'FHII'
     &,t96,'ZETAef',t105,'QHDH')
      if (jeq .eq. 'F') write(luop, 467) 
  467 format(t73,'Fthick',t82,'Treach',t91,'DTCO',t100,'Dton',t109,
     &'DToff')
      if (jeq .eq. 'P') write(luop, 421) 
  421 format(t100,'RECscal',t109,'DToff')
      write(*, 447) fnam
  447 format(/' Output in file : ',a,'~~~~~~~~~~~~~~~~~~~~~~~'/' #'
     &,t7,'Teav',t15,'DLOSav',t24,'DTon',t33,'DToff',t42,'DISav',t51,
     &'dr',t60,'Hden',t69,'FHI')
      close(luop) 
      close(lupf) 
      close(lusl) 
c
      return
c
      end
