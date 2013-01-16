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
      subroutine s3fheader(fnam,filna,filnb,filnc,luop,lupf,lusl
     &                 ,dhpr, xhpr, tepr, tepo,tmag ,banfil)
c
      include 'cblocks.inc'
c
c
      double precision dhpr, xhpr, tepr, tepo,tmag, vshockm
c
      integer*4 luop,lupf,lusl
c
      character filna*20, filnb*20,filnc*20
      character banfil*20,fnam*20
      character fn*20,pfx*8,sfx*4
c
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
c
c    ***DERIVE FILENAME : shckn**.sh3
c
      fn = ' '
      pfx = 'shckn'
      sfx = 'sh3'
      call newfile(pfx,5,sfx,3,fn)
      filna = fn
      fnam = filna
c
c
      fn = ' '
      pfx = 'shapn'
      sfx = 'sh3'
      call newfile(pfx,5,sfx,3,fn)
      filnb = fn
c
c
      fn = ' '
      pfx = 'allion'
      sfx = 'sh3'
      call newfile(pfx,5,sfx,3,fn)
      filnc = fn
c
c
      open(luop, file=filna, status='NEW')
      open(lupf, file=filnb, status='NEW')
      open(lusl, file=filnc, status='NEW')
c
c     main file
c
      write(luop, 2000) fnam, runname, banfil
 2000 format('PLANE PARALLEL SHOCK MODEL: S3',3x,
     & '(separate el. and ion temps)',/,' FILE  ',a16,/,
     & ' RUN: ',a64,/,' INPUT: ',a20,//)
c
c     write main header
c
      write(luop, 2020) dhpr, xhpr, tepr, tmag
 2020 format(/' PRE-SHOCK CONDITIONS :',t28,
     & 'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     & 'FRACTIONAL IONISATION :',1pg10.3/' --------------------',t33,
     & 'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     & 'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write(luop, 2043) zgas
 2043 format(9x,'Abundances of the elements relative',
     &' to hydrogen :','  (Zgas=',f7.4,' Zsun) ',
     &' and the proto-ionisation used.')
c
      call wionabal(luop,propop)
c
 2807 format(/t9,' PHOTOIONISATION SOURCE AT THE SHOCK FRONT :')
      write(luop, 2800) 
 2800 format(/' MOD',t6,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t46,'FQHI',t55,'FQHEI',t65,'FQHEII',t75,'DILUF')

      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii, diluf
c

 2850 format(1x,a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))

      vshockm = vshoc/1.d5
      write(luop, 2100) tepo, vshockm
 2100 format(/' POST-SHOCK CONDITIONS :',t32,
     &'POST-SHOCK TEMPERATURE :',1pg13.6,'K',t76,
     &'COMPUTED SHOCK VELOCITY :',1pg12.5,'KM/SEC'/
     &' ---------------------')
c
c     apn file
c
      write(lupf, 2003) fnam,runname
 2003 format(
     &'PLANE PARALLEL SHOCK MODEL S3 (Separate el. and ion temps)',
     & /,' REF.FILE: ',a16,/,' RUN: ',a64,//)
c
      write(lupf, 2020) dhpr, xhpr, tepr, tmag
      write(lupf, 2043) zgas
      call wionabal(lupf,propop)
c
c     allion file
c
      write(lusl, 2003) fnam,runname
      write(lusl, 2020) dhpr, xhpr, tepr, tmag
      write(lusl, 2043) zgas
      call wionabal(lusl,propop)
c
c
c
c
      close(luop) 
      close(lupf) 
      close(lusl) 
c
      return
c
      end





