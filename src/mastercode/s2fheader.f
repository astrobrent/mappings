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
      subroutine s2fheader(fnam,filna,filnb,filnc,filn,
     &luop,lupf,lusl,luions,dhpr,xhpr, tepr, tepo,tmag, banfil)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision vshockm
      double precision dhpr, xhpr, tepr,tepo, tmag
c
      integer*4 i,ie,j
      integer*4 luop,lupf,lusl,luions(4)
c
      character filna*20, filnb*20,filnc*20,filn(4)*20
      character banfil*20,fnam*20
      character fn*20,pfx*8,sfx*4,tab*4
c
c
c    ***WRITE INITIAL PARAMETERS IN OUTPUT FILE
c
c
c    ***DERIVE FILENAME : shckn**.sh2
c
      fn = ' '
      pfx = 'shckn'
      sfx = 'sh2'
      call newfile(pfx,5,sfx,3,fn)
      filna = fn
      fnam = filna
c
c
      fn = ' '
      pfx = 'shapn'
      sfx = 'sh2'
      call newfile(pfx,5,sfx,3,fn)
      filnb = fn
c
c
      fn = ' '
      pfx = 'allion'
      sfx = 'sh2'
      call newfile(pfx,6,sfx,3,fn)
      filnc = fn
c
c
      do i= 1, ieln
      if (i.eq.1) ie = iel1
      if (i.eq.2) ie = iel2
      if (i.eq.3) ie = iel3
      if (i.eq.4) ie = iel4
      fn = ' '
      pfx = elem(ie)
      sfx = 'sh2'
      if ((mapz(ie).eq.1).or.(mapz(ie).eq.6).or.
     &(mapz(ie).eq.7).or.(mapz(ie).eq.8).or.
     &(mapz(ie).eq.16)) then
         call newfile(pfx,1,sfx,3,fn)
      else
         call newfile(pfx,2,sfx,3,fn)
      endif
      filn(i) = fn
      enddo
c
c
      open(luop, file=filna, status='NEW')
      open(lupf, file=filnb, status='NEW')
      open(lusl, file=filnc, status='NEW')
c
      do i = 1,ieln
         open(luions(i),file=filn(i),status='NEW')
      enddo
c
c     main file
c
      write(luop, 2000) fnam, runname, banfil
 2000 format(' ',t40,'PLANE PARALLEL SHOCK MODEL: S2',3x,
     &'(DIFFUSE FIELD INCLUDED)',t100,'FILE  ',a16/t40,
     &'RUN:',a64,t100,'INPUT ',a12//)
c
c     write main header
c
      write(luop, 2020) dhpr, xhpr, tepr, tmag
 2020 format(/' PRE-SHOCK CONDITIONS :',t28,
     &'NUMBER DENSITY OF HYDROGEN :',1pg11.4,'/CC',t78,
     &'FRACTIONAL IONISATION :',1pg10.3/' --------------------',t33,
     &'PRE-SHOCK TEMPERATURE :',1pg13.6,'K',t74,
     &'TRANSVERSE MAGNETIC FIELD :',1pg9.2,' MICROGAUSS'/)
      write(luop, 2043) zgas
 2043 format(1h ,t9,'Abundances of the elements relative',
     &' to hydrogen :','  (Zgas=',f7.4,' Zsun) '/
     &' and the proto-ionisation used.')
c
      call wionabal(luop,propop)
c
 2807 format(/t9,' PHOTOIONISATION SOURCE AT THE SHOCK FRONT :')
      write(luop, 2800) 
 2800 format(/' MOD',t7,'TEMP.',t16,'ALPHA',t22,'TURN-ON',t30,'CUT-OFF'
     &,t38,'ZSTAR',t47,'FQHI',t56,'FQHEI',t66,'FQHEII',t76,'DILUF')

      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii, diluf
c

 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,4(1pg10.3))

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
 2003 format(' ',t40,
     &'PLANE PARALLEL SHOCK MODEL S2 (DIFFUSE FIELD INCLUDED)',10x,
     &'REF.FILE : ',a16/'RUN: ',a64/
     &t40,'--------------------------'/)
c
      write(lupf, 2020) dhpr, xhpr, tepr, tmag
      write(lupf, 2100) tepo, vshockm
      write(lupf, 2043) zgas
      call wionabal(lupf,propop)
c
c     allion file
c
      write(lusl, 2003) fnam,runname
      write(lusl, 2020) dhpr, xhpr, tepr, tmag
      write(lusl, 2100) tepo, vshockm
      write(lusl, 2043) zgas
      call wionabal(lusl,propop)
c
      tab  = char(9)
c
      write(lusl,*)'It#',tab,'distance',tab,'dr',tab,'Te',      
     &tab,'de',tab,'dh',tab,'nt',tab,'Tloss',tab,'Nloss',tab,'density',
     &tab,'pressure',tab,'internal energy',tab,'sound speed',tab,
     &'mol. wgt'
c
c     
      do i = 1,ieln 
c
      write(luions(i), 2003) fnam,runname
      write(luions(i), 2020) dhpr, xhpr, tepr, tmag
      write(luions(i), 2100) tepo, vshockm
      write(luions(i), 2043) zgas
c
 3000 format(2A4,29A6)

      if (i.eq.1) then
      write(luions(i),3000)'It#','Te ',(rom(j), j = 1,maxion(iel1))
      endif
      if (i.eq.2) then
      write(luions(i),3000)'It#','Te ',(rom(j), j = 1,maxion(iel2))
      endif
      if (i.eq.3) then
      write(luions(i),3000)'It#','Te ',(rom(j), j = 1,maxion(iel3))
      endif
      if (i.eq.4) then
      write(luions(i),3000)'It#','Te ',(rom(j), j = 1,maxion(iel4))
      endif
c
      enddo
c
c
c
      close(luop) 
      close(lupf) 
      close(lusl) 
c
      do i = 1,ieln
         close(luions(i))
      enddo
c
      return
c
      end





