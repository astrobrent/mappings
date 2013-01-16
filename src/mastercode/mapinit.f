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
c     subroutine to initailise mappings for non interactive data
c     ie data files and so on....
c
c
c     RSS 7/90
c
c     if the file map.prefs is found in the current working directory
c     then it will be used in place of dats/ATDAT
c
c     RSS 12/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine mapinit(error)
c
c
      include 'cblocks.inc'
c
c
c           Variables
c

      double precision en
      integer*4 i,ion,j,k,kpre,line,luin,atom
      integer*4 luop,min,ni,npre,ntr,series,trans
c
      character ibell*4,ilgg*4, ibuf(11)*4
c
      logical error,iexi
c
c
c
      double precision nair,lam
c
c	Refractive index of std air for Lam>2000 A
c
      nair(lam) = 1.0d0
     &          +(643.26d0+294981.d0/(146.d0-lam*lam)
     &          + 2554.d0/(41.d0-lam*lam))/1.0d7
c
      error = .false.
c
      luin = 20
      luop = 23
c
      ibell(1:4) = char(7)
c
c
c    ***START READ-IN OF DATA
c
      if (runmode.ne.'batchrun') then
      write(*, 99)
 99   format(///'Start read-in of data :::::::::::::::::::::::'/
     &          ':::::::::::::::::::::::::::::::::::::::::::::'/)
      endif
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
 11   format(I3,1x,a4,1x,I4,1x,f9.5,1x,f8.2,1x,f6.3,1x,I2)
c
c
c
c     set up roman numerals array: 1-30
c
      rom(1) = 'I'
      rom(2) = 'II'
      rom(3) = 'III'
      rom(4) = 'IV'
      rom(5) = 'V'
      rom(6) = 'VI'
      rom(7) = 'VII'
      rom(8) = 'VIII'
      rom(9) = 'IX'
      rom(10) = 'X'
      rom(11) = 'XI'
      rom(12) = 'XII'
      rom(13) = 'XIII'
      rom(14) = 'XIV'
      rom(15) = 'XV'
      rom(16) = 'XVI'
      rom(17) = 'XVII'
      rom(18) = 'XVIII'
      rom(19) = 'XIX'
      rom(20) = 'XX'
      rom(21) = 'XXI'
      rom(22) = 'XXII'
      rom(23) = 'XXIII'
      rom(24) = 'XXIV'
      rom(25) = 'XXV'
      rom(26) = 'XXVI'
      rom(27) = 'XXVII'
      rom(28) = 'XXVIII'
      rom(29) = 'XXIX'
c
      do 100 i= 1,28
         zmap(i) = 0
         mapz(i) = 0
 100  continue
c
c
      inquire(file='map.prefs', exist=iexi)
      if (iexi) then
         open(luin, file='map.prefs', status='OLD')
      else
         open(luin, file='dats/ATDAT', status='OLD')
      endif
c
 110  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 110
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      zen = 0.d0
      read (luin, fmt = *) atypes
      do 120 i = 1,atypes
      read(luin,11) j,elem(i),k,
     &        atwei(i),zion0(i),dion(i),maxion(i)
      zen = zen+zion0(i)
      invdion(i)=1.d0/dion(i)
      zmap(k) = j
      mapz(j) = k
      if (i.ne.j) goto 1000
 120  continue
c
c     read collmode (0 = A&R,1 = Shull L&M)
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      read (luin,*) collmode
c
c     read chargetransfer reaction mode (0 = new,1 = old)
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      read (luin,*) chargemode
c
c      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      read (luin,*) photonmode
c
      photonmode = 1
      crate = 0.d0
c
      hellcoolmode   = 1
      alphacoolmode  = 0
c
      expertmode = 0
c
c     expert mode triggers many more detailed questions,
c     useful for testing and program design.  NOT DOCUMENTED
c     SO THERE.  Contact Ralph Sutherland (ralph@madras.anu.edu.au)
c     For details if you are programming and modifying mappings.
c
c     expertmode = 0 is default user mode.
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      read (luin,*) expertmode
c
      close(luin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Read ionisation and recombination data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readphotdat(luin,error)
      call readphiondat(luin,error)
c
      call readcolldat(luin,error)
c     
      call readiondat(luin,error)
      call readion2(luin,error)
c
      if (oprec) call readoprec(luin,error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Read continuum data 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readcont(luin,error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     LINES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Read X-ray/FUV line data
c
      call readxdata(luin,error)
c
c     Read other Line data
c
      call read2level(luin,error)
      call read9level(luin,error)
      call read6level(luin,error)
      call read5level(luin,error)
      call read3level(luin,error)
      call readfeIIlin(luin,error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Hydrogenic Data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readHHedata(luin,error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DUST
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readdust(luin,error)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Other data
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call readstardat(luin,error)
c
c     End readin
c
      if (runmode.ne.'batchrun') then
      write(*, 19)
 19   format(///'Data read sucessfully :::::::::::::::::::::::'/
     &          ':::::::::::::::::::::::::::::::::::::::::::::'/)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Tidy Up and Precalc some tables, assign bins to lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	Find bin to start phion loops
c
      en = ipotpho(1)
      do j = 1, infph - 1
      if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
                IonStartBin = j
      endif
      enddo
c
c      *FIND CORRESPONDING ENERGY BINS IN EPHOT(N) FOR RESONANCE
c      AND INTERCOMBINATION LINES
c
      do i = 1, xlines
         en = xejk(i)
         xbin(i) = 0
         do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               xbin(i) = j
            endif
         enddo
      enddo
c
c	Hydrogen
c
      do line = 1, 10
         do series = 1, 6
             en = LmeV/hlambda(line,series)
             hbin(line,series) = 0
             do j = 1, infph - 1
             if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
                hbin(line,series) = j
             endif
             enddo
         enddo
      enddo
c
c helium ii
c
      do line = 1, 10
         do series = 1, 6
             en = LmeV/helambda(line,series)
             hebin(line,series) = 0
             do j = 1, infph - 1
             if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               hebin(line,series) = j
             endif
             enddo
         enddo
      enddo
c
c	Heavy Hydrogenic ions
c
      do atom = 1, atypes
      do series = 1, 2
         do line = 1, 10
             xhlambda(line,series,atom) = hlambda(line,series)
     &                                    /(mapz(atom)*mapz(atom))
             en = LmeV/xhlambda(line,series,atom)
             xhbin(line,series,atom) = 0
             do j = 1, infph - 1
             if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
                xhbin(line,series,atom) = j
             endif
             enddo
         enddo
      enddo
      enddo
c
c	Three Level Ions
c
      do ion = 1, nf3ions
         do trans = 1, nf3trans
c             en = 1.985940d-16/f3lam(trans,ion)
             en = (plk*cls/ev)/f3lam(trans,ion)
             f3bin(trans,ion) = 0
             do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               f3bin(trans,ion) = j
            endif
         enddo
         enddo
      enddo
c
c	Five Level Ions
c
      do ion = 1, nfions
         do trans = 1, nftrans
c             en = 1.985940d-16/flam(trans,ion)
             en = (plk*cls/ev)/flam(trans,ion)
             fbin(trans,ion) = 0
             do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               fbin(trans,ion) = j
            endif
         enddo
         enddo
      enddo
c
c	Six Level Ions
c
      do ion = 1, n6ions
         do trans = 1, n6trans
c             en = 1.985940d-16/f6lam(trans,ion)
             en = (plk*cls/ev)/f6lam(trans,ion)
             f6bin(trans,ion) = 0
             do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               f6bin(trans,ion) = j
            endif
         enddo
         enddo
      enddo
c
c	Nine Level Ions
c
      do ion = 1, n9ions
         do trans = 1, n9trans
             en = (plk*cls/ev)/f9lam(trans,ion)
             f9bin(trans,ion) = 0
             do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               f9bin(trans,ion) = j
            endif
         enddo
         enddo
      enddo
c
      do i = 1, xilines
         en = xiejk(i)
         xibin(i) = 0
         do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               xibin(i) = j
            endif
         enddo
      enddo
c
c	He 1
c
      do i = 1, xhelines
         en = xhejk(i)
         xhebin(i) = 0
         do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               xhebin(i) = j
            endif
         enddo
      enddo
c
      do i = 1, nlines
         en = e12r(i) / ev
         lrbin(i) = 0
         do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               lrbin(i) = j
            endif
         enddo
      enddo
c
      do i = 1, nheilines
         en = (LmeV/(heilam(i)*1.d8))
         heibin(i) = 0
         do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               heibin(i) = j
            endif
         enddo
      enddo
c
      do i = 1, mlines
         en = e12fs(i) / ev
         lcbin(i) = 0
         do j = 1, infph - 1
            if ((ephot(j+1).gt.en).and.(ephot(j).le.en)) then
               lcbin(i) = j
            endif
         enddo
      enddo
c
      npre = 0
      kpre = -1
      do 460 i = 1, nlines
      min = infph
      do 470 j = 1, nlines
      if ((lrbin(j) .eq. npre).and.(j.gt.kpre)) goto 480
      if ((lrbin(j).ge.min).or.(lrbin(j).le.npre)) goto 470
      min = lrbin(j)
      k = j
 470  continue
c
      linpos(1,i) = min
      linpos(2,i) = k
      kpre = k
      npre = min

      goto 460
 480  linpos(1,i) = lrbin(j)
      linpos(2,i) = j
      kpre = j
      npre = lrbin(j)
 460  continue
c
c     put branching ratios in xbr, plus xdismul init
c
      do i = 1,xlines
         ni = xiso(i)
         ntr = idint(xtrans(i))
         if (dble(ntr).eq.xtrans(i)) then
            xbr(i) = brr(ni,ntr)
         else
            xbr(i) = 1.d0
         endif
      enddo
c
c      *FINDS IONISATIONS CROSS SECTION NUMBER FOR HEII
c
      do 490 i = 1, ionum
      jhe2p = i
      if ((atpho(i) .eq. 2).and.(ionpho(i).eq.2)) goto 500
 490  continue
 500  continue
c
c
c     calc secondary hydrogen data: anl, dnl 
c
c     As ususal nn,ll is upper level, n,l is lower
c
c     anl, total exit from level nn,ll
c
c      do nn = 1,20
c      do ll = 0,nn-1
c         j = 1+(nn*(nn-1)/2)+ll
c         anl(j) = 0.d0 
c         do n = 1,nn-1
c         do l = ll-1,ll+1
c            if ((l.ge.0).and.(l.ne.n).and.(l.lt.n))then
c               i = 1+(n*(n-1)/2)+l
c               anl(j) = anl(j)+anlnl(i,j)
c            endif
c         enddo
c         enddo
c      enddo
c      enddo
c
c
c     dnl, for l mixing
c
c      do nn = 1,20
c      do ll = 0,nn-1
c         j = 1+(nn*(nn-1)/2)+ll
c         dnl(j) = 6*nn*nn*((nn*nn) - (ll*ll) - (ll) - 1)
c      enddo
c      enddo
c
c
      goto 1010
c     
 1000 error = .true.
c
 1010 continue
c
      return
      end
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
      subroutine readphotdat(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      integer*4 i,j,luin,nentries
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
 7    format(//'ERROR IN FILE : PHOTDAT @#',i3)
c
      error = .false.
c
      open(luin, file='dats/PHOTDAT', status='OLD')
 100  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 100
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
c
      FieldVersion = 0
      IonStartBin = 1
      read(luin, *) nentries
      if (nentries.lt.0) then 
          FieldVersion = nentries
          read(luin, *) nentries
      endif
c
      if (nentries.gt.mxinfph) then
          write(*,*) ' INCOMPATIBLE dats/PHOTDAT, need',
     &                 mxinfph,' bins or less.'
          error = .true.
          stop
      endif
c
      infph = nentries
c
      do j = 1, infph
           read(luin,*) ephot(j)
      enddo
c
      close(luin)
c
      do 110 j = 1, infph - 1
      if (ephot(j).lt.ephot(j + 1)) goto 110
      if (runmode.ne.'batchrun') then
          write(*, 7) j
      endif
      error = .true.
      stop
 110  continue
c
c     initialise skipbin
c
      do i = 1,infph
         skipbin(i)=.true.
      enddo
c
      return
      end
c
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
      subroutine readcolldat(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision aa,bb,cc,dd,el,po
      integer*4 i,j,luin,atom,ion,is,at,io
      integer*4 ni,nions,ns,nsequ,sh,shell
c
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
 6    format(a2,1x,a2,1x,a3,1x,3(i2,1x),g9.3,g9.2,1x,g9.2,g9.2,1x,g9.2)
c
      error = .false.
c
c
      open(luin, file='dats/COLLDAT', status='OLD')
 121  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 121
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read (luin, fmt = *) nsequ
c
c     read in each sequence
c
      do 122 i = 1,nsequ
         read(luin, fmt = 1) (ibuf(j),j = 1,11)
c         write(*, 3) (ibuf(j),j = 1, 11)
         read(luin, *) nions,ns
c
c     for each ion in the sequence
c
         do 123 j = 1 ,nions
            do 124 shell = 1, ns
               atom = 0
               ion = 0
               read (luin,6) el,is,sh,at,io,ni,po,aa,bb,cc,dd
               if (zmap(at).ne.0) then 
c
c     Check to see if the atom type is one chosen in ATDAT
c
                  atom = zmap(at)
                  if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c
                     ion = io
                     ipot(shell,ion,atom) = po
                     aar(shell,ion,atom) = aa
                     bar(shell,ion,atom) = bb
                     car(shell,ion,atom) = cc
                     dar(shell,ion,atom) = dd
                  endif
               endif
 124        continue
            if ((atom.ne.0).and.(ion.ne.0)) nshells(ion,atom) = ns
            if ((at.eq.11).and.(ion.eq.1)) nshells(ion,atom) = 2
 123     continue
 122  continue
      close(luin)
c     
      return
      end
c
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
      subroutine readiondat(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision lt(8)
      double precision alp1,alp2,ep,hhe
      integer*4 at,atom,i,ion,io,j,luin,m,nentries,nions
      double precision sg1,sg2,ve
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
c
      open(luin, file='dats/IONDAT', status='OLD')
 131  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 131
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
c
c     read primary ionisation potentials
c     and convert to ergs
c
      read(luin, fmt = *) nions
      do i = 1,nions
         read(luin, fmt = *) at,io,ep
         if (zmap(at).ne.0) then
c
c     Check to see if the atom type is one chosen in ATDAT
c
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               ion = io
               epot(ion,atom) = ep*ev
            endif
         endif
      enddo
c
c     read in misc ionisation rate stuff (see IONDAT)
c
c
c     H & He rates plus entries for on the spot approx
c
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         if (zmap(1).ne.0) read (luin,*) m,(arad(j,zmap(1)),j = 2,6)
         if (zmap(2).ne.0) read (luin,*) m,(arad(j,zmap(2)),j = 2,6)
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         if (zmap(1).ne.0) read (luin,*) m,(xrad(j,zmap(1)),j = 2,6)
         if (zmap(2).ne.0) read (luin,*) m,(xrad(j,zmap(2)),j = 2,6)
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         if (zmap(1).ne.0) read (luin,*) m,(adi(j,zmap(1)),j = 2,6)
         if (zmap(2).ne.0) read (luin,*) m,(adi(j,zmap(2)),j = 2,6)
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         if (zmap(1).ne.0) read (luin,*) m,(t0(j,zmap(1)),j = 2,6)
         if (zmap(2).ne.0) read (luin,*) m,(t0(j,zmap(2)),j = 2,6)
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         if (zmap(1).ne.0) read (luin,*) m,(bdi(j,zmap(1)),j = 2,6)
         if (zmap(2).ne.0) read (luin,*) m,(bdi(j,zmap(2)),j = 2,6)
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         if (zmap(1).ne.0) read (luin,*) m,(t1(j,zmap(1)),j = 2,6)
         if (zmap(2).ne.0) read (luin,*) m,(t1(j,zmap(2)),j = 2,6)
c
c     low T dielec correction
c
c
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) at ,(lt(j),j = 2,8)
            if (zmap(at).ne.0) then
               atom = zmap(at)
               do ion = 2,8
                  adilt(ion,atom) = lt(ion)
               enddo
            endif
         enddo
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) at ,(lt(j),j = 2,8)
            if (zmap(at).ne.0) then
               atom = zmap(at)
               do ion = 2,8
                  bdilt(ion,atom) = lt(ion)
               enddo
            endif
         enddo
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) at ,(lt(j),j = 2,8)
            if (zmap(at).ne.0) then
               atom = zmap(at)
               do ion = 2,8
                  cdilt(ion,atom) = lt(ion)
               enddo
            endif
         enddo
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) at ,(lt(j),j = 2,8)
            if (zmap(at).ne.0) then
               atom = zmap(at)
               do ion = 2,8
                  ddilt(ion,atom) = lt(ion)
               enddo
            endif
         enddo
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) at ,(lt(j),j = 2,8)
            if (zmap(at).ne.0) then
               atom = zmap(at)
               do ion = 2,8
                  fdilt(ion,atom) = lt(ion)
               enddo
            endif
         enddo
c
c     charge exchange and secondary ionisation data
c
         nchxold = 0
c
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) m,at,io,sg1,sg2,alp1,alp2,hhe
            if (zmap(at).ne.0) then
               atom = zmap(at)
               if (io.le.maxion(atom)) then
                  nchxold = nchxold+1
                  charco(1,nchxold) = dble(atom)
                  charco(2,nchxold) = dble(io)
                  charco(3,nchxold) = sg1
                  charco(4,nchxold) = sg2
                  charco(5,nchxold) = alp1
                  charco(6,nchxold) = alp2
                  charco(7,nchxold) = dble(hhe)
               endif
            endif
         enddo
c
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do i = 1 , nentries
            read (luin,*) m,at,io,ve
            if (zmap(at).ne.0) then
               atom = zmap(at)
               if (io.le.maxion(atom)) then
                  secat(i) = atom
                  secio(i) = io
                  secel(i) = ve
               endif
            endif
         enddo
c
c
         close(luin)
c
      return
      end
c
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
      subroutine readion2(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c

      double precision aa,ac,ad,ar,bd
      double precision ta,tc,ti,to,xr
      integer*4 at,atom,i,io,ion,j,luin,nentries
      integer*4 ni,nions,nsequ
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
c     
c
      open(luin, file='dats/IONDAT2', status='OLD')
 151  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 151
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read (luin, fmt = *) nsequ
c
c     read in each sequence
c
      do 152 i = 1,nsequ
         read(luin, fmt = 1) (ibuf(j),j = 1,11)
         read(luin, *) nions
c
c     for each ion in the sequence
c
         do 153 j = 1 ,nions
            atom = 0
            ion = 0
           read (luin,*) at,io,ni,ac,tc,ar,xr,ad,bd,to,ti,aa,ta
           if (zmap(at).ne.0) then
c     
c     Check to see if the atom type is one chosen in ATDAT
c
               atom = zmap(at)
               if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c     
                  ion = io
                  acol(ion,atom) = ac
                  tcol(ion,atom) = tc
c
                  ion = io+1
                  arad(ion,atom) = ar
                  xrad(ion,atom) = xr
c
c     watch where you're coming from, not where you're going...
c     recom rates belong to the n+1 ionisation stage
c
                  adi(ion,atom) = ad
                  bdi(ion,atom) = bd
                  t0(ion,atom) = to
                  t1(ion,atom) = ti
c
c     yet to implement aau rates...
c
               endif
            endif
 153     continue
 152  continue
c
c     read new charge exchange reactions
c
      nchxr = 0
      nchxi = 0
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, * ) nsequ
      do i = 1,nsequ
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
         read(luin, * ) nentries
         do j = 1, nentries
            atom = 0
            ion = 0
c
c     just lazy, reusing some spare locals...
c
            read (luin,*) ni,at,io,to,ti,ac,tc,ar,xr
            if (zmap(at).ne.0) then
c     
c     Check to see if the atom type is one chosen in ATDAT
c
               atom = zmap(at)
               if (io.le.maxion(atom)) then
c
c     Check to see if the ion is within the limits set in ATDAT
c
                  if (i.lt.3) then
c
c     recombination rates
c
                     nchxr = nchxr +1
                     chxrx(nchxr) = ni
                     chxrat(nchxr) = atom
                     chxrio(nchxr) = io
                     chxrtemp(1,nchxr) = to
                     chxrtemp(2,nchxr) = ti
                     chxrcos(1,nchxr) = ac*1.0d-9
                     chxrcos(2,nchxr) = tc
                     chxrcos(3,nchxr) = ar
                     chxrcos(4,nchxr) = xr
                     chxr(nchxr) = 0.d0
                  else
c     
c     ionisation rates
c
                     nchxi = nchxi +1
                     chxix(nchxi) = ni
                     chxiat(nchxi) = atom
                     chxiio(nchxi) = io
                     chxitemp(1,nchxi) = to
                     chxitemp(2,nchxi) = ti
                     chxicos(1,nchxi) = ac*1.0d-9
                     chxicos(2,nchxi) = tc
                     chxicos(3,nchxi) = ar
                     chxicos(4,nchxi) = xr
                     chxi(nchxi) = 0.d0
                  endif
c
c     end inclusion ifs
c     
               endif
            endif
c
c     end nsequ,nentries
c
         enddo
      enddo
c
      close(luin)
c
c
      return
      end
c
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
      subroutine readoprec(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c

      integer*4 at,atom,i,io,j
      integer*4 k,luin,m,nentries,ns
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
c
c
      open(luin, file='dats/OPRECDAT', status='OLD')
c
 510  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 510
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
c
      nrcspl = 0
c
      read(luin, *) nentries
      do j = 1, nentries
         read(luin,*) at,io
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               nrcspl = nrcspl+1
               recsplat(nrcspl) = atom
               recsplio(nrcspl) = io
               read(luin, *) ns
               recspln(nrcspl) = ns
               do i = 1,ns
               read(luin,*) m,(recspl(nrcspl,i,k),k= 1,5)
               enddo
            endif
         endif
      enddo
c
      close(luin)
c
      return
      end
c
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
      subroutine readcont(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c

      double precision dg20,dgu0,g20,gu0,e2,ez
      double precision zo,zto
      integer*4 at,atom,i,ii,io,ion,isos,j
      integer*4 k,luin,nentries,ni,no,nsequ,ngfbt,m
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
      open(luin, file='dats/CONTDAT', status='OLD')
 155  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 155
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read (luin, fmt = *) nentries
      do i = 1,nentries
         read (luin,*) at,io,isos,no,zto,zo,e2
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               ion = io
               fn0(ion,atom) = no
               zt0(ion,atom) = zto
               zn0(ion,atom) = zo
               epot2(ion,atom) = e2
            endif
         endif
      enddo
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
      read (luin, fmt = *) nsequ
      do i = 1,nsequ
         read(luin, fmt = 1) (ibuf(j),j = 1,11)
c         write(*, fmt = 3) (ibuf(j),j = 1,11)
         read (luin, fmt = *) nentries
         do j = 1,nentries
            read (luin,*) at,io,ez
            if (zmap(at).ne.0) then
               atom = zmap(at)
               if (io.le.maxion(atom)) then
                  ion = io
                  ezz(ion,atom) = ez*ev
               endif
            endif
         enddo
      enddo
c     
c   Read 2S0 effective recombination data 
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
c
c   Case A Hydrogen 
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
      
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) th42(i),(r2s1a(i,k),k = 1,7)
      enddo
c
c   Case B Hydrogen 
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) th42(i),(r2s1b(i,k),k = 1,7)
      enddo
c
c   Case A Helium II 
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) the43(i),(r2s2a(i,k),k = 1,7)
      enddo
c
c   Case B Helium II
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) the43(i),(r2s2b(i,k),k = 1,7)
      enddo
c
c     
c   Read Hydrogenic free bound gaunt factors
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
c
c   Temperatures (T/Z^2)
c
      read(luin,*) ngfbtz
      read(luin,*) (gfbtz(k),k = 1,ngfbtz)
      do  i = 1, ngfbtz
          gfbtz(i) = dlog10(gfbtz(i))
      enddo
c
      read(luin,*) ngfbe
c
c edges (Ryd/Z^2) and lookup tables
c
      ngfbt = ngfbe/2+1
      do  i = 1, ngfbt
          ii = (ngfbt - i) + 1
          gfbe(i) = dlog10(1.d0/dble(ii*ii))
      enddo
c
c   gfbs
c
      do  i = 1, ngfbe
          read(luin,*) (gfb(k,i),k = 1,ngfbtz)
      enddo
c
c     
c   Read free free gaunt factors
c
      read(luin, 1) (ibuf(k),k = 1,11)
c      write(*, fmt = 3) (ibuf(j),j = 1,11)
c
c   gamma2 and u (kT/Ryd), (hnu/kT)
c
      read(luin,*) g20,dg20
      read(luin,*) gu0,dgu0
c
c   gfbs & compute bins
c
      read(luin,*) ngffg2, ngffu
c      write(*,'(i4,1x,i4)') ngffg2, ngffu
c
      do  i = 1, ngffg2
          gffg2(i) = g20 + dg20*dble(i-1)
      enddo
c
      do  j = 1, ngffu
	      read(luin,*) (gff(i,j),i = 1, ngffg2)
          gffu(j) = gu0 + dgu0*dble(j-1)
      enddo
c
c    Convert to logs for interpolation
c
      do  j = 1, ngffu
      do  i = 1, ngffg2
	      gff(i,j) = dlog10(gff(i,j))
      enddo
      enddo
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) m,(gffintspl(i,k),k = 1,5)
      enddo
c
      close(luin)      
c     
      return
      end
c
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
      subroutine readphiondat(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision au,be,den,ep,g
      integer*4 at,atom,i,io,ion,j,luin,nentries,ni
      double precision s,si
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
      open(luin, file='dats/PHIONDAT2', status='OLD')
 100  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 100
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      if (mapmode.eq.1) then
         read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nentries
      ionum  = 0
      den = 0.0d0
      do 110 i = 1, nentries
         atom = 0
         ion = 0
         read (luin,*) at,io,ni,ep,au,si,be,s,g
c     
c     Check to see if the atom type is one chosen in ATDAT
c
         if (zmap(at).ne.0) then
            atom = zmap(at)
c
c     Check to see if the ion is within the limits set in ATDAT
c     
            if (io.le.maxion(atom)) then
c
               ion = io
               ionum = ionum+1
c
               atpho(ionum) = atom
               ionpho(ionum) = ion
               ipotpho(ionum) = ep
               augpho(ionum) = au
               sigpho(ionum) = si*1.0d-18
               betpho(ionum) = be
               spho(ionum) = s
               stwtpho(ionum) = g
               if (ipotpho(ionum).ge.den) goto 120
               write(*,*) ' INCOMPATIBLE dats/PHIONDAT2'
               error = .true.
               stop
 120           den = ipotpho(ionum)
            endif
         endif
 110  continue
c
      close(luin)
c
      return
      end
c
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
      subroutine readxdata(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision br,eg
      double precision rntc,xev
      double precision xf,xl,xt
      integer*4 at,atom,i,io,j,luin,nct,nentries,ni
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
c
      open(luin, file='dats/XLINDAT', status='OLD')
c
 175  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 175 
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      xlines = 0
      read(luin, *) nentries
      do  j = 1, nentries
         read(luin,*) at,io,ni,xl,xev,xf,xt,eg
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               xlines = xlines+1
               xrat(xlines) = atom
               xion(xlines) = io
               xiso(xlines) = ni
               xrlam(xlines) = xl
               xejk(xlines) = xev
               xegj(xlines) = xev*eg
               xfef(xlines) = xf
               xtrans(xlines) = xt
            endif
         endif
      enddo
c
      close(luin)
c
c
c
      open(luin, file='dats/XINTERDAT', status='OLD')
c
 177  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 177 
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      xilines = 0
      read(luin, *) nentries
      do  j = 1, nentries
         read(luin,*) at,io,ni,xl,xev,xf,xt
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               xilines = xilines+1
               xiat(xilines) = atom
               xiion(xilines) = io
               xiiso(xilines) = ni
               xilam(xilines) = xl
               xiejk(xilines) = xev
               xomeg(xilines) = xf
               xitr(xilines) = xt
            endif
         endif
      enddo
c
      close(luin)
c
c
c
      open(luin, file='dats/XHEIFDAT', status='OLD')
c
 178  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 178 
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      xhelines = 0
      read(luin, *) nentries
      do  j = 1, nentries
         read(luin,*) at,io,ni,xl,xev,xf,xt
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               xhelines = xhelines+1
               xheat(xhelines) = atom
               xheion(xhelines) = io
               xhelam(xhelines) = xl
               xhejk(xhelines) = xev
               xhef(xhelines) = xf
            endif
         endif
      enddo
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      if (runmode.ne.'batchrun') then
c      write(*, 3) (ibuf(j),j = 1, 11)
c      endif
      read(luin, *) nentries
      do j = 1, nentries
         read(luin,*) at,xf
         if (zmap(at).ne.0) then
            atom = zmap(at)
            xhebra(atom) = xf
         endif
      enddo
c
      close(luin)
c
c
c
      open(luin, file='dats/LMTDAT', status='OLD')
 176  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 176
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nentries
      do  j = 1, nentries
         read(luin,*) ni,rntc,i,br
         nct = idint(rntc)
         brr(ni,nct) = 1.d0
         if (dble(nct).eq.rntc) then
            clkup(ni,nct) = i
            brr(ni,nct) = br
         endif
      enddo
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      read(luin, *) nentries
      do  j = 1, nentries
         read(luin,*) i,arg(j),brg(j),crg(j),drg(j),erg(j)
                 if (i.ne.j) then
                    if (runmode.ne.'batchrun') then
                   write(*,*) 'error reading LMTDAT'
                   write (*,*) i,j,arg(j),brg(j),crg(j),drg(j),erg(j)
                   endif
                 endif
      enddo
c
      close(luin)
c
      return
      end
c
c
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
      subroutine read2level(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision a12,a21,airref,e12,fab
      double precision om,w1,w2
      integer*4 at,atom,i,io,j,k,luin,m
      integer*4 ml,nentries,nl,ns
      character ilgg*4, ibuf(11)*4
      logical error
c
c           Functions
c
      double precision nair,lam
c
c	Refractive index of std air for Lam>2000 A
c
      nair(lam) = 1.0d0
     &          +(643.26d0+294981.d0/(146.d0-lam*lam)
     &          + 2554.d0/(41.d0-lam*lam))/1.0d7
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
c     mean air refractive index used in old mappings
c
      airref = 1.00025552
c
      open(luin, file='dats/RESONDAT2', status='OLD')
 181  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 181
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      nlines = 0
      read(luin, *) nentries
      do 180 j = 1, nentries
         read(luin,*) at,io,a12,om,fab
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               nlines = nlines+1
               nl = nlines
               ielr(nl) = atom
               ionr(nl) = io
               e12r(nl) = a12
               omres(nl) = om
               fabsr(nl) = fab
               rlam(nl) = plk*cls / e12r(nl)
c
c              approx air refraction applied
c
               if (rlam(nl).gt.2.d-5) then
	               rlam(nl) = rlam(nl)/nair(1.d0/(rlam(nl)*1.d4))
               endif
c               
            endif
         endif
 180  continue
      close(luin)
c
c
c
      open(luin, file='dats/INTERDAT2', status='OLD')
c
 191  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 191
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      mlines = 0
      read(luin, *) nentries
      do 190 j = 1, nentries
         read(luin,*) at,io,e12,w1,w2,a21,om
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               mlines = mlines+1
               ml = mlines
               ielfs(ml) = atom
               ionfs(ml) = io
               e12fs(ml) = e12
               w1fs(ml) = w1
               w2fs(ml) = w2
               a21fs(ml) = a21
               omfs(ml) = om
               fslam(ml) = plk*cls / e12fs(ml)
c
c              approx air refraction applied
c
               if (fslam(ml).gt.2.d-5) then
	               fslam(ml) = fslam(ml)/nair(1.d0/(fslam(ml)*1.d4))
               endif
            endif
         endif
 190  continue
c
      nomg = 0
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      read(luin, *) nentries
      do j = 1, nentries
         read(luin,*) at,io
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
               nomg = nomg+1
               omgspat(nomg) = atom
               omgspio(nomg) = io
               read(luin, *) ns
               omgnsp(nomg) = ns
               do i = 1,ns
               read(luin,*) m,(omgspl(nomg,i,k),k= 1,5)
               enddo
            endif
         endif
      enddo
c
      close(luin)
c
      return
      end
c
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
      subroutine read9level(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision fi
      integer*4 atom,i,ion,fa
      integer*4 j,luin,m,n,nf,nsequ
      character ilgg*4, ibuf(11)*4
      logical error
c
c           Functions
c
      double precision nair,lam
c
c	Refractive index of std air for Lam>2000 A
c
      nair(lam) = 1.0d0
     &          +(643.26d0+294981.d0/(146.d0-lam*lam)
     &          + 2554.d0/(41.d0-lam*lam))/1.0d7
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
      open(luin, file='dats/NINEDAT', status='OLD')
 601  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 601
c
      n9ions = 0
      n9trans = 36
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nsequ
      do m = 1, nsequ
         read(luin, *) ion
         
         if (ion .ne. m) then
               write(*,*) ' INCOMPATIBLE dats/NINEDAT'
               error = .true.
               stop
         endif
         
         read(luin, fmt = *) fa,fi
         if (zmap(fa).ne.0) then
            atom = zmap(fa)
            if (fi.le.maxion(atom)) then
               n9ions = n9ions+1
               nf = n9ions
               f9ion(nf) = fi
               f9atom(nf) = atom
c
               read(luin, *) (wi9(j,nf),j = 1, 9)
               do i = 1, 9
                  read(luin, *) (ei9(j,i,nf),j = 1, 9)
               enddo
c     
               n = 0
               do i = 1, 8
                  do j = i + 1, 9
                     n = n + 1
                     f9lam(n,nf) = plk*cls / ei9(j,i,nf)
c
c              approx air refraction applied
c
		     if (f9lam(n,nf).gt.2.d-5) then
		       f9lam(n,nf) = f9lam(n,nf)
     &                               /nair(1.d0/(f9lam(n,nf)*1.d4))
		     endif
                  enddo
               enddo
c     
               do i = 1, 9
                  read(luin, *) (ai9(j,i,nf),j = 1, 9)
               enddo
c     
               do i = 1, 9
                  read(luin, *) (omi9(j,i,nf),j = 1, 9)
               enddo
c     
               do i = 1, 9
                  read(luin, *) (tdep9(j,i,nf),j = 1, 9)
               enddo
c
            else
               do i = 1,37
                  read(luin, fmt = 1) (ibuf(j),j = 1,11)
               enddo
            endif
         else
            do i = 1,37
               read(luin, fmt = 1) (ibuf(j),j = 1,11)
            enddo
         endif
      enddo
c
      close(luin)
c
      return
      end
c
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
      subroutine read6level(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision fi
      integer*4 atom,i,ion,fa
      integer*4 j,luin,m,n,nf,nsequ
      character ilgg*4, ibuf(11)*4
      logical error
c
c           Functions
c
      double precision nair,lam
c
c	Refractive index of std air for Lam>2000 A
c
      nair(lam) = 1.0d0
     &          +(643.26d0+294981.d0/(146.d0-lam*lam)
     &          + 2554.d0/(41.d0-lam*lam))/1.0d7
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
c
      open(luin, file='dats/SIXDAT', status='OLD')
 602  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 602
c
      n6ions = 0
      n6trans = 15
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nsequ
      do m = 1, nsequ
         read(luin, *) ion
         if (ion .ne. m) then
               write(*,*) ' INCOMPATIBLE dats/SIXDAT'
               error = .true.
               stop
         endif
         read(luin, fmt = *) fa,fi
         if (zmap(fa).ne.0) then
            atom = zmap(fa)
            if (fi.le.maxion(atom)) then
               n6ions = n6ions+1
               nf = n6ions
               f6ion(nf) = fi
               f6atom(nf) = atom
c
               read(luin, *) (wi6(j,nf),j = 1, 6)
               do i = 1, 6
                  read(luin, *) (ei6(j,i,nf),j = 1, 6)
               enddo
c     
               n = 0
               do i = 1, 5
                  do j = i + 1, 6
                     n = n + 1
                     f6lam(n,nf) = plk*cls / ei6(j,i,nf)
c
c              approx air refraction applied
c
		     if (f6lam(n,nf).gt.2.d-5) then
		       f6lam(n,nf) = f6lam(n,nf)
     &                               /nair(1.d0/(f6lam(n,nf)*1.d4))
		     endif
                  enddo
               enddo
c     
               do i = 1, 6
                  read(luin, *) (ai6(j,i,nf),j = 1, 6)
               enddo
c     
               do i = 1, 6
                  read(luin, *) (omi6(j,i,nf),j = 1, 6)
               enddo
c     
               do i = 1, 6
                  read(luin, *) (tdep6(j,i,nf),j = 1, 6)
               enddo
c
            else
               do i = 1,25
                  read(luin, fmt = 1) (ibuf(j),j = 1,11)
               enddo
            endif
         else
            do i = 1,25
               read(luin, fmt = 1) (ibuf(j),j = 1,11)
            enddo
         endif
      enddo
c
      close(luin)
cc
      return
      end
c
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
      subroutine read5level(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision fi
      integer*4 atom,i,ion,j,luin,m,n,nf,nsequ,fa
      character ilgg*4, ibuf(11)*4
      logical error
c
c           Functions
c
      double precision nair,lam
c
c	Refractive index of std air for Lam>2000 A
c
      nair(lam) = 1.0d0
     &          +(643.26d0+294981.d0/(146.d0-lam*lam)
     &          + 2554.d0/(41.d0-lam*lam))/1.0d7
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
      open(luin, file='dats/FIVEDAT', status='OLD')
 201  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 201
c
      nfions = 0
      nftrans = 10
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nsequ
      do m = 1, nsequ
         read(luin, *) ion
         if (ion .ne. m) then
               write(*,*) ' INCOMPATIBLE dats/FIVEDAT'
               error = .true.
               stop
         endif
         read(luin, fmt = *) fa,fi
         if (zmap(fa).ne.0) then
            atom = zmap(fa)
            if (fi.le.maxion(atom)) then
               nfions = nfions+1
               nf = nfions
               fion(nf) = fi
               fatom(nf) = atom
c
               read(luin, *) (wi(j,nf),j = 1, 5)
c
               do i = 1, 5
                  read(luin, *) (ei(j,i,nf),j = 1, 5)
               enddo
c     
               n = 0
               do i = 1, 4
                  do j = i + 1, 5
                     n = n + 1
                     flam(n,nf) = plk*cls / ei(j,i,nf)
c
c              approx air refraction applied
c
		     if (flam(n,nf).gt.2.d-5) then
	               flam(n,nf) = flam(n,nf)
     &                              /nair(1.d0/(flam(n,nf)*1.d4))
		     endif
                  enddo
               enddo
c     
               do i = 1, 5
                  read(luin, *) (ai(j,i,nf),j = 1, 5)
               enddo
c     
               do i = 1, 5
                  read(luin, *) (omi(j,i,nf),j = 1, 5)
               enddo
c     
               do i = 1, 5
                  read(luin, *) (tdep(j,i,nf),j = 1, 5)
               enddo
c
            else
               do i = 1,21
                  read(luin, fmt = 1) (ibuf(j),j = 1,11)
               enddo
            endif
         else
            do i = 1,21
               read(luin, fmt = 1) (ibuf(j),j = 1,11)
            enddo
         endif
      enddo
c
      close(luin)
c
      return
      end
c
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
      subroutine read3level(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision a21,e12,g1,g2
      integer*4 at,atom,i,io,j,luin,m,nentries,ni
      double precision p1,p2,tlam,wl
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
      open(luin, file='dats/THREEDAT', status='OLD')
 202  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 202
c
      nf3ions = 0
      nf3trans = 3
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
c 
c   Type 1 fine structure transitions
c
      read(luin, 1) (ibuf(j),j = 1, 11)
      read(luin, *) nentries
      nentries = nentries / 3
      do i = 1,nentries
         do j = 1,3
         read(luin,*) ni,at,io,wl,e12,tlam,a21,g1,p1,g2,p2,m
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
                 if (j.eq.1) then 
                 nf3ions = nf3ions+1
                 f3type(nf3ions) = ni
                 f3ion(nf3ions) = io
                 f3atom(nf3ions) = atom
                 endif
                 f3lam(j,nf3ions) = tlam
                 wi3(j,nf3ions) = wl
                 ei3(j,nf3ions) = e12
                 ai3(j,nf3ions) = a21
                 gam3(1,j,nf3ions) = g1
                 gam3(2,j,nf3ions) = g2
                 pow3(1,j,nf3ions) = p1
                 pow3(2,j,nf3ions) = p2
            endif
         endif
         enddo
      enddo
c 
c   Now type 2 transitions
c
      read(luin, 1) (ibuf(j),j = 1, 11)
      read(luin, *) nentries
      nentries = nentries / 3
      do i = 1,nentries
         do j = 1,3
         read(luin,*) ni,at,io,wl,e12,tlam,a21,g1,p1,g2,p2,m
         if (zmap(at).ne.0) then
            atom = zmap(at)
            if (io.le.maxion(atom)) then
                 if (j.eq.1) then 
                 nf3ions = nf3ions+1
                 f3type(nf3ions) = ni
                 f3ion(nf3ions) = io
                 f3atom(nf3ions) = atom
                 endif
                 f3lam(j,nf3ions) = tlam
                 wi3(j,nf3ions) = wl
                 ei3(j,nf3ions) = e12
                 ai3(j,nf3ions) = a21
                 gam3(1,j,nf3ions) = g1
                 gam3(2,j,nf3ions) = g2
                 pow3(1,j,nf3ions) = p1
                 pow3(2,j,nf3ions) = p2
            endif
      endif
      enddo
      enddo
c
      close(luin)
c
      return
      end
c
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
      subroutine readfeIIlin(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision tlam
      integer*4 atom, i,j,luin,ml,nl
      character ilgg*4, ibuf(11)*4
      character tfeid*12
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
 12   format(I3,g12.5,1x,a12)
c
      error = .false.
c
c      open(luin, file='dats/FEIIDAT', status='OLD')
c 258  read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      ilgg = ibuf(1)
c      if (ilgg(1:1).eq.'%') goto 258
c     
c      nfetrans = 0
c      if (runmode.ne.'batchrun') then
c      write(*, 3) (ibuf(j),j = 1,11)
c         endif
c      if (zmap(26).ne.0) then
c         atom = zmap(26)
c         if (2.le.maxion(atom)) then
c     
c            read(luin, *) ml,nl
c            read(luin, *) (wfe(j),j = 1, nl)
cc     
c            read(luin, *) ml,nl
c            do i = 1, ml
c               read(luin, *) (omfe(j,i),j = 1, nl)
c            enddo
c     
c            read(luin, *) ml,nl
c            do i = 1, ml
c               read(luin, *) (afe(j,i),j = 1, nl)
c            enddo
c            do i = ml+1,30
c               do j = 1,30
c                  afe(j,i) = 0.d0
c               enddo
c            enddo
c     
c            read(luin, *) ml,nl
c            do i = 1, ml
c               read(luin, *) (efe(j,i),j = 1, nl)
c            enddo
c            do i = ml+1,30
c               do j = 1,30
c                  efe(j,i) = 0.d0
c               enddo
c            enddo
c     
c            n = 0
c            read(luin, *) ml,nl
c            do i = 1, ml
c               read(luin, *) (tfela(j),j = 1, nl)
c               do j = 1,30
c                  if (tfela(j).ne.0.d0) then
c                     n = n+1
c                     felam(n) = tfela(j)*1.d4
c                  endif
c               enddo  
c            enddo
c            nfetrans = n
c      if (runmode.ne.'batchrun') then
c            write(*,*) nfetrans
c      endif
c     
c
c         endif
c      endif    
c     
c      close(luin)
c
c
c
      open(luin, file='dats/FEIIDAT4', status='OLD')
 258  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 258
c     
      nfetrans = 0
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      if (zmap(26).ne.0) then
         atom = zmap(26)
         if (2.le.maxion(atom)) then
c     
            read(luin, *) ml,nl
            read(luin, *) (wfe(j),j = 1, nl)
c     
            read(luin, *) ml,nl
            do i = 1, ml
               read(luin, *) (efe(j,i),j = 1, nl)
            enddo
c
            read(luin, *) ml,nl
            do i = 1, ml
               read(luin, *) (afe(j,i),j = 1, nl)
            enddo
c
            read(luin, *) ml,nl
            do i = 1, ml
               read(luin, *) (omfe(j,i),j = 1, nl)
            enddo
c     
c     read in line ids, different form version 1 FeII
c 
            read(luin, *) nl
            do j = 1,nl
               read(luin, 12) i,tlam,tfeid
               felamap(j) = i
               felam(i) = tlam
               feid(i) = tfeid
            enddo
c
            nfetrans = nl
c     
         endif
      endif 
c     
      close(luin)
c
c      nfetrans = 0
c
      return
      end
c
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
      subroutine readHHedata(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision frac,ti
      integer*4 i,j,k,l,m,series
      integer*4 level,line,luin,nentries,ni
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
 8    format(I2,1x,f9.2,1x,4(1pg14.7,1x))
 10   format(4(1x,I2,1x),8(1pg10.3,1x))
c
      error = .false.
c
      open(luin, file='dats/HRECDAT', status='OLD')
 196  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 196
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
c
c   Read rec rate splines
c
      do  j = 1, 3
         read(luin, 1) (ibuf(k),k = 1,11)
         read(luin,*) ni
         do  i = 1, ni
            read(luin,8) m,(hspline(j,i,k),k = 1,5)
            if (m .ne. i) then
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading HRECDAT'
               endif
               error = .true.
               stop
            endif
         enddo
      enddo
c
c   Read emissivity data
c
      read(luin, 1) (ibuf(k),k = 1,11)
c
c   Densities
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      read(luin,*) (rhhe(k),k = 1,ni)
c
c   Case A Hydrogen Hbeta
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) th42(i),(h42a(i,k),k = 1,7)
      enddo
c
c   Case B Hydrogen Hbeta
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) th42(i),(h42b(i,k),k = 1,7)
      enddo
c
c   Case A Helium II 4686
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) the43(i),(he43a(i,k),k = 1,7)
      enddo
c
c   Case B Helium II
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) the43(i),(he43b(i,k),k = 1,7)
      enddo
c
c   Hydrogen wavelengths by line, series: Hbeta = 2,2
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
         read(luin,*) (hlambda(i,k),k = 1,6)
      enddo
c     
c   Helium II wavelengths by line, series: H4686 = 1,3
c
C   Series 1 is hydrogen ionizing!!
C
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) (helambda(i,k),k = 1,6)
      enddo
c
c
c   Hydrogenic gf values (same for H and He)
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) (hydrogf(i,k),k = 1,6)
      enddo
c
      close(luin)
c
c
      open(luin, file='dats/HRECRATA', status='OLD')
 197  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 197
      if (runmode.ne.'batchrun') then
c      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) nentries
c
      line = 1
      series = 1
      do  j = 1, nentries
          read(luin, 1) (ibuf(k),k = 1,11)
          read(luin,*) ni
          do  i = 1, ni
               read(luin,*) ti,(hylratsa(line,series,i,k),k = 1,7)
          enddo
          line = line+1
          if (line.eq.11) then
              series=series+1
              line=1
          endif
      enddo
cc
      close(luin)
c
      open(luin, file='dats/HRECRATB', status='OLD')
 198  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 198
      if (runmode.ne.'batchrun') then
c      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) nentries
c
      line = 1
      series = 1
      do  j = 1, nentries
          read(luin, 1) (ibuf(k),k = 1,11)
          read(luin,*) ni
          do  i = 1, ni
               read(luin,*) ti,(hylratsb(line,series,i,k),k = 1,7)
          enddo
          line = line+1
          if (line.eq.11) then
              series=series+1
              line=1
          endif
      enddo
cc
      close(luin)
c
      open(luin, file='dats/HERECRATA', status='OLD')
 192  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 192
      if (runmode.ne.'batchrun') then
c      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) nentries
c
      line = 1
      series = 1
      do  j = 1, nentries
          read(luin, 1) (ibuf(k),k = 1,11)
          read(luin,*) ni
          do  i = 1, ni
               read(luin,*) ti,(helratsa(line,series,i,k),k = 1,7)
          enddo
          line = line+1
          if (line.eq.11) then
              series=series+1
              line=1
          endif
      enddo
cc
      close(luin)
c
c
      open(luin, file='dats/HERECRATB', status='OLD')
 193  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 193
      if (runmode.ne.'batchrun') then
c      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) nentries
c
      line = 1
      series = 1
      do  j = 1, nentries
          read(luin, 1) (ibuf(k),k = 1,11)
          read(luin,*) ni
          do  i = 1, ni
               read(luin,*) ti,(helratsb(line,series,i,k),k = 1,7)
          enddo
          line = line+1
          if (line.eq.11) then
              series=series+1
              line=1
          endif
      enddo
cc
      close(luin)
c
c
      open(luin, file='dats/HCOLLDAT', status='OLD')
 194  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 194
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) nentries
      do  j = 1, nentries
         read(luin,10) (colid(j,i),i=1,4),
     &                 (ccoln(j,k),k=1,4),
     &                 (dcoln(j,l),l=1,4)
      enddo
c
c     
c   collisional excitation branching ratios
c
c
      read(luin, 1) (ibuf(k),k = 1,11)
c
c	case A 
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) line,series,level,frac
          collha(line,series,level) = frac
      enddo
c
c	case B 
c
      read(luin, 1) (ibuf(k),k = 1,11)
      read(luin,*) ni
      do  i = 1, ni
          read(luin,*) line, series, level,frac
          collhb(line,series,level) = frac
      enddo
c
      close(luin)
c
c
      open(luin, file='dats/HYDRODAT', status='OLD')
 261  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 261
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, *) (hlam(j),j = 1, 5)
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 260 j = 1, 5
      read(luin, *) m, (balcoe(j,i),i = 1, 5)
      if (m .eq. j) goto 260
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading HYDRODAT'
               endif
               error = .true.
               stop
 260  continue
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 270 j = 1, 5
      read(luin, *) m, (balcoe(j,i),i = 6, 11)
      if (m .eq. j) goto 270
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading HYDRODAT'
               endif
               error = .true.
               stop
 270  continue
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 280 j = 1, 5
      read(luin, *) m, (hbema(j,i),i = 1, 6)
      if (m .eq. j) goto 280
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading HYDRODAT'
               endif
               error = .true.
               stop
 280  continue
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 300 j = 1, 5
      read(luin, *) m, (hbemb(j,i),i = 1, 6)
      if (m .eq. j) goto 300
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading HYDRODAT'
               endif
               error = .true.
               stop
 300  continue
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, *) (heien(i),i = 1, 4)
      heilam(1) = LmeV/(heien(1)/ev)*1.0d-8
      heilam(2) = LmeV/(heien(3)/ev)*1.0d-8
      heilam(3) = LmeV/(heien(4)/ev)*1.0d-8
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, *) (qren(i),i = 1, 4)
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do j = 1, 12
          read(luin, *) (hel0(i,j),i = 1, 5)
          heilam(j+3) = hel0(1,j)
      enddo
      nheilines = 15
      close(luin)
c
      j108m = 0
      k = 6
      do 320 j = 1, 12
      if ((dabs(hel0(1,j)-10830d-8)/hel0(1,j)).lt.3.d-2) j108m = j
      if ((hel0(2,j).gt.0.d0).and.(j108m .ne. j)) then
      hlam(k) = hel0(1,j)
      k = min0(10,k + 1)
      end if
 320  continue
c
      return
      end
c
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
      subroutine readstardat(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c

      integer*4 i,j,luin,m
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
c
      open(luin, file='dats/STARDAT', status='OLD')
 351  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 351
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 360 j = 1, 8
      read(luin, *) m
      read(luin, *) (tmet(j,i),i = 1, 4)
      read(luin, *) (tmet(j,i),i = 5, 8)
      read(luin, *) (tmet(j,i),i = 9, 11)
      if (m .eq. j) goto 360
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading STARDAT'
               endif
               error = .true.
               stop
 360  continue
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 370 j = 1, 16
      read(luin, *) m
      read(luin, *) (tno4(j,i),i = 1, 4)
      read(luin, *) (tno4(j,i),i = 5, 8)
      read(luin, *) (tno4(j,i),i = 9, 11)
      if (m .eq. j) goto 370
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading STARDAT'
               endif
               error = .true.
               stop
 370  continue
c
      read(luin, fmt = 1) (ibuf(j),j = 1, 11)
      do 380 j = 1, 16
      read(luin, *) m
      read(luin, *) (tno5(j,i),i = 1, 4)
      read(luin, *) (tno5(j,i),i = 5, 8)
      read(luin, *) (tno5(j,i),i = 9, 11)
      if (m .eq. j) goto 380
               if (runmode.ne.'batchrun') then
                  write(*,*) 'Error Reading STARDAT'
               endif
               error = .true.
               stop
 380  continue
      close(luin)
c
c
      return
      end
c
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
      subroutine readdust(luin,error)
c
      include 'cblocks.inc'
c
c           Variables
c

      integer spmx
      parameter(spmx=1201)
      double precision x0(spmx),aa(spmx),bb(spmx),cc(spmx),dd(spmx)
      double precision x, ebinc,dq,dx0,dqdx, pahsigma,pahsum
      double precision yield(5,5),yesc,y1,y2,ynu,npahy
      double precision englow, engmax, beta,eC
      integer*4 i,ii,inl,j,k,luin,nentries,drad
      character ilgg*4, ibuf(11)*4
      logical error
c
c     formats for reading files
c
 1    format(11a4)
 3    format(1h ,11a4)
c
      error = .false.
C
C     PAH Data
C
c     PAH crossection  
c	assumes Coronene C(24)H(12), area=4*pi*a^2 from Draine
c	a=3.7A for Coronene from Draine
c     new uses 10A Draine & Li Grain
c
c
c     scattering set to Zero for moment
c     have dropped 4pi from sigma as it seems to fit.                       
c
      pahsigma = ((10.d0)**2)*1.d-16
c 
c   Capicitance
      eC=(eesu**2*pi)/(2*6.6136223D-08)

c
      open(luin, file='dats/DUSTDATpah', status='OLD')
c
 275  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 275 
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
c
      npahy = 0
c
c   PAH ionization potentials
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nentries
      do j = 1, nentries
         read(luin,*) pahion(j),pahIP(j)
      enddo
c
c   PAH Yield
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nentries
      npahy = nentries
      do j = 1, nentries
         read(luin,*) yield(1,j),yield(2,j)
      enddo
c
c     Negative PAHs
c
      do i=1,2
        do j=1,infph-1
         if(ephot(j).gt.pahIP(i)) then
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
c          yesc=1.d0
          pahyield(j,i)=1.d0
          do k = 1,npahy-1
            if (ebinc.lt.(pahIP(i)+yield(1,k+1))) then
              y1 = pahIP(i)+yield(1,k)
              y2 = pahIP(i)+yield(1,k+1)
              ynu = (ebinc-y1)/(y2-y1)
              y1 = yield(2,k)
              y2 = yield(2,k+1)
              pahyield(j,i) = (y1+(ynu*(y2-y1)))
              goto 200
            endif
          enddo
 200      continue
         endif
        enddo
      enddo

c
c     Neutral & positive PAHs
c
      do i=3,4
        englow=-(pahion(i)+1)*eC
        do j=1,infph-1
         if(ephot(j).gt.pahIP(i)) then
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          engmax=ebinc-pahIP(i)
          beta=1.d0/(engmax-englow)**2
          yesc=1.d0-2.d0*beta*englow**2
          pahyield(j,i)=1.d0
          do k = 1,npahy-1
            if (ebinc.lt.(pahIP(i)+yield(1,k+1))) then
              y1 = pahIP(i)+yield(1,k)
              y2 = pahIP(i)+yield(1,k+1)
              ynu = (ebinc-y1)/(y2-y1)
              y1 = yield(2,k)
              y2 = yield(2,k+1)
              pahyield(j,i) = (y1+(ynu*(y2-y1)))*yesc
              goto 210
            endif
          enddo
 210      continue
         endif
        enddo 
      enddo 
c
c  Read PAH emission data
c  Listed as lorentz profiles,with x_0,peak (aa) & sigma (bb)
c  where F_PAH=peak*sigma^2/(sigma^2+(x-x_0)^2)
c  Also uses fixed exponential cut-offs of 1.5eV and 0.01 eV 
c   to prevent overflow
c
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin, *) nentries
      do j = 1, nentries
         read(luin,*) x0(j),aa(j),bb(j)
      enddo
      do inl = 1,infph-1
        ebinc = 0.5d0*(ephot(inl)+ephot(inl+1))
        do j=1,nentries
          pahflux(inl)=pahflux(inl)+
     &          aa(j)*bb(j)*bb(j)/(bb(j)*bb(j)+(ebinc-x0(j))**2)
        enddo
        pahflux(inl)=pahflux(inl)*
     &       exp(-(ebinc/1.5)**3)*exp(-(0.01/ebinc)**3)
      enddo
c
c     Normalise energy flux to 1 erg s-1
c
      pahsum=0.d0
      do i=1,infph-1
         ebinc=(ephot(i+1)-ephot(i))*evplk
         pahsum=pahsum+pahflux(i)*ebinc
      enddo
c
c  convert to photons/s/Hz/Sr
c
      do i=1,infph-1
         ebinc=0.5d0*(ephot(i+1)+ephot(i))*eV
         pahflux(i)=pahflux(i)/(4*pi*pahsum*ebinc)
      enddo
c
c      open(1011,file='pahflux.sou',status='UNKNOWN')
c      do i=1,infph
c         write(1011,865) ephot(i),pahflux(i)
c      enddo
c 865  format(1pg14.7,' ',1pg14.7)
c      close(1011)
c
c  Read data for neutral PAH
c
 276  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 276 
      if (runmode.ne.'batchrun') then
       write(*, 3) (ibuf(j),j = 1, 11)
      endif
c
c PAHabs
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read(luin, *) nentries
      do j = 1, nentries
        read(luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c	Remap to current binning,use logarithmic extrapolation for 
c        ebin < min x
c         
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dQ=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          if (ebinc.lt.x0(1)) then
            pahnabs(j)=pahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))
     &                *dqdx + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
         	pahnabs(j) = pahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)
     &                       +x*dd(ii))))
           else
         	pahnabs(j) = 0.d0
           endif
          endif
        enddo
c
c PAHsca
      read(luin,*)
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read(luin, *) nentries
      do j = 1, nentries
        read(luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c	Remap to current binning,use logarithmic extrapolation for 
c        ebin < min x
c         
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dQ=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          if (ebinc.lt.x0(1)) then
            pahnsca(j)=0.d0
cpahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))
c     &                *dqdx + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
         	pahnsca(j) = 0.d0
cpahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)
c     &                       +x*dd(ii))))
           else
         	pahnsca(j) = 0.d0
           endif
          endif
        enddo
c
c    PAHcos
        read(luin,*)
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
c        write(*, 3) (ibuf(j),j = 1, 11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         

        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))

          if (ebinc.lt.x0(1)) then
           pahncos(j) = aa(1)
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
         	x = ebinc - x0(ii)
         	pahncos(j) = aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
           else
         	pahncos(j) = 0.d0
           endif
          endif
        enddo
c
c  Read data for ionised PAH
c
 277  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 277 
      if (runmode.ne.'batchrun') then
       write(*, 3) (ibuf(j),j = 1, 11)
      endif
c
c PAHabs
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read(luin, *) nentries
      do j = 1, nentries
        read(luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c	Remap to current binning,use logarithmic extrapolation for 
c        ebin < min x
c         
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dQ=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          if (ebinc.lt.x0(1)) then
            pahiabs(j)=pahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))
     &*dqdx
     &              + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
         	pahiabs(j) = pahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)
     &                       +x*dd(ii))))
           else
         	pahiabs(j) = 0.d0
           endif
          endif
        enddo
c
c PAHsca
      read(luin,*)
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
c      write(*, 3) (ibuf(j),j = 1, 11)
      read(luin, *) nentries
      do j = 1, nentries
        read(luin,*) x0(j),aa(j),bb(j),cc(j),dd(j)
      enddo
c
c	Remap to current binning,use logarithmic extrapolation for 
c        ebin < min x
c         
        dx0=dlog10(x0(2))-dlog10(x0(1))
        dQ=dlog10(aa(2))-dlog10(aa(1))
        dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          if (ebinc.lt.x0(1)) then
            pahisca(j)=0.d0
cpahsigma*10.d0**((dlog10(ebinc)-dlog10(x0(1)))
c     &                 *dqdx + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
         	pahisca(j) = 0.d0
cpahsigma*(aa(ii)+x*(bb(ii)+x*(cc(ii)
c     &                       +x*dd(ii))))
           else
         	pahisca(j) = 0.d0
           endif
          endif
        enddo
c
c    PAHcos
        read(luin,*)
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
c        write(*, 3) (ibuf(j),j = 1, 11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         

        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))

          if (ebinc.lt.x0(1)) then
           pahicos(j) = aa(1)
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
         	x = ebinc - x0(ii)
         	pahicos(j) = aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
           else
         	pahicos(j) = 0.d0
           endif
          endif
        enddo
c
c     get PAH extinction
c
      do i=1,infph-1
         pahnext(i)=pahnabs(i)+pahnsca(i)
         pahiext(i)=pahiabs(i)+pahisca(i)
      enddo
      close(luin)

c
c     Dust grain absorption/extinct curves
c
C
C    Splines for graphite grains type=1
C
      open(luin, file='dats/DUSTDATgra', status='OLD')
c
279   read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 279 
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) dustbinmax
C      write(*,*) dustbinmax
      do 1000 drad=1,dustbinmax
        read(luin,*)
        read(luin,*) grainrad(drad)
c        write(*,*) drad,grainrad(drad)
C
C    Spline for Qabs
C
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning,use logarithmic extrapolation for 
c        ebin < min x
c         
          dx0=dlog10(x0(2))-dlog10(x0(1))
          dQ=dlog10(aa(2))-dlog10(aa(1))
          dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          if (ebinc.lt.x0(1)) then
            absorp(j,drad,1)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))
     &                       *dqdx + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
         	absorp(j,drad,1) = aa(ii)+x*(bb(ii)+x*(cc(ii)
     &                             +x*dd(ii)))
           else
         	absorp(j,drad,1) = 0.d0
           endif
          endif
        enddo
C
C    Spline for Qsca
C
        read(luin, *)
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         
          dx0=dlog10(x0(2))-dlog10(x0(1))
          dQ=dlog10(aa(2))-dlog10(aa(1))
          dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))

          if (ebinc.lt.x0(1)) then
            scatter(j,drad,1)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))
     &                       *dqdx + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
          	scatter(j,drad,1) = aa(ii)+x*(bb(ii)+x*(cc(ii)
     &                              +x*dd(ii)))
           else
         	scatter(j,drad,1) = 0.d0
           endif
          endif
        enddo
C
C    Spline for Gcos
C
        read(luin,*)
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         

        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))

          if (ebinc.lt.x0(1)) then
           Gcos(j,drad,1) = aa(1)
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
         	x = ebinc - x0(ii)
         	Gcos(j,drad,1) = aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
           else
         	Gcos(j,drad,1) = 0.d0
           endif
          endif
        enddo
C
C     Determine extinction
C
        do j = 1,infph-1
          extinct(j,drad,1)=absorp(j,drad,1)+scatter(j,drad,1)
        enddo
 1000 continue
      close(luin)
c
C
C    Splines for silicate grains type=2
C
      open(luin, file='dats/DUSTDATsil', status='OLD')
c
 280  read(luin, fmt = 1) (ibuf(j),j = 1,11)
      ilgg = ibuf(1)
      if (ilgg(1:1).eq.'%') goto 280 
      if (runmode.ne.'batchrun') then
      write(*, 3) (ibuf(j),j = 1, 11)
      read(luin, fmt = 1) (ibuf(j),j = 1,11)
      write(*, 3) (ibuf(j),j = 1, 11)
      endif
      read(luin,*) dustbinmax
      do 2000 drad = 1,dustbinmax
        read(luin,*)
        read(luin,*) grainrad(drad)
C
C    Spline for Qabs
C
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         
          dx0=dlog10(x0(2))-dlog10(x0(1))
          dQ=dlog10(aa(2))-dlog10(aa(1))
          dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))

          if (ebinc.lt.x0(1)) then
            absorp(j,drad,2)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))*dqdx
     &              + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
         	absorp(j,drad,2) = aa(ii)+x*(bb(ii)+x*(cc(ii)
     &                             +x*dd(ii)))
           else
          	absorp(j,drad,2) = 0.d0
           endif
          endif
        enddo
C
C    Spline for Qsca
C
        read(luin,*)
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         
          dx0=dlog10(x0(2))-dlog10(x0(1))
          dQ=dlog10(aa(2))-dlog10(aa(1))
          dqdx=dQ/dx0
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))

          if (ebinc.lt.x0(1)) then
            scatter(j,drad,2)=10.d0**((dlog10(ebinc)-dlog10(x0(1)))
     &                        *dqdx + log10(aa(1)))
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
          	x = ebinc - x0(ii)
          	scatter(j,drad,2) = aa(ii)+x*(bb(ii)+x*(cc(ii)
     &                              +x*dd(ii)))
           else
         	scatter(j,drad,2) = 0.d0
           endif
          endif
        enddo
C
C    Spline for Gcos
C
        read(luin,*)
        read(luin, fmt = 1) (ibuf(j),j = 1,11)
        read(luin, *) nentries
        do i = 1, nentries
          read(luin,*) x0(i),aa(i),bb(i),cc(i),dd(i)
        enddo
c
c	Remap to current binning
c         
        do j = 1,infph-1
          ebinc = 0.5d0*(ephot(j)+ephot(j+1))
          if (ebinc.lt.x0(1)) then
           Gcos(j,drad,2) = aa(1)
          else
           ii = 0
           do i = 1, nentries-1
              if ((ebinc.ge.x0(i)).and.(ebinc.lt.x0(i+1))) ii = i
           enddo
           if (ii.gt.0) then
         	x = ebinc - x0(ii)
         	Gcos(j,drad,2) = aa(ii)+x*(bb(ii)+x*(cc(ii)+x*dd(ii)))
           else
         	Gcos(j,drad,2) = 0.d0
           endif
          endif
        enddo
C
C     Determine extinction
C
        do j = 1,infph-1
         extinct(j,drad,2)=absorp(j,drad,2)+scatter(j,drad,2)
        enddo
 2000 continue
C
C     Calculate grain radius bin edges (+-0.025 in dex)
C
      do drad=1,dustbinmax
c        write(*,*) grainrad(drad)
        grainrad(drad)=grainrad(drad)*1d-4 !convert to cm
        gradedge(drad)=grainrad(drad)*9.440608763E-01
      enddo
      gradedge(dustbinmax+1)=grainrad(dustbinmax)*1.059253725
      close(luin)
      
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
