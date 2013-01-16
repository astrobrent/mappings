cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ***CHANGE ELEMENTAL ABUNDANCES ?
c     FINDS ZGAS RELATIVE TO THE SUN
c     Using Anders 1989 for Solar and not including H or He for Zgas.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     
c     
      subroutine abecha()
c     
      include 'cblocks.inc'
c     
      character*4 ilgg,ibuf(11)
      character*20 blanc
      character*80 fnam
      integer*4 atom,at,luf,i,j,nentries
      integer*4 ii, iii, is, it, ibell
      double precision zi,zs(mxelem)
      double precision zmt,rm, aa, rmf
      logical iexi
c     
c     
 1    format(11a4)
 3    format(' ',11a4)
c     
c     
c     read solar abundance file
c     
144   fnam = 'abund/anders.abn'
      inquire(file=fnam, exist=iexi)
c     
      if (iexi) then
c     
c     found the file...
c     
c     if the file is not present the
c     just put zion0 in zs
c     
         luf = 99
         open(unit = luf, file=fnam , status = 'OLD')
 15      read(unit = luf, fmt = 1) (ibuf(j),j = 1,11)
         ilgg = ibuf(1)
         if (ilgg(1:1).eq.'%') goto 15
         read (luf, fmt = *) nentries
         do i = 1,nentries
            read (luf,*) at,zi
            if (zmap(at).ne.0) then
               atom = zmap(at)
               if (zi.lt.0.d0) zi = 10**zi
               zs(atom)=zi
            endif
         enddo
         close(unit = luf)
      else
      do i = 1,atypes
            zs(i) = zion0(i)
      enddo
      endif
c
c
c     
c     Default elemental abundances for Local ISM (Not Solar)
c     
c     H  1.000E+00    He 1.000E-01    C  8.130E-04    N  3.720E-05
c     O  5.010E-04    Ne 7.940E-05    Mg 3.720E-05    Si 4.270E-05
c     S  1.150E-05    Ar 2.630E-06    Ca 1.320E-06    Fe 2.630E-05
cc    
      abnfile = 'Default Abundances'
c     
      do i = 1,atypes
         zion(i) = zion0(i)
      enddo
c     
 299  continue
c
      zmt = 0.d0
c
      do i = 3,atypes
         zmt = zmt + zion(i)/zs(i)
      enddo
c
      if (atypes.gt.2) then
      zmt = zmt/(atypes-2)
      else
      zmt = zs(2)/zion(2)
      endif
c
      rm = 1.d0
c     
      zgas = zmt
c     
 332  if (runmode.ne.'batchrun') then
         write(*, 300)
 300     format(//' Present elemental abundances are:'/)
         write(*, 310) (elem(j), zion(j),j = 1, atypes)
 310     format(4(4x,a2,1pe10.3))
         write(*, 327) zgas
 327     format(' (Metallicity (Zgas) is',f8.4,' that of Sol)'/)
         
         write(*, 330)
 330     format(' Change abundances (y/n) : ',$)
      endif
      
      read(*, 340) ilgg
 340  format(a)
c     
      ilgg = ilgg(1:1)
c     
      if (ilgg .eq. 'y') ilgg = 'Y'
      if (ilgg .eq. 'n') ilgg = 'N'
c     
      if (ilgg .eq. 'N') goto 130
      if (ilgg .ne. 'Y') goto 332
c     
c     
c     Change abundances....
c     
      if (runmode.ne.'batchrun') then
      write (*,400)
 400  format(' Do you want to enter abundances, read a file or cancel.'
     &        ,' (e/f/c) : ')
      endif
      read(*, 410) ilgg
 410  format(a)
      ilgg = ilgg(1:1)
c     
      if (ilgg .eq. 'f') ilgg = 'F'
      if (ilgg .eq. 'e') ilgg = 'E'
c     
      if ((ilgg .ne. 'F').and.(ilgg.ne.'E')) goto 299
c     
      if (ilgg.eq.'E') then
c     
c     enter new values
c     
         abnfile = 'Entered manually'
c     
         if (runmode.ne.'batchrun') then
           write(*, 100)
           write(*, 153)
 100       format(/' When entering abundances, -ve values are taken'/
     &           ' as log10 values: +ve values are normal numbers.'/)
 153       format(/' Enter new values (<CR>=Conserve, Xn=Multiply by n',
     &           ' E=Exit, otherwise  enter new abundance.'/
     &           '  (Use UPPER CASE ) : ')
           write(*, 141) elem(1), zion(1)
         endif
         
 141     format(' ',a2,' :',1pe10.3)
         do j = 2, atypes
            if (runmode.ne.'batchrun') then
               
               write(*, 142) elem(j), zion(j)
 142           format(a2,' :',1pe10.3,'    : ')
               
            endif
            read(*, 134, err=146) blanc
 134        format(a)
c     
c     find tokens if any
c     
            i = index(blanc(1:18),'E')
            ii = index(blanc(1:20),' ')
            iii = index(blanc(1:20),'    ')
            is = index(blanc(1:1),'X')
c     
c     
c     
            if ((i .eq. (is + 1)).and.(iii .ne. (is + 2))) goto 146
            if (i .eq. 1) goto 299
            if ((i .eq. 0).and.(ii .ne. 1)) blanc(ii:ii) = 'E'
            it = i + index(blanc(i + 1:19),' ')
            if (iii.lt.2) it = 2
            if (is .eq. 1) blanc(1:1) = ' '
c     
c     strip out any leading tokens and read into aa
c     
            blanc = blanc(it:20) // blanc(1:it - 1)
            read(blanc, 136, err=146) aa
 136        format(e20.0)
c     
c     allow log input
c     
            if (aa.lt.0.0) aa = 10**aa
c     
c     put new abundance into zion
c     
            rmf = 1.0d0
            if (is .eq. 1) rmf = aa
            if (((i .eq. 0).and.(iii .eq. 1)).or.(is .eq. 1)) 
     &           aa = zion(j)
            if (((i .eq. 0).and.(iii .eq. 2)).and.(is .eq. 1)) 
     &           rmf = rm
            if (is .eq. 1) rm = rmf
            zion(j) = aa * rmf
c     
         enddo
      endif
c     
c     end of manual input
c     
      if(ilgg.eq.'F') then
c     
c     Read an abundance file
c     
 420     fnam = ' '
         if (runmode.ne.'batchrun') then
            write (*,430)
 430        format(/' Enter file name : ')
         endif
c     
         read(*,440) fnam
 440     format(a)
c         fnam = fnam(1:index(fnam," ")-1)
         inquire(file=fnam, exist=iexi)
c     
         if (iexi) then
c     
c     found the file...
c     
c     FORMAT REQUIRED:
c     
c     see example file allen.abn
c     
c     otherwise, simply:
c     
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max 44 chars)
c     3) # of entries in the file (integer)
c     4) two columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Abundance, -ve means log abundance
c     +ve mean number abundance
c     
c     
            abnfile = fnam
            luf = 99
            open(unit = luf, file=fnam , status = 'OLD')
 155        read(unit = luf, fmt = 1) (ibuf(j),j = 1,11)
            ilgg = ibuf(1)
            if (ilgg(1:1).eq.'%') goto 155
            write(*,*) ' Read abundances from :'
            write(*, 3) (ibuf(j),j = 1, 11)
            read (luf, fmt = *) nentries
            do i = 1,nentries
               read (luf,*) at,zi
               if (zmap(at).ne.0) then
                  atom = zmap(at)
                  if (zi.lt.0.d0) zi = 10**zi
                  zion(atom)=zi
               endif
            enddo
            close(unit = luf)
c     
         else
            if (runmode.ne.'batchrun') then
               write (*,450) fnam
 450           format(' File: ',a20,' does not exist....'/
     &              'Try again? or cancel? (a/c) : ')
            endif
            read (*,460) ilgg
 460        format(a)
            ilgg = ilgg(1:1)
c     
            if (ilgg.eq.'a') ilgg = 'A'
c     
            if (ilgg.eq.'A') goto 420
         endif
         
      endif
c     
c     loop back to abundance display
c     
      goto 299
c     
c     
c     
 146  if (runmode.ne.'batchrun') then
         write(*, 148) ibell
 148     format('ERROR *********************************',a1)
      endif
c
      goto 144
c     
c     Happy....
c     
 130  continue
c
      zen = 0.0d0
      do i = 1,atypes
         zion0(i) = zion(i)
         zen = zen + zion(i)
      enddo     
c     
      return
c     
      end
