cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Set Depletion Pattern
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine depcha()
c
      include 'cblocks.inc'
c
      double precision zi,aa,rm,rmf
c
      integer*4 i,ibell,ii,iii,is,it,j,nentries
      integer*4 atom,at,luf
c
      logical iexi
c
      character*4 ilgg,ibuf(11)
      character*20 blanc
      character*80 fnam
c
 1    format(11a4)
 3    format(' ',11a4)
c
      depfile = 'Default Depletions'
c
 299  if (runmode.ne.'batchrun') then
      write(*, 300)
 300  format(//' Present depletion factors are:'/)
      write(*, 310) (elem(j), dion(j),j = 1, atypes)
 310  format(4(4x,a2,1pe10.3))

 332  write(*, 330)
 330  format(/' Change depletions (y/n) : ',$)
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
      if (ilgg .ne. 'Y') goto 130
c     
c
c     Change depletions....
c
 144  continue
      if (runmode.ne.'batchrun') then
 401  write (*,400)
 400  format(' Do you want to Enter depletions, read a File or Cancel.'
     &,' (e/f/c) : ',$)
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
         depfile = 'Entered manually'
c
      if (runmode.ne.'batchrun') then
      write(*, 100)
      write(*, 153)
 100  format(/' When entering depletions, -ve values are taken'/
     &' as log10 values: +ve values are normal numbers.'/)
 153  format(/' Enter new values (<CR>=Conserve, Xn=Multiply by n',
     &' E=Exit, otherwise  enter new depletion.'/
     &'  (Use UPPER CASE ) : ')
      write(*, 141) elem(1), dion(1)
      endif

 141  format(' ',a2,' :',1pe10.3)
      do j = 2, atypes
      if (runmode.ne.'batchrun') then

         write(*, 142) elem(j), dion(j)
 142     format(a2,' :',1pe10.3,'    : ',$)

      endif
 132     read(*, 134, err=146) blanc
 134     format(a)
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
 136     format(e20.0)
c
c     allow log input
c
         if (aa.lt.0.0) aa = 10**aa
c
c     put new depletion into dion
c
         rmf = 1.0d0
         if (is .eq. 1) rmf = aa
         if (((i .eq. 0).and.(iii .eq. 1)).or.(is .eq. 1)) aa = dion(j)
         if (((i .eq. 0).and.(iii .eq. 2)).and.(is .eq. 1)) rmf = rm
         if (is .eq. 1) rm = rmf
         dion(j) = aa * rmf
c
      enddo
      endif
c
c     end of manual input
c
      if(ilgg.eq.'F') then
c
c     read abundance file
c 
 420    fnam = ' '
        if (runmode.ne.'batchrun') then
         write (*,430)
 430     format(/' Enter file name : ',$)
        endif
c
        read(*,440) fnam
c        fnam=fnam(1:index(fnam, " ")-1)
 440    format(a)
        inquire(file=fnam, exist=iexi)
c
        if (iexi) then
c     
c     found the file...
c     
c     FORMAT REQUIRED:
c     
c     see example file deplete.dep
c     
c     otherwise, simply:
c     
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max 44 chars)
c     3) # of entries in the file (integer)
c     4) two columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Depletion Factor, -ve means log depletion
c     +ve mean normal number
c     
c
            depfile = fnam
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
                  dion(atom)=zi
               endif
            enddo
            close(unit = luf)
c     
         else
            if (runmode.ne.'batchrun') then
            write (*,450) fnam
 450        format(' File: ',a20,' does not exist....'/
     &           'Try again? or cancel? (a/c) : ',$)
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
 148  format('ERROR *********************************',a1)
      endif
      goto 144
c     
c     Happy....
c     
 130  continue
c
c     Finally calculate invdion
c
      do j = 1, atypes
         invdion(j)=1.d0/dion(j)
      enddo
c     
c     
      return
c     
      end
