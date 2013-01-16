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
c     General purpose routin to set the ionisation balance in 
c     the population array pop (global)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine popcha(where)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision p,fr,epotmi
c
      integer*4 i1,i2,iflag,in,nentries
      integer*4 atom,at,io,luf,i,j
c
      character blanc*20, fn*20, where*16
      character ilgg*4,ibell*4,fnam*80
      character*4 ibuf(11)
c
      logical iexi
c
c
c    ***CHOICE OF THE ELEMENT
c
      ibell(1:4) = char(7)
c
c
c
 1    format(i2)
 3    format(a)
c
c
c     Default ionisation balance, species with ionisation
c     potetials below 1 Rydberg are singly ionised, those
c     above are neutral (1 part in 10^6 to provide minimal number of
c     electrons in worst case of pure H)
c
      epotmi = ryde
c
      do  i = 1, atypes
         do  j = 3, maxion(i)
            pop(j,i) = 0.0d0
         enddo
         pop(1,i) = 0.999999d0
         if (epot(1,i).lt.epotmi) then
            pop(1,i) = 1.d-6
         endif
         pop(2,i) = 1.0d0-pop(1,i)
      enddo
c
      ionsetup = 'Default ionisation'
c
 710  format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &         '  Setting the ionisation state for ',a16/
     &        ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//
     &' Default: Elements with ionisation potentials'/
     &' above 1 Rydberg are neutral.  All others are'/
     &' singly ionised.'//)
c
      if (runmode.ne.'batchrun') then
      write(*, 710) where
      endif
c
      call copinto(pop, pop0)
c
   30 continue
c
c     Allow an ionisation equilibrium
c     calculation to be a preionisation..
c
      if (runmode.ne.'batchrun') then
      write (*,2)
      endif
 2    format(//' Choose an ionisation balance:'/
     &' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &        '    D : Default ionisation values'/
     &        '    E : Enter ionisation values'/
     &        '    F : read a ionisation balance from a File'/
     &        '    C : Calulate a preionisation balance'/
     &        '    X : eXit with current balance'//
     &        ' :: ',$)
c
      read (*,3)ilgg
      ilgg = ilgg(1:1)
c
      if(ilgg.eq.'D') ilgg = 'd'
      if(ilgg.eq.'E') ilgg = 'e'
      if(ilgg.eq.'C') ilgg = 'c'
      if(ilgg.eq.'F') ilgg = 'f'
      if(ilgg.eq.'X') ilgg = 'x'
c
      if (ilgg.eq.'d') goto 700
c
      if (ilgg.eq.'e') goto 50
c
      if (ilgg.eq.'f') then
c
c     Read balance file
c
 1420    fnam = ' '
         if (runmode.ne.'batchrun') then
            write (*,1430)
         endif
 1430    format(/' Enter file name : ',$)
         read(*,1440,err=1420) fn
 1440    format(a)
         inquire(file=fn, exist=iexi)
c
c
         if (iexi) then
c
c     Found the file...
c
c     FORMAT REQUIRED:
c
c
c     1) header marked by lines beginning with '%' are skipped
c     2) File Title (max 44 chars)
c     3) # of entries in the file (integer)
c     4) three columns of numbers
c     a) Element Z, if Z not allowed in ATDAT then entry skipped
c     b) Ion (1 is neutral etc)
c     c) Fraction (+ve numbers are fractions, -ve are logs)
c
c
c     first zero arrays
c     
            do  i = 1, atypes
               do  j = 1, maxion(i)
                  pop(j,i) = 0.0d0
               enddo
            enddo
c
 153        format(11a4)
 154        format(' ',11a4)
c
            luf = 99
            open(unit = luf, file=fn , status = 'OLD')
 155        read(unit = luf, fmt = 153) (ibuf(j),j = 1,11)
            ilgg = ibuf(1)
            if (ilgg(1:1).eq.'%') goto 155
            if (runmode.ne.'batchrun') then
            write(*,*) ' Read ion balance from :'
            write(*, 154) (ibuf(j),j = 1, 11)
            endif
            ionsetup = fn
            read (luf, fmt = *) nentries
            do i = 1,nentries
               read (luf,*) at,io,fr
               if (zmap(at).ne.0) then
                  atom = zmap(at)
                  if (fr.lt.0.d0) fr = 10**fr
                  pop(io,atom) = fr
               endif
            enddo
            close(unit = luf)
            call copinto(pop, pop0)
c     
         else
            if (runmode.ne.'batchrun') then
               write (*,156) fnam
            endif
 156        format(' File: ',a20,' does not exist....'/
     &           'Try again? or cancel? (a/c) : ',$)
            read (*,157) ilgg
 157        format(a)
            ilgg = ilgg(1:1)
c
            if (ilgg.eq.'a') ilgg = 'A'
c     
            if (ilgg.eq.'A') goto 1420
         endif
c     
         goto 30
      endif
c     
      if (ilgg.eq.'c') then
         call sinsla(where)
         call copinto(pop, pop0)
         ionsetup = 'Calculated at start'
         goto 30
      endif
c
      if (ilgg.eq.'x') goto 1000
c
      goto 30
c
c
c     Original ionisation balance entering code
c
c
 50   if (runmode.ne.'batchrun') then
         write(*, 60) (i, elem(i),i = 1, atypes)
      endif
   60 format(6(i2,1x,a2,5x))
c
      ionsetup = 'Entered manually'
c
  100 if (runmode.ne.'batchrun') then
         write(*, 110) 
      endif
c
  110 format(//' Give element number (0=end, 99=all) : ',$)
      read(*,1) in
      if (in.le.0) goto 1000
      if (in.gt.99) goto 100
c            
      i1 = in
      i2 = in
c
      if (in.eq.99) then
	      i1 = 1
	      i2 = atypes
      endif
c      
      do 800 i = i1, i2
 
      iflag = 0
c
  160 format(//' Species versus ionisation fraction for element ',
     &' : ',a2,'  @@@@@@@')
      if (runmode.ne.'batchrun') then
         write(*, 160) elem(i)
      endif
c      
  182 do  j = 1, maxion(i)
  180 format(' ',a2,a6,1x,':',3x,1pg13.6)
      if (runmode.ne.'batchrun') then
          write(*, 180) elem(i), rom(j), pop(j,i)
      endif
      enddo
c
      if (iflag .eq. 0) goto 230
c      
 190  if (runmode.ne.'batchrun') then
         write(*, 200) 
      endif
 200  format(' ::::::::::    Alter?(y/n) : ',$)
      read(*, 210, err=190) ilgg
 210  format(a)
      ilgg = ilgg(1:1)
      if (ilgg .eq. 'n') ilgg = 'N'
      if (ilgg .eq. 'y') ilgg = 'Y'
c
      if (ilgg .eq. 'N') goto 800
      if (ilgg .ne. 'Y') goto 190
c
c
c    ***ALTERING OF INITIAL VALUES
c
 230  if (runmode.ne.'batchrun') then
         write(*, 240) 
      endif
  240 format(' Enter new values (<CR>=conserve, X=eXit)', 
     &' ////////////////')
c
      iflag = 1
  250 do 600 j = 1, maxion(i)
c
  260 format(' ',a2,a6,1x,':',3x,1pg13.6,2x,':',1x,$)
      if (runmode.ne.'batchrun') then
          write(*, 260) elem(i), rom(j), pop(j,i)
      endif
c
      read(*, 270) blanc
  270 format(a20)
c  
      i1 = index(blanc(1:18),'X')
      i2 = index(blanc(1:20),' ')
c
      if ((i2 .eq. 1).and.(i1 .eq. 0)) goto 600
      if (i1 .eq. 1) goto 350
c
      read(blanc, *, err=900) p
c
      if (p.lt.0.d0) p = 10.d0**p
      if (p.gt.1.d0) goto 900
c
  320 pop(j,i) = p
c
  600 continue
c
c    ***CHECK ON NORMALISATION
c
  350 fr = 0.d0
      do 400 j = 1, maxion(i)
  400 fr = fr+pop(j,i)
      if (dabs(fr-1.d0).gt.0.1d0) goto 420
      do 410 j = 1, maxion(i)
  410 pop(j,i) = pop(j,i)/fr
      goto 500
 420  if (runmode.ne.'batchrun') then
         write(*, 430) fr, ibell
      endif
  430 format(' Normalisation error.',
     &' ERRNO :',1pg9.2,a1)
      goto 230
 500  if (runmode.ne.'batchrun') then
         write(*, 510)  elem(i)
      endif
  510 format(/' New normalised values for element #',3h : ,a2)
c
c
      goto 182
 900  if (runmode.ne.'batchrun') then
         write(*, 910) ibell
      endif
  910 format(' ERROR *******************************',a1)
c
c
      goto 250
c
c
c
  800 continue
c
      call copinto(pop, pop0)
c
      goto 30
c
c     Default ionisation balance, species with ionisation
c     potetials below 1 Rydberg are singly ionised, those
c     above are neutral
c
 700  epotmi = ryde
c
      do  i = 1, atypes
         do  j = 3, maxion(i)
            pop(j,i) = 0.0d0
         enddo
         pop(1,i) = 0.999999d0
         if (epot(1,i).lt.epotmi) then
            pop(1,i) = 1.d-6
         endif
         pop(2,i) = 1.0d0-pop(1,i)
      enddo
c
      ionsetup = 'Default ionisation'
c
      call copinto(pop, pop0)
c
      goto 1000
c
c
c
 1000 return 
c
      end
