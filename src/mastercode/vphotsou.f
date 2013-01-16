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
c       Routine to read in and set up a variable photon
c       source.
c       Inially designed to use a set of .sou files
c
c       requires filenames with the following format:
c       
c            a*XXXX.sou  to  a*YYYY.sou
c
c      where a* a fixed pattern of one or more alphanumeric characters,
c            XXXX is a four digit number with leading zeros: ie 0017
c            YYYY is another four digit number to form an unbroken
c                 sequence in steps of 1
c
c           eg:
c
c           psou0001.sou, psou0002.sou,psou0003.sou
c
c     the program asks for the prefix and the range of numbers and then
c     generates the list of files automaticly.
c
c
c       based on photsou
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c*******TO ADD A PHOTOIONISATION SOURCE IN THE VECTOR SOUPHO(I)
c       EACH POINT CORRESPONDS TO THE INTENSITY OF RADIATION
c       AT THE MIDDLE OF THE ENERGY INTERVAL : EPHOT(I),EPHOT(I+1)
c       UNITS : ERGS.CM-2.SEC-1.STERAD-1.HZ-1    (INU)
c
c       THE NUMBER OF EMERGENT PHOTONS FOR DIFFERENT SPAN IN
c       ENERGY ARE CONTAINED IN QHI,QHEI,QHEII,QHT(=SUM OF THE 3)
c
c       CALL SUBROUTINE FPOL,INTVEC
c
c     **REFERENCES : HUMMER,D.G.,MIHALAS,D.M.(1970)MNRAS,147,P.339
c                    SHIELDS,G.A.,SEARLE,L., (1978)APJ.222,P.830
c
c
c
      subroutine vphotsou()
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision eioh
      double precision souvec(mxinfph)
      double precision vsouvec(4,0:mxinfph)
      double precision scale,tim
c
      integer*4 st,nd,plen,i
      integer*4 iinf,iinto,j,ki,kii,kiii
      integer*4 luin
c
      logical iexi
c
      character caract*80, s2*4, s1*80
      character pref*80, fnam*80, sfx*4
c
      luin = 17
      eioh = Ryd
      infph = infph
c
c    ***ZERO VECTOR SOUVEC AND RELATED QUANTITIES
c
   10 do 20 i = 1, infph-1
   20 souvec(i) = 0.0d0
      do 30 i = 1, 16
      tedge(i,3) = 0.0d0
   30 tedge(i,2) = 0.0d0
      ipho = ipho+1
      zstar = 0.0d0
      alnth = 0.0d0
      teff = 0.0d0
      cut = 0.0d0
      srcfile = 'None'
c
c
c    ***READ SOURCE FILES
c
c
      scale = 1.0000d0
c
 100  pref(1:80) = ' '
      fnam(1:80) = ' '
      sfx = '.sou'
c
      write (*,110)
 110  format(/' Enter file name prefix: ',$)
      read(*,120,err=100) pref
 120  format(a)
c
      i = 0
      plen = 0
 2    continue
      i = i+1
      if (pref(i:i).eq.' ') then
         plen = i-1
         goto 3
      endif
      if (i.eq.80) goto 3
      goto 2
 3    continue
c
      if (plen.eq.0) goto 100
c
      st = 0
      nd = 0
c
      write (*,130)
 130  format(/' Give start and end file numbers: ',$)
      read(*,*) st,nd
c
      do i = st,nd
c
 1    format(i4.4)
c
      s1 = pref(1:plen)
      fnam(1:plen) = s1
c
      s2(1:4) = ' '
      write(s2, 1) i
      fnam(plen+1:plen+4) = s2
c
      fnam(plen+5:plen+9) = sfx
c
      write(*,*) fnam
c
c
      inquire(file=fnam, exist=iexi)
c     
      if (iexi) then
c     
c     found the file...
c     
c 445     write (*,450)
c 450     format(/' Scaling Factor (Y/N): ',$)
c         read(*,440) ilgg
c         ilgg = ilgg(1:1)
c
c         if (ilgg .eq. 'y') ilgg = 'Y'
c         if (ilgg .eq. 'n') ilgg = 'N'
c         if ((ilgg.ne.'Y').and.(ilgg.ne.'N')) goto 445
c
c         if (ilgg.eq.'Y') then
c 455        write(*,460)
c 460        format(/' Enter field scaling factor (>0): ',$)
c            read(*,*) scale
c            if (scale.le.0.d0) goto 455
c         endif
c
         open(luin, file=fnam, status='OLD') 
c
         do j = 1, 3000
            read(luin, 120, err=150)  caract
            ki = index(caract,'PHOTON')
            if (ki.gt.0) then
               kii = index(caract,'SOURCE')
               kiii = index(caract,'FILE')
               if ((kiii.gt.kii).and.(kii.gt.ki)) goto 160
            end if
         enddo
c
 150     write(*, 155) 
 155     format(/' STRING : " PHOTON SOURCE FILE" NOT FOUND IN FILE')
         goto 1000
 160     continue
c     
c     read and display any header..
c     
         read(luin, 120)  caract
         if (caract(2:16).eq.' ELAPSED TIME:') then
         write(*, *) (caract(2:80))
         s2 = caract(17:80)
         read(s2,*) tim
         write(*,*) tim
         vsouvec(1,0) = tim
         endif
         if (caract(1:1).eq.'%') goto 160
c     
c     
         iinto = infph-1
         iinf = 0
         read(luin, *) iinf
         if (iinf .ne. iinto) write(*, 170) iinf, iinto
 170     format(/' Number of bins (',i4,') does not coincide',
     &        ' with the number needed :',i4)
         do j = 1, iinto
            read(luin, *) souvec(j)
            souvec(j) = souvec(j)*scale
         enddo
c
 1000    continue
c
         close(luin)
c
      else
c     
 405     write(*,410) fnam
 410     format(/,' File: ',a,' Not Found !!',// )
c
      endif
c
      enddo
c
      return
c
      end
