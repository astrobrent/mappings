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
c     Subroutine to setup charge exchange reactions..
c     RSS 11/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine chacha()
c
      include 'cblocks.inc'
c
      character ibell*4
      character ilgg*4, il*4, ili*4, ilt*4, ils*4
c      
      double precision chtes
      integer*4 j,i1,i2
c
c	External Functions
c
      double precision densnum
c
c    ***CHOICE OF CHARGE EXCHANGE REACTIONS WITH HYDROGEN
c
      ibell(1:4) = char(7)
 340  format (a)
c
      if (runmode.ne.'batchrun') then
        write(*,*) ' '
        write(*,*) ' Standard charge exchange reactions are being used.'
        write(*,*) ' '
      endif
c
      if (runmode.ne.'batchrun') then
 101    format(' Disable the reactions? (y/n) : ',$)
        write(*,101)
      endif
c
      read(*,340) ilgg
      ilgg = ilgg(1:1)
      if (ilgg.eq.'y') ilgg = 'Y'
      if (ilgg.eq.'n') ilgg = 'N'
      if (ilgg.eq.'Y') chargemode = 2
c
      if (chargemode.eq.1) then
c
       if (runmode.ne.'batchrun') then
 100     format(' Do you wish to modify the reaction list? (y/n) : ',$)
         write(*,100)
       endif
c
       read(*,340) ilgg
       ilgg = ilgg(1:1)
       if (ilgg.eq.'y') ilgg = 'Y'
       if (ilgg.eq.'n') ilgg = 'N'
       if (ilgg.eq.'Y') then
         chtes = ((densnum(1.d0)-zion(1))-zion(2))/(zion(1)+zion(2))
         if (chtes.gt.0.5d0) then
 503      if (runmode.ne.'batchrun') then
c
            write(*, 507) 
 507        format(//'THE USE OF CHARGE EXCHANGE REACTIONS IS NOT ',
     &           'RECOMMENDED'/'WHEN THE ABUNDANCES OF HEAVY ELEMENTS ',
     &           'ARE TO HIGH'/'DISCARD ALL OF THEM (y/n) : ',$)
          endif
c
          read(*, 340) ilgg
          ilgg = ilgg(1:1)
          if ((ilgg(1:1) .ne. 'Y').and.(ilgg .ne. 'N')) goto 503
          if (ilgg(1:1) .eq. 'Y') then
             do j = 1, 26
               charco(1,j) = - abs(charco(1,j))
             enddo
          endif
         endif
c
 200     if (runmode.ne.'batchrun') then
            write(*, 210) 
 210        format(//' Charge exchange reactions to be included :')
            do j = 1, 26
               i1 = charco(1,j)
               if (i1 .eq. 0) goto 230
               ili = 'YES'
               if (i1.lt.0) ili = 'NO '
               i1 = iabs(i1)
               i2 = charco(2,j)
               ils = '(H )'
               if (charco(7,j) .eq. 2.0) ils = '(HE)'
               ilt = '<>'
               if (charco(3,j).le.0.0) ilt = '< '
               if (charco(4,j).le.0.0) ilt = ' >'
               write(*, 220) ils, elem(i1),rom(i2), ilt, 
     &              elem(i1),rom(i2+1), ili
 220           format(2x,a4,6x,a2,a3,1x,a2,5x,a2,a3,1x,':',2x,a3)
 230        enddo
         endif
c
 235     if (runmode.ne.'batchrun') then
            write(*, 240) 
 240        format(' Do you want to alter the list (Y/N) : ',$)
         endif
c
          read(*, 244) ilgg
 244      format(a)
          if (ilgg(1:1).eq.'y') ilgg = 'Y'
          if (ilgg(1:1).eq.'n') ilgg = 'N'
c         
          if (ilgg(1:1) .eq. 'N') goto 290
          if (ilgg(1:1) .ne. 'Y') goto 235
c         
          if (runmode.ne.'batchrun') then
            write(*, 250) 
 250        format(//' Include this reaction (N/CR/Y/E)? : ')
          endif
c
          do 270 j = 1, 26
            i1 = charco(1,j)
            ili = 'YES'
            if (i1.lt.0) ili = 'NO '
            i1 = iabs(i1)
            if (i1 .eq. 0) goto 270
            i2 = charco(2,j)
            ils = '(H )'
            if (charco(7,j) .eq. 2.0) ils = '(HE)'
            ilt = '<>'
            if (charco(3,j).le.0.0) ilt = '< '
            if (charco(4,j).le.0.0) ilt = ' >'
            charco(1,j) = i1
 252        if (runmode.ne.'batchrun') then
              write(*, 254)ils,elem(i1),rom(i2),ilt,
     &              elem(i1),rom(i2+1),ili
 254          format(2x,a4,6x,a2,a3,1x,a2,5x,a2,a3,1x,':',
     &              2x,a3,5x,': ',$)
            endif

            read(*, 340) ilgg
            if (ilgg(1:1) .eq. 'e') ilgg = 'E'
            if (ilgg(1:1) .eq. 'y') ilgg = 'Y'
            if (ilgg(1:1) .eq. 'n') ilgg = 'N'
c           
            il = ilgg
            if ((ilgg.eq.'   ').or.(ilgg(1:1).eq.'E')) ilgg = ili
            if (ilgg(1:1).eq.'Y') ilgg = 'YES'
            if (ilgg(1:1).eq.'N') ilgg = 'NO '
            if ((ilgg.ne.'NO ').and.(ilgg.ne.'YES')) goto 269
            if (ilgg.eq.'NO ') charco(1,j) = - i1
            if (il(1:1).eq.'E') goto 200
            goto 270
 269        if (runmode.ne.'batchrun') then
              write(*, 148) ibell
 148          format(a4)
            endif
            goto 252
 270      enddo
          goto 200
 290    continue
       endif
c
c     end old mappings code
c
      endif
c          
      return
      end
