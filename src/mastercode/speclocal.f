cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c     
c     copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c     Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     
c     this subroutine is meant to be 
c     used to write out various spectral components in a more
c     machine useable form (for plotting packages etc)
c     In general it will make looooonnnnggg lists
c     
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine speclocal(luop
     & ,tl,el,eg,dl
     & ,t,dh,de,fh1
     & ,di,dr
     & ,list,mode)
c     
      include 'cblocks.inc'
c     
c     Variables
c     
      double precision eloslog,fhbtlog,forblog
      double precision tl,el,eg,dl
      double precision t,dh,de,fh1
      double precision di,dr
      double precision chklim
c
c     Storage for sorted output.
c
      double precision linelam(3072)
      double precision linespec(3072)
      integer*4 lineid(3072,2), lineidx(3072), lineacc(3072)
      integer*4 linecount
c
      double precision shortlam(50)
      double precision shortspec(50)
      integer*4 shortlineid(50,2) , shortidx(50), shortacc(50)
      integer*4 shortcount
c     
      integer*4 i,j,line,series
      integer*4 luop
c     
      character mode*4,list*4, Llist*4
c     
c     
c     Formats:
c     
c     
 5    format(' ',1pg12.6,1x,1pg12.6,1x,1pg12.4,1x,a3,a6,1x,I4)
 6    format(' ',a3,a6,1x,1pg12.6,1x,1pg12.6,1x,1pg12.4,1x,I4)
 7    format(/' Local Emission Line Spectrum File:'/
     &        ' ----------------------------------'/)
 8    format( ' RUN:',a64/) 
 9    format(' Local Slab Properties:'/
     &' ---------------------------------------------------'/
     &'  TOTAL LOSS  :',1pg11.3/'  EFF. LOSS   :',1pg11.3/
     &'  EFF.  GAIN  :',1pg11.3/'  FRAC. RESID.:',1pg11.3/
     &'  TEMP.       :',1pg11.4/'  H DENSITY   :',1pg11.3/
     &'  ELECTR. DENS:',1pg11.3/'  FR. NEUT. H :',1pg11.3/
     &'  DISTANCE.   :',1pg11.3/'  SLAB DEPTH  :',1pg11.3/
     &' ---------------------------------------------------')
10    format('Spectrum by MAPPINGS III v',a6/)
c
 2400 format(/' Local Flux Scale (H Beta = 1.000) :'/
     &        ' -----------------------------'/
     &' H-beta  :',f8.4,
     &' (log(ergs/cm^2/s))'/
     &' TOTLOSS :',f8.4)
c
 2401 format(/' Local Flux Scale (H Beta = 1.000) :'/
     &        ' -----------------------------'/
     &' H-beta :',f8.4,
     &' (log(ergs/s))'/
     &' TOTLOSS :',f8.4)
c
 2425 format(/' Local Flux Scale (Absolute units) :'/
     &        ' -----------------------------'/
     &        ' Units (Log) :',f8.4)
c
4000  format(//' Local Line List >',1pg9.2,' :'/
     &  ' ---------------------------'//
     &  ' Lambda(A)     E(eV)         Flux        Species    Accuracy'/
     &  ' -----------------------------------------------------------')
4001  format(//' Flux Sorted Local Line List >',1pg9.2,' :'/
     &  ' ---------------------------------------'//
     &  ' Lambda(A)     E(eV)         Flux        Species    Accuracy'/
     &  ' -----------------------------------------------------------')
4002  format(//' Fifty Strongest Local Lines :'/
     &  ' -----------------------'//
     &  ' Species     Lambda(A)      E(eV)        Flux       Accuracy'/
     &  ' -----------------------------------------------------------')

c
c	Header
c
      write(luop, 7)
      write(luop, 8) runname
      write(luop, 9) tl,el,eg,dl,t,dh,de,fh1,di,dr
      write(luop,10) theVersion
c     
c     chklim is a lower limit for individual line fluxes to print out.
c     1e-4 saves a lot of paper!
c     
c      chklim = epsilon
c     
      chklim = 1.0d-4*(hbeta+epsilon)
      if (mode .eq. 'ABS') then
           chklim = epsilon
      endif
c   
c     Copy list into llist (local list) so it can be modified if
c     necessary.
c
      Llist = list
      if ((Llist.ne.'TTWN').and.
     &    (Llist.ne.'FLUX').and.
     &    (Llist.ne.'LAMB').and.
     &    (Llist.ne.'ALL')) then
          write(*,*) ' Mode error in SPEC2 : ',list
          write(*,*) ' Using lambda sorting.'
          Llist = 'LAMB'
      endif
c     
      linecount = 0
c
      do series = 1, 6
         do line = 1,10 
            if (hydrobri(line,series).gt.chklim) then
               linecount = linecount +1
               linelam(linecount) = hlambda(line,series)
               linespec(linecount) = hydrobri(line,series)
     &                               /(hbeta+epsilon)
               lineid(linecount,1) = 1
               lineid(linecount,2) = 1
               lineacc(linecount)  = 1
               if (series.eq.1) then
                    lineacc(linecount)  = 3
               endif
            endif
         enddo
      enddo
      
      do series = 1, 6
          do line = 1,10 
            if (helibri(line,series).gt.chklim) then
               linecount = linecount +1
               linelam(linecount) = helambda(line,series)
               linespec(linecount) = helibri(line,series)
     &                               /(hbeta+epsilon)
               lineid(linecount,1) = 2
               lineid(linecount,2) = 2
               lineacc(linecount)  = 1
               if (series.eq.1) then
                    lineacc(linecount)  = 3
               endif
            endif
         enddo
      enddo
c     
      do i = 3, atypes
        do series = 1, 2
           do line = 1,10 
               if (xhydrobri(line,series,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount) = xhlambda(line,series,i)
                 linespec(linecount) = xhydrobri(line,series,i)
     &                                 /(hbeta+epsilon)
                 lineid(linecount,1) = i
                 lineid(linecount,2) = mapz(i)
                 lineacc(linecount)  = 2
               endif
            enddo
         enddo
      enddo
c     
      do i = 1,n6ions
         do j = 1, n6trans
            if (f6bri(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount) = (f6lam(j,i)*1.d8)
                 linespec(linecount) = f6bri(j,i)/(hbeta+epsilon)
                 lineid(linecount,1) = f6atom(i)
                 lineid(linecount,2) = f6ion(i)
                 lineacc(linecount)  = 2
            endif
         enddo
      enddo
c     
      do i = 1, nfions
         do  j = 1, nftrans
            if (fbri(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (flam(j,i)*1.d8)
                 linespec(linecount) = fbri(j,i)/(hbeta+epsilon)
                 lineid(linecount,1) = fatom(i)
                 lineid(linecount,2) = fion(i)
                 lineacc(linecount)  = 2
                 if (xion(i).eq.1) then
                     lineacc(linecount)  = 3
                 endif
            endif
         enddo
      enddo
c     
      do i = 1,nf3ions
         do j = 1, nf3trans
            if (f3bri(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (f3lam(j,i)*1.d4)
                 linespec(linecount) = f3bri(j,i)/(hbeta+epsilon)
                 lineid(linecount,1) = f3atom(i)
                 lineid(linecount,2) = f3ion(i)
                 lineacc(linecount)  = 3
            endif
         enddo
      enddo
c     
      do i = 1,n9ions
         do j = 1, n9trans
            if (f9bri(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (f9lam(j,i)*1.d8)
                 linespec(linecount) = f9bri(j,i)/(hbeta+epsilon)
                 lineid(linecount,1) = f9atom(i)
                 lineid(linecount,2) = f9ion(i)
                 lineacc(linecount)  = 2
            endif
         enddo
      enddo
c     
      do i = 1, nfetrans
         if (febri(felamap(i)).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (felam(felamap(i))*1.d4)
                 linespec(linecount) = febri(felamap(i))
     &                                 /(hbeta+epsilon)
                 lineid(linecount,1) = zmap(26)
                 lineid(linecount,2) = 2
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      do i = 1, mlines
         if (fsbri(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (fslam(i)*1.d8)
                 linespec(linecount) = fsbri(i)/(hbeta+epsilon)
                 lineid(linecount,1) = ielfs(i)
                 lineid(linecount,2) = ionfs(i)
                 lineacc(linecount)  = 3
c
c   Special case more accurate lines: CII, NeII, SIV
c
                 if ((xrat(i).eq.zmap(6)).and.(xion(i).eq.2)) then
                     lineacc(linecount)  = 2
                 endif
                 if ((xrat(i).eq.zmap(10)).and.(xion(i).eq.2)) then
                     lineacc(linecount)  = 2
                 endif
                 if ((xrat(i).eq.zmap(16)).and.(xion(i).eq.4)) then
                     lineacc(linecount)  = 2
                 endif
         endif
      enddo
c     
      do i = 1, xilines
         if (xibri(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = xilam(i)
                 linespec(linecount) = xibri(i)/(hbeta+epsilon)
                 lineid(linecount,1) = xiat(i)
                 lineid(linecount,2) = xiion(i)
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      do i = 1, nlines
         if (rbri(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (rlam(i)*1.d8)
                 linespec(linecount) = rbri(i)/(hbeta+epsilon)
                 lineid(linecount,1) = ielr(i)
                 lineid(linecount,2) = ionr(i)
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      do i = 1, xlines
         if (xrbri(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = xrlam(i)
                 linespec(linecount) = xrbri(i)/(hbeta+epsilon)
                 lineid(linecount,1) = xrat(i)
                 lineid(linecount,2) = xion(i)
                 lineacc(linecount)  = 3
c
c   Special case more accurate lines: CIV, NV, OVI
c
                 if ((xrat(i).eq.zmap(6)).and.(xion(i).eq.4)) then
                     lineacc(linecount)  = 2
                 endif
                 if ((xrat(i).eq.zmap(7)).and.(xion(i).eq.5)) then
                     lineacc(linecount)  = 2
                 endif
                 if ((xrat(i).eq.zmap(8)).and.(xion(i).eq.6)) then
                     lineacc(linecount)  = 2
                 endif
         endif
      enddo
c     
      do i = 1, xhelines
         if (xhebri(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = xhelam(i)
                 linespec(linecount) = xhebri(i)/(hbeta+epsilon)
                 lineid(linecount,1) = xheat(i)
                 lineid(linecount,2) = xheion(i)
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      if ((Llist.eq.'LAMB').or.(Llist.eq.'ALL')) then
c
      if (mode .eq. 'REL') then
c
          fhbtlog = dlog10((hbeta*(4.d0*pi))+epsilon)+vunilog
          forblog = dlog10(floss+epsilon)+vunilog
          eloslog = dlog10(eloss+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
              write(luop, 2401) fhbtlog, eloslog
          else
              write(luop, 2400) fhbtlog, eloslog
          endif
c
      endif
c
      if (mode .eq. 'ABS') then
c
      fhbtlog = dlog10(hbeta+epsilon)+vunilog
      write(luop, 2425) fhbtlog
c
      end if
c
      write(luop, 4000) chklim
c     
c
      call heapindexsort(linecount,linelam,lineidx)
c
      do i = 1, linecount
          write(luop,5) linelam(lineidx(i)),
     &                  LmeV/(linelam(lineidx(i))),
     &                  linespec(lineidx(i)),
     &                  elem(lineid(lineidx(i),1)),
     &                  rom(lineid(lineidx(i),2)),
     &                  lineacc(lineidx(i))
      enddo
c
      endif
c     
      if ((Llist.eq.'FLUX').or.(Llist.eq.'ALL')) then
c     
c
      if (mode .eq. 'REL') then
c
          fhbtlog = dlog10((hbeta*(4.d0*pi))+epsilon)+vunilog
          forblog = dlog10(floss+epsilon)+vunilog
          eloslog = dlog10(eloss+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
              write(luop, 2401) fhbtlog, eloslog
          else
              write(luop, 2400) fhbtlog, eloslog
          endif
c
      endif
c
      if (mode .eq. 'ABS') then
c
          fhbtlog = dlog10(hbeta+epsilon)+vunilog
          write(luop, 2425) fhbtlog
c
      end if
c
      write(luop, 4001) chklim
c
      call heapindexsort(linecount,linespec,lineidx)
c
      do i = linecount, 1, -1
          write(luop,5) linelam(lineidx(i)),
     &                  LmeV/(linelam(lineidx(i))),
     &                  linespec(lineidx(i)),
     &                  elem(lineid(lineidx(i),1)),
     &                  rom(lineid(lineidx(i),2)),
     &                  lineacc(lineidx(i))
      enddo
c
      endif
c     
c     
      if (Llist.eq.'TTWN') then
c     
      write(luop, 4002)
c
      call heapindexsort(linecount,linespec,lineidx)
c
      j = 1
      if (linecount.gt.50) j = linecount - 49
c
c     Copy top 50 lines into short array for lambda sorting
c
      shortcount = linecount - j + 1
      do i = j,linecount
           shortlineid(i-j+1,1) = lineid(lineidx(i),1)
           shortlineid(i-j+1,2) = lineid(lineidx(i),2)
           shortlam(i-j+1)      = linelam(lineidx(i))
           shortspec(i-j+1)     = linespec(lineidx(i))
           shortacc(i-j+1)      = lineacc(lineidx(i))
      enddo
c
      call heapindexsort(shortcount, shortlam, shortidx)
c
      do i = 1, shortcount
          write(luop,6) elem(shortlineid(shortidx(i),1)),
     &                  rom(shortlineid(shortidx(i),2)),
     &                  shortlam(shortidx(i)),
     &                  LmeV/(shortlam(shortidx(i))),
     &                  shortspec(shortidx(i)),
     &                  shortacc(shortidx(i))
      enddo
c
      endif
c     
      return 
      end
