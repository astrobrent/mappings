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
      subroutine spec2(luop,list,mode)
c     
      include 'cblocks.inc'
c     
c     Variables
c     
      double precision eloslog,fhbtlog,forblog
      double precision chklim,totallines,total2p
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
c     
 2400 format(/' Flux Scale (H Beta = 1.000) :'/
     &        ' -----------------------------'/
     &' H-beta  :',f8.4, ' (log(ergs/cm^2/s))'/)
c
 2401 format(/' Flux Scale (H Beta = 1.000) :'/
     &        ' -----------------------------'/
     &' H-beta :',f8.4,  ' (log(ergs/s))'/)
c
 2402 format(/' Total Emission Line Flux   :'/
     &        ' ----------------------------'/
     &' F/Hbeta:',f9.4/)
c
 2403 format(/' Total Two-Photon Continuum :'/
     &        ' ----------------------------'/
     &' F/Hbeta:',f9.4/)
c
 2404 format(/' Total OIII Bowen Lines :'/
     &        ' ------------------------'/
     &' F/Hbeta:',f9.4/)
c
 2425 format(/' Flux Scale (Absolute units) :'/
     &        ' -----------------------------'/
     &        ' Units (Log) :',f8.4)
c
4000  format(//' Line List >',1pg9.2,' :'/
     &         ' ---------------------------'//
     &         ' Lambda(A)     E(eV)         Flux        Species    ',
     &         'Accuracy'/ 
     &         ' ---------------------------------------------------',
     &         '--------')     
4001  format(//' Flux Sorted Line List >',1pg9.2,' :'/
     &         ' ---------------------------------------'//
     &         ' Lambda(A)     E(eV)         Flux        Species    ',
     &         'Accuracy'/ 
     &         ' ---------------------------------------------------',
     &         '--------')     
4002  format(//' Fifty Strongest Lines :'/
     &         ' -----------------------'//
     &         ' Species     Lambda(A)      E(eV)        Flux       ',
     &         'Accuracy'/ 
     &         ' ---------------------------------------------------',
     &         '--------')     
c     
c     chklim is a lower limit for individual line fluxes to print out.
c     1e-4 saves a lot of paper!
c     
c     chklim = 0.d0
c     
      chklim = 1.0d-6
      if ((mode .eq. 'ABS').or.(mode.eq.'TOTL')) then
           chklim = epsilon
      endif
c   
c     Copy list into llist (local list) so it can be modified if
c     necessary.
c
      Llist = list
      if ((Llist.ne.'TTWN').and.
     &    (Llist.ne.'TOTL').and.
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
            if (hydroflux(line,series).gt.chklim) then
               linecount = linecount +1
               linelam(linecount) = hlambda(line,series)
               linespec(linecount) = hydroflux(line,series)
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
            if (heliflux(line,series).gt.chklim) then
               linecount = linecount +1
               linelam(linecount) = helambda(line,series)
               linespec(linecount) = heliflux(line,series)
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
               if (xhydroflux(line,series,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount) = xhlambda(line,series,i)
                 linespec(linecount) = xhydroflux(line,series,i)
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
            if (fluxf6(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount) = (f6lam(j,i)*1.d8)
                 linespec(linecount) = fluxf6(j,i)
                 lineid(linecount,1) = f6atom(i)
                 lineid(linecount,2) = f6ion(i)
                 lineacc(linecount)  = 2
            endif
         enddo
      enddo
c     
      do i = 1, nfions
         do  j = 1, nftrans
            if (fluxf(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (flam(j,i)*1.d8)
                 linespec(linecount) = fluxf(j,i)
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
            if (fluxf3(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (f3lam(j,i)*1.d4)
                 linespec(linecount) = fluxf3(j,i)
                 lineid(linecount,1) = f3atom(i)
                 lineid(linecount,2) = f3ion(i)
                 lineacc(linecount)  = 3
            endif
         enddo
      enddo
c     
      do i = 1,n9ions
         do j = 1, n9trans
            if (fluxf9(j,i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (f9lam(j,i)*1.d8)
                 linespec(linecount) = fluxf9(j,i)
                 lineid(linecount,1) = f9atom(i)
                 lineid(linecount,2) = f9ion(i)
                 lineacc(linecount)  = 2
            endif
         enddo
      enddo
c     
      do i = 1, nfetrans
         if (fluxfe(felamap(i)).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (felam(felamap(i))*1.d4)
                 linespec(linecount) = fluxfe(felamap(i))
                 lineid(linecount,1) = zmap(26)
                 lineid(linecount,2) = 2
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      do i = 1, mlines
         if (fluxi(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (fslam(i)*1.d8)
                 linespec(linecount) = fluxi(i)
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
         if (fluxxi(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = xilam(i)
                 linespec(linecount) = fluxxi(i)
                 lineid(linecount,1) = xiat(i)
                 lineid(linecount,2) = xiion(i)
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      do i = 1, nlines
         if (fluxr(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (rlam(i)*1.d8)
                 linespec(linecount) = fluxr(i)
                 lineid(linecount,1) = ielr(i)
                 lineid(linecount,2) = ionr(i)
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      do i = 1, xlines
         if (fluxx(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = xrlam(i)
                 linespec(linecount) = fluxx(i)
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
         if (fheif(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = xhelam(i)
                 linespec(linecount) = fheif(i)
                 lineid(linecount,1) = xheat(i)
                 lineid(linecount,2) = xheion(i)
                 lineacc(linecount)  = 3
         endif
      enddo
c
c Original He I lines
c
      do i = 1, nheilines
         if (fluxhei(i).gt.chklim) then
                 linecount = linecount +1
                 linelam(linecount)  = (heilam(i)*1.d8)
                 linespec(linecount) = fluxhei(i)
                 lineid(linecount,1) = zmap(2)
                 lineid(linecount,2) = 1
                 lineacc(linecount)  = 3
         endif
      enddo
c     
      if ((Llist.eq.'TOTL')) then
c
      if (mode .eq. 'REL') then
c
          fhbtlog = dlog10((fhbeta*(4.d0*pi))+epsilon)+vunilog
          forblog = dlog10(tforbi+epsilon)+vunilog
          eloslog = dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
              write(luop, 2401) fhbtlog
          else
              write(luop, 2400) fhbtlog
          endif
c
      endif
c
      if (mode .eq. 'ABS') then
c
      fhbtlog = dlog10(fhbeta+epsilon)+vunilog
      write(luop, 2425) fhbtlog
c
      end if
c     
      call heapindexsort(linecount,linelam,lineidx)
c
      totallines = 0.d0
      do i = 1, linecount
        totallines = totallines+linespec(lineidx(i))
      enddo
      write(luop, 2402) totallines
c
      total2p = hei2qa
      do i = 1, atypes
        total2p = total2p+h2q(i)
      enddo
      write(luop, 2403) total2p
      write(luop, 2404) heiioiiibfsum
c
      write(luop, 4000) chklim
c     
      endif
c
      if ((Llist.eq.'LAMB').or.(Llist.eq.'ALL')) then
c
      if (mode .eq. 'REL') then
c
          fhbtlog = dlog10((fhbeta*(4.d0*pi))+epsilon)+vunilog
          forblog = dlog10(tforbi+epsilon)+vunilog
          eloslog = dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
              write(luop, 2401) fhbtlog
          else
              write(luop, 2400) fhbtlog
          endif
c
      endif
c
      if (mode .eq. 'ABS') then
c
      fhbtlog = dlog10(fhbeta+epsilon)+vunilog
      write(luop, 2425) fhbtlog
c
      end if
c
      call heapindexsort(linecount,linelam,lineidx)
c
      totallines = 0.d0
      do i = 1, linecount
        totallines = totallines+linespec(lineidx(i))
      enddo
      write(luop, 2402) totallines
c
      total2p = hei2qa
      do i = 1, atypes
        total2p = total2p+h2q(i)
      enddo
      write(luop, 2403) total2p
      write(luop, 2404) heiioiiibfsum
c
      write(luop, 4000) chklim
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
          fhbtlog = dlog10((fhbeta*(4.d0*pi))+epsilon)+vunilog
          forblog = dlog10(tforbi+epsilon)+vunilog
          eloslog = dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
              write(luop, 2401) fhbtlog
          else
              write(luop, 2400) fhbtlog
          endif
c
      endif
c
      if (mode .eq. 'ABS') then
c
          fhbtlog = dlog10(fhbeta+epsilon)+vunilog
          write(luop, 2425) fhbtlog
c
      end if
c
      call heapindexsort(linecount,linespec,lineidx)
c
      totallines = 0.d0
      do i = 1, linecount
        totallines = totallines+linespec(lineidx(i))
      enddo
      write(luop, 2402) totallines
c
      total2p = hei2qa
      do i = 1, atypes
        total2p = total2p+h2q(i)
      enddo
      write(luop, 2403) total2p
      write(luop, 2404) heiioiiibfsum
c
      write(luop, 4001) chklim
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
      if (mode .eq. 'REL') then
c
          fhbtlog = dlog10((fhbeta*(4.d0*pi))+epsilon)+vunilog
          forblog = dlog10(tforbi+epsilon)+vunilog
          eloslog = dlog10(tlosac+epsilon)+vunilog
c
          if (jgeo.eq.'S') then
              write(luop, 2401) fhbtlog
          else
              write(luop, 2400) fhbtlog
          endif
c
      endif
c
      if (mode .eq. 'ABS') then
c
          fhbtlog = dlog10(fhbeta+epsilon)+vunilog
          write(luop, 2425) fhbtlog
c
      end if
c
      call heapindexsort(linecount,linespec,lineidx)
c
      totallines = 0.d0
      do i = 1, linecount
        totallines = totallines+linespec(lineidx(i))
      enddo
      write(luop, 2402) totallines
c
      total2p = hei2qa
      do i = 1, atypes
        total2p = total2p+h2q(i)
      enddo
      write(luop, 2403) total2p
      write(luop, 2404) heiioiiibfsum
c
      write(luop, 4002)
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
c
      subroutine heapindexsort(n,ra,idx)
c
c     heapsort (because of partial ordering already present)
c     return order in arrau idx, and ra is untouched.
c
      implicit none
c
      integer*4 n,idx(3072)
      double precision ra(3072)
c
      integer*4 i,ir,j,l,itmp
      double precision  rra
c
      do i = 1, n
      idx(i) = i
      enddo
c
      if (n.lt.2) return
c
      l  = n/2+1
      ir = n
10    continue
c
        if(l.gt.1)then
          l      = l-1
          itmp   = idx(l)
          rra    = ra(itmp)
        else
          itmp   = idx(ir)
          rra    = ra(itmp)
          idx(ir)= idx(1)
          ir     = ir-1
          if(ir.eq.1)then
            idx(1) = itmp
            return
          endif
        endif
c   
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(idx(j)).lt.ra(idx(j+1))) j = j+1
          endif
          if(rra.lt.ra(idx(j)))then
            idx(i) = idx(j)
            i      = j
            j      = j + j
          else
            j      = ir + 1
          endif
        goto 20
        endif
c       
        idx(i) = itmp
c   
      goto 10
c
      end
     









