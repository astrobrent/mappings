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
c     subroutine to compute cooling curves for collsional ionisation
c     equilibrium.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine coolc
c
c
      include 'cblocks.inc'
c
c
c
      double precision popz(mxion, mxelem),dift,ab
      double precision poplog(mxion, mxelem),t,tfin,tinc
      double precision band1,band2,band3,band4,bandall
      double precision cab,cspd,de,dh,dr,dv,en,epotmi
      double precision fi
      double precision press, rad,rhotot,tf,tl,tnloss
      double precision trea,tstep1,tstep2,tstep3
      double precision ue,wmol,wid
c
      integer*4 l,luop,m,i,j,hesave
      integer*4 np,luions(mxelem)
c
      character imod*4, lmod*4, wmod*4,where*16
      character linemod*4, spmod*4
      character ilgg*4, jprin*4, ill*4,filn(mxelem)*20
      character fn*20,pfx*8,sfx*4,tab*4,caller*4
c
c     External Functions
c
      double precision fpressu, frho
c
c
 1000 format(//,a,$)
 1010 format(a)
c
      luop = 21
c
      do i = 1,atypes
         luions(i) = 21+i
      enddo
c
      write(*,1000) ' Initial temperature (low, <=10 as log):'
      read(*,* ) t 
c
      if (t.le.10.d0) t = 10.d0**t
c
      write(*,1000) ' Final temerature (high <=10 as log):'
      read(*,*) tfin
c
      if (tfin.le.10.d0) tfin = 10.d0**tfin
c
c
      write(*,1000) ' Temperature step factor (<=1 as log):'
      read(*,*) tinc
c
      if (tinc.gt.1.d0) tinc = dlog10(tinc)
c
      t = 10.d0**(dlog10(t)-tinc)
c
      trea = epsilon
      where = 'CIE cooling'
c
      jprin = 'N'
      tab = char(9)
c
c
c
      ieln = 4
      iel1 = zmap(1)
      iel2 = zmap(6)
      iel3 = zmap(8)
      iel4 = zmap(14)
c
      do i = 1,atypes
         j = i
         fn = ' '
         pfx = elem(j)
         sfx = 'cie'
         if ((mapz(j).eq.1).or.(mapz(j).eq.6).or.
     &        (mapz(j).eq.7).or.(mapz(j).eq.8).or.
     &        (mapz(j).eq.16)) then
            call newfile(pfx,1,sfx,3,fn)
         else
            call newfile(pfx,2,sfx,3,fn)
         endif
         filn(i) = fn
      enddo
c
c
      write(*,1000) ' Run/code name for this calculation:'
      read(*,1010) runname
c
c
      do i = 1,atypes
         open(luions(i),file=filn(i),status='NEW')
      enddo
c
c     
      do i = 1,atypes
c
      write(luions(i), *) filn(i),runname
c
 3010 format(2A4,29A6)

      write(luions(i),3010)'It#','Te ',(rom(j), j = 1,maxion(i))
c
      enddo
c
      do i = 1,atypes
         close(luions(i))
      enddo
c
c      write(luop, 841)
c 841  format('Temp dh ntotal de  fflos totloss totloss/(nde) taucool')
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODE
c
      call zer
c
      hesave = hellcoolmode
c     thin in He
      hellcoolmode = 1 
c
      limph = atypes
      jspot = 'NO'
      jcon = 'YES'
      ilgg = 'B'
c
      cab = 0.0d0
      epotmi = ryde
c
      l = 1
c
      cab = 0.0
c
   20 do 40 i = 1, atypes
      do 30 j = 1, maxion(i)
   30 pop(j,i) = 0.0d0
      if ((ilgg .eq. 'C').or.(ilgg .eq. 'B')) pop(1,i) = 1.0d0
      if (((ilgg .eq. 'D').or.(ilgg .eq. 'E')).or.(ilgg .eq. 'F')) 
     &then
      l = 1
      if (epot(1,i).lt.epotmi) l = 2
      pop(l,i) = 1.0d0
      end if
   40 continue
c
c     need a little kick start so there are some electrons
c     to start with....
c
      pop(2,1) = 0.5d0
      pop(1,1) = 0.5d0
c
c
c
      call photsou(where)
c
c
      call copinto(pop, pop0)
c
c
      ill = 'P'
c
c     initial temp -1 step
c
      fn = ' '
      pfx = 'coolcv'
      sfx = 'cie'
      call newfile(pfx,6,sfx,3,fn)
      open(luop, file=fn, status='NEW') 
c
 900  format(' Cooling curve calculation (optically thin):'/
     &' produced by MAPPINGS III v',a4,'  Run:',a64/)
      write(luop,900) theVersion,runname
c
      tab = char(9)
      write(luop,*) 't    ',tab,'de   ',
     & tab,'dh   ',tab,'en   ',tab,'pop(1,1)',
     & tab,'fflos',tab,'tloss',tab,'cspd ',
     & tab,'band1',tab,'band2',tab,'band3',tab,'band4',
     & tab,'bandall'

c      write(luop,*) 'm',tab,'dr',tab,'dr',tab,'t',tab,'de',tab,
c     &'dh',tab,'en',tab,'fflos',tab,'tloss',tab,'tnloss',tab,'fhi',tab,
c     &'press',tab,'ue',tab,'cspd',tab,'wmol'
c
      m = 0
 130  t = 10.d0**(dlog10(t)+tinc)
      if (t.gt.tfin) goto 500
c
      dh = 1.d0
      wdil = 0.d0
      dr = 1.0d0
c
      tf = t
      t = dabs(t)
      de = dh*zion(1)
      tstep1 = 0.0d0
      tstep2 = 0.0d0
      tstep3 = 0.0d0
c
      lmod = 'ALL'
      imod = 'ALL'
      rad = 1.d38
      if (wdil.eq.0.5) rad = 0.d0
      dv = 0.0d0
      fi = 1.d0
c
c
c     Update Radiation Fields
c     
c
      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff(fi, dh)
c
c
      caseab(1) = cab
c
      if (ilgg .eq. 'B') then 
c
         call equion(t, de, dh)
c
         i = 0
c
         trea = 1.d-12
         dift = 0.d0
c
c
 80      continue
c
         call copinto(pop, popz)
c     
         call localem(t,de,dh)
         call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
         call zetaeff(fi, dh)
c
         call equion(t, de, dh)
c
         write(*,'(1x,t7,a2,t19,a2,t31,a2)') 'Te','ne','fhi'
         write(*,'(1x,6(1pg12.5,1x)/)') t, de, pop(1,1), zen*dh,
     &                 fflos,tloss
c
         call difpop(pop, popz, trea, atypes, dift)
c
         call copinto(pop, popz)
c
         i = i+1
c
         if ((dift.ge.0.001).and.(i.le.4)) goto 80
c
      endif
c
      if ((ilgg .eq. 'A').or.(ilgg .eq. 'B')) call cool(t, de, dh)
c
c
c      wmod = 'SCRN'
c      call wmodel(luop,t,de,dh,dr,wmod)
c
c
c
c     total energy in Inu
c
c
      band1 = 0.d0
      band2 = 0.d0
      band3 = 0.d0
      band4 = 0.d0
      bandall = 0.d0
c
      do i = 1,infph-1
       wid = (ephot(i+1)-ephot(i))*evplk
       bandall = bandall + tphot(i)*wid
       if ((ephot(i).ge.ryd).and.(ephot(i).lt.5.d2)) then
          band1 = band1 + tphot(i)*wid
       endif
       if ((ephot(i).ge.5.d2).and.(ephot(i).lt.4.5d3)) then
          band2 = band2 + tphot(i)*wid
       endif
       if ((ephot(i).ge.4.5d3).and.(ephot(i).lt.10.d3)) then
          band3 = band3 + tphot(i)*wid
       endif
       if (ephot(i).ge.10.0d3) then
          band4 = band4 + tphot(i)*wid
       endif
      enddo
c
c     times 4 pi to cover whole sky
c
      band1 = band1*4.d0*pi
      band2 = band2*4.d0*pi
      band3 = band3*4.d0*pi
      band4 = band4*4.d0*pi
      bandall=  bandall*4.d0*pi

      press = fpressu(t,dh,pop)
c
      en = zen*dh
c
      ue = 3/2*(en+de)*rkb*t
      tnloss = tloss/(en*de)
c
      rhotot = frho(de,dh)
      cspd = dsqrt(5/3*press/rhotot)
      wmol = (rhotot/(en+de))/amu
c
      tab = char(9)
      write(luop,1210) t,de,
     & dh,en,pop(1,1),fflos,tloss,cspd,band1,
     & band2,band3,band4,bandall
c
c      write(luop,1210) m,tab,dr,tab,dr,tab,t,tab,de,tab,
c     &dh,tab,en,tab,fflos,tab,tloss,tab,tnloss,tab,pop(1,1),tab,press,
c     &tab,ue,tab,cspd,tab,wmol,tab
c
 1210 format(13(1pg12.5,2x))
c
c
c      wmod = 'FILE'
c      call wmodel(luop,t,de,dh,dr,wmod)
c
c      write(luop, 2040) 
c 2040 format(/' ABUNDANCES OF THE ELEMENTS RELATIVE',
c     &' TO HYDROGEN AND THEIR STATE OF IONISATION :')
c
c      call wionabal(luop,pop)
c
c      write(luop, 725) zetae, qhdh, tstep1, tstep2, tstep3, 
c     &dift, trec
c
c  725 format(/' ZETAE:',1pg10.3,4x,5hQHDH:,1pg10.3,4x,6hTSTEP:
c     &,3(1pg10.3),4x,5hDIFT:,0pf7.3,4x,5hTREC:,1pg9.2/)
c
      if (jprin .ne. 'N') then
      write(luop, 987) caseab(1)
  987 format(//12h CASE A,B : ,f4.2)
      dr = 1.0d0
      call sumdata(t, de, dh, fi, dr, dr, dr, imod)
      call avrdata
      spmod   = 'REL'
      linemod = 'LAMB'
      call spec2(luop, linemod,spmod)
      end if
c
c
c
c     output photon source file
c
c      if (mod(m,10).eq.0) then 
         caller= 'CC'
         pfx = 'psou'
         np = 4
         wmod = 'NORM'
c
         call wpsou(caller,pfx,np,wmod,
     &   t,de,dh,
     &   0.d0,0.d0,0.d0,
     &   0.d0,dr,
     &   0.d0,0.d0,
     &   1.d0,tphot)
c
c      endif
c
c
      do i = 1 ,atypes
      open(luions(i),file=filn(i),status='OLD',access='APPEND')
c
      do j = 1,maxion(i)
C
c     temp remove logs
c
         poplog(j,i) = 0.d0
         if (pop(j,i).gt.0d0) then
            poplog(j,i) = pop(j,i)
c            poplog(j,i) = -1.d0*(dlog10(pop(j,i)))
         endif
c
c     get cie timescale for each ion
c
c         ciepop(j,i) = -1.0d0
c         if (pop(j,i).gt.1.0d-3) then
c            ciepop(j,i) = dlog10(1.d0/((rec(j,i)+col(j,i))))
c         endif
      enddo
c no logs
      tl = T
      ab = zion(i)
c
 3000 format(I3.3,A1,1pg12.5,A1,29(0pg12.5,A1))
      write(luions(i),3000) m,tab,tl,tab,(poplog(j,i),tab, 
     &j = 1,maxion(i))
c 3000 format(I3.3,A1,1pg10.3,A1,1pg11.4,A1,1pg11.4,A1,29(0pf11.4,A1))
c      write(luions(i),3000) m,tab,tl,tab,de,tab,ab,tab,(ciepop(j,i),
c     &tab,j = 1,maxion(i))
c
      close(luions(i))
c
      enddo
c      
c
      m = m+1
      goto 130
c
c
 500  continue
c
c
      close(luop) 
      write(*, 2005) fn 
c
 2005 format(//' OUTPUT CREATED IN',a14,' &&&&&&&&&&&&&&&&&&&&&&')
c
c     restore He mode
c
      hellcoolmode = hesave
c
      return 
      end

