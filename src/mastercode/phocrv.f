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
c     subroutine to compute cooling curves for photo-ionisation
c     equilibrium.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine  phocrv
c
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision band1,band2,band3,band4,bandall,cab,cspd,dr
      double precision dv,en,epotmi,fi
      double precision press,rad,rhotot,tf,tl,tnloss,trea,trec
      double precision tstep,ue,wid,wmol
      double precision ab,t,de,dh
      double precision poplog(mxion, mxelem),rlow,rhigh,rinc
      double precision sv,qlow,qhigh,qinc
      double precision q0,q1,q2,q3
c
      integer*4 l,luop,m,np
      integer*4 luions(mxelem),i,j,hesave
c
      character imod*4, lmod*4, where*16
      character linemod*4, kmod*4, nmod*4
      character ilgg*4, jprin*4, ill*4,filn(mxelem)*20
      character fn*20,pfx*8,sfx*4,tab*4,fmod*4,wmod*4,caller*4
c
c           Functions
c
      double precision fpressu,frectim,frho
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
      fn = ' '
      pfx = 'phocv'
      sfx = 'pie'
      call newfile(pfx,5,sfx,3,fn)
c
      cab = 0.0d0
      epotmi = ryde
      where = 'Photo CV'
c
c
      trea = epsilon
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODES
c
      call zer
c
      limph = atypes
      jspot = 'NO'
      jcon = 'YES'
c
      hesave = hellcoolmode
c     thin in He
      hellcoolmode = 1 
c
 10   write(*, 20) 
 20   format(//' Plane parallel photoionisation equil. series.'/
     &'  ( Optically thin single slabs )'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
      cab = 0.0d0
      epotmi = ryd*ev
c
      open(luop, file=fn, status='NEW') 
c
   99 write(*, 97) 
   97 format(//' Choose a photo-curve model type :'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     A :  Fixed radiation, variable density.'/
     &'     B :  Fixed density, variable radiation.'/
     &'  :: ',$ )
      read(*, 1010) ilgg
      ilgg = ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg = 'A'
      if (ilgg.eq.'b') ilgg = 'B'
c
      fmod = 'dens'
c
      if (ilgg.eq.'A') then
c
c     variable density
c
      fmod = 'dens'
c
      write(*,1000) ' Initial (H) density (low, <=0 as log):'
      read(*,* ) rlow 
c
      if (rlow.le.0.d0) rlow = 10.d0**rlow
c
      write(*,1000) ' Final (H) density (high, <=0 as log):'
      read(*,*) rhigh
c
      if (rhigh.le.0.d0) rhigh = 10.d0**rhigh
c
c
      write(*,1000) ' Density step factor (<=1 as log):'
      read(*,*) rinc
c
      if (rinc.gt.1.d0) rinc = dlog10(rinc)
c
      dh = 10.d0**(dlog10(rlow)-rinc)
c
c     variable density
c
      endif
c
      if (ilgg.eq.'B') then
c
         fmod = 'rads'
c
c     variable radiation
c
      write(*,1000) ' Hydrogen density (<=0 as log):'
      read(*,* ) dh
c
      if (dh.le.0.d0) dh = 10.d0**dh
c
      write(*,1000) ' Initial ionisation parameter'
     &              ,'(Q low, <=10 as log):'
      read(*,*) qlow
c
      if (qlow.le.10.d0) qlow = 10.d0**qlow
c
      write(*,1000) ' Final Ionisation parameter'
     &             ,' (Q high, <=100 as log):'
      read(*,*) qhigh
c
      if (qhigh.le.100.d0) qhigh = 10.d0**qhigh
c
c
      write(*,1000) ' Ionisation parameter step factor (<=1 as log):'
      read(*,*) qinc
c
      if (qinc.gt.1.d0) qinc = dlog10(qinc)
c
      qlow = 10.d0**(dlog10(qlow)-qinc)
c
c     variable radiation
c
      endif
c
c
c
      call  popcha(where)
c
      call photsou(where)
      
      q0 = qht
      q1 = qhi
      q2 = qhei
      q3 = qheii
c
      if (fmod.eq.'rads') then
c
c     adjust source to give desired Q
c
         sv = qlow*dh/q0
c
         do i = 1,infph
            soupho(i) = soupho(i) * sv
         enddo
c
         q1 = q1*sv
         q2 = q2*sv
         q3 = q3*sv
         q0 = q0*sv
c
      endif
c
      call copinto(pop, pop0)
c
c
c
 150  format(//' Give initial conditions:',/
     &':::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    t (K), dr (cm)'/
     &'    (t<10 taken as a log),'/
     &'    (dr<100 taken as log)'//
     &' :: ',$)
 140  write(*,150)
      read(*, *) t, dr
c
      if (dr.le.0.d0) goto 140
c
      if (t.le.10) t = 10**t
c
c     plane ll only
c
      wdil = 0.5d0
c
 120  format(//' ******************************************'/
     &         '  Ionization parameter QHDH:',1pg14.7/
     &         ' ******************************************')
      write(*,120) ((2.0d0*q0*wdil)/dh)
c
      jprin = 'N'
      tab = char(9)
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
         sfx = 'pie'
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
      write(*,*)'  '
c
c
      do i = 1,atypes
         open(luions(i),file=filn(i),status='NEW')
      enddo
c
c     
      do i = 1,atypes
c
      write(luions(i), *) fn,runname
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
c 841  format('Temp dh ntotal de  totloss totloss/(nde) taucool')
c
      l = 1
c
      cab = 0.0
c
c     need a little kick start so there are some electrons
c     to start with....
c
      pop(2,1) = 0.5d0
      pop(1,1) = 0.5d0
c
      ill = 'P'
c
c     initial temp -1 step
c
c      open(luop, file=fn, status='NEW') 
c
 900  format(' Photoionization curve calculation (optically thin):'/
     &' produced by MAPPINGS III v',a4,'  Run:',a64/)
      write(luop,900) theVersion,runname
c
      tab = char(9)
      write(luop,*) 't',tab,'de',tab,
     &'dh',tab,'en',tab,'mu (amu)',tab,'fflos',tab,'tloss',
     &tab,'eloss',tab,'Lambda'
c
c      write(luop,*) 'm',tab,'dr',tab,'dr',tab,'t',tab,'de',tab,
c     &'dh',tab,'en',tab,'Qtot',tab,'Qhi',tab,'Qhei',tab,
c     &'Qheii',tab,'ue',tab,'cspd',tab,'wmol',tab,'Vlg',tab,'Vsg'
c
      m = 0
c
 130  if (fmod.eq.'dens') then
         dh = 10.d0**(dlog10(dh)+rinc)
         if (dh.gt.rhigh) goto 500
      endif
c
      if (fmod.eq.'rads') then
c
         qlow = 10.d0**(dlog10(qlow)+qinc)
c
c     adjust source to give desired Q
c
         sv = qlow*dh/q0
c
         do i = 1,infph
            soupho(i) = soupho(i) * sv
         enddo
c
         q1 = q1*sv
         q2 = q2*sv
         q3 = q3 * sv
         q0 = q0*sv
c
         if (qlow.gt.qhigh) goto 500
c
      endif
c
      tf = t
      t = dabs(t)
c
      tstep = 0.0d0
c
      lmod = 'ALL'
      imod = 'ALL'
c
      rad = 1.d38
      if (wdil.eq.0.5) rad = 0.d0
c
      dv = 0.0d0
      fi = 1.d0
c

      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff(fi, dh)
c
      caseab(1) = 0.0d0
      caseab(2) = 0.0d0
c
      tstep = 0.0d0
      nmod = 'EQUI'
      call teequi(t, tf, de, dh, tstep, nmod)
      trec = frectim(tf,de,dh)
      t = tf
c

      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff(fi, dh)
c
      caseab(1) = 0.0d0
      caseab(2) = 0.0d0
c
      tstep = 0.0d0
      nmod = 'EQUI'
      call teequi(t, tf, de, dh, tstep, nmod)
      trec = frectim(tf,de,dh)
      t = tf
c

      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff(fi, dh)
c
      caseab(1) = 0.0d0
      caseab(2) = 0.0d0
c
      tstep = 0.0d0
      nmod = 'EQUI'
      call teequi(t, tf, de, dh, tstep, nmod)
      trec = frectim(tf,de,dh)
      t = tf
c      
      write(*,120) ((2.0d0*q0*wdil)/dh)
      write(*,'(1x,t7,a2,t19,a2,t31,a2)') 'Te','ne','nH'
      write(*,'(1x,5(1pg12.5,x)/)') t,de,dh,fflos,eloss
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
      do i = 0,infph-1
      wid = ephot(i+1)-ephot(i)
      if (ephot(i).ge.ryd) then
      bandall = bandall + tphot(i)*wid*evplk
      endif
      if ((ephot(i).ge.ryd).and.(ephot(i).lt.2.d2)) then
         band1 = band1 + tphot(i)*wid*evplk
      endif
      if ((ephot(i).ge.2.d2).and.(ephot(i).lt.4.d3)) then
         band2 = band2 + tphot(i)*wid*evplk
      endif
      if ((ephot(i).ge.4.d3).and.(ephot(i).lt.10.d3)) then
         band3 = band3 + tphot(i)*wid*evplk
      endif
      if (ephot(i).ge.4.d3) then
         band4 = band4 + tphot(i)*wid*evplk
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
c
      press = fpressu(t,dh,pop)
c
      en = 0.d0
      do i = 1,atypes
         en = en+zion(i)
      enddo
      en = en*dh
c
      ue = 3/2*(en+de)*rkb*t
      press = (en+de) * rkb*t
      tnloss = eloss/(en*de)
c
      rhotot = frho(de,dh)
      cspd = dsqrt(5/3*press/rhotot)
      wmol = (rhotot/(en+de))/amu
c
      tab = char(9)
      write(luop,1210) t,tab,de,tab,
     &dh,tab,en,tab,wmol,tab,fflos,tab,tloss,tab,eloss,tab,tnloss
c
c
c      qt = (2.0d0*q0*wdil)/dh
c      qh1 = (2.0d0*q1*wdil)/dh
c      qhe1 = (2.0d0*q2*wdil)/dh
c      qhe2 = (2.0d0*q3*wdil)/dh
c
c      write(luop,1210) m,tab,dr,tab,dr,tab,t,tab,de,tab,
c     &dh,tab,en,tab,qt,tab,qh1,tab,qhe1,tab,qhe2,tab,ue,tab,
c     &cspd,tab,wmol,tab,grainpot,tab,grainspot,tab
c
 1210 format(14(1pg12.5,A1))
c
c
c      wmod = 'SCRN'
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
      write(luop, 987) caseab(1),caseab(2)
  987 format(//' CASE A,B (H, He):',f4.2,f4.2)
      dr = 1.0d0
      call sumdata(t, de, dh, fi, dr, dr, dr, imod)
      call avrdata
c      call spectrum(luop, kmod)
      kmod    = 'REL'
      linemod = 'LAMB'
      call spec2(luop, linemod,kmod)
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
     &   1.d0,ffph)
c
c      endif
c
c
c
c      qt = dlog10(qt)
c      qh1 = dlog10(qh1)
c      qhe1 = dlog10(qhe1)
c      qhe2 = dlog10(qhe2)
c
      tl = dlog10(T)
      ab = zion(i)
c
      do i = 1 ,atypes
      open(luions(i),file=filn(i),status='OLD',access='APPEND')
c
      do j = 1,maxion(i)
c
         poplog(j,i) = 0.d0
         if (pop(j,i).gt.0d0) then
c            poplog(j,i) = pop(j,i)
            poplog(j,i) = dlog10(pop(j,i))
         endif
c
      enddo
c
c 3000 format(I3.3,A1,5(1pg12.5,A1),29(0pg12.5,A1))
c      write(luions(i),3000) m,tab,tl,tab,qt,tab,qh1,tab,qhe1,
c     &tab,qhe2,tab,(poplog(j,i),tab,j = 1,maxion(i))
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

