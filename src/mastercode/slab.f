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
c     A single slab model
c     full diffuse field, includes slab depth.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine slab
c
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision cab,de,dh,dr,dv,epotmi,exl,fi
      double precision fron,prescc,rad,t,tex,xhyf
      double precision tf,trea,trec,tsl,tstep1,tstep2,tstep3,wd00
      double precision popz(mxion, mxelem),dift,tst,tstep
c
      integer*4 l,luop,np
      integer*4 i,j
c
      character imod*4, lmod*4, kmod*4, nmod*4,wmod*4,linemod*4
      character ilgg*4, jprin*4, ndee*4, ill*4
      character fn*20,pfx*8,sfx*4,caller*4,where*16
c
c           Functions
c
      double precision feldens,fpressu,frectim
c
      luop = 23
      cab = 0.0d0
      epotmi = ryde
      where = 'Single slab model'
c
      fn = ' '
      pfx = 'slab'
      sfx = 'ss1'
      call newfile(pfx,4,sfx,3,fn)
c
      trea = epsilon
c
      open(luop, file=fn, status='NEW') 
c
 93   continue
c
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODE
c
      call zer
c
      limph = atypes
      jspot = 'NO'
      jcon = 'YES'
c
c     jspec is used here only in mode F, and in evoltem to 
c     write out a full spectrum for each step of evoltem.
c     default to 'NO'
c
      jspec = 'NO'
c
   99 write(*, 97) 
   97 format(//' Choose a single slab model type :'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     A :  Fixed degree of ionisation at a given temp.'/
     &'     B :  Equilibrium ionisation at a fixed temp.'/
     &'     C :  Equilibrium ionisation and equilibrium temp.'/
     &'     D :  Time dependent ionisation at equilibrium temp.'/
     &'     E :  Time dependent ionisation at a fixed temp.'/
     &'     F :  Time dependent ionisation and temperature.'/
     &'  :: ',$ )
      read(*, 1010) ilgg
 1010 format(a)
c
      l = 1
      if ((ilgg .eq. 'CC').or.(ilgg .eq. 'BB')) l = 2
c
      ilgg = ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg = 'A'
      if (ilgg.eq.'b') ilgg = 'B'
      if (ilgg.eq.'c') ilgg = 'C'
      if (ilgg.eq.'d') ilgg = 'D'
      if (ilgg.eq.'e') ilgg = 'E'
      if (ilgg.eq.'f') ilgg = 'F'
c
      if ((ilgg.lt.'A').or.(ilgg.gt.'F')) goto 99
c
      l = 2
  967 continue
   20 do 40 i = 1, atypes
      do 30 j = 1, maxion(i)
   30 pop(j,i) = 0.0d0
      if ((ilgg .eq. 'C').or.(ilgg .eq. 'B')) pop(l,i) = 1.0d0
      if (((ilgg .eq. 'D').or.(ilgg .eq. 'E')).or.(ilgg .eq. 'F')) 
     &then
      l = 1
      if (epot(1,i).lt.epotmi) l = 2
      pop(l,i) = 1.0d0
      end if
   40 continue
c
c
      call  popcha(where)
c
c
      call photsou(where)
c
c
      call copinto(pop, pop0)
c
c
c
 150  format(//' Give initial conditions:',/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'    t (K), dh (N), dr (cm), Geo. dilution, Case A-B:'/
     &'    (t<10 taken as a log, dh<=0 taken as a log),'/
     &'    (dr<100 taken as log, Dilution factor <= 1.0)'/
     &'    (Case A-B (H) taken as 0(A) <-> 1 (B), <0 auto)',//
     &' :: ',$)
 130  write(*,150)
      read(*, *) t, dh, dr, wdil, cab
c
      if (t.le.10) t = 10**t
c
      if (((wdil.gt.1.0d0).or.
     &     (wdil.lt.0.0d0)).or.
     &     (dh.le.0.0d0)) 
     &     goto 130
c     
      if ((t.le.0.0).and.
     &     ((ilgg.lt.'C').or.(ilgg .eq. 'E'))) 
     &     goto 130

      tf = t
c
      de = dh*zion(1)
      if ((ilgg .eq. 'A').and.(ndee .eq. 'N')) de = feldens(dh,pop)
      de = feldens(dh,pop)
c
      tstep1 = 0.0d0
      tstep2 = 0.0d0
      tstep3 = 0.0d0
c
 702  write(*, 700) 
 700  format(//' Spectrum printout required? (y/n) ',$)
      read(*, 1010) jprin
c
      jprin = jprin(1:1)
      if (jprin .eq. 'y') jprin = 'Y'
      if (jprin .eq. 'n') jprin = 'N'
c
      if ((jprin .ne. 'N').and.(jprin .ne. 'Y')) goto 702
c
      if (jprin .eq. 'Y') then
c
         call zerbuf
c
         if (ilgg .eq. 'F') then 
 802     write(*, 800) 
 800     format(//' Spectrum for each iteration? (y/n) ',$)
         read(*, 1010) jspec
c
         jspec = jspec(1:1)
         if (jspec .eq. 'y') jspec = 'Y'
         if (jspec .eq. 'n') jspec = 'N'
c
         if ((jspec .ne. 'N').and.(jspec .ne. 'Y')) goto 802
c
         if (jspec .eq. 'Y') jspec = 'YES'
         if (jspec .eq. 'N') jspec = 'NO'
c
         endif
c
      endif
c
      lmod = 'ALL'
      imod = 'ALL'
c
c
c
c     get runname
c
 300  format (a64)
 310  format(//' Give a name/code for this run: ',$)
      write(*,310)
      read (*,300) runname
c
      rad = 1.d38
      if (wdil.eq.0.5) rad = 0.d0
c
      dv = 0.0d0
      fi = 1.d0
c
      if (cab.ge.0) caseab(1) = cab
      if (cab.ge.0) caseab(2) = cab
c
      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff(fi, dh)
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
 80      continue
c
         call copinto(pop, popz)
c     
c
         if (cab.ge.0) caseab(1) = cab
         if (cab.ge.0) caseab(2) = cab
c
         call localem(t,de,dh)
         call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
         call zetaeff(fi, dh)
         call equion(t, de, dh)
c
         call difpop(pop, popz, trea, atypes, dift)
c
         call copinto(pop, popz)
c
         i = i+1
c
         if ((dift.ge.0.001).and.(i.le.4)) goto 80
c
c
      endif
c
      if ((ilgg .eq. 'A').or.(ilgg .eq. 'B')) call cool(t, de, dh)
c
c
      trec = frectim(t,de,dh)
c
      if (ilgg.gt.'B') then
         if (ilgg .eq. 'C') then
            nmod = 'EQUI'
         else
            if (ilgg .ne. 'F') then
               nmod = 'TIM'
 710           write(*, 717) zetae, qhdh
 717           format(/' Zetae :',1pg10.3,'  QHDH :',1pg10.3/
     &              ' Give time step (dt<100 as log) :',$ )
               read(*, *) tstep1
c     
               if (tstep1.lt.100.d0) tstep = 10**tstep
c     
            else
c     
               nmod = 'TIM'
 407           write(*, 417) zetae, qhdh
 417           format(//' Zetae :',1pg10.3,'  QHDH :',1pg10.3/
     &              ' Give final time and source lifetime:',$ )
               read(*, *, err=407) tst, tsl
               if (tst.le.0.d0) goto 407
               fron = 0.d0
               tstep3 = dmax1(0.d0,tst-tsl)
               tstep1 = fron*dmin1(tsl,tst)
               tstep2 = dmax1(0.0d0,dmin1(tsl,tst)-tstep1)
c     
 480           write(*, 485) 
 485  format(/' Density: Isochoric or Isobaric (c/b) :',$ )
               read(*, 1010) jden
c     
               jden = jden(1:1)
               if (jden.eq.'c') jden = 'C'
               if (jden.eq.'b') jden = 'B'
c     
               if ((jden .ne. 'C').and.(jden .ne. 'B')) goto 480
c     
            end if
            if (ill .eq. 'N') then
               call copinto(pop0, pop)
            else
               call copinto(pop, pop0)
            end if
         end if
         if (ilgg.lt.'E') then
            call teequi(t, tf, de, dh, tstep1, nmod)
            trec = frectim(tf,de,dh)
         else if (ilgg .eq. 'E') then
c
            tst = 0.0d0
            call cool(t, de, dh)
            write(*,'(4(1x,1pg14.7))') tst,de,dh,tloss
c
c     output diffuse field photon source file
c
            pfx = 'sltd'
            np = 4
            caller = 'TE'
            wmod = 'REAL'
c     
            call wpsou(caller,pfx,np,wmod,
     &           t,de,dh,
     &           0.d0,0.d0,0.d0,
     &           0.d0,dr,
     &           0.d0,0.d0,
     &           1.d0,tphot)
c     

            call timion(t, de, dh, xhyf, tstep1)
            tst = tst+tstep1
            call cool(t, de, dh)
            write(*,'(4(1x,1pg14.7))') tst,de,dh,tloss
c
c     output diffuse field photon source file
c
            pfx = 'sltd'
            np = 4
            caller = 'TE'
            wmod = 'REAL'
c     
            call wpsou(caller,pfx,np,wmod,
     &           t,de,dh,
     &           0.d0,0.d0,0.d0,
     &           0.d0,dr,
     &           0.d0,0.d0,
     &           1.d0,tphot)
c     
c     
            dift = 0.d0
c     
 8000       continue
c     
            call copinto(pop, popz)
c     
c
            if (cab.ge.0) caseab(1) = cab
            if (cab.ge.0) caseab(2) = cab
c
            call localem(t,de,dh)
            call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
            call zetaeff(fi, dh)

            call timion(t, de, dh, xhyf, tstep1)
            call cool(t, de, dh)
            tst = tst+tstep1
            write(*,'(4(1x,1pg14.7))') tst,de,dh,tloss
c
c     output diffuse field photon source file
c
            pfx = 'sltd'
            np = 4
            caller = 'TE'
            wmod = 'REAL'
c     
            call wpsou(caller,pfx,np,wmod,
     &           t,de,dh,
     &           0.d0,0.d0,0.d0,
     &           0.d0,dr,
     &           0.d0,0.d0,
     &           1.d0,tphot)
c     
c     
            call difpop(pop, popz, trea, atypes, dift)
c     
            call copinto(pop, popz)
c     
            tstep1 = tstep1*1.05
c
            i = i+1
c     
            if ((dift.ge.0.001)) goto 8000
c     
            trec = frectim(t,de,dh)
         else
            wd00 = 0.0d0
            tex = 100.d0
            exl = 1.d-20
            if (tstep1.gt.0.d0) call teequi(t, tf, de, dh, tstep1, nmod)
            prescc = fpressu(1.d4,dh,pop)
            call evoltem(tf,tf,de,dh,prescc,tstep2,0.0d0,0.0d0,luop)
            trec = frectim(tf,de,dh)
c
            if (cab.ge.0) caseab(1) = cab
            if (cab.ge.0) caseab(2) = cab
c
            call localem(t,de,dh)
            call totphot(t, dh, fi, rad, dr, dv, wd00, lmod)
            call evoltem(tf, tf, de, dh, prescc, tstep3, exl, tex, luop)
         end if
         call difpop(pop, pop0, trea, atypes, dift)
         t = tf
      end if
c     
      wmod = 'SCRN'
      call wmodel(luop,t,de,dh,dr,wmod)
c
 900   format(' Plasma Slab model diffuse field included.'/
     &' produced by MAPPINGS III v',a4,'     Run:',a64/)
      write(luop,900) theVersion,runname
c
      wmod = 'FILE'
      call wmodel(luop,t,de,dh,dr,wmod)
c
      write(luop, 2040) 
 2040 format(/36h ABUNDANCES OF THE ELEMENTS RELATIVE,
     &44h TO HYDROGEN AND THEIR STATE OF IONISATION :)
c
      call wionabal(luop,pop)
c
      write(luop, 725) zetae, qhdh, tstep1, tstep2, tstep3, 
     &dift, trec
c
  725 format(/' ZETAE:',1pg10.3,4x,5hQHDH:,1pg10.3,4x,6hTSTEP:
     &,3(1pg10.3),4x,5hDIFT:,0pf7.3,4x,5hTREC:,1pg9.2/)
      if (jprin .ne. 'N') then
      write(luop, 987) caseab(1), caseab(2)
  987 format(//' CASE A,B (H, He) : ',f4.2,f4.2)
c
      if (ilgg.lt.'E') then
         call sumdata(t, de, dh, fi, dr, dr, dr, imod)
      endif
      call avrdata
c      call spectrum(luop, kmod)
      kmod    = 'REL'
      linemod = 'LAMB'
      call spec2(luop, linemod,kmod)
      end if
c
c      goto 130
c
      close(luop)
c
c     output diffuse field photon source file
c
      pfx = 'slso'
      np = 4
      caller = 'TE'
      wmod = 'REAL'
c
         call wpsou(caller,pfx,np,wmod,
     &   t,de,dh,
     &   0.d0,0.d0,0.d0,
     &   0.d0,dr,
     &   0.d0,0.d0,
     &   1.d0,tphot)
c
      pfx = 'slfn'
      np = 4
      caller = 'TE'
      wmod = 'NFNU'
c
         call wpsou(caller,pfx,np,wmod,
     &   t,de,dh,
     &   0.d0,0.d0,0.d0,
     &   0.d0,dr,
     &   0.d0,0.d0,
     &   1.d0,tphot)
c
c      fn = ' '
c      pfx = 'slab'
c      sfx = 'phot'
c      call newfile(pfx,4,sfx,4,fn)
c
c      open(luop, file=fn, status='NEW')
c
c      do i = 1,infph-1
c         b = ephot(i)
c         db = ephot(i+1)-ephot(i)
c         dbh = db*ev/plk
c        dbl = (LmeV/(b+dd))-(LmeV/b)
c         xl = LmeV/b
c      write (luop,444) b,db,xl,dbl,dbh,cnphot(i),emidif(i),tphot(i)
c444       format(8(1pg14.6))
c      enddo
c      close(luop)
c
c     write out balance file
c
         wmod = 'NONE'
         pfx = 'slbal'
         np = 5
         call  wbal(caller,pfx,np,wmod,
     &        t,de,dh,
     &        0.d0,0.d0,0.d0,
     &        dr,dr,
     &        0.d0,0.d0)
c
 500  continue
      write(*, 2005) 
c
c
 2005 format(//38h OUTPUT CREATED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
