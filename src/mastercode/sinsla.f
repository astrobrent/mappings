cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       single slab ionisation balance routine
c       Without the 'On-The-Spot' approximation
c       and without file output.  Used by popcha
c       for now.
c
c     rss 9.90
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sinsla(where)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision cab,de,dh,dr,dv,epotmi,exl,fi,fron
      double precision prescc,rad,t
      double precision tex,tf,trea,trec,tsl,tst,tstep1,tstep2
      double precision tstep3,wd00,xhyf
      double precision popz(mxion, mxelem),dift
c
      integer*4 i,j,l,luop
c
      character imod*4, lmod*4,  nmod*4
      character ilgg*4, ill*4,where*16
c
c           Functions
c
      double precision fpressu,frectim
c

c
      luop = 23
      cab = 0.0d0
      epotmi = ryde
c
c
      trea = epsilon
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODE
c
      call zer
c
      limph = atypes
      jspot = 'NO'
      jcon = 'YES'
c
c
c
   99 write(*, 97) 
   97 format(//' Choose a single slab model:'/
     &' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &'     B :  Equilibrium ionisation at a fixed temp.'/
     &'     C :  Equilibrium ionisation and equilibrium temp.'/
     &'     D :  Time dependent ionisation at equilibrium temp.'/
     &'     E :  Time dependent ionisation at a fixed temp.'/
     &'     F :  Time dependent ionisation and temperature.'/
     &'     X :  Exit.'/
     &'  :: ',$ )
      read(*, 1010) ilgg
c
 1010 format(a)
      ilgg = ilgg(1:1)
c
      if (ilgg.eq.'b') ilgg = 'B'
      if (ilgg.eq.'c') ilgg = 'C'
      if (ilgg.eq.'d') ilgg = 'D'
      if (ilgg.eq.'e') ilgg = 'E'
      if (ilgg.eq.'f') ilgg = 'F'
      if (ilgg.eq.'x') ilgg = 'X'
c
      l = 1
c      if ((ilgg .eq. 'CC').or.(ilgg .eq. 'BB')) l = 2
      if (ilgg .eq. 'X') goto 500
c
      if ((ilgg.lt.'B').or.(ilgg.gt.'F')) goto 99
c
   20 do 40 i = 1, atypes
      do 30 j = 1, maxion(i)
   30 pop(j,i) = 0.0d0
c
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
      call photsou(where)
c
      call copinto(pop, pop0)
c
c
      if (ilgg.gt.'B') ill = 'N'
c
c      t = 8912.5093813375d0
c
 130  write(*, 150) 
 150  format(/' Give plasma slab conditions:'/
     &' ( T (K) T<10 as log, dr in cm, Dilution <= 1.0 )'/
     &' : T, H density, slab thickness (dr), Geo. Dilution'/
     &' : ',$)
c PJM 21-JUL-08: Changed input to allow comments in mapscript file.
c      read(*, *) t, dh, dr, wdil
c      read(*, *) t
c      read(*, *) dh
c      read(*, *) dr
c      read(*, *) wdil
      read(*, *) t, dh, dr, wdil
      if (t.le.10) t = 10**t
      if (((wdil.gt.1.0d0).or.(wdil.lt.0.0d0)).or.(dh.le.0.0d0)) 
     &goto 130
      if ((t.le.0.0).and.((ilgg.lt.'C').or.(ilgg .eq. 'E'))) 
     &goto 130
      tf = t
      t = dabs(t)
      de = dh*zion(1)
c
      tstep1 = 0.0d0
      tstep2 = 0.0d0
      tstep3 = 0.0d0
c
      lmod = 'SO'
      imod = 'ALL'
      rad = 1.d38
      dr = 1.0d0
      dv = 0.0d0
      fi = 1.d0
c
c
      caseab(1) = cab
      caseab(2) = cab
c
      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
      call zetaeff(fi, dh)
c
      if (ilgg .eq. 'B') then 
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
c         write (*,*) i,dift
c
      endif
c
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
 717           format(/8h ZETAE :,1pg10.3,4x,6hQHDH :,1pg10.3/
     &              18h$GIVE TIME STEP : ,$)
               read(*, *) tstep1
c
               if (tstep1.lt.0.d0) goto 710
            else
               nmod = 'TIM'
 407           write(*, 417) zetae, qhdh
 417           format(/8h ZETAE :,1pg10.3,4x,6hQHDH :,1pg10.3/
     &              43h GIVE AGE OF NEBULA AND SOURCE LIFE-TIME : ,$)
               read(*, *) tst, tsl
               if (tst.le.0.d0) goto 407
               fron = 0.d0
               tstep3 = dmax1(0.d0,tst-tsl)
               tstep1 = fron*dmin1(tsl,tst)
               tstep2 = dmax1(0.0d0,dmin1(tsl,tst)-tstep1)
c
 480           write(*, 485) 
 485  format(/50h DENSITY BEHAVIOR , ISOCHORIC OR ISOBARIC (C/B) : ,$)
      read(*, 1010) jden
      if ((jden .ne. 'C').and.(jden .ne. 'B')) goto 480
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
         call timion(t, de, dh, xhyf, tstep1)
         call cool(t, de, dh)
         trec = frectim(t,de,dh)
      else
         wd00 = 0.0d0
         tex = 100.d0
         exl = 1.d-20
         if (tstep1.gt.0.d0) call teequi(t, tf, de, dh, tstep1, nmod)
         prescc = fpressu(1.d4,dh,pop)
         call evoltem(tf,tf,de,dh,prescc,tstep2,0.0d0,0.0d0,luop)
         trec = frectim(tf,de,dh)
         call localem(t,de,dh)
c         call contcont(t,de,dh)
         call totphot(t, dh, fi, rad, dr, dv, wd00, lmod)
         call evoltem(tf, tf, de, dh, prescc, tstep3, exl, tex, luop)
      end if
      call difpop(pop, pop0, trea, atypes, dift)
      t = tf
      end if
c
c
c
 500  continue
c
c     new population is now in pop
c
      write (*,*) ' T :',t
      write (*,*) ' dh:',dh
      write (*,*) ' de :',de
c
      return 
      end
