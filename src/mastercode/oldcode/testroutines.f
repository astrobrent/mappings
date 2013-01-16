cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******PROGRAM TO TEST SUBROUTINE IONAB
c
c
      subroutine testa(fn)
c
      include 'const.inc'
c
c
      double precision t,t1,t2
      double precision a(6, 7), rec(5), pion(5)
      double precision pionau(4), ab(6), adndt(6)
      integer*4 i,k,nde
c
      character fn*20
c
c
      open(23, file=fn, status='NEW') 
c
      write(23, 7) 
    7 format(41h TESTA : PROGRAM TO TEST SUBROUTINE IONAB/)
   90 write(*, 100) 
      do 80 k = 1, 7
      do 80 i = 1, 6
   80 a(i,k) = 0.0
  100 format(///34h$NUMBER OF ELEMENTS ON DIAGONAL : )
      read(5, *, err=90) nde
  150 write(*, 160) 
  160 format(42h0LIST OF THE ELEMENTS ON THE UPPER DIAG. ://)
c
      read(5, *, err=150) (a(i,i+1),i = 1, nde-1)
  170 write(*, 180) 
  180 format(42h0LIST OF THE ELEMENTS ON THE LOWER DIAG. ://)
c
      read(5, *, err=170) (a(i,i-1),i = 2, nde)
      if (nde.le.2) goto 173
  163 write(*, 167) 
  167 format(45h0LIST OF THE ELEMENTS BELOW THE LOWER DIAG. ://)
      read(5, *, err=163) (a(i,i-2),i = 3, nde)
c
c
  173 continue
c
  190 write(*, 200) 
  200 format(21h0INITIAL ABUNDANCES :/)
      read(5, *, err=190) (a(i,7),i = 1, nde)
      do 205 k = 1, nde
  205 a(k,k) = 0.0
      do 210 k = 1, nde-1
  210 a(k,k) =-a(k+1,k)
      do 220 k = 2, nde
  220 a(k,k) = a(k,k)-a(k-1,k)
      write(23, 260) 
  260 format(27h0MATRIX OF REACTION RATES :/)
      write(23, 270) ((a(k,i),i = 1, 6),k = 1, nde)
  270 format(1h ,t3,6(1pg13.6))
      do 300 k = 1, nde-1
  300 rec(k) = a(k,k+1)
      do 310 k = 1, nde-1
  310 pion(k) = a(k+1,k)
      do 315 k = 1, nde-2
  315 pionau(k) = a(k+2,k)
      do 320 k = 1, 6
  320 ab(k) = a(k,7)
      t1 = 0.0
      t2 = 0.0
      write(23, 325) 
c
c    ***ITERATION WITH TIME STEP T
c
c
  325 format(37h0IONIC ABUNDANCES AFTER TIME STEP T :/)
c
  400 write(*, 420) 
  420 format(/25h$DURATION OF TIME STEP : )
      read(5, *, err=400) t
      write(23, 450) (a(i,7),i = 1, 6), t1, t2
  450 format(1h ,6(1pg13.6),4h  T:,2(1pg14.6)/)
      t1 = t
      t2 = t2+t
      write(*, 455) t2
  455 format(22h TOTAL ELAPSED TIME : ,1pg13.6)
c
c
      if (t.le.0.0) goto 800
c
c
      call ionab(rec, pion, pionau, ab, adndt, nde, t)
c
      do 470 k = 1, 6
  470 a(k,7) = ab(k)
      write(*, 480) (a(i,7),i = 1, nde)
  480 format(1h ,6(1pg13.6))
      goto 400
  800 if (t .eq. 0.) goto 190
      if (t .eq. (-1.)) goto 90
      close(23) 
      write(*, 1005) 
c
c
 1005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******PROGRAM TO TEST SUBROUTINE SDIFEQ
c
c
c
      subroutine testb(fn)
c
      include 'const.inc'
c
      double precision a, b, c, xhy, aa, bb, cc, xh0
      double precision ano,col,delta,dh,dm
      double precision phi,rb,re,rec,tstep,tt,xht
      integer*4 i,it
      character fn*4
c
c
      open(23, file=fn, status='NEW') 
c
      dm = 1.d37
      write(23, 50) 
   50 format(1h ,44hTESTB : PROGRAM TO TEST SUBROUTINE :  SDIFEQ)
c
   80 write(23, 85) 
   85 format(1h0///)
   90 write(*, 100) 
  100 format(/38h$GIVE RATES : PHOT.,COLLION.,RECOM. : )
      read(5, *, err=90) phi, col, rec
  150 write(*, 160) 
  160 format(/41h$GIVE ABUND. (DH) , INITIAL FRACT. (XHY) ,
     &13hAND DELTA  : )
      read(5, *, err=150) dh, xhy, delta
      tt = 0.0
      it = 1
      xh0 = 1.0-xhy
  190 write(*, 200) 
  200 format(/34h$GIVE TIME STEP AND # OF STEPS  : )
c
c
      read(5, *, err=190) tstep, ano
c
      if (tstep .eq. (-1.0)) goto 150
      if (tstep .eq. (-2.0)) goto 80
      if (tstep.lt.0.0) goto 500
      if ((ano.le.0.0).or.(ano.gt.1.d8)) goto 190
      if (it.gt.1) goto 250
      write(23, 220) 
  220 format(1h0,4x,3hPHI,10x,3hCOL,10x,3hREC,11x,2hDH,10x,3hXHY,8x,
     &5hDELTA)
      write(23, 230) phi, col, rec, dh, xhy, delta
      write(*, 235) phi, col, rec, dh, xhy, delta
  230 format(1h ,6g13.6)
  235 format(1h ,6g11.4)
      write(23, 240) 
      write(*, 240) 
  240 format(1h0,4h   *,6x,4hFHII,9x,3hFHI,12x,4hRES.,12x,
     &12hELAPSED TIME)
c
c
  250 continue
c
      re = dh*col
      rb = dh*rec
      a = phi+(delta*re)
      b = (((1.0-delta)*re)-(delta*rb))-phi
      c =-(rb+re)
      aa =-((a+b)+c)
      bb = b+(2.0*c)
c
c    ***CALL SUBROUTINE SDIFEQ NO TIMES
c
c
      cc =-c
c
c
      do 300 i = 1, ano
c
      call sdifeq(a, b, c, xhy, tstep)
c
c
      call sdifeq(aa, bb, cc, xh0, tstep)
c
  300 continue
      if ((tt.lt.dm).or.(tstep.lt.(dm/ano))) tt = tt+(tstep*
     &ano)
c
      xht = 1.d0-(xhy+xh0)
      write(*, 400) it, xhy, xh0, xht, tt
      write(23, 400) it, xhy, xh0, xht, tt
  400 format(1h ,i4,3x,2(1pg13.6),3x,1pg13.6,4x,1pg13.6)
      it = it+1
c
c
      goto 190
c
  500 close(23) 
      write(*, 1005) 
c
c
 1005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******PROGRAM TO TEST SUBROUTINE GAUNT
c
c
c
      subroutine testc(fn)
c
      include 'const.inc'
c
      double precision e,g
      integer*4 it,iz,j,ndel
      double precision r,t,xsi,y
      double precision akte(50)
c
      character fn*20,ilo*4
c
      double precision fgaunt
c
      open(23, file=fn, status='NEW') 
c
      write(23, 50) 
   50 format(44h TESTC : PROGRAM TO TEST SUBROUTINE : GAUNT ///)
   80 continue
   90 write(*, 100) 
  100 format(///38h GIVE ALL VALUES OF KT/E AT WHICH THE ,
     &33hGAUNT FACTOR WILL BE CALCULATED :///)
  105 do 160 j = 1, 50
  110 read(5, *, err=170) akte(j)
      if (akte(j).le.0.0) goto 200
c
c
  160 continue
c
      goto 200
  170 write(*, 175) 
  175 format(/,'ERROR !!!!   REPEAT ENTRY')
c
c
      goto 105
c
  200 write(*, 220) 
  220 format(//39h$GIVE IONIC SPECIES (EX.: 1=NEUTRAL) : )
      read(5, *, err=200) iz
      if (iz .eq. 0.0) goto 80
      if (iz.lt.0.0) goto 600
  777 write(*, 788) 
  788 format(/44h$ANY CHANGE IN THE PRINCIPAL QUANTUM NUMBER ,8h(Y/N) : 
     &)
      read(5, 798, err=777) ilo
  798 format(a1)
      if ((ilo .ne. 'Y').and.(ilo .ne. 'N')) goto 777
      ndel = 0
c
c
      if (ilo .eq. 'Y') ndel = 1
c
      t = 1.d4
      y = 1.0/iz
      it = t
      write(23, 230) it, iz
      write(*, 230) it, iz
  230 format(//5h T : ,i5,6x,5hIZ : ,i3)
      write(*, 300) 
      write(23, 300) 
c
c
  300 format(/13h GAUNT FACTOR,3x,4hKT/E,7x,4hE/KT,5x,4hNDEL/)
c
      do 500 j = 1, 50
      if (akte(j).le.0.0) goto 200
      e = (1.38054d-16*t)/akte(j)
      r = akte(j)
c
c
      xsi = 1.0/r
c
c
      g = fgaunt(iz,ndel,xsi)
c
      write(23, 330) g, r, xsi, ndel
      write(*, 330) g, r, xsi, ndel
  330 format(1h ,3(1pg11.3),i3)
  500 continue
c
c
      goto 200
c
  600 close(23) 
      write(*, 1005) 
c
c
 1005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE COOL
c       COOLING RATE SUBROUTINE
c
c
c
      subroutine testd(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision cab,de,dh,dift,dr,dv,epotmi,exl
      double precision fi,fron
      integer*4 i,j,l,luop
      double precision prescc,rad
      double precision t,tex,tf,trea,trec,tsl,tst,tstep1
      double precision tstep2,tstep3,wd00,xhyf
      character imod*4, lmod*4, nmod*4,fn*20
      character ilgg*4, jprin*4, jj*4, ndee*4, ill*4
      character where*16
c
c           Functions
c
      double precision feldens,fpressu,frectim
c
      where = 'Test D'
      luop = 23
      cab = 0.0d0
      epotmi = ryde
c
c
      trea = epsilon
c
      open(luop, file=fn, status='NEW') 
 702  write(*, 700) 
 700  format(//' PRINTING OF THE SPECTRUM REQUIRED (Y/N) ?' )
      read(*, 1010) jprin

c
      if ((jprin .ne. 'N').and.(jprin .ne. 'Y')) goto 702
      if (jprin .eq. 'N') write(luop, 900) 
 900  format(40h TESTD : PROGRAM TO TEST SUBROUTINE COOL/)
 93   continue
c
c
c    ***ZERO BUFFER ARRAYS AND DEFINE MODE
c
      call zer
c
      do  i = 1, atypes
         do  j = 1, maxion(i)
            pop(j,i) = 0.0d0
         enddo
      enddo
c
      do  i = 1, atypes
         pop(1,i) = 1.0d0
         if (epot(1,i).lt.ryde) pop(1,i) = 0.0d0
         pop(2,i) = 1.d0-pop(1,i)
      enddo
c
      limph = atypes
      jspot = 'YES'
      jcon = 'YES'
   99 write(*, 97) 
   97 format(//30h SETTING OF IONISATION STATE :/1h ,5x,
     &48hA :  FIXED DEGREE OF IONISATION AT A GIVEN TEMP./1h ,5x,
     &42hB :  EQUILIBRIUM IONISATION AT FIXED TEMP./1h ,5x,
     &48hC :  EQUILIBRIUM IONISATION AT EQUILIBRIUM TEMP./1h ,5x,
     &47hD :  TIME DEPENDANT IONISATION AND EQUIL. TEMP./1h ,5x,
     &45hE :  TIME DEPENDANT IONISATION AT FIXED TEMP./1h ,5x,
     &47hF :  TIME DEPENDANT IONISATION AND TEMPERATURE./1h ,5x,
     &27hO :  NORMAL EXIT AND OUTPUT,t45,4h::? )
      read(*, 1010, err=99) ilgg
c
 1010 format(a)
      l = 1
      if ((ilgg .eq. 'CC').or.(ilgg .eq. 'BB')) l = 2
      if (ilgg .eq. 'BB') ilgg = 'B'
c
      if (ilgg .eq. 'CC') ilgg = 'C'
      if (ilgg .eq. 'O') goto 500
c
      if ((ilgg.lt.'A').or.(ilgg.gt.'F')) goto 99
      if (jprin .eq. 'N') goto 967
  947 write(*, 949) 
  949 format(/41h CASE A,B FOR COLLIS. EXCIT. (0<=X<=1) : )
      read(*, *, err=947) cab
      if ((cab.gt.1.0d0).or.(cab.lt.0.0d0)) goto 947
c
  967 continue
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
      if (ilgg .eq. 'A') then
  103 write(*, 100) 
  100 format(//39h GIVE DEGREE OF IONISATION (0<=D<=5) : )
      read(*, *, err=103) l
      if ((l.lt.0).or.(l.gt.5)) goto 103
      do 33 i = 3, atypes
   33 pop(l+1,i) = 1.0
      if (l .eq. 0) pop(1,1) = 1.0d0
      if (l .eq. 0) pop(1,2) = 1.0d0
      if (l.ge.1) pop(2,1) = 1.0d0
      if (l .eq. 1) pop(2,2) = 1.0d0
c
      if (l.gt.1) pop(3,2) = 1.0d0
  133 write(*, 137) 
  137 format(//44h INDIVIDUAL SETTING FOR ANY ELEMENT (Y/N) : )
      read(*, 1010, err=133) jj
c
      if ((jj .ne. 'Y').and.(jj .ne. 'N')) goto 133
c
      if (jj .eq. 'Y') call popcha(where)
c
  173 write(*, 153) 
  153 format(/33h ELECTRONIC DENSITY SET EQUAL TO ,
     &28hTHE HYDROGEN DENSITY (Y/N) :)
      read(*, 1010, err=173) ndee
      if ((ndee .ne. 'Y').and.(ndee .ne. 'N')) goto 173
      end if
c
c
      call photsou(where)
c
c
      call copinto(pop, pop0)
c
c
      if (ilgg.gt.'B') then
  797 write(*, 787) 
  787 format(/45h NEUTRAL IONIC POPULATIONS OR PREVIOUS VALUES,
     &9h (N/P) : )
      read(*, 1010, err=797) ill
      if ((ill .ne. 'N').and.(ill .ne. 'P')) goto 797
      if (ilgg.lt.'E') write(*, 167) 
  167 format(/47h (NB. A NEGATIVE INITIAL TEMPERATURE ALLOWS ONE,
     &23h TO FOLLOW CONVERGENCE))
      end if
c
      t = 8912.5093813375d0
c
  130 write(*, 150) 
  150 format(/29h GIVE TEMPERATURE , H. DENS. ,
     &28hAND DIL. FACTOR (<=0.5)   : )
      read(*, *, err=130) t, dh, wdil
      if (t.le.10) t = 10**t
c
c
c
c 130  t = 10**(dlog10(t)+0.05)
c      if (t.gt.3.2d8) goto 500
c      dh = 100.d0
c      wdil = 0.d0
c      write(*,*) t,dh,wdil
      if (((wdil.gt.0.5d0).or.(wdil.lt.0.0d0)).or.(dh.le.0.0d0)) 
     &goto 130
      if ((t.le.0.0).and.((ilgg.lt.'C').or.(ilgg .eq. 'E'))) 
     &goto 130
      tf = t
      t = dabs(t)
      de = dh*zion(1)
c
      if ((ilgg .eq. 'A').and.(ndee .eq. 'N')) de = feldens(dh,pop)
      tstep1 = 0.0d0
      tstep2 = 0.0d0
      tstep3 = 0.0d0
c
      if (jprin .eq. 'Y') call zerbuf
      lmod = 'SO'
      imod = 'ALL'
      rad = 1.d38
      dr = 0.0d0
      dv = 0.0d0
      fi = 1.0d0
c
c
      call totphot(t, dh, fi, rad, dr, dv, wdil, lmod)
c
      call zetaeff(fi, dh)
c
      caseab(1) = cab
      caseab(2) = cab
c
      if (ilgg.eq.'B') de = feldens(dh,pop)
      if (ilgg .eq. 'B') call equion(t, de, dh)
      if (de.lt.epsilon) stop
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
  710 write(*, 717) zetae, qhdh
  717 format(/8h ZETAE :,1pg10.3,4x,6hQHDH :,1pg10.3/
     &18h$GIVE TIME STEP : )
      read(*, *, err=710) tstep1
c
      if (tstep1.lt.0.d0) goto 710
      else
      nmod = 'TIM'
  407 write(*, 417) zetae, qhdh
  417 format(/8h ZETAE :,1pg10.3,4x,6hQHDH :,1pg10.3/
     &43h GIVE AGE OF NEBULA AND SOURCE LIFE-TIME : )
      read(*, *, err=407) tst, tsl
      if (tst.le.0.d0) goto 407
      fron = 0.d0
      tstep3 = dmax1(0.d0,tst-tsl)
      tstep1 = fron*dmin1(tsl,tst)
      tstep2 = dmax1(0.0d0,dmin1(tsl,tst)-tstep1)
c
  480 write(*, 485) 
  485 format(/50h DENSITY BEHAVIOR , ISOCHORIC OR ISOBARIC (C/B) : )
      read(*, 1010, err=480) jden
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
      call totphot(t, dh, fi, rad, dr, dv, wd00, lmod)
      call evoltem(tf, tf, de, dh, prescc, tstep3, exl, tex, luop)
      end if
      call difpop(pop, pop0, trea, atypes, dift)
      t = tf
      end if
      write(*, 301) tloss, eloss, egain, dlos, t, dh, de, pop(1
     &,1), wdil
c
  301 format(/14h TOTAL LOSS  :,1pg11.3/14h EFF. LOSS   :,1pg11.3/
     &14h EFF. GAIN   :,1pg11.3/14h FRAC. RESID.:,1pg11.3/
     &14h TEMP.       :,1pg12.5/14h DENSITY     :,1pg11.3/
     &14h ELECTR. DENS:,1pg11.3/14h FR. NEUT. H :,1pg11.3/
     &14h DILUT. F.   :,1pg11.3/)
      if (jprin .eq. 'Y') write(luop, 913) 
  913 format(19h SUBROUTINE COOL : /)
      write(luop, 841) 
  841 format(/9h    TEMP.,t13,5hDENS.,t23,6hEL.DEN,t35,4hDLOS,t46,
     &5hELOSS,t57,5hEGAIN,t65,7hTOT.LOS,t77,4hWDIL,t87,5hTSTAR,t93,3hMOD
     &,t97,5hALPHA,t104,7hTURN-ON,t112,7hCUT-OFF)
      write(luop, 851) t, dh, de, dlos, eloss, egain, tloss, 
     &wdil, teff, iso, alnth, turn, cut
c
  851 format
     &(0pf11.0,2(1pg10.3),1pg11.4,5(1pg10.2),1x,a1,1pg9.2,0pf6.2,0pf7.2)
c
c       STOP
c
      write(*, 300) tloss, hloss, xrloss,rloss, fslos,xiloss,
     & floss, cmplos,fflos, colos, pgain, rngain
  300 format(14h TOTAL LOSS  :,1pg12.5/14h COLEXC      :,1pg12.5/
     &14h XRESON      :,1pg12.5/14h RESON       :,1pg12.5/
     &14h INTER       :,1pg12.5/14h XINTER      :,1pg12.5/
     &14h FORBID      :,1pg12.5/14h COMPTON     :,1pg12.5/
     &14h FREFRE      :,1pg12.5/14h COLION      :,1pg12.5/
     &14h PGAIN       :,1pg12.5/14h RNGAIN      :,1pg12.5)
c
c
      write(luop, 840)
 840     format(9h TOT.LOSS,4x,6hCOLEXC,4x,6hXRESON,4x,6hRESON ,
     &4x,6hINTER ,4x,6hXINTER,4x,6hFORBID,
     &4x,6hCOMPTN,4x,6hFREFRE,4x,6hCOLION,
     &4x,6hPGAIN ,4x,6hRNGAIN)
c
      write(luop, 850) tloss, hloss, xrloss,rloss, fslos,xiloss,
     & floss, cmplos,fflos, colos, pgain, rngain
c
  850 format(1h ,12(1pg10.2))
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
  987 format(//' CASE A,B (H, He): ',f4.2, f4.2)
      dr = 1.0d0
      call sumdata(t, de, dh, fi, dr, dr, dr, imod)
      call avrdata
c      call spectrum(luop, kmod)
      call spec2(luop, 'LAMB','REL')
      end if
c
c
 447  write(*, 444) 
 444  format(/' CONTINUE WITH THE SAME MODE (Y/N) : ')
      read(*, 1010, err=447) jj
      jj = jj(1:1)
      if (jj.eq.'y') jj = 'Y'
      if (jj.eq.'n') jj = 'N'
      if ((jj .ne. 'Y').and.(jj .ne. 'N')) jj = 'N'
      if (jj .eq. 'N') goto 93
c
c
c
c
      goto 130
c
  500 close(luop) 
      write(*, 2005) 
c
c
 2005 format(//38h OUTPUT CREATED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE FORBID,RESON,INTER,HYDREC
c	TO EXIT ENTER TEMP.=0.0
c
c
c
      subroutine testf(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision de,dh,ft
      integer*4 i,j,luop
      double precision t,wei
      character mode*4,fn*20
c
      luop = 23
c
      open(luop, file=fn, status='NEW') 
c
      caseab(1) = 1.0d0
      caseab(2) = 1.0d0
      
      do 100 j = 1, atypes
      zion(j) = 1.0d0
      do 100 i = 1, maxion(j)
c
  100 pop(i,j) = 1.0d0
c
      write(*, 80) 
   80 format(//4h NB.,t6,36hALL ABUNDANCES HAVE BEEN SET TO 1.00/1h ,t6,
     &47hH-BETA FLUX IS FROM PURE RECOMBINATION PROCESS ,8h(CASE B))
  120 write(*, 200) 
  200 format(42h$GIVE TEMPERATURE AND ELECTRONIC DENSITY :)
      read(5, *, err=120) t, de
      dh = 1.0d0
c
      if ((t.le.0.0d0).or.(de.le.0.0d0)) goto 500
      call fivelevel(t, de, dh)
      call inter(t, de, dh)
      call reson(t, de, dh)
      call hydrec(t, de, dh)
c
      tloss = (rloss+fslos)+floss
      wei = (fbri(6,5)+fbri(8,5))*dsign(1.0d0,fbri(10,5))
      roiii = wei/(1.d-38+fbri(10,5))
      deoiii = de
      wei = (fbri(6,2)+fbri(8,2))*dsign(1.0d0,fbri(10,2))
      rnii = wei/(1.d-38+fbri(10,2))
c
      denii = de
c
      call findtde
      write(*, 300) tloss
c
  300 format(20h TOTAL LOSS RATE : :,1pg11.4)
      ft = 4*pi
      do 201 j = 1, nlines
  201 fluxr(j) = rbri(j)*ft
      do 210 j = 1, mlines
  210 fluxi(j) = fsbri(j)*ft
      do 220 j = 1, nfions
      do 230 i = 1, nftrans
  230 fluxf(i,j) = fbri(i,j)*ft
  220 continue
      do 240 j = 1, 10
  240 fluxh(j) = hbri(j)*ft
      fhbeta = hbeta*ft
      fheiilr = heiilr*ft
      vunilog = 0.0d0
c
      mode = 'ABS '
      write(*, *) luop, mode
      call spectrum(luop, mode)
c
      write(*, *) luop, mode
      write(luop, 2400) fhbeta
 2400 format(/////38h EMISSION-LINE SPECTRUM  (ERG/ION/SEC),t67,
     &30hABSOLUTE INTENSITY OF H-BETA :,1pg11.4,12h ERG/ION/SEC)
      write(luop, 2420) t, toiii, tnii, de
 2420 format(24h -----------------------,t73,24hELECTRONIC TEMPERATURE :
     &,g13.6,1hK/8h TOIII :,g13.6,10h    TNII :,g13.6,t77,
     &20hELECTRONIC DENSITY :,g13.6,3h/CC//)
c
c
      goto 120
  500 close(luop) 
      write(*, 1005) 
c
 1005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE PHOTSOU
c
c
c
      subroutine testg(fn)
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision de,eioh
      double precision tnu,wei,weito
      double precision  fl(16, 3)
      integer*4 i,j,ja,jb
      integer*4 k,ki,lc,lf,luop,max,ne
      character caract*80,where*16
      character filna*24,jj*4
      character ilgg*4,fn*20
c
c           Functions
c
      double precision fplank
c
      luop = 23
      where = 'Test G'
c
      eioh = Ryd
  110 write(*, 121) 
  121 format(/,"CREATE PHOTON SOURCE FILE OR TEST SUBR. PHOTSOU ",
     &" (C/T) :")
      read(*, 107, err=110) jj
      if ((jj .ne. 'C').and.(jj .ne. 'T')) goto 110
c
c
c
      if (jj .eq. 'C') goto 700
c
      open(luop, file=fn, status='NEW') 
c
      write(luop, 900) 
c
c
  900 format(43h TESTG : PROGRAM TO TEST SUBROUTINE PHOTSOU//)
c
   50 continue
      call zer
      call photsou(where)
      do 52 i = 1, infph
      tphot(i) = soupho(i)
c
c
c
   52 continue
c
      do 200 i = 1, 16
      de = tedge(i,2)/eioh
      tnu = tedge(i,3)*teff
      fl(i,1) = de
      fl(i,2) = tedge(i,3)
      fl(i,3) = 0.0
      if ((de .eq. 0.0).or.(tnu .eq. 0.0)) goto 200
      fl(i,3) = fplank(tnu,de)
  200 continue
      if (iso .ne. 'C') goto 203
      write(*, 220) (fl(i,1), fl(i,2), fl(i,3),i = 1, 16)
  220 format(1h ,1pg9.3,3x,f7.3,3x,1pg9.2)
c
c
  203 continue
c
      write(luop, 2800) 
 2800 format(/4h MOD,t7,5hTEMP.,t16,5hALPHA,t22,7hTURN-ON,t30,7hCUT-OFF
     &,t38,5hZSTAR,t47,4hFQHI,t56,5hFQHEI,t66,6hFQHEII)
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
      if (iso .ne. 'C') goto 241
      write(luop, 230) (fl(i,1),i = 1, 16)
  230 format(3h EN,16(f8.2))
      write(luop, 235) (fl(i,2),i = 1, 16)
  235 format(3h TE,16(f8.3))
      write(luop, 240) (fl(i,3),i = 1, 16)
  240 format(3h FL,16(1pg8.1))
c
c
  241 continue
c
   74 write(*, 100) 
      read(*, 107, err=74) ilgg
  107 format(a1)
c
c
      if (ilgg .eq. 'Y') goto 50
c
  500 close(23) 
      write(*, 2005) 
 2005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
c   ***TO CREATE PHOTON SOURCE FILE
c
c
      return 
c
  700 continue
      call zer
      do 905 i = 1, infph
  905 tphot(i) = 0.0
      weito = 0.0
      k = 0
  915 continue
  920 write(*, 925) k+1
  925 format(/33h$GIVE RELATIVE WEIGHT OF SOURCE *,i3,12h  (0=EXIT) :)
c
      read(*, *, err=920) wei
      if (wei.lt.0.0) goto 920
c
c
      if (wei .eq. 0.0) goto 930
c
c
      call photsou(where)
c
      do 940 i = 1, infph
  940 tphot(i) = tphot(i)+(wei*soupho(i))
      weito = weito+wei
      k = k+1
      goto 915
  930 write(*, 935) 
  935 format(/31h COMPOSITE SOURCE IS NORMALISED)
      if (k .eq. 0) weito = 1.0
      do 910 i = 1, infph
c
c
c
  910 tphot(i) = tphot(i)/weito
c
   22 write(*, 24) 
   24 format(/,"GIVE PHOTON SOURCE FILENAME : ",$)
      read(*, 27, err=22) lf, filna
c
c
   27 format(a4,a24)
c
c
      open(luop, file=filna, status='UNKNOWN', err=22) 
c
  722 write(*, 712) 
  712 format(/,"GIVE COMMENTS AND TYPE EXIT",/)
  710 read(*, 727, err=722) lc, caract
  727 format(a4,a78)
      ki = index(caract(1:lc),'EXIT')
      if (ki.gt.0) caract(ki:ki+3) = '    '
      write(luop, 730) caract
  730 format(1h ,a78)
c
c
      if (ki .eq. 0) goto 710
c
      i = infph-1
      write(luop, 860) i
  860 format(20h PHOTON SOURCE FILE /
     &44h FILE CREATED ORIGINALLY IN SUBROUTINE TESTG/i6)
      ipho = ipho+1
      ne = 8
      max = ((infph+ne)-1)/ne
      do 120 j = 1, max
      ja = ((j-1)*ne)+1
      jb = j*ne
      if (jb.gt.(infph-1)) jb = infph-1
      write(luop, *) (tphot(i),i = ja, jb)
  120 continue
c
c
      close(luop) 
c
  745 write(*, 100) 
  100 format(/17h$CONTINUE (Y/N) :)
      read(*, 107, err=745) ilgg
c
c
      if (ilgg .eq. 'Y') goto 700
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINES HYDRO AND HYDREC
c
c
c
      subroutine testh(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision de,def,dh,dlo
      double precision t,tstep,xhy,xhyf
      double precision  hrt(7, 18)
      double precision cab
      integer*4 i,j,k,l,luop,m
      character ilo(7)*4, irep*4,fn*20,mod*4
c
      luop = 23
      jspot = 'YES'
      ilo(1) = 'LOGT'
      ilo(2) = ' H-A'
      ilo(3) = ' H-B'
      ilo(4) = ' H-G'
      ilo(5) = ' H-D'
      ilo(6) = ' H-E'
c
c
      ilo(7) = 'LOGB'
c
      open(luop, file=fn, status='NEW') 
c
  810 write(*, 820) 
  820 format(///42h$IONISATION FRACTION SET BY TEMP. (Y/N) ? )
      read(5, 830, err=810) irep
  830 format(a)
      if ((irep .ne. 'Y').and.(irep .ne. 'N')) goto 810
      caseab(1) = 0.0
      caseab(2) = 0.0
      if (irep .eq. 'N') goto 222
  947 write(*, 949) 
  949 format(/41h$CASE A,B FOR COLLIS. EXCIT. (0<=X<=1) : )
      read(5, *, err=947) cab
c
      if ((cab.gt.1.0).or.(cab.lt.0.0)) goto 947
      
      caseab(1) = cab
  222 write(luop, 900) caseab
c
  900 format(41h TESTH : PROGRAM TO TEST SUBROUTINE HYDRO/
     &55h LOGB : LOG OF INTENSITY OF H-BETA PER ION PER ELECTRON,
     &14h  (CASE A,B : ,f4.2,1h))
      if (irep .eq. 'Y') goto 845
      write(*, 860) 
      write(luop, 860) 
  860 format(/29h HYDROGEN COMPLETELY IONISED /)
      goto 850
  845 write(*, 800) 
      write(luop, 800) 
  800 format(/43h THE IONISATION FRACTION OF HYDROGEN IS AT ,
     &10hLEAST 0.1%/38h OTHERWISE FIXED BY COLLIS. IONISATION/)
c
c
  850 continue
      tstep = 1.d36
c
c
      xhy = -1
      do 503 m = 1, 7
      de = 10.d0 ** (m-1)
      dlo = dlog10(de)
c
      dh = 1.0
      write(luop, 2900) dlo
c
 2900 format(//11h LOG(NE) : ,f7.3)
      do 200 k = 1, 18
      do 30 i = 1, atypes
      do 20 j = 1, maxion(i)
   20 pop(j,i) = 0.0d0
   30 pop(1,i) = 1.0d0
      t = 1.d4*(2.d0 ** (k-7))
      hrt(1,k) = dlog10(t)
      pop(1,1) = 0.0d0
      pop(2,1) = 1.0d0
c
      if (irep .eq. 'N') goto 193
      mod = 'EQUI'
      call iohyd(dh, xhy, t, tstep, def, xhyf,mod)
      if (pop(2,1).lt.1.d-3) pop(2,1) = 1.d-3
      if (pop(2,1).le.1.d-2) pop(1,1) = 1.d0-pop(2,1)
c
c
  193 xhyf = pop(2,1)
c
c
      call hydro(t, de, dh)
c
      do 25 l = 1, 5
   25 hrt(l+1,k) = hbri(l)/hbeta
      hrt(7,k) = dlog10((12.566371*hbeta)/((de*xhyf)*dh))
c
c
c
  200 continue
c
      write(luop, 3000) (ilo(l), (hrt(l,j),j = 1, 18),l = 1, 7)
 3000 format(1h ,a4,1h ,18f7.3)
c
c
c
  503 continue
c
  500 close(23) 
      write(*, 2005) 
c
c
 2005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******TO FIND TRAVEL DISTANCE OF PHOTONS TROUGHT NEUTRAL GAS
c	SPACE STEPS DERIVED FROM DTAU
c	VARIABLES ENDING BY "LO" ARE DEXPRESSED IN NEPERIAN LOG
c	VARIABLES ENDING BY "LOG" ARE DEXPRESSED IN BASE 10 LOG
c	VARIABLES CONTAINING "UNI" ARE DEXPRESSED IN RUNITS
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine testi(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision de,dh,dr,dtau,dv,dvoluni,fi,fmdi
      double precision qhlo,rav,rdist
      double precision reclo,rin,rnext,rout,rstep,rstsun,runit,t
      double precision popj0(mxion, mxelem)
      double precision vstrlo
      integer*4 luop,m,nmax
      character lmod*4,imod*4, ilo*4, mmod*4,fn*20,where*16
c
c
c           Functions
c
      double precision fdilu
c
      luop = 21
      limph = 6
      where= 'Test I'
c
c
      jspot = 'YES'
  777 continue
c
c
      call zer
c
      call popcha(where)
c
c    ***INPUT OF PHYSICAL CONDITIONS
c
      call copinto(pop, popj0)
c
      call photsou(where)
c
  810 write(*, 820) 
  820 format(//38h GIVE THE PHOTOIONISING SOURCE RADIUS ,
     &14h(RSUN=6.96E10)/34h$(IN SOLAR UNITS (<1.E8) OR IN CM ,
     &10h(>1.E8) : )
      read(*, *, err=810) rstsun
      if (rstsun.le.0.0) goto 810
      if (rstsun.lt.1.d8) then
      rstar = 6.96d10*rstsun
      else
      rstar = rstsun
      rstsun = rstar/6.96d10
c
      end if
c
  830 write(*, 840) 
  840 format(/39h$GIVE FILLING FACTOR AND PEAK HYDROGEN ,10hDENSITY : )
c
c
      read(*, *, err=830) fi, dh
c
  843 write(*, 847) 
  847 format(/33h STEP VALUE OF THE OPTICAL DEPTH /
     &31h$AND THE DISTANCE MULTIPLIER : )
      read(*, *, err=843) dtau, fmdi
      if ((dtau.le.0.d0).or.(dtau.gt.100.d0)) goto 843
      fmdi = dmax1(1.d0,fmdi)
  888 write(*, 883) 
  883 format(/35h$GIVE MODE FOR DIST. DETERMINATION 11h(LIN/DIS) :)
c
      read(*, 881, err=888) mmod
      if ((mmod .ne. 'LIN').and.(mmod .ne. 'DIS')) goto 888
c
c
  881 format(a)
c
  510 write(*, 512) 
  512 format(/19h$NUMBER OF STEPS : )
      read(*, *, err=510) nmax
c
c
      if (nmax.lt.1) goto 510
c
      reclo = dlog(2.6d-13)
      qhlo = dlog(qht)
      qtlog = (dlog10(qht)+1.09921)+(2.0*dlog10(rstar))
      vstrlo = ((((2.531+qhlo)+(2.0*dlog(rstar)))-(2.0*dlog(dh
     &)))-dlog(fi))-reclo
c
c
      rmax = dexp((vstrlo-1.43241)/3.0)
c
  310 write(*, 320) rmax
  320 format(/36h HYDROGEN STROMGREN SPHERE RADIUS : ,1pg10.3/
     &27h GIVE RADIUS OF EMPTY ZONE /
     &54h$(IN CM (>1.D6) OR IN FRACTION OF STROMGREN RADIUS) : )
      read(*, *, err=310) remp
      if (remp.lt.0.0) goto 310
      if (remp.le.1.d6) remp = remp*rmax
      if (remp.lt.rstar) remp = rstar
c
c
      rmax = dexp(((vstrlo+dlog(1.d0+((remp/rmax) ** 3)))-
     &1.43241)/3.0)
      runit = dmax1(1.d-8*rmax,1.d-5*(rmax-remp))
      vunilog = 3.0*dlog10(runit)
      rin = remp
      lmod = 'SO'
      imod = 'ALL'
      dr = 0.0
      dv = 0.0
      rstep = 0.0
c
c    ***WRITE INITIAL PARAMETERS ON OUPUT FILE
c
c
      t = 5.d3
c
      open(luop, file=fn, status='NEW') 
      write(luop, 2000) 
c
c
 2000 format("TESTI : TO TEST ABSORBTION TROUGH NEUTRAL GAS",//)
c
      write(luop, 2800) 
 2800 format(/,"MOD",7x,"TEMP.",16x,"ALPHA",22x,"TURN-ON",30x,"CUT-OFF",
     &,38x,"ZSTAR",48x,"QHI",57x,"QHEI",67x,"QHEII")
      write(luop, 2850) iso, teff, alnth, turn, cut, zstar, qhi
     &, qhei, qheii
c
c
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
c
      write(luop, 2620) 
 2620 format(/,/,"RSTAR",15x,"REMP",27x,"RSTROMG.",39x,"DTAU",51x,
     & "FILLING F.",63x,"PEAK DENS.",77x,"MMOD",87x,"FMDI")
c
      write(luop, 2622) rstar, remp, rmax, dtau, fi, dh, mmod, 
     & fmdi
c
c
 2622 format(1h ,6(1pg10.3,2x),3x,a4,3x,0pf8.3)
c
      write(*, 2610) 
 2610 format(///3h  *,t5,11hTEMPERATURE,t17,5hRDIST,t28,3hFHI,t39,4hFHII
     &,t50,8hEL.DENS.,t61,4hDLOS,t72,5hRSTEP)
      write(luop, 2627) 
c
c    ***MODEL COMPUTATION STARTS HERE
c
c
 2627 format(///3h  *,t5,11hTEMPERATURE,t17,5hRDIST,t28,3hFHI,t39,4hFHII
     &,t50,8hEL.DENS.,t61,4hDLOS,t72,5hRSTEP,t83,5hZETAE,t94,4hQHDH)
c
      do 400 m = 1, nmax
      rout = rin+rstep
      rav = (rin+rout)/2.0
      rdist = rav
      wdil = fdilu(rstar,rav)
c
c
c
      dvoluni = 4.18879*(((rout/runit) ** 3)-(dble(rin/runit
     &) ** 3))
c
      call totphot(t, dh, fi, rav, dr, dv, wdil, lmod)
      call equion(t, de, dh)
      call taudist(dh, fi, dtau, rnext, rin, popj0, mmod)
c
c
c
      call zetaeff(fi, dh)
c
      write(*, 107) m, t, rdist, pop(1,1), pop(2,1), de, 
     &dlos, rstep
c
      write(luop, 107) m, t, rdist, pop(1,1), pop(2,1), de, 
     &dlos, rstep, zetae, qhdh
c
c
  107 format(1h ,i3,0pf9.0,9(1pg11.3))
c
      rin = rout
c
c
      rstep = fmdi*rnext
c
      call copinto(popj0, pop)
c
c
c
      call sumdata(t, de, dh, fi, dvoluni, rstep, rav, imod)
c
  400 continue
  307 format(a)
  223 write(*, 227) 
  227 format(/18h$CONTINUE (Y/N)? :)
      read(*, 307, err=223) ilo
c
c
      if (ilo .eq. 'Y') goto 777
c
      close(luop) 
      write(*, 2005) 
c
c
 2005 format(//38h OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE FINDTDE
c
c
c
      subroutine testj(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision de,dh,dstep,fi
      double precision pval,t
      integer*4 i,j,luop
      character imod*4,fn*20,mode*4
      character jio(4)*4, jf(3)*4
c
      luop = 23
      jio(1) = 'OIII'
      jio(2) = 'NII'
      jio(3) = 'OII'
      jio(4) = 'SII'
      jf(1) = 'R.'
      jf(2) = 'DE'
c
c
      jf(3) = 'TE'
c
c
c
      open(luop, file=fn, status='NEW') 
c
      write(luop, 900) 
c
c
  900 format(43h TESTJ : PROGRAM TO TEST SUBROUTINE FINDTDE//)
c
      do 120 i = 1, atypes
      do 120 j = 1, maxion(i)
      pval = 11.0d0/dble(maxion(i))
c
c
  120 pop(j,i) = 1.0d0/dble(maxion(i))
c
  137 write(*, 127) 
  127 format(//17h$A : KNOWN T & DE,4x,11hB : UNKNOWN,4x,8hC : EXIT,4x,
     &2h: )
      read(5, 121, err=137) mode
  121 format(a2)
      if ((mode .ne. 'A').and.(mode .ne. 'B')) goto 500
  132 write(*, 135) 
  135 format(/42h$GIVE EL. DENS. FOR TEMPERATURE SENSITIVE ,8hLINES : )
c
      read(5, *, err=132) de
      if (de.le.1.d-6) de = 1.d-6
  140 write(*, 145) 
  145 format(/48h$GIVE TEMPERATURE FOR DENSITY SENSITIVE LINES : )
      read(5, *, err=140) t
      if (t.le.500.) t = 500.
      dh = de
      dstep = 1.0
c
c
      fi = 1.0
c
      if (mode .eq. 'B') goto 200
      imod = 'ALL'
      call cool(t, de, dh)
      call sumdata(t, de, dh, fi, dstep, dstep, dstep, imod)
      call avrdata
c
c
      goto 250
c
  200 write(*, 205) 
  205 format(/31h$GIVE OIII(4959+5007)/4363 AND ,
     &30hNII(6548+6583)/6456  RATIOS : )
      read(5, *, err=200) roiii, rnii
      if (roiii.lt.0.1) roiii = 100
c
c
      if (rnii.lt.0.1) rnii = 100
c
  210 write(*, 215) 
  215 format(/38h$GIVE OII 3729/3726 AND SII 6732/6716 ,9hRATIOS : )
c
      read(5, *, err=210) roii, rsii
      if (roii.lt.0.2) roii = 1.0
c
c
      if (rsii.lt.0.2) rsii = 1.0
c
      toii = t
      tsii = t
      deoiii = de
c
c
      denii = de
c
c
      call findtde
c
  250 continue
      write(*, 267) 
  267 format(/4h ION,t8,5hRATIO,t19,7hEL.DENS,t30,11hTEMPERATURE)
      write(*, 350) jio(1), roiii, deoiii, toiii, jio(2), rnii
     &, denii, tnii, jio(3), roii, deoii, toii, jio(4), rsii, desii, 
     &tsii
c
c
  350 format(1h /4(1h ,a4,1x,1pg11.3,1pg11.3,0pf11.0/))
c
      write(luop, 387) (jf(i),i = 1, 3), (jio(j),j = 1, 4)
  387 format(/1h ,3x,12(a2,a4,4x),5h MODE)
      write(luop, 375) roiii, deoiii, toiii, rnii, denii, tnii
     &, roii, deoii, toii, rsii, desii, tsii, mode
  375 format(/1h ,4(2(1pg10.2),0pf10.0),5x,a2)
c
c
      goto 137
c
  500 close(luop) 
      write(*, 2005) 
c
c
 2005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE RESTRANS
c	(RESONANCE LINE TRANSFER PROBLEM)
c
c
c
      subroutine testk(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision de,dh,dmax,dr,dv,en
      integer*4 i,imax,j,j0,line,luop,nl,nt
      double precision pval,t
      double precision elam(60)
      character fn*20
c
c           Functions
c
      double precision fdismul
c
      luop = 23
c
c
      open(luop, file=fn, status='NEW') 
      write(luop, 900) 
c
  900 format(44h TESTK : PROGRAM TO TEST SUBROUTINE RESTRANS/
     &52h GIVES THE CENTER LAMBDA OF THE CORRESPONDING ENERGY,
     &13h BIN IN EPHOT//)
c
      do 123 i = 1, atypes
      do 123 j = 1, maxion(i)
      pval = 1.0d0/dble(maxion(i))
c
  123 pop(j,i) = 1.0d0/dble(maxion(i))
      do 33 i = 1, nlines
      j = lrbin(i)
      if (j .eq. 0) goto 33
      en = 0.5d0*((ephot(j)+ephot(j+1))*ev)
      elam(i) = 1.9863d-8/en
   33 continue
c
  120 write(*, 150) 
  150 format(//40h$GIVE DV (CM/SEC) AND TEMP. (<0=EXIT) : )
      read(5, *, err=120) dv, t
      if (t.le.0.0) goto 500
  130 write(*, 210) 
  210 format(/26h$GIVE DSTEP AND DENSITY : )
      read(5, *, err=130) dr, dh
      if ((dr.le.0.0).or.(dh.le.0.0)) goto 130
c
      de = dh
c
      call localem(t, de, dh)
c
      do line = 1,xlines
         emilin(2,line) = fdismul(t, dh, dr, dv, line)
      enddo
      dmax = 0.0
      imax = 1
      do 167 i = 1, xlines
            if (emilin(2,i).lt.dmax) goto 167
            dmax = emilin(2,i)
            imax = i
c     
c
c
  167 continue
c
      write(luop, 350) dv, t, dr, dh, dmax, elam(imax)
  350 format(/5h DV :,1pg10.3,5x,3hT :,f11.0,5x,7hDSTEP :,1pg10.3,5x,
     &4hDH :,1pg10.3,5x,6hDMAX :,1pg10.3,4h   @,0pf7.0)
      write(*, 365) dmax, elam(imax)
c
c
  365 format(/7h DMAX :,1pg10.3,4h   @,0pf5.0)
c
      nl = (nlines+9)/10
      do 300 i = 1, nl
      nt = 10
      if ((10*i).gt.nlines) nt = nlines-((i-1)*10)
      j0 = (i-1)*10
      write(luop, 310) (elam(j+1),j = j0, (nt+j0)-1)
  310 format(9h LAMBDA :,t11,10(1pg11.4))
      write(luop, 320) (emilin(2,j+1),j = j0, (nt+j0)-1)
  320 format(9h DISMUL :,t11,10(1pg11.4))
  300 continue
c
c
      goto 120
c
  500 close(luop) 
      write(*, 2005) 
c
c
 2005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE : RATEC
c
c
c
      subroutine testl(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision aa,ey,t
      double precision rate
      integer*4 i,j,k,luop,n,nz
      double precision tel(14), rat(7, 14)
c
      character fn*20
c
      luop = 23
c
      open(luop, file=fn, status='NEW') 
c
      write(luop, 900) 
c
c
  900 format(41h TESTL : PROGRAM TO TEST SUBROUTINE RATEC/)
c
      aa = 0.0
      do 20 i = 1, 10
      aa = aa+2000.
   20 tel(i) = aa
      tel(11) = 5.d4
      tel(12) = 1.0d5
      tel(13) = 1.5d5
c
c
      te(14) = 2.0d5
c
      nz = 1
      do 80 k = 1, 14
      do 80 n = 2, 8
c
c
      t = tel(k)
c
      call ratec(nz, n, t, rate, ey)
c
      rat(n-1,k) = rate
   80 continue
      write(luop, 200) (tel(j),j = 1, 14)
  200 format(///5h TEMP,14(f9.0)//)
      do 27 i = 2, 8
   27 write(luop, 220) i, (rat(i-1,j),j = 1, 14)
c
c
  220 format(1h ,i2,2x,14(1pg9.2)/)
c
  500 close(23) 
      write(*, 2005) 
c
c
 2005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE TAUDIST
c
      subroutine testm(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision dh,drta,ff,fi,frta
      double precision rad
      integer*4 i,ibr,j
      character mmm*4,fn*20,where*16,mmod*4,mod*4
c
      where = 'Test M'
      fi = 1.0
      do 27 i = 1, atypes
      pop(maxion(i),i) = 0.0d0
      ff = maxion(i)-1
      do 27 j = 1, maxion(i)-1
      pop(j,i) = 1.d0/ff
c
   27 continue
c
   45 write(*, 34) 
   34 format(/22h$GIVE MODE (DIS/LIN) :)
      read(5, 37, err=45) mmm
   37 format(a)
      mmod = 'DPRI'
      if ((mmm .eq. 'LI').or.(mmm .eq. 'L')) goto 45
      if ((mmm .eq. 'LIN').or.(mmm .eq. 'LIN ')) mmod = 'LPRI'
c
c
  100 continue
c
c
      call photsou(where)
c
   10 write(*, 20) 
   20 format(/36h$GIVE  DENSITY , TAU  AND  RADIUS : )
      read(5, *, err=10) dh, frta, rad
      if (frta.le.0.0) goto 100
      write(*, 121) 
c
c
  121 format(7h   TAUI,t16,6hTAUNEW,t28,4hDRTA,t40,10hGEO.AMPLF.,t52,
     &6hGAMMAX,t64,1hN)
c
      diluf = 1.d-15
c
c
c
      mod = 'SO'
c
      call totphot(1.d4, dh, fi, 0.d0, 0.d0, 0.d0, diluf, mod)
c
c
      call taudist(dh, fi, frta, drta, rad, pop, mmod)
c
  200 write(*, 110) 
  110 format(/42h$1 NEW PHOTSOURCE   2 NEW TAU   3 END   ? )
      read(5, *, err=200) ibr
      if (ibr .eq. 1) goto 100
c
c
      if (ibr .eq. 2) goto 10
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST DIFFUSE FIELD FROM HYDROGEN
c
c
c
      subroutine testn(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision alpha1,alphaa,de,dh,dista,dr,dva,ene
      double precision fhii, fhi
      double precision fi,frdw
      double precision q1,q2,q3,qemi,qtot,rad
      double precision rphtr,t
      integer*4 i,inl,ist,istio,j
      character fn*20,ilo*4,jmod*4,lmod*4
c
c           Functions
c
      double precision feldens
c
      write(*, 440) 
  440 format(//1h /)
  173 write(*, 177) 
  177 format(/49h$ALL ELEMENTS NEUTRAL OR AS IONISED AS H (N/I) : )
      read(*, 171, err=173) ilo
  171 format(a1)
      if ((ilo .ne. 'I').and.(ilo .ne. 'N')) goto 173
      istio = 2
      if (ilo .eq. 'I') then
  463 write(*, 466) 
  466 format(41h$GIVE STAGE THAT IS POPULATED (1<I<=6) : )
      read(*, *, err=463) istio
      if ((istio.lt.2).or.(istio.gt.6)) goto 463
c
c
      end if
c
  300 call zer
      jspot = 'NO'
      dista = 0.0
      limph = atypes
      do 50 i = 2, atypes
      do 60 j = 2, maxion(i)
   60 pop(j,i) = 0.0d0
c
c
   50 pop(1,i) = 1.0d0
c
  180 write(*, 200) 
  200 format(/50h$GIVE NEUTRAL FR. FOR H ,DENS.,FIL.F. AND TEMP. : )
      read(*, *, err=180) fhi, dh, fi, t
      if (((fhi.lt.0.0).or.(fhi.gt.1.0)).or.(t.lt.0.01)) 
     &return 
c
      if ((dh.le.0.0).or.(fi.le.0.0)) goto 180
      write(*, 402) 
c
c
  402 format(t4,5hDIST.,t14,6hALPHA1,t24,4hQEMI,t34,5hRPHTR,t44,4hQTOT
     &,t54,5hELOSS,t64,5hEGAIN)
c
      fhii = 1.d0-fhi
      pop(2,1) = fhii
      pop(1,1) = fhi
      do 63 i = 2, atypes
      ist = min0(maxion(i),istio)
      pop(ist,i) = 0.0d0
      if (ilo .eq. 'I') pop(ist,i) = fhii
   63 pop(1,i) = 1.0d0-pop(ist,i)
c
c
      de = feldens(dh,pop)
c
      rad = 1.d38
      dva = 0.0
      wdil = 0.5
      lmod = 'DW'
c
c
      jmod = 'OUTW'
c
  183 write(*, 187) 
  187 format(22h$GIVE DISTANCE (CM) : )
      read(*, *, err=183) dr
      if (dr.lt.0.0) goto 300
c
c
      dista = dista+dr
c
      call localem(t, de, dh)
      call totphot(t, dh, fi, rad, dr, dva, wdil, lmod)
      call newdif(t, t, dh, fi, rad, dr, dva, 0.0, dva, frdw, jmod)
      call cool(t, de, dh)
c
c
      call intvec(tphot, q1, q2, q3, qtot)
c
      alphaa = ((de*dh)*fhii)*rec(2,1)
      alpha1 = ((de*dh)*fhii)*(rec(2,1)-rec(3,1))
c
c
      rphtr = (dh*fhi)*rphot(1,1)
c
      do 500 inl = 1, infph-1
      ene = 0.5d0*(ephot(inl)+ephot(inl+1))
c
c
  500 updif(inl) = emidif(inl)*ene
c
      call intvec(updif, q1, q2, q3, qemi)
c
c
      qemi = qemi*ev
c
      write(*, 400) dista, alpha1, qemi, rphtr, qtot, eloss
     &, egain
c
  400 format(1h ,7(1pg10.3))
c
c
      goto 183
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SUBROUTINE COLLION
c
c
      subroutine testo(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision abde,de,dh,frab
      double precision rathb,t,tl
      integer*4 i,j,k,l
      integer*4 luop,ml
      character fn*20
c
      luop = 23
c
      caseab(1) = 0.0
      caseab(2) = 0.0
c
      open(luop, file=fn, status='NEW') 
      write(luop, 900) 
c
  900 format(43h TESTO : PROGRAM TO TEST SUBROUTINE COLLION//)
c
      tl = 1.95d0
      t = 10.d0**tl
c
 220  write(*, 230) 
 230  format(//'Give Atomic Number :')
      read(*, *) i
c
      i = zmap(i)
c
      write(luop,*) 'Rates for :',elem(i)
      write(luop,*) ' '
      if (collmode.eq.0) then
         write(luop,*) 'Calculation method: A&R 1985'
      else
         write(luop,*) 'Calculation method: L&M 1990'
      end if
      write(luop,*) ' '
      write(luop,*) ' '
c
 50   continue
c
      tl = tl+0.05d0
      t = 10.0d0**tl
c
      if (t.le.0) goto 500
c
      dh = 1.0d0
      de = 1.0d0
      pop(1,1) =1.0d0
c
      pop(2,1) = 0.0
c
      call cool(t, de, dh)
c
      abde = ((dh*de)*zion(1))*pop(1,1)
      rathb = 0.0d0
      l = 2
      ml = 8
      do 99 k = l, ml-2
      frab = (caseab(1)*hbemb(l,k))+((1.0-caseab(1))*hbema(l,k))
   99 rathb = rathb+(frab*hherat(k+1))
c
      rathb = rathb/abde
c      write(luop, 2020) t, rathb
c 2020 format(//7h TEMP.:,f11.0,10x,22hEFFECTIVE COLLISIONAL ,
c     &51hEXCITATION RATE FOR H-BETA PER ATOM PER ELECTRON : ,1pg10.3,
c     &9h (CASE A))
c      write(luop, 2043) 
c 2043 format(/t10,37hCOLLISIONAL IONISATION RATES PER ATOM,
c     &15h PER ELECTRON :/)
c
c      call wionpop(luop,col)
c
      write(luop,2000) tl,t,(col(j,i),j = 1,maxion(i))
c
c
 2000 format(f11.2,x,f11.0,2x,29(1pg11.4))
c
c     
      if (tl.ge.8.50d0) goto 500
c
      goto 50

  500 close(luop) 
      write(*, 2005) 
c
 2005 format(//38h0OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SETTING UP OF PREIONISATION CONDITIONS IN SHOCKS
c	USING UPSTREAM COMPONENT OF THE DIFFUSE FIELD
c
c
      subroutine testp(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision d,def,dhpr,dr,dva,epotmi,fi
      double precision qtf,qtot,rad,tepr
      double precision tsmax
      integer*4 i,j,luin,luop
      integer*4 lucty
      character imod*4, lmod*4, kmod*4,jprin*4
      character carac*12,fn*20,where*16
c
      luop = 23
      luin = 24
      rad = 1.d38
      epotmi = ryde
      where = 'Test P'
c
  702 write(*, 700) 
  700 format(//43h$PRINTING OF THE SPECTRUM REQUIRED (Y/N) ? )
      read(*, 1010) jprin
 1010 format(a)
c
c
      if ((jprin .ne. 'N').and.(jprin .ne. 'Y')) goto 702
c
c
c
c
      open(luop, file=fn, status='NEW') 
c
 1700 continue
c
c
      call zer
c
      limph = atypes
c
c
      jspot = 'YES'
c
c
   35 call photsou(where)
c
 2123 write(*, 2135) 
 2135 format(/27h$DILUTION FACTOR (<=0.5) : )
      read(*, *, err=2123) diluf
      if ((diluf.le.0.0).or.(diluf.gt.0.5)) goto 2123
c
c
c
 2103 wdil = diluf
c
  200 write(*, 224) 
  224 format(/34h$GIVE DENSITY ,FILLING FACTOR AND ,
     &27hSHOCK VELOCITY (CM/SEC)  : )
      read(*, *, err=200) dhpr, fi, vshoc
      if ((dhpr.le.0.0).or.(vshoc.le.0.0)) goto 200
      if ((fi.gt.1.0).or.(fi.le.0.0)) goto 200
  203 write(*, 207) 
  207 format(/52h$GIVE MAXIMUM TIMESTEP AND ESTIMATE OF TEMPERATURE :)
c
c
      read(*, *, err=203) tsmax, tepr
c
      write(luop, 900) 
c
c
c
  900 format(39h TESTP : PROGRAM TO TEST SETTING UP OF ,
     &25hPREIONISATION CONDITIONS )
c
      dr = 1.d-4
      dva = 0.0
c
      call popcha(where)
c
      do 1118 i = 3, atypes
      arad(2,i) = dabs(arad(2,i))
      if ((pop(2,i) .eq. 1.0).and.(epot(1,i).lt.epotmi)) arad(2,i)
     & =-arad(2,i)
c
c
 1118 continue
c
      do 55 i = 1, atypes
      do 55 j = 1, maxion(i)
      propop(j,i) = pop(j,i)
      popi(j,i) = pop(j,i)
c
c
   55 popf(j,i) = pop(j,i)
c
      lmod = 'SO'
      imod = 'ALL'
      call totphot(tepr, dhpr, fi, rad, dr, dva, wdil, lmod)
      lucty = 5
      call preion(lucty, luop, tsmax, vshoc, fi, dhpr, def, tepr, qtot, 
     &d)
c
c
      call intvec(soupho, qhi, qhei, qheii, qtf)
c
      qhi = qhi*diluf
      qhei = qhei*diluf
      qheii = qheii*diluf
      qtf = qtf*diluf
      write(luop, 2800) 
 2800 format(//4h MOD,t7,5hTEMP.,t14,5hALPHA,t23,7hCUT-OFF,t34,5hZSTAR
     &,t45,3hQHI,t54,4hQHEI,t64,5hQHEII,t74,5hDILUF)
      write(luop, 2850) iso, teff, alnth, cut, zstar, qhi, qhei
     &, qheii, wdil
 2850 format(1h ,a2,2(1pg9.2),6(1pg10.3))
      write(luop, 2040) 
 2040 format(//36h ABUNDANCES OF THE ELEMENTS RELATIVE,
     &44h TO HYDROGEN AND THEIR STATE OF IONISATION :)
 2043 continue
      write(luop, 2050) (elem(i),i = 1, 11)
 2050 format(1h ,t8,11(5x,a2,4x))
      write(luop, 2060) (zion(i),i = 1, 11)
 2060 format(9h ABUND. :,t8,11(1pg11.3)/)
      carac(1:7) = ' I+++++'
      write(luop, 2070) (carac(1:j+1), (pop(j,i),i = 1, 11),j
     & = 1, 6)
c
c
 2070 format(6(a7,t8,11(1pg11.3)/))
c
      if (jprin .eq. 'N') goto 500
      lmod = 'CAB'
      call totphot(tepr, dhpr, fi, rad, d, dva, wdil, lmod)
      write(luop, 987) caseab(1), tepr
  987 format(/12h CASE A,B : ,f4.2,10x,10hEL.TEMP.: ,f9.0)
      call sumdata(tepr, def, dhpr, fi, d, d, d, imod)
      call avrdata
      kmod = 'REL'
c
c
      call spectrum(luop, kmod)
c
  500 close(luop) 
      write(*, 2005) 
c
c
 2005 format(//38h OUTPUT PRINTED &&&&&&&&&&&&&&&&&&&&&&)
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  TO TEST THE TRUNCATION OF A PHOTOIONIZATION MODEL AT A GIVEN
c  STEP NUMBER (% RECOMBINATION OF HYDROGEN).   THIS READS THE FILE SLOTH.DAT
c  THAT IS WRITTEN OUT BY MAPPINGS.
c
c  AUTHOR : STEPHEN MEATHERINGHAM,    NOVEMBER 1986.
c
      subroutine testq(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision de,dh,dr,dvoluni,fi,t
      integer*4 i,icut,ii
      integer*4 j,luop,lusl,m
      integer*4 ia_sjm(11)
      character newfil*20
      character imod*4,fn*20
      data ia_sjm/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /
c
      call zer
      call zerbuf
c
      luop = 26
      lusl = 45
c
      write(*, 6000) 
 6000 format(7h$FILE> )
      read(5, 5000) newfil
 5000 format(a17)
c
      open(luop, file=fn, status='NEW') 
      open(lusl, file=newfil, status='OLD', recl=500) 
c
      write(*, 6010) 
 6010 format(24h$STEP NUMBER TO CUT AT> )
      read(5, 5010) icut
 5010 format(i4)
c
      read(lusl, 65) (zion(i),i = 1, 11)
   65 format(1h ,11(f10.3))
      do 11 i = 1, atypes
   11 write(*, *) zion(i)
c
      do 9999 ii = 1, icut
      read(lusl,66) m,((pop(j,ia_sjm(i)),j= 1,maxion(ia_sjm(
     &i))),i = 1, 11), t, de, dh, fi, dvoluni, dr, disav
c
   66 format(1h ,i3,66(/,e18.5))
c
      imod = 'ALL'
      call cool(t, de, dh)
      call sumdata(t, de, dh, fi, dvoluni, dr, disav, imod)
c
 9999 continue
c
      call avrdata
      call wrsppop(luop)
c
      close(luop) 
      close(lusl) 
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine testr(fn)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision tstep,ve
      double precision t,de,dh,Bmag,ftime,loss
      character where*16,rmod*4,fn*20
c
c           Functions
c
      double precision feldens
c
       call popcha(where)
c
 100  format(//' Give Flow paramters:'/
     &         ' :::::::::::::::::::::::::::'/
     &         '  Vel, T (K), dh, B(muG)'/
     &         '  loss(erg/cm3/s),tstep (s)'/
     &         ' :: ',$)
c
 110  write(*,100)
      read(*,*) ve,t,dh,Bmag,loss,ftime
c
      if (ve.le.0.d0) goto 200
c
      Bmag = Bmag*1.d-6
c
      if (t.lt.10.d0) t = 10.d0**t
c
      if (tstep.lt.100.d0) tstep = 10.d0**tstep
c
      de = feldens(dh,pop)
c
      rmod = 'NEAR'
c
      write(*,*) ' NEAR solution'
      call  rankhug(t, de, dh, ve, Bmag, rmod,loss,ftime) 
c
      write(*,'(4(1pg12.5,x))') vel0,pr0,rho0,bm0
      write(*,'(4(1pg12.5,x))') vel1,pr1,rho1,bm1
c
c      rmod = 'FAR'
c
c      Bmag = bm0
c      write(*,*) ' FAR solution'
c      call  rankhug(t, de, dh, ve, Bmag, rmod,loss,ftime) 
c
c      write(*,'(4(1pg12.5,x))') vel0,pr0,rho0,bm0
c      write(*,'(4(1pg12.5,x))') vel1,pr1,rho1,bm1
c
      goto 110
c
 200  continue
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******TO TEST SPECIFIED SUBROUTINES
c
c
      subroutine testsub()
c
      include 'cblocks.inc'
c
      character ilgg*4, fn*20,pfx*8,sfx*4
c
   50 write(*, 100) 
  100 format(///'WHICH SUBROUTINE TO TEST :')
c
   70 write(*, 120)
  120 format(/' X',t6,'EXIT  (TO MAIN PROGRAM)'/
     &' a',t6,'IONAB  (IONISATION BALANCE CODE)'/
     &' b',t6,'SDIFEQ  (IONISATION BALANCE OF H)'/
     &' c',t6,'FGAUNT  (GAUNT FACTOR)'/
     &' d',t6,'COOL  (THIN SLAB COOLING)'/
     &' f',t6,'FORBID  (EMISSION LINES)'/
     &' g',t6,'PHOTSOU  (UV IONISATION SOURCE)'/
     &' h',t6,'HYDRO  (HYDROGEN RECOMB. LINES)'/
     &' i',t6,'PROPAGATION  IN NEUTRAL GAS'/
     &' j',t6,'FINDTDE  (GET T & DE FROM LINE RATIOS)'/
     &' k',t6,'RESTRANS  (RESON. LINE ESC. PROBABILITY)'/
     &' l',t6,'RATEC  (COLLIS. EXCIT. OF HI,HEII)'/
     &' m',t6,'TAUDIST  (DISTANCE AT A GIVEN ABSORB.)'/
     &' n',t6,'LOCALEM  (DIFFUSE FIELD FOR HYDROGEN)'/
     &' o',t6,'COLLION  (COLLISIONAL IONISATION)'/
     &' p',t6,'PREION  (PREIONISATION IN SHOCKS)'/
     &' q',t6,'CUT A PHOTOIONISATION OFF AFTER DIFFERENT STEPS)'/
     &' r',t6,'Rankine Hugoniot tests'/
     &t45,'::? ',$)
c
      ilgg = '    '
 130  format (a)
      read(*, 130, err=70) ilgg
      ilgg = ilgg(1:1)
c
      if (ilgg .eq. 'x') ilgg = 'X'
      if (ilgg .eq. 'a') ilgg = 'A'
      if (ilgg .eq. 'b') ilgg = 'B'
      if (ilgg .eq. 'c') ilgg = 'C'
      if (ilgg .eq. 'd') ilgg = 'D'
      if (ilgg .eq. 'f') ilgg = 'F'
      if (ilgg .eq. 'g') ilgg = 'G'
      if (ilgg .eq. 'h') ilgg = 'H'
      if (ilgg .eq. 'i') ilgg = 'I'
      if (ilgg .eq. 'j') ilgg = 'J'
      if (ilgg .eq. 'k') ilgg = 'K'
      if (ilgg .eq. 'l') ilgg = 'L'
      if (ilgg .eq. 'm') ilgg = 'M'
      if (ilgg .eq. 'n') ilgg = 'N'
      if (ilgg .eq. 'o') ilgg = 'O'
      if (ilgg .eq. 'p') ilgg = 'P'
      if (ilgg .eq. 'q') ilgg = 'Q'
      if (ilgg .eq. 'r') ilgg = 'R'
c
      pfx = 'test'//ilgg
      sfx = 'doc'
      fn = ' '
      call newfile(pfx,5,sfx,3,fn)
c

      if (ilgg .eq. 'X') return
      if (ilgg .eq. 'A') call testa(fn)
      if (ilgg .eq. 'B') call testb(fn)
      if (ilgg .eq. 'C') call testc(fn)
      if (ilgg .eq. 'D') call testd(fn)
      if (ilgg .eq. 'F') call testf(fn)
      if (ilgg .eq. 'G') call testg(fn)
      if (ilgg .eq. 'H') call testh(fn)
      if (ilgg .eq. 'I') call testi(fn)
      if (ilgg .eq. 'J') call testj(fn)
      if (ilgg .eq. 'K') call testk(fn)
      if (ilgg .eq. 'L') call testl(fn)
      if (ilgg .eq. 'M') call testm(fn)
      if (ilgg .eq. 'N') call testn(fn)
      if (ilgg .eq. 'O') call testo(fn)
      if (ilgg .eq. 'P') call testp(fn)
      if (ilgg .eq. 'Q') call testq(fn)
      if (ilgg .eq. 'R') call testr(fn)
c
      return
c
      end
