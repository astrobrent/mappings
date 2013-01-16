cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******STEADY FLOW PLANE PARALLEL RADIATING SHOCK
c    INCLUDES DIFFUSE FIELD , USES CASE B
c
c      *COMPUTATIONS ARE PERFORMED IN SUBROUTINE COMPSH2
c
c    CALL SUBR.:VELSHOCK,COMPSH2,POPCHA,PHOTSOU
c
c
c
      subroutine shock2()
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision dhpr,epotmi,hmag
      double precision qdl,tepo,tepr
      double precision tmag,vshockm,xhpr
c
      integer*4 i,j
      integer*4 lucty,nstma2
c
      character carac*46
      character ilgg*4, caract*4,ibell*4
      character banfil*12, where*16
c
      lucty = 5
c
      nstma2 = nstmax-50
      jcon = 'YES'
      jspot = 'NO'
      where = 'Shock 2'
      epotmi = ryde
      limph = 6
      do 713 i = limph, atypes
      if (zion(i).gt.0.01) limph = i
  713 continue
c
      ibell(1:4) = char(7)
c
c
 1400 write(*, 1130) 
 1130 format(///' SHOCK MODEL SELECTED : S2'/
     &' PHOTON DIFFUSE FIELD INCLUDED')
c
c     set up ionisation condiditons
c
      where = 'proto-ionisation'
      call popcha(where)
      where = 'Shock 2'
c
      call copinto(pop, propop)
c
      carac(1:33) = '   '
      do 1120 i = 3, atypes
      j = 1+index(carac(1:33),'   ')
      caract = '   '
 1001 format(a2)
      if ((epot(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) then
           write(caract, 1001) elem(i)
      endif
      carac(j:j+2) = caract // '   '
 1120 continue
      i = j+2
c
 1010 format(a2)
 1011 format(a)
c
 1110 write(*, 1112) carac(1:i)
 1112 format(//' Do you allow the elements : ',a/
     &' to recombine freely to neutral state (Y/N) ? ',$)
      read(*, 1011, err=1110) ilgg
      ilgg = ilgg(1:1)
      ilgg = ilgg(1:1)
      if (ilgg .eq. 'y') ilgg = 'Y'
      if (ilgg .eq. 'n') ilgg = 'N'
      if (ilgg .eq. 'Y') goto 1116
      if (ilgg .ne. 'N') goto 1110
      do 1118 i = 3, atypes
      arad(2,i) = dabs(arad(2,i))
      if ((epot(1,i).lt.epotmi).and.(pop(1,i).lt.1d-2)) arad(2,i)
     & =-arad(2,i)
 1118 continue
 1116 continue
c
c
 1018 write(*, 1020) 
 1020 format(//' Select shock conditions :')
 1025 write(*, 1030) (propop(2,1))
 1030 format(//' Give initial hydrogen density and ',
     &'ionised fraction (',1pg14.7,').'/
     &' Fractions less than zero are taken as logs'/
     &' Fractions greater than one cause default to be used'/
     &' :: ',$)
      read(*, *) dhpr, xhpr
      if (dhpr.le.0.d0) dhpr = 10.d0**dhpr
      if (xhpr.gt.1.d0) xhpr = propop(2,1)
      if (xhpr.lt.0.d0) xhpr = 10.d0**xhpr
      if (((xhpr.lt.0.0d0).or.(dhpr.le.0.0d0)).or.(xhpr.gt.1.d0)) 
     &goto 1025
c
c
 1042 write(*, 1040) 
 1040 format(//' Give initial and postshock temperatures : ',$)
      read(*,*) tepr, tepo
      if (tepr.le.10.d0) tepr = 10.d0**tepr
      if (tepo.le.10.d0) tepo = 10.d0**tepo
c
      if ((tepo.le.tepr).or.(tepr.le.100.0d0)) then
      write(*, 29) ibell
   29 format(' UNSUITABLE ANSWER(S) ****************************',a)
      goto 1042
      endif
c
c
 1047 write(*, 1051) 
 1051 format(//' Give magnetic field (micro-Gauss) : ',$)
      read(*, *, err=1047) tmag
      if (tmag.lt.0.0d0) goto 1047
      hmag = tmag*1.0d-6
c
c
c    ***COMPUTES JUMP CONDITION ACROSS SHOCK
c
c
      call velshock(dhpr, xhpr, tepr, tepo, hmag)
      vshockm = vshoc/1.d5
      write(*, 1803) vshockm
 1803 format(//' Shock velocity :',f8.1)
c
c
c    ***PHOTOIONIZATION SOURCES
c
  873 write(*, 877) 
  877 format(//' Any photon source at the shock front',
     &'(Y/N) ? ',$)
      read(*, 1010, err=873) ilgg
      ilgg = ilgg(1:1)
      ilgg = ilgg(1:1)
      if (ilgg .eq. 'y') ilgg = 'Y'
      if (ilgg .eq. 'n') ilgg = 'N'
c
      if ((ilgg .ne. 'Y').and.(ilgg .ne. 'N')) goto 873
      if (ilgg .eq. 'N') goto 2103
c
      call photsou(where)
c
 2123 write(*, 2135) 
 2135 format(/' Dilution factor (<=0.5) : ',$)
      read(*, *, err=2123) diluf
      if ((diluf.lt.0.0).or.(diluf.gt.0.5)) goto 2123
c
c
 2103 continue
c
c
 1300 write(*, 1310) nstma2
 1310 format(//' Give number of steps (MAX:',i4,') '/
     &' and give number of global iterations : ',$)
      read(*, *, err=1300) loop, iter
c
c
      if ((loop.gt.nstma2).or.(loop.lt.8)) goto 1300
      if (iter.le.0) goto 1300
c     
      iel1 = 0
      iel2 = 0
      iel3 = 0
      iel4 = 0
      ieln = 4
 1315  write(*,1320)
 1320 format(//' Choose four elements (atomic numbers)'/
     &' to be followed : ',$)
c
      read (*,*) iel1,iel2,iel3,iel4
      if (iel1.gt.28) goto 1315
      if (iel2.gt.28) goto 1315
      if (iel3.gt.28) goto 1315
      if (iel4.gt.28) goto 1315
c      
c
      if (iel1.le.0) goto 1315
      if (iel2.le.0) goto 1315
      if (iel3.le.0) goto 1315
      if (iel4.le.0) goto 1315
c
      iel1 = zmap(iel1)
      iel2 = zmap(iel2)
      iel3 = zmap(iel3)
      iel4 = zmap(iel4)
c      
      if (iel1.eq.0) goto 1315
      if (iel2.eq.0) goto 1315
      if (iel3.eq.0) goto 1315
      if (iel4.eq.0) goto 1315
c
c    ***SET TYPE OF ENDING FOR THE PROGRAM
c
      xhmin = 0.01d0
      qdl = 1.5d0+((0.2d0*(tepo/1.d5))*(nstma2/(loop)))
      write(*, 1800) vshockm, qdl
 1800 format(//' SHOCK VELOCITY :',f8.1,' KM/SEC',7x,'Q FINE-TUNING :'
     &,f7.2)
 1810 write(*, 1820) xhmin, elem(iel1)
 1820 format(//' CHOOSE TYPE OF ENDING FOR THE PROGRAM :'/t6,
     &' A-  NORMAL ENDING  (FINAL H+/H :',2pf4.1,'% )'/t6,
     &' B-  ACCORDING TO FINAL TEMPERATURE'/t6,
     &' C-  ACCORDING TO ELAPSED TIME'/t6,
     &' D-  ACCORDING TO THE IONISED FRACTION OF : ',a2/t6,
     &' E-  ACCORDING TO LINE INTENSITIES'/t6,' R-  RE-INITIALIZE'/' ' 
     &,t55,':: ',$)
      read(*, 1821, err=1810) jfin
 1821 format(a)
      jfin = jfin(1:1)
c
      if (jfin.eq.'a') jfin = 'A'
      if (jfin.eq.'b') jfin = 'B'
      if (jfin.eq.'c') jfin = 'C'
      if (jfin.eq.'d') jfin = 'D'
      if (jfin.eq.'e') jfin = 'E'
      if (jfin.eq.'r') jfin = 'R'
c
      if ((jfin.lt.'A').or.((jfin.gt.'E').and.(jfin .ne. 'R'))) 
     &goto 1810
      if (jfin .eq. 'R') return 
      if (jfin .eq. 'A') goto 1895
      if (jfin .ne. 'D') goto 1877
c
 1830 write(*, 1840) 
 1840 format(/' WHICH IONISING STAGE : ',$)
      read(*, *, err=1810) istf
      if ((istf.lt.2).or.(istf.gt.maxion(iel1))) goto 1810
c
 1877 continue
c
      if (jfin .ne. 'E') goto 1838
 1863 write(*, 1868) 
 1868 format(/' CHOOSE LINE NUMBER : ' ,
     &' [OIII]   2 [OII]   3 [NII]   4 [SII]  :: ',$)
      read(*, *, err=1810) istf
      if ((istf.lt.1).or.(istf.gt.4)) goto 1810
 1838 continue
c
 1823 write(*, 1825) 
 1825 format(/' GIVE FINAL VALUE : ')
      read(*, *, err=1810) fval
      if (fval.lt.0.d0) fval = 10**fval
c
      if (jfin .eq. 'B') texi = fval
 1895 continue
c
      jpre = 'C'
      if (iter.le.1) goto 1857
 1855 write(*, 1853) 
 1853 format(/' DETERMINATION OF THE SUCCESSIVE ',
     &'PREIONISATION CONDITIONS (U/F/C) :'/t6,
     &' U-  USING UPSTREAM PHOTON FIELD'/t6,
     &' F-  USING CONDITIONS AT THE END OF SHOCK'/t6,
     &' C-  CONSTANT FROM ONE ITERATION TO THE OTHER'/t6,
     &' T-  THERMAL EQUILIBRIUM AT PRESHOCK TEMP.'/,t55,':: ',$)
      read(*, 1011, err=1855) jpre
      jpre = jpre(1:1)
      if (jpre.eq.'u') jpre = 'U'
      if (jpre.eq.'f') jpre = 'F'
      if (jpre.eq.'c') jpre = 'C'
      if (jpre.eq.'t') jpre = 'T'
      if ((((jpre .ne. 'U').and.(jpre .ne. 'F')).and.(jpre .ne. 'C')
     &).and.(jpre .ne. 'T')) goto 1855
c
 1857 continue
c
c     get runname
c
 100  format (a64)
 2000 format(//' Give a name/code for this run: ',$)
      write(*,2000)
      read (*,100) runname
c
c     current UNIX batch system...
c
      if (runmode.ne.'batchstart') then
c
c
c     remains of old vax batch system (not used)
c
c
      banfil = 'Interactive'
c
      call compsh2(lucty, dhpr, xhpr, tepr, tepo, hmag, banfil)
c
      endif
c
      return 
      end
