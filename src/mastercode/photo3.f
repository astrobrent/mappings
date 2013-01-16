cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******PHOTOIONISATION MODEL
c	DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c	SPACE STEPS DERIVED FROM CHANGE IN OPTICAL DEPTH : DTAU
c	CHOICE OF EQUILIBRIUM OR FINITE AGE CONDITIONS
c
c	CALL SUBR.: PHOTSOU,COMPPH3,WRIPH3BA
c
c	NB. COMPUTATIONS PERFORMED IN SUBROUTINE COMPPH3
c
c
c
      subroutine photo3()
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision dhn,dht,dthre,dti,epotmi,fin
      double precision qhlo,qhlog,rcin,rechy,reclo,rstsun,vstrlo
c
      integer*4 i,j
c
      character carac*36, where*16
      character caract*4,ilgg*4
      character banfil*12
c
c           Functions
c
      double precision densnum,fdilu
c
      limph = atypes
      jcon = 'NO'
      jspot = 'NO'
      where = 'Photo 3'
      dtau0 = 0.025d0
c
      write(*, 300) 
c
c    ***INITIAL IONISATION CONDITIONS
c
  300 format(///'Photoionisation model P3 selected:'/
     &'Diffuse field computed : Outward integration')
      epotmi = 2.1775d-11
 1010 format(a)
      carac(1:33) = '   '
      do 1120 i = 3, 11
      j = 1+index(carac(1:33),'   ')
      caract = '   '
      if (epot(1,i).lt.epotmi) write(caract, 1010) elem(i)
      carac(j:j+2) = caract // '   '
 1120 continue
      i = j+2
 1110 write(*, 1112) carac(1:i)
      ilgg = '    '
c
 1112 format(//' Allow the elements : ',a/
     &'to recombine freely (Y/N) ?')
c
      read(*, 1010, err=1110) ilgg
      if (ilgg(1:1) .eq. 'Y') goto 1116
      if (ilgg(1:1) .eq. 'y') goto 1116
      if ((ilgg(1:1).ne.'N').and.(ilgg(1:1).ne.'n')) goto 1110
c
      do 1118 i = 3, 11
      arad(2,i) = dabs(arad(2,i))
      if (epot(1,i).lt.epotmi) arad(2,i) =-arad(2,i)
 1118 continue
 1116 continue
c
c     set up ionisation state
c
      where = 'protoionisation'
      call popcha(where)
      where = 'Photo 3'
c
c    ***DENSITY BEHAVIOR
c
  480 write(*, 485) 
  485 format(/' Density behaviour : isoChoric, isoBaric or ',
     &'F(r) (C/B/F) :')
      read(*, 350, err=480) jden
c
      if ((jden(1:1).eq.'C').or.(jden(1:1).eq.'c')) jden = 'C'
      if ((jden(1:1).eq.'B').or.(jden(1:1).eq.'b')) jden = 'B'
      if ((jden(1:1).eq.'F').or.(jden(1:1).eq.'f')) jden = 'F'

      if ((jden.ne.'C').and.(jden.ne.'B').and.(jden.ne.'F'))
     & goto 480
c
      if (jden .eq. 'F') then
 732     write(*, 735) 
 735     format(/' Density Function '/
     &        ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &        ' n(x) = n0*[fx*exp((x-a)/r0) + fp*((x/a)**b)]+c  '/
     &        ' Give parameters n0 fx fp a b c & r0 (cgs): ',$)
         read(*, *) ofac, xfac, pfac, afac, bfac, cfac, scalen
         dhn = ofac
      endif
c
c    ***FILLING FACTOR AND DENSITY
c
      if (jden .eq. 'B') then
      write(*, 911) 
  911 format(/' A filling factor of 1 is considered ',
     &' invariant (complete relaxation)'/
     &' otherwise the product with density is kept invariant.')
      end if
  490 write(*, 495) 
  495 format(/' Give filling factor and nominal hydrogen density:')
      read(*, *, err=490) fin, dhn
      if (((fin.le.0.0d0).or.(fin.gt.1.0d0)).or.(dhn.le.0.0d0)) 
     &goto 490
c
      dht = densnum(dhn)
c
c    ***CHOOSING PHOTON SOURCE AND GEOMETRICAL PARAMETERS
c
   11 call photsou(where)
c
c
      if (qht.le.0.0d0) goto 11
c
      qhlo = dlog(qht)
c
      rechy = 2.6d-13
      reclo = dlog(rechy)
  333 continue
      abc0 = 0d0
      abc9 = 1.d38
  330 write(*, 335) 
  350 format(a)
  335 format(/' Spherical or plane-parallel geometry (S/P) :' )
      read(*, 350) jgeo
c
      if ((jgeo(1:1).eq.'S').or.(jgeo(1:1).eq.'s')) jgeo = 'S'
      if ((jgeo(1:1).eq.'B').or.(jgeo(1:1).eq.'b')) jgeo = 'P'
c
      if ((jgeo .ne. 'S').and.(jgeo .ne. 'P')) goto 333
c
      if (jgeo .eq. 'S') then
c
c    ***IF SYMETRY IS SPHERICAL :
c
  810 write(*, 820) 
  820 format(/' Give the photoionising source radius ',
     &'(Rsun=6.96e10)'/'(in solar units (<1.e8) or in cm (>1.e8):',$)
      read(*, *, err=810) rstsun
      if (rstsun.le.0.0d0) goto 810
      if (rstsun.lt.1.d8) then
      rstar = 6.96d10*rstsun
      else
      rstar = rstsun
      rstsun = rstar/6.96d10
      end if
c
      qtlog = (dlog10(qht)+1.09921d0)+(2.0d0*dlog10(rstar))
      vstrlo = ((((2.531d0+qhlo)+(2.d0*dlog(rstar)))-(2.d0*dlog(
     &dht)))-dlog(fin))-reclo
      rmax = dexp((vstrlo-1.43241d0)/3.d0)
c      write (*,*) qht,dht
      qhdin = (2.d0*qht)/dht
      qhdav = ((4.d0*qht)/dht)*fdilu(rstar,(rmax+rstar)/2.0d0)
c
  312 write(*, 320) rmax, qhdin, qhdav
  320 format(/' (Hydrogen Stromgren radius of the order of',1pg10.3,
     &' cm'/' Impact parameter QHDHIN :',1pg10.3,4x,'QHDHAV :',
     &1pg10.3,' )'/
     &' Give radius of empty zone surrounding source (in cm (>1.e6)'/
     &' or in fraction of Stromgren radius):',$)
      read(*, *, err=312) remp
c
      if (remp.lt.0.0d0) goto 312
      if (remp.le.1.d6) remp = remp*rmax
      if (remp.lt.rstar) remp = rstar
      rmax = dexp(((vstrlo+dlog(1.d0+((remp/rmax) ** 3)))-
     &1.43241d0)/3.0d0)
c
      wdpl = fdilu(rstar,remp)
      qhdin = ((4.0d0*qht)*wdpl)/dht
      wdpl = fdilu(rstar,(rmax+remp)/2.0d0)
      qhdav = ((4.0d0*qht)*wdpl)/dht
c
   60 write(*, 62) 
   62 format(/' Volume integration of spectrum',
     &'over the whole sphere (Y/N) :',$)
      read(*, 350, err=60) ilgg
c
      if ((ilgg(1:1).eq.'Y').or.(ilgg(1:1).eq.'y')) jeq = 'Y'
      if ((ilgg(1:1).eq.'N').or.(ilgg(1:1).eq.'n')) jeq = 'N'
c
      if ((ilgg .ne. 'Y').and.(ilgg .ne. 'N')) goto 60
c
      if (ilgg .eq. 'N') then
 921     write(*, 923) 
 923     format(/' INTEGRATION THROUGH THE LINE OF SIGHT ',
     &        'FOR A CENTERED RING APERTURE'/
     &        'GIVE INNER AND OUTER RADII (CM) :',$)
         read(*, *, err=921) abc0, abc9
         if ((abc0.lt.0.0d0).or.(abc0.ge.(0.999d0*abc9))) goto 921
      end if
c
      endif
c
c    ***IF GEOMETRY PLANE-PARALLEL :
c
      if (jgeo.eq.'P') then
c
         rstar = 0.0d0
         remp = 0.0d0
c
c
         qtlog = dlog10(qht)
         write(*,*) 'photo3 pp:',qht,qhlo,qhlog,dht,fin,reclo
c
         vstrlo = ((qhlo-(2.0d0*dlog(dht)))-dlog(fin))-reclo
c
         rmax = dexp(vstrlo)
c
         qhdav = (2.0d0*qht)/dht
 400     write(*, 410) qhdav
 410     format(/' (IMPACT PARAMETER , QHDHAV :',1pg10.3,' )'/
     &        'GIVE GEOMETRICAL DILUTION FACTOR (<=0.5) :',$)
         read(*, *, err=400) wdpl
c     
         if ((wdpl.le.0.0d0).or.(wdpl.gt.0.5d0)) goto 400
c
         rmax = rmax*wdpl
         qhdin = ((4.0d0*qht)*wdpl)/dht
         qhdav = qhdin
c
      end if
c
      abc3 = 0.0d0
      tm00 = 0.0d0
c
c    ***EQUILIBRIUM OR FINITE AGE ASSUMPTION
c
  460 write(*, 465) 
  465 format(/'Assumption of Equilibrium, Finite age or Post-Equil.',
     &' (E/F/P) :')
      read(*, 350, err=460) jeq
c
      if ((jeq(1:1).eq.'E').or.(jeq(1:1).eq.'e')) jeq = 'E'
      if ((jeq(1:1).eq.'F').or.(jeq(1:1).eq.'f')) jeq = 'F'
      if ((jeq(1:1).eq.'P').or.(jeq(1:1).eq.'p')) jeq = 'P'
c
      if (((jeq .ne. 'E').and.(jeq .ne. 'F')).and.(jeq .ne. 'P')) 
     &goto 460
c
      if (jeq .ne. 'E') then
      rcin = 0.9
      if (jgeo .eq. 'S') then
      dti =-(dlog(1.0d0-(rcin ** 3))/(rechy*dht))
      else
      dti =-(dlog(1.0d0-rcin)/(rechy*dht))
      end if
      dti = (rmax/3.d10)+dmax1(dti,rmax/3.d10)
      dthre = (1.2/dht)/3d-13
c
      if (jeq .eq. 'F') then
  467 write(*, 470) rmax, dti
  470 format(/43h (HYDROGEN STROMGREN RADIUS OF THE ORDER OF,1pg10.3,
     &4h CM,/35h CHARACTERISTIC TIME FOR IONISING :,1pg10.3,5h SEC)/
     &35h GIVE ELAPSED TIME SINCE TURNED ON /
     &40h$AND LIFE-TIME OF PHOTON SOURCE (SEC) : )
      read(*, *, err=467) telap, tlife
      if ((telap.le.0.0d0).or.(tlife.le.0.0d0)) goto 467
      tm00 = 100.0d0
  777 write(*, 772) 
  772 format(/38h$INITIAL TEMPERATURE OF NEUTRAL GAS : )
      read(*, *, err=777) tm00
      if (tm00.lt.1.0d0) goto 777
      else
      tlife = 0.0d0
  431 write(*, 433) rmax, dthre
  433 format(/43h (HYDROGEN STROMGREN RADIUS OF THE ORDER OF,1pg10.3,
     &4h CM,/42h CHARACTERISTIC TIME SCALE FOR RECOMB. OF ,9hHYDROG. :
     &,1pg7.1,5h SEC)/42h GIVE ELAPSED TIME SINCE TURNED OFF (SEC) /
     &51h$AND FRACTIONAL FINAL LUMINOSITY OF SOURCE (>=0) : )
      read(*, *, err=431) telap, abc3
      if ((telap.le.0.0d0).or.(abc3.lt.0.0d0)) goto 431
      end if
      end if
c
  843 write(*, 847) 
  847 format(/' Step value of the optical depth (<=0.02) :')
      read(*, *, err=843) dtau0
      if ((dtau0.le.0.0d0).or.(dtau0.gt.0.1d0)) goto 843
c
c    ***PRINT SETTING UP
c
      if (jeq .eq. 'E') then
      write(*, 625) dtau0
  625 format(/1h ,t5,'SUMMARY * THERMAL AND IONIC EQUILIBRIUM',5x,
     &'DTAU :',0pf7.4)
      else if (jeq .eq. 'F') then
      write(*, 627) telap, tlife, dtau0
  627 format(/t5,'SUMMARY  ;  AGE :,1pg10.3,4x,13hSOURCE LIFE :'
     &,1pg10.3,4h SEC,4x,6hDTAU :,0pf7.4)
      else
      write(*, 623) telap, dtau0
  623 format(/1h ,t5,24hSUMMARY  ;  EQUIL. PLUS ,19hSWITCHING OFF FOR :
     &,1pg10.3,4h SEC,4x,6hDTAU :,0pf7.4)
      end if
c
      write(*, 620) 
  620 format(1h ,t5,5hRSOU.,t15,5hREMP.,t25,4hRMAX,t35,6hDILUAV,t45,
     &7hFILL.F.,t54,5hHDENS,t64,6hQHDHIN,t74,6hQHDHAV)
      write(*, 615) rstar, remp, rmax, wdpl, fin, dhn, qhdin
     &, qhdav
  615 format(1x,4(1pg10.3),0pf9.6,3(1pg10.3))
  503 continue
c
c    ***TYPE OF EXIT FROM THE PROGRAM
c
      ielen = 1
      jpoen = 1
      tend = 10.0d0
      fren = 0.01d0
      diend = 0.0d0
      tauen = 0.0d0
  505 write(*, 510) fren
c
  510 format(/' CHOOSE TYPE OF EXIT FROM PROGRAM :'/
     &'   A  @@  RADIATION BOUNDED NEBULA , H+/H <',2pf5.2,'%'/
     &'   B  @@  IONISATION BOUNDED NEBULA , X+/X < Y'/
     &'   C  @@  TEMPERATURE BOUNDED NEBULA , TMIN'/
     &'   D  @@  OPTICAL DEPTH BOUNDED NEBULA , TAUX'/
     &'   E  @@  DENSITY BOUNDED NEBULA , DIST.'/
     &'   F  @@  (RE-INITIALISE)'/
     &'   G  @@  (RESET GEOMETRY)'/ 
     &,t45,'::',$)
      read(*, 350, err=505) jend
c
      if ((jend(1:1).eq.'A').or.(jend(1:1).eq.'a')) jend = 'A'
      if ((jend(1:1).eq.'B').or.(jend(1:1).eq.'b')) jend = 'B'
      if ((jend(1:1).eq.'C').or.(jend(1:1).eq.'c')) jend = 'C'
      if ((jend(1:1).eq.'D').or.(jend(1:1).eq.'d')) jend = 'D'
      if ((jend(1:1).eq.'E').or.(jend(1:1).eq.'e')) jend = 'E'
      if ((jend(1:1).eq.'F').or.(jend(1:1).eq.'f')) jend = 'F'
      if ((jend(1:1).eq.'G').or.(jend(1:1).eq.'g')) jend = 'G'
c
      if ((jend .ne. 'A').and.(jend .ne. 'B').and.(jend .ne. 
     &'C').and.(jend .ne. 'D').and.(jend .ne. 'E').and.(jend
     & .ne. 'F').and.(jend .ne. 'G')) goto 503
      if (jend .eq. 'F') return 
      if (jend .eq. 'G') goto 333
      if ((jend .eq. 'B').or.(jend .eq. 'D')) then
  525 write(*, 530) 
  530 format(/31h$APPLIES TO ATOMIC ELEMENT # : )
      read(*, *, err=525) ielen
      if ((ielen.lt.1).or.(ielen.gt.11)) goto 525
      jpoen = 1
      if (arad(2,ielen).le.0.0d0) jpoen = 2
      end if
c
      if (jend .eq. 'B') then
  535 write(*, 540) elem(int(ielen))
  540 format(/35h$GIVE FINAL IONISATION FRACTION OF ,a2,3h : )
      read(*, *, err=535) fren
      if ((fren.lt.0.0d0).or.(fren.gt.1.0d0)) goto 535
      end if
c
      if (jend .eq. 'D') then
  545 write(*, 550) elem(int(ielen))
  550 format(/42h$GIVE FINAL OPTICAL DEPTH AT THRESHOLD OF ,a2,3h : )
      read(*, *, err=545) tauen
      if (tauen.le.0.0d0) goto 545
      end if
c
      if (jend .eq. 'C') then
  555 write(*, 560) 
  560 format(/26h$GIVE FINAL TEMPERATURE : )
      read(*, *, err=555) tend
      if (tend.lt.1.0d0) goto 555
      end if
c
      if (jend .eq. 'E') then
  565 write(*, 570) 
  570 format(/50h GIVE DISTANCE OR RADIUS AT WHICH DENSITY DROPS : /
     &51h$(IN CM (>1.D6) OR IN FRACTION OF STROMGREN RADIUS ,
     &10h(<1D6)) : )
      read(*, *, err=565) diend
      if (diend.lt.1.d6) diend = diend*rmax
      if (diend.le.remp) goto 565
      end if
c
c     get runname
c
 100  format (a64)
 101  format(/' Give a name/code for this run: ',$)
      write(*,101)
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
c
c    *** RUN PROGRAM DIRECTLY
c
      banfil = 'INTERACTIVE'
c
c
      call compph3(dhn, fin, banfil)
c
      endif
c
      return 
      end
