cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS IIII.  An Astrophysical Plasma Modelling Code.
c
c   Modelling And Prediction in PhotoIonised Nebulae and Gasdynamical Shocks
c   M         A   P             P    I       N           G            S
c
c   Developed 1975-1994 mainly at Mt. Stromlo Stromlo and Siding Spring 
c   Observatories, Institute of Advanced Studies, The Australian National
c   University.
c
c   copyright 1994 Ralph S. Sutherland (1), Michael A. Dopita (1)
c                        Luc Binette (2), and Brent Groves (1)
c
c
c   (1) Research School of Astronomy & Astrophysics
c       Mount Stromlo Observatory
c       The Australian National University
c       Canberra, Australia
c
c       email:
c       bgroves@mso.anu.edu.au
c       ralph@mso.anu.edu.au
c       Michael.Dopita@anu.edu.au
c
c   (2) binette@astroscu.unam.mx
c
c   This code may be freely copied for scientific research and educational
c   purposes.  It may not be copied for commercial purposes without prior
c   consent from the authours.
c   
c   No portion of this source code may be used elsewhere without
c   permission and appropriate acknowledgements.
c
c   Private modifications may be made, but do not distribute modified
c   copies.  If you have modifications that should be distributed, then
c   contact the mappings account (mappings@merlin.anu.edu.au), so that
c   new versions can be made for general distribution.
c
c   When MAPPINGS III is used, please refer to the code and version number
c   in all publications.  The references for MAPPINGS III itself are:
c
c   Dopita, M.A., 1976, Ap. J., 209, 395.
c   Binette, L., Dopita, M.A. & Tuohy, I.R., 1985, Ap.J., 297, 476.
c   Sutherland, R.S., & Dopita, M.A. 1993,ApJS, 88, 253.
c
c   If you intend to use MAPPINGS III often then send email to 
c   ralph@mso.anu.edu.au so that a mailing list for update
c   notices and bug reports can be compiled.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Version 1.0.0r of MAPPINGS III (Sun)
c
c     RSS 11/00
c     BAG edition 04/04 added dust properties
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      program mappings
c
      include 'cblocks.inc'
c
c
      character ilgg*4
      logical initerr
      real tarray1(2)
      integer*4 m,i,ion
c
c     real dtarr,dtime
c     integer*4 iargc,nt0,time
c     character argv*16
c
      theVersion = '1.0.0r'
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     UNIX specific argument passing, used to decide between
c     interactive and batch mode operation.  Will also
c     be use to determine whether to start up the machine
c     interface mode as well...
c
c     modes:  'map3' or 'map3 interactive' use normal interactive mode.
c
c             'map3 batchstart' and 'map3 batchrun' for 
c              batch mode operation.  batchstart just does the 
c              interactive prompting and prepares the input file
c              and batch run uses the input file and diables screen
c              output.
c             
c     
c     
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     iargc() is a Sun UNIX fortan library routine.
c     to get command line argument number
c
c
      runmode = 'interactive'
c      m = iargc()
c      if (m.eq.1) then
c         call getarg(1,argv)
c         if ((argv.eq.'batchstart').or.
c     &       (argv.eq.'batchrun')) then
c            runmode = argv
c         endif
c      endif
c
c
      if (runmode.ne.'batchrun') then
c
      write (*,*) ' ************************************************ '
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) '  Welcome to MAPPINGS III version ',theVersion 
      write (*,*) ' '
      write (*,'(a20,a20)') '          Run mode: ',runmode
      write (*,*) ' '
      write (*,*) ' ************************************************ '
c
      endif
c
c      i = ieee_handler ("set", "common", SIGFPE_IGNORE)
c      CALL nonstandard_arithmetic()
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     call the initialsation routine to read datafiles etc...
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
      initerr = .false.
      call mapinit(initerr)
      if (initerr) goto 1000
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (photonmode.eq.0) then
      if (runmode.ne.'batchrun') then

      write (*,*) ' ***********************************************'
      write (*,*) ' *                                             *'
      write (*,*) ' *   WARNING: PHOTON FIELD DISABLED.           *'
      write (*,*) ' *                                             *'
      write (*,*) ' ***********************************************'
c
      endif
      endif
c
      if (alphacoolmode.eq.1) then
      if (runmode.ne.'batchrun') then

      write (*,*) ' ***********************************************'
      write (*,*) ' *                                             *'
      write (*,*) ' *   WARNING: POWERLAW COOLING ENABLED.        *'
      write (*,*) ' *                                             *'
      write (*,*) ' ***********************************************'
c
      endif
      endif
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     begin interactive initialisation and then subroutine selector...
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***CHANGE ELEMENTAL ABUNDANCES ?
c      FINDS ZGAS RELATIVE TO THE SUN
c
 100  call abecha
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     *** Setup charge exchange reactions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      call chacha
 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Use new Ar and S opacity data from OPRECDAT
c  Comment out to leave true
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      oprec=.FALSE.
c      write(*,*) ' '
c 201  format('Use New Ar and S Opacity data (uncertain) ? (y/N) : ',$)
c      write(*,201)
c 202  format (a)
c      read(*,202) ilgg
c      ilgg = ilgg(1:1)
c      if (ilgg.eq.'n') ilgg = 'N'
c      if (ilgg.eq.'N') oprec = .FALSE.
c  
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     *** Setup dust grain properties
c
c     Now available (averages over all grainsizes)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      call grainpar
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***SELECT A PROGRAM
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (runmode.ne.'batchrun') then
c
      write(*, 310) 
 310  format(///' End of initialisation */*/*/*/*/*/*/*/*/*/*/*/'////)
      endif
 320  if (runmode.ne.'batchrun') then
      write(*, 330) 
 330  format(
     &/' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     & ' Choose a model:'/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')

      write(*, 340) 
 340  format(
     & '    SS  :  Single slab models'/
     & '    CC  :  CIE Cooling curves'/
     & '    PP  :  PIE Photoionization curves'/
     & '        :'/
     & '    S2  :  Shock + diffuse field ( v < 150 kms/s )'/
     & '    S3  :  Shock, elect. + ion Temps ( v < 150 kms/s )'/
     & '    S4  :  Shock, adaptive mesh, requires ext. preionisation'/
     & '        :'/
     & '    P3  :  Photoionisation, dtau, simple continuum'/
     & '    P4  :  Photo. Abs. distance step, full cont.'/
     & '    P5  :  Include Radiation Pressure & Dust'/
     & '        :'/
     & '    H1  :  Spherical Flow model'/
     & '        :'/
     & '    O   :  Old models'/
c     & '    T   :  Tests'/
     & '    R   :  Reinitialise'/
     & '  E,X,Q :  Exit'/
     &       ' :: ',$)
      endif
c
      read(*, 350) ilgg
 350  format(a)
      ilgg = ilgg(1:2)
c
c     fiddle lower case entries..
c
      if (ilgg(1:2) .eq.'rm') ilgg = 'RM'
      if (ilgg(1:2) .eq.'ss') ilgg = 'SS'
      if (ilgg(1:2) .eq.'pp') ilgg = 'PP'
      if (ilgg(1:2) .eq.'cc') ilgg = 'CC'
      if (ilgg(1:1) .eq.'s') ilgg(1:1) ='S'
      if (ilgg(1:1) .eq.'p') ilgg(1:1) ='P'
      if (ilgg(1:1) .eq.'h') ilgg(1:1) ='H'
      if (ilgg(1:1) .eq.'o') ilgg(1:1) ='O'
c      if (ilgg(1:1) .eq.'t') ilgg(1:1) ='T'
      if (ilgg(1:1) .eq.'r') ilgg(1:1) ='R'
      if (ilgg(1:1) .eq.'e') ilgg(1:1) ='E'
c
c     X for exit as well
c
      if (ilgg(1:1) .eq.'X') ilgg(1:1) ='E'
      if (ilgg(1:1) .eq.'x') ilgg(1:1) ='E'
c
c     Q for exit as well
c
      if (ilgg(1:1) .eq.'Q') ilgg(1:1) ='E'
      if (ilgg(1:1) .eq.'q') ilgg(1:1) ='E'
c
c    ***ZERO BUFFER ARRAYS AND RESET COUNTERS
c
      call zer
c
c      dtarr = dtime(tarray1)
c
      runname = 'Test Run'
c
      if (ilgg(1:2) .eq. 'RM') call liftoff
      if (ilgg(1:2) .eq. 'S2') call shock2
      if (ilgg(1:2) .eq. 'S3') call shock3
      if (ilgg(1:2) .eq. 'S4') call shock4
      if (ilgg(1:2) .eq. 'P3') call photo3
      if (ilgg(1:2) .eq. 'P4') call photo4
      if (ilgg(1:2) .eq. 'P5') call photo5
      if (ilgg(1:1) .eq. 'E') goto 2000
      if (ilgg(1:2) .eq. 'H1') call sphere1
      if (ilgg(1:1) .eq. 'E') goto 2000
      if (ilgg(1:1) .eq. 'R') goto 100
c      if (ilgg(1:1) .eq. 'O') call oldmodels
c      if (ilgg(1:1) .eq. 'T') call testsub
      if (ilgg(1:2) .eq. 'CC') call coolc
      if (ilgg(1:2) .eq. 'SS') call slab
      if (ilgg(1:2) .eq. 'PP') call phocrv
c
      if (runmode.ne.'batchrun') then
c
c      dtarr = dtime(tarray1)
c
 400  format('Subroutine time (map2,sys): ',2(f8.1,1x),/) 
      write(*,400) tarray1(1),tarray1(2)
c
      endif
c
      goto 320

      if (runmode.ne.'batchrun') then

      write(*, 390) 
 390  format(/' That code is not yet available : : : : : : : : : :')

      endif

      goto 320

c
c error message
c
 1000 continue
      if (runmode.ne.'batchrun') then
      write(*, 410) 
  410 format(' Ion read was not ion expected:')
      write(*, 420) m, i, ion
  420 format(' Parameters M,I,ION have values :',3i3)
      endif
 2000 continue
c
      stop 
c
      end
