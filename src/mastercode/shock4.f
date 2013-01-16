cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     SHOCK4: Shock model with explicit Rankine-Hugoiot solution for
c     each step.  Time steps based on fastest atomic, cooling or 
c     photon absorbion timescales.
c
c     Full continuum and diffuse field, requires precalculated ionisation
c     input.
c
c     RSS 1992
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine Shock4()
c
c
      include 'cblocks.inc'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           Variables
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision dhpr,dr,drdw,drup,dvdw,dvup,frdw,ftime
      double precision tdw,tepo,tepr,tlim,tmag,tup
      double precision vshockm,xhii,xhpr
      double precision t, de, dh, dv, rad, fi, wdilt0
      double precision tt0,l0,l11,l12,l1,cmpf,ue0,ue1
      double precision r0,p0,v0, tnloss,cspd,wmol,mb
      double precision Hfree, Bmag, ve, en, ue,rdis,dva
      double precision pop1(mxion, mxelem),totn,ionbar,abstime
      double precision press, mu, rhotot,cltime,utime 
      double precision ptime, ctime, rtime, htime,stime,eqtime
      double precision atimefrac,cltimefrac,abtimefrac
      double precision timend,teq,dheq,tl0,tl1,tt1
      double precision band1,band2,band3,band4,band5,band6,bandall, wid
c
      integer*4 ie,luop,lusp
      integer*4 count, step, np,i,j,lurt,ludy,lual,m, lupf
      integer*4 luions(4)
c
c     time steps - non-standard so commented out
c
c      integer*4 nt0,nt1,time
c      real tarray(2)
c     real dtime,dtarr
c
      character filn(4)*20
      character specmode*4, rmod*4, imod*4, jmod*4,tsrmod*4,spmod*4
      character fn*20, fl*20, lmod*4, tab*1, ilgg*4, tmod*4, linemod*4
      character fd*20, fr*20, fa*20, fsp*20, allmod*4, ratmod*4,ispo*4
      character pfx*8, sfx*4, caller*4,wmod*4,pollfile*12,vmod*4
      character cht*12,pht*12,clt*12,dynmod*4, mypfx*8,where*16
      character fpf*20
c
      logical iexi,pseudo
c
c           Functions
c
      double precision fcolltim,feldens,fphotim,fpressu,frectim2,frho
c
c     set up logical unit numbers
c
      lual = 20
      luop = 21
      lurt = 22
      ludy = 23
      lusp = 24
      lupf = 25
c
c      nt1 = time()
c      nt0 = nt1
c
      m = 0
c
      ieln = 4
      do i = 1,ieln
         luions(i) = 23+i
      enddo
c
      fl = ' '
      pfx = 'shckn'
      sfx = 'sh4'
      call newfile(pfx,4,sfx,3,fl)
c     
      open(luop, file=fl, status='NEW') 
c
c     zero arrays and define modes: all ions, no on-thespot,
c     full continuum
c     calculations..
c
      call zer
      call zerbuf
      tab = char(9)
      limph = atypes
      jspot = 'NO'
      jcon = 'YES'
      jspec = 'NO'
      ispo = 'SII'
      allmod = 'NO'
      tsrmod = 'NO'
      ratmod = 'NO'
      dynmod = 'NO'
      vmod = 'NONE'
      mypfx = 'psend'
      specmode = 'DW'
      jgeo = 'P'
      jden = 'B'
      pseudo = .FALSE.
      where = 'Shock 4'
 10   format(a)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTRO - Define Input Parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 100  format(///,
     & ' SHOCK 4 model selected:',/,
     & ' Shock code with natural coordinates.',/,
     & ' ****************************************************')
      if (runmode.ne.'batchrun') then
        write(*,100)
      endif
c
c     get ionisation balance to work with...
c
      where = 'pre-ionisation'
      call popcha(where)
      where = 'Shock 4'
c
c     set up current field conditions
c
      call photsou(where)
c
c     Shock model preferences
c
c
 110  format(///,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Setting the shock conditions',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::')
      if (runmode.ne.'batchrun') then
        write(*,110)
      endif
c
 115  format(//,
     & ' Choose Shock Jump Paramter:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    T : Shock jump in terms of Temperature jump',/,
     & '    V : Shock jump in terms of flow Velocity.',//,
     & ' :: ')
c
 116  if (runmode.ne.'batchrun') then
        write(*,115)
      endif
      read(*,10) ilgg
      ilgg = ilgg(1:1)
      if (ilgg.eq.'T') ilgg = 't'
      if (ilgg.eq.'V') ilgg = 'v'
c
      if ((ilgg.ne.'t').and.(ilgg.ne.'v')) goto 116
c
c     Shock in term of flow velocity
c
      if (ilgg.eq.'v') then
 120    format(//,' Give Preshock Conditions:',/,
     &  '  T (K), dh (N), v (cm/s), B, Geo. dilution',/,
     &  ' (T<10 taken as a log, dh <= 0 taken as a log),',/,
     &  ' (B in micro Gauss, Dilution factor <= 0.5)',/,
     &  ' : ')
        if (runmode.ne.'batchrun') then
         write(*,120)
        endif
        read(*, *) t, dh, ve, bmag, wdil
c
        dr = 1.d0
c
        if (t.le.10.d0) t = 10.d0**t
        if (dh.le.0.d0) dh = 10.d0**dh
        if (wdil.gt.0.5d0) wdil = 0.5d0
        wdilt0 = 0.d0
c
        Bmag = Bmag*1.d-6
        vshoc = ve
        de = feldens(dh,pop)
c
c     fill in other default parameters
c
      ftime = 0.d0
      tloss = 0.d0
      rmod = 'NEAR'
      call rankhug(t, de, dh, ve, Bmag, rmod,tloss,ftime)
      tm00 = t
c
      endif
c
c     Shock in terms of temp jump
c
      if (ilgg.eq.'t') then
 125    format(//,' Give Preshock Conditions:',/,
     &  '  (T<10 taken as a log, dh <= 0 taken as a log)',/,
     &  '  (B in micro Gauss, Dilution factor <= 0.5)',/,
     &  ' : T (K), dh (N), B, Geo. dilution',/,
     &  ' : ')
        if (runmode.ne.'batchrun') then
         write(*,125)
        endif
        read(*, *) te0, dh, Bmag, wdil
c
        if (dr.le.0.d0) dr = 1.d0
        if (wdil.gt.0.5d0) wdil = 0.5d0
        Bmag = Bmag*1.d-6
c
 126    format(//,' Give Postshock Conditions:',/,
     &  '  (T<10 taken as a log)',//,
     &  ' : T (K)',/,
     &  ' : ')
        if (runmode.ne.'batchrun') then
         write(*,126)
        endif
        read(*, *) te1
c
        if (te0.le.10.d0) te0 = 10.d0**te0
        if (te1.le.10.d0) te1 = 10.d0**te1
c
        bm0 = Bmag
c
        call velshock(dh, pop(2,1), te0, te1, Bmag)
c
        de   = feldens(dh,pop)
        vel0 = vshoc
        vel1 = vpo
        rho0 = frho(de,dh)
        rho1 = frho(depo,dhpo)
        pr0  = fpressu(te0,dh,pop)
        pr1  = fpressu(te1,dhpo,pop)
        bm1  = bm0*rho1/rho0 
        tm00 = te0
c
      endif
c
      if (alphacoolmode.eq.1) then
 600   format(//,
     &  ' Powerlaw Cooling enabled  Lambda T6 ^ alpha :',/,
     &  ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     &  '    Give index alpha, Lambda at 1e6K',/,
     &  '    (Lambda<0 as log)',/,
     &  ' :: ')
       if (runmode.ne.'batchrun') then
         write(*,600)
       endif
       read(*,*) alphaClaw, alphaC0
       if (alphaC0.lt.0.d0) alphaC0 = 10.d0**alphaC0
      endif
c
c     Diffuse field interaction
c
c
      photonmode = 1
c
 230  format(//,
     & ' Choose Diffuse Field Option :',/,
     & ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    F  : Full diffuse field interaction (Default).',/,
     & '    Z  : Zero diffuse field interaction.',/,
     & '       :',/,
     & '    I  : fInite extent plane shock (approx. rad. transfer)',//,
     & ' :: ')
      if (runmode.ne.'batchrun') then
        write(*,230)
      endif
 231  read(*,10) ilgg
      ilgg = ilgg(1:1)
c
      if (ilgg.eq.'Z') ilgg = 'z'
      if (ilgg.eq.'F') ilgg = 'f'
      if (ilgg.eq.'I') ilgg = 'i'
      if ((ilgg.ne.'z').and.(ilgg.ne.'f').and.(ilgg.ne.'i')) goto 231
c
      if (ilgg.eq.'z') photonmode = 0
      if (ilgg.eq.'f') photonmode = 1
      if (ilgg.eq.'i') then
        photonmode = 1
        jgeo = 'F'
 232    format(/' Give clyindrical shock front radius (cm): ',$)
        read(*,*) rshock
        if (rshock.le.0) rshock = 0.d0
        do i = 1,infph 
           dilradpho(i) = 0.d0
           difftot(i) = 0.d0
        enddo
      endif
c
 111  format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Calculation Settings ',/,
     & ' ,::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      if (runmode.ne.'batchrun') then
       write(*,111)
      endif
c
c     Set up time step limits
c
 127  format(//,' Choose a Time Step Behaviour :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Auto. time steps, based on Atomic timescales.',/,
     & '    M  : Choose a Maximum time step, Auto otherwise.',/,
     & '    N  : Choose a miNimum time step, Auto otherwise.',/,
     & '    S  : Strict preset timestep. (not recomended)',//,
     & ' :: ')
c
 128  if (runmode.ne.'batchrun') then
        write(*,127)
      endif
      read(*,10) tmod
      tmod = tmod(1:1)
      if (tmod.eq.'A') tmod = 'a'
      if (tmod.eq.'M') tmod = 'm'
      if (tmod.eq.'N') tmod = 'n'
      if (tmod.eq.'S') tmod = 's'
c
      if ((tmod.ne.'a').and.(tmod.ne.'m').and.
     & (tmod.ne.'n').and.(tmod.ne.'s')) goto 128
c
      atimefrac = 0.2d0
      cltimefrac = 0.2d0
      abtimefrac = 0.2d0
c
      if (tmod.ne.'s') then
 117     format(//,
     &   ' Give timescale fractions:',/,
     &   ' ( 0 < dtau < 1 )',/,
     &   ' ( recommend: 0.05 < dtau < 0.25 )',/,
     &   ' : atomic, cooling, absorbsion',/,
     &   ' : ')
 118     if (runmode.ne.'batchrun') then
            write(*,117)
         endif
         read(*, *) atimefrac,cltimefrac,abtimefrac
c
         if (atimefrac.le.0.d0) goto 118
         if (cltimefrac.le.0.d0) goto 118
         if (abtimefrac.le.0.d0) goto 118
c
         if (atimefrac.ge.1.d0) goto 118
         if (cltimefrac.ge.1.d0) goto 118
         if (abtimefrac.ge.1.d0) goto 118
c         
      endif
c
      utime = 0.d0
      if (tmod.ne.'a') then
 129    format(//,
     &  ' Give the time step:',/,
     &  ' ( time < 100 taken as a log)',//,
     &  ' : time (sec)',/,
     &  ' : ')
        if (runmode.ne.'batchrun') then
         write(*,129)
        endif
        read(*, *) utime
        if (utime.le.100.d0) utime = 10.d0**utime
      endif
c
c
 112  format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Boundry Conditions ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      if (runmode.ne.'batchrun') then
        write(*,112)
      endif
c
c     Choose ending
c
 136  format(/,
     & ' Choose Ending :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard ending, 1% weighted ionisation.',/,
     & '    B  : Species Ionisation limit.',/,
     & '    C  : Temperature limit.',/,
     & '    D  : Distance limit.',/,
     & '    E  : Time limit.',/,
     & '    F  : Thermal Balance limit.',/,
     & '    G  : Heating limit.',/,
     & ' :: ')
c
 137  if (runmode.ne.'batchrun') then
         write(*,136)
      endif
      read(*,10) jend
      jend = jend(1:1)
c
      if (jend.eq.'A') jend = 'a'
      if (jend.eq.'B') jend = 'b'
      if (jend.eq.'C') jend = 'c'
      if (jend.eq.'D') jend = 'd'
      if (jend.eq.'E') jend = 'e'
      if (jend.eq.'F') jend = 'f'
      if (jend.eq.'G') jend = 'g'
c
      if (
     & (jend.ne.'a').and.
     & (jend.ne.'b').and.
     & (jend.ne.'c').and.
     & (jend.ne.'d').and.
     & (jend.ne.'e').and.
     & (jend.ne.'f').and.
     & (jend.ne.'g')) goto 137
c
c     Secondary info:
c
      if (jend.eq.'b') then
 400    format(//,' Give Atom, Ion and limit fraction: ')
        if (runmode.ne.'batchrun') then
         write(*,400)
        endif
        read(*,*) ielen,jpoen,fren
      endif
      if (jend.eq.'c') then
 410    format(//,' Give final temperature (K > 10, log <= 10): ')
        if (runmode.ne.'batchrun') then
         write(*,410)
        endif
        read(*,*) tend
        if (tend.le.10.d0) tend = 10.d0**tend
      endif
      if (jend.eq.'d') then
 420    format(//,' Give final distance (cm > 100, log<=100): ')
        if (runmode.ne.'batchrun') then
         write(*,420)
        endif
        read(*,*) diend
        if (diend.le.100.d0) diend = 10.d0**diend
      endif
      if (jend.eq.'e') then
 430    format(//,' Give time limit (s > 100, log<=100): ')
        if (runmode.ne.'batchrun') then
         write(*,430)
        endif
        read(*,*) timend
        if (timend.le.100.d0) timend = 10.d0**timend
      endif
c
 113  format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' Output Requirements ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      if (runmode.ne.'batchrun') then
      write(*,113)
      endif
c
c     default monitor elements
c
      iel1 = zmap(1)
      iel2 = zmap(2)
      iel3 = zmap(6)
      iel4 = zmap(8)
      mypfx = 'psend'
c
  440 format(//,' Choose output settings : ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard output. ',/,
     & '    B  : Standard plus ion balance files.',/,
     & '    C  : Standard plus dynamics file.',/,
     & '    D  : Standard plus rates file.',/,
     & '    E  : Standard plus final down stream field.',/,
     & '    F  : Standard plus upstream field at each step.',/,
     & '       :',/,
     & '    G  : B + E',/,
     & '    H  : B + F',/,
     & '       :',/,
     & '    I  : Everything',//,
     & ' :: ')
 2101 if (runmode.ne.'batchrun') then
         write(*, 440)
      endif
      read(*, 10) ilgg
      ilgg = ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg = 'A'
      if (ilgg.eq.'b') ilgg = 'B'
      if (ilgg.eq.'c') ilgg = 'C'
      if (ilgg.eq.'d') ilgg = 'D'
      if (ilgg.eq.'e') ilgg = 'E'
      if (ilgg.eq.'f') ilgg = 'F'
      if (ilgg.eq.'g') ilgg = 'G'
      if (ilgg.eq.'h') ilgg = 'H'
      if (ilgg.eq.'i') ilgg = 'I'
c
      if ((ilgg.ne.'A').and.(ilgg.ne.'B').and.(ilgg.ne.'C').and.
     & (ilgg.ne.'D').and.(ilgg.ne.'E').and.(ilgg.ne.'F').and.
     & (ilgg.ne.'G').and.(ilgg.ne.'H').and.(ilgg.ne.'I')) goto 2101
c
c     Upstream fields
c
      if ((ilgg.eq.'F').or.(ilgg.eq.'H').or.(ilgg.eq.'I')) then
         lmod = 'y'
      endif
c
c     F-Lambda plots, ionising and optical
c
      if ((ilgg.eq.'E').or.(ilgg.eq.'G').or.(ilgg.eq.'I')) then
         jspec= 'YES'
      endif
c
c     timescales and rates file
c
      ratmod = 'n'
      if ((ilgg.eq.'D').or.(ilgg.eq.'I')) then
         ratmod = 'y'
         pfx = 'rates'
         sfx = 'sh4'
         call newfile(pfx,5,sfx,3,fr)
         open(lurt, file=fr, status='NEW')
      endif
c
c     flow dynamics file
c
      dynmod ='n'
      if ((ilgg.eq.'C').or.(ilgg.eq.'I')) then
         dynmod = 'y'
         pfx = 'dyn'
         sfx = 'sh4'
         call newfile(pfx,3,sfx,3,fd)     
         open(ludy, file=fd, status='NEW') 
      endif
c
      pfx = 'spec'
      sfx = 'sh4'
      call newfile(pfx,4,sfx,3,fsp)
      open(lusp, file=fsp, status='NEW') 
c
c     ionbalance files
c
      allmod = 'n'
      tsrmod = 'n'
c
      if ((ilgg.eq.'B').or.(ilgg.eq.'G').or.(ilgg.eq.'H')
     &     .or.(ilgg.eq.'I')) then
c
       if ((ilgg.eq.'B').or.(ilgg.eq.'G').or.(ilgg.eq.'H')) then
c
 234     format(//,' Record all ions file (Y/N)?: ')
         if (runmode.ne.'batchrun') then
          write(*,234)
         endif
         read(*,10) allmod
         allmod = allmod(1:1)
         if (allmod.eq.'Y') allmod = 'y'
       endif
c
       if (ilgg.eq.'I') allmod = 'y'
c
       if (allmod.eq.'y') then
         pfx = 'allion'
         sfx = 'sh4'
         call newfile(pfx,6,sfx,3,fa)
         open(lual, file=fa, status='NEW')      
       endif
c
       if (ilgg.eq.'B') then
 138     format(//,' Record particular ion balances (Y/N)?: ')
         if (runmode.ne.'batchrun') then
           write(*,138)
         endif
         read(*,10) tsrmod
         tsrmod = tsrmod(1:1)     
         if (tsrmod.eq.'Y') tsrmod = 'y'
       endif
c
       if (ilgg.eq.'I') tsrmod = 'y'
c
       if (tsrmod.eq.'y') then
 146     format(//' Give 4 elements by atomic number: ')
         if (runmode.ne.'batchrun') then
           write(*,146)
         endif
         read(*,*) iel1,iel2,iel3,iel4
         iel1 = zmap(iel1)
         iel2 = zmap(iel2)
         iel3 = zmap(iel3)
         iel4 = zmap(iel4)
       endif
      endif
c     
c     get final field file prefix
c
 360  format (a8)
 370  format(//,
     & ' Give a prefix for final field file (5chars): ')
      if (runmode.ne.'batchrun') then
        write(*,370)
      endif
      read (*,360) mypfx
c
c     get screen display mode
c
 500  format(//,' Choose runtime screen display:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard display ',/,
     & '    B  : Detailed Slab display.',/,
     & '    C  : Full Display',/,
     & '       : (Full slab display + timescales).',/,
     & ' :: ')
 510  if (runmode.ne.'batchrun') then
        write(*, 500)
      endif
      read(*, 10) ilgg
      ilgg = ilgg(1:1)
c
      if (ilgg.eq.'a') ilgg = 'A'
      if (ilgg.eq.'b') ilgg = 'B'
      if (ilgg.eq.'c') ilgg = 'C'
c
      if (ilgg.eq.'A') vmod = 'MINI'
      if (ilgg.eq.'B') vmod = 'SLAB'
      if (ilgg.eq.'C') vmod = 'FULL'
c
c     get runname
c
 300  format (a64)
 310  format(//,' Give a name/code for this run: ')
      if (runmode.ne.'batchrun') then
        write(*,310)
      endif
      read (*,300) runname

      if (runmode.ne.'batchstart') then !If 'batchstart', goto end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write Headers
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion balance files if requested
c
      if (tsrmod.eq.'y') then
        do i= 1, ieln 
          if (i.eq.1) ie = iel1
          if (i.eq.2) ie = iel2
          if (i.eq.3) ie = iel3
          if (i.eq.4) ie = iel4
          fn = ' '
          pfx = elem(ie)
          sfx = 'sh4'
          if ((mapz(ie).eq.1).or.(mapz(ie).eq.6).or.
     &       (mapz(ie).eq.7).or.(mapz(ie).eq.8).or.
     &       (mapz(ie).eq.16)) then
            call newfile(pfx,1,sfx,3,fn)
          else
            call newfile(pfx,2,sfx,3,fn)
          endif
          filn(i) = fn
        enddo
c     
        do i = 1,ieln
          open(luions(i),file=filn(i),status='NEW')
        enddo
c     
        do i = 1,ieln    
          write(luions(i), *) fl,runname
          write(luions(i), *) dhpr, xhpr, tepr, tmag
          write(luions(i), *) tepo, vshockm
          write(luions(i), *) zgas
c    
 3010     format(A4,29A6)       
          if (i.eq.1) then
           write(luions(i),3010) 'Te ',(rom(j), j = 1,maxion(iel1))
          endif
          if (i.eq.2) then
           write(luions(i),3010) 'Te ',(rom(j), j = 1,maxion(iel2))
          endif
          if (i.eq.3) then
           write(luions(i),3010) 'Te ',(rom(j), j = 1,maxion(iel3))
          endif
          if (i.eq.4) then
           write(luions(i),3010) 'Te ',(rom(j), j = 1,maxion(iel4))
          endif
        enddo
c     
        do i = 1,ieln
          close(luions(i))
        enddo
c
c     end trsmod = 'y'
c
      endif
c
c     Title
c
 320  format(' SHOCK4: Explicit Rankine-Hugoniot shock code: ',/,
     &       ' --------------------------------------------',/,
     &       ' Diffuse Field, Full Continuum Calculations.',/,
     &       ' Calculated by MAPPINGS III v',a6)
      if (runmode.ne.'batchrun') then
        write (*,320) theVersion
      endif
      write (luop,320) theVersion
      write (lusp,320) theVersion
      if (ratmod.eq.'y') write (lurt,320) theVersion
      if (dynmod.eq.'y') write (ludy,320) theVersion
      if (allmod.eq.'y') write (lual,320) theVersion
 330  format(//,' Run  :',a,/,
     &         ' File :',a)
      write (luop,330) runname,fl
      write (lusp,330) runname,fl
      if (ratmod.eq.'y') write (lurt,330) runname,fl
      if (dynmod.eq.'y') write (ludy,330) runname,fl
      if (allmod.eq.'y') write (lual,330) runname,fl
c
c     Write Model Parameters
c
 1000 format(//,' Model Parameters:',/,
     &        ' -----------------',//,
     &        ' Abundances     : ',a80,/,
     &        ' Pre-ionisation : ',a80,/,
     &        ' Photon Source  : ',a80)
 1005 format(/,' Charge Exchange: ',A12,/,
     &        ' Photon Mode    : ',A12,/,
     &        ' Collision calcs: ',A12)
 1010 format(//,' ',t2,'Jden',t8,'Jgeo',t14,'Jend',t20,'Ielen',t26,
     & 'Jpoen',t35,'Fren',t41,'Tend',t47,'DIend',t57,'TAUen',t67,'Tmod'
     & ,t74,'Teini')
 1020 format('  ',3(a4,2x),2(i2,4x),0pf6.4,0pf6.0,2(1pg10.3),3x,a4
     & ,0pf7.1)
 1030 format(//,' Description of Photon Source :')
 1040 format(/,' MOD',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,
     & 'Cut-off',t38,'Zstar',t47,'FQHI',t56,'FQHEI',t66,'FQHEII')
 1050 format(' ',a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
c
c
c
      write(luop, 1000) abnfile,ionsetup,srcfile
      if (ratmod.eq.'y') write(lurt, 1000) abnfile,ionsetup,srcfile
      if (dynmod.eq.'y') write(ludy, 1000) abnfile,ionsetup,srcfile
      if (allmod.eq.'y') write (lual,1000) abnfile,ionsetup,srcfile
c
      cht = 'New Set'
      if (chargemode.eq.2) cht = 'Disabled'
      if (chargemode.eq.1) cht = 'Old Set'
      pht = 'Normal'
      if (photonmode.eq.0) pht = 'Decoupled'
      clt = 'Younger'
      if (collmode.eq.1) clt = 'Old Methods'
c
      write(luop, 1005) cht,pht,clt
      write(luop, 1010) 
      write(luop, 1020) jden, jgeo, jend, ielen, jpoen, fren, 
     &  tend, diend, tauen, tmod, tm00
      write(luop, 1030) 
      write(luop, 1040) 
      write(luop, 1050) iso, teff, alnth, turn, cut, zstar, qhi
     &  , qhei, qheii
c
      if (ratmod.eq.'y') then 
        write(lurt, 1005) cht,pht,clt
        write(lurt, 1010) 
        write(lurt, 1020) jden, jgeo, jend, ielen, jpoen, fren, 
     &   tend, diend, tauen, tmod, tm00
        write(lurt, 1030) 
        write(lurt, 1040) 
        write(lurt, 1050) iso, teff, alnth, turn, cut, zstar, qhi
     &   , qhei, qheii
      endif
c
      if (dynmod.eq.'y') then
        write(ludy, 1005) cht,pht,clt
        write(ludy, 1010) 
        write(ludy, 1020) jden, jgeo, jend, ielen, jpoen, fren, 
     &   tend, diend, tauen, tmod, tm00
        write(ludy, 1030) 
        write(ludy, 1040) 
        write(ludy, 1050) iso, teff, alnth, turn, cut, zstar, qhi
     &  , qhei, qheii
      endif
c
      if (allmod.eq.'y') then
        write(lual, 1005) cht,pht,clt
        write(lual, 1010)
        write(lual, 1020) jden, jgeo, jend, ielen, jpoen, fren,
     &   tend, diend, tauen, tmod, tm00
        write(lual, 1030)
        write(lual, 1040)
        write(lual, 1050) iso, teff, alnth, turn, cut, zstar, qhi
     &   , qhei, qheii
      endif
c
c     Shock parameters
c
 340  format(//,' Jump Conditions:',/,
     &         ' ----------------')
      write (luop,340)
      write (lusp,340)
      if (ratmod.eq.'y') write (lurt,340)
      if (dynmod.eq.'y') write (ludy,340)
      if (allmod.eq.'y') write (lual,340)
c
 140  format(//,' T1',1pg14.7,' V1',1pg14.7,
     &        ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7,/,
     &         ' T2',1pg14.7,' V2',1pg14.7,
     &        ' RH2',1pg14.7,' P2',1pg14.7,' B2',1pg14.7)
c
      dr = 0.d0
      dv = vel1 - vel0
c
      if (runmode.ne.'batchrun') then
      write(*,140) te0,vel0,rho0,pr0,bm0
     &,te1,vel1,rho1,pr1,bm1
      endif
      write(luop,140) te0,vel0,rho0,pr0,bm0
     &,te1,vel1,rho1,pr1,bm1
      write(lusp,140) te0,vel0,rho0,pr0,bm0
     &,te1,vel1,rho1,pr1,bm1
      if (ratmod.eq.'y') write(lurt,140) te0,vel0,rho0,pr0,bm0
     &,te1,vel1,rho1,pr1,bm1
      if (dynmod.eq.'y') write(ludy,140) te0,vel0,rho0,pr0,bm0
     &,te1,vel1,rho1,pr1,bm1
c
c
      close(lusp)
c
c   Output conditions in the precursor step
c
      wdilt0 = te0
      t = te0
      ve = vel0
      dr = 1.0
      dv = ve/100.d0
      fi = 1.d0
      rad = 1.d38
      if (wdil.eq.0.5) rad = 0.d0
      de = feldens(dh,pop)
c
c     calculate the radiation field and atomic rates in the 
c     precursor step 
c
      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, specmode)
      call zetaeff(fi, dh)
      call cool(t,de,dh)
c
 350  format(//,' Precursor Conditions and Ionisation state',/)
      write (luop,350)
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod = 'SCRN'
        call wmodel(luop,t,de,dh,dr,wmod)
      endif
c
      wmod = 'FILE'
      call wmodel(luop,t,de,dh,dr,wmod)
      if (ratmod.eq.'y') call wmodel(lurt,t,de,dh,dr,wmod)
c
      call wionabal(luop,pop)
      if (allmod.eq.'y') call wionabal(lual,pop)
c
      close(luop)
      if (ratmod.eq.'y') close(lurt)
      if (dynmod.eq.'y') close(ludy)
      if (allmod.eq.'y') close(lual)
c
c     Calculate the conditions at the shock front 
c     (directly post-shock front = step #0)
c
      wdilt0 = te1
      t = te1
      cmpf = rho1/rho0
      dh = dh*cmpf
      de = feldens(dh,pop)
      dr = 1.d0
      ve = vel1
      vpo = vel1
      dv = ve/100.d0
      fi = 1.d0
      rad = 1.d38
      if (wdil.eq.0.5) rad = 0.d0
c
c     calculate the radiation field and atomic rates
c
      call localem(t,de,dh)
      call totphot(t, dh, fi, rad, dr, dv, wdil, specmode)
      call zetaeff(fi, dh)
      call cool(t,de,dh)
c
      en = zen*dh
      ue = 3.d0/2.d0*(en+de)*rkb*t
      tnloss = tloss/(en*de)
      cspd = dsqrt(5.d0/3.d0*pr1/rho1)
      wmol = rho1/(en+de)
c
      step = 0
c
 1200 format(//,A4,A1,11(4X,A6,4X,A1))      
 1210 format(I4.3,A1,11(1pg13.6,A1))    
c  
 1220 format(//' #',t9,'Te',t21,'ne',t33,'nH',t45,'ni',t57,
     &     'B',t69,'FHI'/t8,'Time',t22,'dt',t32,'Dist',t45,
     &     'dr',t57,'v',t67,'dlos'/)
 520  format (i4,6(1x,1pg11.4)/t4,6(1x,1pg11.4))
c
      open(luop, file=fl, status='OLD', access='APPEND')
      
      write(luop,1200) 'Step',tab,'Dist.',tab,'Time',
     & tab,'Te',tab,'de',tab,'dh',tab,'en',tab,'Tloss',tab,'dlos',
     & tab,'FHI',tab,'Bbar',tab,'dr',tab
     
      write(luop,1210) step,tab,dist(step+1),tab,timlps(step+1),
     & tab,t,tab,de,tab,dh,tab,en,tab,tloss,tab,dlos,
     & tab,pop(1,1),tab,bm0,tab,dr,tab
      close(luop)
c
      if (runmode.ne.'batchrun') then
         if (vmod.eq.'MINI') then
            write(*,1220)
            write(*,520) step,t,de,dh,en,bm0,pop(1,1),
     &           timlps(step+1),ftime,dist(step+1),
     &           dr,veloc(step+1),dlos
c
c            dtarr = dtime(tarray)
c     
c 900        format('     Step time (map3, real sec): ',(f8.1,1x,i5),/) 
c            write(*,900) tarray(1),int(tarray(2))
c
         endif
      endif
c
c     now compute free enthalpy per particle & pressure
c
      en = zen*dh
      press = (en+de)*rkb*t
      ue = (3.d0/2.d0)*(en+de)*rkb*t
c
c     Calculate timescales in the shock front zone
c
c     effective cooling timescale, based on net loss
c
      cltime = dabs(ue/tloss)
c
      rhotot = frho(de,dh)
      wmol = rhotot/(en+de)
      mu = en*wmol/((en+de)*(mp+me))
c
      Hfree = Bmag*Bmag/(4.d0*pi*en) + 0.5d0*(wmol*ve*ve)
     &     + 5.d0/2.d0*(rkb*t)/mu
c
      htime = dabs(Hfree/tloss)
c
      abstime = 1.d0/epsilon
c
      ctime = fcolltim(de)
      rtime = frectim2(de)
      ptime = fphotim()
c
      eqtime = (1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
      eqtime = dabs(1.d0/eqtime)
c
c     experiment with timescales
c
      stime = (1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      stime = dabs(1.d0/stime)
c
      ftime = 1.d0/((1.d0/(stime*atimefrac))
     &             +(1.d0/(cltime*cltimefrac)))
c     &             +(1.d0/(abstime*abtimefrac)))
c
      if ((tmod.eq.'m').and.(ftime.gt.utime)) ftime = utime
      if ((tmod.eq.'n').and.(ftime.lt.utime)) ftime = utime
      if (tmod.eq.'s') ftime = utime
c
 130  format ( //' Tot. Internal :',1pg14.7,/
     &     ' Tot.loss rate :',1pg14.7,/' Eff.loss rate :',1pg14.7,/
     &     ' Enthaply time :',1pg14.7,/' Eff.Cool time :',1pg14.7,/
     &     ' Collis.  time :',1pg14.7,/' Recomb.  time :',1pg14.7,/ 
     &     ' Photo.   time :',1pg14.7,/' Atomic   time :',1pg14.7,/
     &     ' Equilib  time :',1pg14.7,/' Absorb.  time :',1pg14.7,/
     &     ' Choice   time :',1pg14.7/)
c
      if (runmode.ne.'batchrun') then
         if (vmod.eq.'FULL') then 
            write(*,130) ue,Tloss,Tloss,htime,cltime,ctime
     &           ,rtime,ptime,stime,eqtime,abstime, ftime
         endif
      endif
c
      if (ratmod.eq.'y') then 
        open(lurt, file=fr, status='OLD', access='APPEND')
        write(lurt,130) ue,Tloss,Tloss,htime,cltime,ctime
     &   ,rtime,ptime,stime,eqtime,abstime,ftime      
        close(lurt)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Setup Initial first step in post shock-front gas(#1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rmod = 'NEAR'
      call rankhug(t, de, dh, ve, Bmag, rmod,tloss,ftime)
c
      dr = ftime*(vel0+vel1)/2.d0
      dv = vel1 - vel0
c
c      write(*,200) tloss,dr,dv
c
      step = 0
c
c     init arrays first step
c
      xh(1) = 0.d0
      xh(2) = pop(2,1)
      veloc(1) = 0.d0
      veloc(2) = vel0
      te(1) = 0.d0
      te(2) = te0
      dhy(1) = 0.d0
      dhy(2) = dh
      dist(1) = 0.0d0
      dist(2) = 0.0d0
      deel(1) = 0.d0
      deel(2) = de
      timlps(1) = 0.0d0
c
c     prepare for time step in ions
c
c
c     record initial conditions, t0 and l0 and pop
c
      call copinto(pop,pop0)
c
      tt0 = te0
      l0 = tloss
      l12 = 0.d0
      l11 = 0.d0
      r0 = rho0
      p0 = pr0
      v0 = vel0
      de0 = de
      dh0 = dh
      count = 0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SHOCK Iteration reentry point:
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 150  continue
c
      count = count +1
c
      call copinto(pop0,pop)
c
      cmpf = rho1/rho0
      dh = (dh0+dh1)/2.d0
      de = feldens(dh,pop)
      t = (te1+te0)/2.d0
      xhii = pop(2,1)
c
      call timion(t,de,dh,xhii,ftime)
c
      call copinto(pop,pop1)
c
c     after timion on the average t,de,dh, now
c     calculate new lossrate at te1,de1,dh1
c
c     calculate the radiation field and atomic rates
c
      call localem(t,de,dh)
      rad = dist(step+1)
      call totphot(t, dh, fi, rad, dr, dv, wdil, specmode)
      call zetaeff(fi, dh)
c
      t = te1
      dh = dh1
      de = feldens(dh,pop)
c
      call cool(t,de,dh)
c
c     record first guess
c
      l12 = l11
      l11 = tloss
      l1 = dabs(1.d0-dabs(l11)/(epsilon+dabs(l12)))
c
c
c    test for pseudo equilibrium
c
c      if ((eqtime.le.stime).and.
c     &   (eqtime.le.cltime).and.
c     &   (eqtime.le.abstime*abtimefrac).and.
      if ((rtime.le.ctime).and.
     &   (ptime.le.ctime).and.
     &   (dlos.lt.0.01d0)) then
         if (count.gt.1) then
            te1 = te0
            dh1 = dh0
            vel1 = vel0
            bm1 = bm0
         endif
        goto 250
      endif
c
      tloss = (l11+l0)/2.d0
c     
c     allow convergence to relax as dlos->0
c
      tlim = (1.d-4)/(dabs(dlos)+epsilon)
c
c     if within limits go and write step
c
      if (l1.le.tlim) then
         write(*,*) '    Converged',l1
         call copinto(pop,pop1)
         call averinto(0.5d0,pop0,pop1,pop)
         t = (te0+te1)/2.d0
         dh = (dh0+dh1)/2.d0
         de = feldens(dh,pop)
         call cool(t,de,dh)
         goto 160 
      endif
c
c     try 5 times then give up
c
      if (count.gt.4) then
c         count = 0
c         ftime = ftime/2.d0
c         write(*,*) '    Failed to converge',l1
         call copinto(pop,pop1)
         call averinto(0.5d0,pop0,pop1,pop)
         t = (te0+te1)/2.d0
         dh = (dh0+dh1)/2.d0
         de = feldens(dh,pop)
         call cool(t,de,dh)
         goto 160 
      endif
c
      en = zen*dh + de
c      
      press = en*rkb*t
      ue = (3.d0/2.d0)*press
c      
      cltime = dabs(ue/tloss)
c
      if (ftime .gt.cltime) then
         ftime = cltime
      endif
c
      t = tt0
      dh = dh0
      de = feldens(dh,pop0)
      ve = v0
      bmag = bm0
c
      call copinto(pop1,pop)
c
      rmod = 'NEAR'
      call rankhug(t, de, dh, ve, Bmag, rmod, tloss,ftime)
c
      dr = ftime*(vel0+vel1)/2.d0
      dv = vel1 - vel0
c
 200  format(t4,3(1pg15.7,1x))
      write(*,200) tloss,dr,dv
c
c     otherwise continue with the normal NEQ
c     mode
c
      goto 150
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pseudo equilibrium calculations for case where
c     heating and cooling very nearly match and extrapolation
c     based on net cooling rate is unstable....
c     = photoionized post-shock zone
c
c     pseudo equilibrium entry point
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 250  continue
c
      write(*,*) 'Equilibrium'
c
      pseudo = .TRUE.
c
      count = count +1
c
      call copinto(pop0,pop)
c
      dh = dh1
      de = feldens(dh,pop)
c
c     initial internal energy /unit mass
c
      en = zen*dh
      ue0 = (3.d0/2.d0)*(en+de)*rkb*te1
c
      tt0 = te1
      t = te1
c
      imod = 'TIM'
      call teequi(tt0,t,de,dh,ftime,imod)
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod = 'SCRN'
        call wmodel(luop,t,de,dh,dr,wmod)
      endif
c
c     compute new internal energy
c     de is updated
c
      teq = t
      dheq = dh
c
      en = zen*dheq
      ue1 = (3.d0/2.d0)*(en+de)*rkb*t
c
c     Net Loss
c
      tloss = (ue0-ue1)/ftime
c
      if (ratmod.eq.'y') then
        open(lurt, file=fr, status='OLD', access='APPEND')
        write(lurt,*) 'pse',ue0,ue1,tloss,tt0,t
        close(lurt)
      endif
c
      call copinto(pop,pop1)
      call copinto(pop0,pop)
      t = tt0
      dh = dh1
      de = feldens(dh,pop)
      ve = vel1
      bmag = bm1
      call copinto(pop1,pop)
c
      rmod = 'NEAR'
      call rankhug(t, de, dh, ve, Bmag, rmod,tloss,ftime)      
c
      if (dabs(1.d0-(teq/te1)).gt.1.d-3) then         
         tl0 = tloss*0.5
         bmag = bm1
         call rankhug(t, de, dh, ve, Bmag, rmod,tl0,ftime)      
         tt1 = te1
         tl1 = tloss*0.5
         bmag = bm1
         call rankhug(t, de, dh, ve, Bmag, rmod,tl1,ftime)      
         tl1 = tl0+(((teq-tt1)/(te1-tt1))*(tl1-tl0))
         bmag = bm1
         call  rankhug(t, de, dh, ve, Bmag, rmod,tl1,ftime) 
         tloss = tl1
      endif
c
      dr = ftime*(vel0+vel1)/2.d0
      dv = vel1 - vel0
c
      call averinto(0.5d0,pop0,pop1,pop)
      t = (te0+te1)/2.d0
      dh = (dh0+dh1)/2.d0
      de = feldens(dh,pop)
c
      call cool(t,de,dh)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Entry point for calculating recording step and end values
c   Both Equilibrium and Non-Equilibrium end up here
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 160  continue
c
      step = step+1
c
c     record step
c
      xh(step+1) = pop(2,1)
      veloc(step+1) = (vel1+vel0)*0.5d0
      te(step+1) = t
      dhy(step+1) = dh
      deel(step+1) = de
      dist(step+1) = dist(step)+dr
      timlps(step+1) = timlps(step)+ftime
c
c
 170  format(i4,1x,i4)
 180  format(1pg14.7,'K',1x,1pg14.7,'cm/s',1x,1pg14.7,'g/cm3',1x,1pg14.7
     &,'dyne/cm2',1x,1pg14.7,'Gauss',1x,1pg14.7,'ergs/cm3/s')
 190  format(2(1pg14.7,' cm',1x),2(1pg14.7,' s',1x))
 195  format(2(1pg14.7,' ergs/cm^3',1x),1pg14.7,' /cm^3',1x,
     &1pg14.7,' cm/s',1x,1pg14.7,' g ')
c
      en = zen*dh
c
      ue = 3.d0/2.d0*(en+de)*rkb*t
      tnloss = tloss/(en*de)
c
      cspd = dsqrt(5.d0/3.d0*(pr1+pr0)/(rho1+rho0))
      wmol = (rho1+rho0)/(2.d0*(en+de))
c
      mb = (bm0+bm1)*0.5
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
      wmod = 'SCRN'
      call wmodel(luop,t,de,dh,dr,wmod)
      endif
c
      if (runmode.ne.'batchrun') then 
c
 525     format('     Case A-B (HI, HeII): ',2(1pg11.3))
         write(*,525) caseab(1), caseab(2)
c     
         if (vmod.eq.'MINI') then
c     
            write(*,520) step,t,de,dh,en,mb,pop(1,1),
     &           timlps(step+1),ftime,dist(step+1),
     &           dr,veloc(step+1),dlos
         endif
c     
c         dtarr = dtime(tarray)
c         nt1   = time()
c         
c     
c         write(*,900) tarray(1),nt1-nt0
c         
c         nt0 = nt1
c
      endif
c
      if (ratmod.eq.'y') then
      open(lurt, file=fr, status='OLD', access='APPEND')
      wmod = 'FILE'
      call wmodel(lurt,t,de,dh,dr,wmod)
      close(lurt)
      endif
c
      if (dynmod.eq.'y') then
      open(ludy, file=fd, status='OLD', access='APPEND')
      write(ludy,170) step,count
      write(ludy,190) dist(step+1),dr,timlps(step+1),ftime
      write(ludy,180) te0,vel0,rho0,pr0,bm0,l0
      write(ludy,180) te1,vel1,rho1,pr1,bm1,l11
      write(ludy,195) ue,Hfree,en,cspd,wmol
      close(ludy)
      endif
c
      open(luop, file=fl, status='OLD', access='APPEND')
      write(luop,1210) step,tab,dist(step+1),tab,timlps(step+1),
     &tab,t,tab,de,tab,dh,tab,en,tab,tloss,tab,dlos,
     &tab,pop(1,1),tab,(0.5*(bm0+bm1)),tab,dr,tab
      close(luop)
c
      if (allmod.eq.'y') then
         open(lual,file=fa,status='OLD',access='APPEND')
         write(lual, 1250)
         write(lual, 1251) step, t, de, dh, vel1, rho1, pr1,
     &                 dist(step+1),timlps(step+1)
 1250 format(//,' Step   Te Ave.(K)   ',
     &                '  ne(cm^-3)   ',
     &                '  nH(cm^-3)   ',
     &                '  V1(cm/s)    ',
     &                ' Rho1(g/cm^3) ', 
     &                ' Pr1(erg/cm^3)',
     &                '   Dist.(cm)  ',
     &                ' Elps. Time(s)')
 1251 format(1x,i4,8(1pg14.7)/)
c
         call wionabal(lual,pop)
         call wionabal(lual,popint)
         close(lual)
      endif
c
      if (tsrmod.eq.'y') then
      do i = 1 ,ieln 
      open(luions(i),file=filn(i),status='OLD',access='APPEND')
c
 3000 format(1pg12.5,A1,29(1pg11.4,A1))
      if (i.eq.1) then
         write(luions(i),3000) te1,tab,(pop(j,iel1),tab, 
     &        j = 1,maxion(iel1))
      endif
      if (i.eq.2) then
         write(luions(i),3000) te1,tab,(pop(j,iel2),tab, 
     &        j = 1,maxion(iel2))
      endif
      if (i.eq.3) then
         write(luions(i),3000) te1,tab,(pop(j,iel3),tab, 
     &        j = 1,maxion(iel3))
      endif
      if (i.eq.4) then
         write(luions(i),3000) te1,tab,(pop(j,iel4),tab, 
     &        j = 1,maxion(iel4))
      endif
c
      close(luions(i))
c
      enddo
c
c     end .tsr files
c
      endif
c
c     poll for photons file
c
      pollfile = 'photons'
      inquire(file=pollfile, exist=iexi)
c
c
      if ((lmod.eq.'y').or.(iexi)) then
c
c     Local photon field
c         
        specmode = 'LOCL'
        call totphot(t, dh, 1.0d0, 1.d38, 1.d0, 0.d0, 0.5d0, specmode)
        call zetaeff(fi, dh)
        caller = 'S4'
        pfx = 'plocl'
        np = 5
        wmod = 'REAL'
        dva = vel1-vel0
        call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        vel1,dva,vshoc,
     &        dist(step+1),dr,
     &        timlps(step+1),ftime,
     &        1.d0,tphot)
c     
        pfx = 'slec'
        sfx = 'sh4'
        call newfile(pfx,4,sfx,3,fpf)
        open(lupf, file=fpf, status='NEW') 
        spmod = 'ABS'
        linemod = 'LAMB'
        call speclocal(lupf,tloss, eloss, egain, dlos, t, de, dh
     &                , pop(1,1), dist(step+1),dr,linemod,spmod)
        close(lupf) 
c     
        open(luop, file=fl, status='OLD', access='APPEND')
c
c     total energy in Inu
c
c
        band1 = 0.d0
        band2 = 0.d0
        band3 = 0.d0
        band4 = 0.d0
        band5 = 0.d0
        band6 = 0.d0
        bandall = 0.d0
c       
        write(*,*) epot(2,2)/ev
        do i = 1,infph-1
          wid = ephot(i+1)-ephot(i)
          bandall = bandall + tphot(i)*wid*evplk
          if ((ephot(i).ge.1.0d2).and.(ephot(i).lt.2.0d2)) then
             band1 = band1 + tphot(i)*wid*evplk
          endif
          if ((ephot(i).ge.2.0d2).and.(ephot(i).lt.3.0d2)) then
             band2 = band2 + tphot(i)*wid*evplk
          endif
          if ((ephot(i).ge.3.0d2).and.(ephot(i).lt.4.0d2)) then
           band3 = band3 + tphot(i)*wid*evplk
          endif
          if ((ephot(i).ge.4.0d2).and.(ephot(i).lt.5.0d2)) then
             band4 = band4 + tphot(i)*wid*evplk
          endif
          if ((ephot(i).ge.5.0d2).and.(ephot(i).lt.2.0d3)) then
             band5 = band5 + tphot(i)*wid*evplk
          endif
          if ((ephot(i).ge.2.0d3).and.(ephot(i).lt.7.0d3)) then
             band6 = band6 + tphot(i)*wid*evplk
          endif
        enddo
c
c     times 4 pi to cover whole sky
c
        band1   =  band1  *4.d0*pi/(wdil*dr)
        band2   =  band2  *4.d0*pi/(wdil*dr)
        band3   =  band3  *4.d0*pi/(wdil*dr)
        band4   =  band4  *4.d0*pi/(wdil*dr)
        band5   =  band5  *4.d0*pi/(wdil*dr)
        band6   =  band6  *4.d0*pi/(wdil*dr)
        bandall =  bandall*4.d0*pi/(wdil*dr)
c     
        tab = char(9)
 1230   format(14(A12,A1)) 
        write(*,1230) "Te",tab,"de",tab,"dh",tab,"en",tab,"FHI",tab
     &   ,"ff",tab,"tloss",
     &   tab,"B1-200eV",tab,"B2-300eV",tab,"B3-400eV",
     &   tab,"B4-500eV",tab,"B0.5-2.0keV",tab,"B2.0-7.0keV",tab,"Ball"
     
 1231   format(14(1pg12.5,A1))
        
c        write(luop,1231) t,tab,de,tab,
c     &   dh,tab,en,tab,pop(1,1),tab,fflos,tab,tloss,tab,band1,
c     &   tab,band2,tab,band3,tab,band4,tab,band5,tab,band6,tab,bandall
        write(*,1231) t,tab,de,tab,
     &   dh,tab,en,tab,pop(1,1),tab,fflos,tab,tloss,tab,band1,
     &   tab,band2,tab,band3,tab,band4,tab,band5,tab,band6,tab,bandall
c       
        close(luop)
c
c     reset mode and tphot
c
        specmode = 'DW'
        rad = dist(step+1)
        call totphot(t, dh, fi, rad, dr, dv, wdil, specmode)      
c
      endif
c
c     get mean ionisation state for step
c
      rdis = (dist(step)+dist(step+1))/2.0d0
      cmpf = rho1/rho0
      fi = 1.0d0
c
c     accumulate spectrum
c
      tdw = te0
      tup = te1
      drdw = dr
      dvdw = vel0-vel1
      drup = dist(step+1)
      dvup = dsqrt(vpo*vel0)
      frdw = 0.5d0
c
      jmod = 'LODW'
      imod = 'ALL'
c
      call localem(t,de,dh)
      call zetaeff(fi, dh)
      call newdif(tdw, tup, dh, fi, rad, drdw, dvdw, drup, dvup, frdw
     &           ,jmod)
c
      call sumdata(t, de, dh, fi, dr, dr, rdis, imod)
c
c     record line ratios
c     
      if (ox3.ne.0) then
         hoiii(step+1) = (fluxf6(7,ox3)+fluxf6(10,ox3))
     &                   /(fluxh(2)+epsilon)
      endif
      if (ox2.ne.0) then
         hoii(step+1) = (fluxf(1,ox2)+fluxf(2,ox2))
     &                   /(fluxh(2)+epsilon)
      endif
      if (ni2.ne.0) then
         hnii(step+1) = (fluxf6(7,ni2)+fluxf6(10,ni2))
     &                   /(fluxh(2)+epsilon)
      endif
      if (su2.ne.0) then
         hsii(step+1) = (fluxf(1,su2)+fluxf(1,su2))/(fluxh(2)+epsilon)
      endif
      if (ox1.ne.0) then
      if (ispo .eq. ' OI') hsii(step+1) = fluxf(3,ox1)
     &                   /(fluxh(1)+epsilon)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set up initial conditions for the next zone 
c     (= end conditions in previous zone)
c     Calculate new cooling times and next step
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copinto(pop1,pop0)
c
      t = te1
      tt0 = te1
      l0 = l11
      l12 = 0.d0
      l11 = 0.d0
      r0 = rho1
      p0 = pr1
      v0 = vel1
      dh0 = dh1
      de0 = de1
      bm0 = bm1
      bmag = bm0
      count = 0
c
      call cool(t, de0, dh0)
c
c      write(*,*) t,de0,dh0,tloss
c
c     now compute free enthalpy per particle
c
      en = zen*dh0
c
      press = (en+de0)*rkb*t
      ue = (3.d0/2.d0)*(en+de0)*rkb*t
c
      cltime = dabs(ue/tloss)
c
      wmol = r0/(en+de0)
      mu = en*wmol/((en+de0)*(mp+me))
c
      Hfree = Bmag*Bmag/(4.d0*pi*en) + 0.5d0*(wmol*ve*ve)
     &     + 5.d0/2.d0*(rkb*t)/mu
c
c     calculate timescales
c
      htime = dabs(Hfree/tloss)
c
      abstime = 1.d0/epsilon
c      if (photonmode.eq.1) then
c      call absdis(t, dh, 1.0d0, abtimefrac, dr, 0.d0, pop)
c      abstime = (dr/vel1)/abtimefrac
c      endif
c
      ctime = fcolltim(de)
      rtime = frectim2(de)
      ptime = fphotim()
c
      eqtime = (1.d0/ctime)+(1.d0/ptime)+(1.d0/rtime)
      eqtime = dabs(1.d0/eqtime)
c
c     experiment with timescales
c
      stime = (1.d0/ctime)+(1.d0/ptime)-(1.d0/rtime)
      stime = dabs(1.d0/stime)
c
      ftime = 1.d0/((1.d0/(stime*atimefrac))
     &             +(1.d0/(cltime*cltimefrac)))
c     &             +(1.d0/(abstime*abtimefrac)))
c
      if ((tmod.eq.'m').and.(ftime.gt.utime)) ftime = utime
      if ((tmod.eq.'n').and.(ftime.lt.utime)) ftime = utime
      if (tmod.eq.'s') ftime = utime
c
      if (runmode.ne.'batchrun') then
         if (vmod.eq.'FULL') then
            write(*,130) ue,Tloss,Tloss,htime,cltime,ctime
     &           ,rtime,ptime,stime,eqtime,abstime,ftime
         endif
      endif
c
      if (ratmod.eq.'y') then
        open(lurt, file=fr, status='OLD', access='APPEND')
        write(lurt,130) ue,Tloss,Tloss,htime,cltime,ctime
     &   ,rtime,ptime,stime,eqtime,abstime,ftime      
        close(lurt)
      endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Program endings
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ionisation, finish when mean ionisation < 1%
c
      if (jend.eq.'a') then
         totn = 0.d0
         totn = zen
         ionbar = 0.d0
         do i = 1, atypes     
            do j = 2 , maxion(i)
               ionbar = ionbar+(pop(j,i)*zion(i))
            enddo
         enddo
         ionbar = ionbar/totn
         if (ionbar.lt.0.01d0) goto 999
      endif
c
c     Specific species ionisation limit
c     
      if (jend.eq.'b') then
         if (pop(jpoen,ielen).lt.fren) goto 999
      endif
c
c     Temperature limit
c
      if ((jend.eq.'c').and.(t.lt.tend)) goto 999
c
c     Distance Limit
c
      if ((jend.eq.'d').and.(dist(step+1).ge.diend)) goto 999
c
c     Time Limit
c
      if ((jend.eq.'e').and.(timlps(step+1).ge.timend)) goto 999
c
c     thermal balance dlos<1e-2
c
      if ((jend.eq.'f').and.(dlos.lt.1.d-2)) goto 999
c
c     cooling function test, finish when tloss goes -ve
c
      if ((jend.eq.'g').and.(tloss.lt.0.d0)) goto 999
c
c     poll for terminate file
c
      pollfile = 'terminate'
      inquire(file=pollfile, exist=iexi)
      if (iexi) goto 999
c
c     otherwise go to normal iteration loop, if interlocks
c     permit
c
      if ((t.gt.100.d0).and.(step.lt.(nstmax-1))) then
         goto 150
      endif
c
 999  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     End model
c     write out spectrum etc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jspec.eq.'YES') then
c     
c     Final Downstream Field
c
        caller = 'S4'
        pfx = 'spdw'
        np = 4
        dva = 0
        wmod = 'REAL'
        call wpsou(caller,pfx,np,wmod,
     &        t,de,dh,
     &        vel1,dva,vshoc,
     &        dist(step+1),dr,
     &        timlps(step+1),ftime,
     &        (2.d0/4.d0),tphot)
      endif
c
c     Upstream photon field
c
      specmode = 'UP'
      rad = dist(step+1)
      call totphot(t, dh, fi, rad, dr, dv, wdil, specmode)      
      caller = 'S4'
      pfx = mypfx
      np = 5
      dva = vel1-vel0
c
      wmod = 'REAL'
      call wpsou(caller,pfx,np,wmod,
     &      te1,de1,dh1,
     &      vel1,dva,vshoc,
     &      dist(step+1),dr,
     &      timlps(step+1),ftime,
     &      (2.d0/4.d0),tphot)
c
      wmod = 'NFNU'
      call wpsou(caller,pfx,np,wmod,
     &      te1,de1,dh1,
     &      vel1,dva,vshoc,
     &      dist(step+1),dr,
     &      timlps(step+1),ftime,
     &      (2.d0/4.d0),tphot)
c     
c     dynamics
c
      open(luop, file=fl, status='OLD', access='APPEND')
c
 1060 format(//,' Model ended.',a1,'Distance:',1pg14.7,a1,
     & 'Time:',1pg14.7,a1,'Temp:',1pg12.5//)
      write(luop,1060) tab,dist(step+1),tab,timlps(step+1),tab,t
c 
      call avrdata
      call wrsppop(luop)
      close(luop)
c  
      open(lusp, file=fsp, status='OLD', access='APPEND')
      spmod = 'REL'
      linemod = 'LAMB'
      call spec2(lusp,linemod,spmod)
      close(lusp) 
c
      endif !     end batchstart loop

c
c     restore photonmode!
c
      photonmode = 1
c
      return
c
      end
      
