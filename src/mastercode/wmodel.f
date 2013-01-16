cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     wmodel:write out current plasma model conditions.
c     
c     wmod = 'SCRN'   for vertical format, luop is not used.
c     wmod = 'FILE'   for side by side format, luop is 
c                     assumed to be open.
c     wmod = 'PROP'   properties only, inline,luop
c                     assumed to be open.
c     wmod = 'LOSS'   losses and gains only in line, luop
c                     assumed to be open.
c
c
c     runmode = 'batchrun' overrides scrn
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wmodel(luop,t,de,dh,dstep,wmod)
c
      include 'cblocks.inc'
c
      double precision t,de,dh,dstep,rl,fsl,sint,aspd
      double precision q2los,cspd,ram,mag
      double precision trec,tcoll,press,en,ue,cts,rhotot,wmol
      integer*4 luop
      character wmod*4
c
c           Functions
c
      double precision fcietim,fpressu,frectim,frho
c
      press = fpressu(t,dh,pop)
      rhotot = frho(de,dh)
c
      trec = frectim(t,de,dh)
      tcoll = fcietim(t,de,dh)
c
      ram = rhotot*vel0*vel0
c
      q2los = h2ql+heii2ql
c
      en = zen*dh
c
      ue = 3/2*(en+de)*rkb*t
      cts = ue/tloss
c
      cspd = dsqrt(5/3*press/rhotot)
c
      mag = bm0
      aspd = dsqrt(mag*mag/(4.d0*pi*rhotot))
c
      wmol = rhotot/(en+de)
c
      rl = xrloss+rloss
      fsl = fslos+xiloss
c
      if(wmod.eq.'LOSS') then
c
 200  format(25(1pg11.4,1x))
c
      write(luop,200) dstep,t,de,dh,en,dlos,eloss,egain,hloss,
     &chgain,xrloss,rloss,fsl,xheloss,
     &floss,feloss,q2los,cmplos,fflos,
     &colos,pgain,rngain,cosgain,
     &gheat-gcool,paheat
c
      endif
      if(wmod.eq.'PROP') then
c
 100  format(15(1pg11.4,1x))
c
      write(luop,100)
     &dstep,t,de,dh,en,pop(1,1), Rhotot,press,vel0,ram,cspd,dstep,
     &aspd,bm0,(avgpot(1)+avgpot(2))/2.d0
c
      endif
c
      if (wmod.eq.'FILE') then 
      write(luop,303) tloss,eloss,tloss,egain,dlos,hloss,chgain,
     &rloss,xrloss
c
 303  format(/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::::: ',
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TOTAL LOSS  :',1pg11.3,' : EFF. LOSS   :',1pg11.3,' : ',
     &': TOTAL LOSS:', 1pg12.5,' :::::::::::::::::::::::::::::'/
     &': EFF. GAIN   :',1pg11.3,' : FRAC. RESID.:',1pg11.3,' : ',
     &': COLEXC H. :', 1pg12.5,' : CHARGE EX.  :', 1pg12.5,' :'/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::::: ',
     &': RESON.    :', 1pg12.5,' : X. RESON.   :', 1pg12.5,' :')
c
      write(luop,304)
     &     t,rhotot,fsl,f3loss,
     &     en,de,floss,feloss,
     &     pop(1,1), press,q2los,cmplos,
     &     vel0,ram,fflos,colos,
     &     cspd,dstep,pgain,rngain,aspd,mag,
     &     cosgain,paheat,avgpot(1),avgpot(2),(gheat-gcool)     
c     
 304  format(
     &': TEMP.       :',1pg11.4,' : DENSITY     :',1pg11.3,' : ',
     &': INTER.    :', 1pg12.5,' : 3 LEV. FINE :', 1pg12.5,' :'/
     &': N IONS      :',1pg11.3,' : N ELECTRONS :',1pg11.3,' : ',
     &': FORBID    :', 1pg12.5,' : FE II       :', 1pg12.5,' :'/
     &': FRAC. N. H  :',1pg11.3,' : PRESSURE    :',1pg11.3,' : ',
     &': 2PHOTON   :', 1pg12.5,' : COMPTON     :', 1pg12.5,' :'/
     &': FLOW VELOC. :',1pg11.3,' : RAM PRESS.  :',1pg11.3,' : ',
     &': FREFRE    :', 1pg12.5,' : COLION      :', 1pg12.5,' :'/
     &': SOUND SPEED :',1pg11.3,' : SLAB DEPTH  :',1pg11.3,' : ',
     &': PGAIN     :', 1pg12.5,' : RNGAIN      :', 1pg12.5,' :'/
     &': ALFEN SPEED :',1pg11.3,' : MAG. FIELD  :',1pg11.3,' : ',
     &': COSMIC    :', 1pg12.5,' : PHOTO.  PAH :', 1pg12.5,' :'/
     &': GRA GRN POT :',1pg11.3,' : SIL GRN POT :',1pg11.3,' : ',
     &': PHOT. GRN :', 1pg12.5,' :'/
     &'::::::::::::::::::::::::::::::::::::::::::::::::::::::: ',
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
c
      endif
c
      if (wmod.eq.'SCRN') then
c
      if (runmode.ne.'batchrun') then
      write(*,301) tloss,eloss,egain,dlos,
     &t,rhotot,en,de,pop(1,1), press,vel0,ram,cspd,dstep
      write(*,302) aspd,mag,avgpot(1),avgpot(2),
     &sint,wdil,teff,alnth,turn,cut
      endif
 301  format(/
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TOTAL LOSS  :',1pg11.3,' : EFF. LOSS   :',1pg11.3,' :'/
     &': EFF. GAIN   :',1pg11.3,' : FRAC. RESID.:',1pg11.3,' :'/
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TEMP.       :',1pg11.4,' : DENSITY     :',1pg11.3,' :'/
     &': N IONS      :',1pg11.3,' : N ELECTRONS :',1pg11.3,' :'/
     &': FRAC. N. H  :',1pg11.3,' : PRESSURE    :',1pg11.3,' :'/
     &': FLOW VELOC. :',1pg11.3,' : RAM PRESS.  :',1pg11.3,' :'/
     &': SOUND SPEED :',1pg11.3,' : SLAB DEPTH  :',1pg11.3,' :')
 302  format(
     &': ALFEN SPEED :',1pg11.3,' : MAG. FIELD  :',1pg11.3,' :'/
     &': GRA GRN POT.:',1pg11.3,' : SIL GRN POT.:',1pg11.3,' :'//
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': SOURCE INT. :',1pg11.3,' : DILUT. F.   :',1pg11.3,' :'/
     &': STAR TEMP.  :',1pg11.3,' : ALPHA       :',1pg11.3,' :'/
     &': TURN-ON     :',1pg11.3,' : CUT-OFF     :',1pg11.3,' :')
c
      if (runmode.ne.'batchrun') then
      write(*,300) tloss,hloss,chgain,xrloss,rloss,fsl,f3loss,
     &floss,feloss,q2los,cmplos,fflos,colos,pgain,rngain,cosgain,
     &paheat,(gheat-gcool)
      endif
 300  format(
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &': TOTAL LOSS:', 1pg12.5,' :::::::::::::::::::::::::::::'/
     &': COLEXC H. :', 1pg12.5,' : CHARGE EX.  :', 1pg12.5,' :'/
     &': RESON.    :', 1pg12.5,' : X. RESON.   :', 1pg12.5,' :'/
     &': INTER SYS.:', 1pg12.5,' : 3 LEV. FINE :', 1pg12.5,' :'/
     &': FORBID    :', 1pg12.5,' : FE II       :', 1pg12.5,' :'/
     &': 2PHOTON   :', 1pg12.5,' : COMPTON     :', 1pg12.5,' :'/
     &': FREFRE    :', 1pg12.5,' : COLION      :', 1pg12.5,' :'/
     &': PGAIN     :', 1pg12.5,' : RNGAIN      :', 1pg12.5,' :'/
     &': COSMIC    :', 1pg12.5,' : PHOTO. PAH  :', 1pg12.5,' :'/
     &': PHOT. GRN :', 1pg12.5,/
     &':::::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
c
      endif
c
      return
c
      end
