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
c     This routine returns the non-relativistic net compton
c     heating.
c
c     ref: Krolik, McKee and Tarter Ap. J. 249:422-442
c
c     Thanks to Wang Chi Lin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compton(t, de, dh)
c
      include 'cblocks.inc'
c
c
      double precision t, de, dh
      double precision sigmat,f1,eav
      double precision sumfe, sumf
      double precision energ,wid
      integer*4 i
c
c
      cmplos   = 0.d0
      cmpcool = 0.d0
      cmpheat = 0.d0
c
c     relativistic (2hv/mc^2) term neglected
c
      sigmat = 6.6524d-25
c
c     get totalflux sumf, and mean photon energy eav
c
      sumf = 0.0d0
      sumfe = 0.0d0
      do i = 1,infph-1
         energ = ev*0.5d0*(ephot(i+1)+ephot(i))
         wid   = (ephot(i+1)-ephot(i))
         f1      = tphot(i)*wid*evplk
         sumfe = sumfe + f1*energ
         sumf  = sumf  + f1
      enddo
c
      if (sumf.gt.0.d0) then
c
      eav = sumfe/sumf
c
c     get the loss rate...
c
      cmplos = ((sigmat*sumf)/(me*cls*cls))*de*((4*rkb*t)-eav)
      cmpcool= ((sigmat*sumf)/(me*cls*cls))*de*((4*rkb*t))
      cmpheat= ((sigmat*sumf)/(me*cls*cls))*de*(-eav)
      endif
c     
      if (expertmode.gt.0) then
         write(*,99) 'Te,Ctot,ht,cl:',t,cmplos,cmpcool,cmpheat
 99      format (a15,1p4e15.4)
      endif
c 
      return 
      end
