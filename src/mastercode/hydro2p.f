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
c*******COMPUTES HI AND HEII COOLING BY TWO QUANTA EMISSION FOR 
c       HI ,HEI AND HEII
c
c	COMPUTE ALSO THE SPECTRUM OF HI ADDING THE RECOMB. COMPONENT
c	USE INTERMEDIATE CASE A,B FOR THE COLLISIONAL EXCITATION
c	COMPONENT OF THE BALMER LINES.  THESE ARE THE OLD METHODS
C   (CF HYDRO.F THE NEW METHODS) BUT THE TWO PHOTON CALCS ARE
C   SO ENTWINED IT IS EASIER TO RETAIN THE OLD COLL AND RECOMBINATION
C   CALCS AND ONLY USE THEM TO DERIVE THE TWO PHOTON VALUES.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hydro2p(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision t, de, dh,qpr,qel,reff,x1,fr,u
      double precision ab,abde,ey,frab,hbcol
      double precision r2q,rate,rt,w
c
      integer*4 k,l,ll,ml,n,nz
c
      qpr(u) = 4.74d-4*(u ** (-0.151d0))
      qel(u) = 0.57d-4*(u ** (-0.373d0))
      reff(u) = 0.838d-13*(u ** (-0.721d0))
      x1(u) = 0.7179d0*(u ** (-0.0584d0))
      fr(u) = 0.588d0*(u ** (-0.234d0))
c
      ml = 8
      u = t/1.d4
      if (u.lt.0.1d0) u = 0.1d0
      if (u.gt.10.d0) u = 10.d0
      w = t/4.d4
      if (w.lt.0.1d0) w = 0.1d0
      if (w.gt.10.d0) w = 10.0d0
c
c    ***COMPUTES COLLISIONAL EXCITATION
c
      do  nz = 1, 2
         abde = ((de*dh)*zion(nz))*pop(nz,nz)
         do n = 2, ml
c
            ll = ((nz-1)*(ml-1))+(n-1)
            call ratec(nz, n, t, rate, ey)
            rt = abde*rate
            hherat(ll) = rt
            hheen(ll) = ey
c
         enddo
      enddo
c
c    ***CALCULATES HI RECOMB. LINES IN HBRI(5)  (CASE B)
c
      call hydrec(t, de, dh)
c
c
c    ***ADDS CONTRIBUTION BY COLLIS. EXCITATION (CASE B)
c
      do  l = 1, 5
      
         hbcol = 0.0d0
         do k = l, ml-2
            frab = (caseab(1)*hbemb(l,k))+((1.0d0-caseab(1))*hbema(l,k))
            hbcol = hbcol+(frab*hherat(k+1))
         enddo
c
c  recorded for comparison with the new hydro lines, not used otherwise
c   
      hbri(l) = hbri(l)+((hbcol*(hheen(l+1)-hheen(1)))/(4.d0*pi))
c     
      
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	HBeta and 4686 are computed in hydro.f
c
c
c      hbeta = hbri(2)
c
c      *HEII 4686
c
      hbcol = ((((0.5164d0*hherat(ml+2))+(0.1876d0*hherat(ml+3)))
     &        +(0.1665d0*hherat(ml+4)))+(0.1580d0*hherat(ml+5)))+(
     &        0.1536d0*hherat(ml+6))
     
      heiilr = heiilr+((hbcol*(hheen(ml+2)-hheen(ml+1)))/(4.d0*pi))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***derives the two photon emission flux from HI (2s level)
c
c
      ab = (dh*zion(1))*pop(2,1)
      abde = ab*de
      r2q = hheen(1)*((fr(u)*hherat(1))+(abde*reff(u)))
      r2q = r2q/(1.0d0+(((ab*qpr(u))+(de*qel(u)))/8.226d0))
c
c      write(*,*) 'old ',r2q
c
      h2ql = r2q/(4.d0*pi)
c
c    ***derives 2 photon emission flux for HeII
c
      ll = ml
      ab = (dh*zion(2))*pop(3,2)
      abde = ab*de
      r2q = hheen(ll)*((fr(w)*hherat(ll))+(abde*(1.d0-x1(w))*rec(5,2)))
c
      ab = (dh*zion(1))*pop(2,1)
      r2q = r2q/(1.0d0+(((ab*qpr(w))+(de*qel(w)))/526.464))
      heii2ql = r2q/(4.d0*pi)
c
c      write(*,*)'Old He 2p',heii2ql
c
      return 
      end
