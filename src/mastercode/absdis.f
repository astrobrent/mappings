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
c*******COMPUTES THE DISTANCE : DRTA  TO OBTAIN
c	A GIVEN TOTAL PHOTON ABSORBTION FRACTION FROM TPHOT
C       TAUAV (MMOD='DIS','LIN') takes account of dust
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine absdis(t, dh, fi, absf, drta, rad, popul)
c
c
      include 'cblocks.inc'
c
c
c
      double precision popul(mxion, mxelem), popt(mxion, mxelem)
      double precision wid,wei,phots
      double precision t, dh, fi, absf, drta, rad,dstmin,pcros
      double precision plos,qto, de, g1, g2, temp
      double precision nfn(mxinfph),r2,r3,cebin,sig,ue
      integer*4 bincount,inl,i,ie,j,k,dtype
c
      double precision ptime,rtime,ctime,abio,crosec
c
c       External Functions
c
      double precision feldens,fphotim,frectim2,fcolltim,fdilu
c
c       Internal Functions
c
      double precision acrs,eph,at,bet,se
c
      acrs(eph,at,bet,se) = (at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
c
c    ***DERIVATION OF TOTAL CROSS SECTIONS AT EACH ENERGY BIN
c
      if (drta.eq.0.d0) drta = 1.0d16
c
c
      call copinto(pop,popt)
      call copinto(popul,pop)
c
      de = feldens(dh, popul)
      ue = 1.5*(de+dh*zen)*rkb*t
c
c      blum = 0.d0
c      do inl = 1, infph-1
c         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
c         blum = blum+tphot(inl)*cebin*evplk
c      enddo
c
      ptime = fphotim()
      rtime = frectim2(de)
      ctime = fcolltim(de)
c
c      write(*,*) ptime,rtime,ctime
c      r4 = rtime/(ptime+rtime)
c      write(*,*) 'Recombination Correction : ',r4
c      absf = absf*r4
c      write(*,*) 'Absorbsion Fraction : ',absf
c
      call copinto(popt,pop)
c
      plos = 0.d0
      xsect = 0.d0
      bincount = 0
      qto = 0.d0
      do inl = 1, infph-1
         xsec(inl)=0.d0
         if (skipbin(inl)) goto 300
         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
         nfn(inl) = tphot(inl)*cebin*evplk
         sig = 0.0d0
         do i = 1, ionum
            ie = atpho(i)
            j = ionpho(i)
            abio = zion(ie)*popul(j,ie)
            eph = cebin/ipotpho(i)
            if (eph.lt.1.0d0) goto 400
            crosec = acrs(eph,sigpho(i),betpho(i),spho(i))
            sig = sig+(abio*crosec)
         enddo
 400     continue
c
         if (sig.lt.0.d0) sig = 0.d0
         sig = (dh*fi)*sig
         xsec(inl) = sig
         xsect = xsect+sig
 300     continue
      enddo
c
c  dust
c      
      if(grainmode)then
c
c   pahs
c
        if (pahmode) then
         do inl=1,infph-1
           sig=0.d0
           do k=1,pahi
             if (pahion(k).gt.0) then
                pcros = pahZ(k)*pahiext(inl)*pahfrac
             else
                pcros = pahZ(k)*pahnext(inl)*pahfrac
             endif
             sig = sig+pcros
           enddo
	   xsec(inl)=xsec(inl)+sig*dh*fi
	   xsect=xsect+sig*dh*fi
        enddo
       endif         
c
c  grains
c 
        dstmin=1
        if(ClinPAH.AND.(.NOT.pahmode)) dstmin=2
c
       do dtype=dstmin,numtypes
         do inl=1,infph-1
c          if (skipbin(inl)) goto 310
          sig=fi*dh*dcrosec(inl,dtype)
          xsec(inl)=xsec(inl)+sig
          xsect=xsect+sig
 310      continue
         enddo
       enddo
      endif
c
c      write(*,*) 'Total X-section : ',xsect
c
      do inl = 1, infph-1
c        if (.NOT.skipbin(inl)) then
         wei = xsec(inl)/xsect
         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
         wid = (ephot(inl+1)-ephot(inl))*evplk
         phots = (tphot(inl))*wid*wei
c         write(*,*) inl,phots,qto
c         write(*,*) tphot(inl),xsec(inl),xsect
         qto = qto+phots
c        endif
      enddo
c
c      write(*,*) 'Weighted total : ',qto
c
      iter = 0
 100  continue
      plos = 0.d0
c
      do inl = 1, infph-1
c        if (.NOT.skipbin(inl)) then
         wei = xsec(inl)/xsect
         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
         wid = (ephot(inl+1)-ephot(inl))*evplk
         phots = (tphot(inl))*wid*wei
c
         plos = plos+(1.d0-dexp(-drta*xsec(inl)))*(phots/qto)
c
      enddo
c
      iter = iter+1
      if (plos.gt.0.d0) then
        r2 = absf/plos
        r3 = dabs(1.d0-r2)
c        write(*,*) 'Change in dr:',r2
        if (r3.gt.0.01d0) then
         drta = drta*r2
         goto 100
        endif
      else
        r2 = 1.d0/epsilon
        drta = r2
      endif
c
c
c     diagnosic printouts, normally commented out.
c
c
c      caller = 'AD'
c      pfx = 'xsect'
c      np = 5
c      wmod = 'REAL'
c      call wpsou(caller,pfx,np,wmod,
c     &        t,de,dh,
c     &        0.d0,0.d0,vshoc,
c     &        0.d0,drta,
c     &        0.d0,0.d0,
c     &        1.d0,xsec)
c      caller = 'AD'
c      pfx = 'nfn'
c      np = 3
c      wmod = 'REAL'
c      call wpsou(caller,pfx,np,wmod,
c     &        t,de,dh,
c     &        0.d0,0.d0,vshoc,
c     &        0.d0,drta,
c     &        0.d0,0.d0,
c     &        1.d0,nfn)
c      caller = 'AD'
c      pfx = 'tphot'
c      np = 5
c      wmod = 'REAL'
c      call wpsou(caller,pfx,np,wmod,
c     &        t,de,dh,
c     &        0.d0,0.d0,vshoc,
c     &        0.d0,drta,
c     &        0.d0,0.d0,
c     &        1.d0,tphot)
c     
c
      if (plos.gt.0.d0) drta = drta*absf/plos
c
c     correct for geometric dilution changes if necessary
c
      if ((jgeo.eq.'S').or.(jgeo.eq.'F')) then 
         g1 = fdilu(rstar,rad)
         r2 = rad + drta
         g2 = fdilu(rstar,r2)
         temp = dsqrt(1.23456789d0*g2/g1)
         r2 = (r2*temp)-rad
         if (r2.lt.0) r2=drta
         r3 = 1.d0/r2 + 1.d0/drta
         drta = 1.d0/r3
c         write(*,*) 'Geometric Correction:',temp
      endif

      if (drta.lt.0.d0) then
       write(*,*) "dr is negative!!",drta
       stop
      endif

c
c      drta = max(1.0e13, drta)
c
c      drta = min(5.0e14,drta)
c
      return 
      end






