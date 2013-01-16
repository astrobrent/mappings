cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c     
c     copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c     Luc Binette and Brent Groves
c     
c     Version 1.0.0r
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Computes HI and HeII cooling by collisional excitation.
c     Computes HI and HeII recombination spectrum.
c     Uses intermediate case A - B for all calculations.
c     
c     RSS1995/6
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine hydro(t, de, dh)
c     
      include 'cblocks.inc'
c     
      double precision t,de,dh,f,omg,rate,aa,totcoll
      integer*4 i,j,k,l,m, nz,clevel,line,series,atom
c     
      double precision tmin, tmax, telec,tz,nz2,z3
      double precision meanden, lmden, tlo, dlo
      double precision logrh0, logrh1, logth0, logth1
      double precision cab,uomg
      double precision d1b,d2b,d3b
      double precision hb1b,hb2b,hb3b,hb4b      
      double precision d1a,d2a,d3a
      double precision hb1a,hb2a,hb3a,hb4a     
      double precision egj,helos,pcoll
      double precision ab,ab0,ab1,abde,r2q,pz
c
      double precision qpr,qel,u
c
c     internal functions
c
      qpr(u) = 4.74d-4*(u ** (-0.151d0))
      qel(u) = 0.57d-4*(u ** (-0.373d0))
c     
c     init
c     
      hloss = 0.0d0
c     
      totcoll = 0.d0
c     
      do atom = 1,atypes
         do j = 1,16
            rateton(atom,j) = 0.0d0
         enddo
         recrate2p(atom) = 0.0d0
         collrate2p(atom) = 0.0d0
      enddo
c     
c      write(*,*) ' H Case A-B (0-1) in hydro2:', caseab(1)
c      write(*,*) ' He Case A-B (0-1) in hydro2:', caseab(2)
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Get total collisional excitation up to first 15 levels H and He II 
c     etc.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      do atom = 1, atypes
c     
         nz = mapz(atom)
         nz2 = dble(nz*nz)
c     
         u = dabs(t+epsilon)/nz2
         uomg = u
c     
         if (uomg.ge.5.0d5) uomg = 5.0d5
c
         f = dsqrt(1.0d0/(t+epsilon))
c     
         ab0 = zion(atom)*pop(nz,atom)
         ab1 = zion(atom)*pop(nz+1,atom)
         abde = ab0*de*dh
c     
         do j = 1,20 
c     
c     ground level excitation      
c     energy of upper level
c     
            egj = epot(nz,atom)*(1.d0-(1.d0/(colid(j,3)*colid(j,3))))
            aa  = egj/(rkb*t)
c     
            if (aa.lt.loghuge) then
c     
               if (u.lt.7.2d4) then
                  omg = ccoln(j,1)+uomg*(ccoln(j,2)+uomg*(ccoln(j,3)
     &                  +uomg*ccoln(j,4)))
               else
                  omg = dcoln(j,1)+uomg*(dcoln(j,2)+uomg*(dcoln(j,3)
     &                  +uomg*dcoln(j,4)))
               endif
c     
               rate    = rka*f*dexp(-aa)*(omg/nz2)
               totcoll = totcoll+rate*egj*abde
c     
               if (j.eq.1) then
c                  write(*,*) '2S excitation rate:',elem(atom), rate
                  collrate2p(atom) = rate
               endif
c     
               rateton(atom,colid(j,3)) = rateton(atom,colid(j,3))
     &                                    +rate*abde
c     
            endif
c     
         enddo
      enddo
c     
c      write(*,'(" T De DH ",3(1pg14.7,1x))'), t, de, dh
c      write(*,'(" H N=2 rate:",1pg14.7)'), rateton(1,2) 
c      write(*,'(" H N=3 rate:",1pg14.7)'), rateton(1,3) 
c      write(*,'(" H N=4 rate:",1pg14.7)'), rateton(1,4) 
c     
      hloss = totcoll
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Calculate H I and He II recombination lines. 
c     
c     Based on Hummer and Storey 1995.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     look up densities
c     
      i = 1
      j = 7
      
 100  k = (i+j)/2
      
      if (de.ge.rhhe(k)) i = k
      if (de.lt.rhhe(k)) j = k
      if ((j-i).ge.2) goto 100 
      
      if (j.eq.1) j = 2
      if (i.eq.7) i = 6
      if (i.eq.j) i = j-1
c     
      logrh0 = dlog(rhhe(i))
      logrh1 = dlog(rhhe(j))-logrh0
      logrh0 = (dlog(de+1.d-10)-logrh0)/(logrh1+epsilon)
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	Check temperatures
c     
      tmin = 100.d0
      tmax = 1.0d6
c     
      telec = t+1.d-3
      hbeta = 0.d0
      heiilr = 0.d0
      do series = 1,6
         do line = 1,10
           hydrobri(line,series) = 0.d0
           helibri(line,series) = 0.d0
         enddo
      enddo	
c
      if ((t.ge.tmin).and.(t.le.tmax)) then 
c      
      tlo = dlog(telec)
      dlo = dlog(de+1.d-10)
      meanden = dsqrt(de*dh)
      lmden = dlog(meanden+1.d-10)
c     
c     
c     ***FINDS ABSOLUTE FLUX FOR H-BETA
c     
      abde = ((de*dh)*zion(1))*pop(2,1)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c   Look up H-Beta temperature
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      l = 1
      m = 10
c     
 101  k = (l+m)/2
c     
      if (telec.ge.th42(k)) l = k
      if (telec.lt.th42(k)) m = k
      if ((m-l).ge.2) goto 101 
c     
      if (m.eq.1) m = 2
      if (l.eq.10) l = 9
      if (l.eq.m) l = m-1
c     
      logth0 = dlog(th42(l))
      logth1 = dlog(th42(m))-logth0
      logth0 = (tlo-logth0)/(logth1+epsilon)
c     
c     Interpolate Hbeta 
c     
      d1b = dlog(h42b(l,i))+(dlog(h42b(m,i))-dlog(h42b(l,i)))*logth0
      d2b = dlog(h42b(l,j))+(dlog(h42b(m,j))-dlog(h42b(l,j)))*logth0
c     
      d1a = dlog(h42a(l,i))+(dlog(h42a(m,i))-dlog(h42a(l,i)))*logth0
      d2a = dlog(h42a(l,j))+(dlog(h42a(m,j))-dlog(h42a(l,j)))*logth0
c     
      d3b = d1b+(d2b-d1b)*logrh0
      d3a = d1a+(d2a-d1a)*logrh0
c     
      hbeta = (caseab(1)*dexp(d3b)+(1.0-caseab(1))*dexp(d3a))
     &        *abde/(4.d0*pi)
c     
c     Do H series
c     
      do series = 1,6
         do line = 1,10
c     
            hb1b = dlog(hylratsb(line,series,l,i)+epsilon)
            hb2b = dlog(hylratsb(line,series,m,i)+epsilon)
            hb3b = dlog(hylratsb(line,series,l,j)+epsilon)
            hb4b = dlog(hylratsb(line,series,m,j)+epsilon)
c     
            d1b = hb1b+(hb2b-hb1b)*logth0
            d2b = hb3b+(hb4b-hb3b)*logth0
c     
            hb1a = dlog(hylratsa(line,series,l,i)+epsilon)
            hb2a = dlog(hylratsa(line,series,m,i)+epsilon)
            hb3a = dlog(hylratsa(line,series,l,j)+epsilon)
            hb4a = dlog(hylratsa(line,series,m,j)+epsilon)
c     
            d1a = hb1a+(hb2a-hb1a)*logth0
            d2a = hb3a+(hb4a-hb3a)*logth0
c     
            d3b = d1b+(d2b-d1b)*logrh0
            d3a = d1a+(d2a-d1a)*logrh0
c     
            hydrobri(line,series) = (caseab(1)*dexp(d3b)+
     &                            (1.d0-caseab(1))*dexp(d3a))*hbeta
c     
         enddo
      enddo	
c     
c     Interpolate 2S0 recombinate rate
c     
      d1b = dlog(r2s1b(l,i))+(dlog(r2s1b(m,i))-dlog(r2s1b(l,i)))*logth0
      d2b = dlog(r2s1b(l,j))+(dlog(r2s1b(m,j))-dlog(r2s1b(l,j)))*logth0
c     
      d1a = dlog(r2s1a(l,i))+(dlog(r2s1a(m,i))-dlog(r2s1a(l,i)))*logth0
      d2a = dlog(r2s1a(l,j))+(dlog(r2s1a(m,j))-dlog(r2s1a(l,j)))*logth0
c     
      d3b = d1b+(d2b-d1b)*logrh0
      d3a = d1a+(d2a-d1a)*logrh0
c     
      recrate2p(1) = (caseab(1)*dexp(d3b)+(1.0-caseab(1))*dexp(d3a))
c
c      write(*,*) 'H rec 2p',recrate2p(1)
c
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Look up He II 4686 temperature and do He II series.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      abde = ((de*dh)*zion(2))*pop(3,2)
c     
      l = 1
      m = 12
c     
 102  k = (l+m)/2
c     
      if (telec.ge.the43(k)) l = k
      if (telec.lt.the43(k)) m = k
      if ((m-l).ge.2) goto 102 
c     
      if (m.eq.1) m = 2
      if (l.eq.12) l = 11
      if (l.eq.m) l = m-1
c     
      logth0 = dlog(the43(l))
      logth1 = dlog(the43(m))-logth0
      logth0 = (tlo-logth0)/(logth1+epsilon)
c     
      d1b = dlog(he43b(l,i))+(dlog(he43b(m,i))-dlog(he43b(l,i)))*logth0
      d2b = dlog(he43b(l,j))+(dlog(he43b(m,j))-dlog(he43b(l,j)))*logth0
c     
      d3b = d1b+(d2b-d1b)*logrh0
c     
      d1a = dlog(he43a(l,i))+(dlog(he43a(m,i))-dlog(he43a(l,i)))*logth0
      d2a = dlog(he43a(l,j))+(dlog(he43a(m,j))-dlog(he43a(l,j)))*logth0
c     
      d3a = d1a+(d2a-d1a)*logrh0
c     
      heiilr = (caseab(2)*dexp(d3b)+(1.0-caseab(2))*dexp(d3a))
     &         *abde/(4.d0*pi)
c     
      do series = 1,6
         do line = 1,10
c     
            hb1b = dlog(helratsb(line,series,l,i)+epsilon)
            hb2b = dlog(helratsb(line,series,m,i)+epsilon)
            hb3b = dlog(helratsb(line,series,l,j)+epsilon)
            hb4b = dlog(helratsb(line,series,m,j)+epsilon)
c     
            d1b = hb1b+(hb2b-hb1b)*logth0
            d2b = hb3b+(hb4b-hb3b)*logth0
c     
            hb1a = dlog(helratsa(line,series,l,i)+epsilon)
            hb2a = dlog(helratsa(line,series,m,i)+epsilon)
            hb3a = dlog(helratsa(line,series,l,j)+epsilon)
            hb4a = dlog(helratsa(line,series,m,j)+epsilon)
c     
            d1a = hb1a+(hb2a-hb1a)*logth0
            d2a = hb3a+(hb4a-hb3a)*logth0
c     
            d3b = d1b+(d2b-d1b)*logrh0
            d3a = d1a+(d2a-d1a)*logrh0
c     
            helibri(line,series) = (caseab(2)*dexp(d3b)
     &                           +(1.0-caseab(2))*dexp(d3a))*heiilr
c     
         enddo
      enddo	
c     
c     Interpolate 2S0 recombinate rate
c     
      d1b = dlog(r2s2b(l,i))+(dlog(r2s2b(m,i))-dlog(r2s2b(l,i)))*logth0
      d2b = dlog(r2s2b(l,j))+(dlog(r2s2b(m,j))-dlog(r2s2b(l,j)))*logth0
c     
      d1a = dlog(r2s2a(l,i))+(dlog(r2s2a(m,i))-dlog(r2s2a(l,i)))*logth0
      d2a = dlog(r2s2a(l,j))+(dlog(r2s2a(m,j))-dlog(r2s2a(l,j)))*logth0
c     
      d3b = d1b+(d2b-d1b)*logrh0
      d3a = d1a+(d2a-d1a)*logrh0
c     
      recrate2p(2) = (caseab(2)*dexp(d3b)+(1.0-caseab(2))*dexp(d3a))
c
c      write(*,*) 'He Rec 2p',recrate2p(2)
c
c     end t limited calcs
c
      endif     
c
      if (atypes.gt.2) then
c
c     Do two series for heavy atoms, based on extrapolation of case A
c      hydrogen
c
      do atom = 3,atypes
        do series = 1,2
         do line = 1,10
           xhbeta(atom) = 0.d0
           xhydrobri(line,series,atom) = 0.d0
         enddo
        enddo	
      enddo	
c      
      do atom = 3,atypes
c      
c     ***FINDS ABSOLUTE FLUX FOR Heavy Z H-BETA
c     
      abde = zion(atom)*pop(mapz(atom)+1,atom)
      
      if (abde.gt.1e-12) then 
c
       abde = abde*de*dh
      
       telec = (t+1.d-3)/(mapz(atom)*mapz(atom))
       z3    = mapz(atom)*mapz(atom)*mapz(atom)
c 
       if ((telec.ge.tmin).and.(telec.le.tmax)) then 
c
        tlo = dlog(telec)
        dlo = dlog(de+1.d-10)
        meanden = dsqrt(de*dh)
        lmden = dlog(meanden+1.d-10)
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Look up equivalent H-Beta temperature
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        l = 1
        m = 10
c     
 103    k = (l+m)/2
c     
        if (telec.ge.th42(k)) l = k
        if (telec.lt.th42(k)) m = k
        if ((m-l).ge.2) goto 103 
c     
        if (m.eq.1) m = 2
        if (l.eq.10) l = 9
        if (l.eq.m) l = m-1
c     
        logth0 = dlog(th42(l))
        logth1 = dlog(th42(m))-logth0
        logth0 = (tlo-logth0)/(logth1+epsilon)
c     
c     Interpolate high Z Hbeta 
c     
        d1b = dlog(h42b(l,i))+(dlog(h42b(m,i))-dlog(h42b(l,i)))*logth0
        d2b = dlog(h42b(l,j))+(dlog(h42b(m,j))-dlog(h42b(l,j)))*logth0
c     
        d1a = dlog(h42a(l,i))+(dlog(h42a(m,i))-dlog(h42a(l,i)))*logth0
        d2a = dlog(h42a(l,j))+(dlog(h42a(m,j))-dlog(h42a(l,j)))*logth0
c     
        d3b = d1b+(d2b-d1b)*logrh0
        d3a = d1a+(d2a-d1a)*logrh0
c     
        xhbeta(atom) = 
     &    z3*(caseab(atom)*dexp(d3b)+(1.0-caseab(atom))*dexp(d3a))
     &    *abde/(4.d0*pi)
c     
c     Do H series by ratio
c     
        do series = 1,2
          do line = 1,10
c     
            hb1b = dlog(hylratsb(line,series,l,i)+epsilon)
            hb2b = dlog(hylratsb(line,series,m,i)+epsilon)
            hb3b = dlog(hylratsb(line,series,l,j)+epsilon)
            hb4b = dlog(hylratsb(line,series,m,j)+epsilon)
c     
            d1b = hb1b+(hb2b-hb1b)*logth0
            d2b = hb3b+(hb4b-hb3b)*logth0
c     
            hb1a = dlog(hylratsa(line,series,l,i)+epsilon)
            hb2a = dlog(hylratsa(line,series,m,i)+epsilon)
            hb3a = dlog(hylratsa(line,series,l,j)+epsilon)
            hb4a = dlog(hylratsa(line,series,m,j)+epsilon)
c     
            d1a = hb1a+(hb2a-hb1a)*logth0
            d2a = hb3a+(hb4a-hb3a)*logth0
c     
            d3b = d1b+(d2b-d1b)*logrh0
            d3a = d1a+(d2a-d1a)*logrh0
c     
            xhydrobri(line,series,atom) = 
     &       (caseab(atom)*dexp(d3b)+(1.d0-caseab(atom))*dexp(d3a))
     &       *xhbeta(atom)
c     
          enddo
        enddo	
       endif           ! End of IF in temperature range.

      endif              ! End of abundance.

      enddo                 ! End of loop over atoms.

      endif                    ! End of have heavies.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Extrapolate high Z 2S rec rates by Z and T/Z^2
c     Based on H data.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (atypes.gt.2) then
         do atom = 3, atypes
c     
            nz = mapz(atom)
            nz2 = dble(nz*nz)
            telec = t+1.d-3
            tz = telec/nz2
            tlo = dlog(tz)
            
            if ((tz.ge.tmin).and.(tz.le.tmax)) then 
c     
            l = 1
            m = 10
c     
 104        k = (l+m)/2
c     
            if (tz.ge.th42(k)) l = k
            if (tz.lt.th42(k)) m = k
            if ((m-l).ge.2) goto 104 
c     
            if (m.eq.1) m = 2
            if (l.eq.10) l = 9
            if (l.eq.m) l = m-1
c     
            logth0 = dlog(th42(l))
            logth1 = dlog(th42(m))-logth0
            logth0 = (tlo-logth0)/(logth1+epsilon)
c     
c     Interpolate 2S0 recombinate rate
c     
            d1a = dlog(r2s1a(l,i))+(dlog(r2s1a(m,i))-dlog(r2s1a(l,i)))
     &            *logth0
            d2a = dlog(r2s1a(l,j))+(dlog(r2s1a(m,j))-dlog(r2s1a(l,j)))
     &            *logth0
c     
            d3a = d1a+(d2a-d1a)*logrh0
c     
c     Assume case A for higher species, this may not hold
c     in metal rich models
c     
            recrate2p(atom) = nz*dexp(d3a)
c     
c            write(*,*) 'extrap',tz,elem(atom),recrate2p(atom)
c     
            endif
c            
         enddo     
      endif
c     
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Add collisional excitation to recombination spectrum of H and He 
c     only.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      do atom = 1, 2
c     
         nz = mapz(atom)
c     
         do series = 1,6
            do line = 1,10
c     
               k = series+line
c     
               egj = ev*Ryd*((1.0/(series*series))-(1.0/(k*k)))*nz*nz
c     
               do clevel = 2,15                  
c     
                  if (clevel.ge.k) then
c     
                     cab = caseab(atom)*collhb(line,series,clevel)
     &                    +(1.0-caseab(atom))*collha(line,series,clevel)
                     
                     if (nz.eq.1) then
                        hydrobri(line,series) = hydrobri(line,series)
     &                       +cab*rateton(atom,clevel)*egj/(4.d0*pi)
                     else
                        helibri(line,series) = helibri(line,series)
     &                       +cab*rateton(atom,clevel)*egj/(4.d0*pi)
                     endif
c     
                  endif
c     
               enddo
c     
            enddo
         enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Adds collisional excitation to recombination spectrum of heavies.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      if (atypes.gt.3) then 
      do atom = 3, atypes
c     
         nz = mapz(atom)
c     
         do series = 1,2
            do line = 1,10
c
            if (xhydrobri(line,series,atom).gt.epsilon) then
c     
               k = series+line
c     
               egj = ev*Ryd*((1.0/(series*series))-(1.0/(k*k)))*nz*nz
c     
               do clevel = 2,15                  
c     
                  if (clevel.ge.k) then
c     
                     cab = 
     &       caseab(atom)*collhb(line,series,clevel)+
     &       (1.0-caseab(atom))*collha(line,series,clevel)
                     
                     xhydrobri(line,series,atom) = 
     &        xhydrobri(line,series,atom)+
     &        cab*rateton(atom,clevel)*egj/(4.d0*pi)
c     
                  endif
c     
               enddo
c
            endif
            enddo
         enddo
c
      enddo
      endif
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Set global H-beta and He II 4686 to include the collisional rates.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      hbeta = hydrobri(2,2)
      heiilr = helibri(1,3)
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Find two-photon emission, calling old hydro spec calculations in
c     the process, but not using them. They are needed for HeI though.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call hydro2p(t, de, dh)
c     
      do atom = 1,atypes 
c     
c     Test new 2 photon calcs
c     
         nz    = mapz(atom)
         pz    = zion(atom)*pop(nz+1,atom)
         ab1   = dh*pz
         ab0   = dh*zion(atom)*pop(nz,atom)
         abde  = ab0*de
c     
         r2q  = 0.75d0*epot(nz,atom)*de*(ab0*collrate2p(atom)+ab1
     &          *recrate2p(atom))
c         
c
c     scale temperature 
c     
         tz   = t/(1.0d4*nz*nz)
c     
c     proton & electron deexcitation collision
c     
         ab    = dh*zion(1)*pop(2,1)
         pcoll = (ab*qpr(tz)+de*qel(tz))*nz
c     
c     A scales as nz^6
c     
         nz   = nz*nz
         nz   = nz*nz*nz
c
         r2q  = r2q/(1.0d0+((ab*qpr(tz)+de*qel(tz))/(8.227*nz)))
c
         h2qbri(atom) = r2q/(4.d0*pi)
         hloss = hloss + (0.75d0*epot(mapz(atom),atom)*de*ab0
     &           *collrate2p(atom))
c         
c         if (atom.eq.2) then
c             write(*,*)'He 2p',h2qbri(atom)
c         endif
c     
      enddo
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Find intensity of important He I lines and the resultant
c     cooling rate (uses results from old hydro spectrum code).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call helioi(t, de, dh, helos)
c     
      hloss = hloss+helos
c     
      return 
c     
      end


