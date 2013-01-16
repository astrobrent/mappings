cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES LOCAL EMISSIVITY IN VECTORS EMIDIF,EMILIN
c	TO BE USED IN SUBROUTINE TOTPHOT AND NEWDIF
c
c	NOTE : UNITS FOR THE VECTORS EMIDIF,EMILIN ARE IN NUMBERS
c	       OF PHOTONS (INSTEAD OF ERGS LIKE ALL OTHER VECTORS
c          (Units are photons/s/cm^3/Hz/Sr (1/4pi))
c
c	CALL SUBR. RESON,INTER,HYDRO,ALLRATES
c	NOTE : CALLING THESE SUBR. COULD CHANGE THE INTENSITIES
c	       OF THE LINES IN THE BUFFER VECTORS,SO IT IS NECESSARY
c	       TO CALL SUBR. COOL BEFORE SUMMING LINE INTENS. (SUMDATA)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine localem(t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision t, de, dh, telc
      double precision xp(mxion, mxelem)
      double precision abio,at,bet,cebin,crosec,eminu,energ,enu
      double precision eph,eto,excerg
      double precision phoi,rkt,se,statwei,w
      double precision wid,xpf,zsqd
c      
      integer*4 i,ie,ii2,inl,j,jn,atom
      integer line,series,ion,trans
      character jjmod*4
c
c     Functions
c
      double precision acrs,densnum,favcha,freem,hir,q2p
c
c    ***DEFINITION OF FUNCTIONS AND CONSTANTS
c
      q2p(w) = (1.9922d-26*w)*(((9.0d0*w)-((9.0d0*w)*w))**0.55d0)
c
      hir(rkt,enu,excerg) = dexp(((-138.4461d0+(2.0d0*dlog(enu)))-
     &                          (1.5d0*dlog(rkt)))-(excerg/rkt))
c
      freem(rkt,energ) = dexp(((-106.366d0-dlog(energ))-
     &                       (0.5d0*dlog(rkt)))-(energ/rkt))
c
      acrs(eph,at,bet,se) = (at*(bet+((1.0d0-bet)/eph)))/(eph**se)
c
      telc = dmax1(t,1.d-2)
      rkt = rkb*telc
c
c
      do  j = 1, atypes
         do  i = 1, maxion(j)-1
            if (j.gt.limph) then
               xp(i,j) = 0.0d0
            else
               xp(i,j) = 1.0d0
            end if
         enddo
      enddo
c
      xpf = 1.0d0
      zsqd = (de*densnum(dh))*(favcha(pop,2) ** 2)
      do i = 1, ionum
        ii2 = i
        if ((atpho(i).eq.2).and.(ionpho(i).eq.2)) goto 75
      enddo
   75 continue
c
      jjmod = 'TEMP'
c
c    ***FINDS LINES BRIGHTNESS
c
      call allrates(t, jjmod)
      call hydro(t, de, dh)
      call inter(t,de,dh)
      call inter2(t,de,dh)
      call fivelevel(t,de,dh)
      call sixlevel(t,de,dh)
      call ninelevel(t,de,dh)
      call fine3(t,de,dh)
      call ironii(t,de,dh)
      call reson(t,de,dh)
      call reson2(t,de,dh)
      call helif(t,de,dh)
c
      jcon = 'YES'
c
      if (jcon.eq.'NO') then
c
c     original continuum calcs for jcon = 'NO'
c
         do inl = 1, infph-1
c
            eminu = 0.0d0
c
            cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
            energ = cebin*ev
            enu = cebin*evplk

c     Adds recapture to ground state for elements up to LIMPH.
c     STATWEI = Ratio of statistical weight of recombination to ion
c     EX. HYDR.: STATWEI 2/1=2
c     
            eto = 0.0d0
            do 290 i = 1, ionum
               ie      = atpho(i)
               j       = ionpho(i)
               statwei = stwtpho(i)
c
               if ((xp(j,ie).gt.0.0d0).and.(statwei.gt.0.0d0)) then
                  phoi     = ipotpho(i)
                  if (cebin.lt.phoi) goto 293
                  excerg   = ev*(cebin-phoi)
                  eph      = cebin/phoi
                  abio     = zion(ie)*pop(j+1,ie)
                  crosec   = acrs(eph,sigpho(i),betpho(i),spho(i))
                  xp(j,ie) = ((abio*crosec)*hir(rkt,enu,excerg))*statwei
                  eto      = eto+xp(j,ie)
               end if
c
 290        continue
 293        continue
c
c     Add free-free emission.
c
            if (xpf .gt. 0.0d0) then
               xpf = zsqd*freem(rkt,energ)
               eminu = eminu+xpf
            end if
c
c     Add recapture to level N=2 of He II.
c
            ie = 2
            j = 2
            jn = j+2
            statwei = 8.0d0
            if (xp(jn,ie).gt.0.0d0) then
               phoi = ipotpho(ii2)/4.0d0
               if (cebin.lt.phoi) goto 297
               excerg = ev*(cebin-phoi)
               eph = cebin/phoi
               abio = zion(ie)*pop(j+1,ie)
               crosec = 2.d0*acrs(eph,sigpho(ii2),betpho(ii2),spho(ii2))
               xp(jn,ie) = ((abio*crosec)*hir(rkt,enu,excerg))*statwei
               eto = eto+xp(jn,ie)
            end if
 297        continue
            eminu = eminu+((de*dh)*eto)
c
c     Add two-photon emission by H I.
c
            w = energ/(0.75d0*epot(1,1))
            if (w.lt.1.0d0) then
               eminu = eminu+(q2p(w)*(h2ql/energ))
            endif
c     
c     Add two-photon emission by He I.
c
            w = energ/heien(2)
            if (w.lt.1.0d0) then
               eminu = eminu+(q2p(w)*(hein(2)/energ))
            end if
c
c     Add two-photon emission by He II.
c
            w = energ/(0.75d0*epot(2,2))
            if (w.lt.1.0d0) then
               eminu = eminu+(q2p(w)*(heii2ql/energ))
            end if

            emidif(inl) = eminu

         end do                 ! End of do inl=1,infph-1 loop.

      else if (jcon .eq. 'YES') then

         do inl=1,infph
            cnphot(inl) = 0.d0
         enddo
c         
         call freefree(t,de,dh)
         call freebound(t,de,dh)
         call twophoton(t,de,dh)
c            
         do inl = 1,infph
c
            cnphot(inl) = ffph(inl)+fbph(inl)+p2ph(inl)
            emidif(inl) = ffph(inl)+fbph(inl)+p2ph(inl)
c
         enddo
c        
      endif   ! End of if (jcon .eq. Y)
c
c  If dust and IR included add dust continuum to continuum emission
c
      if (grainmode) then
       if (IRmode.ne.0) then
        do inl = 1,infph
       	  cnphot(inl) = cnphot(inl) + IRphot(inl)
          emidif(inl) = emidif(inl) + IRphot(inl)
        enddo
       endif
       if ((pahmode).and.(IRmode.ne.0)) then
        do inl = 1,infph
           cnphot(inl) = cnphot(inl)+paheng*pahflux(inl)*pahfrac*dh
           emidif(inl)=emidif(inl)+paheng*pahflux(inl)*pahfrac*dh
        enddo
       endif
      endif
c
c     Now add weak lines directly to vector EMIDIF...
c
c    Five level atoms
c
      do ion = 1,nfions
         do trans = 1,nftrans
             j = fbin(trans,ion)
             if ((j.ne.0).and.(fbri(trans,ion).gt.epsilon)) then 
c                energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
                energ = (plk*cls)/flam(trans,ion)
                wid = evplk*(ephot(j+1)-ephot(j))
                emidif(j) = emidif(j)+(fbri(trans,ion)/(wid*energ))
             endif
         enddo
      enddo
c
c    Six level atoms
c
      do ion = 1,n6ions
         do trans = 1,n6trans
             j = f6bin(trans,ion)
             if ((j.ne.0).and.(f6bri(trans,ion).gt.epsilon)) then 
c                energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
                energ = (plk*cls)/f6lam(trans,ion)
                wid = evplk*(ephot(j+1)-ephot(j))
                emidif(j) = emidif(j)+(f6bri(trans,ion)/(wid*energ))
             endif
         enddo
      enddo
c
c    Nine level atoms
c
      do ion = 1,n9ions
         do trans = 1,n9trans
             j = f9bin(trans,ion)
             if ((j.ne.0).and.(f9bri(trans,ion).gt.epsilon)) then 
c                energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
                energ = (plk*cls)/f9lam(trans,ion)
                wid = evplk*(ephot(j+1)-ephot(j))
                emidif(j) = emidif(j)+(f9bri(trans,ion)/(wid*energ))
             endif
         enddo
      enddo
c
c    Three level atoms
c
      do ion = 1,nf3ions
         do trans = 1,nf3trans
             j = f3bin(trans,ion)
             if ((j.ne.0).and.(f3bri(trans,ion).gt.epsilon)) then 
c                energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
                energ = (plk*cls)/f3lam(trans,ion)
                wid = evplk*(ephot(j+1)-ephot(j))
                emidif(j) = emidif(j)+(f3bri(trans,ion)/(wid*energ))
             endif
         enddo
      enddo
c
c     Add old intercombination lines to vector EMDIF.
c
      do i = 1, mlines
         j = lcbin(i)
         if ((j.ne.0).and.(fsbri(i).gt.epsilon)) then 
c         energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
         energ = e12fs(i)
         wid = evplk*(ephot(j+1)-ephot(j))
         emidif(j) = emidif(j)+(fsbri(i)/(wid*energ))
         endif
      enddo
c
c
      do i = 1, xilines
         j = xibin(i)
         if ((j.gt.0).and.(xibri(i).gt.epsilon)) then
c         energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
            energ = ev*xiejk(i)
            wid   = evplk*(ephot(j+1)-ephot(j))
            emidif(j) = emidif(j)+(xibri(i)/(wid*energ))
         endif
      enddo
c
c old intercombination lines
c     
      do i = 1, nlines
         j = lrbin(i)
         if ((j.ne.0).and.(rbri(i).gt.epsilon)) then 
c         energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
          energ = e12r(i) 
         wid = evplk*(ephot(j+1)-ephot(j))
         emidif(j) = emidif(j)+rbri(i)/(wid*energ)
         endif
      enddo
c
      do i = 1, nheilines
         j = heibin(i)
         if ((j.ne.0).and.(heibri(i).gt.epsilon)) then 
c         energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
         energ = ev*(LmeV/(heilam(i)*1.d8))
         wid = evplk*(ephot(j+1)-ephot(j))
         emidif(j) = emidif(j)+heibri(i)/(wid*energ)
         endif
      enddo
c
c     Add He-like ion intercombination and forbidden lines to EMIDIF.
c
      do  i = 1, xhelines
         j = xhebin(i)
c         energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
         energ = ev*xhejk(i)
         wid = evplk*(ephot(j+1)-ephot(j))
         emidif(j) = emidif(j)+(xhebri(i)/(wid*energ))
      enddo
c
c     adds xlines to emilin(1,line)
c
      do line = 1,xlines
         j = xbin(line)
c         energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
         energ = ev*xejk(line)
         wid = evplk*(ephot(j+1)-ephot(j))
         emilin(1,line) = xrbri(line)/(wid*energ)
      enddo
c
c
c     adds hydrogen and helium lines to hydlin(1,line,series) and 
c      hellin(1,line,series)
c
      do series = 1,6
      do line = 1,10
         j = hbin(line,series)         
         if (j.ne.0) then
c            energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
            energ = ev*LmeV/hlambda(line,series)
            wid = evplk*(ephot(j+1)-ephot(j))
            hydlin(1,line,series) = hydrobri(line,series)/(wid*energ)
         endif
      enddo
      enddo
c
      do series = 1,6
      do line = 1,10
         j = hebin(line,series)
         if (j.ne.0) then
c            energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
            energ = ev*LmeV/helambda(line,series)
            wid = evplk*(ephot(j+1)-ephot(j))
            hellin(1,line,series) = helibri(line,series)/(wid*energ)
         endif
      enddo
      enddo
c
c     Adds heavy hydrogenic series
c
      if (atypes.gt.3) then
      do atom = 3,atypes
c
c      if (xhydrobri(1,1,atom).gt.epsilon) then 
c      only include if xhydlin is zeroed (ie through zerbuf)
c
      do series = 1,2
      do line = 1,10
         j = xhbin(line,series,atom)
         if (j.ne.0) then
c            energ = 0.5d0*(ev*(ephot(j+1)+ephot(j)))
            energ = ev*LmeV/xhlambda(line,series,atom)
            wid = evplk*(ephot(j+1)-ephot(j))
            xhydlin(1,line,series,atom) = xhydrobri(line,series,atom) 
     &                                    /(wid*energ)
         endif
      enddo
      enddo
c      endif
      enddo
      endif
c
      return
c 
      end
