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
c     version 2.2, extended ionisation stages
c     Vector subsampling optimisation
c     rationalised loops
c     new charge exchange
c
c     RSS 4/91
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******DERIVES IONIC ABUNDANCES AFTER TIME : TSTEP
c	FOR THE ATOMIC ELEMENTS HEAVIER THAN HYDROGEN .
c	RETURNS RESULTS IN POP(6,11) AND DERIVATIVES IN DNDT(6,11)
c
c	WHEN ARGUMENT MOD SET TO 'EQUI' ,IT MAKES SURE THAT THE
c	SYSTEM REACHES EQUILIBRIUM (TSTEP IRRELEVANT)
c	THE ARGUMENT NEL DETERMINE IF IT COMPUTES FOR ALL ELEMENTS
c	OR ONLY FOR HELIUM .
c
c	THE AVERAGE ELECTRONIC DENSITY : DE  MUST BE GIVEN AS WELL
c	AS THE FRACTIONAL IONISATION OF HYDROGEN : XHY (OR POP2,1))
c	DEFAULT VALUE (USING POP(2,1)) ASSIGNED TO XHY IF SET TO -1
c	CALL SUBROUTINES SDIFEQ,ALLRATES,SPOTAP,IONAB,IONSEC
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine iobal(mod, nel, de, dh, xhy, t, tstep)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision abio,abne,de2,dep,dey,dpio,drec,dyn
      double precision fhi,fhii
      double precision ph,ph2,sigab,timt,yh
      double precision de, dh, xhy, t, tstep
c
      double precision reco(29), pion(29), ab(29), adndt(29)
      double precision pionau(29),rech(29),pich(29)
      double precision rc(29), pn(29), a(29), adnt(29), pa(29)
c
      integer*4 i,iel,ies,j,jdo,ji,jjj
      integer*4 nom,mnde,nde,mxde,mdelta,ion,at
c
      character jjmod*4, mod*4, nel*4
c
c           Functions
c
      double precision feldens
c
      if (mod .eq. 'EQUI') then
          timt = 1.d37
      else
          timt = tstep
          if (timt.lt.0.0d0) then
              write(*, 200) 
  200         format(/' STOP !!!!!!!!!!!!!!!!!! NEGATIVE TIME STEP'//)
              stop 
          end if
      end if
c
      fhii = pop(2,1)
      fhi = pop(1,1)
      if ((xhy.lt.0.0d0).or.(xhy.gt.1.0d0)) goto 327
      fhii = xhy
      fhi = 1.0d0-xhy
  327 continue
c
c
c
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      jjmod = 'ALL'
      call allrates(t, jjmod)
c
c
c    ***CALCULATE RATES DUE TO SECONDARY IONISATION
c
      call ionsec
c
      dyn = 1.d12
c      dyn = 1.d10
      jjj = atypes
      if (nel .eq. 'HE') jjj = 2
c
      do 100 iel = 2, jjj
      if (zion(iel).le.0.0d0) goto 100
      nde = maxion(iel)
      nom = 0
      jdo = 1
c
c    ***SET INITIAL ABUNDANCES IN AB(6) FOR ELEMENT : IEL
c
      do i = 1, nde
         ab(i) = pop(i,iel)*zion(iel)*dh
      enddo
c
c     reentry point.
c
  111 continue
c
c
c
c
      do i = 1,29
         rech(i) = 0.d0
         pich(i) = 0.d0
         reco(i) = 0.d0
         pion(i) = 0.d0
         pionau(i)=0.d0
      enddo
c
c    ***SET POPULATION RATES FOR EACH IONIC SPECIES OF ELEM : IEL
c
      do i = 2, nde
        reco(i-1) = rec(i,iel)*de
      enddo
c
c
c    ***SET DEPOPULATION RATES FOR EACH IONIC SPECIES OF ELEM : IEL
c
      do i = 1, nde-1
         pion(i) = rphot(i,iel)+(col(i,iel)*de)
      enddo
c
c     Auger rates, off diagnal terms.
c
      do i = 1, nde-2
         pionau(i) = auphot(i,iel)
      enddo
c
c    ***ADD CORRECTION DUE TO SECONDARY IONISATION
c
      do i = 1, min0(3,nde-1)
         pion(i) = pion(i)+rasec(i,iel)
      enddo
c
c    ***ADD CHARGE EXCHANGE RATES WITH HYDROGEN
c
      if (chargemode.eq.1) then 
c
c     old mapping charge rates
c
         do 40 j = 1, nchxold
            if (idint(charte(1,j)) .ne. iel) goto 40
            ji = idint(charte(2,j))
            if (reco(ji).le.0.0d0) goto 40
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c     He III not considered
c
c
            ies = idint(charte(5,j))
            abne = (dh*zion(ies))*pop(1,ies)
            abio = (dh*zion(ies))*pop(2,ies)
            drec = abne*charte(4,j)
            dpio = abio*charte(3,j)
c
            rech(1) = reco(ji)+drec
            pich(1) = pion(ji)+dpio
c
c      *COMPRESS RECH TO FIT INTO DYNAMIC RANGE (IF NECESSARY)
c
            if (reco(ji).ge.(rech(1)/dyn)) goto 45
            if ((mod .eq. 'EQUI').and.(nom.lt.1)) then
               jdo = 2
               goto 40
            end if
            reco(ji) = reco(ji)*dyn
            pion(ji) = pich(1)*(reco(ji)/rech(1))
c     
            goto 40
 45         pion(ji) = pich(1)
            reco(ji) = rech(1)
 40      continue
c
      endif 
c
c     end old reactions
c
      if (chargemode.eq.0) then 
c
c     MAPPINGS 2.2 charge rates
c
c     First, recombination reactions
c     with neutral H and He
c
         do at = 1, nchxr
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c     He III not considered
c
c
            if (chxrat(at).eq.iel) then
               if (chxr(at).gt.0.d0) then
                  ies = chxrx(at)
                  ion = chxrio(at)-1
                  abne = (dh*zion(ies))*pop(1,ies)
                  drec = abne*chxr(at)
c                  if (mapz(iel).eq.14) then
c                     write(*,*) 'drec',elem(iel),rom(ion),drec
c                     write(*,*) 'drec',abne,pop(ion,iel),rech(ion)
c                  endif
                  rech(ion) = rech(ion)+drec
               endif
            endif
c
c     end at
c
         enddo
c
c     Ionising reactions
c     with ionised H and He
c
c
         do at = 1, nchxi
c
c     abne = abund of neutral H or He
c     abio = abund of single ionised H or He
c
c     He III not considered
c
c
            if (chxiat(at).eq.iel) then
               if (chxi(at).gt.0.d0) then
                  ies = chxix(at)
                  ion = chxiio(at)
                  abio = (dh*zion(ies))*pop(2,ies)
                  dpio = abio*chxi(at)
c                  write(*,*) 'dpio',elem(iel),rom(ion),dpio
                  pich(ion) = pich(ion)+dpio
               endif
            endif
c
c     end at
c
         enddo
c
c     renormalise to fit dynamic range
c   
         do ion = 1,maxion(iel)
c
c
            rech(ion) = reco(ion)+rech(ion)
            pich(ion) = pich(ion)+pion(ion)
c
            if (reco(ion).ge.(rech(ion)/dyn)) goto 55
            if ((mod .eq. 'EQUI').and.(nom.lt.1)) then
               jdo = 2
               goto 50
            end if
            reco(ion) = reco(ion)*dyn
            pion(ion) = pich(ion)*(reco(ion)/rech(ion))
c     
            goto 50
 55         pion(ion) = pich(ion)
            reco(ion) = rech(ion)
 50         continue
         enddo
c
c     end char exchange rates
c
      endif
c
c
c    ***ALTER RECOMB. RATES FOR HE WHEN ON THE SPOT APPROX. USED
c
      if ((jspot.eq.'YES').and.(iel.eq.2)) then
c
         call spotap(de, dh, fhi, t, yh, ph, ph2, dey, dep, de2)
         reco(1) = reco(1)-((de*(1.0d0-yh))*(rec(2,2)-rec(4,2)))
         reco(2) = reco(2)-(de*(rec(3,2)-rec(5,2)))
c
      endif
c
c    ***COMPUTES NEW IONIC ABUNDANCES FOR ELEMENT : IEL
c
c     first calculation done in full
c
c      if (nom.eq.0) then
c
c
c
c      call ionab(reco, pion, pionau, ab, adndt, nde, timt)
c
c      else
c
c
c
c      if (mapz(iel).eq.2) then 
c         write(*,*) 'Helium rates (t,de,dh):',t,de,dh 
c         do i = 1,nde
c            write(*,*) rom(i),pion(i),reco(i)
c         enddo
c      endif
c
         mnde = 1
         mxde = 2
c
c     find max and min ionisation stages that matter
c
c
         do i = 1,nde
            if (pop(i,iel).ge.1.d-6) mxde = i
            if (pop(nde-i+1,iel).ge.1.d-6) mnde = nde-i+1
            if (pion(i).gt.0.d0) then
               if (reco(i)/pion(i).lt.1.d3) mxde = i
            endif
            if (pionau(i).gt.0.d0) then
               if (reco(i)/pionau(i).lt.1.d3) mxde = i
            endif
c
c     pull up the minimun, esp for Nickel and Iron with their
c     fairly flat rates
c
c            if (reco(nde-i+1).gt.0.d0) then
c               if (pion(nde-i+1)/reco(nde-i+1).lt.1.d3) mnde = nde-i+1
c            endif
         enddo
c
         mxde = min0(mxde+1,nde)
         mnde = max0(mnde-1,1)
c
c         if  (iel.eq.2) write(*,*) elem(iel),mnde,mxde
c         write(*,*) 'Min rates:',pion(mnde),reco(mnde)
c         write(*,*) 'Max rates:',pion(mxde),reco(mxde)
c         write(*,*) 'min/max pops:',pop(mnde,iel),pop(mxde,iel)
c
         mdelta = mxde-mnde+1
c
c     copy section into ab
c
         sigab = 0.d0
         do i = 1,mdelta
            a(i) = ab(i+mnde-1)
            sigab = sigab+pop(i+mnde-1,iel)
            pn(i) = pion(i+mnde-1)
            rc(i) = reco(i+mnde-1)
            pa(i) = pionau(i+mnde-1)
        enddo
         if (sigab.le.0.99d0) then 
            write (*,*) 'IOBAL VECTOR SAMPLING OPT. FAILED'
            write (*,*) 'Significant population missed:',1.d0-sigab
            write (*,*) elem(iel),'min:',mnde,'Max:',mxde
         endif
c
c     do ionab with shortened vector
c
         call ionab(rc,pn,pa,a,adnt,mdelta,timt)
c
c     copy results back into their correct positions
c
         do i = 1,mdelta
            adndt(i+mnde-1)= adnt(i)
            ab(i+mnde-1) = a(i)
         enddo
c     
         if (mnde.gt.1) then
            do i = 1,mnde-1
               adndt(i)= 0.0d0
               ab(i) = 0.0d0
            enddo
         endif
c     
         if (mxde.lt.maxion(iel)) then
            do i = mxde+1,maxion(iel)
              adndt(i)= 0.d0
               ab(i) = 0.0d0
            enddo
         endif
c     
c     
c     endif
c     
c    ***CHECK IF IT HAS TO RE-DO COMPUTATION TO REACH EQUILIBRIUM
c
      nom = nom+1
      if ((mod .eq. 'EQUI').and.(nom.lt.jdo)) goto 111
c
c    ***COPY FINAL IONIC POPULATIONS AND TIME DERIVATIVE FOR EL. IEL
c
      do 60 i = 1, nde
         dndt(i,iel) = adndt(i)
 60      pop(i,iel) = (ab(i)/(dh*zion(iel)))
c
  100 continue
c
c
c    ***UPDATES ELECTRONIC DENSITY
c
      de = feldens(dh,pop)
c
      return 
      end
