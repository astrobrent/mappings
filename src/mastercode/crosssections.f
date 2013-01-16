ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine crosssections(inl, tauso, sigmat, dustsigmat)
c
c tausso is the integrated optical depth back to the source
c sigmat, and dustsigmat are the local crossections per H atom
c for the current zone
c
c
      include 'cblocks.inc'
c
      double precision tauso, sigmat, dustsigmat
      integer*4 inl

      double precision aph, crosec, den
      double precision acrs, eph, at, bet, se
      integer*4 i, j, ie, dtype, m
c
c     internal function
c
      acrs(eph,at,bet,se) = (at*(bet+((1.0d0-bet)/eph)))*(eph**(-se))
                  
      tauso      = 0.d0
      sigmat     = 0.d0
      dustsigmat = 0.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     
c     ***DERIVES ABSORPTION CROSS-SECTION FOR BIN#INL
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         
      den = 0.5d0*(ephot(inl)+ephot(inl+1))
      do  i = 1, ionum
         ie = atpho(i)
         j = ionpho(i)
         aph = zion(ie)*pop(j,ie)
         eph = den/ipotpho(i)
         if (eph.ge.1.0d0) then
          crosec = acrs(eph,sigpho(i),betpho(i),spho(i))
          tauso = tauso+(popint(j,ie)*crosec)
          sigmat = sigmat+(aph*crosec)
         endif
      enddo
c
c   Dust Cross-section
c
      if (grainmode) then 
c
c   pahs
        if (pahmode) then
           do m=1,pahi
             if (pahion(m).gt.0) then
               crosec = pahiext(inl)*pahfrac
             else
               crosec = pahnext(inl)*pahfrac
             endif
             tauso = tauso+pahint(m)*crosec
             sigmat = sigmat + pahZ(m)*crosec
           enddo
        endif 
c
c   dust grains
        if (.NOT.ClinPAH) then !C grains not linked to PAH
          do dtype=1,numtypes
            tauso = tauso+dustint*dcrosec(inl,dtype)
            sigmat = sigmat+dcrosec(inl,dtype)
            dustsigmat=dustsigmat+dcrosec(inl,dtype)
          enddo
        else !C linked to PAH
          if (pahmode) then !and PAHs exist
            do m=1,pahi
               tauso = tauso+pahint(m)*dcrosec(inl,1)
            enddo
            sigmat = sigmat+dcrosec(inl,1)
            dustsigmat=dustsigmat+dcrosec(inl,1)
         endif
            do dtype=2,numtypes
             tauso = tauso+dustint*dcrosec(inl,dtype)
             sigmat = sigmat+dcrosec(inl,dtype)
             dustsigmat=dustsigmat+dcrosec(inl,dtype)
            enddo
        endif	  
      endif
c
      return
  
      end
