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
c
c     Subroutine freefree calculates the freefree continuum contribution
c     to the local diffuse field..
c     adds the emission to cnphot (in photons/cm3/s/Hz/sr)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine freefree(t,de,dh)
c
      include 'cblocks.inc'
c
      double precision cebin,bev,et
c
      double precision t,de,dh,zpop,abde,ffconst
      double precision u, gffm, g2, zn2,t12,rkt
      double precision fgfflin,fgfflog,phots
      integer*4 atom,ion,i,j,m,inl
c
c     external functions
c
      integer*4 GetIndexUpwards
c
c
      t12    = 1.d0/dsqrt(t)
      rkt    = rkb*t
c
      do inl = 1,infph
         ffph(inl) = 0.d0
      enddo
c
c     loop through all ions present
c
      do atom = 1 , atypes
         do ion = 1,maxion(atom)-1
c
c     only consider abundant ions
c
         zpop = zion(atom)*pop(ion+1,atom)
         abde = zpop*de*dh
c         
         if (zpop.ge.pzlimit) then
c
c     calc gamma^2 : zn0 from data file: CONTDAT
c
         zn2     = zn0(ion,atom)*zn0(ion,atom)
         g2      = dlog10(ryde*zn2/(rkt))
         ffconst = ffk*abde*zn2*t12

         if (g2.lt.-4.d0)then
c
c        lower limit of -4 is a v. good assymptote (unlike +4)
c
              g2 = -4.0d0               
c
         endif
c
         do inl = 1,infph-1
c
c     check energy compared to kt
c
         cebin = 0.5d0*(ephot(inl)+ephot(inl+1))
         bev   = cebin*ev
c
c     get scaled energy of bin wrt te
c
         et  = (bev)/(rkt)
         u   = dlog10(et)
c
c
c prevent over/underflow for exp -et
c
         if (dabs(et).lt.loghuge) then
c
c     get free free gaunt factor
c
c
c     Lookup table, extrapolate if necessary in u - and g2 if it is above 4
c
            if (u.lt.-4.d0) then
              
                  if (g2.le.4.d0) then
c
c    linear extrapolation in log lin space at g2, common and accurate
c
                      m = 2
c
                      i = GetIndexUpwards(g2, ngffg2, gffg2)
                      j = 1
                      gffm = fgfflin(m,g2,u,i,j)
c
c                      write(*,'("Ext2:",3(1pg14.7))'),g2,u,gffm
c
                  else
c
c    quadratic extrapolation in log lin space for u and g2 (scary!) should be rare
c
                      m = 3
c
                      i = ngffg2
                      j = 1
                      gffm = fgfflin(m,g2,u,i,j)
c
c                     SANITY!
c
                      if (gffm.lt.0.d0) gffm = 0.d0
c
c                      write(*,'("WARNING ff Ext3:",3(1pg14.7))'),g2,u,gffm
c
                  endif
c
            else if (u.gt.4.d0) then
c
c    Cannot happen unless contcont is changed.  u is limited to log(loghuge) (ie about 2.4)
c    since answer is mulitplied by exp(-u) anyway.
c
                  if (g2.le.4.d0) then
c
c    powerlaw extrapolation in u log log space at g2
c
c    potentially inaccurate at very low T
c
                      m = 2
c
                      i = GetIndexUpwards(g2, ngffg2, gffg2)
                      j = ngffu
                      gffm = fgfflog(m,g2,u,i,j)
c
c                      write(*,'("Ext4:",3(1pg14.7))'),g2,u,gffm
c
                  else
c
c    powerlaw extrapolation in log log space for u and g2
c
                      m = 2
c
                      i = ngffg2
                      j = ngffu
                      gffm = fgfflog(m,g2,u,i,j)
c
c                     SANITY!
c
                      if (gffm.lt.0.d0) gffm = 0.d0
c
c                      write(*,'("WARNING ff Ext5:",3(1pg14.7))'),g2,u,gffm
c
                  endif
c
            else
c              
c   Normal evaluations here
c
                  if (g2.le.4.d0) then
c
                    i = GetIndexUpwards(g2, ngffg2, gffg2)
                    j = GetIndexUpwards(u, ngffu, gffu)
c
c  Cubic interpolation in log log space
c    
                    m = 4
c
                    gffm = fgfflog(m,g2,u,i,j)
c               
                  else
c
c    quad extrapolation in u log log space for g2 - OK
c
                      m = 3
c
                      i = ngffg2
                      j = GetIndexUpwards(u, ngffu, gffu)
                      gffm = fgfflog(m,g2,u,i,j)
c
c                     SANITY!
c
                      if (gffm.lt.1.d0) gffm = 1.d0
c
c                      write(*,'("WARNING ff Ext6:",3(1pg14.7))'),g2,u,gffm
c
                  endif
c
            endif
c               
c    not collecting Mewe style gaunt factors normally - place after phots 
c    line when using
c
c            ffgau(inl)  = ffgau(inl)  + bev*phots/(ffk*de*de*t12*dexp(-et))
c
            phots       = ffconst*gffm*dexp(-et)/bev
            ffph(inl)   = ffph(inl) + phots
            cnphot(inl) = cnphot(inl) + phots
c
c     end et limit
c
         endif
c
c     end inl/infph loop
c  
         enddo
c
c     end abundant ions if loop
c
         endif
c
c     end ion and atom loops
c
         enddo
      enddo
c
c
      return
c
      end
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
c*******TO OBTAIN FREE-FREE COOLING RATE ; FFLOS (ERG.CM-3.S-1)
c
c     **REFERENCES : ALLEN,C.W. 1973 AP.QTIES. ED.3,103
c     **REFERENCES : Sutherland 1997
c     **REFERENCES : Karzas & Latter 1961
c
c
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine frefre(t, de, dh)
c
      include 'cblocks.inc'
c
c old
c
      double precision t, de, dh 
      double precision chsq,u
c
c new
c
      double precision zpop,abde,abdftz2
      double precision zn2,t12,rkt,gffint
      double precision lx,lt,lr,x1,g2
      integer*4 j,jj,icount
      integer*4 ion,atom
c
c     double precision g2bar,gion
c
c     internal function: gef
c
      double precision gef,densnum,favcha
c
      gef(u)=1.0d0+0.44d0*dexp(-((dabs(dlog10(u)-5.5d0)/2.4d0)**1.5d0))
c
      fflos = 0.d0
      chsq = favcha(pop,2)
      chsq = chsq*chsq
c
      fflos =1.435d-27*gef(t/chsq)*dsqrt(t)*de*densnum(dh)*chsq
c
      t12    = dsqrt(t)
      rkt    = rkb*t
      gffint = 0.d0
c
c     loop through all ions present
c
c      g2bar = 0.d0
c      icount = 0
c      gion = 0.d0
c
      do atom = 1 , atypes
         do ion = 1 , maxion(atom)-1
c
c     only consider abundant ions
c
         zpop = zion(atom)*pop(ion+1,atom)
         abde = zpop*de*dh
c         
         if (zpop.ge.pzlimit) then
c
c     calc gamma^2 : zn0 from data file: CONTDAT
c
         zn2     = zn0(ion,atom)*zn0(ion,atom)
         g2      = epot(ion, atom)/(rkt)
c
c         if (atom.gt.2) then
c         g2bar = g2bar + g2*zpop
c         gion  = gion  + zpop
c         endif
c
         icount = icount+1
c        
         abdftz2 = abde*zn2*t12
c     
c     Evaluate Spline at log(g2) and convert back to normal units
c     
c     
c     Get spline nodes, g2 will be different for each ion
c     
         lt = dlog10(g2)
         x1 = gffintspl(1,1)
         jj = 1
         do j = 1,ngffg2
            if (lt.ge.gffintspl(j,1)) then
               x1 = gffintspl(j,1)
               jj = j
            endif
         enddo
c     
c     interp and evaluate in RPN
c     
         lx = lt - gffintspl(jj,1)
         lr = lx*gffintspl(jj,5)
         do j = 4,3,-1
            lr = lx*(gffintspl(jj,j) + lr)
         enddo
         lr = lr+gffintspl(jj,2)
c
         gffint = gffint + lr*abdftz2
c
         endif
c
         enddo
      enddo
c
c      write(*,*) g2bar/gion
c
c      write(*,'(3(1pg14.7))') fflos, Fk*gffint, Fk*gffint/fflos
c
      fflos = Fk*gffint
c
      return 
      end
