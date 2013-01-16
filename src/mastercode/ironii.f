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
c*******COMPUTES  16-LEVEL ATOM COOLING FOR THE IRON II ION
c
c     RSS 12/91
c     RSS 94
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ironii(t, de, dh)
c
      include 'cblocks.inc'
c
c
c
c           Variables
c
      double precision t, de, dh
      double precision w(16),e(16, 16),a(16,16),om(16,16) 
      double precision rd(16,16), re(16,16)
      double precision alph(18,17),x(17)
      double precision da
      double precision f,po,zi,aa,aion,bion
c
      integer*4 j,jc,jl,il,ion,atom,nl
      integer*4 jll,m
c
      feloss = 0.0d0
      f = dsqrt(1.0d0/t)
c
      nl = 16
      ion = 2
      atom = zmap(26)
c
c     first check if we have iron at all
c
      if ((atom.ne.0).and.(nfetrans.ge.1) ) then

c
         do j = 1,nfetrans
            febri(j) = 0.d0
         enddo
c
         po = pop(ion,atom)
         zi = zion(atom)
c
      if (po*zi.ge.1.d-16) then
c
c
         do j = 1, nl
            w(j) = wfe(j)
            do m = 1, nl
               e(m,j) = efe(m,j)
               a(m,j) = afe(m,j)
               om(m,j) = omfe(m,j)
            enddo
         enddo

c
c    ***SET UP COLL. EXCIT. AND DEEXCIT. MATRIX
c
         do jl = 1, nl
            do jc = 1, nl
            aa = e(jc,jl)/(rkb*t)
            if (aa.ge.650.d0) aa = 650.d0
            rd(jc,jl) = 0.0d0
            if (jc.gt.jl) rd(jc,jl) = ((rka*f)*om(jc,jl))/w(jc)
            re(jc,jl) = 0.0d0
            if (jc.lt.jl) re(jc,jl) = (((rka*f)*om(jc,jl))
     &                                *dexp(-aa))/w(jc)
            enddo
         enddo
c
c    ***SET UP MATRIX ELEMENTS
c
         do 40 jl = 1, nl
            do 41 jc = 1, nl
               if (jl .eq. 1) goto 113
               if (jl .eq. jc) goto 114
               alph(jc,jl) = (de*(rd(jc,jl)+re(jc,jl)))+a(jc,jl)
               goto 41
 113           alph(jc,jl) = 1.0d0
               goto 41
 114           alph(jc,jl) = 0.0d0
               do 42 jll = 1, nl
                  if (jc .eq. jll) goto 42
                  alph(jc,jl) = alph(jc,jl)-(de*(rd(jc,jll)+re(jc,jll))
     &                 +a(jc,jll))
 42            continue
 41         continue
 40      continue
c
         do  jl = 1, nl
            alph(nl+1,jl) = 0.0d0
         enddo
         alph(nl+1,1) = 1.0d0
c
c    ***SOLVE FOR LEVEL POPULATIONS  X
c
c    ***CHECK ON NORMALISATION
c
         call mdiag16(alph, x)
         da = 0.0d0
         do 200 jc = 1, nl
            if (x(jc).lt.0.0d0) x(jc) = 0.0d0
            if (x(jc).gt.0.0d0) da = da+x(jc)
 200     continue
         do 205 jc = 1, nl
 205        x(jc) = x(jc)/da
c     
c    ***GET ION COOLING
c
         aion = zi*po*dh
         il = 0
         do jl = 1, nl-1
            do jc = jl+1, nl
               bion = x(jc)*e(jc,jl)*aion*a(jc,jl)
               if (afe(jc,jl).ne.0.d0) then
                  il = il+1
                  febri(il) = bion/(4.0d0*pi)
                  feloss = feloss+bion
              endif
           enddo
        enddo
c
c     end abundant ion
c
      endif
c
c     end fe
c
      endif
c
      return 
      end
