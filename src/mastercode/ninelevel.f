cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES  9-LEVEL ATOM COOLING FOR THE 
c	FOLLOWING IONS : Fe VII
c
      subroutine ninelevel(t, de, dh)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision t, de, dh
c
      double precision w(9), e(9, 9), a(9, 9),om(9,9),re(9,9) 
      double precision rd(9, 9), alph(11, 10), x(10)
      double precision f,t4,po,zi
      double precision aa,aion,bion,da
c
      integer*4 il,jc,jl,jll
      integer*4 i,j,m,nl
c
      nl = 9
c
      f9loss = 0.0d0
      f = dsqrt(1.0d0/t)
      t4 = t/(1.d4)
c      
      do i = 1, n9ions
c
c     only calculate abundant ions
c
         do j = 1,n9trans
            f9bri(j,i) = 0.d0
         enddo
c
         po = pop(f9ion(i),f9atom(i))
         zi = zion(f9atom(i))
c
      if (po*zi.ge.1.d-16) then
c
c Now have t dependence for all 9 level atoms, may be 
c 0.000 for cases with no tdep.
c
         do j = 1, nl
            w(j) = wi9(j,i)
            do m = 1, nl
               e(m,j) = ei9(m,j,i)
               a(m,j) = ai9(m,j,i)
               om(m,j) = omi9(m,j,i)*(t4**tdep9(m,j,i))
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
     &                                   *dexp(- aa))/w(jc)
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

         do  jl = 1, nl
            alph(nl+1,jl) = 0.0d0
         enddo
         alph(nl+1,1) = 1.0d0
c
c
c    ***SOLVE FOR LEVEL POPULATIONS  X
c
c    ***CHECK ON NORMALISATION
c
         call mdiag9(alph, x)
         
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
c
         aion = zi*po*dh
c
         il = 1
         do jl = 1, nl-1
            do jc = jl+1, nl
               bion = x(jc)*e(jc,jl)*aion*a(jc,jl)
               f9bri(il,i) = bion/(4.d0*pi)
               il = il+1
               f9loss = f9loss+bion
            enddo
         enddo
c     
c     end abundant ion
c
      endif
c
c     next ion..
c
      enddo
c
      return 
      end
