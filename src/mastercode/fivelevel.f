cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES  5-LEVEL ATOM COOLING FOR THE 
c
c	RETURNS : FORBIDDEN LINE COOLING RATE FLOSS (ERG.CM-3.S-1)
c	AND BRIGHTNESS OF EACH LINE FBRI(10,16) (ERG.CM-3.S-1.STER.-1)
c	AVAILABLE IN COMMON BLOCK /FLINE/
c	CALL SUBROUTINE MDIAG
c
c
c
      subroutine fivelevel(t, de, dh)
c
      include 'cblocks.inc'
c
      double precision t, de, dh
      double precision w(5), e(5, 5), a(5, 5),om(5,5),re(5,5) 
      double precision rd(5, 5), alph(7, 6), x(6)
      double precision f,t4,po,zi,filoss
      double precision aa,aion,bion,da
c
      integer*4 il,jc,jl,jll
      integer*4 i,j,m
c
      floss = 0.0d0
      f = dsqrt(1.0d0/t)
      t4 = t/(1.d4)
c
      do 21 i = 1, nfions
c
c     only calculate abundant ions
c
         do j = 1,10
            fbri(j,i) = 0.d0
         enddo
c
         po = pop(fion(i),fatom(i))
         zi = zion(fatom(i))
c
c
      if (po*zi.ge.1.d-16) then
c
c Now have t dependence for all 5 level atoms, may be 
c 0.000 for cases with no tdep.
c
         do j = 1, 5
            w(j) = wi(j,i)
            do m = 1, 5
               e(m,j) = ei(m,j,i)
               a(m,j) = ai(m,j,i)
               om(m,j) = omi(m,j,i)*(t4**tdep(m,j,i))
            enddo
         enddo
c
c    ***SET UP COLL. EXCIT. AND DEEXCIT. MATRIX
c
         do jl = 1, 5
            do jc = 1, 5
               aa = e(jc,jl)/(rkb*t)
               if (aa.ge.650.d0) aa = 650.d0
               rd(jc,jl) = 0.0d0
               if (jc.gt.jl) rd(jc,jl) = ((rka*f)*om(jc,jl))/w(jc)
               re(jc,jl) = 0.0d0
               if (jc.lt.jl) re(jc,jl) = (((rka*f)*om(jc,jl))
     &                                   *dexp(- aa))/w(jc)
c
            enddo
         enddo
c
c    ***SET UP MATRIX ELEMENTS
c
         do 40 jl = 1, 5
            do 41 jc = 1, 5
               if (jl .eq. 1) goto 113
               if (jl .eq. jc) goto 114
               alph(jc,jl) = (de*(rd(jc,jl)+re(jc,jl)))+a(jc,jl)
               goto 41
 113           alph(jc,jl) = 1.0d0
               goto 41
 114           alph(jc,jl) = 0.0d0
               do 42 jll = 1, 5
                  if (jc .eq. jll) goto 42
                  alph(jc,jl) = alph(jc,jl)-(de*(rd(jc,jll)+re(jc,jll))
     &                 +a(jc,jll))
 42            continue
 41         continue
 40      continue

         do  jl = 1, 5
            alph(6,jl) = 0.0d0
         enddo
         alph(6,1) = 1.0d0
c
c    ***SOLVE FOR LEVEL POPULATIONS  X
c
         jl = 5
c
c    ***CHECK ON NORMALISATION
c
         call mdiag(jl, alph, x)
c
c         if ((mapz(fatom(i)).eq.6).and.(fion(i).eq.3)) then
c            write(*,'(5(1pg12.5,1x))') (x(m),m=1,5)
c         endif
c
         da = 0.0d0
c
         do jc = 1, 5
            if (x(jc).lt.0.0d0) x(jc) = 0.0d0
            if (x(jc).gt.0.0d0) da = da+x(jc)
         enddo
c
         do jc = 1, 5
            x(jc) = x(jc)/da
         enddo
c     
c    ***GET ION COOLING
c
         aion = zi*po*dh
c
         filoss = 0.d0
c
         il = 1
         do  jl = 1, 4
            do jc = jl+1, 5
               bion = x(jc)*e(jc,jl)*aion*a(jc,jl)
               fbri(il,i) = bion/(4.d0*pi)
               il = il+1
               floss = floss+bion
               filoss = filoss+bion
            enddo
         enddo
c
c         if (mapz(fatom(i)).eq.6) then
c            write(*,*) elem(fatom(i)), rom(fion(i)), po
c            write(*,'(1pg14.7)') filoss
c         endif
c
c     end abundant ion
c
      endif
c
c     next ion..
c
   21 continue
c
      return 
      end
