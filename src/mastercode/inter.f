cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES INTERCOMBINATION LINE COOLING
c       NB.  SUBR. HYDRO SHOULD BE CALLED PREVIOUSLY
c
c	RETURNS DATA IN COMMON BLOCK /CLINE/
c	FSLOS (ERG.CM-3.S-1)
c
c     NEW: uses detailed fit to CII 158 mu
c
      subroutine inter(t, de, dh)
c
      include 'cblocks.inc'
c
c
c
      double precision t, de, dh
      double precision rp,q21,q12,aa,pn2,fbr
      double precision lx,lt,lr,x1,omg
      double precision f, zi, po,ba,t2,hrat
      integer*4 jel,jjj,m,j,jj
c
c
c
c
      fslos = 0.0d0
c
c    ***COMPUTES RATES Q12,Q21 , HENCE POPULATION OF LEVEL 2
c
      f = dsqrt(1.0d0/(t+epsilon))
      t2 = t/100.d0
      do 10 m = 1, mlines
         fbr = 0.0d0
         jel = ielfs(m)
         jjj = ionfs(m)
         fsbri(m) = 0.d0
         zi = zion(jel)
         po = pop(jjj,jel)
c
         if (zi*po.ge.pzlimit) then

            aa = e12fs(m)/(rkb*t)
            ba = dexp(-aa)
            omg = omfs(m)
c
c     use fit for NeII 12.8 mu
c
            if(  (jel.eq.zmap(10))
     &           .and.(jjj.eq.2)
     &           .and.(e12fs(m).eq.1.5501E-13)) then
                 omg = 0.293+(0.8e-6*t)
            endif
c
c     use fit for CII 158 mu
c
            if(  (jel.eq.zmap(6))
     &           .and.(jjj.eq.2)
     &           .and.(e12fs(m).eq.1.2598E-14)) then
c
c     very important, assumes CII 158 is first spline
cc
c     CII spline in linear te
c
               lt = t
               x1 = omgspl(1,1,1)
               jj = 1
               do j = 1,omgnsp(1)
                  if (lt.ge.omgspl(1,j,1)) then
                     x1 = omgspl(1,j,1)
                     jj = j
                  endif
               enddo
c
               lx = lt - omgspl(1,jj,1)
               lr = lx*omgspl(1,jj,5)
               do j = 4,3,-1
                  lr = lx*(omgspl(1,jj,j) + lr)
               enddo
               lr = lr+omgspl(1,jj,2)
c
               if (lr.lt.0.d0) lr = 0.d0
               if (lr.gt.3.d0) lr = 3.d0
c
               omg = lr
c
            endif
c
c
c     use fit for S IV 10.5 mu
c
            if(  (jel.eq.zmap(16))
     &           .and.(jjj.eq.4)
     &           .and.(e12fs(m).eq.1.8876E-13)) then
c
c     very important, assumes S IV is second spline
c
c     SIV spline in log Te
c
               lt = dlog10(t)
               x1 = omgspl(2,1,1)
               jj = 1
               do j = 1,omgnsp(2)
                  if (lt.ge.omgspl(2,j,1)) then
                     x1 = omgspl(2,j,1)
                     jj = j
                  endif
               enddo
c
               lx = lt - omgspl(2,jj,1)
               lr = lx*omgspl(2,jj,5)
               do j = 4,3,-1
                  lr = lx*(omgspl(2,jj,j) + lr)
               enddo
               lr = lr+omgspl(2,jj,2)
c
               if (lr.lt.0.d0) lr = 0.d0
c               write(*,*) 'Siv :', lr
c
               omg = lr
c
            endif
c
            q12 = (rka*f)*omg*ba/w1fs(m)
            q21 = (rka*f)*omg/w2fs(m)
            if(((jel.eq.zmap(6)).or.(jel.eq.zmap(14)))
     &           .and.(jjj.eq.2)
     &           .and.((e12fs(m).eq.1.2598E-14)
     &           .or.(e12fs(m).eq.5.7076E-14))) then
c     
c     add hydrogen rate to CII and SiII
c     
               hrat = ((dh/de)*pop(1,1)*8.d-10*(t2**0.07))
               q21 = q21 + hrat
               q12 = q12 + (hrat*ba*w2fs(m)/w1fs(m)) 
c     
            endif
c     
            if((jel.eq.zmap(10))
     &           .and.(jjj.eq.2)
     &           .and.(e12fs(m).eq.1.5501E-13)) then
c     
c     add hydrogen rate to Ne II (different t depend)
c     
               hrat = (dh/de)*pop(1,1)*1.3d-9
               q21 = q21 + hrat
               q12 = q12 + (hrat*ba*w2fs(m)/w1fs(m)) 
c     
            endif
            rp = (de*q12)/(a21fs(m)+(de*q21))
            pn2 = dh*zi*po*(rp/(1.0d0+rp))
            fbr = (pn2*a21fs(m))*e12fs(m)
            fslos = fslos+fbr
            fsbri(m) = fbr/(4.d0*pi)
            goto 10
 100        continue
c 
c	Now handled in heibri as part of helio
c    
c            fsbri(m) = hein(1)
c
         endif
c
 10   continue
c
      return 
      end
