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
c*******TO EXTRACT A ROOT OF A POLYNOMIAL
c	USING BERNOUILLI'S METHOD
c
c	A IS A VECTOR OF COEFF. IN ORDER OF
c	DECREASING POWER ; TOTAL*OF COEFF.=NCOEF (<9)
c	ROOT=SOLUTION ROOT
c	FEED A GUESS FOR ROOT FIRST
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine roots(a, ncoef, root, rmod)
c
      include 'const.inc'
c
      double precision root, unew,anew,roota,fract
      double precision a(10), u(10)
      integer*4 ncoef,j,iter,m,mcoef
      character rmod*4
c
c      write(*,*) 'ROOTS: init guess:',root
c
      u(1) = 1.d0
      u(2) = 1.d0/root
      do 2 j = 3, ncoef
    2 u(j) = 0.d0
      iter = 0
   51 unew = 0.d0
      iter = iter+1
      mcoef = ncoef-1
      do 3 j = 1, mcoef
    3 unew = unew-(a(j+1)*u(j))
      anew = dabs(unew)
      if (anew.lt.1.0d20) goto 50
c
c    ***RENORMALISATION OF COEFF. IF TOO LARGE
c
      do 4 j = 1, ncoef
    4 u(j) = u(j)/unew
      unew = 1.d0
   50 continue
c
      do 5 j = 2, ncoef
      m = (ncoef-j)+1
    5 u(m+1) = u(m)
c
c    ***EXTRACTION OF ROOT ; TEST FOR CONVERGENCE
c
      u(1) = unew
      roota = root
      root = u(1)/u(2)
      fract = dabs((root-roota)/root)
      if (iter.gt.150) return 
      if (fract.gt.1.0d-7) goto 51
c
c      write(*,*) 'ROOTS:',ncoef,'coeffs'
c      do i = 1,ncoef
c         write(*,100) (a(j),j=1,ncoef)
c         write(*,100) (u(j),j=1,ncoef)
c         write(*,200) iter,root,roota
c 100     format(10(1pg14.7))
c 200     format(I3.3,2(1pg14.7))
c      enddo
      return 
      end
