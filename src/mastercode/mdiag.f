cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c	BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c	BACK SUBSTITUTION
c
c	JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c	ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c	ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c	R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c	X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag(jl, alph, x)
c
      include 'const.inc'
c
      double precision aa,aaa,dx
      double precision ra,rab
      double precision alph(7, 6), b(7, 6), a(7, 6), x(6), r(6)
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      integer*4 jl
c
      jc = jl+1
      jca = jc+1
      do 72 il = 1, jl
      do 73 ic = 1, jca
   73 a(ic,il) = alph(ic,il)
   72 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa = 0.0d0
      do 1 il = 1, jl
      aa = aa+dabs(a(il,il))
      x(il) = 0.0d0
    1 continue
      do 2 il = 1, jl
      do 9 ic = 1, jc
      a(ic,il) = a(ic,il)/aa
      b(ic,il) = a(ic,il)
    9 continue
   52 a(jca,il) = a(jc,il)
      b(jca,il) = a(jc,il)
    2 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 4 il = 1, jl
      do 10 kl = il, jl
      rab = a(il,kl)
      ra = dabs(rab)
      if (ra.lt.epsilon) goto 10
      do 11 kc = 1, jc
      a(kc,kl) = a(kc,kl)/rab
   11 continue
   10 continue
c
      do 16 kl = il, jl
      r(kl) = a(il,kl)
   16 continue
      do 12 kl = 1, jl
      if (kl .eq. il) goto 12
      ra = dabs(r(kl))
      if (ra.lt.epsilon) goto 12
      do 13 kc = 1, jc
   13 a(kc,kl) = a(kc,kl)-a(kc,il)
   12 continue
    4 continue
c
c    ***BACK SUBSTITUTION
c
      do 17 ila = 1, jl
      il = (jl+1)-ila
      ics = il+1
      aaa = 0.0d0
      if (ics.gt.jl) goto 50
      do 18 ic = ics, jl
   18 aaa = aaa+(a(ic,il)*x(ic))
   50 dx = (a(jc,il)-aaa)/a(il,il)
      x(il) = x(il)+dx
   17 continue
   57 continue
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c	BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c	BACK SUBSTITUTION
c
c	JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c	ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c	ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c	R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c	X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag3(jl, alph, x)
c
      include 'const.inc'
c
      double precision aa,aaa,dx
      double precision ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      double precision alph(5, 4), b(5, 4), a(5, 4), x(4), r(4)
      integer*4 jl
c
      jc = jl+1
      jca = jc+1
      do 72 il = 1, jl
      do 73 ic = 1, jca
   73 a(ic,il) = alph(ic,il)
   72 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa = 0.0d0
      do 1 il = 1, jl
      aa = aa+dabs(a(il,il))
      x(il) = 0.0d0
    1 continue
      do 2 il = 1, jl
      do 9 ic = 1, jc
      a(ic,il) = a(ic,il)/aa
      b(ic,il) = a(ic,il)
    9 continue
   52 a(jca,il) = a(jc,il)
      b(jca,il) = a(jc,il)
    2 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 4 il = 1, jl
      do 10 kl = il, jl
      rab = a(il,kl)
      ra = dabs(rab)
      if (ra.lt.epsilon) goto 10
      do 11 kc = 1, jc
      a(kc,kl) = a(kc,kl)/rab
   11 continue
   10 continue
c
      do 16 kl = il, jl
      r(kl) = a(il,kl)
   16 continue
      do 12 kl = 1, jl
      if (kl .eq. il) goto 12
      ra = dabs(r(kl))
      if (ra.lt.epsilon) goto 12
      do 13 kc = 1, jc
   13 a(kc,kl) = a(kc,kl)-a(kc,il)
   12 continue
    4 continue
c
c    ***BACK SUBSTITUTION
c
      do 17 ila = 1, jl
      il = (jl+1)-ila
      ics = il+1
      aaa = 0.0d0
      if (ics.gt.jl) goto 50
      do 18 ic = ics, jl
   18 aaa = aaa+(a(ic,il)*x(ic))
   50 dx = (a(jc,il)-aaa)/a(il,il)
      x(il) = x(il)+dx
   17 continue
   57 continue
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c	BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c	BACK SUBSTITUTION
c
c	JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c	ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c	ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c	R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c	X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag6(alph, x)
c
      include 'const.inc'
c
      double precision aa,aaa,dx
      double precision ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      double precision alph(8, 7), b(8, 7), a(8, 7), x(7), r(7)
      integer*4 jl
c
      jl = 6
c
      jc = jl+1
      jca = jc+1
      do 72 il = 1, jl
      do 73 ic = 1, jca
   73 a(ic,il) = alph(ic,il)
   72 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa = 0.0d0
      do 1 il = 1, jl
      aa = aa+dabs(a(il,il))
      x(il) = 0.0d0
    1 continue
      do 2 il = 1, jl
      do 9 ic = 1, jc
      a(ic,il) = a(ic,il)/aa
      b(ic,il) = a(ic,il)
    9 continue
   52 a(jca,il) = a(jc,il)
      b(jca,il) = a(jc,il)
    2 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 4 il = 1, jl
      do 10 kl = il, jl
      rab = a(il,kl)
      ra = dabs(rab)
      if (ra.lt.epsilon) goto 10
      do 11 kc = 1, jc
      a(kc,kl) = a(kc,kl)/rab
   11 continue
   10 continue
c
      do 16 kl = il, jl
      r(kl) = a(il,kl)
   16 continue
      do 12 kl = 1, jl
      if (kl .eq. il) goto 12
      ra = dabs(r(kl))
      if (ra.lt.epsilon) goto 12
      do 13 kc = 1, jc
   13 a(kc,kl) = a(kc,kl)-a(kc,il)
   12 continue
    4 continue
c
c    ***BACK SUBSTITUTION
c
      do 17 ila = 1, jl
      il = (jl+1)-ila
      ics = il+1
      aaa = 0.0d0
      if (ics.gt.jl) goto 50
      do 18 ic = ics, jl
   18 aaa = aaa+(a(ic,il)*x(ic))
   50 dx = (a(jc,il)-aaa)/a(il,il)
      x(il) = x(il)+dx
   17 continue
   57 continue
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c	BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c	BACK SUBSTITUTION
c
c	JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c	ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c	ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c	R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c	X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
      subroutine mdiag9(alph, x)
c
      include 'const.inc'
c
      double precision aa,aaa,dx
      double precision ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      double precision alph(11, 10), b(11, 10), a(11, 10), x(10), r(10)
      integer*4 jl
c
      jl = 9
c
      jc = jl+1
      jca = jc+1
      do 72 il = 1, jl
      do 73 ic = 1, jca
   73 a(ic,il) = alph(ic,il)
   72 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa = 0.0d0
      do 1 il = 1, jl
      aa = aa+dabs(a(il,il))
      x(il) = 0.0d0
    1 continue
      do 2 il = 1, jl
      do 9 ic = 1, jc
      a(ic,il) = a(ic,il)/aa
      b(ic,il) = a(ic,il)
    9 continue
   52 a(jca,il) = a(jc,il)
      b(jca,il) = a(jc,il)
    2 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 4 il = 1, jl
      do 10 kl = il, jl
      rab = a(il,kl)
      ra = dabs(rab)
      if (ra.lt.epsilon) goto 10
      do 11 kc = 1, jc
      a(kc,kl) = a(kc,kl)/rab
   11 continue
   10 continue
c
      do 16 kl = il, jl
      r(kl) = a(il,kl)
   16 continue
      do 12 kl = 1, jl
      if (kl .eq. il) goto 12
      ra = dabs(r(kl))
      if (ra.lt.epsilon) goto 12
      do 13 kc = 1, jc
   13 a(kc,kl) = a(kc,kl)-a(kc,il)
   12 continue
    4 continue
c
c    ***BACK SUBSTITUTION
c
      do 17 ila = 1, jl
      il = (jl+1)-ila
      ics = il+1
      aaa = 0.0d0
      if (ics.gt.jl) goto 50
      do 18 ic = ics, jl
   18 aaa = aaa+(a(ic,il)*x(ic))
   50 dx = (a(jc,il)-aaa)/a(il,il)
      x(il) = x(il)+dx
   17 continue
   57 continue
c
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c*******COMPUTES SOLUTION TO A SYSTEM OF LINEAR EQUATIONS
c	BY DIAGONALISATION TO TRIANGULAR FORM AND THEN
c	BACK SUBSTITUTION
c
c	JL = NUMBER OF VARIABLES ; A AND B : WORKING MATRIX
c	ALL MATRICES , MINIMUM DIMENSION : (JL+2) BY (JL+1)
c	ALPH = INPUT MATRIX (REMAIN UNCHANGED)
c	R = WORKSPACE VECTOR (LENGTH JLMAX+1)
c	X = SOLUTION VECTOR (LENGTH JLMAX+1)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mdiag16( alph, x)
c
      include 'const.inc'
c
      double precision aa,aaa,dx
      double precision ra,rab
      integer*4 ic,ics,il,ila,jc
      integer*4 jca,kc,kl
      double precision alph(18,17),b(18,17),a(18,17),x(17),r(17)
      integer*4 jl
c
      jl = 16
      jc = jl+1
      jca = jc+1
      do 72 il = 1, jl
      do 73 ic = 1, jca
   73 a(ic,il) = alph(ic,il)
   72 continue
c
c    ***RENORMALISATION OF MATRIX
c
      aa = 0.0d0
      do 1 il = 1, jl
      aa = aa+dabs(a(il,il))
      x(il) = 0.0d0
    1 continue
      do 2 il = 1, jl
      do 9 ic = 1, jc
      a(ic,il) = a(ic,il)/aa
      b(ic,il) = a(ic,il)
    9 continue
   52 a(jca,il) = a(jc,il)
      b(jca,il) = a(jc,il)
    2 continue
c
c    ***GAUSSIAN REDUCTION TO TRIANGULAR FORM
c
      do 4 il = 1, jl
      do 10 kl = il, jl
      rab = a(il,kl)
      ra = dabs(rab)
      if (ra.lt.epsilon) goto 10
      do 11 kc = 1, jc
      a(kc,kl) = a(kc,kl)/rab
   11 continue
   10 continue
c
      do 16 kl = il, jl
      r(kl) = a(il,kl)
   16 continue
      do 12 kl = 1, jl
      if (kl .eq. il) goto 12
      ra = dabs(r(kl))
      if (ra.lt.epsilon) goto 12
      do 13 kc = 1, jc
   13 a(kc,kl) = a(kc,kl)-a(kc,il)
   12 continue
    4 continue
c
c    ***BACK SUBSTITUTION
c
      do 17 ila = 1, jl
      il = (jl+1)-ila
      ics = il+1
      aaa = 0.0d0
      if (ics.gt.jl) goto 50
      do 18 ic = ics, jl
   18 aaa = aaa+(a(ic,il)*x(ic))
   50 dx = (a(jc,il)-aaa)/a(il,il)
      x(il) = x(il)+dx
   17 continue
   57 continue
c
      return 
      end
