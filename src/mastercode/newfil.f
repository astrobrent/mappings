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
c     This utility subroutine will return a filename that is
c     unique in the current directory.
c     It accepts a prefix and suffix.  Prefixes can be up to 8 chars
c     and suffixes can be 4.  The final filemane will be a character*20
c     string.  Multiple of four are used for SPARC and RISC optimisation.
c
c     RSS 10/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine newfile(pref,p,suff,s,filena)
c
      include 'const.inc'
c
      character pref*8,suff*4,filena*20,s1*8,s2*4
      integer*4 p,s,i,j
      logical iexi
c
c
c
 1    format(i4.4)
c
c
c
      filena = ' '
      j = p + 4 + 1 + s
      filena(1:j) = ' '
c
      i = 0
 100  i  = i+1
      s2 = ' '
      write(s2, 1) i
      j = p+1
      filena(j:j+4) = s2
c
      s1 = ' '
      if ((p.gt.0).and.(p.lt.9)) s1 = pref(1:p)
      j = p+1
      filena(1:p) = s1
c
      s2 = ' '
      if ((s.gt.0).and.(s.lt.5)) s2 = suff(1:s)
      j = p+6
      filena(j:j+s) = s2
      j = j-1
      filena(j:j) = '.'
c
      inquire(file=filena, exist=iexi) 
c
      if (iexi) goto 100
c
      return
c
      end
