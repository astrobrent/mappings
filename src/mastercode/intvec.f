cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*******INTEGRATES ARRAY BUFPHO THAT CONTAINS THE LOCAL MEAN INTENS.
c   JNU TO DERIVE THE AVAILABLE NUMBER OF PHOTONS TO IONISE H,HE
c   (HENCE INTEGRATES 4PI*JNU/HNU)
c   NB.  BUFPHO MUST BE DEXPRESSED IN  ERGS.CM-2.SEC-1.HZ-1.STER-1 (1/4pi)
c   CALL SUBR. INTERPOL
c
c
c
      subroutine intvec(bufpho, qahi, qahei, qaheii, qatot)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision bufpho(mxinfph)
      double precision qahi, qahei, qaheii, qatot
      double precision piplklo
      double precision ain,alh,alo,alpha,bet,cflo,e1,e2
      double precision etr,rnutr,rr1,rr2,sa
      double precision phq(2, 2)
c
      integer*4 i, j, inl
c
c           Functions
c
      double precision dfuint
c
c     ln(4pi/h)
c
      piplklo = 62.8098090594892475d0

      do  i = 1, 2
      do  j = 1, 2
          phq(j,i) = 0.d0
      enddo
      enddo
c
      do 2500 inl = 1, infph-1
      if (bufpho(inl).le.epsilon) goto 2500
      e1 = ephot(inl)
      e2 = ephot(inl+1)
c
       call interpol(inl, cflo, alpha, bufpho)
c
       alh = alpha-1.0d0
c
       do 2400 i = 1, 2
        do 2300 j = 1, maxion(i)-1
         etr = epot(j,i)/ev
         rr2 = e2/etr
         if (rr2.le.1.0001d0) goto 2300
         rr1 = e1/etr
         if (rr1.lt.1.0d0) rr1 = 1.0d0
         bet = 1.0d0
         sa = alh
         rnutr = etr*evplk
         alo = (piplklo+cflo)+((1.0d0+alh)*dlog(rnutr))
         ain = dfuint(bet,sa,rr1,rr2,alo)
         phq(j,i) = phq(j,i)+ain
 2300   continue
 2400  continue
 2500 continue

      qaheii = phq(2,2)
      qahei = dmax1(0.d0,phq(1,2)-phq(2,2))
      qahi = dmax1(0.d0,phq(1,1)-phq(1,2))
c
      qatot = phq(1,1)
c
      return 
      end
