cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FINDS THE EQUILIBRIUM TEMPERATURE AND THE
c     CORRESPONDING IONISING STATE OF THE GAS
c     AT EQUILIBRIUM IONISATION OF AFTER TIME STEP :STEP
c     TEI : INITIAL GUESS FOR TEMPERATURE
c     TEF,EDENS: FINAL TEMPERATURE AND ELECTRONIC DENSITY
c
c     NMOD = 'EQUI'  :  EQUILIBRIUM IONISATION
c     NMOD = 'TIM'   :  INITIAL IONIC POP. EVOLVED BY TSTEP SEC.
c
c     CALL SUBR. COOL,EQUION
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine teequi(tei, tef, edens, hdens, tstep, nmod)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision a,b,c1,c2,cc,csub,delexi,delmin
      double precision dl1,dlodef,dlpr,dt12,dvv,fadl,ftu,ftu0
      double precision t2,tdef,timin,tl1,tl2,tl2a,tm1,tmin,tsub
      double precision tti,u,varf,xhyf
      double precision tei, teinit,tef, edens, hdens, tstep
      double precision popul(mxion, mxelem)
c
      integer*4 iconsis,isis,jok,luop,luty,n,nf
c
      character nmod*4, ntest*4
c
      double precision arctanh
c
      arctanh(u) = 0.5d0*dlog((1.d0+u)/(1.d0-u))
c
      luty = 6
      luop = 23
c
      if ((nmod .ne. 'EQUI').and.(nmod .ne. 'TIM')) then
      write(luty, 212) nmod
  212 format(/' MODE IMPROPERLY SET FOR SUBR. TEEQUI  : ',a4)
      stop 
      end if
c
c    Fudge which sets everything constant once T =30 K
c
      if (tei.eq.30.d0) then
         tef=tei
         goto 220   !goes to end - 
      endif
      
c
c
      ntest = 'N'
      if (tef.lt.0.0d0) ntest = 'Y'
      delmin = 4.d-4
      nf = 3
c      tmin = 100.d0
c      timin = 100.d0
      tmin = 10.d0
      timin = 30.0d0
      ftu0 = 0.0d0
      isis = 0
      if (nmod .eq. 'TIM') then
         call copinto(pop, popul)
      end if
c
c    ***GUESS AN IMPROVED VALUE FOR THE TEMPERATURE
c
c
c      te0 = dmax1(timin,tei,dsign(2.d4,-tei))
      teinit = dmax1(timin,tei)
      tef = teinit
      tm1 = teinit
c     
      tl1 = dlog(tm1)
      if (nmod .eq. 'TIM') then
         call timion(tm1, edens, hdens, xhyf, tstep)
      else
         call equion(tm1, edens, hdens)
      end if
c     
      call cool(tm1, edens, hdens)
c     
      if (dabs(dlos).lt.delmin) goto 500
      if (dabs(dabs(dlos)-1.0d0).lt.1.d-6)
     &     dlos = dsign(1.0d0-1.d-6,dlos)
      if (dabs(dlos).gt.0.06d0) nf = nf+2
      if (dabs(dlos).gt.0.20d0) nf = nf+3
      c1 = arctanh(dlos)
      tl2 = tl1-dlos
      t2 = dmax1(dexp(tl2),4.0d0*tmin)
      tl2 = dlog(t2)
      tti = t2
      dl1 = dlos
      tdef = tm1
      dlodef = dl1
      varf = 0.02d0+(dabs(dl1) ** 0.70d0)
      jok = 0
c     
c     
c     ***USING 2 TEMP. (TM1,T2)  AND 2 FRACT. RESID. (DLOS1,DLOS2)
c     IT FITS THE FUNCTION : ARCTANH(DLOS)=B+A*LN(T)
c     NEXT VALUE OF T=DEXP(-B/A)
c
      do 100 n = 1, nf
         dlpr = dlos
         if (nmod .eq. 'TIM') then
            call copinto(popul, pop)
            call timion(t2, edens, hdens, xhyf, tstep)
         else
            call equion(t2, edens, hdens)
         end if
         call cool(t2, edens, hdens)
         if (dabs(dabs(dlos)-1.0d0).lt.1.d-6)
     &        dlos = dsign(1.0d0-1.d-6,dlos)
         if (dabs(dlos).le.dabs(dlodef)) then
            dlodef = dlos
            tdef = t2
         end if
         c2 = arctanh(dlos)
         csub = c2-c1
         tsub = tl2-tl1
         if (dabs(tsub).lt.1.0d-36) tsub = dsign(1.0d-36,tsub)
         if (tsub .eq. 0.0d0) tsub = dsign(1.0d-36,csub)
c     **DOES NOT ALLOW  A  TO BE NEGATIVE
         a = dabs(csub/tsub)
         b = c1-(a*tl1)
         cc = c1
         c1 = c2
         tl1 = tl2
         tm1 = t2
         tl2a =-(b/(dmax1(dabs(a),dabs(b)/20.0d0)*dsign(1.d0,a)))
         dt12 = tl2a-tl1
         tl2 = tl1+dsign(dmin1(varf,dabs(dt12)),dt12)
         t2 = dexp(tl2)
         ftu = (t2-tm1)/dmax1(t2,tm1,tmin/10.d0)
         iconsis = idint((1.1d0*ftu0)*ftu)/(1.d-20+dabs(ftu0*ftu))
         ftu0 = ftu
c     
         dvv = dabs(tl2-tl1)
         if (iconsis.lt.0) then
            isis = isis+1
            varf = varf/(2.0d0+(1.d0/(isis*isis)))
            tl2 = tl1+dsign(dmin1(varf,dabs(dt12)),dt12)
            t2 = dexp(tl2)
         else if (iconsis.gt.0) then
            varf = dmin1(1.1d0,dmax1(5.0d-4,(1.3d0*varf)*(dmax1(1.5d-1,
     &           dmin1(1.d0,dvv/varf)) ** 0.2d0)))
         end if
c     
c     *PRINT OUT WHEN TESTED BY SUBR. TESTI
      if (ntest .ne. 'Y') goto 50
      if (n .ne. 1) goto 30
      write(luop, 43) 
      write(luty, 43) 
 43   format(1h0,t4,3hTM1,t14,4hDLOS,t27,2hT2,t39,1hA,t49,1hB,t58,2hC1
     &        ,t68,2hC2)
      write(luty, 33) teinit, dl1, tti
      write(luop, 33) teinit, dl1, tti
 33   format(1h ,1pg12.5,g10.3,g12.5)
 30   write(luty, 37) tm1, dlos, t2, a, b, cc, c2, isis
      write(luop, 37) tm1, dlos, t2, a, b, cc, c2, isis, varf
 37   format(1h ,1pg12.5,g10.3,g12.5,4g10.3,i4,1pg10.3)
      write(luop, 850) tloss, hloss, rloss, fslos, floss, fflos
     &, colos, pgain, rngain
 850  format(1h ,9(1pg10.2))
 50   continue
c     
      fadl = dmin1(100.0d0,dmax1(0.5d0,dabs(dlpr)/(dabs(dlos)+1.d-4)))
      delexi = 7.0d-4+(2.0d-4*dmin1(5.d0,6.d-1*fadl))
      if (dabs(dlos).lt.2.d-3) jok = jok+1
      if (((dabs(dlos).lt.delexi).or.(t2.lt.tmin)).or.((jok.gt.
     &       3).and.(dabs(dlos).lt.(dble(jok)*1.3d-4)))) goto 200
 100  continue
c     
c     ***IF CONVERGENCE POOR , USES TEMP. WITH MIN. ABS(DLOS)
c     
 200  if (dabs(dlos).lt.0.02d0) goto 210
      t2 = tdef
      write(luty, 120) dlodef, tdef
      if (ntest .eq. 'Y') write(luop, 120) dlodef, tdef
 120  format(44h CONVERGENCE FOR EQUIL. TEMP. IS TOO SLOW   ,3hDL:
     &     ,1pg9.2,5h  TE:,0pf9.0)
 210  continue
c     
c     ***USES FINAL TEMPERATURE
c     
      if (dabs(dlos).gt.(2.0d0*dabs(dlodef))) t2 = tdef
      tef = dmax1(t2,30.d0)      
 220  if (nmod .eq. 'TIM') then
         call copinto(popul, pop)
         call timion(tef, edens, hdens, xhyf, tstep)
      else
         call equion(tef, edens, hdens)
      end if
      if (tef.lt.200.d0) write(*,*) 
     &     "WARNING: temperature may be unrealistic"      
c     
c      WRITE(*,*) N,DL1,TEinit,DLODEF,TDEF,DLOS,TEF
c
      call cool(tef, edens, hdens)
 500  continue
c
c
c
      return 
c
      end
      
