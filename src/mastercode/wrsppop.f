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
c*******TO OUTPUT ON MEDIUM*LUOP  THE EMISSION SPECTRUM AND
c    THE AVERAGE TEMPERATURES AND IONIC POPULATIONS
c    CALL SUBROUTINE SPECTRUM
c
c
c     modified for varable numbers of elemnts and ionisation stages
c     
c     RSS 8/90
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wrsppop(luop)
c
      include 'cblocks.inc'
c
c           Variables
c
      double precision achh,dhas,dhtt
      integer*4 luop
      character mode*4
c
c     FOR population total calculation
c      double precision popintot(mxion, mxelem)
c      integer*4 i,j
c
c           Functions
c
      double precision densnum,favcha
c
      achh = favcha(pam,1)
      dhtt = densnum(1.d0)
      dhas = deav/((dhtt*achh)+1.d-38)
      write(luop, 2043) zgas
c
c
 2043 format(/' ',t9,'Abundances of the elements relative',
     &' to Hydrogen and their average state of ionisation :','  (Zgas='
     &,f7.4,' Zsun)'/)
c
      call wionabal(luop,pam)
c
c
c
      write(luop, 2123) 
 2123 format(//' ',t9,'Average distances for a given state ',
     &'of ionisation :'/)
c
      call wionpop(luop,rdisa)
c
c
c
      write(luop, 2120) 
 2120 format(//' ',t9,'Average ionic temperatures :'/)
c
      call wionpop(luop,teav)
c
c
c
      write(luop, 2130) 
 2130 format(//' ',t9,'Integrated column densities:'/)
c
c skip maxion+1 bin
c       do j = 1, atypes
c          do  i = 1, maxion(j)
c            popintot(i,j)=popint(i,j)
c          enddo
c       enddo
      call wionpop(luop,popint)
c
c
c
      write(luop, 2150) deav, dhas, achh, dhtt
 2150 format(//,' ',t9,'Average ionic electron dens. :',2x,'<De>='
     & ,1pg10.3,3x,'<DH>=',1pg10.3,3x,'<CHARGE>=',1pg10.3,3x,'Dnum/DH :'
     & ,1pg11.4,/)
c
      call wionpop(luop,deam)
c
c
c
      write(luop, 2080) 
 2080 format(/' ',t9,'Line weighted averaged quantities :'/)
      write(luop, 2085) 
 2085 format(' ',t11,'<T>NII',t21,'<De>NII',t31,'<T>OIII',t41,'<De>OIII'
     &,t51,'<DE>SII',t61,'<T>SII',t71,'<DE>OII',t81,'<T>OII'/)
      write(luop, 2065) tnii, denii, toiii, deoiii, desii, tsii
     &, deoii, toii
c
c
c
 2065 format(' ',t8,2(0pf10.0,1pg10.3),2(1pg10.3,0pf10.0)/)
      write(luop, 2233) 
 2233 format(' ',t9,'Photon field averaged quantities :'/)
      write(luop, 2237) 
 2237 format(' ',t11,'<ZETAE>',t21,'<QHDH>'/)
      write(luop, 2243) zetaeav, qhdha
 2243 format(' ',t8,2(1pg10.3)/)
c
c    ***OUTPUT EMISSION SPECTRUM RELATIVE TO H-BETA
c

      mode = 'REL'
c
      call spectrum(luop, mode)
c
      return 
      end
