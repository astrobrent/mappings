c*******************************************
c
c   include file for global constants. 

c   (cf August 1996 Physics Today, 75th CRC Handbook 1995)
c
c   v3.0.0q MAPPINGS III
c
c*******************************************
c
      implicit none
c
      double precision pi,ftpi
      double precision rt2,euler
      double precision fpi,ifpi
      double precision epsilon, logepsilon
      double precision huge, loghuge
      double precision plk,cls,ev,eesu
      double precision me,mp,amu,mh,mh1
      double precision bohr,stefan,nav,rgas
      double precision ryd,rydlambda,LmeV
      double precision evplk,ryde,pzlimit
      double precision pars,rka,rkb,ffk,fk
      double precision kbb
      integer*4 mxinfph,mxelem,mxion
      integer*4 mxdust,mxdtype,pahi
      integer*4 mxTempbin,dinfph
      integer*4 nstmax,nionedges
c
c*******************************************
c
c        Universal constants
c
c*******************************************
c
c   planck's constant (erg seconds) 
c
      parameter  (plk = 6.6260755d-27)
c
c   speed of light in vacuum (cm/s)
c
      parameter  (cls = 2.99792458d10)
c
c     1 electron Volt in ergs
c
      parameter (ev = 1.60217733d-12)
c
c     1 electron Volt/h
c
      parameter (evplk = ev/plk)
c
c   electron charge in esu units
c
      parameter (eesu = 4.8032067991d-10)
c
c     rest mass of an electron me (g)
c
      parameter (me = 9.1093897d-28 )
c
c     rest mass of a proton mp (g)
c
      parameter (mp = 1.6726231d-24 )
c
c     Nav Avogadro's number 
c
      parameter (nav = 6.0221367d23 )
c
c     mean mass of a H atom : 1.00797 amu
c
      parameter (mh = 1.673774705394e-24 )
c
c     mass of a 1H atom : 1.007825 amu
c
      parameter (mh1 = 1.673533927065e-24 )
c
c   Ionisation potential of hydrogen (ev) Primary Def.
c
      parameter (ryd = 13.59844d0)
c
c
c   Ionisation potential of hydrogen (ergs)
c
      parameter (ryde = ryd*ev)
c
c
c   Wavelength of Hydrogen Limit (Angstoms) Secondary
c
      parameter (rydlambda = 911.7534400)
c
c
c   Wavelength (Angstoms) to Energy (eV) Primary Def.
c
      parameter (LmeV = 12398.424449)
c
c
c   First Bohr radius  (cm)
c
c   a0 = h^2/(4 pi^2 me e^2
c
c
      parameter (bohr = 5.29177249d-9)
c
c     atomic mass unit assuming C12 == 12
c
      parameter (amu = 1.6605402d-24 )
c
c   Stefan-Boltzmann constant
c
c   (ergs cm^-2 s^-1 K^-1) 
c
c      2 pi^5 k^4
c   = ------------
c      15 c^2 h^3
c
      parameter (stefan = 5.6705085385d-5)
c
c   collection of constants for coll. rate = sqrt((2*pi*hbar^4)/(k me^3))
c
c   in ergs and cm units (SI use d-12)
c
      parameter (rka = 8.62915d-6)
c
c   boltzman's constant k, ergs per kelvin
c
      parameter (rkb = 1.380658d-16)
c
c     continuum emission constant
c
c    16/3 (pi/6)^0.5 (eesu^6 / (c^3 me^2)) (me/k)^0.5 
c  = 16/3 (pi/(6 k me^3))^0.5  (e^2 / c )^3 
c
      parameter(ffk = 5.44435784954d-39)
c
c
c     total FreeFree emission constant
c
c    (16/(3 hbar me)) ((2 pi  k )/(3 me))^(1/2)  ((eesu^2)/c)^3
c
c
      parameter(Fk = 1.42553656926202256d-27)
      
c
c constant for plank function 1/((cls**2)*(plk**3))
c
      parameter ( kbb = 3.824630234294027d57)

c
c   universal gas constant J K-1 M-1 using K and nav above
c
      parameter ( rgas = 8.3145112d+7 )
c
c   1 parsec (cm)
c
      parameter ( pars = 3.08568d18)
c
c*******************************************
c
c        Mathematical constants
c
c*******************************************
c
c   pi and four thirds pi
c
      parameter ( pi = 3.1415926535897d0, ftpi = (4.d0/3.d0)*pi )
      parameter (fpi = 4.d0*pi          , ifpi =  1.d0/(4.d0*pi))
c
c   Euler' constant to 20 places
c
      parameter ( euler = 0.57721566490153286061d0 ) 
c
c   square root of 2
c
      parameter ( rt2 = 1.414213562373095 ) 
c
c*******************************************
c
c        Program constants
c
c*******************************************
c
c   epsilon is a small number to prevent division by zero
c   and nasty logs. it is also used to test quantities to see
c   if they are worth calculating (see intvec), 1.0d-38 is chosen
c   to match 1E-38 in single precision on a VAX, note:
c   smallest single precision number on sun is ~1.17e-38 
c
c
c   Changed to 1e-99 for double precision
c
      parameter (epsilon = 1.0d-99, logepsilon = -2.2795592420641d2)
c
c
      parameter (huge = 1.0d99, loghuge = 2.2795592420641d2)
c
c     limit for abundant ions: ab * pop
c
      parameter (pzlimit = 1.0d-12)
c
c     maximum number of atoms (elements)
c
      parameter (mxelem = 16)
c
c     maximum number of ionization stages
c
      parameter (mxion = 29)
c
c     radiation vector size
c
      parameter (mxinfph = 1802)
c
c     maximum number of dust sizes (36=MRN, 81=FULL)
c     need to change DUSTDAT files as well (.mrn <==>.full)
c
      parameter (mxdust = 81)
c
c     maximum number of dust types (ie graphite,silicate etc)
c
      parameter (mxdtype = 2)
c
c     Number of PAH ions     
c
      parameter(pahi=5)
c
c     ionsation edges vector size
c
      parameter (nionedges = 310)
c
c     model stucture step array limit
c
      parameter (nstmax = 1024)
c
c
      parameter (mxTempbin=401)      !Max No. of Temp. bins
c
      parameter (dinfph=401)         !Max No of dust photon bins
c
