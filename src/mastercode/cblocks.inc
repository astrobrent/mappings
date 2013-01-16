ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Mappings III 1.0.0r COMMON BLOCK GLOBAL VARIABLES
c     
c     common blocks are all in this one include file
c     this ensures identical declarations for all files.
c     
c     Need documentation and more organisation....
c     will document each variable as they are edited.
c     
c     RSS 8/90-96
c
c     Version 1.0.0r
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'const.inc'
c     
c     Need some global enviroment variables
c     
      character*8 theVersion, launchstate
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   declare ATDAT block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   zion: gas phase abundance fractions
c   zion0: total (undepleted) abundance fractions
c   dion: depletion (dust) factors (mulipliers)
c   invdion: inverse of depletion factors
c   zen: sum of all species abundance fractions
c           ie zen*dh = total number of ions
c   atwei: atomic weights
c   maxion: maximum ionisation stage for a given atom type
c   zmap&mapz : convert atomic numbers to array entries and v v 
c   elem:   Atomic symbols
c   rom :   Ion roman numerals 
c
      double precision zion(mxelem),zion0(mxelem)
      double precision atwei(mxelem),dion(mxelem)
      double precision invdion(mxelem),zen
      integer*4 maxion(mxelem)
      integer*4 zmap(mxion),mapz(mxion)
      character elem(mxelem)*4,rom(mxion)*8
c
c      declare ARCOL block variables
c   arrays by shells,ion,elem + nshells to control loops
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Arnaud and Rothenflug collisional ionisation coeffs
c   see collion.f and COLLDAT
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   col: collisional ionisation rates by ion and atom
c
      double precision ipot(3,mxion,mxelem)
      double precision aar(3,mxion,mxelem),bar(3,mxion,mxelem)
      double precision car(3,mxion,mxelem),dar(3,mxion,mxelem)
      integer*4 nshells(mxion,mxelem)
c
c   Landini& Monsignori Fosse & Shull coll. ion. fits
c
      double precision acol(mxion, mxelem),tcol(mxion, mxelem)
c
      double precision col(mxion, mxelem)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Original mappings collisional ionisation coeffs
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   (also uses epot from photio)
c
c      double precision beta(6,11),effn(6,11)
c
c      declare ELABU block variables
c
c   pop: ion fractional populations
c   dndt: rate of change of pop with time step
c   popint: integrated column densities
c
      double precision pop(mxion, mxelem)
      double precision dndt(mxion, mxelem)
      double precision popint(mxion, mxelem)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare PHOTIO block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   epot: ion valence ionisation potential (ergs)
c   epot2: ion ionisation potetial from no+1 level (ev)
c       (used in ffbgau.f)
c       atpho,ionpho: atom and ion entry id for photoionisation
c                     data arrays
c       ipotpho: ionisation edge potential
c       augpho: Auger electron potential
c       sigpho: threshold crossection
c       betpho,spho: beta and s for crossection as f(e) fit
c       stwtpho: stat weight for ion ground state (2J+1) RS coupling
c   secel: secondary ionisation data array
c   xsec : mean photo crossection for each bin
c   xsect: total crossection over ionisisin band
c   .
c
      double precision epot(mxion, mxelem)
      double precision epot2(mxion, mxelem)
      integer*4 atpho(nionedges),ionpho(nionedges)
      integer*4 secat(26),secio(26),secel(28)
      double precision ipotpho(nionedges)
      double precision augpho(nionedges)
      double precision sigpho(nionedges)
      double precision betpho(nionedges)
      double precision spho(nionedges)
      double precision stwtpho(nionedges)
      double precision secra(26)
      double precision heaph(5,mxion, mxelem), anr(2,mxion, mxelem)
      double precision wnr(2,mxion, mxelem),rasec(mxion, mxelem)
      double precision rphot(mxion, mxelem),auphot(mxion, mxelem)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare RECOMB block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   arad,xrad: radiative recombination rate coeffs by ion and atom
c   adi,bdi,t0,t1: dielectronic recombination coeffs
c   a-f dilt: low T dielectronic rate coeffs
c   rec: total recombination rates by ion and atom
c
c   recspl: spline coeffs
c   nrcspl: number of ions spline recomb data
c   recsplat: the atoms
c   recsplio: the ions
c   recspln: number of spline points
c
c   oprec: True if the new data in OPRECDAT for Ar + S to be included 
c
      double precision arad(mxion, mxelem),xrad(mxion, mxelem)
      double precision adi(mxion, mxelem),t0(mxion, mxelem)
      double precision bdi(mxion, mxelem),t1(mxion, mxelem)
      double precision adilt(mxion, mxelem),bdilt(mxion, mxelem)
      double precision cdilt(mxion, mxelem),ddilt(mxion, mxelem)
      double precision fdilt(mxion, mxelem)
      double precision rec(mxion, mxelem)
      double precision recspl(10,16,5)
      integer nrcspl,recsplat(10),recsplio(10),recspln(10)
	  logical oprec
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare PHOTDAT block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   ephot: energies of photon bin edges in eV: from file PHOTDAT
c   soupho: source photon vector
c   tphot: total local photon vector (ergs)
c   optphot: total local optical photon vector (ergs)
c   emidif,updif,dwdif: diffuse field photon vectors 
c       (local,upstream and downstream)
c   emlin,uplin,dwnlin: resonance emission line arrays with distance 
c       multipliers (local,upstream and downstream)
c   skipbin: logical array to optimise calculation of photionisation
c       to skip bin if it has too few photons to bother with.
c   dilradpho effective source radius for finite extent plane parallel 
c   models as a funtion of energy, goes with rshock in s4
c   infph: number of bins in vectors 
c   linpos: line bin position information.
c
c       ffphot,fbphot,p2phot: continuum spectra (in ergs)
c                            (free-free,fre-bnd,2photon)
c   dustphot: photon field absorbed by dust locally
c
      double precision ephot(mxinfph),tphot(mxinfph)
      double precision soupho(mxinfph), wdil
      double precision difftot(mxinfph),dilradpho(mxinfph),rshock
      double precision emidif(mxinfph),updif(mxinfph),dwdif(mxinfph)
      double precision xsec(mxinfph), xsect
      double precision dustphot(mxinfph)
      integer*4 infph,linpos(2,60)
      integer*4 FieldVersion, IonStartBin
      logical skipbin(mxinfph)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   declare gaunt factor calculation globals
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   arg,brg,crg,drg,erg: resonace gaunt factor parameters
c   
c
        integer*4 clkup(19,19)
        double precision arg(38),brg(38),crg(38),drg(38),erg(38)
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare CHARXC block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   charge exchange reaction rate variables.....
c   charco and charte:old mappings: chargemode = 1
      integer*4 nchxold
      double precision charco(7,26),charte(5,26)
c
c   ref Arnaud and Rothenflug 1985, A.A.Suppl.Ser. 60:425
c
c   recombination reactions
c
c   nchxr: number of reactions used
c   chxrat: element #
c   chxrion: ion stage
c   chxrx: reactant (H or He)
c   chxrtemp: suggested temp range
c   chxrcos: coefficients
c       1:a
c       2:b
c       3:c
c       4:d
c
c   chxi series equivalent
c   chxicos: coefficients
c       1:a
c       2:b
c       3:c
c       4:delta E (eV)
c
      integer*4 nchxr,chxrat(55),chxrio(55),chxrx(55)
      double precision chxrtemp(2,55),chxrcos(4,55),chxr(55)
c
c   ionisation reactions
c
      integer*4 nchxi,chxiat(15),chxiio(15),chxix(15)
      double precision chxitemp(2,15),chxicos(4,15),chxi(15)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       declare COOLPR block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   cooling rate quantities (ergs)
c   alphaCloss, alphaClaw, alphaC0: artificial powerlaw losses , 
c      index and normalisation used if alphacoolmode = 1
c   colos: collisional ionisation losses
c   fflos: free-free radiation losses
c   cploss: compton losses := cmpcool, cmpheat separate heating 
c      and cooling processes
c   chgain: charge exchange heating
c   pgain: photionisation heating
c   rngain: recombination heating
c   gheat: photoelectric heating by large grains
c   tloss :total absolute loss
c   eloss: effective loss
c   egain: effective heating
c   wgain: heating due to pdv work
c   dlos: measure of equilibration
c   tem: current electron teperature
c   fixtemp: fixed electron temp for jthm = 'T' models
c   cosgain:cosmic ray hydrogen heating
c   crate: cosmic ray primary event rate (~10^-17 s^-1)
c   cosphi: cosmic ray secondary ionisations
c
      double precision alphaCloss, alphaClaw, alphaC0
      double precision colos,fflos,chgain,pgain,rngain,tem
      double precision cosgain,crate,cosphi,cmplos, cmpcool, cmpheat
      double precision tloss,eloss,egain,dlos, fixtemp, wgain
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Dust Grains
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Bgrain: threshold parameter(large}
c   Yinf: yield parameter (larges)
c   haccom: accomodation parameter
c   segrain,spgrain: sticking for electrons and protons
c   siggrain: area of large grains per H atom
c   Avgrpot: Average grain potential in volts
c   grainpot: grain potential for each size
c   amin,amax: grain size limits (microns)
c   galpha: distribution index, default to n^-3.5
c   gheat,gcool:heating and cooling due to grains
c   paheat:  heating due to PAH molecules
c   pahfrac: # of PAHs per H
c   pahCfrac: % of C in PAHs 
c   pahyield: the PE yield curve for all PAH ion states
c   pahion: ionisation state of PAH
c   pahIP: Ionisation potential of that PAH state
c   pahZ: fraction of PAHs in that ion state
c   paheng: energy absorbed per sec by PAH (=energy emitted)
c   pahn/iabs: PAH absorbtion curve for positive and neutral/negative 
c        states
c   pahn/isca: PAH scattering as above
c   pahn/icos: PAH scattering angle (as above)
c   pahint:integrated column of PAH
c    ionising vector, using same bins.
c   extinct(i,k,j) dust extinction curve for grains, size k and typej, 
c    (Q ext values)
c   scatter(i,k,j) dust scattering curves for grains, size k and typej,
c    (Qsca values)
c   absorp(i,k,j) dust absorption curve for grains, size k and typej,
c    (Q abs values)
c   Gcos(i,k,j) dust scattering angle curve for grains, size k and type
c    (<cos(theta)> value)
c   
c   dustint:integrated column of Hydrogen for dust extinction
c   dustbinmax: number of dust radius bins, 
c   mindust,maxdust: bins corresponding to min and max grainsizes
c
c   npahx, npahy : number of points in the curves
c
c   grainrad,gradedge: grain radius bin centres and edges
c   dustnum, dustsig: number and crossection distribution of dust
c   dcrosec(inl): extinction cross-sections by dust at wvl inl
c   dfscacros,dbscacros:forward and back scattering crossections 
c
c   IRphot: infrared flux
c
c   ClinPAH: Condition for Carbon to be linked to PAH existence
c

      double precision grainrad(mxdust),gradedge(mxdust+1),
     &           dustmass(mxdtype)
      double precision dustnum(mxdust,mxdtype),dustsig(mxdust,mxdtype)
      double precision amin(mxdtype),amax(mxdtype),graindens(mxdtype)
      double precision Bgrain,Yinf,haccom
      double precision gheat,gcool
      double precision segrain,spgrain
      double precision pahint(pahi),dustint
      double precision grainpot(mxdust,mxdtype),avgpot(mxdtype),galpha
      double precision extinct(mxinfph,mxdust,mxdtype), 
     &           absorp(mxinfph,mxdust,mxdtype)
      double precision scatter(mxinfph,mxdust,mxdtype), 
     &		 Gcos(mxinfph,mxdust,mxdtype)
      double precision dcrosec(mxinfph,mxdtype),
     &           dfscacros(mxinfph,mxdtype),dbscacros(mxinfph,mxdtype)
      double precision IRphot(mxinfph),dustAv,A_vend
      integer*4 dustbinmax,mindust(mxdtype),maxdust(mxdtype),numtypes
      logical IRtemp, ClinPAH
c
c PAHs
c
      double precision paheat,pahfrac,pahCfrac
      double precision pahIP(pahi),pahZ(pahi),pahion(pahi),pahlimit
      double precision pahyield(mxinfph,pahi),paheng,pahflux(mxinfph)
      double precision pahnabs(mxinfph),pahnsca(mxinfph),
     &           pahncos(mxinfph),pahnext(mxinfph)
      double precision pahiabs(mxinfph),pahisca(mxinfph),
     &           pahicos(mxinfph),pahiext(mxinfph)
      character pahend
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision engy(mxinfph)
      double precision phot(mxinfph)
      double precision IRFlux(mxinfph)
      double precision Cabs(mxinfph)
      double precision photabs(mxinfph)
      double precision dustsigt(mxinfph)
      double precision collheat(mxinfph)
      
      double precision edphot(dinfph)
      double precision dengy(dinfph)     
      double precision ddE(dinfph)
      double precision v_e(dinfph)
      double precision dabsphot(dinfph),delE
     
      double precision BBemiss(mxinfph,mxTempbin)
      double precision T_grid(mxTempbin)
      double precision T_edge(mxTempbin)
      double precision H_grid(mxTempbin) 
      double precision Hmin(mxTempbin)
      double precision Hmax(mxTempbin) 
      double precision deltaH(mxTempbin)
      double precision invH(mxTempbin, mxTempbin)
      double precision Tr_matrix(mxTempbin, mxTempbin)
      double precision T_Prob(mxTempbin)
      
      double precision grainvol(mxdust)
      double precision atom_no(mxdtype,mxdust)
      
      double precision flxdist(mxinfph,mxdust)
      double precision d_rad(mxdust)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      ** LINE DATA
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare XLINE block vaiables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   new resonance (permitted) line data
c   xlam: wavelengths (units?)
c   xrbri: line brightnesses (ergs)
c   rloss: resonance line heat losses (in rline)
c   920 ytransitions included from M&L 1990
c   
c
      integer*4 xlines,xrat(920),xion(920),xiso(920),xbin(920)
      double precision xrlam(920),xrbri(920),xbr(920),xrloss
      double precision xejk(920),xfef(920),xtrans(920),brr(19,19)
      double precision xegj(920),fluxx(920)
      double precision emilin(2,920),uplin(2,920),dwlin(2,920)
c
      integer*4 xilines,xiat(32),xiion(32),xiiso(32),xibin(32)
      double precision xilam(32),xibri(32),xiloss
      double precision xiejk(32),xiegj(32),xomeg(32),xitr(32)
      double precision xilin(2,32),xiuplin(2,32),xidwlin(2,32)
      double precision fluxxi(32)
c
      integer*4 xhelines,xheat(26),xheion(26),xhebin(26)
      double precision xhelam(26),xhebri(26),xheloss
      double precision xhejk(26),xhef(26),xhebra(16)
      double precision xhelin(2,32),xheuplin(2,32), xhedwlin(2,32)
      double precision fheif(26)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare RLINE block vaiables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   resonance (permitted) line data
c   rlam: wavelengths (units?)
c   rbri: line brightnesses (ergs)
c   rloss: resonance line heat losses
c   lrbin : location in diffuse field bins for each line
c    rlin(2,60),ruplin(2,60), rdwlin(2,60) lines with distance 
c    multipliers like emilin etc
c    fluxr: integrated fluxes
c
      double precision rlam(60),rbri(60),rloss
      integer*4 lrbin(60)
      integer*4 nlines,ielr(60),ionr(60)
      double precision rlin(2,60),ruplin(2,60),rdwlin(2,60)
      double precision fluxr(60)
c
c   resonance line strengths data
c
      double precision e12r(60),omres(60),fabsr(60)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare CLINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   semi-forbidden intercombination line information
c
c   fslam: wavelengths
c   fsbri: brightnesses
c   fslos :semi forbidden cooling rate (ergs)
c   intam: ion/atom identification 
c       (will be changed to numerical arrays
c       and use the global ion & roman numeral character arrays)
c   lcbin : location in diffuse field for each line
c
      double precision fslam(80),fsbri(80),fslos
      integer*4 lcbin(80)
      integer*4 mlines,ielfs(80),ionfs(80)
      character intam(2,80)*4
      double precision fslin(2,80),fsuplin(2,80),fsdwlin(2,80)
      double precision fluxi(80)
c
c      semi-forbidden line strengths data
c
      double precision e12fs(80),w1fs(80),w2fs(80),a21fs(80),omfs(80)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare F3LINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   three level atom line information
c   f3lam: wavelengths
c   f3bri: brightnesses
c   f3loss: forbidden cooling rate (ergs)
c   f3ion,f3atom: forbiden line ion identification,
c   read from threedat
c   f3type: type of three level lines, 1 = finestructure
c   2 = forbidden + others
c
      double precision f3lam(3,60),f3bri(3,60),f3loss
      integer*4 f3ion(60),f3atom(60),nf3ions,nf3trans
      integer*4 f3type(60),f3bin(3,60)
      double precision f3lin(2,3,60),f3uplin(2,3,60),f3dwlin(2,3,60)
      double precision fluxf3(3,60)
c
c   three level line strengths data
c
      double precision wi3(3,60),ei3(3,60),ai3(3,60),gam3(2,3,60)
      double precision pow3(2,3,60)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare FLINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Five Level line information
c   flam: wavelengths
c   fbri: brightnesses
c   floss: forbidden cooling rate (ergs)
c   fion,fatom: forbiden line ion identification,
c   read from forbidat
c   fluxf: integrated lines
c
      double precision flam(10,40),fbri(10,40),floss
      integer*4 fion(40),fatom(40),nfions,nftrans,fbin(10,40)
      double precision flin(2,10,40),fuplin(2,10,40),fdwlin(2,10,40)
      double precision fluxf(10,40)
c
c   Five Level line strengths data
c
      double precision wi(5,40),ei(5,5,40)
      double precision ai(5,5,40),omi(5,5,40),tdep(5,5,40)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare F6LINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Six level line information (room for 5 species)
c   f6lam: wavelengths
c   f6bri: brightnesses
c   f6loss: six level cooling rate (ergs)
c   f6ion,f6atom: six level line ion identification,
c
      double precision f6lam(15,5),f6bri(15,5),f6loss
      integer*4 f6ion(5),f6atom(5),n6ions,n6trans,f6bin(15,5)
      double precision f6lin(2,15,5),f6uplin(2,15,5),f6dwlin(2,15,5)
      double precision fluxf6(15,5)
c
c   Six Level line strengths data (room for 5)
c
      double precision wi6(6,5),ei6(6,6,5)
      double precision ai6(6,6,5),omi6(6,6,5),tdep6(6,6,5)
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare F9LINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Nine level line information
c   f9lam: wavelengths
c   f9bri: brightnesses
c   f9loss: nine level cooling rate (ergs)(Fe VII only for now, room 
c      for 3)
c   f9ion,f9atom: nine level line ion identification,
c
      double precision f9lam(36,3),f9bri(36,3),f9loss
      integer*4 f9ion(3),f9atom(3),n9ions,n9trans,f9bin(36,3)
      double precision f9lin(2,36,3),f9uplin(2,36,3), f9dwlin(2,36,3)
      double precision fluxf9(36,3)
c
c   Nine Level line strengths data (room for 3 species)
c
      double precision wi9(9,3),ei9(9,9,3)
      double precision ai9(9,9,3),omi9(9,9,3),tdep9(9,9,3)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare FELINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Iron II multiplet line information
c
c   felam: wavelengths
c   febri: brightnesses
c   feloss: forbidden cooling rate (ergs)
c   read from feiidat
c
      double precision felam(160),febri(160),feloss
      integer*4 nfetrans,felamap(160),febin(160)
      character*12 feid(160)
      double precision fluxfe(160)
      double precision felin(2,160),feuplin(2,160), fedwlin(2,160)
c
c   Fe II line strengths data
c
      double precision wfe(24),efe(24,24),afe(24,24),omfe(24,24)
c
c   data structures for fitted collision strengths
c   ist ion/line to use this is CII 158mu
c   make room for 10 ions, 10points, 5 values
c   will ususally use log(t) for x value.
c
c   data lives in INTERDAT2 for now as only one 
c   fs line uses this now. This data is intened to 
c   be generally used for aspects of different lines.
c   other data such as wavelenght and A values are 
c   stored separately with the other transition data.
c
c   omgspl: spline coeffs
c   nomg: number of ions    
c   omgspat: the atoms
c   omgspio: the ions
c   omgnsp: number of spline points
c
      double precision omgspl(10,16,5)
      integer nomg,omgspat(10),omgspio(10),omgnsp(10)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Hydrogen data.....
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   New for hydrec2:
c
c   hspline: cubic spline fits to accurate hydrogenic recombination 
c    rate for total: hspline(1,x,coeffs), sum of 210 levels: 
c    hspline(2,x,coeffs), and n=1 , l = 0: hspline(3,x,coeffs)
c    x = log(t), result = log(rate)
c
c   hcoln: n -> n' collisional rate fits
c   dnl: level term for nl->nl' collisions
c
c   Line series:
c
c   hlambda(line,series) ,ie Hbeta = hlambda(2,2) (Balmer line 2)
c
c
c    th42(10),the43(12)   temperatures for hydrogen & Helium spec calc 
c      data points
c    rhhe(7)  densities (n) for  hydrogen & Helium spec calc data 
c      points (same for h and he)
c    h42a(10,7),h42b(10,7) lookup tables for hydrogen hbeta emission
c    he43a(12,7),he43b(12,7) lookup tables for helium 4686 emission
c    Use powerlaw interpoalations and extrapolations
c    hylratsa(line,series,temps,rho) id = same order as lambdas:
c      ie lyman series then balmer etc
c    hylratsb(line,series,temps,rho)      ten lines each for six series
c    helratsa(line,series,temps,rho)
c    helratsb(line,series,temps,rho)
c
c    collrate2p(mxelem): collision rate to 2S level of H-like atoms 
c       from ground only
c    rateton(mxelem,n): total collision rate to level n of H-like atoms
c
c
c    Use lyman and Balmers series only for heavier elements, scaled
c    from hydrogen hylrats(a,b) data, (unlike helium which is done 
c    properly)
c    Generally the lines will be assumed to be in case A until this 
c    causes a problem.
c       
c      double precision xhlambda(10,2,mxelem)
c      double precision xhydroflux(10,2,mxelem)
c      double precision xhydrobri(10,2,mxelem)
c      double precision xhbin(10,2,mxelem)
c      double precision xhydlin(2,10,2,mxelem)
c      double precision xhyduplin(2,10,2,mxelem)
c      double precision xhyddwlin(2,10,2,mxelem)
c
      double precision balcoe(5,16),hbema(5,6),hbemb(5,6),hherat(14)
      double precision hheen(14),hlam(10),hbri(10),h2ql
      double precision h2qav, hei2qa, heii2qa,fluxh(10),fhbeta
c
      double precision collrate2p(mxelem), rateton(mxelem,16)
      double precision recrate2p(mxelem), h2q(mxelem),h2qbri(mxelem)
      double precision r2s1a(10,7),r2s1b(10,7)
      double precision r2s2a(12,7),r2s2b(12,7)
c      
      double precision hspline(3,18,5)
      double precision hlambda(10,6),helambda(10,6)
      double precision th42(10),the43(12),rhhe(7)
      double precision h42a(10,7),h42b(10,7)
      double precision he43a(12,7),he43b(12,7)
      double precision hylratsa(10,6,10,7), hylratsb(10,6,10,7)
      double precision helratsa(10,6,12,7), helratsb(10,6,12,7)
      double precision hydroflux(10,6),heliflux(10,6)
      double precision heiioiiibf, heiioiiibfsum
      double precision hydrobri(10,6),helibri(10,6)
      double precision hydrogf(10,6),hbin(10,6),hebin(10,6)
      double precision hydlin(2,10,6), hyduplin(2,10,6), 
     &           hyddwlin(2,10,6)
      double precision hellin(2,10,6), heluplin(2,10,6), 
     &           heldwlin(2,10,6)
c
      double precision xhlambda(10,2,mxelem)
      double precision xhydroflux(10,2,mxelem)
      double precision xhydrobri(10,2,mxelem)
      double precision xhbeta(mxelem)
      double precision xhbin(10,2,mxelem)
      double precision xhydlin(2,10,2,mxelem)
      double precision xhyduplin(2,10,2,mxelem)
      double precision xhyddwlin(2,10,2,mxelem)
c
      double precision hbeta, hloss, caseab(mxelem)
c
c   Collisional Excitation
c
c   ccoln(56,4),dcoln(56,4): collision strength fits for hydrogenic 
c     lines
c   colid(56,4): collisional level ids for ccoln and dcoln
c   collha(line,series,coll.level), collhb(line,series,coll.level):
c     Contribution to line in series from collisions to coll.level
c     for case a and b, assuming 'n method'
c
c
      double precision ccoln(56,4),dcoln(56,4)
      integer*4 colid(56,4)
      double precision collha(10,6,15),collhb(10,6,15)

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Helium data......
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision heii2ql,heien(4),qren(4)
      double precision hein(4),heiilr,hel0(6,12)
      double precision heilam(16), heibri(16), fluxhei(16)
      integer*4  heibin(16), nheilines

c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Continuum (free-free, f-b and two photon) emmision data
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   zn0 = effective charge of ion
c   ezz = 2 photon 1s-2s energy gaps (from grotian diags) (eV)
c   fn0 = quantum nunber of the ground state
c   zt0 = # of vacancies available
c   ffph,fbph,p2ph: arrays of continuum component photons
c   cnphot: total continuum array (photons), uses ephot bins
c
c   gff(41,81),gffg2(41),gffu(81) -free free gaunt factor lookup table
c   gfb(20,50),gfbe(50),gfbtz(20) -free bound gaunt factor lookup table
c
c   gffintspl: (41,5) cubic spline fits to accurate integrated 
c     freefree emission gfactor 
c   x = log(gamma^2), result = <gff(g^2)>
c   coeffs: 1 = log g2, 2 = a0, 3 = a1, 4=a2, 5 = a3
c
      double precision zn0(mxion, mxelem),ezz(mxion, mxelem)
      integer*4 fn0(mxion, mxelem),zt0(mxion, mxelem)
      double precision ffph(mxinfph),fbph(mxinfph),p2ph(mxinfph)
      double precision cnphot(mxinfph)
c
      double precision gff(41,81),gffg2(41),gffu(81)
      double precision gfb(20,50),gfbe(50),gfbtz(20)
      double precision gffintspl(41,5)
      integer*4 ngffg2,ngffu,ngfbe,ngfbtz
c
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Summed (integerated) model quantities
c
c      declare SLINE block variables
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision deav,teav(mxion, mxelem),rdisa(mxion, mxelem),
     &           pam(mxion, mxelem)
c
c      declare TONLIN block variables 
c
      double precision vunilog
c
c      declare DENLIN block variables
c
      double precision roiii,rnii,weoiii
      double precision wenii,deoiii,denii,toiii,tnii
      double precision deam(mxion, mxelem),dhav,roii
      double precision rsii,weoii,wesii,toii,tsii
      double precision deoii,desii,qtosoh,qhdh,zetae,qhdha,zetaeav
c     
c   declare PHOBAT block variables
c
c   abc! what the hell is abc3?  Anyway  I am changing them
c   as I figure them out. 
c   ofac, xfac, pfac, vfac, afac,bfac,cfac, efac, ffacm & scalen are  
c   density law parameters:
c   f(x) = ofac*[xfac*exp((x-afac)/scalen) + pfac*((x/afac)**bfac)
c                vfac*(((r0-x)/efac)**-ffac)]+cfac
c   tofac,txfac, tpfac,tafac,tbfac,tcfac are the same temperature 
c      profile
c   aef: accretion mass-luminosity efficiency approx 0.1
c   mdot: accretion rate
c
      double precision dhavv,viofr,aef
      double precision pop0(mxion, mxelem),wdpl,telap,tlife,remp,rmax
      double precision fren,tauen,tend,diend,dtau0,qhdin
      double precision qhdav,tm00,abc3,mdot
      double precision ofac,xfac,pfac,afac,bfac,cfac,scalen
      double precision vfac, efac,ffac
      double precision tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
      double precision abc4,abc5,abc6,abc7,abc8
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      **SHOCK FRONT DATA :
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare VELSHO block variables
c
      double precision abc9,abc0
      double precision dhpo,vshoc,vpo,vfin,hmago
      double precision fm,pres,ww,qst,qfin,depr
c
c
c   Declare flow and stucture variables
c   pr0,1: pressures at point a and b.
c   rho0,1: densities (gcm^-3)
c   te0,1: electron temps
c   vel0,1: flow velocities (cm s-1)
c   dh0,1: Hydrogen number densities
c   de0,1: electron number densities
c   bm0,1: magnetic field (Gauss)
c
      double precision pr0,rho0,te0,vel0,dh0,de0,bm0
      double precision pr1,rho1,te1,vel1,dh1,de1,bm1
      double precision dis0,dis1,disav,rad0,rad1
      double precision gamma,GM,flowconst
c
c      declare SHOARR block variables
c
c   arrays to store shock models, nstmax = 1024 (defined in const.inc)
c
      double precision depo
      double precision dhy(nstmax),xh(nstmax),dist(nstmax),
     &           veloc(nstmax)
      double precision deel(nstmax),te(nstmax), abn(mxelem,nstmax)
      double precision popi(mxion, mxelem),popf(mxion, mxelem)
      double precision popini(mxion, mxelem)
      double precision propop(mxion, mxelem),timlps(nstmax),
     &           distim(2,0:nstmax)
      double precision hoiii(nstmax),hoii(nstmax)
      double precision hnii(nstmax),hsii(nstmax)
      double precision cav(5),diluf,xhmin,fval,texi
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      ** CHARACTER  and INTEGER type common blocks
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      declare JBLOCK block variables
c
c   atypes: number of different atoms used (up to mxelem)
c   collmode: 0 or 1 : collisional ionisation calculation mode(def 0)
c   dielcmode: 0 or 1 :dielectr. recomb. calculation mode(def 0)
c   mapmode: 0 or 1 :overall mode, new=0, old = 1(default 0)
c   photonmode: 0 or 1: no coupling=0, normal coupling = 1 
c        (should be 1)
c   chargemode: 0 or 1: new charex = 0, oldcharex = 1 (default 0)
c   grainmode: T or F: nograins = F, include grains = T (default 0)
c   pahmode: T or F: if grains exist T, if not, F default (F)
c   IRmode: 0,1 or 2: nograins = 0,if grainsand not want IR or 
c        dusttemp=0, want IR =1, want IR and dusttempdist=2 
c        (default = 0)
c   hellcoolmode: 0 or 1: thick/thin Helium resonance cooling 
c        (default 0)
c   alphacoolmode: 0 or 1: physical =0, or powerlaw=1 cooling 
c        (default 0)
c   expertmode: 0 or greater, expert diagnostic mode (default 0)
c
c   ??????????????????????????
c   iphom,ipho
c   jhe2p,j108m
c   ionum,limph
c   ielen,jpoen
c   istf,loop,iter,iel1,iel2
c   nsteps: total number of steps in structure arrays (ie deel(nstmax))
c   jspot,jpre,jfin,jend
c   jden,jgeo,jeq,iso: structure settings
c     jden - density; C=isoChoric, B=isoBaric, F=function, I=input file
c     jgeo - geometry; S=spherical, P=plane, F=??
c     jeq  - equil.; E=equilibrium, F=finite source, P=post equil decay, 
c            (C=CIE)
c     iso  - ??
c   jthm: thermal structure - should be 'S' for self consistent 
c         (or T for isothermal)
c
c   jcon : flag to enable full continuum calculation
c         (used in photo4,sh2,sh3)
c   jspec: flag to enable optical flambda spectrum integration
c           (in newdif,totphot,photo4)
c   jall: flag to enable production of allion file
c   jiel: flag to enable monitoring of 4 elements
c   jsou: flag to enable production of final source file
c       jpfx: final source file prefix
c
c   jbal: flag to enable production of balance files
c   jbfx: balance file prefix 
c   gdist:  grain distribution
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   general run info:
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   runname: code for a given run, included in all output files
c   ionsetup: Ionisation balance setup
c   abnfile: name of abundance file or 'default'
c   srcfile: name of photon source file or 'none'
c   depfile: depletion pattern file or 'default'
c   structfile: radius - dh file or 'none'
c   runmode: 'interactive' or 'batchstart' or 'batchrun'    
c
      integer*4 atypes,collmode,dielcmode,mapmode,photonmode
      integer*4 iphom,ipho,jhe2p,j108m,chargemode,expertmode
      integer*4 ionum,limph,ielen,jpoen,ieln,hellcoolmode, alphacoolmode
      integer*4 istf,loop,iter,iel1,iel2,iel3,iel4,nsteps
      integer*4 ox1,ox2,ox3,ni2,su2,lygline
      integer*4 IRmode
      logical grainmode,pahmode
      character jspot*4,jpre*4,jfin*4,jend*4,jsou*4,jall*4
      character jbal*4, jbfx*8, jiel*4, jthm*4, gdist*4
      character jden*4,jgeo*4,jeq*4,iso*4,jcon*4,jspec*4,jpfx*8
      character runname*64,ionsetup*80,abnfile*80,srcfile*80
      character runmode*64,depfile*80, structfile*80
c     
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      stellar atmosphere data
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   star blocks 16 is nothing to do with atoms
c
      double precision tno4(16,11),tno5(16,11),tmet(8,11)
      double precision tedge(16,3),turn
      double precision alnth,teff,cut,zgas,zstar,rstar
c
c   source total photon fluxes (not just stars)
c
      double precision qht,qhi,qhei,qheii,qtlog
c
c   ** Globals for avrdata,summdata,and z* routines
c   ** accounting stuff
c
c   declare ACCUM variables
c
      double precision tlosac,tforbi,fheiilr,vas111
      double precision vas222,vas333,vas444,vas555
c
c
c
      common /environ/ theVersion,launchstate
c
c
c      **COMMON BLOCKS ,FIRST : ATOMIC DATA
c
c
      common /atdat1/ zion,zion0,atwei,dion,invdion,zen
     &                ,maxion,zmap,mapz
      common /atdatchar/ elem,rom
c
      common /arcol1/ ipot,aar,bar
     &                ,car,dar ,nshells  ,acol,tcol ,col
c
      common /elabu1/ pop,dndt,popint
c     
      common /phonn1/ atpho,ionpho  ,ipotpho,augpho,sigpho
     &                ,betpho,spho,stwtpho  ,epot,epot2,secra
     &                ,secat,secio,secel,xsec ,heaph,anr,xsect
     &                ,wnr,rasec ,rphot,auphot
c
      common /recom1/ arad,xrad,adi,t0,bdi,t1
     &                ,adilt,bdilt,cdilt,ddilt,fdilt ,rec,recspl
     &                ,nrcspl,recsplat,recsplio,recspln,oprec
c
      common /phovec1/ ephot,tphot,wdil ,emidif,updif,dwdif
     &                ,emilin,uplin,dwlin ,soupho,dilradpho,difftot
     &                ,dustphot,infph,linpos,FieldVersion
     &                ,skipbin,IonStartBin 
c
      common /resga1/  arg,brg,crg,drg,erg,clkup
c
      common /charex1/ charco,charte ,nchxr,chxrat,chxrio,chxrx
     &                ,chxrtemp,chxrcos,chxr ,nchxi,chxiat,chxiio,chxix
     &                ,chxitemp,chxicos,chxi
c
c   loss rates
c
      common /coolr1/ colos,fflos,chgain,pgain,rngain,tem,cmplos,wgain
     &                ,gheat,gcool,paheat,pahfrac,pahCfrac
     &                ,tloss,eloss,egain,dlos,cosgain, crate, cosphi
     &                ,cmpcool, cmpheat, alphaCloss, alphaClaw, alphaC0
c
c   dust grains
c
      common /grain1/ grainrad,gradedge,dustmass,amin,amax,dustnum
     &                ,Bgrain,Yinf,haccom,galpha
     &                ,segrain,spgrain,grainpot,avgpot,dustsig,graindens
     &                ,dustint,extinct,absorp,scatter,Gcos
     &                ,dustbinmax,mindust,maxdust,numtypes
     &                ,IRphot,dcrosec,dfscacros,dbscacros
     &		      ,dustAv,A_vend,IRtemp, ClinPAH
c
      common /PAHs/   pahIP,pahZ,pahion,pahlimit
     &                ,pahyield,pahflux,paheng,pahint
     &		      ,pahnabs,pahnsca,pahncos,pahnext
     &		      ,pahiabs,pahisca,pahicos,pahiext
      common /PAHchar/   pahend
c
c NewIR dust parameters
c
      common /dustpar1/  engy, phot,Cabs, photabs, dustsigt,collheat
      common /dustpar2/  edphot, dengy, ddE, v_e, dabsphot,delE
      common /dustpar3/  BBemiss,T_grid,T_edge,Tr_matrix,T_Prob
      common /dustpar4/  H_grid, Hmin, Hmax, deltaH, invH
      common /dustpar5/  grainvol, atom_no
      common /dustpar6/  flxdist, d_rad, IRflux
c
c   line calculation globals
c
      common /xline1/ xrlam,xrbri,xbr,xrloss
     &                ,xejk,xfef,xtrans,xegj
     &                ,xlines,xrat,xion,xiso,xbin
c
      common /xilin1/ xilam,xibri,xiloss,brr
     &                ,xiejk,xomeg,xitr,xiegj
     &                ,xilin,xiuplin,xidwlin
     &                ,xilines,xiat,xiion,xiiso,xibin

      common /xhlin1/ xhelam,xhebri,xheloss
     &                ,xhejk,xhef,xhebra
     &                ,xhelin,xheuplin, xhedwlin,fheif
     &                ,xhelines,xheat,xheion,xhebin
c
      common /rline1/ rlam,rbri,rloss
     &                ,rlin,ruplin, rdwlin
     &                ,e12r,omres,fabsr
     &                ,nlines,ielr,ionr,lrbin
c
      common /cline1/ fslam,fsbri,fslos
     &                ,fslin,fsuplin, fsdwlin
     &                ,e12fs,w1fs,w2fs,a21fs,omfs
     &                ,mlines,ielfs,ionfs,lcbin
      common /clinechar/ intam
c
      common /fline1/  flam,fbri,floss
     &                ,wi6,ei6,ai6,omi6,tdep6
     &                ,fion,fatom,nfions,nftrans
     &                ,flin,fuplin, fdwlin,fbin
     &                ,wi,ei,ai,omi,tdep
c
      common /f6line1/ f6lam,f6bri,f6loss
     &                ,f6lin,f6uplin, f6dwlin,f6bin
     &                ,f6ion,f6atom,n6ions,n6trans
c
      common /f9line1/ f9lam,f9bri,f9loss,f9bin
     &                ,wi9,ei9,ai9,omi9,tdep9
     &                ,f9lin,f9uplin, f9dwlin
     &                ,f9ion,f9atom,n9ions,n9trans
c
      common /f3line1/ f3ion,f3atom,f3type,nf3ions,nf3trans
     &                ,f3lam,f3bri,f3loss,f3bin,f3lin,f3uplin, f3dwlin
     &                ,wi3,ei3,ai3,gam3,pow3
c
      common /feline/ felam,febri,feloss
     &                ,felin,feuplin, fedwlin
     &                ,wfe,efe,afe,omfe
     &                ,felamap,febin, nfetrans
      common /felinechar/ feid
c     
      common /linest/ omgspl,nomg,omgspat,omgspio,omgnsp
c     old 
      common /hline1/  balcoe,hbema,hbemb,hherat,hheen,hlam,hbri
     &                ,h2ql,hbeta,caseab,hloss
     &                ,hspline,hlambda,helambda,th42,the43,rhhe
     &                ,h42a,h42b,he43a,he43b,hbin,hebin
     &                ,hylratsa, hylratsb, helratsa, helratsb
     &                ,hydroflux,heliflux,hydrobri,helibri,hydrogf
     &                ,hydlin,hyduplin, hyddwlin
     &                ,hellin,heluplin, heldwlin
     &                ,heiioiiibf, heiioiiibfsum
     &                ,collrate2p, rateton, recrate2p,h2q,h2qbri
     &                ,r2s1a, r2s1b, r2s2a, r2s2b
     &                ,xhlambda,xhydroflux,xhydrobri,xhbeta
     &                ,xhbin,xhydlin,xhyduplin,xhyddwlin
     &                ,ccoln,dcoln,collha,collhb,colid
c     
      common /helio1/  heii2ql,heien,qren,hein,heiilr,hel0
     &                ,heilam,heibri,fluxhei
     &                ,heibin,nheilines
c
c
c   continuum calculations
c
c
      common /cndat1/ zn0,ezz
     &                ,fn0,zt0
     &                ,cnphot,ffph,fbph,p2ph
     &                ,gff,gffu,gffg2,gfb,gfbe,gfbtz,gffintspl
     &                ,ngffg2,ngffu,ngfbe,ngfbtz
c
c     
c     SUMMED LINES
c
      common /sline1/ fluxr,fluxi,fluxf,h2qav,hei2qa,heii2qa
     &                ,fluxh,fhbeta,pam,deav,teav
     &                ,rdisa,vunilog,fluxf3,fluxf6,fluxf9
     &                ,fluxx,fluxxi,fluxfe
      common /tonlin/ roiii,rnii,weoiii,wenii,deoiii,denii,toiii,tnii
      common /denlin/ deam,dhav,roii,rsii,weoii,wesii,toii,tsii
     &                ,deoii,desii,qtosoh,qhdh,zetae
     &                ,qhdha,zetaeav,dhavv,viofr
c
c   PHOBAT  
c
      common /phobat/ pop0,wdpl,telap,tlife,remp,rmax,distim
     &                ,fren,tauen,tend,diend,dtau0,mdot,fixtemp
     &                ,qhdin,qhdav,tm00,aef
     &                ,abc3,abc4,abc5,abc6,abc7,abc8,abc9,abc0
     &                ,ofac,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac
     &	              ,ffac,tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
c
c   shock wave blocks
c
      common /velsho/ dhpo,vshoc,vpo,vfin,hmago,fm,pres
     &                ,ww,qst,qfin,depr,depo, rshock
      common /flow1/ pr0,rho0,te0,vel0,dh0,de0,bm0
     &                ,pr1,rho1,te1,vel1,dh1,de1,bm1
     &                ,gamma,GM,flowconst
c
      common /struc1/ dis0,dis1,disav,rad0,rad1
c
      common /shoarr/ dhy,xh,dist,veloc,deel,te,timlps,abn
     &                ,popi,popf,popini,propop
     &                ,hoiii,hoii,hnii,hsii
     &                ,cav,diluf,xhmin,fval,texi
c
c   Character  and integer global variables
c
      common /jints1/ atypes,collmode,dielcmode,mapmode,photonmode
     &                ,iphom,ipho,jhe2p,j108m,chargemode,nchxold
     &                ,ionum,limph,ielen,jpoen,ieln,nsteps
     &                ,istf,loop,iter,iel1,iel2,iel3,iel4
     &                ,ox1,ox2,ox3,ni2,su2,lygline
     &                ,alphacoolmode,expertmode,hellcoolmode
     &                ,grainmode,IRmode,pahmode
      common /jchars/ jspot,jpre,jfin,jend,jsou,jall,jbal,jbfx
     &                ,jden,jgeo,jeq,iso,jcon,jspec,jpfx,jiel, jthm
     &                ,gdist, runname,ionsetup,abnfile,srcfile,runmode
     &                ,depfile,structfile
c
c   star blocks
c
      common /star1/ tno4,tno5,tmet
     &                ,tedge,turn
     &                ,alnth,teff,cut,zgas,zstar,rstar
     &                ,qht,qhi,qhei,qheii,qtlog
c
      common /accum1/ tlosac,tforbi,fheiilr,vas111
     &                ,vas222,vas333,vas444,vas555
c