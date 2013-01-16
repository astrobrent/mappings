cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******CALCULATES IONISATION RATE DUE TO SECONDARY ELECTRONS
c	FOR THE ELEMENTS INCLUDED IN MATRIX SECEL
c	OUTPUT RATES IN RASEC(3,11)
c
c	REF.: SHULL,M.J. APJ.234,P761 (1979)
c	      BERGERON,J.,SOUFFRIN,S. A&A 14,P167
c
c
      subroutine ionsec()
c
c
      include 'cblocks.inc'
c
c
c           Variables
c
      double precision a1,a2,abr,abwe,aij,ath,atx,b1
      double precision b2,deno,eah,eau1,eau2,eax,eij,emin
      double precision eps,epsh,epsij,fhiieff,frio
      double precision poef,poij,ps,valel,wei,wij,wth,wtx
      double precision x,zshull,ztot
      double precision sig,beff
c
      integer*4 i,j,k
      integer*4 atom,ion
c
      sig(eps) = dmin1(1.0d0-(1.0d0/eps),1.0d0/eps)+epsilon
      beff(x) = (0.0352d0*(1.d0-((x)**(-(0.1d0*dlog10(x))))))
     &/dmax1(1.d-20,1.d0-x)
c
      eau1 = ryd
      eau2 = 30.0d0
      emin = ryd
c
      zshull = 1.1d0
      do 51 i = 1, atypes
      do 50 j = 1, 3
         rasec(j,i) = 0.0d0
   50 continue
   51 continue
c
      if ((anr(1,1,1).le.0.0d0).and.(wnr(1,1,1).le.0.0d0)) goto 420
      wth = 0.0d0
      wtx = 0.0d0
      ath = 0.0d0
      atx = 0.0d0
c
      do 101 i = 1, atypes
      do 100 j = 1, maxion(i)-1
         abr = zion(i)*pop(j,i)
         ath = ath+(abr*anr(1,j,i))
         atx = atx+(abr*anr(2,j,i))
         wth = wth+(abr*wnr(1,j,i))
         wtx = wtx+(abr*wnr(2,j,i))
 100  continue
 101  continue
c
      if (ath.le.0.0d0) goto 420
      eah = (ath/(wth+1.d-38))/ev
      eax = (atx/(wtx+1.d-38))/ev
      deno = dlog(eau2)-dlog(eau1)
      b1 = (dlog(atx+1.d-38)-dlog(ath+1.d-38))/deno
      a1 = dlog(ath+1.d-38)-(b1*dlog(eau1))
      b2 = (dlog(eax+1.d-38)-dlog(eah+1.d-38))/deno
      a2 = dlog(eah+1.d-38)-(b2*dlog(eau1))
      ztot = 0.0d0
      ps = 0.0d0
      wei = 0.0d0
c
c
      do 300 k = 1, 25
      atom = secat(k)
      if (atom.le.0) goto 300
      ion = secio(k)
      poij = epot(ion,atom)/ev
      poef = dmax1(poij,emin)
      eij = dexp(a2+(b2*dlog(poef)))+poef
      aij = dexp(a1+(b1*dlog(poef)))
      epsij = eij/poij
      epsh = eij/eau1
      valel = dble(secel(k))
      wij = ((aij*(valel*((eau1/poij) ** 2)))*sig(epsij))/sig(epsh)
      secra(k) = wij
      ztot = ztot+zion(atom)
      abwe = zion(atom)*wij
      ps = ps+(abwe*pop(ion,atom))
      wei = wei+abwe
  300 continue
c
      fhiieff = dmax1(1.d-20,1.0d0-(ps/(wei+1.d-36)))
      frio = ((zshull*beff(fhiieff))/ztot)/ev
c
      do 400 k = 1, 25
      atom = secat(k)
      if (atom.le.0) goto 400
      ion = secio(k)
      rasec(ion,atom) = frio*secra(k)
  400 continue
c
  420 return 
      end

