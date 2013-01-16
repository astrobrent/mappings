cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     functions to lookup and interpolate structure arrays using 
c     bisections
c     assume the index array is monotonic - ie radius or time
c     from small to large. Linear interps.
c     using many functions so the usage is clear and easy to understand.
c     and to prevent sillines like get temperature from electon density
c
c     functions: fGetdhFromRadius   (get the density of h atoms - use
c                                      frho to get density from this)
c                fGetdhFromTime     (from time coordinate instead)
c                fGetdeFromRadius   (get electron denisty)
c                fGetdeFromTime     (get electron denisty, from time)
c                fGetTeFromRadius   (temperature)
c                fGetTeFromTime     ( from time)
c                fGetVelFromRadius   (velocity)
c                fGetVelFromTime     ( from time)
c                fGetxhFromRadius   (hydrogen ion fraction)
c                fGetxhFromTime     ( from time)
c                fGetAbundFromRadius   (abundance of a species from 
c                 radius - ie in novae)
c                fGetAbundFromTime     ( from time)
c
c
      integer*4 function GetIndexUpwards(value, np, data)
c      
      include 'const.inc'
c
      double precision value,data(nstmax)
      integer*4 np, index,lowerbound,upperbound
c      
      index = 1
c
      upperbound = np
      lowerbound = 1
c
c     bisection search
c
 100  if((upperbound-lowerbound).gt.1) then
c
            index = (upperbound+lowerbound)/2
c
c            write(*,*) index,lowerbound,upperbound
c     
            if (value.ge.data(index)) then
               lowerbound = index
	        else
               upperbound = index
            endif

            goto 100
         
      endif
c     
      index = lowerbound
c
      GetIndexUpwards = index
c  
      end
      
      integer*4 function GetIndexDownwards(value, np, data)
c      
      include 'const.inc'
c
      double precision value,data(nstmax)
      integer*4 np, index,lowerbound,upperbound
c
      index = 2
c      
      if (np.gt.2) then
c
         upperbound = np
         index = 1 + np/2  
         lowerbound = 1
c
c     bisection search
c
 100     continue
c     
         if (value.lt.data(index)) then
            upperbound = index
            index = lowerbound+(upperbound-lowerbound)/2
         endif
         
         if (value.ge.data(index)) then
            lowerbound = index
            index = lowerbound+(upperbound-lowerbound)/2
         endif
         
         if (((value.lt.data(index)).or.(value.ge.data(index)))
     &        .and.(index.gt.1).and.(index.le.np)) goto 100
c     
         if ((index.lt.2).or.(index.ge.np)) then
            write(*,*) 'Error in interpolation: Off the profile array!'
     &      ,index
            stop
         endif
c     
        endif
c
        GetIndexDownwards = index
  
      end

      double precision function fGetdhFromRadius( r )
c
      include 'cblocks.inc'
c
      double precision r,dhl
      integer*4 index,GetIndexUpwards
c
c     assumes radius info in dist(nstmax) and nsteps are up to date
c
c
      dhl = 0.d0
c
      index = GetIndexUpwards(r,nsteps,dist)
c          
      dhl = (r-dist(index-1))/(dist(index)-dist(index-1))
      dhl = (dhl*(dhy(index)-dhy(index-1)))+dhy(index-1)
c
      fGetdhFromRadius = dhl
c
      return
c
      end

c
c
      double precision function fGetdhFromTime( t )
c
      include 'cblocks.inc'
c
      double precision t,dhl
      integer*4 index,GetIndexUpwards
c
c     assumes time info in timlps(nstmax) and nsteps are up to date
c
c
      dhl = 0.d0
c
      index = GetIndexUpwards(t,nsteps,timlps)
c     
      dhl = (t-timlps(index-1))/(timlps(index)-timlps(index-1))
      dhl = (dhl*(dhy(index)-dhy(index-1)))+dhy(index-1)
c
      fGetdhFromTime = dhl
c
      return
c
      end

c
      double precision function fGetdeFromRadius( r )
c
      include 'cblocks.inc'
c
      double precision r,del
      integer*4 index,GetIndexUpwards
c
c     assumes radius info in dist(nstmax) and nsteps are up to date
c
      del = 0.d0
c
      index = GetIndexUpwards(r,nsteps,dist)
c     
      del = (r-dist(index-1))/(dist(index)-dist(index-1))
      del = (del*(deel(index)-deel(index-1)))+deel(index-1)
c
      fGetdeFromRadius = del
c
      return
c
      end
c
c
      double precision function fGetdeFromTime( t )
c
      include 'cblocks.inc'
c
      double precision t,del
      integer*4 index,GetIndexUpwards
c
c     assumes time info in timlps(nstmax) and nsteps are up to date
c
c
      del = 0.d0
c
      index = GetIndexUpwards(t,nsteps,timlps)
c     
      del = (t-timlps(index-1))/(timlps(index)-timlps(index-1))
      del = (del*(deel(index)-deel(index-1)))+deel(index-1)
c
      fGetdeFromTime = del
c
      return
c
      end

c
      double precision function fGetTeFromRadius( r )
c
      include 'cblocks.inc'
c
      double precision r,Tel
      integer*4 index,GetIndexUpwards
c
c     assumes radius info in dist(nstmax) and nsteps are up to date
c
c
      Tel = 0.d0
c
      index = GetIndexUpwards(r,nsteps,dist)
c     
      Tel = (r-dist(index-1))/(dist(index)-dist(index-1))
      Tel = (Tel*(te(index)-te(index-1)))+te(index-1)
c
      fGetTeFromRadius = Tel
c
      return
c
      end

c
c
      double precision function fGetTeFromTime( t )
c
      include 'cblocks.inc'
c
      double precision t,Tel
      integer*4 index,GetIndexUpwards
c
c     assumes time info in timlps(nstmax) and nsteps are up to date
c
c
      Tel = 0.d0
c
      index = GetIndexUpwards(t,nsteps,timlps)
c     
      Tel = (t-timlps(index-1))/(timlps(index)-timlps(index-1))
      Tel = (Tel*(te(index)-te(index-1)))+te(index-1)
c
      fGetTeFromTime = Tel
c
      return
c
      end

c
      double precision function fGetVelFromRadius( r )
c
      include 'cblocks.inc'
c
      double precision r,Vell
      integer*4 index,GetIndexUpwards
c
c     assumes radius info in dist(nstmax) and nsteps are up to date
c     small to large radius
c
      Vell = 0.d0
c
      index = GetIndexUpwards(r,nsteps,dist)
c     
      Vell = (r-dist(index-1))/(dist(index)-dist(index-1))
      Vell = (Vell*(veloc(index)-veloc(index-1)))+veloc(index-1)
c      write(*,*) ' Velocity interp: '
c      write(*,*) index, index-1
c      write(*,'(3(1pg14.7,1x))') Vel1, veloc(index), veloc(index-1)
c      write(*,'(3(1pg14.7,1x))') r, dist(index), dist(index-1)
c     
      fGetVelFromRadius = Vell
c
      return
c
      end

c
c
      double precision function fGetVelFromTime( t )
c
      include 'cblocks.inc'
c
      double precision t,Vell
      integer*4 index,GetIndexUpwards
c
c     assumes time info in timlps(nstmax) and nsteps are up to date
c
c
      Vell = 0.d0
c
      index = GetIndexUpwards(t,nsteps,timlps)
c     
      Vell = (t-timlps(index-1))/(timlps(index)-timlps(index-1))
      Vell = (Vell*(veloc(index)-veloc(index-1)))+veloc(index-1)
c     
      fGetVelFromTime = vell
c
      return
c
      end
c
      double precision function fGetXhFromRadius( r )
c
      include 'cblocks.inc'
c
      double precision r,Xhl
      integer*4 index,GetIndexUpwards
c
c     assumes radius info in dist(nstmax) and nsteps are up to date
c
c
      Xhl = 0.d0
c
      index = GetIndexUpwards(r,nsteps,dist)
c     
      Xhl = (r-dist(index-1))/(dist(index)-dist(index-1))
      Xhl = (Xhl*(xh(index)-xh(index-1)))+xh(index-1)
c
      fGetXhFromRadius = Xhl
c
      return
c
      end
c
c
      double precision function fGetXhFromTime( t )
c
      include 'cblocks.inc'
c
      double precision t,Xhl
      integer*4 index, GetIndexUpwards
c
c     assumes time info in timlps(nstmax) and nsteps are up to date
c
c
      Xhl = 0.d0
c
      index = GetIndexUpwards(t,nsteps,timlps)
c     
      Xhl = (t-timlps(index-1))/(timlps(index)-timlps(index-1))
      Xhl = (Xhl*(xh(index)-xh(index-1)))+xh(index-1)
c
      fGetXhFromTime = Xhl
c
      return
c
      end
c
      double precision function fGetAbundFromRadius( species, r )
c
      include 'cblocks.inc'
c
      double precision r,Abnl
      integer*4 index,species,GetIndexUpwards
c
c     assumes radius info in dist(nstmax) and nsteps are up to date
c
c
      Abnl = 0.d0
c
      index = GetIndexUpwards(r,nsteps,dist)
c     
      Abnl = (r-dist(index-1))/(dist(index)-dist(index-1))
      Abnl = (Abnl*(abn(species,index)-abn(species,index-1)))
     &       +abn(species,index-1)
c
      fGetAbundFromRadius = Abnl
c
      return
c
      end
c
c
      double precision function fGetAbundFromTime( species, t )
c
      include 'cblocks.inc'
c
      double precision t,Abnl
      integer*4 index,species,GetIndexUpwards
c
c     assumes time info in timlps(nstmax) and nsteps are up to date
c
c
      Abnl = 0.d0
c
      index = GetIndexUpwards(t,nsteps,timlps)
c     
      Abnl = (t-timlps(index-1))/(timlps(index)-timlps(index-1))
      Abnl = (Abnl*(abn(species,index)-abn(species,index-1)))
     &       +abn(species,index-1)
c
      fGetAbundFromTime = Abnl
c
      return
c
      end



