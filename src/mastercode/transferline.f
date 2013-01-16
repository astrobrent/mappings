ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function tauline(dismul, tau0)
c
c tau * dismul subject to an upper limit of 1.d10
c
      include 'const.inc'
c
      double precision dismul,  tau0
      double precision tau_eff
c
      tau_eff  = tau0
c      
      if ( dismul .gt. 1.d0) then 
           tau_eff  =  tau0*dismul
      endif
c
      tauline = min(tau_eff,  1.d10)
c
      return
c  
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function transferout(dismul, tau)
c
c  Simply the attenuation of radition across a depth tau
c  with no sources along the path.
c  The fraction retained in the zone (ie absorbed) is 
c  1-transferout.
c

      include 'const.inc'
c
      double precision dismul, tau, tau_eff
      double precision escape, tauline
c                  
      tau_eff = tauline(dismul, tau)
c
      escape  = dexp(-tau_eff)
c
      transferout = escape
      return
c  
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       double precision function localout(dismul, tau)
c
c  A zone producing a uniform field.  Emission at one side
c sees a different optical depth to escape than at the other for a given direction
c The local production of energy escaping sees an average of exp(-tau)
c and the fraction that escapes is (1-exp(-tau))/tau.  The fraction
c retained in the zone is 1-escape
c
      include 'const.inc'
c
      double precision dismul,  tau , tau_eff
      double precision escape, tauline
c                  
      tau_eff = tauline(dismul, tau)
c
      escape = 1.d0
      if (tau_eff .ge. 1.0d-6) then
            escape = (1.d0-dexp(-tau_eff))/tau_eff
      endif
c      
      localout = escape
c      
      return
c  
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function meanfield(dismul, tau)
c
c A field crossing a zone depth tau, attenuated to exp(-tau)
c at the end has a mean field in the zone of (1-exp(tau))/tau
c Exactly the same as the one-side escape fraction of a locally produced
c field!  Called meanfield just to distinguish the case from local 
c production.
c
c The difference is that the amount escaping is exp(-tau), not
c 1 - meanfield.  This function averages a crossing field with
c no local sources.
c
      include 'const.inc'
c
      double precision dismul,  tau , tau_eff
      double precision mean, tauline
c                  
      tau_eff =tauline(dismul, tau)
c
      mean = 1.d0
      if (tau_eff .ge. 1.0d-6) then
            mean = (1.d0-dexp(-tau_eff))/tau_eff
      endif
c      
      meanfield = mean
c
      return
c  
      end
c
