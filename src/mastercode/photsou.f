cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   MAPPINGS III.  An Astrophysical Plasma Modelling Code.
c
c   copyright 1997 Ralph S. Sutherland, Michael A. Dopita
c                        Luc Binette and Brent Groves
c       Version 1.0.0r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c*******TO ADD A PHOTOIONISATION SOURCE IN THE VECTOR SOUPHO(I)
c     EACH POINT CORRESPONDS TO THE INTENSITY OF RADIATION
c     AT THE MIDDLE OF THE ENERGY INTERVAL : EPHOT(I),EPHOT(I+1)
c     UNITS : ERGS.CM-2.SEC-1.STERAD-1.HZ-1    (INU )
c
c     Inu in 2D 1/pi units *not* 3D 1/(4pi) like the diffuse field vectors
c     This is so that geometric dilution is simpler etc
c
c     External files are all 3D 1/(4pi) because they are generally written
c     from diffuse field vectors.  If external stellar atmosphere files are
c     converted for use in MAPPINGS, the units will have to be checked and
c     either converted to 3D units or scaled down by 4 after readin.
c
c     THE NUMBER OF EMERGENT PHOTONS FOR DIFFERENT SPAN IN
c     ENERGY ARE CONTAINED IN QHI,QHEI,QHEII,QHT(=SUM OF THE 3)
c     
c     CALL SUBROUTINE FPOL,INTVEC
c     
c     **REFERENCES : HUMMER,D.G.,MIHALAS,D.M.(1970)MNRAS,147,P.339
c     SHIELDS,G.A.,SEARLE,L., (1978)APJ.222,P.830
c     
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine photsou(where)
c     
      include 'cblocks.inc'
c     
c
c           Variables
c
      double precision bv,del,dt0,dt1,dtot,etr,flu
      double precision xx(5), yy(5), al(5), cc(5), souvec(mxinfph)
      double precision scale,wid,den
      double precision blum,ilum,sv,tm1,pt0,gr
      double precision q1,q2,q3, fe1, fk1, fk2, fe2, alpha1, alpha2
      double precision eng1,eng2
      double precision q4,rnuh,tnu,yc,zmod,ze
c
      integer*4 fv,iinf,iinto,j,jg,jmax,ki,kii
      integer*4 kiii,luin,luop,n,n1,i
c
      character where*16, caract*80
      character fnam*80, ilgg*4
      character kstr*4, nf*4
c
      logical iexi
c
c     old (Draine 79 UV)
c     double precision fev,a1,a2,a3,
c     integer*4 fuvmin
c
c           Functions
c
      double precision fplank, fheavyside
c     
      luin = 17
      luop = 18
c     
c     ***ZERO VECTOR SOUVEC AND RELATED QUANTITIES
c     
      do i = 1, infph-1
         souvec(i) = 0.0d0
         skipbin(i) = .TRUE.
      enddo
      do  i = 1, 16
         tedge(i,3) = 0.0d0
         tedge(i,2) = 0.0d0
      enddo
c     
      ipho = ipho+1
      zstar = 0.0d0
      alnth = 0.0d0
      teff = 0.0d0
      cut = 0.0d0
      srcfile = 'None'
      scale = 1.0000d0
      blum = 0.d0
      ilum = 0.d0
      etr = ephot(infph)/ryd
c     
c     ***CHOICE OF THE TYPE OF SOURCE
c     
c     
      if (runmode.ne.'batchrun') then
      write(*,80) where
 80   format(//' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' Setting the radiation field for ',a16/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'//)
 85   format(/' Initial field is zeroed, successive source types'/
     &     ' adds to total field until exit.  The field may be',/
     &     ' re-initialised with the (Z) zero option'//)
      endif
c     
c     total energy in Inu
c     
c     
 100  blum = 0.d0
      ilum = 0.d0
c     
      do i = 1,infph-1
         wid = ephot(i+1)-ephot(i)
         blum = blum + souvec(i)*wid*evplk
         if (ephot(i).ge.Ryd) ilum = ilum + souvec(i)*wid*evplk
      enddo
c
c     Plane Parallel x pi
c
      blum = pi*blum
      ilum = pi*ilum
c
c
c    ***INTEGRATE NUMBER OF EMEGENT PHOTONS TO IONISE H,HE (=PI*JNU)
c    Assuming plane parallel geometry and dilution of 0.5 and intvec
c    assumes 1/4pi units and souvec is 1/pi plane parallel units.
c
      q1 = 0.0d0
      q2 = 0.0d0
      q3 = 0.0d0
      q4 = 0.0d0
c
      call intvec(souvec, q1, q2, q3, q4)
c
      qhi   = q1*0.25d0
      qhei  = q2*0.25d0
      qheii = q3*0.25d0
      qht   = q4*0.25d0
c    
      if (runmode.ne.'batchrun') then
         write(*, 110) blum,ilum,qhi,qhei,qheii,qht,infph,infph
      endif
 110  format(//' Multiple component photoionisation sources'/
     &     ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     ' : Total intensity   : ',1pg12.5,' (ergs/s/cm^2)      :'/
     &     ' : Ion.  intensity   : ',1pg12.5,' (ergs/s/cm^2)      :'/
     &     ' : FQHI  (1-1.8Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     ' : FQHeI (1.8-4Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     ' : FQHeII   (>4Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     ' : FQ tot   (>1Ryd)  : ',1pg12.5,' (phots/cm^2/s)     :'/
     &     ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Non-thermal component (c*nu**Alpha)'/
     &     '    B  :   Black body component (T>2500)'/
     &     '    C  :   4G stellar models (26000-56000K)'/
     &     '    E  :   Bremsstrahlung component (c*exp(-nu/kT))'/
     &     '    F  :   Local ISRF (Mathis etal (1993)'//
     &     '    G  :   AGN (Bland-Hawthorn et al 1997) '/
     &     '    H  :   AGN with BBB (Nagao et al 2001) '//
     &     '    I  :   Input data file    (',i4,' POINTS)'/
     &     '    O  :   Output data file    (',i4,' POINTS)'//
     &     '    Z  :   Zero (or Nil ) radiation'/
     &     '    S  :   Apply overall Scaling factor.'/
     &     '    N  :   Normalise fluxes.'//
     &     '    X  :   eXit with current source.'//
     &     ' :: ',$)
      iso = '    '
      read(*, 107) iso
 107  format(a)
      iso = iso(1:1)
      if (iso .eq. 'a') iso = 'A'
      if (iso .eq. 'A') goto 1000
      if (iso .eq. 'b') iso = 'B'
      if (iso .eq. 'B') goto 900
      if (iso .eq. 'c') iso = 'C'
      if (iso .eq. 'C') goto 200
      if (iso .eq. 'e') iso = 'E'
      if (iso .eq. 'E') goto 800
      if (iso .eq. 'f') iso = 'F'
      if (iso .eq. 'F') goto 500
      if (iso .eq. 'g') iso = 'G'
      if (iso .eq. 'G') goto 1100
      if (iso .eq. 'h') iso = 'H'
      if (iso .eq. 'H') goto 1300
      if (iso .eq. 'i') iso = 'I'
      if (iso .eq. 'I') goto 400
      if (iso .eq. 'o') iso = 'O'
      if (iso .eq. 'O') goto 300
      if (iso .eq. 'n') iso = 'N'
      if (iso .eq. 'N') goto 1200
      if (iso .eq. 'z') iso = 'Z'
      if (iso .eq. 'Z') goto 600
      if (iso .eq. 's') iso = 'S'
      if (iso .eq. 'S') goto 700
      if (iso .eq. 'x') iso = 'X'
      if (iso .eq. 'X') goto 2000
c     
      goto 100
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Zero field
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
 600  do i = 1,infph
         souvec(i) = 0.0d0
      enddo
      do  i = 1, 16
         tedge(i,3) = 0.0d0
         tedge(i,2) = 0.0d0
      enddo
c     
      ipho = ipho+1
      zstar = 0.0d0
      alnth = 0.0d0
      teff = 0.0d0
      cut = 0.0d0
      srcfile = 'None'
      scale = 1.0000d0
      blum = 0.d0
      ilum = 0.d0
      etr = ephot(infph)/ryd
c     
      goto 100
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     scale factor
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
 700  if (runmode.ne.'batchrun') then
      write(*,760)
 760  format(/' Give overall source scaling factor (>0): ',$)
      endif
      read(*,*) scale
      if (scale.le.0.d0) goto 700
c     
      do j = 1,infph
         souvec(j) = souvec(j)*scale
      enddo
c     
      goto 100
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Normalise factors
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
1200  if (runmode.ne.'batchrun') then
      write(*,1210)
1210  format(//' Multiple component photoionisation sources'/
     &     '  Choose the normalisation range :'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Normalise over H   I  up (1 Ryd +)'/
     &     '    B  :   Normalise over H   I (1 - 1.8 Ryd)'/
     &     '    C  :   Normalise over He  I (1.8 - 4 Ryd)'/
     &     '    D  :   Normalise over He II (4 Ryd +)'/
     &     ' :: ',$)
      endif
      read(*,'(a)') ilgg
      ilgg = ilgg(1:1)
c     
      if (ilgg .eq. 'a') ilgg = 'A'
      if (ilgg .eq. 'b') ilgg = 'B'
      if (ilgg .eq. 'c') ilgg = 'C'
      if (ilgg .eq. 'd') ilgg = 'D'
c
c    ***INTEGRATE NUMBER OF EMEGENT PHOTONS TO IONISE H,HE (=PI*JNU)
c    Assuming plane parallel geometry and dilution of 0.5 and intvec
c    assumes 1/4pi units and souvec is 1/pi plane parallel units.
c
      q1 = 0.0d0
      q2 = 0.0d0
      q3 = 0.0d0
      q4 = 0.0d0
c
      call intvec(souvec, q1, q2, q3, q4)
c
      qhi   = q1*0.25d0
      qhei  = q2*0.25d0
      qheii = q3*0.25d0
      qht   = q4*0.25d0
c
1220  if (runmode.ne.'batchrun') then
      write(*,1230)
1230  format(//' Multiple component photoionisation sources'/
     &     '  Normalise Flux to :'/
     &     ' ::::::::::::::::::::::::::::::::::::::::::::::::::::'/
     &     '    A  :   Normalise to 1 photon/cm^2/s'/
     &     '    B  :   Normalise to U = 1 (at 1 cm^-3)'/
     &     ' :: ',$)
      endif
      read(*,'(a)') nf
      nf = nf(1:1)
c
      if (nf .eq. 'a') nf = 'A'
      if (nf .eq. 'b') nf = 'B'
c
      scale = 1.d0
      if (nf.eq.'B') scale = cls
      if (ilgg.eq.'A') scale = scale/qht
      if (ilgg.eq.'B') scale = scale/qhi
      if (ilgg.eq.'C') scale = scale/qhei
      if (ilgg.eq.'D') scale = scale/qheii
c
      do j = 1,infph
         souvec(j) = souvec(j)*scale
      enddo
c     
      goto 100
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ***READ SOURCE FROM EXTERNAL FILE
c
c     file NEEDS "PHOTON SOURCE FILE" in header to be recognizable
c     
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
 400  continue
      
 420  fnam = ' '
      if (runmode.ne.'batchrun') then
         write (*,430)
 430     format(/' Enter file name : ',$)
      endif
      read(*,440) fnam
 440  format(a)
      inquire(file=fnam, exist=iexi)
c     
      if (iexi) then
c     
c     found the file...
c     
c     get turn on!
c
      etr = ephot(1)
 447  if (runmode.ne.'batchrun') then
      write(*, 443) etr
 443  format(/' Give the turn-on energy (Min.:',f6.2,' eV.) : ',$)
      endif
      read(*,*) yc
      if (yc.lt.0.0d0) goto 447
      turn = yc
c
c     get turn off!
c
      etr = ephot(infph)
 415  if (runmode.ne.'batchrun') then
      write(*, 425) etr
 425  format(/' Give high energy cut-off (Max.:',f9.2,'eV) : ',$)
      endif
c
      read(*, *) yc
      if (yc.lt.0.0d0) goto 415
      cut = yc
c
      open(luin, file=fnam, status='OLD') 
 43   format(a)
      do 40 j = 1, 3000
         read(luin, 43)  caract
         ki = index(caract,'PHOTON')
         if (ki.gt.0) then
            kii = index(caract,'SOURCE')
            kiii = index(caract,'FILE')
            if ((kiii.gt.kii).and.(kii.gt.ki)) goto 48
         end if
 40   continue
 198  if (runmode.ne.'batchrun') then
         write(*, 45)
 45      format(/'STRING : " PHOTON SOURCE FILE" NOT FOUND IN FILE')
      goto 100
c
      else
c
      if (runmode.ne.'batchrun') then
         write(*, 223) (fnam)
      endif
 223  format(/'FILE :',(a),' NOT FOUND')
      goto 400
c
      endif
c
 48   continue
c     
c     read and display any header..
c     
      read(luin, 43)  caract
      if (runmode.ne.'batchrun') then
         write(*, 79) (caract(2:80))
      endif
 79   format(a/)
      if (caract(1:1).eq.'%') goto 48
c     
c
c   Check for v2.0.0+ Field Verions ID (<0)
c
      read(luin, *) fv
      if (fv.ne.FieldVersion) then
c
c   Special case old files if old style compatibility PHOTDAT.old is used.
c
         if ((fv.eq.359).and.(FieldVersion.eq.0)) then
c  
               iinto = infph-1
               iinf = fv
               
               if (iinf .ne. iinto) then
               write(*, 1244) iinf, iinto
1244           format(/,' Invalid number of data points (',i6,').',
     &         ' Require :',i5,' points')
               endif
c
c     Multiply by 4 to convert 1/4pi sr in files to 
c     1/pi for Inu source field.
c     only read one column for old files
c
c
               do j = 1, iinto
                   read(luin, *) sv
                   den = 0.5d0*(ephot(iinto+1)-ephot(iinto))
                   if ((den.ge.turn).and.(den.lt.cut)) then
                       souvec(j) = souvec(j)+sv
                   endif
               enddo
c
         else 
         
            write(*, 243) fv, FieldVersion
 243  format(/' Invalid source vector version',i4,
     &        ' needed ',i4,'!!!!')
     
     
         endif
         
      else
c
c     Input matched the current field ID, assume the 
c     bins are OK and ignore energies, only reading
c     fluxes.  This avoids failing to compare energies when the 
c     input file may just not have enough decimals.
c
c     WARNING: All source vector with a given binning *must*
c     have a unique Field Version.  
c                                   Old = 0, only works with fv read as 359
c                                   Std = -2
c                                 Short = -3
c                              reserved = -4...-1024
c                          free for use   -1025... -max long int
c
      iinto = infph-1
      iinf = 0
      read(luin, *) iinf
c
c     error
c     
 241  if (iinf .ne. iinto) then
         if (runmode.ne.'batchrun') then
            write(*, 244) iinf, iinto
 244  format(/' Invalid number of data points (',i6,').',
     &     ' Require :',i5,' points')
         endif
      endif
c
c     Multiply by 4 to convert 1/4pi sr in files to 
c     1/pi for Inu source field.
c
      do j = 1, iinto
         read(luin, *) bv, sv
         den = 0.5d0*(ephot(iinto+1)-ephot(iinto))
         if ((den.ge.turn).and.(den.lt.cut)) then
            souvec(j) = souvec(j)+sv*4.d0
         endif
      enddo
c
      endif
c
      close(luin)
c
      srcfile = fnam
c
      goto 100
c         
      else
c     
 405     if (runmode.ne.'batchrun') then
            write(*,410)
         endif
 410     format(/,' File Not Found !! Try again or cancel (a/c): ',$)
         read(*,440) ilgg
         ilgg = ilgg(1:1)
         if (ilgg.eq.'a') ilgg = 'A'
         if (ilgg.eq.'c') ilgg = 'C'
         if (ilgg.eq.'A') goto 420
         if (ilgg.eq.'C') goto 100
         goto 405
c     
      endif 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ***WRITE A SOURCE FILE OF CURRENT FIELD
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
 300  continue
      
 320  fnam = ' '
      if (runmode.ne.'batchrun') then
         write (*,330)
 330     format(/' Enter new source file name : ',$)
      endif
      read(*,340) fnam
 340  format(a)
      inquire(file=fnam, exist=iexi)
c     
      if (iexi) then
c     
c     
 305  if (runmode.ne.'batchrun') then
         write(*,310)
      endif
 310  format(//,' File Already Exists!! Try again or cancel(a/c): ',$)
         read(*,340) ilgg
         ilgg = ilgg(1:1)
         if (ilgg.eq.'a') ilgg = 'A'
         if (ilgg.eq.'c') ilgg = 'C'
         if (ilgg.eq.'A') goto 320
         if (ilgg.eq.'C') goto 100
         goto 305
c
      else
c     
c     Valid file name...
c
      if (runmode.ne.'batchrun') then
         write (*,350)
 350     format(/' Give code/identifier for this source file : ',$)
      endif
      read(*,340) srcfile
      
      open(luop, file=fnam, status='NEW')
c
c
c   write header
c
10    format('% PHOTON SOURCE FILE'/
     &       '% Multi-component source file.'/
     &       '% UNITS: eV vs Fnu (ergs/cm2/s/Hz/sr)'/
     &       '% Fnu in 3D units (1/(4pi)) sr not '/
     &       '% 2D source (1/pi) units. '/
     &       '% bin energies are lower edge of bins.'/
     &       '% fluxes are average over bin.'/
     &       '% Field identification: ',a40)
      write(luop,10) srcfile
20    format('Field produced by MAPPINGS III v',a4)
      write(luop,20) theVersion
c
      srcfile = fnam
c
c
c   Field Version Number
c
      write(luop,*) FieldVersion
c
c   number of points
c
      write(luop,*) infph-1
c
c   write souvec.  Divide by 4 to convert Inu to 3D units
c
 355  format(1pg14.7,' ',1pg14.7)
      do j = 1,infph-1
          if(souvec(j)/4.d0.gt.epsilon)then
            write(luop,355) ephot(j),souvec(j)/4.d0
          else
            write(luop,355) ephot(j),epsilon
          endif
      enddo
c
      close(luop) 
c
      goto 100
c
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***STAR WITH BLACK BODY DISTRIBUTION
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 900  if (runmode.ne.'batchrun') then
         write(*, 905) 
 905  format(//' Give Blackbody temperature : ',$)
      endif
c
      srcfile = 'Black Body'
c
      read(*, *) teff
      if (teff.le.10) teff = 10**teff
      if (teff.lt.1.0d2) goto 900
c     
c     get turn on!
c
      etr = ephot(1)/ryd
 947  if (runmode.ne.'batchrun') then
      write(*, 943) etr
 943  format(/' Give the turn-on energy (Min.:',f6.2,' Ryd.) : ',$)
      endif
      read(*,*) yc
      if (yc.lt.0.0d0) goto 947
      turn = yc
c
c     get turn off!
c
      etr = ephot(infph)/ryd
 915  if (runmode.ne.'batchrun') then
         write(*, 925) etr
 925  format(/' Give high energy cut-off (Max.:',f6.2,'Ryd.) : ',$)
      endif
c
      read(*, *) yc
      if (yc.lt.0.0d0) yc = 0.d0
      cut = yc
c      
c
c     
 930  if (runmode.ne.'batchrun') then
      write(*,960)
 960  format(/' Give component scaling factor (>0): ',$)
      endif
      read(*,*) scale
      if (scale.le.0.d0) goto 930
c
c
c
      do i = 1, infph-1
          rnuh = 0.5d0*(ephot(i)+ephot(i+1))/(ryd)
          if ((rnuh.lt.cut).and.(rnuh.ge.turn)) then
             sv = fplank(teff,rnuh)
             souvec(i) = souvec(i)+sv*scale
             skipbin(i) = .FALSE.
          endif
       enddo
c
      goto 100
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***X-Ray Thermal Bremsstrahlung ; A exp (-h nu/kT)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
 800  if (runmode.ne.'batchrun') then
         write(*, 810)
      endif
      
 810  format(//' Give effective Temperature'/
     &     '(<=10 as log(K), <=1e5 as eV, >1e5 as K )',$)
      read(*, *) teff
c     
      srcfile = 'Bremsstrah'
c     
      if (teff.le.10.d0) then 
         teff = 10.d0**teff
         goto 830
      endif
c     
      if (teff.le.1.d5) then
         teff = teff*(ev/rkb)
         goto 830
      endif
c     
 830  if (runmode.ne.'batchrun') then
         write(*, 840)
 840     format(/' Give total flux, ( ergs/cm^2/s ) ',$ )
      endif
      read(*, *) yc
      
      if (yc.lt.0.0d0) goto 830
c     
c     get normalization constant
c     
      flu = (yc/(1.41d-6))*(877*plk)*dsqrt(rkb*teff/(1.d3*ev))
c
c     to get Inu...
c
      flu = flu/pi
c
c     get turn on!
c
 850  etr = ephot(1)
c
 847  if (runmode.ne.'batchrun') then
         write(*, 843) etr
 843  format(/' Give the turn-on energy  (Min.:',f6.2,' eV.) : ',$)
      endif
      read(*,*) yc
      if (yc.lt.0.0d0) goto 847
      turn = yc
c
c     get turn off!
c
      etr = ephot(infph)
c
 815  if (runmode.ne.'batchrun') then
         write(*, 825) etr
 825  format(/' Give high energy cut-off (Max.:',f9.2,'eV.) : ',$)
      endif
c
      read(*, *) yc
      if (yc.lt.0.0d0) goto 815
      cut = yc
c
c
c     
 855  if (runmode.ne.'batchrun') then
      write(*,856)
 856  format(/' Give component scaling factor (>0): ',$)
      endif
      read(*,*) scale
      if (scale.le.0.d0) goto 855
c     
      do i = 1, infph-1
         den = 0.5d0*(ephot(i)+ephot(i+1))
         wid = (ephot(i+1)-ephot(i))*evplk
         sv = flu*dexp(-(den*eV)/(teff*rkb))
         skipbin(i) = .TRUE.
         if (((den).ge.turn).and.((den).lt.cut)) then
            souvec(i) = souvec(i) + sv*scale
            skipbin(i) = .FALSE.
         endif
      enddo
c
      goto 100
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***NON-THERMAL SOURCE ; INTENSITY PROPORTIONAL TO NU**ALPHA
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
 1000 if (runmode.ne.'batchrun') then
         write(*, 1010)
      endif
 1010 format(//' Give spectral index, Alpha : ',$)
      read(*, *) alnth
c
      srcfile = 'Non-Thermal'
c
 1030 if (runmode.ne.'batchrun') then
         write(*, 1040)
 1040 format(/' Give zero point intensity I(E), E :'/
     &        '   I(E) in Inu 1/pi units (<0 as log),'/
     &        '   E in eV  (<0 as -Rydbergs)'/
     &        ':::: ',$ )
      endif
      read(*, *) yc,ze
c     
      if (yc.lt.0.d0) yc = 10.d0**yc
      flu = yc
c
      if (ze.lt.0.d0) ze = -Ryd*ze
c
c     get turn on!
c
 1042 etr = ephot(1)
 1047 if (runmode.ne.'batchrun') then
         write(*, 1043) etr
 1043 format(/' Give the turn-on energy  (Min.:',f6.2,' eV.) : ',$)
      endif
      read(*,*) yc
      if (yc.lt.0.0d0) goto 1047
      turn = yc
c
c     get turn off!
c
      etr = ephot(infph)
 1015 if (runmode.ne.'batchrun') then
         write(*, 1025) etr
 1025 format(/' Give high energy cut-off (Max.:',f9.2,' eV.) : ',$)
      endif
c
      read(*, *) yc
      if (yc.lt.0.0d0) goto 1015
      cut = yc
c     
c     
 1050 if (runmode.ne.'batchrun') then
      write(*,1060)
 1060 format(/' Give component scaling factor (>0): ',$)
      endif
      read(*,*) scale
      if (scale.le.0.d0) goto 1050
c
      do i = 1, infph-1
         den = 0.5d0*(ephot(i)+ephot(i+1))
         sv = flu*((den/ze) ** alnth)
         skipbin(i) = .TRUE.
         if ((den.ge.turn).and.(den.lt.cut)) then
            souvec(i) = souvec(i) + sv*scale
            skipbin(i) = .FALSE.
         endif
      enddo
c
      if (runmode.ne.'batchrun') then
      write(*, 1307) ze, flu
 1307 format(/' @ ',1pg10.3, ' eV : Intensity :',1pg10.3)
      endif
c
      goto 100
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***AGN QUASAR SOURCE ; INTENSITY Big Blue Bump + Powerlaw and cutoffs
c
c	L = k1 e^(-alpha1) exp(-e/e1) + k2 e^(-alpha2) exp(-e/e2) H(e - e1)
c
c	H is the heavyside operator  Spectrum in photons
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c

 1100 if (runmode.ne.'batchrun') then
         write(*, 1110)
 1110    format(//' Give warm component spectral index, Alpha 1: ',$)
      endif
      read(*, *) alpha1
      
      if (runmode.ne.'batchrun') then
         write(*, 1141)
 1141    format(/' Give warm component fraction, k1, (<0 as log): ',$ )
      endif
      read(*, *) fk1
      if (runmode.ne.'batchrun') then
         write(*, 1142)
 1142    format(/' Give warm component cutoff energy, e1,'
     &         ,' (eV, <0 as log): ',$ )
      endif
      read(*, *) fe1

      if (runmode.ne.'batchrun') then
         write(*, 1120)
 1120    format(//' Give hot component spectral index, Alpha 2 : ',$)
      endif
      read(*, *) alpha2
      if (runmode.ne.'batchrun') then
         write(*, 1140)
 1140    format(/' Give hot component fraction, k2, (<0 as log): ',$ )
      endif
      read(*, *) fk2
      if (runmode.ne.'batchrun') then
         write(*, 1143)
 1143    format(/,' Give hot component cutoff energy, e2,'
     &         ,' (eV, <0 as log): ',$ )
      endif
      read(*, *) fe2
c
      srcfile = 'JBH AGN'
c
c
c     get turn on!
c
      etr = ephot(1)
 1147 if (runmode.ne.'batchrun') then
         write(*, 1144) etr
 1144 format(/' Give the turn-on energy  (Min.:',f6.2,' eV) : ',$)
      endif
      read(*,*) yc
      if (yc.lt.0.0d0) goto 1047
      turn = yc
c
c     get turn off!
c
      etr = ephot(infph)
 1115 if (runmode.ne.'batchrun') then
         write(*, 1125) etr
 1125 format(/' Give high energy cut-off (Max.:',f9.2,'eV) : ',$)
      endif
c
      read(*, *) yc
      if (yc.lt.0.0d0) goto 1015
      cut = yc
c     
c     
 1150 if (runmode.ne.'batchrun') then
      write(*,1160)
 1160 format(/' Give component scaling factor (>0): ',$)
      endif
      read(*,*) scale
      if (scale.le.0.d0) goto 1050
c
      do i = 1, infph-1
         den = 0.5d0*(ephot(i)+ephot(i+1))
         
         sv = fk1 * (den ** alpha1) *dexp(-den/fe1)
         sv = sv  + (fk2 * (den ** alpha2) *dexp(-den/fe2))
     &              *fheavyside(ephot(i) - fe1)
c         sv = sv  + (fk2 * (den ** (alnth)) *dexp(-den/fe2))
         
         skipbin(i) = .TRUE.
         if ((den.ge.turn).and.(den.lt.cut)) then
            souvec(i) = souvec(i) + sv*scale
            skipbin(i) = .FALSE.
         endif
      enddo
c
      goto 100
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***AGN QUASAR SOURCE ; INTENSITY Big Blue Bump + Powerlaw and 
c                            exponential cutoffs
c
c	fv = k1*eng^(alpha1) exp(-eng/eng1) exp(-eng2/eng)  
c                + k2*eng(alpha2)exp(0.1eV/eng) exp(-eng2/10000eV)
c
c	Spectrum in photons
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 1300 if (runmode.ne.'batchrun') then
         write(*, 1310)
 1310    format(//,' Give index of low energy side of BBB, alpha1'
     &            ,'(~-0.5) : ',$)
         read(*, *) alpha1
         write(*, 1315)
 1315    format(//,' Give energy of UV cutoff of BBB, eng1'
     &            ,' (~42eV) : ',$)
         read(*, *) eng1
         write(*, 1320)
 1320    format(//,' Give energy of IR cutoff of BBB, eng2'
     &            ,' (~0.1eV) : ',$)
         read(*, *) eng2
         write(*, 1325)
 1325    format(//,' Give fraction of BBB, k1 : ',$)
         read(*, *) fk1
         write(*, 1330)
 1330    format(//,' Give index of X-ray power-law continuum, alpha2'
     &            ,'(~-0.85) : ',$)
         read(*, *) alpha2
         write(*, 1335)
 1335    format(//,' Give relative fraction of x-ray continuum, k1',/
     &            ,' (k2/k1 determines relative UV-Xray slope'
     &            ,' (~0.1)) : ',$)
         read(*, *) fk2
c
c
c     get turn on!
c
         etr = ephot(1)
 1339    write(*, 1340) etr
 1340    format(/' Give the turn-on energy  (Min.:',f6.2,' eV) : ',$)
         read(*,*) yc
         if (yc.lt.0.0d0) goto 1339
         turn = yc
c
c     get turn off!
c
         etr = ephot(infph)
 1344    write(*, 1345) etr
 1345    format(/' Give high energy cut-off (Max.:',f9.2,'eV) : ',$)
         read(*, *) yc
         if (yc.lt.0.0d0) goto 1344
         cut = yc
c     
c     
c     
     
 1350    write(*,1355)
 1355    format(/' Give Flux at 1 Ryd (>0): ',$)
         read(*,*) scale
         if (scale.le.0.d0) goto 1350
         sv=fk1*(ryd**alpha1)*dexp(-ryd/eng1)*dexp(-eng2/ryd)
         sv = sv  + (fk2 * (ryd**alpha2) *dexp(-ryd/1.0d5)
     &              *dexp(-eng2/ryd))
         scale=scale/sv

      do i = 1, infph-1
         den = 0.5d0*(ephot(i)+ephot(i+1))
         
         sv = fk1*(den**alpha1)*dexp(-den/eng1)*dexp(-eng2/den)
         sv = sv  + (fk2*(den**alpha2)*dexp(-den/1.0d5)
     &              *dexp(-eng2/den))
         
         skipbin(i) = .TRUE.
         if ((den.ge.turn).and.(den.lt.cut)) then
            souvec(i) = souvec(i) + sv*scale
            skipbin(i) = .FALSE.
         endif
      enddo
c
      endif
      goto 100



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***STELLAR ATMOSPHERE MODELS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

 200  if (runmode.ne.'batchrun') then
         write(*, 210)
      endif
 210  format(//' Give effective temperature : (26000<Teff<56000) : ',$)
      read(*, *) teff
c
      if ((teff.lt.2.6d4).or.(teff.gt.5.6d4)) goto 200
c
      
 217  if (runmode.ne.'batchrun') then
         write(*, 220) zgas
 220  format(/' Give the metalicity of the star ',
     &     /'(Zgas =',1pg9.3,'Zsun)   Zstar : ',$ )
      endif
c
      read(*, *) zstar
c
      if ((zstar.lt.0.d0).or.(zstar.gt.1.d1)) goto 217
      zmod = zstar/0.3216d0
c
c    ***CALCULATES RATIOS OF TEMPERATURE (AT THE EDGES)
c       TO EFFECTIVE TEMPERATURE
c
      gr = 4.0d0
      kstr = 'TNO4'
      if (gr .eq. 5.d0) kstr = 'TNO5'
c
      do i = 1, 16
         call fpol(kstr, tm1, i)
         tedge(i,1) = tno4(i,1)
         tedge(i,2) = tno4(i,2)
         tedge(i,3) = tm1
      enddo
c
      kstr = 'TMET'
c
      i = 1
      call fpol(kstr, pt0, i)
      do i = 2, 8
c
          j = tmet(i,1)
c
          dt1 = tedge(j,3)-pt0
c
          if (dt1.lt.1.0d-17) dt1 = 1.d-17
c
          call fpol(kstr, del, i)
          if (del.lt.0.0d0) del = 0.0d0
          dt0 = dt1+del
          dtot = dt0*((1.d0+((((dt0/dt1)**2)-1.d0)*zmod))
     &           **(-0.5d0))
          tedge(j,3) = pt0+dtot
c
      enddo
c
c     get turn on!
c
      etr = ephot(1)/ryd
 247  if (runmode.ne.'batchrun') then
         write(*, 245) etr
 245     format(/' Give the turn-on energy  (Min.:',f6.2,' Ryd.) : ',$)
      endif
      read(*,*) yc
      if (yc.lt.0.0d0) goto 247
      turn = yc
c
c     get turn off!
c
      etr = ephot(infph)/ryd
 215  if (runmode.ne.'batchrun') then
         write(*, 225) etr
 225  format(/' Give high energy cut-off (Max.:',f6.2,'Ryd.) : ',$)
      endif
c
      read(*, *) yc
      if (yc.lt.0.0d0) goto 215
      cut = yc
c
c
c     
 265  if (runmode.ne.'batchrun') then
      write(*,260)
 260  format(/' Give component scaling factor (>0): ',$)
      endif
      read(*,*) scale
      if (scale.le.0.d0) goto 265
c
c    ***INTERPOLATES BETWEEN EDGES TO GET TNU
c
      j = 1
      jmax = 16
      jg = 0
c
      do i = 1, infph-1
          den = 0.5d0*(ephot(i)+ephot(i+1))
          rnuh = den/ryd
c
  270     if (tedge(j+1,2).lt.den) j = j+1
          if (j.gt.jmax-1) goto 100
          if (tedge(j+1,2).lt.den) goto 270
          if (jg .eq. j) goto 275
c
          do n = 1, 2
             n1 = (j+n)-1
             xx(n) = dlog(tedge(n1,2))
             yy(n) = dlog(tedge(n1,3))
          enddo
c
          al(1) = (yy(2)-yy(1))/(xx(2)-xx(1))
          cc(1) = ((yy(1)*xx(2))-(yy(2)*xx(1)))/(xx(2)-xx(1))
          jg = j
c
  275     tnu = teff*dexp(cc(1)+(al(1)*dlog(den)))
          sv = fplank(tnu,rnuh)
c
          skipbin(i) = .TRUE.
          if ((rnuh.ge.turn).and.(rnuh.lt.cut)) then
             souvec(i) = souvec(i) + sv*scale
             skipbin(i) = .FALSE.
          endif
c
       enddo
c
      goto 100
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Local Interstellar Radiation Field 
c     	Mathis, Mezger, Panagia 1983
c
c     in photons cm-2 sr-1 eV-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 500    write(*,510)
 510	format(/' Give component scaling factor (>0): ',$)
	read(*,*) scale
      	if (scale.le.0.d0) goto 500
	i= 1
 	rnuh = 0.5d0*(ephot(i)+ephot(i+1))*eV/Ryd
 	do while (rnuh.lt.0.3706314) !lambda > 0.246 mum
	  sv = 1.0D-14*fplank(7.5d3,rnuh)+ 1.65D-13*fplank(4.0d3,rnuh)
     &	     + 4.0D-13*fplank(3.0d3,rnuh)
          souvec(i) = souvec(i)+sv
          skipbin(i) = .FALSE.
          i=i+1
	  rnuh=0.5d0*(ephot(i)+ephot(i+1))*eV/Ryd
	enddo
	  
	den=cls*1e4/(0.5d0*(ephot(i)+ephot(i+1))*eVplk) !den=wave in mum
	do while (den.lt.0.134) !0.246>L_mum>0.134
	  sv=7.115D-4*den**(0.3322) !factor of den**2/cls for J_Lam to J_nu
	  souvec(i)=souvec(i) + sv/(cls*1e4)  !cls in mum
          skipbin(i) = .FALSE.
	  i=i+1
	  den=cls*1e4/(0.5d0*(ephot(i)+ephot(i+1))*eVplk)
	enddo
	
	do while (den.lt.0.110)
	  sv=2.045D-2*den**2
	  souvec(i)=souvec(i) + sv/(cls*1e4)
          skipbin(i) = .FALSE.
	  i=i+1
	  den=cls*1e4/(0.5d0*(ephot(i)+ephot(i+1))*eVplk)	
	enddo
	
	do while (den.lt.0.0912) !0.110>L_mum>0.0912- upto 13.6 eV
	  sv= 38.57*den**(5.4172)!factor of den**2/cls for J_Lam to J_nu
	  souvec(i)=souvec(i) + sv/(cls*1e4)  !cls in mum
          skipbin(i) = .FALSE.
	  i=i+1
	  den=cls*1e4/(0.5d0*(ephot(i)+ephot(i+1))*eVplk)
	enddo

	do j=i, infph-1
	  skipbin(i)=.TRUE.
	enddo
	
	goto 100
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     'Standard' local ambient UV field.  Based on cubic 
c     fit  by Draine 1979. between 4 and 14 eV.
c
c     in photons cm-2 sr-1 eV-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c 500  srcfile = 'Soft UV'
c
c
c     get turn off!
c
c      etr = 1.d0
c 510  if (runmode.ne.'batchrun') then
c         write(*, 525) etr
c 525     format(/' Give high energy cut-off (Max.:',f6.2,'Ryd.) : ',$)
c      endif
c
c      read(*, *) yc
c      if (yc.lt.0.0d0) yc = 0.d0
c      if (yc.gt.1.0d0) yc = 1.d0
c      yc = yc*ryd
c     
c
c
c      scale = 1.d0
c
c 550  if (runmode.ne.'batchrun') then
c         write(*,560)
c 560     format(/' Give component scaling factor (>0): ',$)
c      endif
c      read(*,*) scale
c      if (scale.le.0.d0) goto 550
c
c     fit coeffs
c     
c      a1 = 1.658d6
c      a2 = -2.152d5
c      a3 = 6.919d3
c     
c     still loop through whole vector in case vector
c     changes one day.  only want bins below 1 Ryd.
c     
c      i=1
c      do 570 while (ephot(i).lt.6.d0) !5eV FUV min
c         i=i+1
c 570  continue
c      fuvmin=i
c      do i = fuvmin,infph
c
c
c         fev = 0.d0
c
c         if (ephot(i).lt.yc) then
c
c            den = 0.5d0*(ephot(i)+ephot(i+1))
c            wid = (ephot(i+1)-ephot(i))
c
c            fev = den*(a1+den*(a2+den*a3))
c
c     convert to ergs per eV, mult by energy per photon
c
c            fev = fev * den *ev
c
c     convert to per Hz over bin width (ev)
c     wid cancels, but this is clearer.
c
c     total bin energy
c
c            fev = fev * wid
c
c     divide by bin width in Hz
c
c            fev = fev/(wid*evplk)
c
c         endif
c
c         souvec(i) = souvec(i) + fev*scale
c
c      enddo
c
c      goto 100
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     End of sources, tidy up vector:
c
c    ***DETERMINATION OF HIGH ENERGY CUT-OFF
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 2000 continue
c
c    ***FIND actual CUT-OFF FREQUENCY
c
 2098 cut = 0.0d0
      do 2017 i = infph-1, 1, -1
      j = i
      if (souvec(i).gt.epsilon) goto 2013
 2017 continue
      goto 2011
 2013 cut = ephot(j+1)/ryd
 2011 continue
c
c    ***FIND actual TURN-ON FREQUENCY
c
 2073 turn = 0.0d0
      do 2087 i = 1, infph-1
      j = i
      if (souvec(i).gt.epsilon) goto 2083
 2087 continue
      goto 2071
 2083 turn = ephot(j)/ryd
 2071 continue
c
c    ***INTEGRATE NUMBER OF EMEGENT PHOTONS TO IONISE H,HE (=PI*JNU)
c    assumes 1/4pi units and souvec is 1/pi plane parallel units.
c
      q1 = 0.0d0
      q2 = 0.0d0
      q3 = 0.0d0
      q4 = 0.0d0
c
      call intvec(souvec, q1, q2, q3, q4)
c
      qhi   = q1*0.25d0
      qhei  = q2*0.25d0
      qheii = q3*0.25d0
      qht   = q4*0.25d0
c
      if (runmode.ne.'batchrun') then
      write(*, 2800) 
 2800 format(/' Mod',t7,'Temp.',t16,'Alpha',t22,'Turn-on',t30,'Cut-off'
     &,t38,'Zstar',t47,'FQHI',t56,'FQHeI',t66,'FQHeII')
      write(*, 2850) iso,teff,alnth,turn,cut,zstar,qhi, 
     &qhei, qheii
 2850 format(1h ,a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
      endif
c
c    ***PUT PHOTON SOURCE IN VECTOR SOUPHO IN /PHOTDAT/
c
      do 2851 i = 1, infph-1
      soupho(i) = souvec(i)
c      write (*,*) i,' skpb:',skipbin(i),' sv:',soupho(i)
 2851 continue
c
c     cosmic ray event rate
c
      crate = 0.d0
c
 3000 if (runmode.ne.'batchrun') then
         write(*,3010)
 3010    format(/' Include cosmic ray heating? (Y/N): ',$)
      endif
      read(*,440) ilgg
      ilgg = ilgg(1:1)
c     
      if (ilgg .eq. 'y') ilgg = 'Y'
      if (ilgg .eq. 'n') ilgg = 'N'
      if ((ilgg.ne.'Y').and.(ilgg.ne.'N')) goto 3000
c     
      if (ilgg.eq.'Y') then 
c
 3020    if (runmode.ne.'batchrun') then
         write(*,3030)
 3030    format(/' Give event rate per second (~1e-17): ',$)
         endif
         read(*,*) crate
         if (crate.lt.0.d0) goto 3020
      endif
c     
      return 
c
      end
