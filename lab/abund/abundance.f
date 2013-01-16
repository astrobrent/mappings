	program abundance
	implicit none
	integer num, atom(16), i
	real*8	abund(16), fracsolar, oxygen
	character*20 title, filename

 101	format('What Multiple of solar? :')
 	write(*,101)
 	read(*,*) fracsolar
 102	format(f4.2,'xSolar chosen')	
	write(*,102) fracsolar
 103	format('asplundx',f4.2,'.abn')
	write(filename,103) fracsolar
	open(1,file='asplund.abn')
	read(1,*) title
	do while (title(1:1).eq.'%')
	  read(1,*) title
	enddo
	write(*,*) title
	read(1,*) num
	
	do i=1, num
	  read(1,*) atom(i), abund(i)
	  if (atom(i).eq.8) oxygen=abund(i)+log10(fracsolar)
	enddo
	
	close(1)
	
	do i=2, num  !skip H
  	  if (atom(i).eq.2) then
  	    abund(i)=log10(0.0737+0.0293*fracsolar)
	  else if (atom(i).eq.7) then
	    abund(i)=log10(10**(oxygen)*(10**(-1.6)+10**(2.33+oxygen)))
	  else
	    abund(i)=abund(i)+log10(fracsolar)
	  endif
	enddo
	
	open(1,file=filename)
 105	format('%',/
     &		'%	SOLAR ABUNDANCES:ref Grevesse & Sauval 1998 Space Science Revies 85, 161'/
     &		'%	meteoritic abundances with C,N,O, Si, Fe from Asplund papers (2000,01, 02)'/
     &		'%	All approx +/- 0.05'/
     &		'%'/
     &		'%	FORMAT'/
     &		'%	Title'/
     &		'%	# entries'/
     &		'%	Z  Abundance'	/
     &		'%	.'/
     &		'%	.'/
     &		'%	.'/
     &		'%'/
     &		'%'/
     &		'%	Note: negative abundances are logarithmic relaive to H = 0.0'/
     &		'%'/
     &		'%	Elements not included in the data file ATDAT will not be read'/
     &		'%'/
     &		'%	Elements not present in this file will have no effect on the' /
     &		'%	default abundances.'/
     &		'%'/
     &		'%	The Header lines are not needed but they are useful'/
     &		'%')
     	write(1,105)
 106	format(f4.2,'xSolar abundance from Grevesse & Sauval, and Asplund (2000)')
     	write(1,106) fracsolar
	write(1,'(i2)') num
 107	format(i2,2x,f6.3)
      	do i=1, num
      	  write(1,107) atom(i), abund(i)
      	enddo
	close(1)
	end