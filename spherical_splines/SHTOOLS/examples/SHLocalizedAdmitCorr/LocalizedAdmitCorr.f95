Program LocalizedAdmitCorr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program demonstrates how to calucate localized admittance and
!	correlation functions
!
!	Written by Mark Wieczorek (April, 2005)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	use SHTOOLS
	use PlanetsConstants
	
	implicit none
	
	integer, parameter ::	degmax = 360
	character*80 :: topography_file, potential_file, outfile
	real*8 ::	topo(2, degmax+1, degmax+1), header(8), mpr, pot(2, degmax+1, degmax+1), r0_pot, &
			gm, mass, lat, lon, pi, theta0, alpha, sn, r
	real*8, allocatable ::	grav(:,:,:), admit(:), corr(:), admit_error(:), corr_error(:), tapers(:,:), &
				eigenvalues(:)
	integer ::	lmax_topo, lmax_pot, lmax, option1, l, lwin, lmaxwin, astat(5), K
	integer, allocatable ::	taper_order(:)
	
	pi = acos(-1.0d0)

	topography_file = "../ExampleDataFiles/Mars2000.shape"
	potential_file = "../ExampleDataFiles/jgm85h02.sh"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Read topography and gravity fields, and from the header information determine the 
	!	the mean planetary radius and mass. Verify that these values are correct as	
	!	it is possible the input files could be using a different set of base	
	!	units (i.e., km instead of meters).	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	print*, "Reading = ", topography_file	
	call SHRead(topography_file, topo, lmax_topo, header=header(1:8))	
	
	print*, "Lmax of topography file = ", lmax_topo
	
	mpr = topo(1,1,1)	
	print*, "Mean planetary radius (km) = ", mpr/1.d3	
		
	print*, "Reading = ", potential_file	
	call SHRead(potential_file, pot, lmax_pot, header=header(1:2))	
	
	print*, "Lmax of potential file = ", lmax_pot	
	
	r0_pot = header(1)	
	gm = header(2)	
	mass = gm/Grav_constant	
	
	print*, "Reference radius of potential coefficients (km) = ", r0_pot/1.d3	
	
	print*, "Mass of planet (kg) = ", mass		
	
	lmax = min(lmax_topo, lmax_pot)	
	print*, "Maximum spherical harmonic degree to be used in calculations = ", lmax		
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!	
	!	Get localization parameters from user.	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
	print*, "Number of tapers to use > "
	read(*,*) K

	print*, "Latitude of feature of interest (degrees) > "	
	read(*,*) lat		
	print*, "Longitude (degrees) > "	
	read(*,*) lon			
	print*, "Angular radius of localization window (degrees) > "	
	read(*,*) theta0	
	theta0 = theta0*pi/180.0d0
			
	print*, "Create localization window using"	
	print*, "(1) Desired concentration factor (alpha^2)"	
	print*, "(2) Desired (approximate) Shannon number; (Lwin+1) Theta0 / pi)"	
	print*, "(3) Desired Spectral bandwidth (Lwin)"	
	read(*,*) option1		
	
	if (option1==1) then		
		print*, "Input desired concentration factor of the Kth taper > "		
		read(*,*) alpha		
		lwin = SHFindLWin(theta0, 0, alpha, K)		
		print*, "Corresponding spherical harmonic bandwidth = ", lwin		
		print*, "Corresponding approximate Shannon number = ", (lwin+1)*theta0/pi	
	elseif (option1==2) then		
		print*, "Input Shannon number > "		
		read(*,*) sn		
		lwin = nint(sn*pi/theta0) - 1		
		print*, "Corresponding spherical harmonic bandwidth = ", lwin	
	elseif (option1==3) then		
		print*, "Input Lwin > "		
		read(*,*) lwin		
		print*, "Corresponding approximate Shannon number = ", (lwin+1)*theta0/pi	
	else		
		stop	
	endif		
	
	allocate(tapers(lwin+1, (lwin+1)**2), stat = astat(1))
	allocate(taper_order((lwin+1)**2), stat = astat(2))
	allocate(eigenvalues((lwin+1)**2), stat = astat(3))
	
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 ) then
		print*, "Problem allocatig arrays for tapers, taper_order, or eigenvalues", astat(1), astat(2), astat(3)
		stop
	endif
	
	call SHReturnTapers(theta0, lwin, tapers, eigenvalues, taper_order)
	
	print*, "Concentration factor of first taper = ", eigenvalues(1)
	print*, "Taper order = ", taper_order(1)
	
	lmaxwin = lmax+lwin		
	
	r = r0_pot
	
	print*, "Gravity field evaluated at R = (km) ", r0_pot/1.d3
		
	print*, "Name of output admittance and corrence file > "	
	read(*,*) outfile	
	open(12, file=outfile)		
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!	
	!	Allocate memory for arrays based on Lmax and Lwin. 	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	


	allocate(grav(2,lmax+1,lmax+1), stat=astat(1))	
	allocate(admit(lmax+lwin+1), stat=astat(2))	
	allocate(corr(lmax+lwin+1), stat=astat(3))	
	allocate(admit_error(lmax+lwin+1), stat=astat(4))
	allocate(corr_error(lmax+lwin+1), stat=astat(5))
	
	grav = 0.0d0
	admit = 0.0d0
	corr = 0.0d0
	admit_error = 0.0d0
	corr_error = 0.0d0		

	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. astat(4) /= 0 &
		.or. astat(5) /=0) then		
		print*, "Problem allocating memory"		
		stop	
	endif		
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!
	!	Create gravity coefficients in units of mGals. Convert Topography	
	!	coefficients to km.	
	!	
	!	Cilm(gravity) = 1.d5*G*M*Cilm(potential)*(l+1)*(r0_pot/r)**(l+2)/r0_pot**2	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

	do l=0, lmax		
		grav(1:2,l+1,1:l+1) = pot(1:2,l+1,1:l+1) * dble(l+1) * (r0_pot/r)**(l+2)	
	enddo		
	
	grav = grav * 1.0d5 * gm / (r0_pot**2)	
	topo = topo / 1.0d3	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!	
	!	Compute localized admittances, and write data to output file.	
	!	Units for the admittance are mgals/km.	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

	call SHLocalizedAdmitCorr(tapers, taper_order, lwin, lat, lon, grav, topo, lmax, admit, &		
		corr, K, admit_error=admit_error, corr_error=corr_error, mtdef=1)		
		
	if (K==1) then	
		do l = 0, lmax-lwin		
			write(12,*) l, admit(l+1), corr(l+1)
		enddo	
	else
		do l = 0, lmax - lwin		
			write(12,*) l, admit(l+1), admit_error(l+1), corr(l+1), corr_error(l+1)
		enddo	
	endif

	close(12)		
	
	deallocate(grav)	
	deallocate(admit)	
	deallocate(admit_error)	
	deallocate(corr)
	deallocate(corr_error)
	deallocate(tapers)
	deallocate(taper_order)
	deallocate(eigenvalues)
		
end program LocalizedAdmitCorr
