program MakeMag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read in a set of magnetic spherical harmonic coefficients
!	and will expand these into gridded raster ascii files corresponding
!	to the phi, theta, radial and total field. The field is calculated
!	on a speroid with radius r and flattening f.
!
!	The included spherical harmonic file FSU_mars90.sh is the martian magnetic
!	field model of Cain et al., 2003.
!
!	Dependencies:	SHRead
!			MakeMagGrid2D
!
!	Written by Mark Wieczorek, February 2004
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS

	implicit none
	integer, parameter ::	maxdeg = 100, gridmax = 1
	character*80 ::		infile, radf, thetaf, phif, totalf
	real*8 ::		glm(2,maxdeg+1, maxdeg+1), header(4), interval, r0, r, &
				rad(181*gridmax, 361*gridmax), phi(181*gridmax, 361*gridmax), &
				theta(181*gridmax, 361*gridmax), total(181*gridmax, 361*gridmax), &
				f, mpr, z, timein, timeout
	integer ::		lmax, nlong, nlat, i, j
	
	
	f = 1.0d0/169.864881d0	! Mars flattening = (R_eq - R_p)/R_eq
	mpr = 3389.508d3	! Mean radius of mars
	z = 145.d3		! mean altitude to calculate field at.
	
	infile = "../ExampleDataFiles/FSU_mars90.sh"
	
	call SHRead(infile, glm, lmax, header=header(1:4), skip=1)
	r0 = header(1)*1.d3
	print*, "R0 (km) = ", r0/1.d3
	print*, "Lmax of file = ", lmax
	
	r = mpr + z
	print*, "R (km) = ", r/1.d3
		
	print*, "Interval for gridded data (degrees) > "
	read(*,*) interval
	
	if (interval < gridmax) then
		print*, gridmax, interval
		stop
	endif
	
	radf = "radial_145f.dat"
	thetaf = "theta_145f.dat"
	phif = "phi_145f.dat"
	totalf = "total_145f.dat"
	
	open(12,file=radf)
	open(13,file=phif)
	open(14,file=thetaf)
	open(15,file=totalf)
	
	call cpu_time(timein)
	
	call MakeMagGrid2D(rad, phi, theta, total, glm, r0, r, f, lmax, interval, nlat, nlong)
	
	call cpu_time(timeout)
	
	print*, "Elapsed time (sec) = ", timeout-timein
	
	print*, "Maximum and minimum intensity (nT) = ", maxval(total(1:nlat,1:nlong)), minval(total(1:nlat,1:nlong))
	
	print*, nlat, nlong
	
	! write(12,*) nlat, nlong
	! write(13,*) nlat, nlong
	! write(14,*) nlat, nlong
	! write(15,*) nlat, nlong
	
	do i=1, nlat
		do j=1, nlong
			write(12,*) rad(i,j)
			write(13,*) phi(i,j)
			write(14,*) theta(i,j)
			write(15,*) total(i,j)
		enddo
	enddo
	
	close(12)
	close(13)
	close(14)
	close(15)
	
end program