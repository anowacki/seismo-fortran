subroutine MakeMagGrid2D(rad, phi, theta, total, cilm, r0, a, f, lmax, interval, nlat, nlong, &
	north, south, east, west)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given the Schmidt normalized Spherical Harmonic coefficients cilm, this subroutine
!	will compute a 2D grid with equal latitude and longitude spacings for the radial, 
!	phi and theta directions, as well as the total field.  The field is calculated 
!	on a flattended ellipsoid with semi-major axis A and flattening F.
!	
!	Note that since this is NOT done using FFTs, this routine is relatively SLOW!
!
!	The value at grid(1,1) correspons to 90 degrees latitude and 0 degress
!	longitude. Longitude points are calculated from 0 to 360, with both points 
!	being included. If the optional parameters NORTH, SOUTH, EAST and WEST are specified, 
!	the upper-left and lower right coordinates of the output grid are (NORTH, WEST) and 
!	(SOUTH, EAST), respectively.
!
!	Calling Parameters:
!		IN
!			cilm:		Schmidt Normalized spherical harmonic coefficients with 
!					dimensions (2, lmax+1, lmax+1).
!			r0:		Reference radius of cilm coeficients.
!			a:		Semi-major axis of the reference ellipsoid.
!			f:		Flattening of reference radius: (R_eq - R_p)/R_eq, where R_eq = A.
!			lmax:		Maximum spherical harmonic degree of expansions to be performed.
!			interval:	Spacing of output grid in DEGREES.
!
!		OUT
!			rad:	 	Grid of radial magnetic field.
!			phi:		Grid of magnetic field in phi direction.
!			theta:		Grid of magnetic field in theta direction.
!			total:		Grid of magnetic intensity.
!			nlat:		Number of latitude points for the grid.
!			nlong:		Number of longitude points for the grid.
!		OPTIONAL (IN)
!			north		Maximum latitude to compute, in degrees.
!			south		Minimum latitude to compute, in degrees.
!			east		Maximum longitude to compute, in degrees.
!			west		Minimum latitude to compute, in degrees.
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, then this file will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!		2. 	Latitude is geocentric latitude.
!
!		Dependencies:	PlmSchmidt_d1, PlSchmidt, PlmIndex
!
!	Written by Mark Wieczorek (2003)
!
!	August 2005 - 	sines and cosines are now precomputed, leading to a savings in time
!			by a factor of up to 15 (if the code is compiled with normal
!			optimizations.
!	May 31, 2006 -  Modified so that the Condon-Shortley phase is never used, regardless of 
!			what the default value of CSPHASE_DEFAULT is.
!	June 21, 2007 - Modified routine such that the field is calculated on a flattened
!			ellisoid, using exact equations for the ellipsoid. Previously, only
!			the degree-2 Legendre expansion was utilized. It is now necessary to
!		 	input the semi-major axis of the ellispoid instead of the mean ellipsoid
!			radius.
!	August 19, 2007	Added option to calculate subdomains of a raster grid.
!
!	Copyright (c) 2005-2007, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PlmSchmidt_d1, PlSchmidt, PlmIndex

	implicit none
	
	real*8, intent(in) :: 	cilm(:,:,:), interval, r0, a, f
	real*8, intent(out) :: 	rad(:,:), phi(:,:), theta(:,:), total(:,:)
	integer, intent(in) :: 	lmax
	integer, intent(out) :: nlat, nlong
	real*8, intent(in), optional :: north, south, east, west
	integer :: 		l, m, j, k, index, l1, m1, lmax_comp, temp, astat(4)
	real*8 :: 		pi, latmax, latmin, longmin, longmax, lat, &
				longitude, x, intervalrad, r_ex, coef
	real*8, allocatable ::	pl(:), dpl(:), cosm(:,:), sinm(:,:)
		
    	
    	temp = 0
	if (present(north)) temp = temp + 1
	if (present(south)) temp = temp + 1
	if (present(east)) temp = temp + 1
	if (present(west)) temp = temp + 1
	
	if (temp /=0 .and. temp /=4) then
		print*, "Error --- MakeMagGrid2d"
		print*, "The optional parameters NORTH, SOUTH, EAST, and WEST must all be specified", &
			present(north), present(south), present(east), present(west)
		stop
	endif
	
	if (temp == 4) then
		latmax = north
		latmin = south
		longmin = west
		longmax = east
		
		if (latmax < latmin) then
			print*, "Error --- MakeMagGrid2d"
			print*, "NORTH must be larger than SOUTH."
			print*, "NORTH = ", latmax
			print*, "SOUTH = ", latmin
			stop
		endif
		
		if (latmax > 90.0d0 .or. latmin < -90.0d0) then
			print*, "Error --- MakeMagGrid2d"
			print*, "NORTH and SOUTH must lie between 90 and -90."
			print*, "NORTH = ", latmax
			print*, "SOUTH = ", latmin
			stop
		endif
		
		if (longmin > longmax) longmax = longmax + 360.0d0
		
	else
		latmax = 90.0d0
		latmin = -90.0d0
		longmin = 0.0d0
		longmax = 360.d0
	endif
	
	nlat = (latmax - latmin) / interval + 1
	nlong = (longmax - longmin) / interval + 1
    	
    	
 	if (size(rad(:,1)) < nlat .or. size(rad(1,:)) < nlong ) then
		print*, "Error --- MakeMagGrid2D"
		print*, "RAD must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
		print*, "INTERVAL = ", interval
		print*, "LATMAX = ", latmax
		print*, "LATMIN = ", latmin
		print*, "LONGMIN = ", longmin
		print*, "LONGMAX = ", longmax
		print*, "Input array is dimensioned ", size(rad(:,1)), size(rad(1,:))
		stop
 	elseif (size(phi(:,1)) < nlat .or. size(phi(1,:)) < nlong ) then
		print*, "Error --- MakeMagGrid2D"
		print*, "PHI must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
		print*, "INTERVAL = ", interval
		print*, "LATMAX = ", latmax
		print*, "LATMIN = ", latmin
		print*, "LONGMIN = ", longmin
		print*, "LONGMAX = ", longmax
		print*, "Input array is dimensioned ", size(phi(:,1)), size(phi(1,:))
		stop
 	elseif (size(theta(:,1)) < nlat .or. size(theta(1,:)) < nlong ) then
		print*, "Error --- MakeMagGrid2D"
		print*, "THETA must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
		print*, "INTERVAL = ", interval
		print*, "LATMAX = ", latmax
		print*, "LATMIN = ", latmin
		print*, "LONGMIN = ", longmin
		print*, "LONGMAX = ", longmax
		print*, "Input array is dimensioned ", size(theta(:,1)), size(theta(1,:))
		stop
 	elseif (size(total(:,1)) < nlat .or. size(total(1,:)) < nlong ) then
		print*, "Error --- MakeMagGrid2D"
		print*, "TOTAL must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
		print*, "INTERVAL = ", interval
		print*, "LATMAX = ", latmax
		print*, "LATMIN = ", latmin
		print*, "LONGMIN = ", longmin
		print*, "LONGMAX = ", longmax
		print*, "Input array is dimensioned ", size(total(:,1)), size(total(1,:))
		stop
	elseif (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- MakeMagGrid2D"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	allocate(pl((lmax+1)*(lmax+2)/2), stat = astat(1))
	allocate(cosm(int(360./interval + 1), lmax+1), stat = astat(2))
	allocate(sinm(int(360./interval + 1), lmax+1), stat = astat(3))
	allocate(dpl((lmax+1)*(lmax+2)/2), stat = astat(4))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /=0 .or. astat(4) /= 0) then
		print*, "Error --- MakeMagGrid2D"
		print*, "Problem allocating arrays PL, COSM, SINM, and DPL", astat(1), astat(2), astat(3), astat(4)
		stop
	endif
	
	pi = acos(-1.0d0)
	rad = 0.0d0
	phi = 0.0d0
	theta = 0.0d0
	total = 0.0d0

	intervalrad = interval*pi/180.0d0
	
	lmax_comp = min(lmax, size(cilm(1,1,:))-1)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Loop over latitude. Precompute sin and cosine terms.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do k=1, nlong
	
		longitude = longmin*pi/180.0d0 + dble(k-1)*intervalrad
		do m=0, lmax
			cosm(k, m+1) = cos(m*longitude)
			sinm(k, m+1) = sin(m*longitude)
		enddo
		
	enddo

	do j=1, nlat
	
		lat = latmax - dble(j-1)*interval
		x = sin(lat*pi/180.0d0)
		
		if (lat == 90.0d0 .or. lat == -90.0d0) then
		
			call PlSchmidt(pl, lmax_comp, x)
	
			r_ex = a * (1.0d0 - f)	! Reference ellipsoid to calculate field			

			coef = (r0/r_ex)**2
			
			do l = 1, lmax_comp, 1
				l1 = l+1
				coef = coef * (r0/r_ex)
				rad(j,1) = rad(j,1) + dble(l+1) * cilm(1,l1,1) * pl(l1) * coef
			enddo
	
			phi(j,1:nlong) = 0.0d0 
			rad(j,2:nlong) = rad(j,1)
			theta(j,1:nlong) = 0.0d0
			
		else
		
			r_ex = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
				(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
			r_ex = sqrt(r_ex)

			call PlmSchmidt_d1(pl, dpl, lmax_comp, x, csphase = 1)
			dpl = - dpl * cos(lat*pi/180.0d0)	! convert derivative with respect to z to
								! a derivative with respect to theta.
			do k = 1, nlong
					
				coef = (r0/r_ex)**2
					
				do l = 1, lmax_comp, 1
					l1 = l+1
					coef = coef * (r0/r_ex)
					do m = 0, l
						index = PlmIndex(l,m)
						m1 = m+1
						rad(j,k) = rad(j,k) + dble(l+1) * (cilm(1,l1,m1)*cosm(k,m1) + &
							cilm(2,l1,m1)*sinm(k, m1) ) * pl(index) * coef
						phi(j,k) = phi(j,k) + dble(m)*(cilm(1,l1,m1)*sinm(k,m+1) - &
							cilm(2,l1,m1)*cosm(k, m+1) ) * pl(index) * coef &
							/ cos(lat*pi/180.0d0)
						theta(j,k) = theta(j,k) - (cilm(1,l1,m1)*cosm(k,m1) + &
							cilm(2,l1,m1)*sinm(k,m1) ) * dpl(index)  &
							* coef
						
					enddo
				enddo
			
			enddo
			
		endif
		
	enddo
	

	total(1:nlat, 1:nlong) = sqrt( rad(1:nlat,1:nlong)**2 + phi(1:nlat,1:nlong)**2 + theta(1:nlat,1:nlong)**2)
	
	! deallocate memory
	call PlmSchmidt_d1(pl, dpl, -1, x, csphase = 1)
	deallocate(pl)
	deallocate(cosm)
	deallocate(sinm)
	deallocate(dpl)

end subroutine MakeMagGrid2D

