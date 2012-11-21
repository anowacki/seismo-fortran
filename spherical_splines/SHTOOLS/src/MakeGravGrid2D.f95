subroutine MakeGravGrid2D(rad, cilm, lmax, r0, a, f, gm, gravpot, interval, nlat, nlong, &
	theta, phi, total, omega, north, south, east, west, normal_gravity)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given the gravitational spherical harmonic coefficients CILM, this subroutine
!	will compute a 2D grid with equal latitude and longitude spacings of either the 
!	radial gravity or gravitational potential. These quantities are calculated on a 
!	flattened ellipsoid, with flattening F and semimajor axis A. In order to calculate 
!	the entire gravitational acceleration, it is necessary that the degree-0 term be set equal to 1. 
!	Radial gravity is assumed to be positive when directed UPWARDS. If the optional parameter 
!	OMEGA is specified, the gravitational acceleration will be calculated in the reference frame 
!	of a rotating body.
!
!	Note that (1) since this routine does not use FFTs, it is relatively slow, (2) the calculations
!	are only strictly exact when all radii on the ellispoid are less than the maximum radius
!	of the planet, (3) the radial gravity is in the radial direction, not normal to the ellipsoid,
!	amd (4) lattiude is geocentric latitude.  
!	
!	The value at rad(1,1) corresponds to 90 degrees latitude and 0 degrees longitude, and 
!	the longitude spacing goes from 0 to 360, with both points being calculated. If the optional
!	parameters NORTH, SOUTH, EAST and WEST are specified, the upper-left and lower right coordinates
!	of the output grid are (NORTH, WEST) and (SOUTH, EAST), respectively. The output units are in
!	m/s**2.
!
!	Calling Parameters:
!		IN
!			cilm:		Gravitational spherical harmonic coefficients.
!			lmax:		Maximum degree of expansions to be performed.
!			r0:		Reference radius of potential coefficients.
!			a:		The semimajor axis of the flattened ellipsoid (i.e., the 
!					mean equatorial radius).
!			f:		Flattening of the planet.
!			gravpot:	If "U" create potential. If "G" creat radial gravity.
!			interval:	Spacing of output grid in DEGREES.
!		IN, OPTIONAL
!			omega:		Angular rotation rate of the planet.
!			north		Maximum latitude to compute, in degrees.
!			south		Minimum latitude to compute, in degrees.
!			east		Maximum longitude to compute, in degrees.
!			west		Minimum latitude to compute, in degrees.
!		OUT
!			rad:		Gridded expansion of the radial component of the 
!					gravitational field or potential.
!			nlat:		Number of latitude points for the grid.
!			nlong:		Number of longitude points for the grid.
!
!		OUT, OPTIONAL
!			theta		Gridded expansaion of the theta component of the 
!					gravitational field.
!			phi		Gridded expansaion of the phi component of the 
!					gravitational field.
!			total		Gridded expansaion of the the magnitude of the 
!					gravitational field.
!			normal_gravity	If 1, the magnitude of the normal gravity on the ellipsoid 
!					will be removed from the magnitude of the total gravity vector.
!					This is the "gravity disturbance." 
!
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, then this file will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!		2. 	Latitude is geocentric latitude.
!
!	Dependencies:	PlmBar, PlBar
!
!	Written by Mark Wieczorek (2007)
!
!	June 21, 2007 -	Modified routine such that the field is calculated on a flattened
!			ellisoid, using exact equations for the ellipsoid. Previously, only
!			the degree-2 Legendre expansion was utilized. It is now necessary to
!			input the semi-major axis of the ellispoid instead of the mean ellipsoid
!			radius.
!	April 20, 2008	Added options to calculate the theta and phi components of the gravity field,
!			as well as the magnitue. THE RADIAL GRAVITY IS NOW DEFINED POSITIVE WHEN
!			DIRECTED UPWARDS. OUTPUT UNITS ARE SI (I.E., M/S**2).
!	August 20, 2009	Add the optional argument NORMAL_GRAVITY. If 1, the magnitude of the normal 
!			gravity (on the ellipsoid) will be removed from the magnitude of the total
!			gravity.
!
!	Copyright (c) 2007-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PlmBar, PlBar, PlmBar_d1, PlBar_d1

	implicit none
	
	interface
		real*8 function NormalGravity(geocentric_lat, GM, omega, a, b)
			real*8, intent(in) ::	geocentric_lat, gm, omega, a, b
		end function NormalGravity
	end interface

	real*8, intent(in) :: 		cilm(:,:,:), interval, f, r0, a, gm
	real*8, intent(out) :: 		rad(:,:)
	integer, intent(in) :: 		lmax
	integer, intent(out) :: 	nlat, nlong
	character*1, intent(in) ::	gravpot
	real*8, intent(in), optional :: omega, north, south, east, west
	real*8, intent(out), optional :: theta(:,:), phi(:,:), total(:,:)
	integer, intent(in), optional ::	normal_gravity
	integer :: 			l, m, j, k, index, l1, m1, lmax_comp, temp, astat(4)
	real*8 :: 			pi, latmax, latmin, longmin, longmax, lat, longitude, &
					x, intervalrad, r_ex, coef, b
	real*8, allocatable ::		pl(:), dpl(:), cosm(:, :), sinm(:, :), cilm2(:,:,:)
	
	
	temp = 0
	if (present(north)) temp = temp + 1
	if (present(south)) temp = temp + 1
	if (present(east)) temp = temp + 1
	if (present(west)) temp = temp + 1
	
	if (temp /=0 .and. temp /=4) then
		print*, "Error --- MakeGravGrid2d"
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
			print*, "Error --- MakeGravGrid2d"
			print*, "NORTH must be larger than SOUTH."
			print*, "NORTH = ", latmax
			print*, "SOUTH = ", latmin
			stop
		endif
		
		if (latmax > 90.0d0 .or. latmin < -90.0d0) then
			print*, "Error --- MakeGravGrid2d"
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
	
	if (present(normal_gravity)) then
		if (normal_gravity /= 0 .and. normal_gravity /= 1) then
			print*, "Error --- MakeGravGrid2d"
			print*, "NORMAL_GRAVITY must be either 1 (remove normal gravity)"
			print*, "or 0 (do not remove normal gravity)."
			print*, "Input value of NORMAL_GRAVITY is ", normal_gravity
			stop
		elseif (.not.present(total) .and. normal_gravity == 1) then
			print*, "Error --- MakeGravGrid2d"
			print*, "TOTAL must be specified when removing the normal gravity."
			stop
		elseif (.not.present(omega) .and. normal_gravity == 1) then
			print*, "Error --- MakeGravGrid2d"
			print*, "OMEGA must be specified when removing the normal gravity."
			stop
		endif
	endif
	
	temp = 0
	if (present(theta)) temp = temp + 1
	if (present(phi)) temp = temp + 1
	if (present(total)) temp = temp + 1
	
	if (temp /= 0 .and. temp /= 3) then
		print*, "Error --- MakeGravGrid2d"
		print*, "The optional parameters THETA, PHI, and TOTAL must all be specified", &
			present(theta), present(phi), present(total)
		stop
	endif

 	if (gravpot /= "U" .and. gravpot /= "G" .and. gravpot /= "u" .and. gravpot /= "g") then
 		print*, "Error --- MakeGravGrid2D"
 		print*, "GRAVPOT must be equal to either U (potential) or G (gravity)."
 		print*, "Input value is ", gravpot
 		stop
 	endif
 	
 	if (gravpot == "U" .or. gravpot == "u") then
 		if (present(phi) .or. present(theta) .or. present(total) .or. present(normal_gravity)) then
 			print*, "Error --- MakeGravGrid2d"
 			print*, "PHI, THETA, TOTAL and NORMAL_GRAVITY can not be present when calculating the gravitational potential."
 			stop
 		endif
 	endif
    	
 	if (size(rad(:,1)) < nlat .or. size(rad(1,:)) < nlong ) then
		print*, "Error --- MakeGravGrid2D"
		print*, "RAD must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
		print*, "INTERVAL = ", interval
		print*, "LATMAX = ", latmax
		print*, "LATMIN = ", latmin
		print*, "LONGMIN = ", longmin
		print*, "LONGMAX = ", longmax
		print*, "Input array is dimensioned ", size(rad(:,1)), size(rad(1,:))
		stop
	elseif (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- MakeGravGrid2D"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if (present(phi) .and. present(theta) .and. present(total)) then
		if (size(theta(:,1)) < nlat .or. size(theta(1,:)) < nlong &
			.or. size(phi(:,1)) < nlat .or. size(phi(1,:)) < nlong &
			.or. size(total(:,1)) < nlat .or. size(total(1,:)) < nlong) then
			print*, "Error --- MakeGravGrid2D"
			print*, "THETA, PHI, and TOTAL must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
			print*, "INTERVAL = ", interval
			print*, "LATMAX = ", latmax
			print*, "LATMIN = ", latmin
			print*, "LONGMIN = ", longmin
			print*, "LONGMAX = ", longmax
			print*, "Input arrays are dimensioned "
			print*, "THETA: ", size(theta(:,1)), size(theta(1,:))
			print*, "PHI: ", size(phi(:,1)), size(phi(1,:))
			print*, "TOTAL: ", size(total(:,1)), size(total(1,:))
			stop
		endif
	endif
	
	if (cilm(1,1,1) /= 1.0d0 .and. f /= 0.0d0) then
		print*, "Warning --- MakeGravGrid2D"
		print*, "The degree-0 term of the spherical harmonic coefficients is not equal to 1."
		print*, "The variation in gravity due to variations in radius of the flattened"
		print*, "ellipsoid will not be taken into accout."
		print*, "C00 = ", cilm(1,1,1)
		print*, "F = ", f
	endif
	
	allocate(pl((lmax+1)*(lmax+2)/2), stat = astat(1))
	allocate(cosm(int(360./interval + 1), lmax+1), stat = astat(2))
	allocate(sinm(int(360./interval + 1), lmax+1), stat = astat(3))
	allocate(cilm2(2,lmax+1,lmax+1), stat = astat(4))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /=0 .or. astat(4) /= 0) then
		print*, "Error --- MakeGravGrid2D"
		print*, "Problem allocating arrays PL, COSM, SINM, and CILM2", astat(1), astat(2), astat(3), astat(4)
		stop
	endif
	
	if (present(total)) then
		allocate(dpl((lmax+1)*(lmax+2)/2), stat = astat(1))
		if (astat(1) /= 0) then
			print*, "Error --- MakeGravGrid2D"
			print*, "Problem allocating arrays DPL", astat(1)
			stop
		endif
		theta = 0.0d0
		phi = 0.0d0
		total = 0.0d0
	endif
	
	pi = acos(-1.0d0)
	rad = 0.0d0
	
	lmax_comp = min(lmax, size(cilm(1,1,:))-1)
	
	intervalrad = interval*pi/180.0d0
	
	cilm2(1:2, 1:lmax_comp+1, 1:lmax_comp+1) = cilm(1:2, 1:lmax_comp+1, 1:lmax_comp+1)
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Convert coefficients to radial gravity if necessary 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (gravpot == "G" .or. gravpot == "g") then
		do l=0, lmax_comp, 1
			cilm2(1:2,l+1, 1:l+1) = -dble(l+1.0d0) * cilm2(1:2, l+1, 1:l+1)
		enddo
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Precomputing sines and cosines leads to an increase in speed by a factor of 
	!	almost 4 with no optimization, and by a factor of about 15 with normal optimizations.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do k=1, nlong
	
		longitude = longmin*pi/180.0d0 + dble(k-1)*intervalrad
		do m=0, lmax
			cosm(k, m+1) = cos(m*longitude)
			sinm(k, m+1) = sin(m*longitude)
		enddo
		
	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate default case of ONLY the potential or radial gravity
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if ( .not. present(total) ) then
	
		do j=1, nlat
	
			lat = latmax - dble(j-1)*interval
			x = sin(lat*pi/180.0d0)
		
			if (lat == 90.0d0 .or. lat == -90.0d0) then
		
				call PlBar(pl, lmax_comp, x)
	
				r_ex = a * (1.0d0 - f)	! Reference ellipsoid to calculate field
		
				coef = gm/r_ex
				if (gravpot == "g" .or. gravpot == "G") coef = coef / r_ex 
	
				rad(j,1) = cilm2(1,1,1) * pl(1) * coef	! degree 0
				do l = 1, lmax_comp, 1
					l1 = l+1
					coef = coef * (r0/r_ex)
					rad(j,1) = rad(j,1) + cilm2(1,l1,1) * pl(l1) * coef
				enddo
			
				rad(j, 2:nlong) = rad(j,1)
		
			else
		
				r_ex = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
					(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
				r_ex = sqrt(r_ex)
		
				call PlmBar(pl, lmax_comp, x, csphase = 1)
			
				coef = gm / r_ex
				if (gravpot == "g" .or. gravpot == "G") coef = coef / r_ex
		
				do l = 0, lmax_comp, 1
					if (l /= 0) coef = coef * (r0/r_ex) 
					l1 = l+1
					index = (l+1)*l/2 + 1
			
					rad(j,1:nlong) = rad(j,1:nlong) + cilm2(1,l1,1) * pl(index) * coef
			
					do k = 1, nlong
						
						do m = 1, l, 1
							m1 = m+1
							index = (l+1)*l/2 + m + 1
							rad(j,k) = rad(j,k) + ( cilm2(1,l1,m1)*cosm(k,m1) + &
								cilm2(2,l1,m1)*sinm(k,m1) ) * pl(index) * coef
						enddo
					enddo
				enddo
			
			endif
		
		enddo
	
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate all three components of the gravity field
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (present(total) .and. (gravpot == "G" .or. gravpot == "g") ) then
	
		do j=1, nlat
	
			lat = latmax - dble(j-1)*interval
			x = sin(lat*pi/180.0d0)
		
			if (lat == 90.0d0 .or. lat == -90.0d0) then
		
				call PlBar(pl, lmax_comp, x)
	
				r_ex = a * (1.0d0 - f)	! Reference ellipsoid to calculate field
		
				coef = gm/r_ex**2 
	
				rad(j,1) = cilm2(1,1,1) * pl(1) * coef	! degree 0
				do l = 1, lmax_comp, 1
					l1 = l+1
					coef = coef * (r0/r_ex)
					rad(j,1) = rad(j,1) + cilm2(1,l1,1) * pl(l1) * coef
				enddo
			
				rad(j, 2:nlong) = rad(j,1)
				theta(j, 1:nlong) = 0.0d0
				phi(j, 1:nlong) = 0.0d0
		
			else
		
				r_ex = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
					(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
				r_ex = sqrt(r_ex)
		
				call PlmBar_d1(pl, dpl, lmax_comp, x, csphase = 1)
				dpl = - dpl * cos(lat*pi/180.0d0)	! convert derivative with respect to z to
									! a derivative with respect to theta.
			
				coef = gm / r_ex**2
		
				do l = 0, lmax_comp, 1
					if (l /= 0) coef = coef * (r0/r_ex) 
					l1 = l+1
					index = (l+1)*l/2 + 1
			
					rad(j,1:nlong) = rad(j,1:nlong) + cilm2(1,l1,1) * pl(index) * coef
					theta(j,1:nlong) = theta(j,1:nlong) + cilm(1,l1,1) * dpl(index) * coef
			
					do k = 1, nlong
						
						do m = 1, l, 1
							m1 = m+1
							index = (l+1)*l/2 + m + 1
							rad(j,k) = rad(j,k) + ( cilm2(1,l1,m1)*cosm(k,m1) + &
								cilm2(2,l1,m1)*sinm(k,m1) ) * pl(index) * coef
							theta(j,k) = theta(j,k) + ( cilm(1,l1,m1)*cosm(k,m1) + &
								cilm(2,l1,m1)*sinm(k,m1) ) * dpl(index) * coef
							phi(j,k) = phi(j,k) + ( -cilm(1,l1,m1)*sinm(k,m1) + &
								cilm(2,l1,m1)*cosm(k,m1) ) * pl(index) * coef * dble(m) &
								/ cos(lat*pi/180.0d0)
						enddo
					enddo
				enddo
			
			endif
		
		enddo
			
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Add rotational effects
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (present(omega)) then
		
		do j=1, nlat
			lat = latmax - dble(j-1)*interval
			x = sin(lat*pi/180.0d0)
			r_ex = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
				(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
			r_ex = sqrt(r_ex)
			
			if (gravpot == "G" .or. gravpot == "g") then
				rad(j,1:nlong) = rad(j,1:nlong)  &
					+ r_ex * ( sin( (90.0d0-lat)*pi/180.0d0 ) * omega )**2 
				
				if (present(theta) ) then
					theta(j,1:nlong) = theta(j,1:nlong) &
						+ sin( (90.0d0-lat)*pi/180.0d0 ) * &
						cos( (90.0d0-lat)*pi/180.0d0 ) * r_ex * omega**2  
				endif
				
			else
				rad(j,1:nlong) = rad(j,1:nlong)  &
					+ 0.50d0 * ( r_ex * sin((90.0d0-lat)*pi/180.0d0) * omega )**2
			endif
		enddo
		
	endif
	
	if (present(total) .and. (gravpot == "G" .or. gravpot == "g")) then
	
		total(1:nlat, 1:nlong) = sqrt( rad(1:nlat,1:nlong)**2 + phi(1:nlat,1:nlong)**2 + theta(1:nlat,1:nlong)**2)
		
	endif
	
	! remove normal gravity from total gravitational acceleration
	if (present(normal_gravity)) then
		b = a * ( 1.0d0 - f)
		do j=1, nlat
			lat = latmax - dble(j-1)*interval	! geocentric latitude in degrees!
			total(j,1:nlong) = total(j,1:nlong) - NormalGravity(lat, GM, omega, a, b)
		enddo
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Deallocate memory
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (present(total) .and. (gravpot == "G" .or. gravpot == "g")) then
		 call PlmBar_d1(pl, dpl, -1, x, csphase = 1)
		 deallocate(dpl)
		 deallocate(pl)
	else
		call PlmBar(pl, -1, x, csphase = 1)
		deallocate(pl)
	endif
	
	deallocate(cosm)
	deallocate(sinm)
	deallocate(cilm2)
	
end subroutine MakeGravGrid2D


real*8 function NormalGravity(geocentric_lat, GM, omega, a, b)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Compute the total predicted gravity normal to the elipsoid at a given
!	GEOCENTRIC latitude (input in DEGREES), using Somigliana's formula.
!	This is taken from Physical Geodesy (Hofmann-Wellenhof and Moritz, sec. ed.),
!	sections 2.7 and 2.8.
!
!	Written by Mark Wieczorek (August 2009).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) ::	geocentric_lat, gm, omega, a, b
	real*8 ::		geodetic_lat, pi, ga, gb, m, ep, bigE, q0, q0p

	if (a<b) then
		print*, "Warning --- NormalGravity"
		print*, "The semimajor axis A should be greater than the semiminor axis B."
	endif
	
	pi = acos(-1.0d0)
	
	m = omega**2 * a**2 * b / GM
	
	bigE = sqrt(a**2 - b**2) ! linear eccentricity
	
	ep = bigE / b ! second eccentricity
	
	q0 = 0.50d0 * ( (1.0d0 + 3.0d0 * (b/bigE)**2) * atan(bigE/b) - 3.00 * b/bigE )
	
	q0p = 3.0d0 * (1.0d0 + (b/bigE)**2) * (1.0d0 - b/bigE * atan(bigE/b) ) - 1.0d0
	
	ga = GM / (a*b) * (1.0d0 - m - m * ep * q0p / 6.0d0 / q0) ! gravity on equator
		
	gb = GM / a**2 * (1.0d0 + m * ep * q0p / 3.0d0 / q0) ! gravity at poles
	
	geodetic_lat = atan((a/b)**2 * tan(geocentric_lat*pi/180.0d0))
	
	NormalGravity = a * ga * cos(geodetic_lat)**2 + b * gb * sin(geodetic_lat)**2
			
	NormalGravity = NormalGravity / sqrt(a**2 * cos(geodetic_lat)**2 + &
					b**2 * sin(geodetic_lat)**2 ) 
	
end function NormalGravity

