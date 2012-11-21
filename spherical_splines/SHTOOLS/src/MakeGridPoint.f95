real*8 function MakeGridPoint(cilm, lmax, lat, longitude, norm, csphase)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will determine the value at a given latitude and 
!	longitude corresponding to the given set of spherical harmonics.
!	Latitude and Longitude are assumed to be in DEGREES!
!
!	Calling Parameters:
!		IN
!			cilm:	Spherical harmonic coefficients, with dimensions
!				(2, lmax+1, lmax+1).
!			lmax:	Maximum degree used in the expansion.
!			lat:	latitude (degrees).
!			long:	longitude(degrees).
!		OPTIONAL (IN)
!			norm	Spherical harmonic normalization:
!					(1) "geodesy" (default)
!					(2) Schmidt
!					(3) unnormalized
!					(4) orthonormalized
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!
!
!	Dependencies:		PlmBar, PlBar, PlmSchmidt, PlmON, CSPHASE_DEFAULT
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, then this file will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!
!	Written by Mark Wieczorek (June 2004)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PlmBar, PLegendreA, PlmSchmidt, PlmON, CSPHASE_DEFAULT

	implicit none
	
	real*8, intent(in):: 	cilm(:,:,:), lat, longitude
	integer, intent(in) ::	lmax
	integer, intent(in), optional ::	norm, csphase
	real*8 ::		pi, x, expand, lon, mlong
	integer ::		index, l, m, l1, m1, lmax_comp, phase, astat
	real*8, allocatable :: 	 pl(:)
	
	
	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- MakeGridPoint"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
		
	if (present(norm)) then
		if (norm >4 .or. norm < 1) then
			print*, "Error - MakeGridPoint"
			print*, "Parameter NORM must be 1, 2, 3, or 4"
			stop
		endif
	endif
	
    	if (present(csphase)) then
     		if (csphase == -1) then
     			phase = -1.0d0
     		elseif (csphase == 1) then
     			phase = 1.0d0
     		else
     			print*, "Error --- MakeGridPoint"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		endif
     	else
     		phase = dble(CSPHASE_DEFAULT)
     	endif
     	
     	allocate(pl((lmax+1)*(lmax+2)/2), stat = astat)
     	if (astat /= 0) then
		print*, "Error --- MakeGridPoint"
		print*, "Problem allocating array PL", astat
		stop
	endif
 

	pi = acos(-1.0d0)
	x = sin(lat*pi/180.0d0)
	lon = longitude * pi/180.0d0
	
	lmax_comp = min(lmax, size(cilm(1,1,:))-1)
	
	if (present(norm)) then
		if (norm == 1) call PlmBar(pl, lmax_comp, x, csphase = phase)
		if (norm == 2) call PlmSchmidt(pl, lmax_comp, x, csphase = phase)
		if (norm == 3) call PLegendreA(pl, lmax_comp, x, csphase = phase)
		if (norm == 4) call PlmON(pl, lmax_comp, x, csphase = phase)
	else
		call PlmBar(pl, lmax_comp, x, csphase = phase)
	endif
	
	expand = 0.0d0

	do l = lmax_comp, 0, -1
	
		l1 = l+1
		index = (l+1)*l/2 + 1
		expand = expand + cilm(1,l1,1) * pl(index)

		do m = 1, l, 1
			m1 = m+1
			mlong = m*lon
			index = (l+1)*l/2 + m + 1
			expand = expand + ( cilm(1,l1,m1) * cos(mlong) + &
				cilm(2,l1,m1) * sin(mlong) ) * pl(index)
		enddo
	enddo
	
	MakeGridPoint = expand
	
	! deallocate memory
	if (present(norm)) then
		if (norm == 1) call PlmBar(pl, -1, x, csphase = phase)
		if (norm == 2) call PlmSchmidt(pl, -1, x, csphase = phase)
		if (norm == 4) call PlmON(pl, -1, x, csphase = phase)
	else
		call PlmBar(pl, -1, x, csphase = phase)
	endif
	deallocate(pl)
	
end function MakeGridPoint

