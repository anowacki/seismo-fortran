subroutine SHExpandDH(grid, n, cilm, lmax, norm, sampling, csphase, lmax_calc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will expand a grid containing N samples in both longitude 
!	and latitude (or N x 2N, see below) into spherical harmonics. This routine 
!	makes use of the sampling theorem presented in Driscoll and Healy (1994) 
!	and employs FFTs when calculating the sin and cos terms. The number of 
!	samples, N, must be EVEN for this routine to work, and the spherical 
!	harmonic expansion is exact if the function is bandlimited to degree N/2-1.
!	Legendre functions are computed on the fly using the scaling methodolgy 
!	presented in Holmes and Featherston (2002). When NORM is 1,2 or 4, these are 
!	accurate to about degree 2800. When NORM is 3, the routine is only stable 
!	to about degree 15. If the optional parameter LMAX_CALC is specified, the 
!	spherical harmonic coefficients will only be calculated up to this degree. 
!	
!	If SAMPLING is 1 (default), the input grid contains N samples in latitude from 
!	90 to -90+interval, and N samples in longitude from 0 to 360-2*interval, where 
!	interval is the latitudinal sampling interval 180/N. Note that the datum at 
!	90 degees North latitude is ultimately downweighted to zero, so this point 
!	does not contribute to the spherical harmonic coefficients. If SAMPLING is 2,
!	the input grid must contain N samples in latitude and 2N samples in longitude.
!	In this case, the sampling intervals in latitude and longitude are 180/N and 
!	360/N respectively. When performing the FFTs in longitude, the frequencies greater
!	than N/2-1 are simply discarded to prevent aliasing.
!
!	Calling Parameters:
!		IN
!			grid		Equally sampled grid in latitude and longitude of dimension
!					(1:N, 1:N) or and equally spaced grid of dimension (1:N,2N).
!			N		Number of samples in latitude and longitude (for SAMPLING=1),
!					or the number of samples in latitude (for SAMPLING=2).
!		OUT
!			cilm		Array of spherical harmonic coefficients with dimension
!					(2, LMAX+1, LMAX+1), or, if LMAX_CALC is present 
!					(2, LMAX_CALC+1, LMAX_CALC+1).
!			lmax		Spherical harmonic bandwidth of the grid. This corresponds
!					to the maximum spherical harmonic degree of the expansion
!					if the optional parameter LMAX_CALC is not specified.
!		OPTIONAL (IN)
!			norm		Normalization to be used when calculating Legendre functions
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!			sampling	(1) Grid is N latitudes by N longitudes (default).
!					(2) Grid is N by 2N. The higher frequencies resulting
!					from this oversampling are discarded, and hence not
!					aliased into lower frequencies.
!			csphase		1: Do not include the Condon-Shortley phase factor of (-1)^m.
!					-1: Apply the Condon-Shortley phase factor of (-1)^m.
!			lmax_calc	The maximum spherical harmonic degree calculated in the 
!					spherical harmonic expansion.
!
!	Notes:
!		1.	This routine does not use the fast legendre transforms that
!		   	are presented in Driscoll and Heally (1994).
!		2. 	Use of a N by 2N grid is implemented because many geographic grids are 
!			sampled this way. When taking the Fourier transforms in longitude, all of the 
!			higher frequencies are ultimately discarded. If, instead, every other column of the 
!			grid were discarded to form a NxN grid, higher frequencies could be aliased 
!			into lower frequencies.
!
!	Dependencies:		DHaj, FFTW3, CSPHASE_DEFAULT
!	
!
!	Written by Mark Wieczorek (May 2004).
!
!	September 3, 2005. 	The final triple sum was reordered, and the array plm, which
!				was orginally (n, ((n/2)*(n/2+1))/2), was changed to 
!				((n/2)*(n/2+1))/2 in order to save memory.
!	September 26, 2005. 	Added optional argument NORM.
!	May 17, 2006. 		Added option for using Nx2N grids.
!	May 30, 2006. 		Modified routine to take into account the symmetry of the Legendre 
!				functions about the equator, i.e., (-1)^(l+m).
!	October 17, 2006. 	The Legendre functions are now computed within this program
!				during the summations over l and m. This leads to an increase in 
!				speed of about a factor of 2.
!	July 23, 2007. 		Added optional parameter LMAX_CALC.
!	August 18, 2009.	Fixed memory issue that could have given rise to errors when
!				called with LMAX_CALC = 0.
!	November 21, 2011	Fixed problem where saved variables used in Plm recursion were not recalculated
!				if NORM changed from one call to the next (with the same value of N).
!
!	Copyright (c) 2006-2011, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use FFTW3
	use SHTOOLS, only: DHaj, CSPHASE_DEFAULT
		
	implicit none
	
	real*8, intent(in) ::	grid(:,:)
	real*8, intent(out) ::	cilm(:,:,:)
	integer, intent(in) ::	n
	integer, intent(out) ::	lmax
	integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
	complex*16 ::		cc(n+1)
	integer ::		l, m, k, kstart, i, l1, m1, i_eq, i_s, phase, lnorm, astat(5), &
				lmax_comp, nlong
	integer*8 ::		plan
	real*8 ::		pi, gridl(2*n), aj(n), fcoef1(2, n/2+1), fcoef2(2, n/2+1), &
				theta, prod, scalef, rescalem, u, p, pmm, pm1, pm2, z
	real*8, save, allocatable ::	f1(:), f2(:), sqr(:), symsign(:)
	integer, save ::	lmax_old=0, norm_old = 0
	
	lmax = n/2 - 1	
	
	if (present(lmax_calc)) then
		lmax_comp = min(lmax, lmax_calc)
	else
		lmax_comp = lmax
	endif
	
	if (present(lmax_calc)) then
		if (lmax_calc > lmax) then
			print*, "Error --- SHExpandDH"
			print*, "LMAX_CALC must be less than or equal to the effective bandwidth of GRID."
			print*, "Effective bandwidth of GRID = N/2 - 1 = ", lmax
			print*, "LMAX_CALC = ", lmax_calc
			stop
		endif
	endif
			
	if (present(sampling)) then
		if (sampling /= 1 .and. sampling /=2) then
			print*, "Error --- SHExpandDH"
			print*, "Optional parameter SAMPLING must be 1 (N by N) or 2 (N by 2N)."
			print*, "Input value is ", sampling
			stop
		endif
	endif
	
	if (mod(n,2) /=0) then
		print*, "Error --- SHExpandDH"
		print*, "The number of samples in latitude and longitude, n, must be even."
		print*, "Input value is ", n
		stop
	elseif (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) <  lmax_comp+1 .or. &
			size(cilm(1,1,:)) < lmax_comp+1) then
		print*, "Error --- SHExpandDH"
		print*, "CILM must be dimensioned as (2, LMAX_COMP+1, LMAX_COMP+1) where"
		print*, "LMAX_COMP = MIN(N/2, LMAX_CALC+1)"
		print*, "N = ", n
		if (present(lmax_calc)) print*, "LMAX_CALC = ", lmax_calc
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if (present(sampling)) then
		if (sampling == 1) then
			if (size(grid(:,1)) < n .or. size(grid(1,:)) < n) then
				print*, "Error --- SHExpandDH"
				print*, "GRIDDH must be dimensioned as (N, N) where N is ", n
				print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
				stop
			endif
		elseif (sampling ==2) then
			if (size(grid(:,1)) < n .or. size(grid(1,:)) < 2*n) then
				print*, "Error --- SHExpandDH"
				print*, "GRIDDH must be dimensioned as (N, 2*N) where N is ", n
				print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
				stop
			endif
		endif
	else
		
		if (size(grid(:,1)) < n .or. size(grid(1,:)) < n) then
			print*, "Error --- SHExpandDH"
			print*, "GRIDDH must be dimensioned as (N, N) where N is ", n
			print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
			stop
		endif

	endif
			
		
	if (present(csphase)) then
     		if (csphase /= -1 .and. csphase /= 1) then
     			print*, "SHExpandDH --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif
     	
     	
    	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error --- SHExpandDH"
			print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
			print*, "Input value is ", norm
			stop
		endif
		
		lnorm = norm
		
	else
		lnorm = 1
	endif
	

	pi = acos(-1.0d0)
	
	cilm = 0.0d0
	
	scalef = 1.0d-280
	
	call DHaj(n, aj)
	aj(1:n) = aj(1:n)*sqrt(4.0d0*pi) 	! Driscoll and Heally use unity normalized spherical harmonics
	
	if (present(sampling)) then
		if (sampling == 1) then
			nlong = n
		else
			nlong = 2*n
		endif
	else
		nlong = n
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate recursion constants used in computing the Legendre polynomials
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (lmax_comp /= lmax_old .or. lnorm /= norm_old) then

		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		if (allocated(symsign)) deallocate(symsign)
		
		allocate(sqr(2*lmax_comp+1), stat=astat(1))
		allocate(f1((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(2))
		allocate(f2((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(3))
		allocate(symsign((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(4))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. astat(4) /= 0) then
			print*, "Error --- SHExpandDH"
			print*, "Problem allocating arrays SQR, F1, F2, or SYMSIGN", astat(1), astat(2), astat(3), astat(4)
			stop
		endif
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Calculate signs used for symmetry of Legendre functions about equator
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		k = 0
		do l = 0, lmax_comp, 1
			do m = 0, l, 1
				k = k + 1
				if (mod(l-m,2) == 0) then
					symsign(k) = 1.0d0
				else
					symsign(k) = -1.0d0
				endif
			enddo
		enddo
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax_comp+1
			sqr(l) = sqrt(dble(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		P(l,m) = x*f1(l,m)*P(l-1,m) - P(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l as a different 
		!	recursion is used. Furthermore, for m=l-1, Plmbar(l-2,m) is assumed to be zero.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		k = 1
		
		select case(lnorm)
		
			case(1,4)
	
				if (lmax_comp /= 0) then
					k = k + 1
					f1(k) = sqr(3)
					f2(k) = 0.0d0
					k = k + 1
				endif
				
				do l=2, lmax_comp, 1
					k = k + 1
					f1(k) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
					f2(k) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)
					do m=1, l-2, 1
						k = k+1
						f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                				f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                  				 	/ sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
					enddo
					k = k+1
					f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                			f2(k) = 0.0d0
					k = k + 1
				enddo
			
			case(2)
			
				if (lmax_comp /= 0) then
					k = k + 1
					f1(k) = 1.0d0
					f2(k) = 0.0d0
					k = k + 1
				endif
				
				do l=2, lmax_comp, 1
					k = k + 1
					f1(k) = dble(2*l-1) /dble(l)
					f2(k) = dble(l-1) /dble(l)
					do m=1, l-2, 1
						k = k+1
						f1(k) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  				f2(k) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
					enddo
					k = k+1
					f1(k) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  			f2(k) = 0.0d0
					k = k + 1
				enddo
			
			case(3)
		
				do l=1, lmax_comp, 1
					k = k + 1
					f1(k) = dble(2*l-1) /dble(l)
					f2(k) = dble(l-1) /dble(l)
					do m=1, l-1, 1
						k = k+1
						f1(k) = dble(2*l-1) / dble(l-m)
                  				f2(k) = dble(l+m-1) / dble(l-m)
					enddo
					k = k + 1
				enddo

		end select
	
		lmax_old = lmax_comp
		norm_old = lnorm
	
	endif	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Create generic plan for gridl
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call dfftw_plan_dft_r2c_1d_(plan, nlong, gridl(1:nlong), cc, FFTW_MEASURE)

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		
	! 	Integrate over all latitudes. Take into account symmetry of the 
	!	Plms about the equator.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	i_eq = n/2 + 1	! Index correspondong to the equator
	
	do i=2, i_eq - 1, 1
	
		theta = pi * dble(i-1)/dble(n)
		z = cos(theta)
		u = sqrt( (1.0d0-z) * (1.0d0+z) )
		
		gridl(1:nlong) = grid(i,1:nlong)
		call dfftw_execute_(plan)	! take fourier transform
		fcoef1(1,1:n/2) = sqrt(2.0d0*pi) * aj(i) * dble(cc(1:n/2)) / dble(nlong)
		fcoef1(2,1:n/2) = -sqrt(2.0d0*pi) * aj(i) * dimag(cc(1:n/2)) / dble(nlong)
		
		i_s = 2*i_eq -i
		
		gridl(1:nlong) = grid(i_s,1:nlong)
		call dfftw_execute_(plan)	! take fourier transform
		fcoef2(1,1:n/2) = sqrt(2.0d0*pi) * aj(i_s) * dble(cc(1:n/2)) / dble(nlong)
		fcoef2(2,1:n/2) = -sqrt(2.0d0*pi) * aj(i_s) * dimag(cc(1:n/2)) / dble(nlong)
		
		select case(lnorm)
			case(1,2,3);	pm2 = 1.0d0
			case(4);	pm2 = 1.0d0 / sqrt(4.0d0*pi)
		end select

		cilm(1,1,1) = cilm(1,1,1) + pm2 * (fcoef1(1,1) + fcoef2(1,1) )
		! symsign(1) = 1
		
		if (lmax_comp == 0) cycle
				
		k = 2
		pm1 =  f1(k) * z * pm2
		cilm(1,2,1) = cilm(1,2,1) + pm1 * ( fcoef1(1,1) - fcoef2(1,1) )
		! symsign(2) = -1
				
		do l=2, lmax_comp, 1
			l1 = l+1
			k = k+l
			p = f1(k) * z * pm1 - f2(k) * pm2
			pm2 = pm1
			pm1 = p
			cilm(1,l1,1) = cilm(1,l1,1) + p * ( fcoef1(1,1) + fcoef2(1,1) * symsign(k) )
		enddo
				
		select case(lnorm)
			case(1,2);	pmm = sqr(2) * scalef
			case(3);	pmm = scalef
			case(4);	pmm = sqr(2) * scalef / sqrt(4.0d0*pi)
		end select
				
		rescalem = 1.0d0/scalef
		kstart = 1
	
		do m = 1, lmax_comp-1, 1
				
			m1 = m+1
			rescalem = rescalem * u
			kstart = kstart+m+1
					
			select case(lnorm)
				case(1,4)
					pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm
				case(2)
					pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm / sqr(2*m+1)
				case(3)
					pmm = phase * pmm * dble(2*m-1)
					pm2 = pmm
			end select
					
			fcoef1(1:2,m1) = fcoef1(1:2,m1) * rescalem
			fcoef2(1:2,m1) = fcoef2(1:2,m1) * rescalem
					
			cilm(1:2,m1,m1) = cilm(1:2,m1,m1) + pm2 * &
					( fcoef1(1:2,m1) + fcoef2(1:2,m1) )
			! symsign(kstart) = 1
					
			k = kstart+m+1
	   		pm1 = z * f1(k) * pm2
	   			
	   		cilm(1:2,m1+1,m1) = cilm(1:2,m1+1,m1) + pm1 * &
               				( fcoef1(1:2,m1) - fcoef2(1:2,m1) )
               		! symsign(k) = -1
					
			do l = m+2, lmax_comp, 1
				l1 = l+1
				k = k+l
                  		p = z * f1(k) * pm1-f2(k) * pm2
                  		pm2 = pm1
                  		pm1 = p
				cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + p * &
						( fcoef1(1:2,m1) + fcoef2(1:2,m1) * symsign(k) )
			enddo
               				
		enddo
				
		rescalem = rescalem * u
			
            	select case(lnorm)
            		case(1,4);	pmm = phase * pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            		case(2);	pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
            		case(3);	pmm = phase * pmm * dble(2*lmax_comp-1) * rescalem
        	end select
        			
        	cilm(1:2,lmax_comp+1,lmax_comp+1) = cilm(1:2,lmax_comp+1,lmax_comp+1) + pmm * &
        				( fcoef1(1:2,lmax_comp+1) + fcoef2(1:2,lmax_comp+1) )	
        				! symsign = 1
	enddo
	
	! Finally, do equator
	
	i = i_eq
	
	z = 0.0d0
	u = 1.0d0
	
	gridl(1:nlong) = grid(i,1:nlong)
	call dfftw_execute_(plan)	! take fourier transform
	fcoef1(1,1:n/2) = sqrt(2.0d0*pi) * aj(i) * dble(cc(1:n/2)) / dble(nlong)
	fcoef1(2,1:n/2) = -sqrt(2.0d0*pi) * aj(i) * dimag(cc(1:n/2)) / dble(nlong)


	select case(lnorm)
		case(1,2,3);	pm2 = 1.0d0
		case(4);	pm2 = 1.0d0 / sqrt(4.0d0*pi)
	end select

	cilm(1,1,1) = cilm(1,1,1) + pm2 * fcoef1(1,1)
	
	if (lmax_comp /= 0) then
				
		k = 2
				
		do l=2, lmax_comp, 2
			l1 = l+1
			k = k+l
			p = - f2(k) * pm2
			pm2 = p
			cilm(1,l1,1) = cilm(1,l1,1) + p * fcoef1(1,1)
			k = k + l + 1
		enddo
				
		select case(lnorm)
			case(1,2);	pmm = sqr(2) * scalef
			case(3);	pmm = scalef
			case(4);	pmm = sqr(2) * scalef / sqrt(4.0d0*pi)
		end select
				
		rescalem = 1.0d0/scalef
		kstart = 1
	
		do m = 1, lmax_comp-1, 1
				
			m1 = m+1
			kstart = kstart+m+1
					
			select case(lnorm)
				case(1,4)
					pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm
				case(2)
					pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm / sqr(2*m+1)
				case(3)
					pmm = phase * pmm * dble(2*m-1)
					pm2 = pmm
			end select

			fcoef1(1:2,m1) = fcoef1(1:2,m1) * rescalem
					
			cilm(1:2,m1,m1) = cilm(1:2,m1,m1) + pm2 * fcoef1(1:2,m1)
					
			k = kstart+m+1
	   						
			do l = m+2, lmax_comp, 2
				l1 = l+1
				k = k+l
                	  	p = - f2(k) * pm2
                  		pm2 = p
				cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + p * fcoef1(1:2,m1) 
				k = k + l + 1
			enddo

		enddo		
       
       		select case(lnorm)
            		case(1,4);	pmm = phase * pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            		case(2);	pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
            		case(3);	pmm = phase * pmm * dble(2*lmax_comp-1) * rescalem
       		end select
       			
       		cilm(1:2,lmax_comp+1,lmax_comp+1) = cilm(1:2,lmax_comp+1,lmax_comp+1) + pmm * fcoef1(1:2,lmax_comp+1) 
       		
       	endif

	call dfftw_destroy_plan_(plan) 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Divide by integral of Ylm*Ylm 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	select case(lnorm)
	
		case(1)
	
			do l=0, lmax_comp, 1
				cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) / (4.0d0*pi)
			enddo
		
		case(2)
	
			do l=0, lmax_comp, 1
				cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) * dble(2*l+1) / (4.0d0*pi)
			enddo
	
		case(3)
		 
			do l=0, lmax_comp, 1
				prod = 4.0d0*pi/dble(2*l+1)
				cilm(1,l+1,1) = cilm(1,l+1,1) / prod
				prod = prod / 2.0d0
				do m=1, l-1, 1
					prod = prod * dble(l+m) * dble(l-m+1)
					cilm(1:2,l+1,m+1) = cilm(1:2,l+1,m+1) / prod
				enddo
				!do m=l case
				if (l /= 0) cilm(1:2,l+1,l+1) = cilm(1:2,l+1, l+1)/(prod*dble(2*l))
			enddo
			
	end select
		
end subroutine SHExpandDH

