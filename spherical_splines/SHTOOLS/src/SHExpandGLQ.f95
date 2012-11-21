subroutine SHExpandGLQ(cilm, lmax, gridglq, w, plx, zero, norm, csphase, lmax_calc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
!	This program will expand a grid of data (regulary spaced in longitude,
!	and irregurarly in latitude according to the Guass-Legendre quadrature points) 
!	into spherical harmonics by utilizing an FFT of each latitudinal band, and a 
!	Guass-Legendre Quadrature in latitute. Note that the array PLX is
!	optional, and should not be precomputed when memory is an issue (i.e., LMAX>360).
!	It is implicitly assumed that the gridded data are bandlimited to degree LMAX.
!	If the optional parameter LMAX_CALC is specified, the spherical harmonic coefficients
!	will be calculated only up to and including this degree.
!
!	If PLX is not present, the Legendre functions are computed on the fly
!	using the scaling methodolgy presented in Holmes and Featherston (2002). 
!	When NORM is 1,2 or 4, these are accurate to degree 2800. When NORM is 3, the
!	routine is only stable to about degree 15.
!
!	Calling Parameters:
!		IN
!			gridglq		Gridded data used in the expansion generated
!					from a call to MakeGridGLQ, with dimensions
!					(LMAX+1, 2*LMAX+1).
!			lmax		Spherical harmonic bandwidth of the input grid. 
!					If LMAX_CALC is not specified, this corresponds
!					to the maximum spherical harmonic degree of the expansion.
!			w		Gauss-Legendre points used in the integrations
!					(determined from a call to PreCompute).
!		OUT
!			cilm 		Spherical harmonic coefficients of expansion with 
!					dimensions (2, LMAX+1, LMAX+1), or if LMAX_CAL is
!					present (2, LMAX_CALC+1, LMAX_CALC+1).
!		OPTIONAL (IN)
!			plx		Input array of Associated Legendre Polnomials computed
!					at the Gauss-Legendre points (determined from a call to
!					PreCompute). If this is not included, then the optional
!					array ZERO MUST be inlcuded.
!			zero		Array of dimension LMAX+1 that contains the latitudinal
!					gridpoints used in the Gauss-Legendre quadrature integration
!					scheme, calculated from a call to PreCompute. This is only 
!					needed if PLX is not given.
!			norm		Normalization to be used when calculating the Legendre functions
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!			csphase		1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			lmax_calc	The maximum spherical harmonic degree calculated in the 
!					spherical harmonic expansion.
!
!	Dependencies:	FFTW3, NGLQSH, CSPHASE_DEFAULT
!	
!	Written by Mark Wieczorek (October 2003)
!
!	September 3, 2005. 	Modifed so that the array plx is now optional.
!	September 26, 2005. 	Added optional argument NORM.
!	October 16, 2006. 	The Legendre functions are now computed within this program
!				during the summations over l and m. This leads to an increase in speed
!				of about a factor of 2.
!	July 23, 2007.		Added optional argument LMAX_CALC.
!	August 18, 2009.	Fixed bug that could have given incorrect results when LMAX_CAlC = 0
!	November 21, 2011	Fixed problem where saved variables used in Plm recursion were not recalculated
!				if NORM changed from one call to the next (with the same value of N).
!
!	Copyright (c) 2005-2011, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use FFTW3
	use SHTOOLS, only: NGLQSH, CSPHASE_DEFAULT
	
	implicit none
	real*8, intent(in) ::	w(:), gridglq(:,:)
	real*8, intent(in), optional ::	plx(:,:), zero(:)
	real*8, intent(out) ::	cilm(:,:,:)
	integer, intent(in) ::	lmax
	integer, intent(in), optional :: norm, csphase, lmax_calc
	integer :: 		nlong, nlat, i, l, m, k, kstart, l1, m1, phase, i_s, astat(4), &
				lnorm, lmax_comp
	real*8 :: 		pi, gridl(2*lmax+1), prod, scalef, rescalem, u, p, pmm, &
				pm1, pm2, z, fcoef1(2, lmax+1), fcoef2(2, lmax+1)
	complex*16 ::		cc(lmax+1)
	integer*8 :: 		plan
	real*8, save, allocatable ::	f1(:), f2(:), sqr(:), symsign(:)
	integer, save ::	lmax_old=0, norm_old = 0

	if (present(lmax_calc)) then
		lmax_comp = min(lmax, lmax_calc)
	else
		lmax_comp = lmax
	endif
	
	if (present(lmax_calc)) then
		if (lmax_calc > lmax) then
			print*, "Error --- SHExpandGLQ"
			print*, "LMAX_CALC must be less than or equal to the effective bandwidth of GRID."
			print*, "Effective bandwidth of GRID = N/2 - 1 = ", lmax
			print*, "LMAX_CALC = ", lmax_calc
			stop
		endif
	endif
	
	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax_comp+1 .or. size(cilm(1,1,:)) < lmax_comp+1) then
		print*, "Error --- SHExpandGLQ"
		print*, "CILM must be dimensioned as (2, LMAX_COMP+1, LMAX_COMP+1) where"
		print*, "LMAX_COMP = MIN(LMAX+1, LMAX_CALC+1)"
		print*, "LMAX = ", lmax
		if (present(lmax_calc)) print*, "LMAX_CALC = ", lmax_calc
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	elseif (size(gridglq(1,:)) < 2*lmax+1 .or. size(gridglq(:,1)) < lmax+1 ) then
		print*, "Error --- SHExpandGLQ"
		print*, "GRIDGLQ must be dimensioned as (LMAX+1, 2*LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(gridglq(:,1)), size(gridglq(1,:))
		stop
	elseif (size(w) < lmax+1) then
		print*, "Error --- SHExpandGLQ"
		print*, "W must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(w)
		stop
	endif

	
	if (present(plx)) then
		if (size(plx(:,1)) < lmax+1 .or. size(plx(1,:)) < (lmax+1)*(lmax+2)/2) then
			print*, "Error --- SHExpandGLQ"
			print*, "PLX must be dimensioned as (LMAX+1, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
			print*, "Input array is dimensioned as ", size(plx(:,1)), size(plx(1,:))
			stop
		endif
	elseif (present(zero)) then
		if (size(zero) < lmax+1) then
			print*, "Error --- SHExpandGLQ"
			print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
			print*, "Input array is dimensioned ", size(zero)
			stop
		endif
	else
		print*, "Error --- SHExpandGLQ"
		print*, "Either PLX or ZERO must be specified."
		stop
	endif
	
	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error - SHExpandGLQ"
			print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
			print*, "Input value is ", norm
			stop
		endif
		
		lnorm = norm
		
	else
		lnorm = 1
	endif
	
	if (present(csphase)) then
     		if (csphase /= -1 .and. csphase /= 1) then
     			print*, "Error --- SHExpandGLQ"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif
     	
 	
	nlong = 2*(lmax+1) -1	! This is the number of points (and period) used
				! in the FFT in order that m=0 to lmax.
	
	nlat = NGLQSH(lmax)	! nlat = lmax+1
	
	pi = acos(-1.0d0)
	
	cilm = 0.0d0
	
	scalef = 1.0d-280
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate recursion constants used in computing the Legendre polynomials
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if ( (lmax_comp /= lmax_old .or. lnorm /= norm_old) .and. .not. present(plx) ) then

		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		if (allocated(symsign)) deallocate(symsign)
		
		allocate(sqr(2*lmax_comp+1), stat=astat(1))
		allocate(f1((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(2))
		allocate(f2((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(3))
		allocate(symsign((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(4))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. astat(4) /= 0) then
			print*, "Error --- SHExpandGLQ"
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
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Create generic plan for gridl
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call dfftw_plan_dft_r2c_1d(plan, nlong, gridl(1:nlong), cc, FFTW_MEASURE)
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine Cilms, one l at a time I by integrating over all 
	!	latitudes using Gauss-Legendre Quadrature. When PLX is not 
	!	present, the Legendre functions are computed on the fly 
	!	during the summations over l and m. These are scaled using
	!	the methodology of Holmesand Featherstone (2002), with the
	!	exception of the m=0 terms that do not need to be scaled 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	if (present(plx)) then
	
		do i=1, nlat
		
			gridl(1:nlong) = gridglq(i,1:nlong)
			call dfftw_execute(plan)				! take fourier transform
			fcoef1(1,1:lmax+1) = 2.0d0 * w(i) * pi * dble(cc(1:lmax+1)) / dble(nlong)
			fcoef1(2,1:lmax+1) = -2.0d0 * w(i) * pi * dimag(cc(1:lmax+1)) / dble(nlong)
			
			k = 0
			do l = 0, lmax_comp, 1
				l1 = l + 1
				do m = 0, l, 1
					m1 = m+1
					k = k + 1
					cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + plx(i,k) * fcoef1(1:2,m1)
				enddo
			enddo
		enddo
		
	else
				
		do i=1, (nlat+1)/2
		
			if (i == (nlat+1)/2 .and. mod(nlat,2) /= 0) then		! This latitude is the equator; z=0, u=1
			
				gridl(1:nlong) = gridglq(i,1:nlong)
				call dfftw_execute(plan)				! take fourier transform
				fcoef1(1,1:lmax+1) = 2.0d0 * w(i) * pi * dble(cc(1:lmax+1)) / dble(nlong)
				fcoef1(2,1:lmax+1) = -2.0d0 * w(i) * pi * dimag(cc(1:lmax+1)) / dble(nlong)
	
				u = 1.0d0
				
				select case(lnorm)
					case(1,2,3);	pm2 = 1.0d0
					case(4);	pm2 = 1.0d0 / sqrt(4.0d0*pi)
				end select
				
				cilm(1,1,1) = cilm(1,1,1) + pm2 * fcoef1(1,1)
				
				if (lmax_comp == 0) cycle
				
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
				
					m1 = m + 1
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
        			
        			cilm(1:2,lmax_comp+1,lmax_comp+1) = cilm(1:2,lmax_comp+1,lmax_comp+1) &
        					+ pmm * fcoef1(1:2,lmax_comp+1) 

				
			else			
			
				z = zero(i)
				u = sqrt( (1.0d0-z) * (1.0d0+z) )
				
				gridl(1:nlong) = gridglq(i,1:nlong)
				call dfftw_execute(plan)			
				fcoef1(1,1:lmax+1) = 2.0d0 * w(i) * pi * dble(cc(1:lmax+1)) / dble(nlong)
				fcoef1(2,1:lmax+1) = -2.0d0 * w(i) * pi * dimag(cc(1:lmax+1)) / dble(nlong)
				
				i_s = nlat + 1 - i	! point symmetric about the equator
				
				gridl(1:nlong) = gridglq(i_s,1:nlong)
				call dfftw_execute(plan)			
				fcoef2(1,1:lmax+1) = 2.0d0 * w(i_s) * pi * dble(cc(1:lmax+1)) / dble(nlong)
				fcoef2(2,1:lmax+1) = -2.0d0 * w(i_s) * pi * dimag(cc(1:lmax+1)) / dble(nlong)

				
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
			
			endif
		enddo
	
	endif
	
	call dfftw_destroy_plan(plan) 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Divide by integral of Ylm*Ylm 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	select case(lnorm)
	
		case(1)
		
			do l=0, lmax_comp, 1
				cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) / (4.0d0*pi)
			enddo	
		
		case(2)
	
			do l=0, lmax_comp, 1
				cilm(1:2,l+1, 1:l+1) = dble(2*l+1) * cilm(1:2,l+1, 1:l+1) / (4.0d0*pi)
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
	
end subroutine SHExpandGLQ

