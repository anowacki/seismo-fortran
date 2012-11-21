!===============================================================================
!  spherical splines contains a few functions for fitting and evaluating 
!  spherical splines, as given by 
!     Wang & Dahlen (1995)  Spherical-spline parameterization of three-dimensional
!     Earth models.  GRL (22) 3099-3102.
!
!  This module uses LAPACK, and so must be compiled on a system with a version
!  available.  On Mac OS X, one can use the optimised 'vecLib' framework
!  by compiling with the options '-framework vecLib -llapack -lblas
!
!-------------------------------------------------------------------------------
!
!  Andy Nowacki, 
!  University of Bristol
!
!-------------------------------------------------------------------------------
!  HISTORY:
!     2011-01-30: Inception date.  First working version (of sorts).
!     2012-02-24: Added alternative spline basis function, the Abel-Poisson
!                 kernel (sph_slines_func_AP), and subroutine to calculate the
!                 most appropriate value for its parameter, h, for a desired
!                 spline width.  Change splines_set type to keep track of which
!                 kernel we're using.
!
!-------------------------------------------------------------------------------
!  TODO:
!     + Allow for strongly irregular knot spacing by calculating f using a
!       local value for dbar, the mean knot spacing (in that region).
!===============================================================================
module spherical_splines
!===============================================================================

   implicit none

!  ** size constants
      integer, parameter, private :: i4 = selected_int_kind(9)       ! long int
      integer, parameter, private :: r4 = selected_real_kind(6,37)   ! SP
      integer, parameter, private :: r8 = selected_real_kind(15,307) ! DP
      
!  ** precision selector
      integer, parameter, private :: rs = r8
      
!  ** maths constants and other useful things
      real(rs), parameter, private :: pi = 3.141592653589793238462643_rs ;
      real(rs), parameter, private :: pi2 = pi/2._rs
      real(rs), parameter, private :: ZERO = 0._rs, ONE = 1._rs, TWO = 2._rs, &
         HALF = 0.5_rs, FOURTH = 0.25_rs, THIRD = 1._rs/3._rs
      
!  ** logical units
      integer, parameter, private :: lu_stdin  = 5
      integer, parameter, private :: lu_stdout = 6
      integer, parameter, private :: lu_stderr = 0
      integer, parameter, private :: lu_file = 10  ! LU to be used for reading and writing
      
!  ** Names and number of different kernels available
      integer, parameter :: NUM_KERNELS = 2  ! Number of kernel types available
      integer, parameter :: KERNEL_NAME_LEN = 2, &
                            KERNEL_LONG_NAME_LEN = 30
      character(len=KERNEL_NAME_LEN),dimension(NUM_KERNELS) :: &
                                          KERNEL_NAME = (/ 'B ','AP' /)
      character(len=KERNEL_LONG_NAME_LEN),dimension(NUM_KERNELS) :: &
        KERNEL_LONG_NAME = (/ 'B-spline                      ', &
                              'Abel-Poisson                  ' /)

!  ** Hide the helper functions
      private :: sph_splines_delta
      private :: sph_splines_fatal
      private :: sph_splines_info
      
!  ** derived types used by the module, containing knot/data points, spline
!     coefficients and so on
      type sph_splines_set
         real(rs),dimension(:),allocatable :: phi,theta, & ! Lon & colat
                                              c, &  ! Value at control points
                                              a     ! Spline coefficient
         real(rs) :: dbar = -1._rs, & ! Mean knot spacing globally (for B-splines) &
                     h = -1._rs       !  Abel-Poisson coefficient (-1 if unset)
         integer :: n                 ! Number of control knots
         logical :: degrees = .true.  ! Knot units; defaults to degrees
         character(len=KERNEL_NAME_LEN) :: kernel ! Type of kernel being used
      end type

!-------------------------------------------------------------------------------
CONTAINS

!===============================================================================
subroutine sph_splines_fit(s,quiet)
!===============================================================================
!  Calculate the coefficients of the spherical spline representation of a surface
!  field, with knot points as given in the data.
!  
!  INPUT: sph_splines_set structure.  The following are unchanged
!         phi  \  Column vectors of implicit length which contain coordinates
!         theta } of the knot points and the surface values there.
!         c    /
!         degress Logical variable specifying whether input coordinates
!                 are in degrees or not.  (Default .true.)
!         
!  OUTPUT: sph_splines_set structure.  The following are changed
!         a      Column vector of coefficients which can be used to calculate
!                the fitted function.  (Re-)Allocated within the subroutine.
!         dbar   Mean distance between knots points globally.
!
!  Constructs the set of knot coefficients by inverting the linear equations
!     c(theta_i,phi_i) = f_j(Delta_ij) * a(i)
!  ->                a = F^-1 * c
!  For this, I use LAPACK

   implicit none
   
   type(sph_splines_set) :: s
   logical,optional,intent(in) :: quiet
   logical :: silent
   integer :: n  ! number of knots; length of phi,theta,r
   integer :: i,j,ier,nd,npoints
   integer,allocatable :: ipiv(:)  ! Pivot rows used by LAPACK
   real(rs),allocatable :: F(:,:)   ! Matrix of distances between point i and j
   real(rs) :: dmin

!  Get the level of verbosity.  Default to noisy
   silent = .false.
   if (present(quiet)) silent = quiet

!  Check for the right things being present
   if (.not.(allocated(s%phi) .and. allocated(s%theta) .and. allocated(s%c))) &
      call sph_splines_fatal('sph_splines_fit: location or values of knot points not allocated')
   
!  Get length of input arrays and make sure they're the same length
   n = size(s%phi)
   if (s%n /= n .or. size(s%theta) /= n .or. size(s%c) /= n) &
      call sph_splines_fatal('sph_splines_fit: structure arrays must all be same length.')
      
!  (Re-)Allocate s%a if already allocated to a length != n
   if (allocated(s%a)) then
      if (size(s%a) /= n) then
         deallocate(s%a);  allocate(s%a(n))
         if (.not.silent) call sph_splines_info( &
               'sph_splines_fit: allocating memory for spline coefficients.')
      endif
   else
      allocate(s%a(n))
      if (.not.silent) call sph_splines_info( &
            'sph_splines_fit: allocating memory for spline coefficients.')
   endif
   
!  Construct matrix of distances between knots, F_ij = Delta_ij
   if (.not.silent) call sph_splines_info( &
      'sph_splines_fit: Calculating inter-knot spacings.')
   allocate(F(n,n))
   do i=1,n
      do j=i,n
         if (i==j) then
            F(i,j) = 0._rs
         else
            F(i,j) = sph_splines_delta(s%phi(i),s%theta(i),s%phi(j),s%theta(j),degrees=s%degrees)
         endif
         F(j,i) = F(i,j)
      enddo
   enddo
   
!  Work out average distance between adjacent points globally, dbar
   if (.not.silent) call sph_splines_info( &
      'sph_splines_fit: Calculating mean knot spacing.')
   s%dbar = 0._rs
   nd = 0
   do i=1,n
!  Take the nearest 6 or so points and take the average distance of these
      dmin = 0._rs
      npoints = 0
      do while (npoints <= 6)
         dmin = minval(F(i,:), F(i,:) > dmin) 
         s%dbar = s%dbar + sum(F(i,:), F(i,:) == dmin)
         npoints = npoints + count(F(i,:) == dmin)
         nd = nd + count(F(i,:) == dmin)
      enddo
   enddo
   s%dbar = s%dbar/nd
   
!  Update F to contain the value of the B-spline function, not just the distance
!  F_ij = f_j(Delta_ij)
   if (.not.silent) call sph_splines_info('sph_splines_fit: Constructing F matrix.')
   do i=1,n
      do j=i,n
         F(i,j) = sph_splines_func(F(i,j),s%dbar)
!         F(j,i) = F(i,j)  ! Only need the upper half in passing to ?posv
      enddo
   enddo
   
!  Solve the set of equations F*a = c  (In LAPACK terms, A*X = B).
!  The knot coefficients, a, are placed in B, so we copy c into a for this.
!  F is replaced by the L and U factors from the factorisation.
   s%a = s%c
   
!  Create array for pivot rows
   allocate(ipiv(n))
!  Call the appropriate LAPACK routine for single or double precision
!  The matrix F is real, symmetric, and the diagonal is 0.  All eigenvalues are
!  positive, so it is positive definite.  (Please tell me if this is wrong!)
!  Hence solving by Cholesky decomposition should be faster than calling the more
!  general ?gesv, which uses LU decomposition.
!  If the input values do not give a positive definite matrix, LAPACK will report
!  ier > 0 and the users has probably entered the wrong knot locations or values.
   if (.not.silent) write(lu_stderr,'(a,i0,a,i0,a)') &
      'spherical_splines: sph_splines_fit: Finding spline coefficients by inverting ',&
      n,' x ',n,' matrix.'
   if (rs == r8) then
!      call dgesv(n,1,F,n,ipiv,s%a,n,ier)  ! Works for all non-singular matrices
      call dposv('U',n,1,F,n,s%a,n,ier)   ! Works for positive definite matrices
   else if (rs == r4) then
!      call sgesv(n,1,F,n,ipiv,s%a,n,ier)
      call sposv('U',n,1,F,n,s%a,n,ier)
   else
      call sph_splines_fatal('module compiled with neither real or double precision')
   endif
   
!  Check status
!   if (ier > 0) call sph_splines_fatal( &
!      'sph_splines_fit: trying to solve singular matrix.  Check knots and values.')
   if (ier > 0) call sph_splines_fatal( &
      'sph_splines_fit: F matrix is not positive definite.  Check knots and values.')
   if (ier < 0) call sph_splines_fatal( &
      'sph_splines_fit: illegal value somewhere in call to LAPACK.')
   if (ier == 0 .and..not.silent) call sph_splines_info( &
      'sph_splines_fit: Successfully solved for knot coefficients.')

!  Deallocate things we won't use outside the subroutine
   deallocate(F,ipiv)
   
!  Mark down which kernel we've fitted
   s%kernel = 'B '

end subroutine sph_splines_fit
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sph_splines_fit_AP(s,half_width,quiet)
!===============================================================================
!  Calculate the spline coefficients for the Abel-Poisson kernel representation
!  h is the desired scaling parameter.  Using sph_splines_AP_width2h can give you
!  a value of h for the desired half-width in degrees of the kernel.
   implicit none
   type(sph_splines_set),intent(inout) :: s
   real(rs),intent(in),optional :: half_width  ! Degrees
   logical,intent(in),optional :: quiet
   logical :: silent
   integer :: n  ! number of knots; length of phi,theta,r
   integer :: i,j,ier
   real(rs),allocatable :: F(:,:)   ! Matrix of distances between point i and j
!   real(rs) :: d,dmin

!  Get the level of verbosity.  Default to noisy
   silent = .false.
   if (present(quiet)) silent = quiet

!  Check for conflict between supplied h or requested half_width
   if (s%h > 0._rs .and. present(half_width) .and. .not.silent) &
      call sph_splines_info('sph_splines_fit_AP: overwriting value of h with that ' &
                         // 'generated by call to sph_splines_AP_width2h')

!  Calculate the desired value of h if we don't know it already
   if (present(half_width)) then
      s%h = sph_splines_AP_width2h(half_width,degrees=.true.)
      if (.not.silent) write(lu_stderr,'(2a,f0.3)') 'spherical_splines: ', &
            'sph_splines_fit_AP: Using h-value of ',s%h
   endif

!  Check for the right things being present
   if (.not.(allocated(s%phi) .and. allocated(s%theta) .and. allocated(s%c))) &
      call sph_splines_fatal('sph_splines_fit_AP: location or values of knot points not allocated')
   
!  Get length of input arrays and make sure they're the same length
   n = size(s%phi)
   if (s%n /= n .or. size(s%theta) /= n .or. size(s%c) /= n) &
      call sph_splines_fatal('sph_splines_fit_AP: structure arrays must all be same length.')
      
!  (Re-)Allocate s%a if already allocated to a length != n
   if (allocated(s%a)) then
      if (size(s%a) /= n) then
         deallocate(s%a);  allocate(s%a(n))
         if (.not.silent) call sph_splines_info( &
               'sph_splines_fit_AP: allocating memory for spline coefficients.')
      endif
   else
      allocate(s%a(n))
      if (.not.silent) call sph_splines_info( &
            'sph_splines_fit_AP: allocating memory for spline coefficients.')
   endif
   
!  Construct matrix of distances between knots, F_ij = Delta_ij
   if (.not.silent) call sph_splines_info( &
      'sph_splines_fit_AP: Calculating inter-knot spacings.')
   allocate(F(n,n))
   do i=1,n
      do j=i,n
         if (i==j) then
            F(i,j) = 0._rs
         else
            F(i,j) = sph_splines_delta(s%phi(i),s%theta(i),s%phi(j),s%theta(j),degrees=s%degrees)
         endif
         F(j,i) = F(i,j)
      enddo
   enddo
   
!  Work out average distance between adjacent points globally, dbar
!   if (.not.silent) call sph_splines_info( &
!      'sph_splines_fit: Calculating mean knot spacing.')
!      
!   s%dbar = 0._rs
!   nd = 0
!   do i=1,n
!  Take the nearest 6 or so points and take the average distance of these
!      dmin = 0._rs
!      npoints = 0
!      do while (npoints <= 6)
!         dmin = minval(F(i,:), F(i,:) > dmin) 
!         s%dbar = s%dbar + sum(F(i,:), F(i,:) == dmin)
!         npoints = npoints + count(F(i,:) == dmin)
!         nd = nd + count(F(i,:) == dmin)
!      enddo
!   enddo
!   s%dbar = s%dbar/nd
   
!  If necessary, convert distances into radians, as assumed by sph_splines_func_AP
   if (s%degrees) F = F*pi/180_rs
   
!  Update F to contain the value of the AP-spline function, not just the distance
!  F_ij = f_j(Delta_ij)
   if (.not.silent) call sph_splines_info('sph_splines_fit_AP: Constructing F matrix.')
   do i=1,n
      do j=i,n
         F(i,j) = sph_splines_func_AP(F(i,j),s%h)
!         F(j,i) = F(i,j)  ! Only need the upper half in passing to ?posv
      enddo
   enddo
   
!  Solve the set of equations F*a = c  (In LAPACK terms, A*X = B).
!  The knot coefficients, a, are placed in B, so we copy c into a for this.
!  F is replaced by the L and U factors from the factorisation.
   s%a = s%c
   
!  Call the appropriate LAPACK routine for single or double precision
!  The matrix F is real, symmetric, and the diagonal is 0.  All eigenvalues are
!  positive, so it is positive definite.  (Please tell me if this is wrong!)
!  Hence solving by Cholesky decomposition should be faster than calling the more
!  general ?gesv, which uses LU decomposition.
!  If the input values do not give a positive definite matrix, LAPACK will report
!  ier > 0 and the users has probably entered the wrong knot locations or values.
   if (.not.silent) write(lu_stderr,'(a,i0,a,i0,a)') &
      'spherical_splines: sph_splines_fit_AP: Finding spline coefficients by inverting ',&
      n,' x ',n,' matrix.'
   if (rs == r8) then
      call dposv('U',n,1,F,n,s%a,n,ier)   ! Works for positive definite matrices
   else if (rs == r4) then
      call sposv('U',n,1,F,n,s%a,n,ier)
   else
      call sph_splines_fatal('module compiled with neither real or double precision')
   endif
   
!  Check status
!   if (ier > 0) call sph_splines_fatal( &
!      'sph_splines_fit: trying to solve singular matrix.  Check knots and values.')
   if (ier > 0) call sph_splines_fatal( &
      'sph_splines_fit_AP: F matrix is not positive definite.  Check knots and values.')
   if (ier < 0) call sph_splines_fatal( &
      'sph_splines_fit_AP: illegal value somewhere in call to LAPACK.')
   if (ier == 0 .and..not.silent) call sph_splines_info( &
      'sph_splines_fit_AP: Successfully solved for knot coefficients.')

!  Mark down which kernel we've fitted
   s%kernel = 'AP'

end subroutine sph_splines_fit_AP
!-------------------------------------------------------------------------------
   

!===============================================================================
subroutine sph_splines_eval(phi,theta,s,c,degrees)
!===============================================================================
!  Evaluate a set of splines with coefficients at a set of arbitrary points
!  on a sphere
   implicit none
   type(sph_splines_set),intent(in) :: s
   real(rs),intent(in),dimension(:) :: phi,theta
   real(rs),dimension(:),intent(out) :: c
   logical,optional,intent(in) :: degrees
   real(rs) :: conversion,d,torad
   integer :: n  ! length of phi and theta
   integer :: i,j
   
!  Get size of arrays and check a few things
   n = size(phi)
   if (size(theta) /= n .or. size(c) /= n) call sph_splines_fatal( &
      'sph_splines_eval: evaluation coordinates and values must be in arrays of same length')
   if (.not.(allocated(s%phi) .and. allocated(s%theta) .and. allocated(s%a) &
       .and. allocated(s%c))) call sph_splines_fatal( &
      'sph_splines_eval: one or more arrays in the supplied spline set not allocated.')
   if (.not.any(KERNEL_NAME == s%kernel)) call sph_splines_fatal( &
      'sph_splines_eval: kernel type "' // s%kernel // '" not supported.')
   
!  If necessary, convert input phi and theta to the same units as the
!  splines are given in.  Assume they are the same
   conversion = 1._rs
   if (present(degrees)) then
      if (degrees .eqv. s%degrees) then
         conversion = 1._rs
      else if (degrees) then      ! Must mean s%degrees = .false.
         conversion = pi/180._rs
      else if (.not.degrees) then ! Must mean s%degrees = .true.
         conversion = 180._rs/pi
      else
         call sph_splines_fatal('sph_splines_eval: Andy''s coding is wrong.')
      endif
   else  ! Assume input is in degrees
      if (s%degrees) conversion = 1._rs
      if (.not. s%degrees) conversion = pi/180._rs
   endif
   
!  We need to pass radians to the AP spline function, so one final conversion.
!  The B-sline function merely needs d and dbar in the same units, so no need to 
!  change units for this function.
   torad = 1._rs
   if (s%degrees .and. s%kernel=="AP") torad = pi/180._rs
         
!  Evaluate spline at each point
   do i=1,n
!      c(i) = sum(s%a*sph_splines_func( &
!               sph_splines_delta(phi(i),theta(i),s%phi,s%theta,degrees=s%degrees),&
!               s%dbar)))
      c(i) = 0._rs
      do j=1,s%n
         d = torad * sph_splines_delta(conversion*phi(i),conversion*theta(i),&
                               s%phi(j),s%theta(j),degrees=s%degrees)
         if (s%kernel == 'B ') then
            c(i) = c(i) + sph_splines_func(d,s%dbar)*s%a(j)
         else if (s%kernel == 'AP') then
            c(i) = c(i) + sph_splines_func_AP(d,s%h)*s%a(j)
         else
            call sph_splines_fatal('sph_splines_fit: kernel "' // s%kernel // &
               '" is not recognised.')
         endif
      enddo
   enddo
   
end subroutine sph_splines_eval
!-------------------------------------------------------------------------------

!===============================================================================
function sph_splines_func(d,dbar) result(f)
!===============================================================================
!  Compute the spherical B-splines function:
!         { (3/4)d^3/(dbar^3) - (6/4)d^2/(dbar^2) + 1,         d <= dbar
!     f = {-(1/4)d1^3 + (3/4)d1^2 - (3/4)d1 + (1/4)  , dbar <= d <= 2dbar
!         { 0                                        ,         d >  2dbar
!  where
!     d1 = (d - dbar)/dbar.

   implicit none
   real(rs),intent(in)  :: d    ! Distance between this point and knot point
   real(rs),intent(in)  :: dbar ! Mean knot spacing
   real(rs) :: f
   real(rs) :: d1
   
!  Where the distance is zero, the value of the function is 1.
   if (d < 0._rs) then
      call sph_splines_fatal('sph_splines_func: distance must be > 0.')
   else if (d == 0.) then
      f = 1._rs
   else if (d <= dbar) then
      f = (3._rs/4._rs)*d**3/dbar**3 - (6._rs/4._rs)*d**2/dbar**2 + 1._rs
   else if (d > dbar .and. d <= 2._rs*dbar) then
      d1 = (d - dbar)/dbar
      f = -(1._rs/4._rs)*d1**3 + (3._rs/4._rs)*d1**2 - &
          (3._rs/4._rs)*d1 + (1._rs/4._rs)
   else  ! Outside 2dbar, the value and its derivatives are zero.
      f = 0._rs
   endif
      
end function sph_splines_func
!-------------------------------------------------------------------------------

!===============================================================================
function sph_splines_func_AP(d,h) result(K)
!===============================================================================
!  Compute the Abel-Poisson spherical splines kernel:
!     Kh(x,y) = (1/4*pi)*(1 - h**2)/(1 + h**2 - 2*h*dot_product(x,y))**(3/2),
!  where 
!     x and y are points on the unit sphere.  (Hence x.y is the cosine of the 
!  angular distance (d) between the two points.)
!     h is the scaling factor, 0<=h<=1.  Values of h=0.96 correspond to a
!  half-width of maximum intensity at ~1-h (radians)

   implicit none
   real(rs),intent(in) :: d,h  ! Angle between points (radians); and scaling
   real(rs) :: K
   
   if (d < 0._rs) &
      call sph_splines_fatal('sph_splines_func_AP: distance must be > 0.')
   K = (1._rs/(4._rs*pi))*(1._rs - h**2)/(1._rs + h**2 - 2._rs*h*cos(d))**(1.5_rs)
   
end function sph_splines_func_AP
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sph_splines_write_coefs(s,outfile)
!===============================================================================
!  Write out the coefficients for a spherical splines set.
!  
!  Format is (ascii):
!     [1]      n kernel g     ! n = num coefficients 
!                             ! kernel = type of kernel
!                             ! g is dbar or h, depending on kernel
!     [2..n+1] lon lat c      ! lon, colat (degrees), value and coefficient

   implicit none
   type(sph_splines_set),intent(in) :: s
   character(len=*), intent(in) :: outfile
   real(rs) :: temp_val,conversion
   logical :: opened
   integer :: unit,i,iostatus
   
   call sph_splines_fatal('sph_splines_write_ceofs: Subroutine not tested yet!')
   
!  Look for available logical unit
   unit = 10
   opened = .true.
   do while(opened)
      inquire(unit,opened=opened)
      if (.not.opened) exit
      if (unit == 100) then
         call sph_splines_fatal( &
            'sph_splines_write_coefs: Unable to find available unit for writing.')
      endif
      unit = unit + 1
   enddo
   
!  Open file
   open(unit,file=trim(outfile),iostat=iostatus)
   if (iostatus /= 0) call sph_splines_fatal( &
      'sph_splines_write_coefs: Problem opening file "' // trim(outfile) // '"')
   
!  Write header line
   if (s%kernel == 'B') then
      temp_val = s%dbar
      if (.not.s%degrees) temp_val = s%dbar * 180._rs/pi
   elseif (s%kernel == 'AP') then
      temp_val = s%h
   else
      call sph_splines_fatal( &
         'sph_splines_write_coefs: Unknown kernel "' // s%kernel // '"')
   endif
   
!  If necessary, convert phi and theta to degrees for output
   conversion = 1._rs
   if (.not.s%degrees) conversion = 180._rs/pi
!  Write header
   write(unit,*) s%n, s%kernel, temp_val
!  Write knot points and coefficients
   do i=1,s%n
      write(unit,*) conversion*s%phi(i), conversion*s%theta(i), s%c(i), s%a(i)
   enddo
   
   close(unit)

end subroutine sph_splines_write_coefs
!-------------------------------------------------------------------------------   

!===============================================================================
function sph_splines_read_coefs(infile) result(s)
!===============================================================================
!  Read the coefficients for a spherical splines set
!  For format of files, see sph_splines_write_coefs

   implicit none
   type(sph_splines_set) :: s
   character(len=*),intent(in) :: infile
   real(rs) :: temp_val
   integer :: i,iostatus,unit
   logical :: opened
   
   call sph_splines_fatal('sph_splines_read_coefs: Subroutine not tested yet!')
   
!  Look for available logical unit
   unit = 10
   opened = .true.
   do while(opened)
      inquire(unit,opened=opened)
      if (.not.opened) exit
      if (unit == 100) then
         call sph_splines_fatal( &
            'sph_splines_read_coefs: Unable to find available unit for reading.')
      endif
      unit = unit + 1
   enddo

!  Open file
   open(unit,file=trim(infile),iostat=iostatus)
   if (iostatus /= 0) call sph_splines_fatal( &
      'sph_splines_read_coefs: Problem opening file "' // trim(infile) // '"')
   
!  Read header
   read(unit,*) s%n, s%kernel, temp_val
   if (s%kernel == 'B') then
      s%dbar = temp_val
   elseif (s%kernel == 'AP') then
      s%h = temp_val
   else
      call sph_splines_fatal( &
         'sph_splines_read_coefs: Unknown kernel "' // trim(s%kernel) // '"')
   endif
   
!  Allocate space for coordinates and coefficients
   call sph_splines_realloc(s%n,s)

!  Read in data
   do i=1,s%n
      read(unit,*) s%phi(i), s%theta(i), s%c(i), s%a(i)
   enddo
   
   close(unit)
   
!  Set the units of the coefficients to be degrees
   s%degrees = .true.
   
end function sph_splines_read_coefs
!-------------------------------------------------------------------------------

!===============================================================================
function sph_splines_AP_width2h(x,degrees) result(h)
!===============================================================================
!  Using an empirical power law, return the value of h which produces a
!  Abel-Poisson kernel which has the requested half-width.
!  Optionally, specify desired half-width in degrees; default is radians.
   implicit none
   real(rs),intent(in) :: x
   logical,intent(in),optional :: degrees
   real(rs) :: h,r
   real(rs),parameter :: a = -1.31743_rs, b = 0.89017_rs, c = -0.315546_rs, &
                         d = 0.0468784_rs
   
   r = x  ! Input in radians
   if (present(degrees)) then
      if (degrees) r = pi*x/180._rs  ! Input in degrees
   endif
   
   h = 1._rs + a*r + b*r**2 + c*r**3 + d*r**4

end function sph_splines_AP_width2h
!-------------------------------------------------------------------------------

!===============================================================================
function sph_splines_delta(phi1,theta1,phi2,theta2,degrees) result(delta)
!===============================================================================
!  delta returns the angular distance between two points on a sphere given 
!  the colat and lon of each using the Haversine formula

   implicit none
   real(rs),intent(in) :: phi1,theta1,phi2,theta2
   real(rs) :: delta
   real(rs) :: lon1,lon2,lat1,lat2,conversion
   logical,optional,intent(in) :: degrees
   
!   real(rs) :: cos_lat1,cos_lat2,sin_lat1,sin_lat2,&
!               cos_lon1,cos_lon2,sin_lon1,sin_lon2
      
!  Default to degrees input
   conversion = pi/180._rs
   if (present(degrees)) then
      if (.not.degrees) conversion = 1._rs
   endif
   
!  Convert to radians and from colat to lat for calculation
   lon1 = phi1*conversion         ;  lon2 = phi2*conversion
   lat1 = pi2 - theta1*conversion ;  lat2 = pi2 - theta2*conversion
   
   delta = atan2( sqrt( (cos(lat2)*sin(lon2-lon1))**2 + (cos(lat1)*sin(lat2) - &
           sin(lat1)*cos(lat2)*cos(lon2-lon1))**2) , &
           sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1))

 
!  Convert back to desired units
   delta = delta/conversion
  
end function sph_splines_delta
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sph_splines_fatal(string)
!===============================================================================
!  Utility to flag up a fatal error which stops the program.  Input should be:
!  <name of subroutine>: <error message>
   implicit none
   character(len=*) :: string
   
   write(lu_stderr,'(2a)') "spherical_splines: ",trim(string)
!  Assuming compiler compatibility, exit with non-0 status
   stop 1
end subroutine sph_splines_fatal
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sph_splines_info(string)
!===============================================================================
!  Utility to write out some information to stderr.  Input should be:
!  <name of subroutine>: <info message>
   implicit none
   character(len=*) :: string
   
   write(lu_stderr,'(2a)') "spherical_splines: ",trim(string)
end subroutine sph_splines_info
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sph_splines_realloc(n,s)
!===============================================================================
!  (Re)allocate the right amount of space for the arrays in the sph_splines_set
!  structure.
!  Defaults to verbose, so supply quiet=.true. to hush the subroutine up.

   implicit none
   type(sph_splines_set),intent(inout) :: s
   integer,intent(in) :: n
   
   if (allocated(s%phi)) then
      if (size(s%phi) /= n) then
         deallocate(s%phi);  allocate(s%phi(n))
      endif
   else
      allocate(s%phi(n))
   endif

   if (allocated(s%theta)) then
      if (size(s%theta) /= n) then
         deallocate(s%theta);  allocate(s%theta(n))
      endif
   else
      allocate(s%theta(n))
   endif

   if (allocated(s%a)) then
      if (size(s%a) /= n) then
         deallocate(s%a);  allocate(s%a(n))
      endif
   else
      allocate(s%a(n))
   endif

   if (allocated(s%c)) then
      if (size(s%c) /= n) then
         deallocate(s%c);  allocate(s%c(n))
      endif
   else
      allocate(s%c(n))
   endif

end subroutine sph_splines_realloc
!-------------------------------------------------------------------------------

   
!_______________________________________________________________________________
end module spherical_splines
