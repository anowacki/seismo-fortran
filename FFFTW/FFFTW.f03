!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
!
! FFFTW - Fastest (Fortran) Fourier Transform in the West
!
!===============================================================================
!
! FFFTW is a module interface to the FFTW3 C DFT library.  It provides
! a simple wrapper to call FFTW's routines from Fortran programs.  It requires
! a F2003-aware compiler in order to use the ISO_C_BINDING intrinsic module,
! but there are no requirements other than that.
!
! The module requires FFTW3 to be installed and probably compiled with the same
! type of compiler (e.g., GNU, Intel, etc.) as used to compile the module.
!
! Module procedures allow one to call just a forward or reverse FFT: just use
! the array you want to be transformed and the module will automatically choose
! the correct FFTW routine.
!
! E.g., to transform a 1D single-precision real array:
!     real :: time_domain_array(100)
!     complex, allocatable :: freq_domain_array(:)
!     call FFFTW_fwd(time_domain_array, freq_domain_array)
!
! To then transform back:
!     call FFFTW_rev(freq_domain_array, time_domain_array)
!
! (N.B.: freq_domain_array is destroyed in the process; to save it, copy the
! array to another variable before the reverse transform.)
!
! FFFTW computes the required size of the FFTed complex trace automatically and
! (re)allocates it as necessary.  Hence the FFT array must be declared
! allocatable.  The compiler should flag any attempt to do otherwise, but note
! that the error message may not say explicitly so--instead, you may get an
! undefined reference to the generic overloaded procedure FFFTW_{fwd,rev}.
!===============================================================================
module FFFTW
!===============================================================================

   use, intrinsic :: iso_c_binding

   implicit none

   ! Only allow explicitly public variables and routines to be visible
!    private

   include 'fftw3.f03'

   ! Precision selectors
   integer, parameter, private :: r4 = selected_real_kind(6,37), &
                                  r8 = selected_real_kind(15,307)
   integer, parameter, private :: c4 = r4, &
                                  c8 = r8

   ! Current state of FFFTW: The first time we're called, make sure that the
   ! Fortran and C types match.  Otherwise we're doomed.
   logical, private :: FFFTW_initialised = .false.

   ! Can call each type of transform with the same call regardless of using
   ! real(4) or real(8) precision.
   interface FFFTW_fwd
      module procedure :: FFFTW_fwd_r2c_1d_real4, &
                          FFFTW_fwd_r2c_1d_real8, &
                          FFFTW_fwd_c2c_1d_complex4, &
                          FFFTW_fwd_c2c_1d_complex8
   end interface FFFTW_fwd

   interface FFFTW_rev
      module procedure :: FFFTW_rev_c2r_1d_real4, &
                          FFFTW_rev_c2r_1d_real8, &
                          FFFTW_rev_c2c_1d_complex4, &
                          FFFTW_rev_c2c_1d_complex8
   end interface FFFTW_rev

   interface FFFTW_allocate
      module procedure :: FFFTW_allocate_real4, &
                          FFFTW_allocate_real8, &
                          FFFTW_allocate_complex4, &
                          FFFTW_allocate_complex8
   end interface FFFTW_allocate

   private :: FFFTW_allocate_real4, &
              FFFTW_allocate_real8, &
              FFFTW_allocate_complex4, &
              FFFTW_allocate_complex8, &
              FFFTW_err

contains

!===============================================================================
subroutine FFFTW_fwd_r2c_1d_real4(real_array, complex_fft, plan_type)
!===============================================================================
!  FFFTW_fwd_*() computes the FT for a 1d array.  This is done out-of-place, so
!  the input array is preserved.  complex_fft
   implicit none
   real(r4), intent(in), dimension(:) :: real_array
   complex(c4), intent(inout), allocatable, dimension(:) :: complex_fft
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc

   call FFFTW_initialise

   ! Size of real and complex FFT traces
   nr = size(real_array)
   nc = nr/2 + 1

   ! Allocate memory for the FFT trace if necessary
   call FFFTW_allocate(complex_fft, nc)

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call sfftw_plan_dft_r2c_1d(plan, &
                              int(nr, kind=C_INT), &
                              real(real_array, kind=r4), &
                              complex_fft, &
                              fftw_plan_type)
   call sfftw_execute_dft_r2c(plan, real_array, complex_fft)
   call sfftw_destroy_plan(plan)

end subroutine FFFTW_fwd_r2c_1d_real4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_fwd_r2c_1d_real8(real_array, complex_fft, plan_type)
!===============================================================================
   implicit none
   real(r8), intent(in), dimension(:) :: real_array
   complex(c8), intent(inout), allocatable, dimension(:) :: complex_fft
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc

   call FFFTW_initialise

   ! Size of real and complex FFT traces
   nr = size(real_array)
   nc = nr/2 + 1

   ! Allocate memory for the FFT trace if necessary
   call FFFTW_allocate(complex_fft, nc)

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call dfftw_plan_dft_r2c_1d(plan, &
                              int(nr, kind=C_INT), &
                              real(real_array, kind=r8), &
                              complex_fft, &
                              fftw_plan_type)
   call dfftw_execute_dft_r2c(plan, real_array, complex_fft)
   call dfftw_destroy_plan(plan)

end subroutine FFFTW_fwd_r2c_1d_real8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_fwd_c2c_1d_complex4(complex_array, complex_fft, plan_type)
!===============================================================================
   implicit none
   complex(c4), intent(in), dimension(:) :: complex_array
   complex(c4), intent(inout), allocatable, dimension(:) :: complex_fft
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc

   call FFFTW_initialise

   ! Size of real and complex FFT traces
   nr = size(complex_array)
   nc = nr

   ! Allocate memory for the FFT trace if necessary
   call FFFTW_allocate(complex_fft, nc)

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call sfftw_plan_dft_1d(plan, &
                          int(nr, kind=C_INT), &
                          complex_array, &
                          complex_fft, &
                          FFTW_FORWARD, &
                          fftw_plan_type)
   call sfftw_execute_dft(plan, complex_array, complex_fft)
   call sfftw_destroy_plan(plan)

end subroutine FFFTW_fwd_c2c_1d_complex4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_fwd_c2c_1d_complex8(complex_array, complex_fft, plan_type)
!===============================================================================
   implicit none
   complex(c8), intent(in), dimension(:) :: complex_array
   complex(c8), intent(inout), allocatable, dimension(:) :: complex_fft
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc

   call FFFTW_initialise

   ! Size of real and complex FFT traces
   nr = size(complex_array)
   nc = nr

   ! Allocate memory for the FFT trace if necessary
   call FFFTW_allocate(complex_fft, nc)

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call dfftw_plan_dft_1d(plan, &
                          int(nr, kind=C_INT), &
                          complex_array, &
                          complex_fft, &
                          FFTW_FORWARD, &
                          fftw_plan_type)
   call dfftw_execute_dft(plan, complex_array, complex_fft)
   call dfftw_destroy_plan(plan)

end subroutine FFFTW_fwd_c2c_1d_complex8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_rev_c2r_1d_real4(complex_fft, real_array, plan_type)
!===============================================================================
   implicit none
   complex(c4), intent(inout), allocatable, dimension(:) :: complex_fft
   real(r4), intent(inout), dimension(:) :: real_array
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc

   call FFFTW_initialise

   ! Check array sizes
   nc = size(complex_fft)
   nr = size(real_array)
   if (nc /= nr/2 + 1) &
      call FFFTW_err('FFFTW_rev_c2r_1d_real4: Array sizes of FFT and T-domain' &
         // ' trace not as expected')

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call sfftw_plan_dft_c2r_1d(plan, &
                              int(nr, kind=C_INT), &
                              complex_fft, &
                              real_array, &
                              fftw_plan_type)
   call sfftw_execute_dft_c2r(plan, complex_fft, real_array)
   call sfftw_destroy_plan(plan)

   ! Normalise values
   real_array = real_array/real(nr, kind=r4)

end subroutine FFFTW_rev_c2r_1d_real4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_rev_c2r_1d_real8(complex_fft, real_array, plan_type)
!===============================================================================
   implicit none
   complex(c8), intent(inout), allocatable, dimension(:) :: complex_fft
   real(r8), intent(inout), dimension(:) :: real_array
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc

   call FFFTW_initialise

   ! Check array sizes
   nc = size(complex_fft)
   nr = size(real_array)
   if (nc /= nr/2 + 1) &
      call FFFTW_err('FFFTW_rev_c2r_1d_real8: Array sizes of FFT and T-domain' &
         // ' trace not as expected')

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call dfftw_plan_dft_c2r_1d(plan, &
                              int(nr, kind=C_INT), &
                              complex_fft, &
                              real_array, &
                              fftw_plan_type)
   call dfftw_execute_dft_c2r(plan, complex_fft, real_array)
   call dfftw_destroy_plan(plan)

   ! Normalise values
   real_array = real_array/real(nr, kind=r4)

end subroutine FFFTW_rev_c2r_1d_real8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_rev_c2c_1d_complex4(complex_fft, complex_array, plan_type)
!===============================================================================
   implicit none
   complex(c4), intent(inout), allocatable, dimension(:) :: complex_fft
   complex(c4), intent(inout), dimension(:) :: complex_array
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc
   complex(c4) :: norm

   call FFFTW_initialise

   ! Check array sizes
   nc = size(complex_fft)
   nr = size(complex_array)
   if (nc /= nr) &
      call FFFTW_err('FFFTW_rev_c2c_1d_complex4: Array sizes of FFT and T-domain' &
         // ' trace not as expected')

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call sfftw_plan_dft_1d(plan, &
                          int(nr, kind=C_INT), &
                          complex_fft, &
                          complex_array, &
                          FFTW_BACKWARD, &
                          fftw_plan_type)
   call sfftw_execute_dft(plan, complex_fft, complex_array)
   call sfftw_destroy_plan(plan)

   ! Normalise values
   norm = real(nr, kind=r4)
   complex_array = complex_array/norm

end subroutine FFFTW_rev_c2c_1d_complex4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_rev_c2c_1d_complex8(complex_fft, complex_array, plan_type)
!===============================================================================
   implicit none
   complex(c8), intent(inout), allocatable, dimension(:) :: complex_fft
   complex(c8), intent(inout), dimension(:) :: complex_array
   integer(C_INT), intent(in), optional :: plan_type
   integer(C_INT) :: fftw_plan_type
   type(C_PTR) :: plan
   integer :: nr, nc
   complex(c8) :: norm

   call FFFTW_initialise

   ! Check array sizes
   nc = size(complex_fft)
   nr = size(complex_array)
   if (nc /= nr) &
      call FFFTW_err('FFFTW_rev_c2c_1d_complex4: Array sizes of FFT and T-domain' &
         // ' trace not as expected')

   fftw_plan_type = FFTW_ESTIMATE
   if (present(plan_type)) fftw_plan_type = int(plan_type, kind=C_INT)

   call dfftw_plan_dft_1d(plan, &
                          int(nr, kind=C_INT), &
                          complex_fft, &
                          complex_array, &
                          FFTW_BACKWARD, &
                          fftw_plan_type)
   call dfftw_execute_dft(plan, complex_fft, complex_array)
   call dfftw_destroy_plan(plan)

   ! Normalise values
   norm = real(nr, kind=r8)
   complex_array = complex_array/norm

end subroutine FFFTW_rev_c2c_1d_complex8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_allocate_real4(a,n)
!===============================================================================
!  FFFTW_allocate_*() allocates memory to an array.  If that array is already
!  allocated but is of different length to n, it is reallocated.  If it is the
!  same length as n, it is not touched.
   implicit none
   real(r4), dimension(:), allocatable, intent(inout) :: a
   integer, intent(in) :: n
   if (allocated(a)) then
      if (size(a) /= n) allocate(a(n))
   else
      allocate(a(n))
   endif
end subroutine FFFTW_allocate_real4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_allocate_real8(a,n)
!===============================================================================
   implicit none
   real(r8), dimension(:), allocatable, intent(inout) :: a
   integer, intent(in) :: n
   if (allocated(a)) then
      if (size(a) /= n) allocate(a(n))
   else
      allocate(a(n))
   endif
end subroutine FFFTW_allocate_real8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_allocate_complex4(a,n)
!===============================================================================
   implicit none
   complex(c4), dimension(:), allocatable, intent(inout) :: a
   integer, intent(in) :: n
   if (allocated(a)) then
      if (size(a) /= n) allocate(a(n))
   else
      allocate(a(n))
   endif
end subroutine FFFTW_allocate_complex4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_allocate_complex8(a,n)
!===============================================================================
   implicit none
   complex(c8), dimension(:), allocatable, intent(inout) :: a
   integer, intent(in) :: n
   if (allocated(a)) then
      if (size(a) /= n) allocate(a(n))
   else
      allocate(a(n))
   endif
end subroutine FFFTW_allocate_complex8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_initialise()
!===============================================================================
!  Perform initialisation procedures.
   implicit none
   if (.not.FFFTW_initialised) then
      ! Test for equivalence of Fortran selected kinds and C kinds.
      ! The Fortran kinds are used by the compiler to pick the correct
      ! subroutines when calling the generic procedures (FFFTW_fwd(), etc.),
      ! therefore as long as they are equal to the C kinds we can pick the
      ! correct subroutine and call the FFTW functions safely.
      if (r4 /= C_FLOAT .or. r8 /= C_DOUBLE .or. c4 /= C_FLOAT_COMPLEX .or. &
          c8 /= C_DOUBLE_COMPLEX) &
          call FFFTW_err('FFFTW_initialise: Fortran kinds do not match C types')
      FFFTW_initialised = .true.
   endif
end subroutine FFFTW_initialise
!-------------------------------------------------------------------------------

!===============================================================================
subroutine FFFTW_err(string)
!===============================================================================
!  Private utility routine to write a message and stop the program
   implicit none
   character(len=*), intent(in) :: string
   write(0,'(a)') 'FFFTW: Error: ' // trim(string)
   stop
end subroutine FFFTW_err
!-------------------------------------------------------------------------------

end module FFFTW
!-------------------------------------------------------------------------------
