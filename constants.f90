!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
!  Module containting some physical and precision constants.
!===============================================================================
   module constants
!===============================================================================
      implicit none

!  ** size constants
      integer, parameter :: i4 = selected_int_kind(9) ; ! long int
      integer, parameter :: r4 = selected_real_kind(6,37) ; ! SP
      integer, parameter :: r8 = selected_real_kind(15,307) ; ! DP

!  ** precision selector
      integer, parameter, public :: rs = r8
      
!  ** maths constants and other useful things
      real(r8), parameter :: pi = 3.141592653589793238462643d0 ;
      real(r8), parameter :: to_rad = 1.74532925199433D-002 ;  
      real(r8), parameter :: to_deg = 57.2957795130823d0 ;  
      real(r8), parameter :: to_km = 111.194926644559 ;      

      real(r8), parameter :: big_number = 10.d36 ;      
      real(r8), parameter :: Earth_radius = 6371.d0 ;
      
!  ** IO direction
      integer, parameter :: lu_stderr = 0
      integer, parameter :: lu_stdin  = 5
      integer, parameter :: lu_stdout = 6
           
   end module constants
!===============================================================================


