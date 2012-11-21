!===============================================================================
!     Fortran 90/95 Source Code File
!===============================================================================
!
!     PROGRAM : useful_constants
!     FILE    : useful_constants.f90
!     AUTHOR  : James Wookey
!     PLACE   : School of Earth Sciences, University of Leeds
!     DATE    : July 2004
!     PURPOSE : 
!     VERSION : 0.1
!     COMPLETE: No
!     COMMENTS: 
!
!-------------------------------------------------------------------------------
!     This software is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!-------------------------------------------------------------------------------
!
!     A simple module of useful constants
!
!-------------------------------------------------------------------------------
!     Changes log
!-------------------------------------------------------------------------------
!     2004-07-21     * Incept date
!     2010-12-06     * Added Earth radius


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


