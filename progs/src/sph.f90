!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program sph
!  Evenly sample a sphere.  Supply the spacing between points as the first argument,
!  otherwise defaults to 15 degrees.  Values nearly divisible into 90 are best,
!  or small ones.

   use spherical_geometry, only: sphere_sample

   implicit none

   real(8), allocatable, dimension(:) :: lon, lat
   real(8) :: d
   integer :: i, iostat, n
   character(len=80) :: arg

   if (command_argument_count() == 1) then
      call getarg(1,arg)
      read(arg,*,iostat=iostat) d
      if (iostat /= 0) call usage
   else if (command_argument_count() == 0) then
      d = 15._8
   else
      call usage
   endif

   call sphere_sample(d, lon, lat, n)
   do i=1,n
      write(*,*) lon(i), lat(i)
   enddo

contains
   !============================================================================
   subroutine usage
   !============================================================================
      write(0,'(a)') 'Usage: sph (point spacing in degrees)', &
                     'Writes a set of evenly-spaced (lon,lat) points on a sphere to stdout'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------
end program sph