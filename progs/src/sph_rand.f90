!===============================================================================
program sph_rand
!===============================================================================
! Randomly sample a sphere.
! This program takes samples one at a time and writes them out as it goes,
! rather than allocating memory for all points at once, to save memory in case
! of very large numbers of points.
! Output is lon, lat in degrees.

   use spherical_geometry

   implicit none
   integer :: n = 100
   real(8) :: lon, lat
   character(len=250) :: arg
   integer :: i, iostat

   ! Check args and get n if necessary
   if (command_argument_count() > 1) call usage
   if (command_argument_count() == 1) then
      call get_command_argument(1, arg)
      if (arg(1:2) == '-h') call usage
      read(arg, *, iostat=iostat) n
      if (iostat /= 0) then
         write(0,'(a)') 'sph_rand: Error: Cannot get n from first argument "' &
            //trim(arg)//'"'
         call usage
      endif
   endif

   do i = 1, n
      call sg_random_point_geog(lon, lat, degrees=.true.)
      write(*,*) lon, lat
   enddo

contains

   subroutine usage
      write(0,'(a)') &
         'Usage: sph_rand (n)', &
         'Randomly take points on a sphere with uniform probability density.', &
         'If specified, first argument defines the number of points.', &
         'Output is lon, lat (degrees).'
      stop
   end subroutine usage
end program sph_rand
!-------------------------------------------------------------------------------
