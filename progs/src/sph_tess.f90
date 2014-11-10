!===============================================================================
program sph_tess
!===============================================================================
! Dump the points for a level-n tesselation of a sphere, starting with an
! icosahedron.

   use sphere_tesselate
   use spherical_geometry

   implicit none
   type(tesselation) :: t
   integer :: i, n
   logical :: output_triangles = .false., output_geog = .false.

   call get_args

   t = st_icosahedron()
   do i = 1, n
      call st_iterate_level(t)
   enddo
   if (.not.output_triangles) then
      call st_dump_points(t, geog=output_geog)
   else
      call st_dump_triangles(t, geog=output_geog)
   endif

contains

   subroutine usage
      write(0,'(a)') &
         'Usage: sph_tess (options) [n]', &
         'Dump points for a level n tesselation of a sphere, starting with an', &
         'icosahedron.  x, y, z values on the unit sphere are output to stdout.', &
         'n must be >= 0', &
         'Options:', &
         '   -g : Write values as lon, lat instead of x, y, z', &
         '   -t : Dump the triangles as multisegment lines'
      stop
   end subroutine usage

   subroutine get_args
      character(len=250) :: arg
      integer :: iarg, narg
      narg = command_argument_count()
      if(narg == 0) call usage
      ! Get options
      iarg = 1
      do while (iarg <= narg - 1)
         call get_command_argument(iarg, arg)
         select case(arg)
            case('-g');     output_geog = .true.;      iarg = iarg + 1
            case('-t');     output_triangles = .true.; iarg = iarg + 1
            case default;   call usage
         end select
      enddo
      ! Get tesselation level
      call get_command_argument(narg, arg)
      read(arg,*) n
      if (n < 0) call usage
   end subroutine get_args

end program sph_tess
!-------------------------------------------------------------------------------
