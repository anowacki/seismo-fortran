!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program sph_tess
!===============================================================================
! Dump the points for a level-n tesselation of a sphere, starting with an
! icosahedron.
! This version will create a tesselation of higher order if it is not available
! as a pre-saved version on disk.

   use sphere_tesselate
   use spherical_geometry

   implicit none
   type(tesselation) :: t
   integer :: i, n
   logical :: output_triangles = .false., output_geog = .false., &
      cache_exists = .false.
   ! Path to directory containing pre-computed tesselations
   character(len=250) :: file

   call get_args

   ! Try to load the file from the global cache
   call st_load_cache(t, cache_exists)
   ! Otherwise, create the tesselation and save the cache file for later use
   if (.not.cache_exists) then
      write(0,'(a)') 'Generating tesselation and saving cache file...'
      t = st_icosahedron(pole=.true.)
      do i = 1, n
         call st_iterate_level(t)
      enddo
      call st_save_cache(t)
   endif
   ! Output
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
         '   -g        : Write values as lon, lat instead of x, y, z', &
         '   -p [path] : Set directory path to precomputed tesselations', &
         '   -t        : Dump the triangles as multisegment lines'
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
            case('-g');   output_geog = .true.;       iarg = iarg + 1
            case('-p');   call get_command_argument(iarg+1, st_cache_dir)
                                                      iarg = iarg + 2
            case('-t');   output_triangles = .true.;  iarg = iarg + 1
            case default; call usage
         end select
      enddo
      ! Get tesselation level
      call get_command_argument(narg, arg)
      read(arg,*) n
      t%level = n
      if (n < 0) call usage
   end subroutine get_args

end program sph_tess
!-------------------------------------------------------------------------------
