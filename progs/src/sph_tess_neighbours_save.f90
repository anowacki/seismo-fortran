!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
! Written by Jack Walpole
!===============================================================================
program sph_tess_neighbours_save
!===============================================================================
! Save the neighbour point indices for a level-n tesselation of a sphere, starting
! with an icosahedron, to default cache location.

   use sphere_tesselate

   implicit none
   type(tesselation) :: t
   type(list_neighbours) :: lnbrs
   integer :: i, n

   call get_args

   ! Create the tesselation and save the cache file for later use
   write(0,'(a)') 'Generating tesselation...'
   t = st_icosahedron(pole=.true.)
   do i = 1, n
     call st_iterate_level(t)
   enddo
   write(0,'(a)') 'Generating neighbours cache...'
   call st_generate_list_neighbours(t, lnbrs)
   call st_save_neighbours_cache(lnbrs)


contains

   subroutine usage
      write(0,'(a)') &
         'Usage: sph_tess_neighbours_save [n]', &
         'Create cache of neighbour (point indices) for a level n tesselation of a &
         &sphere, starting with an icosahedron. ', &
         'n must be >= 0'
      stop
   end subroutine usage

   subroutine get_args
      character(len=250) :: arg
      integer :: iarg, narg
      narg = command_argument_count()
      if(narg == 0) call usage
      ! Get tesselation level
      call get_command_argument(narg, arg)
      read(arg,*) n
      t%level = n
      if (n < 0) call usage
   end subroutine get_args

end program sph_tess_neighbours_save
!-------------------------------------------------------------------------------
